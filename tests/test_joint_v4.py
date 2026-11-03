import csv
import json
import sys
from types import SimpleNamespace
from pathlib import Path

import pytest
from ecd_epitope_scanner.integrations import JsonWorkerRunner, cache_key
from ecd_epitope_scanner.integrations._adapter import request_for
from ecd_epitope_scanner.integrations.base import ToolError
from ecd_epitope_scanner.integrations.spurs_adapter import SPURSAdapter
from ecd_epitope_scanner.joint_design_v4 import run_joint_v4_ml
from ecd_epitope_scanner.joint_design_v4.contexts import write_mapping_sidecar
from ecd_epitope_scanner.joint_design_v4.engine import (
    _adapter,
    _foldx_shortlist,
    _formal_candidate_reasons,
)
from ecd_epitope_scanner.joint_design_v4.evidence import attach_epitope_islands, epitope_evidence
from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence,
    EpitopeEvidence,
    PredictionResult,
    ResidueKey,
    SurfaceEvidence,
    V4Candidate,
    V4MutationOption,
    V4Sector,
    V4SurfaceHypothesis,
)
from ecd_epitope_scanner.joint_design_v4.geometry_evidence import annotate_hypothesis_geometry
from ecd_epitope_scanner.joint_design_v4.hypotheses import hypotheses_from_surface_v2
from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
from ecd_epitope_scanner.joint_design_v4.sectors import build_epitope_sectors
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord
from ecd_epitope_scanner.util import sha256_file


@pytest.fixture(autouse=True)
def supplied_repeat_fixture(monkeypatch):
    from ecd_epitope_scanner.joint_design_v4 import engine
    # These frozen shortlist fixtures exercise the historical five-slot protocol.
    monkeypatch.setitem(engine.config.JOINT_V4_ML_PARAMETERS, "execution_protocol", {"name": "legacy"})
    original, saved = engine._prediction, {}
    def predict(adapter, structure, chain, sequence, output, structure_type, parameters, supplied, key=None, mutations=None):
        source = supplied.get("foldx") if supplied else None
        if key == "foldx" and source and source.status == "success":
            saved.update({tuple(sorted(r["mutations"])): r for r in source.records})
        if key == "foldx" and saved and (source or "foldx_refinement" in str(output)):
            rows = []
            for n, group in enumerate(mutations or []):
                row = dict(saved[tuple(sorted(group))])
                start = int(parameters.get("sampling_slot_start", 0))
                row["paired_models"] = [{"sampling_slot": i, "ddg_kcal_mol": row["ddg"],
                    "batch_path": str(output), "group": n, "run": i}
                    for i in range(start,start+int(parameters.get("number_of_runs", 5)))]
                rows.append(row)
            return PredictionResult("success", records=rows)
        return original(adapter,structure,chain,sequence,output,structure_type,parameters,supplied,key,mutations)
    monkeypatch.setattr(engine,"_prediction",predict)


def _fixture_records(count: int = 12):
    sequence = "DEKRFWYAILMSTNQ"[:count]
    return [
        ResidueRecord(
            uniprot_pos=index,
            aa=aa,
            chain_id="A",
            auth_seq_id=index,
            ca=(float(index * 3.4), 0.0, 0.0),
            sasa=100.0,
            sidechain_sasa=40.0,
            rsasa=0.5,
            surface_accessible=True,
        )
        for index, aa in enumerate(sequence, 1)
    ]


def _protein(records):
    sequence = "".join(item.aa for item in records)
    return ProteinContext(
        accession="TEST",
        protein_name="test",
        organism="test",
        sequence=sequence,
        ecd_start=1,
        ecd_end=len(sequence),
        membrane_anchor_position=len(sequence),
        topology="test",
    )


def test_external_engine_disable_ignores_configured_and_environment_commands(tmp_path, monkeypatch):
    structure = tmp_path / "input.pdb"
    structure.write_text("PDB fixture", encoding="utf-8")
    monkeypatch.setenv("ECD_V4_SPURS_COMMAND", json.dumps(["must-not-run"]))
    adapter = _adapter(
        "spurs_single",
        {"backend": "wsl_gpu", "wsl_python": "/env/bin/python", "command": ["configured"]},
        tmp_path / "cache",
        tool_name="spurs_single",
        disabled=True,
    )

    result = adapter.predict(request_for(
        adapter, structure, "A", "ACD", tmp_path, parameters={"device": "cuda"}
    ))

    assert adapter.command == []
    assert adapter.inventory().available is False
    assert adapter.inventory().failure_reason == "disabled by request"
    assert result.status == "failed"
    assert result.errors[0]["code"] == "TOOL_NOT_CONFIGURED"


def test_json_worker_contract_rejects_stdout_pollution(tmp_path: Path):
    request_code = (
        "import json,sys; print('log'); print(json.dumps({'status':'success','records':[]}))"
    )
    runner = JsonWorkerRunner([sys.executable, "-c", request_code], timeout_seconds=5)
    # The runner must not silently accept progress text on stdout.
    from ecd_epitope_scanner.joint_design_v4.models import PredictionRequest

    try:
        runner.run(
            PredictionRequest(
                schema_version="1.0",
                request_id="x",
                structure_path=__file__,
                chain_id="A",
                sequence="A",
            )
        )
    except Exception as exc:
        assert getattr(exc, "code", None) == "OUTPUT_SCHEMA_INVALID"
    else:
        raise AssertionError("stdout pollution was accepted")


def test_json_worker_reports_timeout_and_crash():
    from ecd_epitope_scanner.joint_design_v4.models import PredictionRequest

    request = PredictionRequest(
        schema_version="1.0", request_id="x", structure_path=__file__, chain_id="A", sequence="A"
    )
    timeout = JsonWorkerRunner(
        [sys.executable, "-c", "import time; time.sleep(1)"], timeout_seconds=0.01
    )
    with pytest.raises(ToolError, match="exceeded") as timeout_error:
        timeout.run(request)
    assert timeout_error.value.code == "TOOL_TIMEOUT"
    crash = JsonWorkerRunner([sys.executable, "-c", "raise SystemExit(3)"], timeout_seconds=5)
    with pytest.raises(ToolError, match="code 3") as crash_error:
        crash.run(request)
    assert crash_error.value.code == "TOOL_CRASHED"


def test_cache_key_includes_structure_content_and_model_revision(tmp_path: Path):
    structure = tmp_path / "same-name.pdb"
    structure.write_text("ATOM A", encoding="ascii")
    first = cache_key(
        structure,
        chain_id="A",
        sequence="A",
        structure_type="solved",
        tool_version="1",
        model_revision="a",
        relevant_parameters={},
    )
    second = cache_key(
        structure,
        chain_id="A",
        sequence="A",
        structure_type="solved",
        tool_version="1",
        model_revision="b",
        relevant_parameters={},
    )
    assert first != second


def test_gpu_request_device_and_multi_mutations_are_cache_scoped(tmp_path: Path):
    structure = tmp_path / "target.pdb"
    structure.write_text("ATOM", encoding="ascii")
    worker = (
        "import json,sys; request=json.load(sys.stdin); "
        "print(json.dumps({'status':'success','records':["
        "{'mutations':item,'ddg':0.1} for item in request['mutations']]}))"
    )
    adapter = SPURSAdapter(
        command=[sys.executable, "-c", worker], timeout_seconds=5, cache_dir=tmp_path / "cache"
    )
    first_request = request_for(
        adapter,
        structure,
        "A",
        "AA",
        tmp_path,
        mutations=[["A1C"]],
        parameters={"device": "cuda"},
    )
    second_request = request_for(
        adapter,
        structure,
        "A",
        "AA",
        tmp_path,
        mutations=[["A2C"]],
        parameters={"device": "cuda"},
    )
    assert first_request.device == "cuda"
    assert adapter.predict(first_request).records[0]["mutations"] == ["A1C"]
    second = adapter.predict(second_request)
    assert second.records[0]["mutations"] == ["A2C"]
    assert second.provenance["cache_hit"] is False
    assert second.provenance["cache_disabled_reason"] == "model_digest_unavailable"
    assert not (tmp_path / "cache").exists()


def test_v4_foldx_repair_cache_is_scoped_by_binary_digest(tmp_path: Path, monkeypatch):
    from subprocess import CompletedProcess

    from ecd_epitope_scanner import stability

    structure = tmp_path / "target.pdb"
    structure.write_text("ATOM", encoding="ascii")
    output = tmp_path / "run" / "output"
    output.mkdir(parents=True)
    calls = []

    def fake_run(command, workdir, timeout):
        calls.append(Path(workdir))
        (Path(workdir) / "structure_Repair.pdb").write_text("REPAIRED", encoding="ascii")
        return CompletedProcess(command, 0)

    monkeypatch.setattr(stability, "_run_engine", fake_run)
    first = stability.repair_with_foldx("foldx", structure, output, binary_sha256="a" * 64)
    assert stability.repair_with_foldx("foldx", structure, output, binary_sha256="a" * 64) == first
    second = stability.repair_with_foldx("foldx", structure, output, binary_sha256="b" * 64)
    assert first != second
    assert len(calls) == 2


def test_spurs_adapter_rejects_incomplete_single_and_multi_outputs(tmp_path: Path):
    from ecd_epitope_scanner.joint_design_v4.models import PredictionRequest

    adapter = SPURSAdapter()
    request = PredictionRequest("1.0", "single", __file__, "A", sequence="AC")
    single = PredictionResult(status="success", records=[{"position": 1, "wt": "A", "scores": {"C": 0.1}}])
    assert adapter.checked_result(request, single).errors[0]["code"] == "RESIDUE_COUNT_MISMATCH"
    multi_request = PredictionRequest(
        "1.0", "multi", __file__, "A", sequence="AC",
        mutations=[["A1D"], ["C2D"]],
    )
    multi = PredictionResult(status="success", records=[{"mutations": ["A1D"], "ddg": 0.1}])
    assert adapter.checked_result(multi_request, multi).errors[0]["code"] == "OUTPUT_MISSING"


def test_discotope_keeps_raw_and_calibrated_scores_and_uses_surface_rank():
    records = _fixture_records(3)
    keys = {
        item.uniprot_pos: ResidueKey("0", "A", item.auth_seq_id, uniprot_pos=item.uniprot_pos, one_letter=item.aa)
        for item in records
    }
    graph = {1: {2: 3.4}, 2: {1: 3.4}}
    rows = [
        {"position": 1, "score": 2.0, "calibrated_score": 0.8, "epitope": True},
        {"position": 2, "score": 1.0, "calibrated_score": 0.2, "epitope": False},
        {"position": 3, "score": 100.0, "calibrated_score": 1.0, "epitope": True},
    ]
    evidence = epitope_evidence(records, keys, rows, graph)
    assert evidence[0].raw_score == 2.0
    assert evidence[0].calibrated_score == 0.8
    assert evidence[0].surface_percentile == 1.0
    assert evidence[0].local_peak
    assert evidence[0].predicted_epitope_label is True
    assert evidence[2].surface_percentile is None


def test_epitope_islands_ignore_scores_outside_surface_graph():
    records = _fixture_records(3)
    evidence = []
    for item in records:
        evidence.append(
            {
                "position": item.uniprot_pos,
                "score": 0.0,
            }
        )
    from ecd_epitope_scanner.joint_design_v4.models import EpitopeEvidence, ResidueKey

    typed = [
        EpitopeEvidence(
            residue=ResidueKey(
                "1", "A", row["position"], one_letter="A", uniprot_pos=row["position"]
            ),
            surface_percentile=1.0,
        )
        for row in evidence
    ]
    attach_epitope_islands(typed, {1: {2: 1.0}, 2: {1: 1.0}}, 0.8, 0.65)
    assert typed[0].island_id is not None
    assert typed[2].island_id is None


def test_v4_preview_is_explicitly_not_formal_and_respects_hard_protection(tmp_path: Path):
    records = _fixture_records()
    records[4].disulfide = True
    records[4].hard_exclusions = ["disulfide-cysteine"]
    success = PredictionResult(status="success", records=[])
    discotope = PredictionResult(
        status="success",
        records=[
            {"position": item.uniprot_pos, "score": 100.0 if item.uniprot_pos == 5 else 0.0}
            for item in records
        ],
    )
    result = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path,
        parameters={"mode": "preview", "max_primary_constructs": 2},
        assay_format="cell_surface",
        prediction_results={"procaliper": success, "spurs": success, "discotope": discotope},
    )
    assert result.recommendation_status == "PREVIEW_NOT_FOR_FORMAL_RECOMMENDATION"
    assert not any(c.impact_evidence.get("policy") == "regional_mechanism_pareto_v1" for c in result.candidates)
    assert any("accessibility is unassessed" in warning for warning in result.warnings)
    protected = [item for item in result.mutation_options if item.residue.uniprot_pos == 5]
    assert protected and all(not item.eligible for item in protected)
    assert all(
        item.stability is None or item.stability.foldx_ddg_kcal_mol is None
        for item in result.candidates
    )
    assert result.legacy_designs
    primary_designs = [
        item for item in result.legacy_designs if item.design_role == "STAGE1_PROBE"
    ]
    stage2_designs = [
        item for item in result.legacy_designs if item.design_role.startswith("STAGE2")
    ]
    assert primary_designs and all(item.status == "PREVIEW_CANDIDATE" for item in primary_designs)
    assert stage2_designs
    assert all(item.status in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"} for item in stage2_designs)
    assert all(
        item.parent_design_id in {primary.id for primary in primary_designs}
        for item in stage2_designs
    )
    assert all("contact" in item.rationale.lower() for item in stage2_designs)
    assert all(
        next(candidate for candidate in result.candidates if candidate.candidate_id == candidate_id).order < 4
        for candidate_id in result.panel.selected_candidate_ids
    )
    assert result.candidates
    assert all(
        {"spurs_multi_ddg_kcal_mol", "foldx_ddg_kcal_mol", "hypothesis_near_glycan_fraction"}
        <= set(candidate.risk_components)
        for candidate in result.candidates
    )
    assert set().union(*(set(item.residue_positions) for item in result.sectors)) == {
        item.uniprot_pos for item in records
    }
    graph_path = tmp_path / "surface_graph_v4.json"
    graph_snapshot = json.loads(graph_path.read_text(encoding="utf-8"))
    manifest = json.loads((tmp_path / "v4_manifest.json").read_text(encoding="utf-8"))
    manifest_v41 = json.loads((tmp_path / "v4_manifest_v41.json").read_text(encoding="utf-8"))
    assert graph_snapshot["position_system"] == "uniprot_pos"
    assert manifest["provenance"]["surface_graph"]["snapshot_sha256"] == sha256_file(graph_path)
    assert manifest["schema_version"] == "4.0"
    assert manifest_v41["schema_version"] == "4.1"
    assert manifest_v41["base_schema_version"] == "4.0"
    assert (tmp_path / "panel_identifiability_v41.json").is_file()


def test_v4_strict_blocks_without_external_evidence(tmp_path: Path):
    records = _fixture_records()
    result = run_joint_v4_ml(records, _protein(records), tmp_path, parameters={"mode": "strict"})
    assert result.recommendation_status == "BLOCKED_EXTERNAL_EVIDENCE"
    assert "Missing external evidence" in " ".join(result.warnings)


def test_v4_combination_disagreement_is_recorded(tmp_path: Path):
    records = _fixture_records(10)
    evidence = _complete_per_residue_evidence(records)
    seed = run_joint_v4_ml(
        records, _protein(records), tmp_path / "seed",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    combinations, foldx_combinations = _combination_rows(seed)
    foldx_combinations = [dict(row, ddg=3.0) for row in foldx_combinations]
    evidence["spurs_multi"] = PredictionResult(status="success", records=combinations)
    evidence["foldx"] = PredictionResult(status="success", records=foldx_combinations)
    result = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path,
        parameters={
            "mode": "strict",
            "regional_recommendation_policy": None,
            "candidate_search": {"orders": [2]},
        },
        prediction_results=evidence,
    )
    assert result.candidates
    assert any(
        item.stability and item.stability.consensus_status == "CONSENSUS_REVIEW"
        for item in result.candidates
    )
    assert result.recommendation_status == "BLOCKED_EXTERNAL_EVIDENCE"
    assert not result.panel.selected_candidate_ids
    assert all(item.status != "RECOMMENDED_LOW_PREDICTED_RISK" for item in result.legacy_designs)


def _complete_per_residue_evidence(records):
    return {
        "procaliper": PredictionResult(
            status="success",
            records=[
                {"position": item.uniprot_pos, "sasa": 100.0, "rsasa": 0.5} for item in records
            ],
        ),
        "spurs": PredictionResult(
            status="success",
            records=[
                {
                    "position": item.uniprot_pos,
                    "wt": item.aa,
                    "scores": {
                        aa: 0.0 for aa in "ACDEFGHIKLMNPQRSTVWY" if aa != item.aa
                    },
                }
                for item in records
            ],
        ),
        "discotope": PredictionResult(
            status="success",
            records=[
                {"position": item.uniprot_pos, "score": float(item.uniprot_pos)} for item in records
            ],
        ),
    }


def _combination_rows(seed):
    # Frozen legacy selection fixtures; the new policy has its own request-driven integration checks.
    def rows(candidates):
        return [
            {"mutations": [option.notation for option in item.mutations], "ddg": 0.1}
            for item in candidates
        ]

    return rows(seed.candidates), rows(_foldx_shortlist(seed.candidates, 42))


def test_v4_strict_recommends_only_candidate_level_consensus_pass(tmp_path: Path):
    records = _fixture_records(10)
    evidence = _complete_per_residue_evidence(records)
    seed = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path / "seed",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    assert seed.candidates
    assert {item.recommendation_class for item in seed.candidates} == {
        "EPITOPE_PRIORITIZED",
        "STABILITY_FALLBACK",
    }
    combinations, foldx_combinations = _combination_rows(seed)
    strict_evidence = dict(evidence)
    strict_evidence["spurs_multi"] = PredictionResult(status="success", records=combinations)
    strict_evidence["foldx"] = PredictionResult(status="success", records=foldx_combinations)
    result = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path / "strict",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=strict_evidence,
    )
    assert result.recommendation_status == "FORMAL_RECOMMENDATION"
    assert result.panel.selected_candidate_ids
    selected = {
        item.candidate_id: item
        for item in result.candidates
        if item.candidate_id in result.panel.selected_candidate_ids
    }
    assert selected and all(item.formal_eligible for item in selected.values())
    assert all(
        item.stability and item.stability.consensus_status == "CONSENSUS_PASS"
        for item in selected.values()
    )
    primary_designs = [
        item for item in result.legacy_designs if item.design_role == "STAGE1_PROBE"
    ]
    stage2_designs = [
        item for item in result.legacy_designs if item.design_role.startswith("STAGE2")
    ]
    assert primary_designs
    assert {item.status for item in primary_designs} == {"RECOMMENDED_LOW_PREDICTED_RISK"}
    assert stage2_designs
    assert all(item.status in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"} for item in stage2_designs)
    manifest = json.loads((tmp_path / "strict" / "v4_manifest.json").read_text(encoding="utf-8"))
    assert manifest["worker_runs"]["spurs_multi"]["status"] == "success"
    assert manifest["software"]["source_code_sha256"]
    assert manifest["normalized_structure_sha256"] is None


@pytest.mark.parametrize("failed_tool", ["spurs_multi", "foldx"])
def test_v4_strict_blocks_failed_combination_worker_even_with_scores(tmp_path: Path, failed_tool: str):
    records = _fixture_records(10)
    evidence = _complete_per_residue_evidence(records)
    seed = run_joint_v4_ml(
        records, _protein(records), tmp_path / "seed",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    combinations, foldx_combinations = _combination_rows(seed)
    evidence["spurs_multi"] = PredictionResult(status="success", records=combinations)
    evidence["foldx"] = PredictionResult(status="success", records=foldx_combinations)
    evidence[failed_tool] = PredictionResult(
        status="failed", records=(foldx_combinations if failed_tool == "foldx" else combinations),
        errors=[{"code": "TOOL_CRASHED", "message": "simulated failure"}],
    )
    result = run_joint_v4_ml(
        records, _protein(records), tmp_path / "failed",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    assert result.recommendation_status == "BLOCKED_EXTERNAL_EVIDENCE"
    assert not result.panel.selected_candidate_ids
    assert result.errors and result.errors[0]["code"] == "TOOL_CRASHED"


def test_v4_panel_coverage_is_union_over_all_sectors():
    left = _series_candidate("VS001", 0, "EPITOPE_PRIORITIZED")
    right = _series_candidate("VS002", 1, "STABILITY_FALLBACK")
    left.covered_positions_by_scenario = {"geodesic_12": [1, 2]}
    right.covered_positions_by_scenario = {"geodesic_12": [3]}
    sectors = [
        V4Sector("VS001", [1, 2], 1, 20.0, 20.0, 0.9),
        V4Sector("VS002", [3, 4], 1, 20.0, 20.0, 0.1),
    ]
    panel, selected = select_panel(
        [left, right], max_primary=2, sectors=sectors,
        surface_areas={1: 10.0, 2: 10.0, 3: 10.0, 4: 10.0},
    )
    assert len(selected) == 2
    assert panel.coverage_denominator_residues == 4
    assert panel.covered_residues_by_scenario["geodesic_12"] == 3
    assert panel.achieved_coverage["geodesic_12"] == 0.75
    assert panel.achieved_area_coverage["geodesic_12"] == 0.75


def _probe(candidate_id: str, order: int, positions: tuple[int, ...], response: dict[str, float], exploratory=False):
    mutations = [
        V4MutationOption(
            residue=ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D"),
            mutant_aa="A", notation=f"D{position}A", chemistry_class="charge_removal",
            spurs_ddg_kcal_mol=0.0, score=1.0, eligible=True,
        )
        for position in positions
    ]
    return V4Candidate(
        candidate_id=candidate_id, sector_id="VS001", anchor_region_id="VS001",
        mutations=mutations, order=order, hypothesis_response=response,
        hypothesis_perturbation=dict(response),
        search_role="EXPLORATORY" if exploratory else "PRIMARY_POOL",
    )


def test_v41_identifiability_panel_never_selects_exploratory_four_site_probes():
    hypotheses = [
        V4SurfaceHypothesis("H1", {1: 1.0}, discotope_prior=0.8),
        V4SurfaceHypothesis("H2", {2: 1.0}, discotope_prior=0.2),
    ]
    probes = [
        _probe("double", 2, (1, 2), {"H1": 0.45, "H2": 0.45}),
        _probe("triple", 3, (3, 4, 5), {"H1": 1.0, "H2": 0.0}),
        _probe("z_duplicate", 3, (10, 11, 12), {"H1": 1.0, "H2": 0.0}),
        _probe("exploratory-four", 4, (6, 7, 8, 9), {"H1": 0.0, "H2": 1.0}, True),
    ]
    panel, selected = select_panel(probes, max_primary=3, hypotheses=hypotheses)
    assert "triple" in panel.selected_candidate_ids
    assert len({"triple", "z_duplicate"} & set(panel.selected_candidate_ids)) == 1
    assert all(item.order < 4 and item.search_role != "EXPLORATORY" for item in selected)
    assert panel.prior_weighted_hypothesis_response > panel.mean_hypothesis_response
    assert set(panel.selected_candidate_ids) == set(panel.marginal_information_gain)

    only_exploratory = [probes[-1]]
    panel, selected = select_panel(only_exploratory, max_primary=3, hypotheses=hypotheses)
    assert selected == []
    assert panel.selected_candidate_ids == []


def test_v41_double_and_triple_search_budgets_are_independent():
    positions = list(range(1, 9))
    sector = V4Sector("VS001", positions, 1, 800.0, 800.0, 0.5)
    options = [
        V4MutationOption(
            residue=ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D"),
            mutant_aa="A", notation=f"D{position}A", chemistry_class="charge_removal",
            spurs_ddg_kcal_mol=0.0, score=1.0, eligible=True,
            accessibility_score=1.0, chemical_disruption_score=1.0,
        )
        for position in positions
    ]
    graph = {
        position: {
            neighbor: 3.0 for neighbor in (position - 1, position + 1)
            if neighbor in positions
        }
        for position in positions
    }
    candidates = generate_candidates(
        [sector], options, graph,
        {
            "orders": [2, 3], "exploratory_orders": [],
            "site_pool_per_sector": 8, "shortlist_per_sector": 12,
            "beam_width_by_order": {2: 1, 3: 100},
            "shortlist_by_order": {2: 2, 3: 3},
            "maximum_geodesic_span_a": 24.0,
        },
    )
    by_order = {order: sum(item.order == order for item in candidates) for order in (2, 3)}
    assert by_order == {2: 1, 3: 3}
    assert all(item.anchor_region_id == "VS001" for item in candidates)


def test_v41_glycan_evidence_separates_mutation_prohibition_from_nearby_shielding():
    records = [
        ResidueRecord(1, "N", "A", 1, glycosylation_site=True),
        ResidueRecord(2, "A", "A", 2, glycan_shielded=True),
    ]
    footprint = SimpleNamespace(
        footprint_id="H1", center_vertex_id=0, center_residue_id=1,
        generation_mode="adaptive_area", radius=None, target_area=600.0,
        residues=[1, 2], residue_surface_area={1: 10.0, 2: 30.0},
        residue_surface_fraction={1: 1.0, 2: 1.0}, total_surface_area=40.0,
        quality_flags=[], region_id="R1", mutable_residues=[1, 2],
        protected_residues=[],
    )
    hypotheses = hypotheses_from_surface_v2(
        SimpleNamespace(raw_footprints=[], footprints=[footprint]), records
    )
    glycan = hypotheses[0].glycan_evidence
    assert glycan["mutation_prohibition_positions"] == [1]
    assert glycan["near_resolved_glycan_positions"] == [2]
    assert glycan["direct_glycosylation_area_fraction"] == 0.25
    assert glycan["near_resolved_glycan_area_fraction"] == 0.75


def test_v41_geometry_keeps_relative_anchor_axis_distinct_from_membrane_plane():
    records = [
        ResidueRecord(1, "A", "A", 1, ca=(0.0, 0.0, 5.0)),
        ResidueRecord(2, "A", "A", 2, ca=(0.0, 0.0, 10.0)),
        ResidueRecord(3, "A", "A", 3, ca=(0.0, 0.0, 100.0)),
    ]
    relative = V4SurfaceHypothesis("H1", {1: 1.0, 2: 1.0}, mean_normal=(0.0, 0.0, 1.0))
    annotate_hypothesis_geometry(
        [relative], records,
        assay_orientation_geometry={
            "status": "relative_axis_only", "axis_source": "C-terminal anchor surrogate",
            "anchor_xyz": [0.0, 0.0, 0.0], "extracellular_axis": [0.0, 0.0, 1.0],
        },
    )
    assert relative.membrane_geometry["relative_axis_projection_from_anchor_a"] == 7.5
    assert relative.membrane_geometry["weighted_plane_distance_a"] is None
    assert relative.membrane_geometry["weighted_accessibility_proxy"] is None
    assert relative.membrane_geometry["physical_membrane_plane_inferred"] is False

    plane = V4SurfaceHypothesis("H2", {1: 1.0, 2: 1.0}, mean_normal=(0.0, 0.0, 1.0))
    annotate_hypothesis_geometry([plane], records, membrane_anchor=(0.0, 0.0, 0.0), membrane_normal=(0.0, 0.0, 1.0))
    assert plane.membrane_geometry["weighted_plane_distance_a"] == 7.5
    assert plane.membrane_geometry["membrane_plane_surrogate_available"] is True
    assert plane.membrane_geometry["physical_membrane_plane_inferred"] is False
    assert plane.membrane_geometry["weighted_membrane_position_factor_proxy"] is not None
    assert "not molecular or assay accessibility" in plane.membrane_geometry["accessibility_proxy_interpretation"]


def test_v4_formal_candidate_rejects_surface_and_sector_qc_failures():
    candidate = _series_candidate("VS001", 0, "STABILITY_FALLBACK")
    reasons = _formal_candidate_reasons(
        candidate, allow_spurs_only=False, minimum_coverage_gain=0.05,
        lower_order_coverage=0.0, sasa_conflict_positions={1},
        sasa_context_sensitive_positions={2},
        sasa_method_unverified_positions={1},
        sector_quality_flags={"DISCONNECTED_SECTOR"},
    )
    assert "SASA_SURFACE_CLASSIFICATION_DISAGREEMENT" not in reasons
    assert "SASA_STRUCTURE_CONTEXT_CLASSIFICATION_CHANGED" in reasons
    assert "SASA_METHOD_SETTINGS_UNVERIFIED" not in reasons
    assert "SECTOR_GEOMETRY_QC_FAILED" in reasons


@pytest.mark.parametrize("gain,expected", [
    (None, "HIGHER_ORDER_COVERAGE_BASELINE_MISSING"),
    (float("nan"), "HIGHER_ORDER_COVERAGE_BASELINE_MISSING"),
    (0.035, "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN"),
    (0.086, None),
])
def test_formal_higher_order_uses_regional_gain_instead_of_global_response(gain, expected):
    candidate = _series_candidate("VS001", 0, "STABILITY_FALLBACK")
    candidate.order = 3
    reasons = _formal_candidate_reasons(candidate, allow_spurs_only=False,
        minimum_coverage_gain=0.05, lower_order_coverage=1.0,
        higher_order_information_gain=0.001,
        regional_coverage_gain=gain, regional_gain_required=True)
    if expected:
        assert expected in reasons
    else:
        assert "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN" not in reasons
        assert "HIGHER_ORDER_COVERAGE_BASELINE_MISSING" not in reasons


@pytest.mark.parametrize("model", ["spurs_multi_ddg_kcal_mol", "foldx_ddg_kcal_mol"])
@pytest.mark.parametrize("missing", [None, float("nan"), float("inf")])
def test_formal_panel_requires_both_finite_models_despite_cached_pass(model, missing):
    candidate = _series_candidate("VS001", 0, "STABILITY_FALLBACK")
    candidate.stability.consensus_status = "CONSENSUS_PASS"
    setattr(candidate.stability, model, missing)
    reasons = _formal_candidate_reasons(candidate, allow_spurs_only=True,
        minimum_coverage_gain=0.05, lower_order_coverage=0.0)
    name = "FOLDX" if model == "foldx_ddg_kcal_mol" else "SPURS"
    assert f"COMBINATION_{name}_MISSING" in reasons


@pytest.mark.parametrize("mutant", ["C", "P", "G"])
def test_authoritative_mutation_policy_blocks_new_and_cached_options(mutant):
    from ecd_epitope_scanner.joint_design_v4.evidence import build_single_evidence
    from ecd_epitope_scanner.joint_design_v4.contexts import residue_keys
    records = _fixture_records(1)
    evidence = build_single_evidence(records, residue_keys(records),
        [{"position": 1, "mutant": mutant, "ddg": -2.0}], {})
    option = next(row for row in evidence if row.mutant_aa == mutant)
    assert not option.eligible
    assert "auto-excluded-mutant" in option.hard_exclusions
    candidate = _series_candidate("VS001", 0, "STABILITY_FALLBACK")
    candidate.mutations[0].mutant_aa = mutant
    candidate.mutations[0].hard_exclusions = []  # Stale cache claims eligibility.
    reasons = _formal_candidate_reasons(candidate, allow_spurs_only=True,
        minimum_coverage_gain=0.05, lower_order_coverage=0.0)
    assert "AUTO_EXCLUDED_MUTANT" in reasons


def test_v4_procaliper_is_diagnostic_and_does_not_block_projection(tmp_path: Path):
    records = _fixture_records(10)
    evidence = _complete_per_residue_evidence(records)
    seed = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path / "seed",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    combinations, foldx_combinations = _combination_rows(seed)
    evidence["procaliper"] = PredictionResult(status="failed")
    evidence["spurs_multi"] = PredictionResult(status="success", records=combinations)
    evidence["foldx"] = PredictionResult(status="success", records=foldx_combinations)
    result = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path / "blocked",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    assert result.recommendation_status == "FORMAL_RECOMMENDATION"
    assert result.legacy_designs
    for row in result.panel.sector_recommendations:
        for recommendation in [row, *row.get("alternatives", [])]:
            if recommendation.get("candidate_id"):
                assert recommendation["validation_status"]["foldx"] != "UNKNOWN"
    primary_designs = [
        item for item in result.legacy_designs if item.design_role == "STAGE1_PROBE"
    ]
    stage2_designs = [
        item for item in result.legacy_designs if item.design_role.startswith("STAGE2")
    ]
    assert primary_designs
    assert all(item.status == "RECOMMENDED_LOW_PREDICTED_RISK" for item in primary_designs)
    assert any("ProCaliper SASA cross-check is unavailable" in item for item in result.warnings)
    assert all(item.status in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"} for item in stage2_designs)


def test_mapping_sidecar_has_explicit_normalized_position(tmp_path: Path):
    records = _fixture_records(3)
    keys = {
        item.uniprot_pos: ResidueKey(
            model_id="0",
            chain_id=item.chain_id,
            auth_seq_id=item.auth_seq_id,
            uniprot_pos=item.uniprot_pos,
            one_letter=item.aa,
        )
        for item in records
    }
    path = write_mapping_sidecar(tmp_path / "mapping.csv", records, keys)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [int(row["normalized_position"]) for row in rows] == [1, 2, 3]


def test_foldx_shortlist_is_capped_and_sector_balanced():
    candidates = [
        V4Candidate(
            candidate_id=f"{sector}-{index}",
            sector_id=sector,
            mutations=[],
            order=2,
            epitope_score=10.0 - index,
            structural_risk=float(index),
        )
        for sector in ("VS001", "VS002", "VS003")
        for index in range(3)
    ]
    selected = _foldx_shortlist(candidates, 5)
    assert len(selected) == 5
    assert [item.sector_id for item in selected[:3]] == ["VS001", "VS002", "VS003"]
    assert [item.candidate_id for item in selected] == [
        "VS001-0",
        "VS002-0",
        "VS003-0",
        "VS001-1",
        "VS002-1",
    ]


def test_foldx_shortlist_supports_configured_budget_of_60():
    candidates = [
        V4Candidate(
            candidate_id=f"VS{index // 10 + 1:03d}-{index}",
            sector_id=f"VS{index // 10 + 1:03d}",
            mutations=[],
            order=2,
            structural_risk=float(index),
        )
        for index in range(75)
    ]
    selected = _foldx_shortlist(candidates, 60)
    assert len(selected) == 60
    assert len({candidate.candidate_id for candidate in selected}) == 60


def _series_candidate(sector: str, index: int, recommendation_class: str) -> V4Candidate:
    positions = (index * 2 + 1, index * 2 + 2)
    options = [
        V4MutationOption(
            residue=ResidueKey(
                "0", "A", position, uniprot_pos=position, one_letter="D"
            ),
            mutant_aa="A",
            notation=f"D{position}A",
            chemistry_class="charge_removal",
            spurs_ddg_kcal_mol=0.0,
            score=1.0,
            eligible=True,
            epitope_prior=0.9 if recommendation_class == "EPITOPE_PRIORITIZED" else 0.1,
            epitope_supported=recommendation_class == "EPITOPE_PRIORITIZED",
        )
        for position in positions
    ]
    return V4Candidate(
        candidate_id=f"{sector}-{recommendation_class}-{index}",
        sector_id=sector,
        mutations=options,
        order=2,
        epitope_score=sum(item.epitope_prior for item in options),
        structural_risk=float(index) / 10.0,
        coverage_by_scenario={"geodesic_12": 0.5},
        stability=CombinationStabilityEvidence(
            mutations=tuple(item.notation for item in options),
            spurs_multi_ddg_kcal_mol=0.1,
            foldx_ddg_kcal_mol=0.2,
            consensus_status="CONSENSUS_PASS",
        ),
        formal_eligible=True,
        epitope_supported_sites=(2 if recommendation_class == "EPITOPE_PRIORITIZED" else 0),
        recommendation_class=recommendation_class,
    )


def test_v4_panel_returns_epitope_and_stability_series_per_sector():
    candidates = [
        _series_candidate(sector, index, recommendation_class)
        for sector in ("VS001", "VS002")
        for recommendation_class in ("EPITOPE_PRIORITIZED", "STABILITY_FALLBACK")
        for index in range(4)
    ]
    panel, selected = select_panel(
        candidates,
        max_primary=12,
        recommendations_per_sector=6,
        epitope_priority_per_sector=3,
        stability_fallback_per_sector=3,
    )
    assert len(selected) == 12
    for sector in ("VS001", "VS002"):
        series = [item for item in selected if item.sector_id == sector]
        assert len(series) == 6
        assert sum(item.recommendation_class == "EPITOPE_PRIORITIZED" for item in series) == 3
        assert sum(item.recommendation_class == "STABILITY_FALLBACK" for item in series) == 3
    assert len(panel.candidates_by_role["EPITOPE_PRIORITIZED"]) == 6
    assert len(panel.candidates_by_role["STABILITY_FALLBACK"]) == 6


def test_v4_discotope_failure_allows_audited_stability_fallback(tmp_path: Path):
    records = _fixture_records(10)
    evidence = _complete_per_residue_evidence(records)
    evidence["discotope"] = PredictionResult(status="failed")
    seed = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path / "fallback_seed",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    assert seed.candidates
    combinations, foldx_combinations = _combination_rows(seed)
    evidence["spurs_multi"] = PredictionResult(status="success", records=combinations)
    evidence["foldx"] = PredictionResult(status="success", records=foldx_combinations)
    result = run_joint_v4_ml(
        records,
        _protein(records),
        tmp_path / "fallback_final",
        parameters={"mode": "strict", "regional_recommendation_policy": None, "candidate_search": {"orders": [2]}},
        prediction_results=evidence,
    )
    assert result.recommendation_status == "FORMAL_RECOMMENDATION"
    assert result.panel.selected_candidate_ids
    selected_ids = set(result.panel.selected_candidate_ids)
    assert all(
        item.recommendation_class == "STABILITY_FALLBACK"
        for item in result.candidates
        if item.candidate_id in selected_ids
    )
    assert any("DiscoTope evidence is unavailable" in item for item in result.warnings)


def test_v4_sector_k_center_avoids_one_residue_island_sectors():
    records = _fixture_records(14)
    keys = {
        item.uniprot_pos: ResidueKey(
            "0", "A", item.auth_seq_id, uniprot_pos=item.uniprot_pos, one_letter=item.aa
        )
        for item in records
    }
    graph = {position: {} for position in keys}
    for position in range(1, 14):
        graph[position][position + 1] = 3.4
        graph[position + 1][position] = 3.4
    surface = [
        SurfaceEvidence(residue=keys[position], mapped_surface_area_a2=80.0)
        for position in keys
    ]
    epitope = [
        EpitopeEvidence(
            residue=keys[position],
            surface_percentile=1.0 if position in {1, 2} else 0.1,
            island_id=f"EI{position:03d}" if position in {1, 2} else None,
        )
        for position in keys
    ]
    sectors = build_epitope_sectors(
        surface,
        epitope,
        graph,
        {item.uniprot_pos: item for item in records},
        max_sectors=3,
        parameters={"minimum_sector_residues": 4, 'policy': 'connected_bounded_v3',
                    'maximum_geodesic_diameter_a': 26.0, 'maximum_euclidean_diameter_a': 18.0},
    )
    assert len(sectors) == 3
    assert min(len(sector.residue_positions) for sector in sectors) >= 4
    assert set().union(*(set(sector.residue_positions) for sector in sectors)) == set(keys)
