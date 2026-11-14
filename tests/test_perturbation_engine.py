import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.joint_design_v4 import engine
from ecd_epitope_scanner.joint_design_v4.models import PredictionResult, V4SurfaceHypothesis
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord
from ecd_epitope_scanner.util import safe_name


def test_default_execution_protocol_matches_published_site_first_v2():
    from ecd_epitope_scanner import config

    assert (
        config.JOINT_V4_ML_PARAMETERS["execution_protocol"]["name"]
        == "ordinary_site_first_v2"
    )


@pytest.mark.parametrize("offset,target", [(0, "antigen"), (100, "antigen"), (100, "binder")])
@pytest.mark.parametrize("protocol", ["ordinary_three_principle_adaptive_v1", "ordinary_site_first_v2", "ordinary_site_first_perturbation_v1"])
def test_adaptive_engine_maps_worker_positions_and_exports_only_qualified_winners(tmp_path, monkeypatch, offset, target, protocol):
    sequence = "DEKRFWYAIL"
    records = [ResidueRecord(uniprot_pos=offset+i, aa=aa, chain_id="A", auth_seq_id=offset+i,
        ca=(i*3.4, 0., 0.), sasa=100., sidechain_sasa=40., rsasa=.5, surface_accessible=True,
        sidechain_exposure={'status': 'KNOWN', 'value': .4, 'numerator_a2': 40.,
                            'denominator_a2': 100., 'method': 'fixture_sidechain_sasa'})
        for i, aa in enumerate(sequence, 1)]
    protein = ProteinContext("TEST", "test", "test", "A"*offset+sequence,
        offset+1, offset+len(sequence), offset+len(sequence), "test")
    calls = []

    def predict(adapter, structure, chain, seq, output, structure_type, parameters, supplied, key=None, mutations=None):
        key = key or adapter.name
        if key in {"spurs_multi", "foldx"}:
            calls.append((key, mutations, Path(output)))
            assert all(1 <= int(value[1:-1]) <= len(sequence) for group in mutations for value in group)
            rows = [{"mutations": group, "ddg": .1} for group in mutations]
            if key == "foldx":
                Path(output).mkdir(parents=True, exist_ok=True)
                for row in rows:
                    row["paired_models"] = [{"sampling_slot": i, "ddg_kcal_mol": .1,
                        "batch_path": str(output), "group": 0, "run": i}
                        for i in range(parameters["sampling_slot_start"],
                                       parameters["sampling_slot_start"]+parameters["number_of_runs"])]
                return PredictionResult("success", rows,
                    provenance={"repair_seconds": 1.25,
                                "buildmodel_seconds": 2.5,
                                "repair_executed": True},
                    runtime_seconds=3.75)
            return PredictionResult("success", rows)
        if key == "spurs":
            return PredictionResult("success", [{"position": i, "wt": aa, "scores": {
                mutant: .1 if mutant == ("A" if aa != "A" else "D") else 3.
                for mutant in "ACDEFGHIKLMNPQRSTVWY" if mutant != aa}}
                for i, aa in enumerate(sequence, 1)])
        if key == "discotope":
            return PredictionResult("success", [{"position": i, "score": float(i)} for i in range(1,len(sequence)+1)])
        return PredictionResult("success", [{"position": i, "sasa": 100., "rsasa": .5}
            for i in range(1,len(sequence)+1)])

    monkeypatch.setattr(engine, "_prediction", predict)
    surface_v2 = None
    if protocol in {"ordinary_site_first_v2", "ordinary_site_first_perturbation_v1"}:
        surface_v2 = SimpleNamespace(surface_regions=[], raw_footprints=[], footprints=[])
        monkeypatch.setattr(engine, "hypotheses_from_surface_v2", lambda *_: [
            V4SurfaceHypothesis("fixture", {record.uniprot_pos: 1. for record in records})
        ])
    result = engine.run_joint_v4_ml(records, protein, tmp_path,
        source_provenance={"analysis_target": target}, disable_external_engines=True,
        parameters={"execution_protocol": protocol, "use_hypothesis_hierarchy": False},
        surface_v2_result=surface_v2)
    assert calls and any(key == "foldx" for key, _, _ in calls)
    assert result.panel.selected_candidate_ids
    assert len(result.panel.selected_candidate_ids) <= 7
    partition = json.loads((tmp_path / "partition_diagnostics_v4.json").read_text(encoding="utf-8"))
    assert "panel_budget_after_sectors" not in partition
    selected = [c for c in result.candidates if c.candidate_id in result.panel.selected_candidate_ids]
    assert all(c.formal_eligible and c.stability.consensus_status == "CONSENSUS_PASS" for c in selected)
    if protocol in {"ordinary_site_first_v2", "ordinary_site_first_perturbation_v1"}:
        assert result.panel.mean_hypothesis_response > 0
        assert any(value > 0 for signature in result.panel.hypothesis_signatures.values() for value in signature)
    assert all(int(value[1:-1]) > offset for c in selected for value in c.stability.mutations)
    ledger = json.loads((tmp_path / "ordinary_candidate_ledger.json").read_text())
    assert ledger["foldx"]["attempted_slots"] <= 300
    foldx_calls = [output for key, _, output in calls if key == "foldx"]
    scored_ids = {safe_name(row["candidate_id"]) for row in ledger["foldx_results"]}
    assert foldx_calls
    if protocol in {"ordinary_site_first_v2", "ordinary_site_first_perturbation_v1"}:
        assert result.execution_protocol['name'] == protocol
        assert result.execution_protocol['representative_policy'] == 'role_ladder_3_6_7_v1'
        assert result.execution_protocol['foldx_site_exploration_policy'] == 'uncovered_residue_before_retry_v1'
        assert result.execution_protocol['foldx_site_exploration_reserve'] == 1
        assert all(path.parent.name == 'fx_batches' for path in foldx_calls)
        assert not (tmp_path/'ordinary_candidate_index.sqlite3').exists()
        assert any(len(groups) > 1 for key, groups, _ in calls if key == 'foldx')
        manifest = json.loads((tmp_path / "v4_manifest.json").read_text(encoding="utf-8"))
        foldx_worker = manifest["worker_runs"]["foldx"]
        assert foldx_worker["runtime_seconds"] == pytest.approx(3.75 * len(foldx_calls))
        assert foldx_worker["provenance"]["worker_metrics"] == {
            "batches": len(foldx_calls),
            "repair_executions": len(foldx_calls),
            "repair_seconds": pytest.approx(1.25 * len(foldx_calls)),
            "buildmodel_seconds": pytest.approx(2.5 * len(foldx_calls)),
        }
        from ecd_epitope_scanner.joint_design_v4.impact_completeness import write_card_acceptance
        analysis = {'execution_protocol': result.execution_protocol,
            'design_summary': {'sector_recommendations': result.panel.sector_recommendations}}
        receipt = write_card_acceptance(tmp_path, analysis)
        missing = sorted({field for row in receipt['cards'] for field in row['missing_fields']})
        assert receipt['status'] == 'PASS', missing
    else:
        assert all(path.parent.parent.name == "foldx_adaptive" and path.parent.name in scored_ids
                   for path in foldx_calls)
    assert len({path.parent for path in foldx_calls}) == len({path.parent.name for path in foldx_calls})
    recommendations = result.panel.sector_recommendations
    manifest = json.loads((tmp_path / "v4_manifest.json").read_text(encoding="utf-8"))
    exported = json.loads((tmp_path / "sector_recommendations_v4.json").read_text(encoding="utf-8"))
    assert exported == recommendations
    for row in recommendations:
        for card in [row, *row.get("alternatives", [])]:
            provenance = card.get("provenance")
            if not provenance:
                continue
            source = provenance["source"]
            assert source["ordinary_adaptive_ledger_ref"] == "ordinary_candidate_ledger.json"
            assert "ordinary_adaptive_ledger" not in source
            for name, model in provenance["models"].items():
                assert model["full_provenance_ref"] == f"v4_manifest.json#/worker_runs/{name}"
                worker = manifest["worker_runs"][name]
                assert model["status"] == worker["status"]
                assert model["provenance"] == {key: value for key, value in worker["provenance"].items()
                                                if key != "adaptive_ledger"}
    recommended_principles = {principle for row in recommendations
                              for recommendation in [row, *row.get("alternatives", [])]
                              for principle in recommendation.get("principles", [])}
    expected = {"mechanism_pareto", "sector_coverage"}
    if target == "binder":
        assert all(q["state"] == "not_applicable" for q in ledger["queues"] if q["principle"] == "epitope_support")
        assert recommended_principles == expected
        for row in recommendations:
            assert any(item["principle"] == "epitope_support" and item["state"] == "not_applicable"
                       for item in row["unavailable_principles"])
    else:
        expected.add("epitope_support")
        assert recommended_principles == expected
    queue_winners = {(str(q["sector_id"]), q["principle"]): next(
        (item["candidate_id"] for item in q["completed"]
         if item["outcome"] == "complete" or (isinstance(item["outcome"], dict)
                                                and item["outcome"].get("outcome") == "complete")), None)
        for q in ledger["queues"]}
    for row in recommendations:
        for recommendation in [row, *row.get("alternatives", [])]:
            if recommendation.get("candidate_id") is not None:
                assert recommendation["site_details"]
                assert recommendation["gate_checks"]
                assert recommendation["regional_geometry_support"]
            for principle in recommendation.get("principles", []):
                assert recommendation["candidate_id"] == queue_winners[(str(row["sector_id"]), principle)]
