from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.integrations import wsl_bridge
from ecd_epitope_scanner.integrations.spurs_adapter import SPURSAdapter
from ecd_epitope_scanner.integrations.workers import procaliper_worker
from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
from ecd_epitope_scanner.joint_design_v4.diagnostics import recommendation_diagnostics
from ecd_epitope_scanner.joint_design_v4.engine import _prediction
from ecd_epitope_scanner.joint_design_v4.models import (
    PredictionResult, ResidueKey, V4MutationOption, V4Sector, V4SurfaceHypothesis,
)


def test_real_ctla4_overlapping_footprints_have_global_construct_identity():
    fixture = json.loads((Path(__file__).parent / "fixtures" / "ctla4_v41_candidate_overlap.json").read_text(encoding="utf-8"))
    graph = {node["uniprot_pos"]: {} for node in fixture["graph_nodes"]}
    for edge in fixture["graph_edges"]:
        left, right = edge["source"], edge["target"]
        graph[left][right] = graph[right][left] = edge["distance_a"]
    sectors = [V4Sector(**row) for row in fixture["sectors"]]
    hypotheses = []
    for row in fixture["hypotheses"]:
        row["residue_weights"] = {int(k): v for k, v in row["residue_weights"].items()}
        hypotheses.append(V4SurfaceHypothesis(**row))
    options = []
    for row in fixture["mutation_options"]:
        row["residue"] = ResidueKey(**row["residue"])
        options.append(V4MutationOption(**row))
    candidates = generate_candidates(sectors, options, graph, {
        "orders": [2, 3], "shortlist_per_order": 10,
    }, hypotheses=hypotheses)
    signatures = [tuple(sorted(option.notation for option in c.mutations)) for c in candidates]
    assert len(signatures) == len(set(signatures))
    triple = next(c for c in candidates if c.order == 3)
    assert {option.notation for option in triple.mutations} == {"K65W", "A66W", "N93C"}
    assert set(triple.targeted_hypotheses) == {"C074_A0600", "C033_A0600"}


def test_inference_boundary_scores_each_physical_combination_once(tmp_path):
    captured = []

    def checked(request, result):
        captured.extend(request.mutations)
        return SPURSAdapter().checked_result(request, result)

    adapter = SimpleNamespace(name="spurs", checked_result=checked)
    result = _prediction(
        adapter, None, "A", "AC", tmp_path, "solved", {},
        {"spurs_multi": PredictionResult("success", records=[{"mutations": ["A1D", "C2D"], "ddg": 0.1}])},
        "spurs_multi", [["A1D", "C2D"], ["C2D", "A1D"]],
    )
    assert captured == [["A1D", "C2D"]]
    assert result.status == "success"


def test_wsl_access_denied_on_utf16_stdout_is_actionable(monkeypatch):
    def denied(*args, **kwargs):
        return subprocess.CompletedProcess(args[0], 4294967295,
            "拒绝访问。错误代码: Wsl/E_ACCESSDENIED\r\n".encode("utf-16-le"), b"")

    monkeypatch.setattr(wsl_bridge.subprocess, "run", denied)
    with pytest.raises(wsl_bridge.WSLBridgeError) as error:
        wsl_bridge._wsl_path(r"E:\user\New project\input.pdb", "Ubuntu-22.04")
    assert error.value.code == "WSL_ACCESS_DENIED"
    assert "Wsl/E_ACCESSDENIED" in str(error.value)
    assert "integrity label" in str(error.value)
    assert error.value.details["stdout"].startswith("拒绝访问")


def test_wsl_exec_preserves_spaces_and_rejects_empty_paths(monkeypatch):
    calls = []

    def converted(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, b"/mnt/e/New project/input.pdb\n", b"")

    monkeypatch.setattr(wsl_bridge.subprocess, "run", converted)
    path = r"E:\New project\input.pdb"
    assert wsl_bridge._wsl_path(path, "Ubuntu-22.04") == "/mnt/e/New project/input.pdb"
    assert calls[0][-1] == path
    assert "--exec" in calls[0]
    monkeypatch.setattr(wsl_bridge.subprocess, "run", lambda *a, **kw: subprocess.CompletedProcess(a[0], 0, b"", b""))
    with pytest.raises(RuntimeError, match="invalid Linux path"):
        wsl_bridge._wsl_path(path, "Ubuntu-22.04")


def test_wsl_mixed_stderr_preserves_worker_traceback():
    raw = "wsl: 代理警告\r\n".encode("utf-16-le") + b"Traceback: model weights unavailable\n"
    assert wsl_bridge._decode_output(raw) == "wsl: 代理警告\nTraceback: model weights unavailable"


def test_procaliper_worker_uses_available_biopython_select():
    from Bio.PDB import Select

    assert issubclass(procaliper_worker._ProteinChainSelect, Select)


def test_empty_panel_distinguishes_missing_single_evidence_from_rejection():
    failed = PredictionResult("failed", errors=[{"code": "WSL_ACCESS_DENIED", "message": "identity rejected"}])
    diagnostics = recommendation_diagnostics([], [], {"spurs_single": failed}, mutation_options=[])
    assert diagnostics["outcome"] == "SINGLE_STABILITY_EVIDENCE_UNAVAILABLE"
    assert diagnostics["worker_failures"]["spurs_single"][0]["code"] == "WSL_ACCESS_DENIED"
    assert diagnostics["formal_exclusion_counts"] == {}


def test_stability_rejection_is_not_misreported_as_a_foldx_crash():
    candidate = SimpleNamespace(
        mutations=[SimpleNamespace(notation="A1D")], formal_eligible=False,
        formal_exclusion_reasons=["COMBINATION_REJECT"],
    )
    diagnostics = recommendation_diagnostics([candidate], [], {
        "spurs_single": PredictionResult("success"),
        "spurs_multi": PredictionResult("success"),
        "foldx": PredictionResult("failed", errors=[{
            "code": "NO_CONFIRMABLE_CANDIDATES", "message": "No SPURS pass candidates",
        }]),
    }, mutation_options=[])
    assert diagnostics["outcome"] == "ALL_CANDIDATES_FAILED_FORMAL_GATES"
    assert diagnostics["worker_failures"] == {}
    assert "foldx" in diagnostics["workers_not_run"]


def test_job_completion_retains_empty_panel_reason(tmp_path, monkeypatch):
    from ecd_epitope_scanner import jobs
    from ecd_epitope_scanner.pipeline import PipelineRequest

    store = jobs.JobStore(tmp_path / "jobs")
    job = store.create(PipelineRequest(input_type="uniprot", identifier="P16410"))
    monkeypatch.setattr(jobs, "run_pipeline", lambda request, progress: {
        "design_summary": {
            "recommendation_status": "BLOCKED_EXTERNAL_EVIDENCE",
            "recommendation_diagnostics": {"outcome": "SINGLE_STABILITY_EVIDENCE_UNAVAILABLE"},
        },
    })
    runner = jobs.JobRunner(store)
    try:
        runner._run(job["id"])
    finally:
        runner.executor.shutdown(wait=True)
    completed = store.get(job["id"])
    assert completed["status"] == "COMPLETE"
    assert "no formal panel" in completed["message"]
    assert "SINGLE_STABILITY_EVIDENCE_UNAVAILABLE" in completed["message"]


def test_v41_diagnostics_and_hypothesis_artifacts_are_delivered_in_zip(tmp_path, monkeypatch):
    import zipfile
    from ecd_epitope_scanner import exporting
    from ecd_epitope_scanner.models import ProteinContext
    from ecd_epitope_scanner.util import sha256_file

    names = {
        "v4_manifest_v41.json", "recommendation_diagnostics_v41.json",
        "epitope_hypotheses_v41.json", "hypothesis_membership_v41.csv",
        "probe_hypothesis_matrix_v41.csv", "panel_identifiability_v41.json",
        "surface_graph_v41.json", "surface_display_regions_v41.json",
    }
    for name in names:
        (tmp_path / name).write_text("{}\n", encoding="utf-8")
    (tmp_path / "structure_for_energy.pdb").write_text("END\n", encoding="ascii")
    monkeypatch.setattr(exporting, "build_offline_html", lambda *args: "<html></html>")
    protein = ProteinContext("P16410", "CTLA4", "human", "AC", 1, 2, 2, "single_pass")
    files = exporting.write_scientific_outputs(tmp_path, {}, protein, [], [], [])
    assert names <= files.keys()
    manifest = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert all(manifest["files"][name] == sha256_file(tmp_path / name) for name in names)
    with zipfile.ZipFile(tmp_path / "ecd_epitope_design.zip") as archive:
        assert names <= set(archive.namelist())
    # The pipeline adds its final file map after the initial export. The bundle
    # must contain that final document with matching manifest checksums.
    final_analysis = {"files": files, "validation": "final"}
    (tmp_path / "analysis.json").write_text(json.dumps(final_analysis), encoding="utf-8")
    exporting.finalize_scientific_bundle(tmp_path)
    manifest = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["files"]["analysis.json"] == sha256_file(tmp_path / "analysis.json")
    with zipfile.ZipFile(tmp_path / "ecd_epitope_design.zip") as archive:
        assert json.loads(archive.read("analysis.json")) == final_analysis
