import json
from pathlib import Path

import pytest
from ecd_epitope_scanner.benchmark_panels import evaluate_manifest


def _analysis(engine: str, positions: list[int], *, formal: bool = True) -> dict:
    return {
        "parameters": {"patch_engine": engine},
        "protein": {"sequence_length": 10, "ecd_start": 1, "ecd_end": 10, "ecd_sequence": "ACDEFGHIKL"},
        "source": {"structure_sha256": "same-structure"},
        "request": {"assay_format": "cell_surface"},
        "joint_v4_ml": {"recommendation_status": "FORMAL_RECOMMENDATION" if formal else "BLOCKED_EXTERNAL_EVIDENCE"},
        "designs": [{
            "design_role": "STAGE1_PROBE", "status": "RECOMMENDED_LOW_PREDICTED_RISK",
            "mutations": [{"position": position, "mutant": "A"} for position in positions],
        }],
    }


def test_benchmark_uses_equal_budget_and_keeps_blocked_runs(tmp_path: Path):
    cases = []
    for antigen, truth, selections, formal in (
        ("A", [5, 6], {"surface_v2": [5], "joint_v3": [1], "joint_v4_ml": [5, 6]}, True),
        ("B", [8], {"surface_v2": [1], "joint_v3": [8], "joint_v4_ml": [8]}, False),
    ):
        runs = {}
        for engine, positions in selections.items():
            path = tmp_path / f"{antigen}_{engine}.json"
            path.write_text(
                json.dumps(_analysis(engine, positions, formal=formal)), encoding="utf-8"
            )
            runs[engine] = path.name
        cases.append({
            "antigen_id": antigen, "antibody_id": f"antibody-{antigen}",
            "split": "frozen_test", "true_contact_positions": truth, "runs": runs,
        })
    manifest = tmp_path / "benchmark.json"
    manifest.write_text(json.dumps({"construct_budget": 1, "cases": cases}), encoding="utf-8")
    result = evaluate_manifest(manifest, bootstrap_draws=20)
    assert result["case_count"] == 2
    assert result["summary"]["joint_v4_ml"]["hit_at_least_one"]["mean"] == 0.5
    assert result["summary"]["joint_v4_ml"]["hit_at_least_two"]["mean"] == 0.5
    blocked = next(
        row for row in result["case_rows"]
        if row["antigen_id"] == "B" and row["engine"] == "joint_v4_ml"
    )
    assert blocked["construct_count"] == 0
    assert blocked["budget_shortfall"] == 1


def test_benchmark_rejects_different_structure_between_engines(tmp_path: Path):
    runs = {}
    for engine in ("surface_v2", "joint_v3", "joint_v4_ml"):
        analysis = _analysis(engine, [5])
        if engine == "joint_v4_ml":
            analysis["source"]["structure_sha256"] = "different"
        path = tmp_path / f"{engine}.json"
        path.write_text(json.dumps(analysis), encoding="utf-8")
        runs[engine] = path.name
    manifest = tmp_path / "benchmark.json"
    manifest.write_text(json.dumps({
        "construct_budget": 1,
        "cases": [{
            "antigen_id": "A", "antibody_id": "antibody-A", "split": "frozen_test",
            "true_contact_positions": [5], "runs": runs,
        }],
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="do not share input"):
        evaluate_manifest(manifest, bootstrap_draws=0)


def test_benchmark_rejects_missing_structure_identity(tmp_path: Path):
    runs = {}
    for engine in ("surface_v2", "joint_v3", "joint_v4_ml"):
        analysis = _analysis(engine, [5])
        analysis["source"].pop("structure_sha256")
        path = tmp_path / f"{engine}.json"
        path.write_text(json.dumps(analysis), encoding="utf-8")
        runs[engine] = path.name
    manifest = tmp_path / "benchmark.json"
    manifest.write_text(json.dumps({
        "construct_budget": 1,
        "cases": [{
            "antigen_id": "A", "antibody_id": "antibody-A", "split": "frozen_test",
            "true_contact_positions": [5], "runs": runs,
        }],
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="lacks sequence"):
        evaluate_manifest(manifest, bootstrap_draws=0)
