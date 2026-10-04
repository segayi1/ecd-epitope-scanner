"""Optional live fixture: user-provided CTLA4 activation experiment, never copied to repo."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FLOWKIT_PYTHON = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"
FIXTURE = Path(os.environ.get("FACS_260903_DIR", r"D:\MFL lab\FACS\260903"))
WORKSPACE = FIXTURE / "03-Sep-2026.wsp"
PLATES = (
    "260903PLATE1-l32-49",
    "260903PLATE2-Y159-W202",
    "260903PLATE3-W212-CD80",
)
PLATE_REAGENTS = {
    PLATES[0]: {"L32", "L49"},
    PLATES[1]: {"Y159", "W202"},
    PLATES[2]: {"W212", "CD80"},
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run_worker(request: dict[str, object]) -> dict[str, object]:
    process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(request),
        capture_output=True,
        text=True,
        timeout=180,
        check=True,
    )
    return json.loads(process.stdout)


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file()
    or not WORKSPACE.is_file()
    or any(not (FIXTURE / plate).is_dir() for plate in PLATES),
    reason="local 260903 FCS/FlowKit fixture unavailable",
)
def test_real_activation_wsp_count_mismatches_are_informational_for_consistent_batches():
    fcs_files = sorted(FIXTURE.rglob("*.fcs"))
    assert len(fcs_files) == 162
    before = {path: _sha256(path) for path in [WORKSPACE, *fcs_files]}
    all_replayed_ids: list[tuple[str, str]] = []
    qc_status_counts: Counter[str] = Counter()
    qc_finding_counts: Counter[str] = Counter()

    for plate in PLATES:
        files = sorted((FIXTURE / plate).glob("*.fcs"))
        assert len(files) == 54
        response = _run_worker(
            {
                "schema_version": "1.0",
                "request_id": f"260903-{plate}",
                "operation": "apply_wsp_gates",
                "path": str(WORKSPACE),
                "fcs_paths": [str(path) for path in files],
            }
        )

        assert response["status"] == "COMPLETE", response.get("errors")
        records = response["records"]
        assert records["status"] == "COMPLETE"
        assert len(records["samples"]) == 54
        assert len(records["gate_counts"]) == 54 * 3
        verification = records["gate_count_verification"]
        assert verification["status"] == "MISMATCH"
        assert verification["compared_gate_count"] == 54 * 3
        assert verification["mismatched_gate_count"] > 0
        assert verification["unmatched_gate_count"] == 0
        assert len(records["gate_count_diagnostics"]) == 54 * 3
        assert records["unresolved_nodes"] == []
        assert records["gate_method_consistency"]["status"] == "CONSISTENT"
        assert records["gate_method_consistency"]["samples_checked"] == 54
        assert records["gate_method_consistency"]["method_variants"] == 1
        assert records["gate_method_consistency"]["gate_coordinates_compared"] is False
        for sample in records["samples"]:
            all_replayed_ids.append((plate, sample["sample_id"]))
            assert {gate["name"] for gate in sample["gates"]} == {
                "Live",
                "GFP-A, SSC-A subset",
                "mCherry-A, SSC-A subset",
            }
            assert sample["qc"]["events_examined"] == sample["event_count"]
            assert sample["qc"]["events_removed"] == 0
            qc_status_counts[sample["qc"]["status"]] += 1
            qc_finding_counts.update(item["code"] for item in sample["qc"]["findings"])

        assert {
            row["sample_id"] for row in records["gate_count_diagnostics"]
        } == {sample["sample_id"] for sample in records["samples"]}
        assert {
            row["sample_id"] for row in records["gate_counts"]
        } == {sample["sample_id"] for sample in records["samples"]}
        assert all(
            row["flowjo_cached_count"] != row["flowkit_recomputed_count"]
            for row in records["gate_count_diagnostics"]
            if not row["matches"]
        )
        if plate == PLATES[0]:
            b3_positive = next(
                row
                for row in records["gate_count_diagnostics"]
                if row["sample_id"] == "B3.fcs"
                and row["gate_name"] == "mCherry-A, SSC-A subset"
            )
            assert b3_positive["flowjo_cached_count"] == 906
            assert b3_positive["flowkit_recomputed_count"] == 755

    assert len(all_replayed_ids) == 162
    assert len(set(all_replayed_ids)) == 162
    assert {well for _, well in all_replayed_ids} == {
        path.name.rsplit("_", 1)[-1]
        for path in (FIXTURE / PLATES[0]).glob("*.fcs")
    }
    assert qc_status_counts == {"REVIEW_REQUIRED": 162}
    assert qc_finding_counts["DETECTOR_UPPER_SATURATION"] > 0
    assert qc_finding_counts["LOW_GATE_EVENT_COUNT"] > 0
    assert {path: _sha256(path) for path in [WORKSPACE, *fcs_files]} == before


@pytest.mark.skipif(
    not (FIXTURE / "03-Sep-2026 FlowJo table-1.xls").is_file(),
    reason="local 260903 mapping workbook unavailable",
)
def test_real_activation_excel_maps_plate_reagents_and_replicates():
    xlrd = pytest.importorskip("xlrd")
    base = FIXTURE
    summary = xlrd.open_workbook(
        str(base / "03-Sep-2026 FlowJo table-1.xls")
    ).sheet_by_name("Sheet1")
    assert PLATE_REAGENTS == {
        "260903PLATE1-l32-49": {"L32", "L49"},
        "260903PLATE2-Y159-W202": {"Y159", "W202"},
        "260903PLATE3-W212-CD80": {"W212", "CD80"},
    }

    panels = (
        ("L32", "03-Sep-2026 FlowJo table-1.xls", 13, "left"),
        ("L49", "03-Sep-2026 FlowJo table-1.xls", 26, "right"),
        ("Y159", "03-Sep-2026 FlowJo table-2.xls", 37, "left"),
        ("W202", "03-Sep-2026 FlowJo table-2.xls", 48, "right"),
        ("W212", "03-Sep-2026 FlowJo table-3.xls", 59, "left"),
        ("CD80", "03-Sep-2026 FlowJo table-3.xls", 72, "right"),
    )
    expected_constructs = ["CTLA4-WT", *(f"CTLA4-E{i}" for i in range(1, 8)), "NC"]
    for reagent, workbook_name, start_row, side in panels:
        raw_sheet = xlrd.open_workbook(str(base / workbook_name)).sheet_by_index(0)
        raw_values = {
            str(raw_sheet.cell_value(row, 0)): float(
                str(raw_sheet.cell_value(row, 1)).replace("%", "").strip()
            ) / 100
            for row in range(1, raw_sheet.nrows)
            if str(raw_sheet.cell_value(row, 0)).endswith(".fcs")
        }
        block_rows = range(start_row - 1, start_row - 1 + 9)
        assert {
            str(summary.cell_value(row, 8)).strip()
            for row in block_rows
            if summary.cell_value(row, 8)
        } == {reagent}
        wells = ("B", "C", "D") if side == "left" else ("E", "F", "G")
        for condition_index, row in enumerate(block_rows):
            construct = str(summary.cell_value(row, 9)).strip()
            assert construct == expected_constructs[condition_index]
            summary_values = [float(summary.cell_value(row, col)) for col in (10, 11, 12)]
            raw_values_for_wells = [
                raw_values[f"{well}{condition_index + 3}.fcs"] for well in wells
            ]
            # Compare each replicate to its own plate export. The workbook now
            # contains the corrected W212 first replicate as well.
            assert summary_values == pytest.approx(raw_values_for_wells, abs=0.00051)


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file()
    or not WORKSPACE.is_file()
    or any(not (FIXTURE / plate).is_dir() for plate in PLATES),
    reason="local 260903 FCS/FlowKit fixture unavailable",
)
def test_real_activation_wsp_rejects_multiple_plates_with_duplicate_flowjo_sample_ids():
    files = sorted(FIXTURE.rglob("*.fcs"))
    response = _run_worker(
        {
            "schema_version": "1.0",
            "request_id": "260903-duplicate-flowjo-ids",
            "operation": "apply_wsp_gates",
            "path": str(WORKSPACE),
            "fcs_paths": [str(path) for path in files],
        }
    )

    assert response["status"] == "FAILED"
    error = response["errors"][0]["message"]
    assert "duplicate FlowJo sample IDs" in error
    assert "Import one sample group at a time" in error
    assert "54 duplicated IDs" in error


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file()
    or not WORKSPACE.is_file()
    or any(not (FIXTURE / plate).is_dir() for plate in PLATES),
    reason="local 260903 FCS/FlowKit fixture unavailable",
)
def test_real_activation_blocks_a_mixed_gate_channel_transform_within_one_batch(tmp_path: Path):
    workspace_tree = ET.parse(WORKSPACE)
    changed = 0

    def local_name(value: str) -> str:
        return value.rsplit("}", 1)[-1]

    for sample_el in workspace_tree.getroot().iter():
        if local_name(sample_el.tag) != "Sample":
            continue
        keywords = {
            child.attrib.get("name", "").lstrip("$").casefold(): child.attrib.get("value", "")
            for keyword_list in list(sample_el)
            if local_name(keyword_list.tag) == "Keywords"
            for child in list(keyword_list)
            if local_name(child.tag) == "Keyword"
        }
        if str(keywords.get("fil", "")).casefold() != "b3.fcs":
            continue
        sample_node = next(
            (child for child in list(sample_el) if local_name(child.tag) == "SampleNode"), None
        )
        if sample_node is None:
            continue
        transformations = next(
            (child for child in list(sample_el) if local_name(child.tag) == "Transformations"),
            None,
        )
        if transformations is None:
            continue
        for transform in list(transformations):
            parameters = [
                child for child in list(transform) if local_name(child.tag) == "parameter"
            ]
            if not any(
                any(local_name(key) == "name" and value == "FSC-A"
                    for key, value in parameter.attrib.items())
                for parameter in parameters
            ):
                continue
            for attribute, value in list(transform.attrib.items()):
                if local_name(attribute) == "maxRange":
                    transform.set(attribute, "1048576")
                    changed += 1
                    break
            break
    assert changed > 0
    inconsistent_workspace = tmp_path / "inconsistent-method.wsp"
    workspace_tree.write(inconsistent_workspace, encoding="utf-8", xml_declaration=True)

    plate = FIXTURE / PLATES[0]
    files = [next(plate.glob("*_B3.fcs")), next(plate.glob("*_B4.fcs"))]
    response = _run_worker({
        "schema_version": "1.0",
        "request_id": "260903-inconsistent-method",
        "operation": "apply_wsp_gates",
        "path": str(inconsistent_workspace),
        "fcs_paths": [str(path) for path in files],
    })
    records = response["records"]
    assert response["status"] == "WAITING_FOR_USER"
    assert records["gate_counts"] == []
    assert records["gate_method_consistency"]["status"] == "INCONSISTENT"
    assert records["gate_method_consistency"]["method_variants"] == 2
    assert "FLOWJO_GATE_METHOD_INCONSISTENT" in {
        node["reason_code"] for node in records["unresolved_nodes"]
    }
