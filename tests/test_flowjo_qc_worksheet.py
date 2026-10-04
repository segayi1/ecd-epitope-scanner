from __future__ import annotations

from io import BytesIO

import pytest
from openpyxl import load_workbook

from ecd_epitope_scanner.flowjo_qc_worksheet import (
    QC_WORKSHEET_EDITABLE_HEADERS,
    QC_WORKSHEET_HEADERS,
    QCWorksheetError,
    build_qc_review_worksheet,
    parse_qc_review_worksheet,
)


def _row(sample_id: str = "sample-1") -> dict[str, str]:
    values = {header: "" for header in QC_WORKSHEET_HEADERS}
    values.update({
        "schema_version": "1.0",
        "project_id": "project-1",
        "export_project_revision": "16",
        "import_id": "import-1",
        "batch_id": "batch-1",
        "source_sha256": "a" * 64,
        "sample_id": sample_id,
        "sample_artifact_id": "artifact-1",
        "well_id": "A1",
        "panel_label": "L32",
        "condition_label": "CTLA4-WT",
        "antigen_label": "CTLA4-WT",
        "replicate_type": "TECHNICAL",
        "replicate_index": "1",
        "displayed_value": "70,0%",
        "percentage": "70",
    })
    return values


def test_qc_review_worksheet_round_trips_text_values_and_decisions():
    source = _row()
    source["decision"] = "PASS_WITH_OVERRIDE"
    source["reviewer"] = "审核人 A"
    source["reason"] = "保留原始显示值 70,0%"
    payload = build_qc_review_worksheet([source])

    workbook = load_workbook(BytesIO(payload), data_only=False)
    sheet = workbook["QC Review"]
    headers = {cell.value: cell for cell in sheet[1]}
    assert tuple(headers) == QC_WORKSHEET_HEADERS
    assert headers["decision"].fill.fgColor.rgb.endswith("FFF2CC")
    for cell in sheet[2]:
        assert cell.data_type == "s"
    assert workbook["Instructions"][2][1].value.startswith("Review the imported")
    workbook.close()

    assert parse_qc_review_worksheet(payload) == [source]
    assert QC_WORKSHEET_EDITABLE_HEADERS == {"decision", "reviewer", "reason"}


@pytest.mark.parametrize("mutation", ["formula", "header", "numeric_cell", "duplicate_sample"])
def test_qc_review_worksheet_rejects_unsafe_or_ambiguous_edits(mutation):
    payload = build_qc_review_worksheet([_row(), _row("sample-2")])
    workbook = load_workbook(BytesIO(payload), data_only=False)
    sheet = workbook["QC Review"]
    if mutation == "formula":
        sheet["U2"] = "=1+1"
    elif mutation == "header":
        sheet["A1"] = "renamed_schema"
    elif mutation == "numeric_cell":
        sheet["P2"] = 70
    else:
        sheet["G3"] = "sample-1"
    edited = BytesIO()
    workbook.save(edited)
    workbook.close()

    with pytest.raises(QCWorksheetError):
        parse_qc_review_worksheet(edited.getvalue())


def test_qc_review_worksheet_rejects_empty_and_oversized_files():
    with pytest.raises(QCWorksheetError):
        parse_qc_review_worksheet(b"")
    with pytest.raises(QCWorksheetError):
        parse_qc_review_worksheet(b"x" * (2 * 1024 * 1024 + 1))
