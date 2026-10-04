"""Strict XLSX transfer format for human QC of imported FlowJo percentage tables."""

from __future__ import annotations

from io import BytesIO
from typing import Any

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill


class QCWorksheetError(ValueError):
    """The uploaded review worksheet is malformed or contains unsupported edits."""


QC_WORKSHEET_HEADERS = (
    "schema_version",
    "project_id",
    "export_project_revision",
    "import_id",
    "batch_id",
    "source_sha256",
    "sample_id",
    "sample_artifact_id",
    "well_id",
    "panel_label",
    "condition_label",
    "antigen_label",
    "replicate_type",
    "replicate_index",
    "displayed_value",
    "percentage",
    "current_qc_review_id",
    "current_status",
    "current_reviewer",
    "current_reason",
    "decision",
    "reviewer",
    "reason",
)

QC_WORKSHEET_EDITABLE_HEADERS = frozenset({"decision", "reviewer", "reason"})
QC_WORKSHEET_STATUSES = frozenset({"PASS", "PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"})
QC_WORKSHEET_MEDIA_TYPE = (
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
)
_MAX_WORKSHEET_BYTES = 2 * 1024 * 1024
_MAX_WORKSHEET_ROWS = 384


def _text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return str(value)


def build_qc_review_worksheet(rows: list[dict[str, Any]]) -> bytes:
    """Build a spreadsheet whose only editable columns are human QC decisions."""
    if not rows or len(rows) > _MAX_WORKSHEET_ROWS:
        raise QCWorksheetError("QC worksheet requires between 1 and 384 observed wells")

    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "QC Review"
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = f"A1:W{len(rows) + 1}"
    sheet.sheet_view.showGridLines = False

    header_fill = PatternFill("solid", fgColor="244A5A")
    editable_fill = PatternFill("solid", fgColor="FFF2CC")
    for column, header in enumerate(QC_WORKSHEET_HEADERS, start=1):
        cell = sheet.cell(row=1, column=column, value=header)
        cell.font = Font(color="FFFFFF", bold=True)
        cell.fill = editable_fill if header in QC_WORKSHEET_EDITABLE_HEADERS else header_fill
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    sheet.row_dimensions[1].height = 30

    for row_number, row in enumerate(rows, start=2):
        unknown = set(row) - set(QC_WORKSHEET_HEADERS)
        missing = set(QC_WORKSHEET_HEADERS) - set(row)
        if unknown or missing:
            raise QCWorksheetError("QC worksheet row does not match the supported schema")
        for column, header in enumerate(QC_WORKSHEET_HEADERS, start=1):
            cell = sheet.cell(row=row_number, column=column, value=_text(row[header]))
            # Store every field as text so Excel cannot reinterpret IDs, values, or reasons.
            cell.number_format = "@"
            cell.alignment = Alignment(vertical="top", wrap_text=header in {
                "displayed_value", "current_reason", "reason",
            })
            if header in QC_WORKSHEET_EDITABLE_HEADERS:
                cell.fill = editable_fill

    widths = {
        "A": 15, "B": 22, "C": 18, "D": 22, "E": 22, "F": 68,
        "G": 34, "H": 34, "I": 12, "J": 20, "K": 24, "L": 24,
        "M": 18, "N": 16, "O": 24, "P": 16, "Q": 34, "R": 22,
        "S": 24, "T": 32, "U": 24, "V": 24, "W": 36,
    }
    for column, width in widths.items():
        sheet.column_dimensions[column].width = width

    instructions = workbook.create_sheet("Instructions")
    instructions.append(["FlowJo percentage QC review worksheet", "Instructions"])
    instructions.append(["Scope", "Review the imported, already-gated percentage observations only."])
    instructions.append(["Editable columns", "decision, reviewer, reason (yellow columns on QC Review)."])
    instructions.append(["Decision values", "PASS; PASS_WITH_OVERRIDE; HOLD; EXCLUDE"])
    instructions.append(["Reason", "Required for PASS_WITH_OVERRIDE, HOLD, and EXCLUDE."])
    instructions.append(["Reviewer", "Required for every submitted decision."])
    instructions.append(["Source columns", "Do not edit any other column or add/remove observation rows."])
    instructions.append(["Import behavior", "The application validates all rows before saving any decision."])
    instructions.append(["Scientific scope", "Percentages do not provide event counts or denominators."])
    instructions.append(["Technical replicates", "Wells are not independent biological replicates."])
    instructions.column_dimensions["A"].width = 24
    instructions.column_dimensions["B"].width = 110
    for cell in instructions[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = header_fill

    output = BytesIO()
    workbook.save(output)
    payload = output.getvalue()
    if len(payload) > _MAX_WORKSHEET_BYTES:
        raise QCWorksheetError("Generated QC worksheet exceeds the supported size")
    return payload


def parse_qc_review_worksheet(payload: bytes) -> list[dict[str, str]]:
    """Read the strict review schema and reject formulas or changed workbook layout."""
    if not payload or len(payload) > _MAX_WORKSHEET_BYTES:
        raise QCWorksheetError("QC worksheet is empty or exceeds the 2 MB limit")
    try:
        workbook = load_workbook(BytesIO(payload), read_only=True, data_only=False)
    except Exception as exc:
        raise QCWorksheetError("Upload a valid .xlsx QC review worksheet") from exc
    try:
        if workbook.sheetnames != ["QC Review", "Instructions"]:
            raise QCWorksheetError("QC worksheet sheets or order were changed")
        sheet = workbook["QC Review"]
        if sheet.max_column != len(QC_WORKSHEET_HEADERS):
            raise QCWorksheetError("QC worksheet columns were added or removed")
        header_cells = next(sheet.iter_rows(min_row=1, max_row=1), ())
        if tuple(cell.value for cell in header_cells) != QC_WORKSHEET_HEADERS:
            raise QCWorksheetError("QC worksheet headers or order were changed")

        rows: list[dict[str, str]] = []
        seen_samples: set[str] = set()
        for cells in sheet.iter_rows(min_row=2, max_col=len(QC_WORKSHEET_HEADERS)):
            if any(cell.data_type == "f" for cell in cells):
                raise QCWorksheetError("Formulas are not allowed in the QC worksheet")
            values = [cell.value for cell in cells]
            if all(value is None or value == "" for value in values):
                continue
            if any(value is not None and not isinstance(value, str) for value in values):
                raise QCWorksheetError("QC worksheet values must remain text cells")
            row = {header: value or "" for header, value in zip(QC_WORKSHEET_HEADERS, values, strict=True)}
            sample_id = row["sample_id"].strip()
            if not sample_id or sample_id in seen_samples:
                raise QCWorksheetError("QC worksheet has a missing or duplicate sample ID")
            seen_samples.add(sample_id)
            rows.append(row)
            if len(rows) > _MAX_WORKSHEET_ROWS:
                raise QCWorksheetError("QC worksheet has too many observation rows")
        if not rows:
            raise QCWorksheetError("QC worksheet contains no observation rows")
        return rows
    finally:
        workbook.close()
