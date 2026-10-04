"""Read already-gated FlowJo percentage tables without inventing event counts."""

from __future__ import annotations

import csv
import io
import math
import re
from typing import Any, Literal


class FlowJoTableError(ValueError):
    """An exported FlowJo result table cannot be mapped unambiguously."""


_WELL = re.compile(r"^([A-P])\s*0?([1-9]|1[0-9]|2[0-4])(?:\.FCS)?$", re.IGNORECASE)


def parse_flowjo_percentage_table(
    content: bytes,
    filename: str,
    *,
    numeric_scale: Literal["AUTO", "PERCENT", "FRACTION"] = "AUTO",
    decimal_separator: Literal[".", ","] = ".",
) -> dict[str, Any]:
    """Parse a two-column, per-well result export.

    The first column must contain plate well IDs and the second the per-well
    percentage. Excel formatting and literal percent signs are honored. Bare
    numeric values in CSV files require an explicit scale because 0.2 could
    mean either 0.2% or 20%. Text decimal separators are explicit so a comma
    decimal such as ``2,1%`` cannot be silently changed into ``21%``.
    """
    if decimal_separator not in {".", ","}:
        raise FlowJoTableError("文本小数分隔符必须选择点 (.) 或逗号 (,)")
    suffix = filename.lower().rsplit(".", 1)[-1] if "." in filename else ""
    if suffix == "xls":
        rows, sheet_name = _read_xls(content)
    elif suffix == "xlsx":
        rows, sheet_name = _read_xlsx(content)
    elif suffix in {"csv", "tsv", "txt"}:
        rows, sheet_name = _read_delimited(content, "\t" if suffix in {"tsv", "txt"} else None)
    else:
        raise FlowJoTableError("上传已圈门的 .xls、.xlsx、.csv 或 .tsv 结果表")

    if not rows:
        raise FlowJoTableError("结果表为空")
    header_row = rows[0]
    gate_label = _cell_text(header_row[1]) if len(header_row) > 1 else ""
    observations: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=1):
        if not row or not row[0]:
            continue
        well_match = _WELL.fullmatch(_cell_text(row[0]))
        if well_match is None:
            continue
        if len(row) < 2 or row[1] in (None, ""):
            raise FlowJoTableError(f"第 {row_number} 行孔位 {well_match.group(1)}{well_match.group(2)} 缺少百分比")
        well_id = f"{well_match.group(1).upper()}{int(well_match.group(2))}"
        if well_id in seen:
            raise FlowJoTableError(f"孔位 {well_id} 在结果表中重复")
        seen.add(well_id)
        fraction, raw_value, decimal_places = _parse_percentage(
            row[1], numeric_scale, decimal_separator
        )
        observations.append({
            "well_id": well_id,
            "fraction": fraction,
            "raw_value": raw_value,
            "display_precision_decimal_places": decimal_places,
            "source_row": row_number,
        })
    if not observations:
        raise FlowJoTableError("未找到形如 B3、D12 或 B3.fcs 的逐孔结果")
    if gate_label and len(gate_label) > 500:
        raise FlowJoTableError("门路径表头过长")
    return {
        "source_filename": filename,
        "sheet_name": sheet_name,
        "gate_label": gate_label or "FlowJo per-well percentage",
        "metric": "FREQUENCY_OF_PARENT" if "freq" in gate_label.casefold() and "parent" in gate_label.casefold()
        else "USER_SUPPLIED_PERCENTAGE",
        "percentage_unit": "percent",
        "denominator_available": False,
        "numeric_scale": numeric_scale,
        "decimal_separator": decimal_separator,
        "observations": observations,
    }


def _read_xls(content: bytes) -> tuple[list[list[Any]], str]:
    try:
        import xlrd
    except ImportError as exc:  # pragma: no cover - declared runtime dependency
        raise FlowJoTableError("读取 .xls 需要安装 xlrd") from exc
    try:
        book = xlrd.open_workbook(file_contents=content, formatting_info=True)
        sheet = book.sheet_by_index(0)
    except Exception as exc:
        raise FlowJoTableError(f"无法读取 .xls 工作簿：{exc}") from exc
    rows: list[list[Any]] = []
    for row_index in range(sheet.nrows):
        cells = []
        for column_index in range(sheet.ncols):
            cell = sheet.cell(row_index, column_index)
            value: Any = cell.value
            if cell.ctype == xlrd.XL_CELL_EMPTY:
                value = None
            elif cell.ctype == xlrd.XL_CELL_NUMBER:
                format_key = book.xf_list[cell.xf_index].format_key
                format_string = book.format_map[format_key].format_str
                value = _NumericCell(float(value), format_string)
            cells.append(value)
        rows.append(cells)
    return rows, sheet.name


def _read_xlsx(content: bytes) -> tuple[list[list[Any]], str]:
    try:
        from openpyxl import load_workbook
    except ImportError as exc:  # pragma: no cover - declared runtime dependency
        raise FlowJoTableError("读取 .xlsx 需要安装 openpyxl") from exc
    try:
        book = load_workbook(io.BytesIO(content), read_only=True, data_only=True)
        sheet = book.worksheets[0]
        rows = []
        for row in sheet.iter_rows(values_only=False):
            cells = []
            for cell in row:
                value = cell.value
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    value = _NumericCell(float(value), cell.number_format or "")
                cells.append(value)
            rows.append(cells)
        name = sheet.title
        book.close()
        return rows, name
    except Exception as exc:
        raise FlowJoTableError(f"无法读取 .xlsx 工作簿：{exc}") from exc


def _read_delimited(content: bytes, delimiter: str | None) -> tuple[list[list[Any]], str]:
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise FlowJoTableError("CSV/TSV 请使用 UTF-8 编码导出") from exc
    if delimiter is None:
        try:
            delimiter = csv.Sniffer().sniff(text[:8192], delimiters=",;\t").delimiter
        except csv.Error:
            delimiter = ","
    try:
        return [list(row) for row in csv.reader(io.StringIO(text), delimiter=delimiter)], ""
    except csv.Error as exc:
        raise FlowJoTableError(f"无法读取分隔文本：{exc}") from exc


class _NumericCell(float):
    def __new__(cls, value: float, number_format: str):
        instance = float.__new__(cls, value)
        instance.number_format = number_format
        return instance


def _cell_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, _NumericCell):
        return str(float(value))
    return str(value).strip()


def _parse_percentage(
    value: Any,
    numeric_scale: Literal["AUTO", "PERCENT", "FRACTION"],
    decimal_separator: Literal[".", ","],
) -> tuple[float, str, int | None]:
    if isinstance(value, _NumericCell) or isinstance(value, (int, float)) and not isinstance(value, bool):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise FlowJoTableError("百分比必须是有限数值")
        number_format = getattr(value, "number_format", "")
        percent_formatted = "%" in number_format
        if percent_formatted:
            fraction = numeric
            decimal_places = _format_decimal_places(number_format)
            raw_text = f"{100 * numeric:.{decimal_places}f} %"
        elif numeric_scale == "PERCENT":
            fraction = numeric / 100
            decimal_places = _number_decimal_places(numeric)
            raw_text = str(numeric)
        elif numeric_scale == "FRACTION":
            fraction = numeric
            decimal_places = _number_decimal_places(numeric)
            raw_text = str(numeric)
        else:
            raise FlowJoTableError("裸数值无法判断是百分数还是比例，请选择输入尺度")
    else:
        raw_text = _cell_text(value)
        text = raw_text
        has_percent = text.endswith("%")
        numeric_text = (text[:-1] if has_percent else text).strip()
        if any(character.isspace() for character in numeric_text):
            raise FlowJoTableError("文本百分比数值内部不能含空白字符")
        if decimal_separator == ",":
            if "." in numeric_text or numeric_text.count(",") > 1:
                raise FlowJoTableError("文本数值的小数分隔符混用或不明确")
            numeric_text = numeric_text.replace(",", ".")
        elif "," in numeric_text:
            raise FlowJoTableError(
                "文本数值包含逗号，无法判断小数/千位分隔符；请明确选择逗号小数格式"
            )
        if has_percent:
            numeric = _to_float(numeric_text)
            fraction = numeric / 100
            decimal_places = _number_decimal_places(numeric_text)
        else:
            numeric = _to_float(numeric_text)
            if numeric_scale == "PERCENT":
                fraction = numeric / 100
            elif numeric_scale == "FRACTION":
                fraction = numeric
            else:
                raise FlowJoTableError("无百分号的文本数值无法判断尺度，请指定 percent 或 fraction")
            decimal_places = _number_decimal_places(numeric_text)
    if not math.isfinite(fraction) or not 0 <= fraction <= 1:
        raise FlowJoTableError(f"百分比超出 0–100% 范围：{raw_text}")
    return fraction, raw_text, decimal_places


def _to_float(value: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise FlowJoTableError(f"无效百分比数值：{value}") from exc


def _number_decimal_places(value: float | str) -> int:
    rendered = str(value)
    if "e" in rendered.casefold():
        rendered = f"{float(value):.12f}".rstrip("0").rstrip(".")
    return len(rendered.partition(".")[2])


def _format_decimal_places(number_format: str) -> int:
    normalized = number_format.split(";")[0]
    if "." not in normalized:
        return 0
    return len(normalized.partition(".")[2].split("%", 1)[0].replace("#", "0").replace(",", ""))
