from __future__ import annotations

import io

import pytest
from ecd_epitope_scanner.experiment_projects import _flowjo_well_id
from ecd_epitope_scanner.flowjo_result_tables import (
    FlowJoTableError,
    parse_flowjo_percentage_table,
)
from openpyxl import Workbook


def test_csv_flowjo_percentages_keep_raw_text_and_do_not_infer_denominators():
    content = (
        b",Live/GFP-A subset/mCherry-A subset | Freq. of Parent\n"
        b"B3.fcs,8.72 %\nC3.fcs,0.064 %\nMean,4.392 %\n"
    )

    result = parse_flowjo_percentage_table(content, "flowjo.csv")

    assert result["gate_label"].endswith("Freq. of Parent")
    assert result["metric"] == "FREQUENCY_OF_PARENT"
    assert result["denominator_available"] is False
    assert [item["well_id"] for item in result["observations"]] == ["B3", "C3"]
    assert result["observations"][0]["fraction"] == pytest.approx(0.0872)
    assert result["observations"][1]["fraction"] == pytest.approx(0.00064)
    assert result["observations"][1]["raw_value"] == "0.064 %"
    assert result["observations"][1]["display_precision_decimal_places"] == 3


def test_xlsx_percent_cell_uses_excel_number_format():
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "FlowJo"
    sheet.append(["Sample", "Live/positive | Freq. of Parent"])
    sheet.append(["B3.fcs", 0.217])
    sheet["B2"].number_format = "0.0%"
    output = io.BytesIO()
    workbook.save(output)

    result = parse_flowjo_percentage_table(output.getvalue(), "flowjo.xlsx")

    assert result["sheet_name"] == "FlowJo"
    assert result["observations"][0]["fraction"] == pytest.approx(0.217)
    assert result["observations"][0]["raw_value"] == "21.7 %"


@pytest.mark.parametrize(
    "content, message",
    [
        (b"sample,rate\nB3,20 %\nB3,21 %\n", "重复"),
        (b"sample,rate\nB3,0.21\n", "尺度"),
        (b"sample,rate\nB3,120 %\n", "超出"),
    ],
)
def test_ambiguous_or_invalid_flowjo_rows_are_rejected(content, message):
    with pytest.raises(FlowJoTableError, match=message):
        parse_flowjo_percentage_table(content, "flowjo.csv")


def test_bare_csv_values_require_an_explicit_scale():
    percent = parse_flowjo_percentage_table(
        b"sample,rate\nB3,21\n", "flowjo.csv", numeric_scale="PERCENT"
    )
    fraction = parse_flowjo_percentage_table(
        b"sample,rate\nB3,0.21\n", "flowjo.csv", numeric_scale="FRACTION"
    )

    assert percent["observations"][0]["fraction"] == pytest.approx(0.21)
    assert fraction["observations"][0]["fraction"] == pytest.approx(0.21)


def test_comma_decimal_text_requires_explicit_separator_and_is_preserved():
    content = b"sample;rate\nB3;2,1 %\nC3;0,064 %\n"

    with pytest.raises(FlowJoTableError, match="无法判断小数/千位分隔符"):
        parse_flowjo_percentage_table(content, "flowjo.csv")

    result = parse_flowjo_percentage_table(
        content, "flowjo.csv", decimal_separator=","
    )

    assert result["decimal_separator"] == ","
    assert result["observations"][0]["fraction"] == pytest.approx(0.021)
    assert result["observations"][0]["raw_value"] == "2,1 %"
    assert result["observations"][0]["display_precision_decimal_places"] == 1
    assert result["observations"][1]["fraction"] == pytest.approx(0.00064)
    assert result["observations"][1]["raw_value"] == "0,064 %"


def test_mixed_decimal_separators_are_rejected():
    with pytest.raises(FlowJoTableError, match="混用或不明确"):
        parse_flowjo_percentage_table(
            b"sample;rate\nB3;2.1 %\n",
            "flowjo.csv",
            decimal_separator=",",
        )


def test_internal_whitespace_is_rejected_without_changing_valid_percent_spacing():
    with pytest.raises(FlowJoTableError, match="数值内部不能含空白字符"):
        parse_flowjo_percentage_table(b"sample,rate\nB3,2 1%\n", "flowjo.csv")

    result = parse_flowjo_percentage_table(
        b"sample,rate\nB3,2.1 %\n", "flowjo.csv"
    )
    assert result["observations"][0]["fraction"] == pytest.approx(0.021)
    assert result["observations"][0]["raw_value"] == "2.1 %"


@pytest.mark.parametrize(
    "source, expected",
    [
        ("B5.fcs", "B5"),
        ("c010.FCS", "C10"),
        ("20260903-panel_E11.fcs", "E11"),
        ("well-not-registered", None),
    ],
)
def test_flowjo_sample_ids_normalize_to_registered_plate_wells(source, expected):
    assert _flowjo_well_id(source) == expected
