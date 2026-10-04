import pytest
from ecd_epitope_scanner.experimental_effects import (
    ProbeBatchEffect,
    WellCount,
    WellRate,
    calculate_antibody_variant_batch_effect,
    calculate_batch_effect,
    calculate_factorial_interaction,
    project_probe_effects,
    summarize_antibody_variant_effects,
    summarize_factorial_interactions,
    visual_fade,
)


def _well(condition, positive, *, batch="B1", qc="PASS"):
    return WellCount(batch, "FabA:Jurkat:gate1", condition, positive, 100, qc)


def _effect(mutant_percent, *, batch="B1", nc=10, reference=70):
    return calculate_batch_effect(
        [
            _well("M", mutant_percent, batch=batch),
            _well("WT", reference, batch=batch),
            _well("NC", nc, batch=batch),
        ],
        mutant="M", reference="WT", mutant_nc="NC",
    )


def test_batch_effect_uses_matching_nc_and_keeps_raw_difference():
    result = _effect(35)
    assert result.difference_pp == pytest.approx(-35)
    assert result.mutant_control_window_pp == pytest.approx(25)
    assert result.reference_control_window_pp == pytest.approx(60)
    assert result.retention == pytest.approx(25 / 60)
    assert result.decline_pct == pytest.approx(100 * (1 - 25 / 60))
    assert _effect(5, nc=10, reference=10).decline_pct is None


def test_technical_wells_equal_weight_and_overrides_are_visible():
    result = calculate_batch_effect(
        [
            WellCount("B1", "fixed", "M", 1, 10, "PASS_WITH_OVERRIDE"),
            WellCount("B1", "fixed", "M", 90, 100, "PASS"),
            WellCount("B1", "fixed", "WT", 80, 1000, "EXCLUDE"),
            WellCount("B1", "fixed", "WT", 80, 100, "PASS"),
            WellCount("B1", "fixed", "NC", 0, 100, "PASS"),
        ],
        mutant="M", reference="WT", mutant_nc="NC",
    )
    assert result.mutant_fraction == pytest.approx(0.5)
    assert result.reference_fraction == pytest.approx(0.8)
    assert "PASS_WITH_OVERRIDE" in result.qc_flags


def test_percentage_only_wells_use_equal_weight_without_fabricated_counts():
    result = calculate_batch_effect(
        [
            WellRate("B1", "FlowJo table", "M", 0.1),
            WellRate("B1", "FlowJo table", "M", 0.5),
            WellRate("B1", "FlowJo table", "M", 0.7),
            WellRate("B1", "FlowJo table", "WT", 0.6),
            WellRate("B1", "FlowJo table", "WT", 0.8),
            WellRate("B1", "FlowJo table", "WT", 0.9),
            WellRate("B1", "FlowJo table", "NC", 0.05),
            WellRate("B1", "FlowJo table", "NC", 0.15),
            WellRate("B1", "FlowJo table", "NC", 0.10),
        ],
        mutant="M", reference="WT", mutant_nc="NC",
    )

    assert result.mutant_fraction == pytest.approx(0.5)
    assert result.reference_fraction == pytest.approx(0.8)
    assert result.mutant_nc_fraction == pytest.approx(0.1)
    assert result.difference_pp == pytest.approx(-30)
    assert result.decline_pct == pytest.approx(100 * (1 - 0.4 / 0.7))


def test_percentage_only_well_rejects_values_outside_probability_range():
    with pytest.raises(ValueError, match="between 0 and 1"):
        WellRate("B1", "FlowJo table", "M", 1.01)


def test_projection_never_splits_combination_effects():
    rows = [
        ProbeBatchEffect("single", (10,), _effect(35)),
        ProbeBatchEffect("same", (11, 12), _effect(40)),
        ProbeBatchEffect("cross", (12, 20), _effect(20)),
        ProbeBatchEffect("unmapped", (10, 999), _effect(20)),
    ]
    projection = project_probe_effects(rows, {"S1": {10, 11, 12}, "S2": {20}})
    assert projection.single_site_probe_ids == {10: ("single",)}
    assert projection.sector_probe_ids == {"S1": ("same", "single")}
    assert projection.sector_decline_pct["S1"] == pytest.approx(
        (rows[0].effect.decline_pct + rows[1].effect.decline_pct) / 2
    )
    by_id = {probe.probe_id: probe for probe in projection.probes}
    assert by_id["cross"].attribution == "CROSS_SECTOR_UNASSIGNED"
    assert by_id["unmapped"].attribution == "UNMAPPED_OR_NOT_SURFACE"


def test_independent_batches_are_median_summarized_without_event_pooling():
    rows = [
        ProbeBatchEffect("P1", (10,), _effect(58, batch="B1")),
        ProbeBatchEffect("P1", (10,), _effect(100, batch="B2")),
        ProbeBatchEffect("P1", (10,), _effect(70, batch="B3")),
        ProbeBatchEffect("P1", (10,), _effect(5, batch="B4", nc=10, reference=10)),
    ]
    summary = project_probe_effects(rows, {"S1": {10}}, seed=1).probes[0]
    assert summary.decline_pct == pytest.approx(0)
    assert summary.batch_count == 4
    assert summary.valid_decline_batch_count == 3
    assert summary.positive_decline_batch_count == 1
    assert summary.negative_decline_batch_count == 1
    assert summary.zero_decline_batch_count == 1
    assert summary.decline_interval_pct is not None
    assert summary.decline_interval_pct[0] <= 0 <= summary.decline_interval_pct[1]


def test_visual_fade_stays_on_connected_surface():
    graph = {10: {11: 3.0}, 11: {10: 3.0, 12: 10.0}, 12: {11: 10.0}, 20: {}}
    fade = visual_fade(graph, 10)
    assert fade[10] == 1.0
    assert fade[11] == pytest.approx(0.8352702114)
    assert 12 not in fade  # geodesic distance 13 A exceeds the cutoff
    assert 20 not in fade  # disconnected component


def test_factorial_interaction_keeps_four_cells_and_bootstraps_independent_batches():
    parent = [
        calculate_batch_effect(
            [_well("M", 40, batch=f"B{i}"), _well("WT", 70, batch=f"B{i}"),
             _well("NC", 5, batch=f"B{i}")],
            mutant="M", reference="WT", mutant_nc="NC",
        )
        for i in range(1, 4)
    ]
    variant = [
        calculate_batch_effect(
            [_well("M", 60, batch=f"B{i}"), _well("WT", 80, batch=f"B{i}"),
             _well("NC", 5, batch=f"B{i}")],
            mutant="M", reference="WT", mutant_nc="NC",
        )
        for i in range(1, 4)
    ]
    rows = [calculate_factorial_interaction(p, v) for p, v in zip(parent, variant, strict=True)]
    assert rows[0]["parent_antigen_effect_pp"] == pytest.approx(-30)
    assert rows[0]["variant_antigen_effect_pp"] == pytest.approx(-20)
    assert rows[0]["interaction_pp"] == pytest.approx(10)
    assert rows[0]["four_cells"]["parent_antibody"]["mutant_antigen"] == pytest.approx(0.4)
    summary = summarize_factorial_interactions(rows, seed=11, bootstrap_iterations=500)
    assert summary["median_interaction_pp"] == pytest.approx(10)
    assert summary["interaction_interval_pp"] == pytest.approx((10, 10))
    assert summary["positive_batch_count"] == 3
    assert summary["interpretation"] == "DESCRIPTIVE_COUPLING_CLUE"


def test_factorial_interaction_keeps_raw_difference_when_control_window_is_invalid():
    parent = calculate_batch_effect(
        [_well("M", 10), _well("WT", 10), _well("NC", 20)],
        mutant="M", reference="WT", mutant_nc="NC",
    )
    variant = calculate_batch_effect(
        [_well("M", 35), _well("WT", 45), _well("NC", 20)],
        mutant="M", reference="WT", mutant_nc="NC",
    )
    interaction = calculate_factorial_interaction(parent, variant)
    assert interaction["interaction_pp"] == pytest.approx(-10)
    assert "PARENT_REFERENCE_WINDOW_INVALID" in interaction["qc_flags"]
    assert parent.retention is None


def test_antibody_variant_effect_uses_each_antibody_own_nc_and_summarizes_batches():
    rows = []
    for batch, parent_rate, variant_rate in (("B1", 0.75, 0.55), ("B2", 0.70, 0.60), ("B3", 0.80, 0.50)):
        parent = [
            WellCount(batch, "parent", "REFERENCE", int(parent_rate * 100), 100),
            WellCount(batch, "parent", "REFERENCE_NC", 5, 100),
        ]
        variant = [
            WellCount(batch, "variant", "REFERENCE", int(variant_rate * 100), 100),
            WellCount(batch, "variant", "REFERENCE_NC", 10, 100),
        ]
        rows.append(calculate_antibody_variant_batch_effect(parent, variant))
    assert rows[0]["parent_control_window_pp"] == pytest.approx(70)
    assert rows[0]["variant_control_window_pp"] == pytest.approx(45)
    assert rows[0]["variant_minus_parent_pp"] == pytest.approx(-25)
    summary = summarize_antibody_variant_effects(rows, seed=7, bootstrap_iterations=500)
    assert summary["batch_count"] == summary["valid_batch_count"] == 3
    assert summary["median_variant_minus_parent_pp"] == pytest.approx(-25)
    assert summary["interval_pp"] == pytest.approx((-35, -15))
    assert summary["negative_batch_count"] == 3
    assert summary["interpretation"] == "DESCRIPTIVE_ANTIBODY_EFFECT"


def test_antibody_variant_effect_retains_raw_contrast_but_excludes_invalid_windows():
    parent = [
        WellCount("B1", "parent", "REFERENCE", 10, 100),
        WellCount("B1", "parent", "REFERENCE_NC", 20, 100),
    ]
    variant = [
        WellCount("B1", "variant", "REFERENCE", 60, 100),
        WellCount("B1", "variant", "REFERENCE_NC", 10, 100),
    ]
    row = calculate_antibody_variant_batch_effect(parent, variant)
    assert row["variant_minus_parent_pp"] == pytest.approx(60)
    assert row["valid"] is False
    assert "PARENT_REFERENCE_WINDOW_INVALID" in row["qc_flags"]
    summary = summarize_antibody_variant_effects([row])
    assert summary["batch_count"] == 1
    assert summary["valid_batch_count"] == 0
    assert summary["median_variant_minus_parent_pp"] is None
    assert summary["interval_pp"] is None
