import csv
import json

from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence, ResidueKey, V4Candidate, V4MutationOption, V4Sector,
)
from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations
from ecd_epitope_scanner.joint_design_v4.spatial_support import annotate_spatial_support


def candidate(identifier, positions, passed=True, ddg=0.0):
    mutations = [V4MutationOption(ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D"),
        "A", f"D{position}A", "charge", 0.0, 1.0, True) for position in positions]
    return V4Candidate(identifier, "VS999", mutations, order=len(positions), formal_eligible=passed,
        stability=CombinationStabilityEvidence(tuple(), ddg, ddg, len(positions)),
        formal_exclusion_reasons=[] if passed else ["CONSENSUS_REVIEW"])


def sectors():
    return [V4Sector("VS001", [1, 2, 3], 1, 30, 30, 0), V4Sector("VS002", [4, 5, 6], 1, 30, 30, 0)]


def graph():
    return {p: {q: 6.0 for q in (p - 1, p + 1) if 1 <= q <= 6} for p in range(1, 7)}


def build(pool, regions=None, parameters=None, epitope_evidence=()):
    areas = dict.fromkeys(range(1, 7), 10)
    annotate_spatial_support(pool, graph(), areas)
    return build_sector_recommendations(pool, regions or sectors(), areas, parameters or {}, graph=graph(), epitope_evidence=epitope_evidence)


def test_exclusive_primary_and_cross_alternative_have_independent_evidence():
    pool = [candidate("exclusive", [1, 2], passed=False, ddg=.7), candidate("cross", [2, 5])]
    rows = build(pool)
    assert rows[0]["candidate_id"] == "exclusive"
    assert rows[0]["status"] == "NEAR_THRESHOLD"
    alt = rows[0]["alternatives"][0]
    assert alt["candidate_id"] == "cross" and alt["status"] == "PASSED"
    assert alt["ownership"] == "CROSS_REGION"
    assert all(value is not None for value in alt["impact_evidence"]["axes"].values())
    assert "policy" not in alt["impact_evidence"]
    assert len(alt["impact_evidence"]["mutation_site_details"]) == 2
    assert set(alt["regional_geometry_support"]["scenarios"]) == {
        "geodesic_9", "geodesic_12", "geodesic_15"}
    assert len(alt["impact_evidence"]["local_interaction_diagnostics"]["pairs"]) == 1
    assert rows[1]["status"] == "NO_COMBINATION"
    assert rows[1]["alternatives"][0]["candidate_id"] == "cross"
    assert rows[0]["gate_details"] and not alt["gate_details"]


def test_primary_merges_repeated_principles_and_excludes_unsafe_fallback():
    one = candidate("one", [1, 2])
    unsafe = candidate("unsafe", [4, 5])
    unsafe.mutations[0].hard_exclusions = ["PROTECTED"]
    rows = build([one, unsafe])
    assert rows[0]["principles"] == ["sector_coverage", "lower_structural_risk", "epitope_support"]
    assert not rows[0]["alternatives"]
    assert rows[1]["candidate_id"] is None and not rows[1]["alternatives"]


def test_input_and_sector_order_do_not_change_metrics():
    pool = [candidate("a", [1, 2]), candidate("b", [4, 5]), candidate("shared", [3, 4])]
    first = build(pool)
    second = build(list(reversed(pool)), list(reversed(sectors())))
    assert first == second
    assert first[0]["coverage_gain_area_a2"] == 20
    assert first[1]["coverage_gain_area_a2"] == 20


def test_complete_details_include_passes_unknown_values_and_local_metrics():
    row = build([candidate("a", [1, 2])])[0]
    assert next(check for check in row["gate_checks"] if check["gate"] == "combination_spurs")["passed"]
    assert row["site_details"][0]["sasa_a2"] is None
    assert row["site_details"][0]["auth_seq_id"] == 1
    assert row["pair_distances"][0]["geodesic_distance_a"] == 6
    scenario = row["regional_geometry_support"]["scenarios"]["geodesic_12"]
    assert "leave_one_out_gain_a2" in scenario and "outside_covered_area_a2" in scenario
    assert row["stability_evidence"]["foldx_ddg_kcal_mol"] == 0


def test_lower_risk_prefers_known_model_values_and_epitope_alternative():
    high = candidate("high", [1, 2], ddg=.4)
    low = candidate("low", [2, 3], ddg=-.2)
    high.mutations[0].epitope_supported = True
    high.mutations[0].epitope_prior = .9
    rows = build([high, low])
    items = [rows[0], *rows[0]["alternatives"]]
    assert next(item for item in items if "lower_structural_risk" in item["principles"])["candidate_id"] == "low"
    assert next(item for item in items if "epitope_support" in item["principles"])["candidate_id"] == "high"


def test_preview_keeps_legacy_cross_primary():
    row = build([candidate("cross", [2, 5])], parameters={"optimization_profile": "coverage_scan_v3_preview"})[0]
    assert row["candidate_id"] == "cross"
    assert "alternatives" not in row


def test_passed_triple_discloses_actual_gain_and_minimum_gate():
    triple = candidate("triple", [1, 2, 3])
    triple.risk_components["higher_order_information_gain"] = .08
    row = build([triple], parameters={"panel_optimizer": {"minimum_coverage_gain": .07}})[0]
    gate = next(check for check in row["gate_checks"] if check["gate"] == "higher_order_gain")
    assert gate["observed"] == .08 and gate["threshold"] == .07 and gate["passed"]
    assert gate["comparison"] == "minimum" and gate["unit"] is None
    assert gate["measurement_basis"] == "information_gain"


def test_old_spurs_only_option_cannot_establish_new_formal_qualification():
    one = candidate("one", [1, 2])
    one.stability.foldx_ddg_kcal_mol = None
    one.stability.consensus_status = "SPURS_ONLY_REVIEW_REQUIRED"
    row = build([one], parameters={"allow_spurs_only_final": True})[0]
    gate = next(check for check in row["gate_checks"] if check["gate"] == "combination_foldx")
    assert gate["required"] and gate["passed"] is None and gate["evidence_missing"]
    assert not row["formal_eligible"] and row["status"] == "NEAR_THRESHOLD"
    assert next(check for check in row["gate_checks"] if check["gate"] == "formal_eligibility")["passed"] is False
    assert one.formal_eligible  # Reading cached evidence must not rewrite it.
    strict = build([one])[0]
    assert next(check for check in strict["gate_checks"] if check["gate"] == "combination_foldx")["required"]


def test_local_nmr_header_and_unknown_predicted_structure_are_distinguished(tmp_path):
    from ecd_epitope_scanner.joint_design_v4.sector_recommendations import read_structure_evidence
    cif = tmp_path / "experimental.cif"
    cif.write_text("data_1AH1\n_entry.id 1AH1\n_exptl.method 'Solution NMR'\n", encoding="utf-8")
    evidence = read_structure_evidence(cif)
    assert evidence["experimental"] and evidence["entry_id"] == "1AH1"
    assert evidence["method"] == "Solution NMR" and evidence["resolution_a"] is None
    assert evidence["plddt_applicable"] is False and evidence["local_precision_status"] == "NOT_EVALUATED"
    unknown = tmp_path / "unknown.cif"
    unknown.write_text("data_unknown\n_entry.id unknown\n", encoding="utf-8")
    assert read_structure_evidence(unknown)["plddt_applicable"] is None
    predicted = read_structure_evidence(unknown, predicted=True)
    assert predicted["plddt_applicable"] is True and not predicted["experimental"]
    assert "plddt" not in predicted  # Never derive model confidence from B factors.
    cif.write_text("data_model\n_exptl.method 'AB INITIO MODEL'\n", encoding="utf-8")
    assert not read_structure_evidence(cif)["experimental"]
    pdb = tmp_path / "nmr.pdb"
    pdb.write_text("HEADER    IMMUNE SYSTEM                           01-JAN-00   1AH1\nEXPDTA    SOLUTION NMR\nEND\n", encoding="ascii")
    assert read_structure_evidence(pdb)["experimental"]


def test_site_confidence_and_risk_have_interpretable_evaluated_evidence():
    from types import SimpleNamespace
    one = candidate("one", [1, 2])
    one.mutations[0].warnings = ["SASA_METHOD_SETTINGS_UNVERIFIED", "SASA_CONTEXT_SENSITIVE"]
    record = SimpleNamespace(uniprot_pos=1, chain_id="A", auth_seq_id=1, insertion_code="", aa="D",
        ca=(1, 2, 3), sasa=42, sidechain_sasa=21, rsasa=.5, plddt=None, risk_flags=["glycan-shielded"])
    source = {"structure_evidence": {"entry_id": "1AH1", "method": "Solution NMR", "experimental": True, "plddt_applicable": False}}
    row = build_sector_recommendations([one], sectors(), dict.fromkeys(range(1, 7), 10), {}, records=[record], source_provenance=source)[0]
    site = row["site_details"][0]
    assert site["coordinate_confidence"]["coordinate_available"] and site["coordinate_confidence"]["mapping_matches"]
    assert "pLDDT不适用" in site["coordinate_confidence"]["conclusion"]
    assert all("尚未解释" not in note["explanation"] for note in site["risk_assessment"]["notes"])
    rules = [check for check in row["gate_checks"] if check.get("check_type") == "rule"]
    assert all(isinstance(check["passed"], bool) for check in rules)
    assert "\ufffd" not in json.dumps(row, ensure_ascii=False)


def test_missing_foldx_reason_distinguishes_scheduling_from_failure():
    one = candidate("one", [1, 2], passed=False, ddg=.8)
    one.stability.foldx_ddg_kcal_mol = None
    row = build_sector_recommendations([one], sectors(), dict.fromkeys(range(1, 7), 10), {},
        source_provenance={"resource_audits": {"VS001": {"scheduled_candidate_ids": []}}})[0]
    missing = next(gate for gate in row["gate_details"] if gate["gate"] == "combination_foldx")
    assert row["validation_status"]["foldx"] == "NOT_SCHEDULED"
    assert "未安排" in missing["description"] and "0.800" in missing["description"]
    unknown = build_sector_recommendations([one], sectors(), dict.fromkeys(range(1, 7), 10), {})[0]
    missing = next(gate for gate in unknown["gate_details"] if gate["gate"] == "combination_foldx")
    assert "不能推断模型失败" in missing["description"]


def test_real_secondary_structure_and_engineered_neighborhood_notes_are_explained():
    one = candidate("one", [1, 2])
    one.mutations[0].warnings = ["beta-strand", "engineered-neighborhood", "salt-bridge"]
    risk = build([one])[0]["site_details"][0]["risk_assessment"]
    notes = {note["code"]: note["explanation"] for note in risk["notes"]}
    assert "β链" in notes["beta-strand"]
    assert "参考序列不一致" in notes["engineered-neighborhood"]
    assert "盐桥" in notes["salt-bridge"] and "不等于已证实" in notes["salt-bridge"]


def test_source_and_combination_model_provenance_are_preserved():
    one = candidate("one", [1, 2])
    one.stability.provenance = {"spurs_multi": {"weights_sha256": "saved-hash", "version": "saved-version"}}
    source = {"models": {"foldx": {"version": "actual-version"}}}
    row = build_sector_recommendations([one], sectors(), dict.fromkeys(range(1, 7), 10), {}, source_provenance=source)[0]
    assert row["provenance"]["models"] == source["models"]
    assert row["provenance"]["combination_models"] == one.stability.provenance


def test_alternative_quantitative_tradeoffs_are_exact_and_missing_baseline_unknown():
    from ecd_epitope_scanner.joint_design_v4.models import EpitopeEvidence
    primary = candidate("primary", [1, 3], ddg=.4)
    alternative = candidate("alternative", [1, 2], ddg=.1)
    alternative.mutations[0].epitope_supported = True
    alternative.mutations[0].epitope_prior = .6
    epitope = [EpitopeEvidence(option.residue, surface_percentile=option.epitope_prior, status="SCORED")
               for option in [*primary.mutations, *alternative.mutations]]
    rows = build([primary, alternative, candidate("cross", [3, 4])], epitope_evidence=epitope)
    baseline = rows[0]
    alternate = next(item for item in baseline["alternatives"] if item["candidate_id"] == "alternative")
    comparison = alternate["primary_comparison"]
    assert comparison["baseline_candidate_id"] == "primary"
    assert abs(comparison["worst_threshold_utilization_difference"] - (.1 / .5 - .4 / .5)) < 1e-12
    assert comparison["epitope_supported_fraction_difference"] == .5
    assert comparison["mean_epitope_prior_difference"] == .3
    assert comparison["coverage_difference"] == alternate["regional_geometry_support"]["mean_coverage_fraction"] - baseline["regional_geometry_support"]["mean_coverage_fraction"]
    absent = rows[1]["alternatives"][0]["primary_comparison"]
    assert absent["baseline_candidate_id"] is None
    assert absent["coverage_difference"] is None and absent["worst_threshold_utilization_difference"] is None


def test_csv_flattens_alternatives_without_increasing_primary_selection(tmp_path):
    from ecd_epitope_scanner.joint_design_v4.models import JointV4Result, V4Panel
    from ecd_epitope_scanner.joint_design_v4.reporting import write_v4_outputs
    pool = [candidate("a", [1, 2]), candidate("cross", [2, 5]), candidate("review", [4, 5], passed=False, ddg=2.5)]
    rows = build(pool)
    panel = V4Panel(["a"], {}, {}, [], "test", sector_recommendations=rows)
    result = JointV4Result("4.0", "FORMAL_RECOMMENDATION", {}, [], [], [], sectors(), pool, panel)
    write_v4_outputs(tmp_path, result, parameters={}, contexts={})
    data = list(csv.DictReader((tmp_path / "sector_recommendations_v4.csv").open(encoding="utf-8")))
    assert len(data) == len(rows) + sum(len(row["alternatives"]) for row in rows)
    assert any(row["recommendation_role"] == "ALTERNATIVE" for row in data)
    exported_alt = next(row for row in data if row["recommendation_role"] == "ALTERNATIVE")
    assert json.loads(exported_alt["primary_comparison"]) == rows[0]["alternatives"][0]["primary_comparison"]
    exported_gap = next(row for row in data if row["sector_id"] == "VS002" and row["candidate_id"] == "review")
    assert json.loads(exported_gap["unavailable_principles"]) == rows[1]["unavailable_principles"]
    assert any(item["gate_details"] for item in rows[1]["unavailable_principles"])
    assert json.loads((tmp_path / "selected_panel_v4.json").read_text(encoding="utf-8"))["selected_candidate_ids"] == ["a"]
