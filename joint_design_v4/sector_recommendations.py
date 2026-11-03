"""Auditable display-region suggestions without changing scientific eligibility."""
from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

from ..stability import AUTO_EXCLUDED_MUTANTS
from .site_first_display import confirmed_repeat_pass


def _direct_positions(candidate):
    return set(candidate.covered_positions_by_scenario.get(
        "geodesic_0", [option.residue.uniprot_pos for option in candidate.mutations]
    ))


def _unsafe_reasons(candidate):
    reasons = set(candidate.formal_exclusion_reasons) & {
        "AUTO_EXCLUDED_MUTANT", "AUTO_PROTECTED_NATIVE_RESIDUE",
        "DISCONNECTED_CANDIDATE", "SECTOR_GEOMETRY_QC_FAILED",
    }
    if not candidate.connected:
        reasons.add("DISCONNECTED_CANDIDATE")
    if candidate.order not in {2, 3}:
        reasons.add("NOT_PRIMARY_DOUBLE_OR_TRIPLE")
    for option in candidate.mutations:
        if option.mutant_aa in AUTO_EXCLUDED_MUTANTS:
            reasons.add("AUTO_EXCLUDED_MUTANT")
        if option.residue.one_letter in "CGP":
            reasons.add("AUTO_PROTECTED_NATIVE_RESIDUE")
        reasons.update(option.hard_exclusions)
    return sorted(reasons)


def _gate_details(candidate, parameters):
    if candidate.formal_eligible:
        return []
    reasons = set(candidate.formal_exclusion_reasons)
    spurs_params = parameters.get("spurs_multi", {})
    foldx_params = parameters.get("foldx_confirmation", {})
    stability = candidate.stability
    spurs = stability.spurs_multi_ddg_kcal_mol if stability else None
    foldx = stability.foldx_ddg_kcal_mol if stability else None
    contract = (stability.provenance.get("stability_contract") if stability else None)
    contract = contract or spurs_params.get("stability_contract", "legacy")
    spurs_primary = contract == "spurs_primary_foldx_sanity_v3"
    details = {}

    def add(gate, reason, observed=None, threshold=None, *, missing=False, lower=False):
        finite = observed is not None and math.isfinite(float(observed))
        excess = None
        if finite and threshold is not None:
            excess = max(0.0, threshold - observed if lower else observed - threshold)
        details[gate] = {
            "gate": gate, "reason": reason, "observed": observed if finite else None,
            "threshold": threshold, "excess": excess,
            "evidence_missing": missing or (observed is not None and not finite),
        }

    stability_reasons = reasons & {
        "COMBINATION_REJECT", "CONSENSUS_REVIEW", "MODEL_DISAGREEMENT",
        "STRUCTURAL_ENERGY_CONFLICT", "SPURS_ONLY_REVIEW_REQUIRED",
        "FOLDX_ONLY_LEGACY", "NOT_SCORED", "COMBINATION_NOT_SCORED",
        "COMBINATION_SPURS_MISSING", "COMBINATION_FOLDX_MISSING",
    }
    spurs_limit = 1.0 if spurs_primary else float(spurs_params.get("preferred_max_ddg", 0.5))
    foldx_limit = float(foldx_params.get("combination_gate_kcal_mol", 2.0))
    if stability_reasons:
        if not (spurs_primary and candidate.order == 3):
            if spurs is None or not math.isfinite(float(spurs)):
                add("combination_spurs", "COMBINATION_SPURS_MISSING", threshold=spurs_limit, missing=True)
            elif spurs > spurs_limit:
                add("combination_spurs", "COMBINATION_SPURS_ABOVE_THRESHOLD", spurs, spurs_limit)
        if foldx is None or not math.isfinite(float(foldx)):
            add("combination_foldx", "COMBINATION_FOLDX_MISSING", threshold=foldx_limit, missing=True)
        elif foldx > foldx_limit:
            add("combination_foldx", "COMBINATION_FOLDX_ABOVE_THRESHOLD", foldx, foldx_limit)
        if "MODEL_DISAGREEMENT" in reasons:
            disagreement = abs(spurs - foldx) if spurs is not None and foldx is not None else None
            add("model_agreement", "MODEL_DISAGREEMENT", disagreement,
                float(spurs_params.get("disagreement_gate_kcal_mol", 1.5)), missing=disagreement is None)
    if "HIGHER_ORDER_REQUIRES_FOLDX" in reasons:
        add("combination_foldx", "HIGHER_ORDER_REQUIRES_FOLDX", threshold=foldx_limit, missing=True)
    if "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN" in reasons:
        regional = "regional_gain_basis" in candidate.risk_components
        gain = candidate.risk_components.get("regional_higher_order_gain") if regional else candidate.risk_components.get("higher_order_information_gain")
        if gain is None and not regional:
            gain = candidate.risk_components.get("higher_order_coverage_gain")
        add("higher_order_gain", "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN", gain,
            float(parameters.get("panel_optimizer", {}).get("minimum_coverage_gain", 0.05)),
            missing=gain is None, lower=True)
        if regional:
            details["higher_order_gain"]["measurement_basis"] = candidate.risk_components["regional_gain_basis"]
            details["higher_order_gain"]["baseline_candidate_id"] = candidate.risk_components.get("regional_gain_baseline_candidate_id")
    if "HIGHER_ORDER_COVERAGE_BASELINE_MISSING" in reasons:
        add("higher_order_gain", "HIGHER_ORDER_COVERAGE_BASELINE_MISSING", missing=True)
    if "SINGLE_EVIDENCE_MISSING_OR_EXCLUDED" in reasons:
        add("single_evidence", "SINGLE_EVIDENCE_MISSING_OR_EXCLUDED", missing=True)
    mapped_reasons = stability_reasons | {
        "HIGHER_ORDER_REQUIRES_FOLDX", "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN", "HIGHER_ORDER_COVERAGE_BASELINE_MISSING",
        "SINGLE_EVIDENCE_MISSING_OR_EXCLUDED",
    }
    for reason in sorted(reasons - mapped_reasons):
        add(reason.lower(), reason, missing="REVIEW" in reason or "UNVERIFIED" in reason)
    if not details:
        # Preserve an authoritative rejection whose status has no numeric gate.
        for reason in sorted(reasons or {"FORMAL_ELIGIBILITY_NOT_ESTABLISHED"}):
            add(reason.lower(), reason, missing=True)
    return [details[key] for key in sorted(details)]


def _legacy_sector_recommendations(candidates, sectors, surface_areas, parameters, selected=()):
    """Give every region a qualified combination or a clearly labelled review suggestion."""
    candidates = [candidate for candidate in candidates if getattr(candidate, "display_eligible", True)]
    selected_ids = {item if isinstance(item, str) else item.candidate_id for item in selected}
    areas = {int(position): max(0.0, float(area)) for position, area in surface_areas.items()}
    denominator = set(areas)
    covered = set()
    output = []
    for sector in sorted(sectors, key=lambda item: item.sector_id):
        membership = set(sector.residue_positions)
        associated = [candidate for candidate in candidates if membership & {
            option.residue.uniprot_pos for option in candidate.mutations
        }]
        safe = [candidate for candidate in associated if not _unsafe_reasons(candidate)]
        sector_area = sum(areas.get(position, 0.0) for position in membership)
        row = {
            "sector_id": sector.sector_id, "candidate_id": None, "status": "NO_COMBINATION",
            "formal_eligible": False, "failed_gates": [], "gate_details": [],
            "gate_distance": {"unmet_gate_count": 0, "normalized_excess": 0.0, "missing_evidence_count": 0},
            "coverage_gain": 0.0, "coverage_gain_area_a2": 0.0,
            "mutation_positions_in_sector": [], "mutation_notations": [], "selected_in_panel": False,
            "ranking_key": [], "geometry_support": {}, "regional_geometry_support": {},
        }
        if not safe:
            row["failed_gates"] = sorted({reason for candidate in associated for reason in _unsafe_reasons(candidate)}) or ["NO_CANDIDATE_WITH_MUTATION_IN_SECTOR"]
            output.append(row)
            continue
        passed = [candidate for candidate in safe if candidate.formal_eligible]
        ranked = []
        for candidate in passed or safe:
            direct = _direct_positions(candidate) & denominator
            new = direct - covered
            local_gain = sum(areas.get(position, 0.0) for position in new & membership)
            global_gain = sum(areas.get(position, 0.0) for position in new)
            details = _gate_details(candidate, parameters)
            distance = {
                "unmet_gate_count": len(details),
                "normalized_excess": sum((detail["excess"] or 0.0) / max(abs(detail["threshold"] or 0.0), 1e-12) for detail in details),
                "missing_evidence_count": sum(detail["evidence_missing"] for detail in details),
            }
            geometry = getattr(candidate, "geometry_support", {}) or {}
            marginal = geometry.get("mean_marginal_area_a2")
            overlap = geometry.get("mean_overlap_fraction")
            trusted = (geometry.get("area_basis") == "eligible_surface_area_a2"
                       and isinstance(marginal, (int, float)) and math.isfinite(marginal)
                       and isinstance(overlap, (int, float)) and math.isfinite(overlap))
            scenarios = geometry.get("scenarios", {})
            regional = {}
            if geometry.get("area_basis") == "eligible_surface_area_a2" and all(
                isinstance(scenarios.get(f"geodesic_{radius}", {}).get("covered_positions"), list)
                for radius in (9, 12, 15)
            ):
                regional = {"basis": "surface_geodesic_neighborhood_proxy",
                            "area_basis": "eligible_surface_area_a2", "sector_area_a2": sector_area,
                            "scenarios": {}}
                for radius in (9, 12, 15):
                    name = f"geodesic_{radius}"
                    positions = set(scenarios[name]["covered_positions"]) & membership
                    area = sum(areas.get(position, 0.0) for position in positions)
                    regional["scenarios"][name] = {
                        "covered_area_a2": area,
                        "coverage_fraction": area / sector_area if sector_area else None,
                        "uncovered_positions": sorted(membership - positions),
                    }
                regional["mean_coverage_fraction"] = (
                    sum(item["coverage_fraction"] for item in regional["scenarios"].values()) / 3
                    if sector_area else None)
            fraction = regional.get("mean_coverage_fraction")
            key = (() if passed else (
                distance["unmet_gate_count"], distance["missing_evidence_count"], distance["normalized_excess"],
            )) + ((-1, -fraction) if fraction is not None else (0, 0.0)) + (
                (-1, -marginal, overlap) if trusted else (0, 0.0, 0.0)) + (
                -local_gain, -global_gain, candidate.structural_risk, candidate.candidate_id)
            ranked.append((key, candidate, direct, local_gain, details, distance, regional))
        key, candidate, direct, local_gain, details, distance, regional = min(ranked, key=lambda item: item[0])
        row.update({
            "candidate_id": candidate.candidate_id, "status": "PASSED" if candidate.formal_eligible else "NEAR_THRESHOLD",
            "formal_eligible": candidate.formal_eligible,
            "failed_gates": list(candidate.formal_exclusion_reasons), "gate_details": details,
            "gate_distance": distance, "coverage_gain": local_gain / sector_area if sector_area else 0.0,
            "coverage_gain_area_a2": local_gain,
            "mutation_positions_in_sector": sorted(membership & {option.residue.uniprot_pos for option in candidate.mutations}),
            "mutation_notations": sorted(option.notation for option in candidate.mutations),
            "selected_in_panel": candidate.candidate_id in selected_ids, "ranking_key": list(key),
            "geometry_support": getattr(candidate, "geometry_support", {}) or {},
            "regional_geometry_support": regional,
        })
        # Review suggestions cannot establish qualified coverage for later regions.
        if candidate.formal_eligible:
            covered.update(direct)
        output.append(row)
    return output


def _finite(value):
    return value is not None and math.isfinite(float(value))


def read_structure_evidence(path, *, predicted=False):
    """Read deposited method metadata; an uploaded file alone proves no quality."""
    result = {"entry_id": None, "method": None, "resolution_a": None,
              "experimental": False, "structure_kind": "unknown", "metadata_source": None,
              "metadata_status": "UNAVAILABLE", "local_precision_status": "NOT_EVALUATED"}
    if path is not None and Path(path).is_file():
        path = Path(path)
        try:
            if path.suffix.lower() in {".cif", ".mmcif"}:
                from Bio.PDB.MMCIF2Dict import MMCIF2Dict
                metadata = MMCIF2Dict(str(path))
                def first(key):
                    value = metadata.get(key, [None])[0]
                    return None if value in {None, "?", "."} else value
                result.update(entry_id=first("_entry.id"), method=first("_exptl.method"))
                resolution = first("_refine.ls_d_res_high") or first("_em_3d_reconstruction.resolution")
                result["resolution_a"] = float(resolution) if resolution else None
            else:
                from Bio.PDB import parse_pdb_header
                metadata = parse_pdb_header(str(path))
                result.update(entry_id=metadata.get("idcode") or None,
                              method=metadata.get("structure_method") or None,
                              resolution_a=metadata.get("resolution"))
            result.update(metadata_source=path.name, metadata_status="READ_FROM_STRUCTURE_HEADER")
        except (OSError, ValueError, KeyError, TypeError):
            result["metadata_status"] = "METADATA_PARSE_FAILED"
    method = (result["method"] or "").lower()
    experimental_methods = ("nmr", "x-ray diffraction", "electron microscopy", "electron crystallography",
                            "electron diffraction", "neutron diffraction", "fiber diffraction", "powder diffraction", "solution scattering")
    result["experimental"] = bool(any(token in method for token in experimental_methods) and not predicted)
    if not _finite(result["resolution_a"]) or result["resolution_a"] <= 0:
        result["resolution_a"] = None
    result["structure_kind"] = "experimental" if result["experimental"] else "predicted" if predicted else "unknown"
    result["plddt_applicable"] = False if result["experimental"] else True if predicted else None
    return result


def _checks(candidate, parameters):
    """Numeric checks disclose the active contract; formal eligibility stays authoritative."""
    evidence = candidate.stability
    spurs = evidence.spurs_multi_ddg_kcal_mol if evidence else None
    foldx = evidence.foldx_ddg_kcal_mol if evidence else None
    params = parameters.get("spurs_multi", {})
    contract = (evidence.provenance.get("stability_contract") if evidence else None) or params.get("stability_contract", "legacy")
    primary = contract == "spurs_primary_foldx_sanity_v3"
    ordinary = contract == "ordinary_v4_directional_v1"
    checks = []
    def add(gate, value, threshold, required=True, *, lower=False, unit="kcal/mol"):
        known = _finite(value)
        checks.append({"gate": gate, "observed": value if known else None, "check_type": "numeric",
                       "threshold": threshold, "required": required,
                       "passed": (value >= threshold if lower else value <= threshold) if known and required else None,
                       "evidence_missing": not known,
                       "excess": max(0.0, threshold - value if lower else value - threshold) if known else None,
                       "comparison": "minimum" if lower else "maximum",
                       "unit": unit, "contract": contract})
    add("combination_spurs", spurs, 1.0 if primary or ordinary else float(params.get("preferred_max_ddg", .5)),
        not (primary and candidate.order == 3))
    add("combination_foldx", foldx, 2.0 if primary else float(parameters.get("foldx_confirmation", {}).get("combination_gate_kcal_mol", 2.0)))
    add("model_agreement", abs(spurs - foldx) if _finite(spurs) and _finite(foldx) else None,
        float(params.get("disagreement_gate_kcal_mol", 1.5)), not (primary or ordinary))
    if candidate.order == 3:
        regional = "regional_gain_basis" in candidate.risk_components
        gain = candidate.risk_components.get("regional_higher_order_gain") if regional else candidate.risk_components.get("higher_order_information_gain")
        basis = candidate.risk_components["regional_gain_basis"] if regional else "information_gain"
        if gain is None and not regional:
            gain = candidate.risk_components.get("higher_order_coverage_gain")
            basis = "coverage_gain"
        add("higher_order_gain", gain, float(parameters.get("panel_optimizer", {}).get("minimum_coverage_gain", .05)), lower=True, unit=None)
        checks[-1]["measurement_basis"] = basis
        if regional:
            checks[-1]["baseline_candidate_id"] = candidate.risk_components.get("regional_gain_baseline_candidate_id")
            checks[-1]["region_ids"] = candidate.risk_components.get("regional_gain_region_ids", [])
            checks[-1]["baseline_coverage_fraction"] = candidate.risk_components.get("regional_lower_order_coverage")
            checks[-1]["candidate_coverage_fraction"] = candidate.risk_components.get("regional_candidate_coverage")
    for gate, passed, description in [
        ("model_evidence_complete", _finite(spurs) and _finite(foldx), "正式推荐必须同时有SPURS组合与FoldX组合的有限评分。"),
        ("candidate_connectivity", candidate.connected, "组合必须满足连通要求。"),
        ("protected_residue_policy", not any(option.residue.one_letter in "CGP" or option.mutant_aa in AUTO_EXCLUDED_MUTANTS for option in candidate.mutations), "保护残基与替换类型按现有科学规则审查。"),
        ("single_evidence", all(_finite(option.spurs_ddg_kcal_mol) and not option.hard_exclusions for option in candidate.mutations), "每个位点须有单点模型证据且没有硬排除项。"),
        ("structural_confidence", candidate.qualification != "REVIEW_REQUIRED" and not any(option.review_required for option in candidate.mutations), "结构置信度审查与预测折叠损伤分开判断。"),
    ]:
        checks.append({"gate": gate, "passed": passed, "required": True, "check_type": "rule",
            "observed": None, "threshold": None, "evidence_missing": False,
            "unit": None, "contract": contract, "description": description})
    for detail in _gate_details(candidate, parameters):
        if detail["gate"] not in {check["gate"] for check in checks}:
            checks.append({**detail, "required": True, "passed": False, "unit": None, "contract": contract,
                           "check_type": "numeric" if detail.get("threshold") is not None else "rule"})
    checks.append({"gate": "formal_eligibility", "passed": candidate.formal_eligible, "check_type": "rule",
                   "required": True, "observed": None, "threshold": None,
                   "evidence_missing": False, "unit": None, "contract": contract})
    return checks


def _ordinary_adaptive_sector_recommendation(sector, rows, ledger, surface_areas, selected):
    """Project only completed adaptive queue winners onto full evidence rows."""
    from .ordinary_adaptive import PRINCIPLES

    queue = {(str(item.get("sector_id")), item.get("principle")): item
             for item in (ledger or {}).get("queues", [])}
    by_id = {candidate.candidate_id: (candidate, row) for candidate, row in rows}
    winners = {}
    missing = []
    for principle in PRINCIPLES:
        item = queue.get((str(sector.sector_id), principle))
        state = item.get("state", "unknown") if item else "missing_queue"
        completed = [entry for entry in (item or {}).get("completed", [])
                     if entry.get("outcome") == "complete"
                     or (isinstance(entry.get("outcome"), dict)
                         and entry["outcome"].get("outcome") == "complete")]
        winner_id = completed[0].get("candidate_id") if state == "complete" and completed else None
        winner = by_id.get(winner_id) if winner_id is not None else None
        if state == "complete" and not completed:
            state = "winner_not_recorded"
        if winner is not None and winner[1].get("ownership") == "EXCLUSIVE":
            winners[principle] = winner
        else:
            if winner_id is not None:
                state = "winner_row_unavailable"
            missing.append({"principle": principle, "state": state,
                            "candidate_id": winner_id,
                            "completed_candidate_count": len(completed),
                            "pending_count": (item or {}).get("pending_count", 0)})

    # Fixed protocol ordering defines the first display row; labels are merged
    # when the same physical candidate won more than one independent queue.
    primary_pair = next(iter(winners.values()), None)
    if primary_pair:
        primary = primary_pair[1]
        primary["recommendation_role"] = "PRIMARY"
        primary["principles"] = []
        primary["recommendation_reasons"] = []
    else:
        primary = _legacy_sector_recommendations([], [sector], surface_areas, {}, selected)[0]
        primary.update({"recommendation_role": "PRIMARY", "principles": [],
                        "recommendation_reasons": []})
    alternatives = []
    by_candidate = {primary.get("candidate_id"): primary} if primary.get("candidate_id") else {}
    reasons = {
        "epitope_support": "该组合是本区表位支持原则队列完成验证的首个通过项。",
        "mechanism_pareto": "该组合是本区机制帕累托原则队列完成验证的首个通过项。",
        "sector_coverage": "该组合是本区覆盖原则队列完成验证的首个通过项。",
    }
    for principle, (candidate, row) in winners.items():
        target = by_candidate.get(candidate.candidate_id)
        if target is None:
            target = row
            target["recommendation_role"] = "ALTERNATIVE"
            target["principles"] = []
            target["recommendation_reasons"] = []
            alternatives.append(target)
            by_candidate[candidate.candidate_id] = target
        target["principles"].append(principle)
        target["recommendation_reasons"].append(reasons[principle])
    primary["recommendation_policy"] = ledger.get("policy", "ordinary_three_principle_adaptive_v1")
    primary["alternatives"] = alternatives
    primary["formal_candidate_count"] = sum(c.formal_eligible for c, _ in rows if _["ownership"] == "EXCLUSIVE")
    primary["confirmed_alternative_count"] = sum(row.get("status") == "PASSED" for row in alternatives)
    primary["all_formal_local_candidates"] = [deepcopy(row) for c, row in rows
                                               if row["ownership"] == "EXCLUSIVE" and c.formal_eligible]
    primary["unavailable_principles"] = missing
    primary["sector_summary"] = {"residue_count": len(set(sector.residue_positions)),
        "surface_area_a2": sum(surface_areas.get(position, 0) for position in sector.residue_positions),
        "mutable_site_count": len({option.residue.uniprot_pos for c, _ in rows
                                    for option in c.mutations if option.residue.uniprot_pos in sector.residue_positions}),
        "safe_exclusive_candidate_count": sum(row["ownership"] == "EXCLUSIVE" for _, row in rows)}
    return primary


def build_sector_recommendations(candidates, sectors, surface_areas, parameters, selected=(), *,
                                 graph=None, records=None, source_provenance=None, epitope_evidence=(),
                                 confirmation_planning=False):
    """Exclusive sector primary plus independently qualified alternatives."""
    candidates = [candidate for candidate in candidates if getattr(candidate, "display_eligible", True)]
    if parameters.get("optimization_profile") in {"coverage_scan_v2_preview", "coverage_scan_v3_preview"}:
        return _legacy_sector_recommendations(candidates, sectors, surface_areas, parameters, selected)
    from dataclasses import replace
    from .regional_perturbation import POLICY, candidate_impact, rank_impacts, dominates, complete_candidate_impact
    # Enrich display metadata without rewriting scored mutation records.
    candidates = [replace(candidate, mutations=deepcopy(candidate.mutations)) for candidate in candidates]
    complete_candidate_impact(candidates, records, graph)
    if graph is not None and surface_areas is not None:
        from .spatial_support import annotate_spatial_support
        annotate_spatial_support(candidates, graph, surface_areas)
    mechanism_policy = parameters.get("regional_recommendation_policy") == POLICY
    effective = []
    for candidate in candidates:
        missing = [reason for field, reason in (("spurs_multi_ddg_kcal_mol", "COMBINATION_SPURS_MISSING"),
                                                ("foldx_ddg_kcal_mol", "COMBINATION_FOLDX_MISSING"))
                   if not _finite(getattr(candidate.stability, field, None))]
        effective.append(replace(candidate, formal_eligible=False,
            formal_exclusion_reasons=sorted(set(candidate.formal_exclusion_reasons + missing)))
            if candidate.formal_eligible and missing else candidate)
    candidates = effective
    sectors = sorted(sectors, key=lambda sector: sector.sector_id)
    owners = {position: sector.sector_id for sector in sectors for position in sector.residue_positions}
    records = {record.uniprot_pos: record for record in records} if records and not isinstance(records, dict) else (records or {})
    epitope_by_position = {item.residue.uniprot_pos: item for item in epitope_evidence}
    selected_ids = {item if isinstance(item, str) else item.candidate_id for item in selected}
    panel_positions = {f"geodesic_{radius}": set() for radius in (9, 12, 15)}
    for candidate in candidates:
        if candidate.formal_eligible and candidate.candidate_id in selected_ids:
            for name in panel_positions:
                panel_positions[name].update(candidate.geometry_support.get("scenarios", {}).get(name, {}).get("covered_positions", []))
    distances = None
    if graph is not None:
        from .candidate_search import _distances
        distances = {position: _distances(graph, position, set(graph)) for candidate in candidates
                     for position in {option.residue.uniprot_pos for option in candidate.mutations}}
    output = []
    for sector in sectors:
        membership = set(sector.residue_positions)
        rows = []
        safe = [candidate for candidate in candidates if not _unsafe_reasons(candidate)
                and membership & {option.residue.uniprot_pos for option in candidate.mutations}]
        for candidate in safe:
            row = _legacy_sector_recommendations([candidate], [sector], surface_areas, parameters, selected)[0]
            row["impact_evidence"] = candidate.impact_evidence
            row["sensitivity"] = candidate.impact_evidence.get("sensitivity", {})
            positions = {option.residue.uniprot_pos for option in candidate.mutations}
            from .spatial_support import candidate_region_membership
            assignment = candidate_region_membership(candidate, sectors)
            row["ownership"] = assignment["assignment"]
            if mechanism_policy:
                row["recommendation_policy"] = POLICY
                row["impact_evidence"] = candidate.impact_evidence
            row["mutation_sector_ids"] = assignment["region_ids"]
            row["gate_checks"] = _checks(candidate, parameters)
            row["site_details"] = []
            for option in sorted(candidate.mutations, key=lambda option: option.residue.uniprot_pos):
                residue = option.residue
                record = records.get(residue.uniprot_pos)
                epitope = epitope_by_position.get(residue.uniprot_pos)
                epitope_known = epitope is not None and _finite(epitope.surface_percentile)
                epitope_summary = {"basis": "surface_percentile_not_binding_probability",
                    "source_status": epitope.status if epitope else "NOT_PROVIDED",
                    "surface_percentile": epitope.surface_percentile if epitope_known else None,
                    "local_peak": epitope.local_peak if epitope else None, "supported_threshold": .80,
                    "conclusion": "位点达到局部峰值或表面排序前20%的高先验支持条件；这不是结合概率。" if epitope and option.epitope_supported
                    else "表位先验已计算，但未达到局部峰值或表面排序前20%的高支持条件。" if epitope_known
                    else "未提供可核验的表位先验记录；不能把缺失解释为没有表位或零结合概率。"}
                structure = (source_provenance or {}).get("structure_evidence", {})
                coordinate_available = record is not None and getattr(record, "ca", None) is not None and len(record.ca) == 3 and all(_finite(value) for value in record.ca)
                mapping_matches = (record is not None and record.chain_id == residue.chain_id
                    and record.auth_seq_id == residue.auth_seq_id and record.insertion_code == residue.insertion_code
                    and record.aa == residue.one_letter)
                confidence = {**structure, "coordinate_available": coordinate_available,
                    "mapping_matches": mapping_matches, "plddt": getattr(record, "plddt", None),
                    "review_required": option.review_required,
                    "conclusion": "实验结构提供坐标；pLDDT不适用。逐残基实验精度未评估，不能据此宣称局部高置信度。" if structure.get("experimental")
                    else "结构坐标与残基映射分别核对；缺少逐残基置信度时不推断高置信度。"}
                codes = sorted(set(option.hard_exclusions + option.warnings + option.glycan_risk + option.structural_risk_flags
                                   + (getattr(record, "risk_flags", []) if record else [])))
                meanings = {
                    "GLYCAN_SHIELDED": "位点存在糖链遮挡提示，可能降低抗体可及性；不等于突变导致结构不稳定。",
                    "SASA_CONTEXT_SENSITIVE": "位点暴露估计对结构上下文敏感；是否影响准入须看原生归一化上下文的表面分类。",
                    "SASA_METHOD_SETTINGS_UNVERIFIED": "外部暴露度交叉检查的设置尚未完全核验；这是方法诊断，不能据此判定结构损伤。",
                    "SASA_DISAGREEMENT": "两种暴露度估计不一致，应复核计算条件；不直接等于结构不稳定。",
                    "RSASA_BOUNDARY_020_025": "位点相对暴露度处于分类边界附近，表面可及性判断对条件变化敏感。",
                    "glycan-shielded": "位点有糖链遮挡提示，可能限制抗体可及性；这不等于已证实结构损伤。",
                    "native-interface": "位点位于输入结构记录的界面附近，可能影响相关接触；尚无实验功能影响结论。",
                    "domain-interface": "位点邻近结构域界面，需留意局部接触影响；不能仅凭此提示判定突变不稳定。",
                    "beta-strand": "位点位于β链二级结构，应留意局部构象影响；该注释本身不证明突变会破坏结构。",
                    "salt-bridge": "输入结构提示该位点参与盐桥，电荷替换可能改变局部联系；这不等于已证实结构失稳或结合变化。",
                    "engineered-neighborhood": "位点邻近输入结构与参考序列不一致的位置，结构代表性可能受影响；这不是已证实的突变损伤。",
                }
                note_codes = sorted({"GLYCAN_SHIELDED" if code == "glycan-shielded" else code for code in codes})
                risk = {"status": "CAUTION" if codes else "NO_RECORDED_WARNING", "codes": codes,
                        "conclusion": "存在需复核的保护、可及性或方法提示；组合稳定性与完整资格另行审查。" if codes else "未记录额外保护或风险提示；这不代表已证明没有风险。",
                        "notes": [{"code": code, "explanation": meanings.get(code, "尚未解释的保护或风险项，请查看原始证据。")} for code in note_codes]}
                row["site_details"].append({"notation": option.notation, "position": residue.uniprot_pos,
                    "native_aa": residue.one_letter, "mutant_aa": option.mutant_aa,
                    "sector_id": owners.get(residue.uniprot_pos), "chain_id": residue.chain_id,
                    "auth_seq_id": residue.auth_seq_id, "insertion_code": residue.insertion_code,
                    "label_seq_id": residue.label_seq_id, "model_id": residue.model_id,
                    "sasa_a2": getattr(record, "sasa", None), "sidechain_sasa_a2": getattr(record, "sidechain_sasa", None),
                    "sidechain_exposure": option.impact_evidence.get("sidechain_exposure", {}),
                    "chemical_descriptors": option.impact_evidence.get("chemical_descriptors", {}),
                    "rsasa": getattr(record, "rsasa", None), "plddt": getattr(record, "plddt", None),
                    "coordinate_confidence": confidence, "risk_assessment": risk,
                    "accessibility_score": option.accessibility_score,
                    "single_spurs_ddg_kcal_mol": option.spurs_ddg_kcal_mol,
                    "epitope_supported": option.epitope_supported if epitope_evidence and epitope else None,
                    "epitope_prior": option.epitope_prior if epitope_known else None, "epitope_assessment": epitope_summary,
                    "hard_exclusions": option.hard_exclusions, "warnings": option.warnings,
                    "glycan_risk": option.glycan_risk, "structural_risk_flags": option.structural_risk_flags})
            evidence = candidate.stability
            row["selected_in_panel"] = row["selected_in_panel"] and candidate.formal_eligible
            for check in row["gate_checks"]:
                if check["gate"] == "structural_confidence":
                    source = (source_provenance or {}).get("structure_evidence", {})
                    check["evaluated_evidence"] = {"structure_kind": source.get("structure_kind", "unknown"),
                        "method": source.get("method"), "plddt_applicable": source.get("plddt_applicable"),
                        "review_required_sites": [site["position"] for site in row["site_details"] if site["coordinate_confidence"]["review_required"]],
                        "mapped_sites": sum(site["coordinate_confidence"]["mapping_matches"] for site in row["site_details"]),
                        "coordinate_available_sites": sum(site["coordinate_confidence"]["coordinate_available"] for site in row["site_details"]),
                        "mutation_count": candidate.order}
                    check["description"] = ("当前没有触发结构置信度复审；不代表已证明局部坐标精度高。" if check["passed"]
                        else "当前触发结构置信度复审，须与预测折叠损伤分开判断。")
                    if source.get("experimental"):
                        check["description"] += " 实验结构不适用pLDDT，逐残基实验精度尚未评估。"
                elif check["gate"] == "protected_residue_policy":
                    check["evaluated_evidence"] = {"native_residues": [option.residue.one_letter for option in candidate.mutations],
                        "mutant_residues": [option.mutant_aa for option in candidate.mutations],
                        "protected_native_classes": ["C", "G", "P"], "excluded_mutant_classes": sorted(AUTO_EXCLUDED_MUTANTS)}
                elif check["gate"] == "single_evidence":
                    check["evaluated_evidence"] = {"single_scores_present": sum(_finite(option.spurs_ddg_kcal_mol) for option in candidate.mutations),
                                                   "mutation_count": candidate.order,
                                                   "hard_exclusions": [reason for option in candidate.mutations for reason in option.hard_exclusions]}
                elif check["gate"] == "candidate_connectivity":
                    check["evaluated_evidence"] = {"connected": candidate.connected}
            row["stability_evidence"] = ({"spurs_multi_ddg_kcal_mol": evidence.spurs_multi_ddg_kcal_mol,
                "foldx_ddg_kcal_mol": evidence.foldx_ddg_kcal_mol, "disagreement": evidence.disagreement,
                "epistasis_estimate": evidence.epistasis_estimate, "consensus_status": evidence.consensus_status,
                "provenance": evidence.provenance} if evidence else {})
            repeat_diagnostic = candidate.impact_evidence.get("stability_repeat_audit", {})
            row["stability_repeat_diagnostic"] = deepcopy(repeat_diagnostic)
            row["recommendation_risk_label"] = repeat_diagnostic.get("recommendation_risk_label", "REPEAT_EVIDENCE_INCOMPLETE")
            row["recommendation_risk_explanation"] = (
                "FoldX已完成十槽复算且中位数通过；仍存在跨门采样波动，详见逐槽值。"
                if confirmed_repeat_pass(repeat_diagnostic) else
                "FoldX重复值跨越配置门限；中位数资格按现行合同保留，但结果对采样敏感，建议复核。"
                if repeat_diagnostic.get("sampling_decision_sensitive") is True else
                "FoldX重复证据缺失或含无效值，采样敏感性未知。"
                if repeat_diagnostic.get("sampling_decision_diagnostic_status") != "COMPLETE" else
                "本批FoldX重复值未跨越配置门限；请结合门限两侧槽数与中位数查看。"
            )
            row["provenance"] = {"source": source_provenance or {}, "generation_anchor": candidate.sector_id,
                                  "targeted_hypotheses": candidate.targeted_hypotheses,
                                  "residue_numbering": "UniProt positions; author chain/residue mapping retained"}
            row["provenance"]["models"] = (source_provenance or {}).get("models", {})
            row["provenance"]["combination_models"] = evidence.provenance if evidence else {}
            row["provenance"]["structure_evidence"] = (source_provenance or {}).get("structure_evidence", {})
            audits = (source_provenance or {}).get("resource_audits", {})
            scheduling_known = bool(audits) and all("scheduled_candidate_ids" in audit for audit in audits.values())
            scheduled_ids = {identifier for audit in audits.values() for identifier in audit.get("scheduled_candidate_ids", [])}
            foldx_known = evidence is not None and _finite(evidence.foldx_ddg_kcal_mol)
            row["validation_status"] = {"foldx": "SCORED" if foldx_known else
                "SCHEDULED_WITHOUT_FINITE_RESULT" if candidate.candidate_id in scheduled_ids else
                "NOT_SCHEDULED" if scheduling_known else "UNKNOWN"}
            for detail in row["gate_details"]:
                if detail["gate"] == "combination_foldx" and detail["evidence_missing"]:
                    status = row["validation_status"]["foldx"]
                    detail["description"] = ("未安排此组合的FoldX验证；没有组合稳定性确认评分。" if status == "NOT_SCHEDULED"
                        else "已安排FoldX，但没有保存此组合的有效评分；不能据此判定模型运行失败。" if status == "SCHEDULED_WITHOUT_FINITE_RESULT"
                        else "缺少FoldX组合评分；现有记录未说明是否安排验证，不能推断模型失败。")
                    spurs_check = next(check for check in row["gate_checks"] if check["gate"] == "combination_spurs")
                    if status == "NOT_SCHEDULED" and spurs_check["passed"] is False:
                        detail["description"] += f" 当前SPURS预测值{spurs_check['observed']:.3f} kcal/mol，高于门槛{spurs_check['threshold']:.3f} kcal/mol。"
            regional = row["regional_geometry_support"]
            row["direct_area_a2"] = sum(surface_areas.get(position, 0) for position in positions & membership)
            local_12 = regional.get("scenarios", {}).get("geodesic_12", {})
            row["local_12_covered_area_a2"] = local_12.get("covered_area_a2")
            row["local_12_coverage_fraction"] = local_12.get("coverage_fraction")
            row["local_12_uncovered_positions"] = local_12.get("uncovered_positions", [])
            row["local_12_uncovered_area_a2"] = (sum(surface_areas.get(position, 0) for position in membership)
                - local_12["covered_area_a2"] if _finite(local_12.get("covered_area_a2")) else None)
            for name, scenario in regional.get("scenarios", {}).items():
                union = set(candidate.geometry_support["scenarios"][name]["covered_positions"])
                outside = union - membership
                scenario["outside_covered_positions"] = sorted(outside)
                scenario["outside_covered_area_a2"] = sum(surface_areas.get(position, 0) for position in outside)
                scenario["additional_panel_area_a2"] = sum(surface_areas.get(position, 0) for position in union - panel_positions[name])
            if distances is not None:
                from .spatial_support import support_metrics
                local_areas = {position: surface_areas[position] for position in membership if position in surface_areas}
                local = support_metrics(positions, distances, local_areas, membership, sector.sector_id)
                for name, values in local["scenarios"].items():
                    if name in regional.get("scenarios", {}):
                        regional["scenarios"][name].update({key: values[key] for key in (
                            "marginal_area_a2", "leave_one_out_gain_a2", "mean_pair_overlap_fraction", "max_pair_overlap_fraction")})
                regional["mean_marginal_area_a2"] = local["mean_marginal_area_a2"]
                regional["mean_overlap_fraction"] = local["mean_overlap_fraction"]
                from itertools import combinations
                row["pair_distances"] = [{"positions": [left, right], "geodesic_distance_a":
                    distances[left].get(right) if _finite(distances[left].get(right)) else None}
                    for left, right in combinations(sorted(positions), 2)]
            else:
                row["pair_distances"] = []
            distance = row["gate_distance"]
            fraction = regional.get("mean_coverage_fraction")
            marginal = regional.get("mean_marginal_area_a2", candidate.geometry_support.get("mean_marginal_area_a2"))
            overlap = regional.get("mean_overlap_fraction", candidate.geometry_support.get("mean_overlap_fraction"))
            row["ranking_key"] = ([] if candidate.formal_eligible else [distance["unmet_gate_count"],
                distance["missing_evidence_count"], distance["normalized_excess"]]) + (
                [-1, -fraction] if fraction is not None else [0, 0.0]) + (
                [-1, -marginal, overlap] if _finite(marginal) and _finite(overlap) else [0, 0.0, 0.0]) + [
                -row["coverage_gain_area_a2"], candidate.structural_risk, candidate.candidate_id]
            row["principles"] = []
            rows.append((candidate, row))
        adaptive_ledger = (source_provenance or {}).get("ordinary_adaptive_ledger")
        if adaptive_ledger is not None:
            output.append(_ordinary_adaptive_sector_recommendation(
                sector, rows, adaptive_ledger, surface_areas or {}, selected))
            continue
        exclusive = [(candidate, row) for candidate, row in rows if row["ownership"] == "EXCLUSIVE"]
        exclusive.sort(key=lambda item: item[0].candidate_id)
        def coverage_key(item):
            candidate, row = item
            return (not candidate.formal_eligible, *row["ranking_key"])
        coverage_choice = min(exclusive, key=coverage_key) if exclusive else None
        primary = coverage_choice
        if mechanism_policy:
            regional_pool = [item for item in exclusive if item[0].impact_evidence.get("policy") == POLICY]
            passed = [item for item in regional_pool if item[0].formal_eligible]
            pool = passed or regional_pool
            if not passed and pool:
                nearest = min(tuple(item[1]["gate_distance"][key] for key in (
                    "unmet_gate_count", "missing_evidence_count", "normalized_excess")) for item in pool)
                pool = [item for item in pool if tuple(item[1]["gate_distance"][key] for key in (
                    "unmet_gate_count", "missing_evidence_count", "normalized_excess")) == nearest]
            doubles = [item[0] for item in pool if item[0].order == 2]
            triples = [item[0] for item in pool if item[0].order == 3]
            ranked_doubles = rank_impacts(doubles)
            ranked_triples = rank_impacts(triples)
            winner = ranked_doubles[0] if ranked_doubles else None
            # Cross-order replacement requires raw-axis dominance and the existing gain gate.
            if ranked_triples:
                gain_limit = float(parameters.get("panel_optimizer", {}).get("minimum_coverage_gain", .05))
                upgrades = [item for item in ranked_triples if item.risk_components.get("regional_higher_order_gain", -math.inf) is not None
                    and item.risk_components.get("regional_higher_order_gain", -math.inf) >= gain_limit
                    and (winner is None or dominates(candidate_impact(item), candidate_impact(winner)))]
                if upgrades:
                    winner = upgrades[0]
            primary = next((item for item in pool if item[0] is winner), None)
            for candidate, row in regional_pool:
                row["impact_evidence"] = candidate.impact_evidence
                row["sensitivity"] = candidate.impact_evidence.get("sensitivity", {})
        if primary:
            primary_row = primary[1]
            primary_row["principles"] = ["mechanism_pareto" if mechanism_policy else "sector_coverage"]
            primary_row["recommendation_reasons"] = [(
                "本区合格候选中优先双点，按侧链暴露和六项化学变化的帕累托前沿选取；各化学机制组取最强轴效用，组间折中。表位倾向仅作注释，排序不是亲和力预测。" if mechanism_policy else
                "优先选择通过计算门槛的本区专属组合，并按本区三尺度平均估计覆盖排序。") if primary[0].formal_eligible else "本区机制主代表仍需审查；本区其他正式合格组合请展开候选列表查看。"]
            if mechanism_policy and primary[0].formal_eligible and primary[0].impact_evidence.get("missing_axes"):
                primary_row["status"] = "EVIDENCE_LIMITED"
        else:
            primary_row = _legacy_sector_recommendations([], [sector], surface_areas, parameters, selected)[0]
            primary_row.update({"failed_gates": ["NO_SAFE_EXCLUSIVE_COMBINATION"], "ownership": "EXCLUSIVE",
                "principles": ["sector_coverage"], "recommendation_reasons": ["本区没有安全的专属双点或三点组合；跨区组合仅列为补充建议。"]})
        if mechanism_policy:
            primary_row["recommendation_policy"] = POLICY
        primary_row["recommendation_role"] = "PRIMARY"
        alternatives = []
        def choose(pool, metric):
            passed = [item for item in pool if item[0].formal_eligible]
            # Unscored principles still need confirmation priority before final qualification.
            eligible = passed or pool
            return min(eligible, key=lambda item: (() if passed else (
                item[1]["gate_distance"]["unmet_gate_count"], item[1]["gate_distance"]["missing_evidence_count"],
                item[1]["gate_distance"]["normalized_excess"])) + tuple(metric(item)) + tuple(item[1]["ranking_key"])) if eligible else None
        def risk(item):
            checks = [check for check in item[1]["gate_checks"] if check.get("required")
                      and check["gate"] in {"combination_spurs", "combination_foldx"}]
            known = all(check.get("observed") is not None for check in checks)
            return (not known, max((check["observed"] / max(abs(check["threshold"]), 1e-12) for check in checks if check.get("observed") is not None), default=math.inf))
        choices = [
            ("lower_structural_risk", choose(exclusive, risk), "在本区专属组合中优先选择组合模型相对门槛占用较低者。"),
            ("epitope_support", choose(exclusive, lambda item: (-sum(option.epitope_supported for option in item[0].mutations) / item[0].order,
                -sum(option.epitope_prior for option in item[0].mutations) / item[0].order)), "在本区专属组合中优先选择表位先验支持位点比例和平均先验较高者。"),
            ("cross_region_coverage", choose([item for item in rows if item[1]["ownership"] == "CROSS_REGION"], lambda item: (
                -sum(scenario.get("additional_panel_area_a2", 0) for scenario in item[1]["regional_geometry_support"].get("scenarios", {}).values()),)),
                "跨区组合优先补充当前合格主面板尚未覆盖的表面邻域；不增加主面板构建数量。")]
        if mechanism_policy:
            choices = [item for item in choices if item[0] != "epitope_support"]
            choices.insert(0, ("sector_coverage", coverage_choice, "保留原覆盖优先方案作为面积对照；几何邻域不等于实际抗体受影响面积。"))
            for principle, axes in (("charge_extreme", ("charge_delta",)),
                                    ("volume_extreme", ("volume_delta_a3",)),
                                    ("hydropathy_extreme", ("hydropathy_delta", "aromatic_delta")),
                                    ("hbond_extreme", ("donor_delta", "acceptor_delta"))):
                for axis in axes:
                    known = [item for item in exclusive if item[0].impact_evidence.get("policy") == POLICY
                             and _finite(item[0].impact_evidence.get("axes", {}).get(axis))]
                    choices.append((principle, choose(known, lambda item, axis=axis: (
                        -item[0].impact_evidence["axes"][axis],)), f"保留{axis}变化较大的机制备选；不能据此推断结合自由能变化。"))
        missing = []
        by_id = {primary_row["candidate_id"]: primary_row} if primary else {}
        for principle, choice, reason in choices:
            if choice is None or (not choice[0].formal_eligible and not confirmation_planning):
                missing.append({"principle": principle, "reason": "NO_QUALIFIED_PRINCIPLE_CANDIDATE",
                    "candidate_id": choice[0].candidate_id if choice else None,
                    "gate_details": deepcopy(choice[1]["gate_details"]) if choice else []})
                continue
            row = choice[1]
            if row["candidate_id"] in by_id:
                row = by_id[row["candidate_id"]]
            else:
                row["recommendation_role"] = "ALTERNATIVE"
                row["recommendation_reasons"] = []
                alternatives.append(row)
                by_id[row["candidate_id"]] = row
            if principle not in row["principles"]:
                row["principles"].append(principle)
            row["recommendation_reasons"].append(reason)
        primary_row["alternatives"] = alternatives
        def comparison_metrics(candidate, row):
            checks = [check for check in row["gate_checks"] if check.get("required")
                      and check["gate"] in {"combination_spurs", "combination_foldx"}]
            energy = (max(max(0.0, check["observed"]) / max(abs(check["threshold"]), 1e-12) for check in checks)
                      if checks and all(_finite(check.get("observed")) for check in checks) else None)
            priors = [site["epitope_prior"] for site in row["site_details"]]
            supported = [site["epitope_supported"] for site in row["site_details"]]
            return {"mean_coverage_fraction": row["regional_geometry_support"].get("mean_coverage_fraction"),
                    "worst_threshold_utilization": energy,
                    "epitope_supported_fraction": sum(supported) / candidate.order if all(isinstance(value, bool) for value in supported) else None,
                    "mean_epitope_prior": sum(priors) / candidate.order if all(_finite(value) for value in priors) else None}
        baseline = comparison_metrics(*primary) if primary else {}
        candidates_by_id = {candidate.candidate_id: candidate for candidate, _ in rows}
        for alternative in alternatives:
            metrics = comparison_metrics(candidates_by_id[alternative["candidate_id"]], alternative)
            comparison = {"baseline_candidate_id": primary_row["candidate_id"], "baseline_metrics": baseline,
                          "candidate_metrics": metrics, "explanation": "与本区专属首选比较；正覆盖差代表更多估计覆盖，正风险占用差代表更接近或超过门槛。" if primary else "本区没有专属首选，无法计算相对差；保留本组合自己的指标与资格。"}
            for name, field in [("mean_coverage_fraction", "coverage_difference"),
                                ("worst_threshold_utilization", "worst_threshold_utilization_difference"),
                                ("epitope_supported_fraction", "epitope_supported_fraction_difference"),
                                ("mean_epitope_prior", "mean_epitope_prior_difference")]:
                comparison[field] = metrics[name] - baseline[name] if _finite(metrics[name]) and _finite(baseline.get(name)) else None
            alternative["primary_comparison"] = comparison
        for recommendation in [primary_row, *alternatives]:
            if recommendation.get("candidate_id") is None:
                continue
            explanation = recommendation.get("recommendation_risk_explanation")
            if explanation and explanation not in recommendation.setdefault("recommendation_reasons", []):
                recommendation["recommendation_reasons"].append(explanation)
        primary_row["formal_candidate_count"] = sum(c.formal_eligible for c, row in exclusive)
        primary_row["confirmed_alternative_count"] = sum(row.get("status") == "PASSED" for row in alternatives)
        primary_row["all_formal_local_candidates"] = [deepcopy(row) for c, row in exclusive if c.formal_eligible]
        available = {principle for row in [primary_row, *alternatives] if row.get("formal_eligible")
                     for principle in row["principles"]}
        primary_row["unavailable_principles"] = list({row["principle"]: row for row in missing
                                                    if row["principle"] not in available}.values())
        primary_row["sector_summary"] = {"residue_count": len(membership),
            "surface_area_a2": sum(surface_areas.get(position, 0) for position in membership),
            "mutable_site_count": len({option.residue.uniprot_pos for candidate in safe for option in candidate.mutations if option.residue.uniprot_pos in membership}),
            "safe_exclusive_candidate_count": len(exclusive)}
        output.append(primary_row)
    return output
