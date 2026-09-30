"""Independent exposure, coordinate and mutation qualification audit."""
from __future__ import annotations

from typing import Any


def apply_coordinate_qualification(records, reliability):
    rows = {row["residue_id"]: row for row in (reliability or {}).get("residues", [])}
    for record in records:
        row = rows.get(record.uniprot_pos)
        record.coordinate_eligible = bool(record.ca is not None and (row or {}).get("use_for_3d_footprint", True))
        record.coordinate_status = (row or {}).get("state", "NOT_ASSESSED")
        if not record.coordinate_eligible:
            record.mutation_qualification = "EXCLUDED_GEOMETRY"
            record.folding_screen_pass = False


def qualification_audit(records, threshold=0.25) -> dict[str, Any]:
    baseline = {r.uniprot_pos for r in records if r.rsasa >= threshold}
    return {
        "policy": "total_rsasa_only; sidechain_sasa_descriptive; coordinate_and_mutation_independent",
        "rsasa_threshold": threshold,
        "residues": [{
            "position": r.uniprot_pos, "aa": r.aa, "rsasa": r.rsasa,
            "sidechain_sasa_a2": r.sidechain_sasa,
            "surface_exposed": r.surface_accessible,
            "coordinate_eligible": r.coordinate_eligible,
            "coordinate_status": r.coordinate_status,
            "mutation_qualification": r.mutation_qualification,
            "hard_exclusions": r.hard_exclusions, "review_reasons": r.review_reasons,
            "risk_flags": r.risk_flags,
        } for r in records],
        "exposure_sensitivity": [{
            "threshold": value,
            "surface_positions": [r.uniprot_pos for r in records if r.rsasa >= value],
            "added_vs_default": [r.uniprot_pos for r in records if r.rsasa >= value and r.uniprot_pos not in baseline],
            "coordinate_eligible_count": sum(r.rsasa >= value and r.coordinate_eligible for r in records),
            "interpretation": "Exposure counts only; candidate/panel comparison requires separate pipeline replay.",
        } for value in (0.15, 0.20, 0.25)],
    }
