"""Deterministic geometry adapters; retain every raw footprint for audit."""
from __future__ import annotations

import math
from dataclasses import fields
from types import SimpleNamespace

from .hypotheses import hypotheses_from_surface_v2
from .models import FootprintFamily, InferenceHypothesis, V4SurfaceHypothesis
from .candidate_search import _distances, _connected

DEFAULT_HIERARCHY = {
    "weighted_jaccard": 0.75, "maximum_center_geodesic_a": 6.0,
    "maximum_normal_angle_deg": 60.0, "maximum_nominal_geodesic_a": 30.0,
    "maximum_nominal_ca_a": 24.0,
    "core_weight": 1.0, "nominal_weight": 0.7, "extended_weight": 0.3,
}


def build_hierarchy(surface_result, records, evidence, graph, parameters=None, *, confidence=None, predicted=False):
    if surface_result is None:
        return [], [], []
    params = {**DEFAULT_HIERARCHY, **(parameters or {})}
    raw = list(getattr(surface_result, "raw_footprints", []) or surface_result.footprints)
    raw_hypotheses = hypotheses_from_surface_v2(SimpleNamespace(footprints=raw), records, evidence)
    source = {h.hypothesis_id: h for h in raw_hypotheses}
    records_by_position = {r.uniprot_pos: r for r in records}
    from ..structure_analysis import median_region_pae
    groups = {}
    for footprint in raw:
        if footprint.footprint_id in source:
            groups.setdefault(footprint.center_vertex_id, []).append(footprint)
    families = []
    for center, members in sorted(groups.items()):
        members.sort(key=lambda f: f.footprint_id)
        def scale(radius):
            return next((f for f in members if f.radius is not None and abs(f.radius - radius) < 1e-6), None)
        nominal = scale(12) or next((f for f in members if f.generation_mode == "adaptive_area"), members[0])
        eligible = set(graph)
        nominal_positions = set(nominal.residues) & eligible
        if not nominal_positions:
            continue
        core = set((scale(9) or nominal).residues) & nominal_positions
        extended = (set((scale(15) or nominal).residues) | nominal_positions) & eligible
        membership = {p: params["core_weight"] if p in core else params["nominal_weight"] if p in nominal_positions else params["extended_weight"] for p in extended}
        # Areas remain physically defined by Surface V2; scale weights are
        # engineering membership factors, never response probabilities.
        area = {p: max(f.residue_surface_area.get(p, 0.0) for f in members) for p in extended}
        domains = sorted({records_by_position[p].domain for p in nominal_positions if p in records_by_position})
        states = sorted({records_by_position[p].coordinate_status for p in nominal_positions if p in records_by_position})
        local_pae = median_region_pae(sorted(nominal_positions), records, confidence) if predicted and confidence else None
        families.append(FootprintFamily(
            family_id=f"HF_C{center:06d}", center_vertex_id=center,
            member_footprint_ids=[f.footprint_id for f in members],
            core_positions=sorted(core), nominal_positions=sorted(nominal_positions),
            extended_positions=sorted(extended), membership=membership,
            residue_weights={p: area[p] * membership[p] for p in extended},
            representative_footprint_id=nominal.footprint_id,
            center_position=nominal.center_residue_id,
            component_id=nominal.surface_component_id, domains=domains,
            reliability_states=states, mean_normal=nominal.mean_normal,
            structural_uncertainty={
                "contains_review_positions": any(records_by_position[p].review_reasons for p in nominal_positions if p in records_by_position),
                "pae_band": ("HIGH" if local_pae >= 15 else "LOW") if local_pae is not None else "UNAVAILABLE",
            },
        ))
    allowed = set(graph)
    center_distances = {f.family_id: _distances(graph, f.center_position, allowed) for f in families if f.center_position in graph}
    nominal_distance_cache = {}
    def bounded(cluster):
        nominal = set().union(*(set(f.nominal_positions) for f in cluster))
        if not _connected(nominal, graph, nominal):
            return False
        points = [records_by_position[p].ca for p in nominal if p in records_by_position and records_by_position[p].ca is not None]
        if len(points) != len(nominal):
            return False
        if any(math.dist(a, b) > params["maximum_nominal_ca_a"] for i, a in enumerate(points) for b in points[i+1:]):
            return False
        key = frozenset(nominal)
        if key not in nominal_distance_cache:
            nominal_distance_cache[key] = max((max(_distances(graph, p, nominal).values(), default=0) for p in nominal), default=0)
        return nominal_distance_cache[key] <= params["maximum_nominal_geodesic_a"]
    def compatible(a, b):
        if a.component_id is None or a.component_id != b.component_id:
            return False
        if len(a.domains) != 1 or a.domains != b.domains or a.reliability_states != b.reliability_states:
            return False
        if a.structural_uncertainty != b.structural_uncertainty:
            return False
        if predicted and a.structural_uncertainty.get("pae_band") == "UNAVAILABLE":
            return False
        if not a.reliability_states or "NOT_ASSESSED" in a.reliability_states:
            return False
        if a.structural_uncertainty.get("pae_band") == "HIGH":
            return False
        if center_distances.get(a.family_id, {}).get(b.center_position, math.inf) > params["maximum_center_geodesic_a"]:
            return False
        if not a.mean_normal or not b.mean_normal:
            return False
        norm = math.sqrt(sum(v*v for v in a.mean_normal) * sum(v*v for v in b.mean_normal))
        if not norm or sum(x*y for x, y in zip(a.mean_normal, b.mean_normal))/norm < math.cos(math.radians(params["maximum_normal_angle_deg"])):
            return False
        positions = set(a.residue_weights) | set(b.residue_weights)
        denominator = sum(max(a.residue_weights.get(p, 0), b.residue_weights.get(p, 0)) for p in positions)
        numerator = sum(min(a.residue_weights.get(p, 0), b.residue_weights.get(p, 0)) for p in positions)
        return denominator > 0 and numerator/denominator >= params["weighted_jaccard"]
    # Complete-link insertion prevents a chain of overlapping pairs from
    # silently joining two mutually incompatible ends.
    clusters = []
    for family in families:
        target = next((cluster for cluster in clusters if all(compatible(family, other) for other in cluster) and bounded(cluster + [family])), None)
        if target is None:
            clusters.append([family])
        else:
            target.append(family)
    hypotheses = []
    for index, cluster in enumerate(clusters, 1):
        representative = cluster[0]
        original = source[representative.representative_footprint_id]
        data = {field.name: getattr(original, field.name) for field in fields(V4SurfaceHypothesis)}
        support = set().union(*(set(f.extended_positions) for f in cluster))
        nominal = set().union(*(set(f.nominal_positions) for f in cluster))
        core = set().union(*(set(f.core_positions) for f in cluster))
        weights = {p: max(f.residue_weights.get(p, 0) for f in cluster) for p in support}
        prior_by_position = {e.residue.uniprot_pos: e.surface_percentile for e in evidence if e.surface_percentile is not None}
        scored = {p: weights[p] for p in support if p in prior_by_position}
        reviewed = [records_by_position[p] for p in support]
        glycan_positions = [r.uniprot_pos for r in reviewed if r.glycan_shielded]
        glycan_evidence = {
            "mutation_prohibition_positions": sorted(r.uniprot_pos for r in reviewed if r.glycosylation_site or r.glycosylation_sequon),
            "direct_glycosylation_positions": sorted(r.uniprot_pos for r in reviewed if r.glycosylation_site),
            "sequon_positions": sorted(r.uniprot_pos for r in reviewed if r.glycosylation_sequon),
            "near_resolved_glycan_positions": sorted(glycan_positions),
        }
        measured = [r for r in reviewed if r.plddt is not None]
        data.update(hypothesis_id=f"HI{index:04d}", residue_weights=weights,
                    core_positions=sorted(core), boundary_positions=sorted(support-core),
                    membership_fraction={p: max(f.membership.get(p, 0) for f in cluster) for p in support},
                    mutable_positions=sorted(p for p in nominal if not records_by_position[p].hard_exclusions),
                    protected_positions=sorted(p for p in nominal if records_by_position[p].hard_exclusions),
                    realized_area_a2=None, geodesic_diameter_a=None,
                    discotope_prior=sum(prior_by_position[p]*w for p,w in scored.items())/sum(scored.values()) if scored else None,
                    glycan_fraction=sum(weights[p] for p in glycan_positions)/sum(weights.values()) if weights else 0.,
                    glycan_evidence=glycan_evidence,
                    mean_plddt=sum(r.plddt*weights[r.uniprot_pos] for r in measured)/sum(weights[r.uniprot_pos] for r in measured) if measured else None,
                    minimum_plddt=min((r.plddt for r in measured), default=None),
                    low_confidence_surface_fraction=sum(weights[r.uniprot_pos] for r in measured if r.plddt <70)/sum(weights.values()) if measured else None,
                    structural_uncertainty=representative.structural_uncertainty.copy(),
                    quality_flags=sorted({flag for f in cluster for fid in f.member_footprint_ids for flag in source[fid].quality_flags}))
        hypotheses.append(InferenceHypothesis(**data,
            family_ids=[f.family_id for f in cluster], nominal_positions=sorted(nominal),
            source_footprint_ids=sorted({fid for f in cluster for fid in f.member_footprint_ids})))
    return raw_hypotheses, families, hypotheses
