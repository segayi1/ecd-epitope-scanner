"""Independent candidate search controls; unmeasured structures stay unmeasured.

All controls use the same substitution menu and graph. These diagnostics do
not infer a physical footprint or waive stability/model/null qualification.
"""
from __future__ import annotations

import time
from ..joint_design_v4.models import V4SurfaceHypothesis
from ..joint_design_v4.candidate_search import generate_candidates, _distances


def signature(candidate):
    return tuple(sorted((o.residue.uniprot_pos,o.mutant_aa) for o in candidate.mutations))


def candidate_hypotheses(route,options,graph,parameters,hypotheses):
    if route=='hard_partition':return None
    if route=='overlapping_windows':return hypotheses or None
    if route!='unpartitioned_local':raise ValueError('unknown physical candidate search route')
    qualified=[o for o in options if o.eligible and not o.review_required]
    radius=float(parameters.get('maximum_geodesic_span_a',24.))/2
    global_allowed=set(graph);windows=[]
    for position in sorted({o.residue.uniprot_pos for o in qualified} & global_allowed):
        neighborhood={p:1. for p,d in _distances(graph,position,global_allowed).items() if d<=radius}
        windows.append(V4SurfaceHypothesis(hypothesis_id=f'GLOBAL_LOCAL_{position}',
            residue_weights=neighborhood,center_position=position,
            generation_mode='unpartitioned_local_control',geodesic_radius_a=radius))
    return windows


def search_controls(sectors,options,graph,parameters,hypotheses,measured=()):
    """Three actual search enumerations, with explicit evidence intersection.

    The unpartitioned control uses overlapping geodesic neighborhoods centered
    at every qualified site. Graph connectivity/local span is preserved; display
    ownership never limits a window. It does not enumerate unlimited long-range
    combinations. Search CPU cost is reported separately from physical workers.
    """
    qualified=[o for o in options if o.eligible and not o.review_required]
    known={signature(c) for c in measured};report={}
    for route in ('hard_partition','unpartitioned_local','overlapping_windows'):
        route_hypotheses=candidate_hypotheses(route,qualified,graph,parameters,hypotheses)
        started=time.perf_counter();diagnostics={}
        candidates=generate_candidates(sectors,qualified,graph,parameters,
                                       hypotheses=route_hypotheses,diagnostics=diagnostics)
        signatures=sorted({signature(c) for c in candidates})
        report[route]={'candidate_signatures':[list(s) for s in signatures],
            'unique_candidates':len(signatures),'measured_candidates':sum(s in known for s in signatures),
            'unmeasured_candidates':sum(s not in known for s in signatures),
            'generation_seconds':time.perf_counter()-started,'search_diagnostics':diagnostics,
            'qualification':'search enumeration only; physical evidence required for all new combinations'}
    return report
