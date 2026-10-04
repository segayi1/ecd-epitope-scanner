"""Prioritize further real verification at connected blind surface patches."""
from __future__ import annotations
import numpy as np
from scipy.sparse.csgraph import connected_components


def blindspot_priority(reference,covered,candidates,proxies):
    area=np.asarray(reference['area']);covered=np.asarray(covered,bool)
    ids=np.flatnonzero(~covered&(area>0));largest=np.zeros(len(area),bool)
    if len(ids):
        _,labels=connected_components(reference['graph'][ids][:,ids],directed=False)
        sizes=np.bincount(labels,weights=area[ids]);largest[ids[labels==int(np.argmax(sizes))]]=True
    diagnostics={}
    for candidate in candidates:
        proxy=np.asarray(proxies[candidate.candidate_id],bool)
        if proxy.shape!=area.shape:raise ValueError('blindspot priority reference mismatch')
        diagnostics[candidate.candidate_id]={'largest_blind_proxy_area_a2':float(area[proxy&largest].sum()),
            'new_proxy_area_a2':float(area[proxy&~covered].sum()),
            'interpretation':'priority only; actual combination repeats and matched null still required'}
    ordered=sorted(candidates,key=lambda c:(-diagnostics[c.candidate_id]['largest_blind_proxy_area_a2'],
        -diagnostics[c.candidate_id]['new_proxy_area_a2'],c.structural_risk,c.candidate_id))
    return ordered,diagnostics
