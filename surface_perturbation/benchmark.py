"""Cheap geometric controls on the same frozen SES denominator and pool.

Their masks are proxy scanning opportunities, not validated physical changes.
Contact labels, when supplied, are evaluated AFTER design only.
"""
from __future__ import annotations
import time
import numpy as np
from scipy.spatial import cKDTree
from .coverage import scan_support,nested_panel,largest_blind


def controls(candidates,records,reference,maximum=12):
    started=time.perf_counter();points=reference['points'];owner=reference['owner']
    mapped={r.uniprot_pos:r for r in records}; matrices={k:[] for k in ('euclidean_12a','geodesic_12a','v4_decay_proxy','sasa_geometry_proxy')}
    for candidate in candidates:
        positions={o.residue.uniprot_pos for o in candidate.mutations}
        seeds=np.isin(owner,list(positions))
        xyz=np.array([mapped[p].ca for p in positions if mapped[p].ca is not None])
        euclidean=np.zeros(len(points),bool)
        if len(xyz): euclidean=cKDTree(xyz).query(points)[0]<=12.
        matrices['euclidean_12a'].append(euclidean)
        matrices['geodesic_12a'].append(scan_support(reference['graph'],seeds,owner,radius=12.))
        matrices['v4_decay_proxy'].append(scan_support(reference['graph'],seeds,owner,radius=min(15.,-6*np.log(.25))))
        exposed=np.array([mapped.get(int(p)) is not None and mapped[int(p)].rsasa>=.25 for p in owner])
        matrices['sasa_geometry_proxy'].append(euclidean&exposed)
    report={};ids=[c.candidate_id for c in candidates]
    for method,masks in matrices.items():
        matrix=np.array(masks,bool).reshape(len(candidates),len(points))
        order=nested_panel(matrix,reference['area'],[c.structural_risk for c in candidates],reference['graph'],ids,maximum)
        covered=np.zeros(len(points),bool);counts=np.zeros(len(points),int);curve=[]
        for k in range(1,13):
            if k<=len(order): covered|=matrix[order[k-1]];counts+=matrix[order[k-1]]
            area=reference['area'];denom=float(area.sum())
            curve.append({'K':k,'actual_constructs':min(k,len(order)),
                          'proxy_center_coverage_all':float(area[covered].sum()/denom) if denom else None,
                          'largest_blind_patch_a2':largest_blind(reference['graph'],covered,area),
                          'area_weighted_redundancy':float(np.sum(area*np.maximum(counts-1,0))/denom) if denom else None})
        report[method]={'curve':curve,'nested_ids':[ids[i] for i in order],
                        'interpretation':'geometric proxy; not verified hard physical coverage'}
    return {'controls':report,'runtime_seconds':time.perf_counter()-started,
            'contact_label_sanity_check':'UNAVAILABLE_WITHOUT_FROZEN_HIDDEN_CONTACT_LABELS',
            'comparison_scope':'same frozen reference and formally gated measured candidate pool'}
