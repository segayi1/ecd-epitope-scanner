"""Surface quadrature and graph sampling audit, separate from hard fields."""
from __future__ import annotations
import numpy as np
from .reference import canonical_reference
from .coverage import scan_support


def sampling_audit(coarse_mesh,fine_mesh,atom_map,records,reliable,designable,*,maximum_anchors=12):
    anchors=sorted(reliable)
    if len(anchors)>maximum_anchors:
        anchors=[anchors[i] for i in np.linspace(0,len(anchors)-1,maximum_anchors,dtype=int)]
    results=[]
    for name,mesh,spacing in (('canonical',coarse_mesh,2.),('dense_independent_centers',coarse_mesh,.5),
                               ('finer_SES',fine_mesh,2.)):
        ref=canonical_reference(mesh,atom_map,records,reliable,designable,spacing=spacing)
        all_area=float(ref['area'].sum());curves={}
        for position in anchors:
            seeds=ref['owner']==position
            footprint=scan_support(ref['graph'],seeds,ref['owner'],radius=12.)
            curves[str(position)]=float(ref['area'][footprint].sum())
        results.append({'setting':name,'spacing_a':spacing,'samples':len(ref['area']),
            'target_area_a2':all_area,'designable_area_a2':float(ref['area'][ref['designable']].sum()),
            'mesh_area_a2':float(mesh.area.sum()),'anchor_geodesic_proxy_area_a2':curves})
    base,dense,fine=results
    denominator=max(base['target_area_a2'],np.finfo(float).eps)
    return {'settings':results,'anchors':anchors,
        'relative_target_area_resampling_difference':abs(fine['target_area_a2']-base['target_area_a2'])/denominator,
        'max_dense_graph_proxy_difference_fraction_all':max((abs(dense['anchor_geodesic_proxy_area_a2'][str(p)]-
            base['anchor_geodesic_proxy_area_a2'][str(p)])/denominator for p in anchors),default=None),
        'status':'MEASURED_SENSITIVITY_NOT_AUTOMATIC_CONVERGENCE',
        'scope':'independent dense centers validate quadrature and geodesic proxy only; no invented dense physical hard fields'}
