"""Candidate-local area completeness, independent of global any-valid flags."""
import numpy as np


def local_scope(calibration, required_mask, area):
    required=np.asarray(required_mask,bool);area=np.asarray(area,float)
    if required.shape != area.shape:raise ValueError('local calibration scope shape mismatch')
    from .directions import VECTOR_GROUPS,resolve_vector
    vector_names={name for names in VECTOR_GROUPS.values() for name in names}
    known=np.zeros(required.shape,bool);channels={};components={}
    for name,value in calibration.items():
        valid=np.asarray(value.get('valid',np.zeros(required.shape,bool)),bool)
        if valid.shape!=required.shape:raise ValueError('local calibration reference mismatch')
        denominator=float(area[required].sum())
        fraction=float(area[required&valid].sum())/denominator if denominator>0 else None
        diagnostic={'calibrated_local_area_fraction':fraction,
            'complete':bool(required.any() and denominator>0 and valid[required].all())}
        if name in vector_names:
            components[name]=diagnostic
        else:
            channels[name]=diagnostic
            known |= valid
    # Directional channels use joint vector null magnitudes. Independent
    # component stability cannot establish convergence of the joint field.
    for label,names in VECTOR_GROUPS.items():
        if not any(name in calibration for name in names):continue
        protocol=calibration.get(names[0],{}).get('protocol_id')
        _,_,unknown,_=resolve_vector(np.zeros((5,len(area),3)),calibration,names,protocol)
        valid=~unknown
        known |= valid
        denominator=float(area[required].sum())
        fraction=float(area[required&valid].sum())/denominator if denominator>0 else None
        channels[label]={'calibrated_local_area_fraction':fraction,
            'complete':bool(required.any() and denominator>0 and valid[required].all())}
    denominator=float(area[required].sum())
    fraction=float(area[required&known].sum())/denominator if denominator>0 else None
    return {'complete':bool(channels) and all(v['complete'] for v in channels.values()),
        'calibrated_local_area_fraction':fraction,
        'unresolved_local_area_fraction':1.-fraction if fraction is not None else None,
        'channels':channels,
        'component_diagnostics':components,
        'scope_contract':'scalar_and_joint_vector_local_scope_2',
        'interpretation':'complete requires every declared scalar and joint-vector hard channel at every required direct-support sample; coordinate components are diagnostics, not independent physical channels; union fraction alone cannot grant completeness'}
