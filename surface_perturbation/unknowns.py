"""Distinguish improvable quantile uncertainty from deterministic unknowns."""
import numpy as np


def calibration_unknowns(calibration, required, *, runs=None):
    required = np.asarray(required, bool)
    rows = {}
    for name, value in calibration.items():
        blocks = value.get('null_blocks')
        error = value.get('numerical_error')
        valid = value.get('valid')
        if blocks is None or error is None or valid is None:
            rows[name] = {'reason': 'QUERY_OR_CALIBRATION_MISSING', 'appendable': False}
            continue
        blocks, error, valid = np.asarray(blocks), np.asarray(error), np.asarray(valid, bool)
        finite = np.isfinite(blocks).all(axis=0) & np.isfinite(error) & (error >= 0)
        unresolved = required & ~valid
        deterministic = unresolved & ~finite
        statistical = unresolved & finite
        rows[name] = {'deterministic_unknown_points': int(deterministic.sum()),
            'nonfinite_null_points':int((required&~np.isfinite(blocks).all(axis=0)).sum()),
            'missing_or_invalid_numerical_bound_points':int((required&(~np.isfinite(error)|(error<0))).sum()),
            'quantile_unstable_points': int(statistical.sum()),
            'reason': 'NUMERICAL_OR_FIELD_UNKNOWN' if deterministic.any() else 'QUANTILE_UNSTABLE' if statistical.any() else 'RESOLVED',
            'appendable': bool(statistical.any() and not deterministic.any())}
    observed={}
    if runs is not None:
        observed['topology_rejected_runs']=sum(bool(r.get('topology_changed',False)) for r in runs)
        observed['first_five_topology_rejected_runs']=sum(bool(r.get('topology_changed',False)) for r in runs[:5])
        missing=[];protonation=[];directions=[]
        for run in runs:
            chemistry=run.get('wt_chemistry')
            if chemistry is None:continue
            chemistry=np.asarray(chemistry)
            queried=np.isfinite(chemistry[:,1])
            scalar_unknown=queried&~np.isfinite(chemistry[:,2:4]).all(axis=1)
            vector_unknown=queried&~scalar_unknown&~np.isfinite(chemistry[:,4:10]).all(axis=1)
            missing.append(required&~queried);protonation.append(required&scalar_unknown)
            directions.append(required&vector_unknown)
        for name,masks in (('required_query_missing_points',missing),
            ('wt_histidine_protonation_unknown_points',protonation),('wt_direction_unknown_points',directions)):
            observed[name]=int(np.logical_or.reduce(masks).sum()) if masks else None
        observed['scope']='Observed WT field masks and explicit topology flags; mutant or tool failure causes are not inferred'
    impossible_four=(len(runs)>=5 and observed['first_five_topology_rejected_runs']>=2) if runs is not None else False
    return {'channels': rows, 'observed_failure_modes':observed,
            'appendable': any(r['appendable'] for r in rows.values()) and not impossible_four,
            'extension_stop_reason':'FIRST_FIVE_TOPOLOGY_CANNOT_SUPPORT_FOUR_REPEATS' if impossible_four else None,
            'interpretation': 'Additional repeats do not resolve missing numerical bounds or protonation/direction definitions'}
