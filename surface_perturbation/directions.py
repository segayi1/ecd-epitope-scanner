"""Joint vector noise and repeat direction checks in canonical coordinates."""
from __future__ import annotations
from itertools import combinations
import numpy as np
from .calibration import calibrate_channel

VECTOR_GROUPS={
    'surface_normal_delta':('surface_normal_x','surface_normal_y','surface_normal_z'),
    'donor_delta':('donor_x','donor_y','donor_z'),
    'acceptor_delta':('acceptor_x','acceptor_y','acceptor_z'),
    'aromatic_delta':('aromatic_x','aromatic_y','aromatic_z')}


def consistent_four(vectors,eligible,minimum_cosine=0.):
    if not 0<=minimum_cosine<1:raise ValueError('direction cosine must be in [0,1)')
    values=np.asarray(vectors,float)[:5];eligible=np.asarray(eligible,bool)[:5]
    if len(values)<5:return np.zeros(values.shape[1],bool)
    norms=np.linalg.norm(values,axis=2)
    valid=np.isfinite(values).all(axis=2)&(norms>1e-12)&eligible
    unit=np.divide(values,norms[:,:,None],out=np.zeros_like(values),where=norms[:,:,None]>1e-12)
    cosine=np.einsum('rnc,snc->nrs',unit,unit)
    hit=np.zeros(values.shape[1],bool)
    for group in combinations(range(5),4):
        agreement=np.ones(len(hit),bool)
        for left,right in combinations(group,2):agreement&=cosine[:,left,right]>minimum_cosine
        hit|=valid[list(group)].all(axis=0)&agreement
    return hit


def resolve_vector(values,calibration,names,protocol,*,minimum_cosine=0.):
    values=np.asarray(values,float)
    n=values.shape[1];empty=np.zeros(n,bool);unknown=np.ones(n,bool)
    entries=[calibration.get(name,{}) for name in names]
    if len(values)<5 or any(e.get('protocol_id')!=protocol or 'null_blocks' not in e for e in entries):
        return empty,empty,unknown,{}
    # The same shared-WT block index is retained across components. Norms are
    # Euclidean vector magnitudes, not an arbitrary sum of independent scores.
    null=np.stack([e['null_blocks'] for e in entries],axis=2)
    error=np.linalg.norm(np.stack([e['numerical_error'] for e in entries],axis=1),axis=1)
    null_norm=np.linalg.norm(null,axis=2)
    derived=calibrate_channel(null_norm,error,protocol,'joint_vector_magnitude')
    if 'valid' not in derived:return empty,empty,unknown,{}
    amplitude=np.linalg.norm(values[:5],axis=2)
    known=derived['valid']&np.isfinite(values[:5]).all(axis=(0,2))
    lower_bound = all(e.get('positive_lower_bound') for e in entries)
    if lower_bound:
        holdout = np.stack([e['holdout_blocks'] for e in entries],axis=2)
        failed = consistent_four(holdout, np.linalg.norm(holdout,axis=2)>derived['threshold'], minimum_cosine)
        derived['valid'] &= np.isfinite(holdout).all(axis=(0,2)) & ~failed
        finite = np.isfinite(values[:5]).all(axis=2)
        positive = consistent_four(values[:5], amplitude>derived['threshold'], minimum_cosine)
        possible = np.sum((amplitude>derived['threshold']) | ~finite, axis=0) >= 4
        known = derived['valid'] & (positive | ~possible | finite.all(axis=0))
    hit=consistent_four(values[:5],amplitude>derived['threshold'],minimum_cosine)&known
    q99=consistent_four(values[:5],amplitude>np.maximum(derived['q99'],error),minimum_cosine)&known
    evidence={key:v for key,v in derived.items() if isinstance(v,np.ndarray)}
    evidence['delta_norm']=amplitude
    return hit,q99,~known,evidence
