"""Exact fixed-cohort quorum witnesses, including site-specific support."""
from itertools import combinations
import numpy as np

from .null_store import FIELD_CHANNELS
from .directions import VECTOR_GROUPS
from .calibration import calibrate_channel


def witnesses(runs,site_support,calibration):
    result={};vector_names={c for names in VECTOR_GROUPS.values() for c in names}
    for site, support in enumerate(np.asarray(site_support,bool).transpose(1,0,2)):
        for name in FIELD_CHANNELS:
            if name in vector_names:continue
            entry=calibration[name];values=np.array([r[name] for r in runs[:5]])
            valid=np.asarray(entry['valid']);threshold=entry['threshold']
            bits=np.zeros(values.shape[1],np.uint8)
            for group in combinations(range(5),4):
                q=list(group)
                hit=support[q].all(axis=0)&valid&((values[q]>threshold).all(axis=0)|(values[q]<-threshold).all(axis=0))
                bits[(bits==0)&hit]=sum(1<<i for i in group)
            result[f'site{site}__{name}']=bits
        for label,names in VECTOR_GROUPS.items():
            entries=[calibration[c] for c in names]
            values=np.stack([[r[c] for c in names] for r in runs[:5]],axis=0).transpose(0,2,1)
            amplitude=np.linalg.norm(values,axis=2)
            null=np.stack([e['null_blocks'] for e in entries],axis=2)
            error=np.linalg.norm(np.stack([e['numerical_error'] for e in entries],axis=1),axis=1)
            derived=calibrate_channel(np.linalg.norm(null,axis=2),error,entries[0]['protocol_id'],label)
            from .directions import resolve_vector
            supported=np.where(support[...,None],values,0.)
            hit,_,_,_=resolve_vector(supported,calibration,names,entries[0]['protocol_id'])
            unit=np.divide(values,amplitude[...,None],out=np.zeros_like(values),where=amplitude[...,None]>1e-12)
            bits=np.zeros(len(hit),np.uint8)
            for group in combinations(range(5),4):
                q=list(group);selected=hit & support[q].all(axis=0) & (amplitude[q]>derived['threshold']).all(axis=0)
                for a,b in combinations(group,2):selected &= np.sum(unit[a]*unit[b],axis=1)>0
                bits[(bits==0)&selected]=sum(1<<i for i in group)
            result[f'site{site}__{label}']=bits
    return result
