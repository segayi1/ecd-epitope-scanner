"""Shared paired-field classification for design and independent dense audits."""
import numpy as np
from .null_store import FIELD_CHANNELS
from .chemistry import CHANNELS
from .directions import VECTOR_GROUPS,resolve_vector
from .calibration import resolved_changes
from .coverage import scan_support


def classify_fields(runs,pairs,native,records,reference,calibration,protocol,*,minimum_cosine=0.,support=None,run_support=None):
    from .engine import direct_support,direct_site_masks,supported_run_values
    if support is None:support=direct_support(pairs[0],native,records,reference)
    if run_support is None:run_support=np.array([direct_support(p,native,records,reference) for p in pairs[:len(runs)]])
    resolved=np.zeros(len(support),bool);known=np.zeros(len(support),bool)
    channel_hits={};channel_known={};resolved_q99=np.zeros(len(support),bool)
    vector_names={name for names in VECTOR_GROUPS.values() for name in names}
    for channel in FIELD_CHANNELS:
        if channel in vector_names:continue
        values=supported_run_values(np.array([r[channel] for r in runs]),run_support)
        hit,missed=resolved_changes(values,calibration.get(channel,{}),protocol)
        missed|=~run_support[:5].any(axis=0)
        resolved|=hit;known|=~missed;channel_hits[channel]=hit;channel_known[channel]=~missed
        value=calibration.get(channel,{})
        if 'q99' in value:
            stricter={**value,'threshold':np.maximum(value['q99'],value['numerical_error'])}
            q99_hit,_=resolved_changes(values,stricter,protocol);resolved_q99|=q99_hit
    vector_evidence={}
    for label,names in VECTOR_GROUPS.items():
        values=np.stack([[r[name] for name in names] for r in runs],axis=0).transpose(0,2,1)
        values=supported_run_values(values,run_support)
        hit,q99_hit,missed,evidence=resolve_vector(values,calibration,names,protocol,minimum_cosine=minimum_cosine)
        missed|=~run_support[:5].any(axis=0)
        resolved|=hit;resolved_q99|=q99_hit;known|=~missed
        for name in names:channel_hits[name]=hit;channel_known[name]=~missed
        vector_evidence.update({label+'__'+key:value for key,value in evidence.items()})
    repeat_support=np.sum(run_support[:5],axis=0)>=4
    changed=resolved&repeat_support
    repeated_sites=[direct_site_masks(p,native,records,reference) for p in pairs[:min(5,len(runs))]]
    if len({len(row) for row in repeated_sites})!=1:raise ValueError('direct site support is ambiguous across repeated models')
    atom_names=('atom_occupancy_gain_a','atom_occupancy_loss_a')
    definitions=(('atom_occupancy_shape',atom_names),('shape_chemistry',(*atom_names,*CHANNELS)),('shape_chemistry_ses',FIELD_CHANNELS))
    group_counts={label:np.zeros(len(changed),int) for label,_ in definitions}
    site_counts=np.zeros(len(changed),int);q99_counts=np.zeros(len(changed),int);site_changes=[]
    for site_index in range(len(repeated_sites[0]) if len(repeated_sites)>=5 else 0):
        site_runs=np.array([row[site_index] for row in repeated_sites])
        local={};local_q99=np.zeros(len(changed),bool)
        for channel in FIELD_CHANNELS:
            if channel in vector_names:continue
            values=supported_run_values(np.array([r[channel] for r in runs[:5]]),site_runs)
            value=calibration.get(channel,{})
            local[channel],_=resolved_changes(values,value,protocol)
            if 'q99' in value:
                stricter={**value,'threshold':np.maximum(value['q99'],value['numerical_error'])}
                q99_hit,_=resolved_changes(values,stricter,protocol);local_q99|=q99_hit
        for names in VECTOR_GROUPS.values():
            values=np.stack([[r[name] for name in names] for r in runs[:5]],axis=0).transpose(0,2,1)
            values=supported_run_values(values,site_runs)
            hit,q99_hit,_,_=resolve_vector(values,calibration,names,protocol,minimum_cosine=minimum_cosine)
            for name in names:local[name]=hit
            local_q99|=q99_hit
        # One residue must support the SAME four modeled changes. Intersecting
        # two independently four-of-five aggregate masks can share only three.
        seed=np.any(list(local.values()),axis=0)&changed
        site_changes.append(seed)
        site_counts+=scan_support(reference['graph'],seed,np.zeros(len(changed),int))
        q99_seed=local_q99&resolved_q99&repeat_support
        q99_counts+=scan_support(reference['graph'],q99_seed,np.zeros(len(changed),int))
        for label,names in definitions:
            baseline=np.any([channel_hits[name] for name in names],axis=0)&repeat_support
            seeds=np.any([local[name] for name in names],axis=0)&baseline
            group_counts[label]+=scan_support(reference['graph'],seeds,np.zeros(len(changed),int))
    groups={label:counts>=1 for label,counts in group_counts.items()}
    all_known = np.logical_and.reduce(list(channel_known.values())) if channel_known else np.zeros(len(support),bool)
    lower_bound = bool(calibration) and all(v.get('positive_lower_bound') for v in calibration.values())
    if lower_bound:
        known = resolved | all_known
    return {'support':support,'run_support':run_support,'known':known,'repeat_support':repeat_support,
        'direct':changed,'hard':site_counts>=1,'q99_hard':q99_counts>=1,'two_direct_sites':site_counts>=2,
        'groups':groups,'vector_evidence':vector_evidence,'channel_known':channel_known,
        'site_changes':np.array(site_changes,bool).reshape(-1,len(changed))}
