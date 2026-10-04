"""Target feasibility is distinct from closure of a physical candidate search."""
from __future__ import annotations
import numpy as np

CONTRACT = 'combo_2_3_hard80_v1'


def enabled(settings):
    return settings.get('panel_contract') == CONTRACT


def waiting_confirmation(candidate):
    reasons=set(candidate.formal_exclusion_reasons)
    return (candidate.order in (2,3) and bool(reasons)
            and reasons <= {'SPURS_ONLY_REVIEW_REQUIRED','HIGHER_ORDER_REQUIRES_FOLDX'}
            and candidate.stability is not None
            and candidate.stability.consensus_status=='SPURS_ONLY_REVIEW_REQUIRED')


def actual_pair_controls(triple, resolved_controls):
    """Require all exact real pair counterfactuals, including resolved zero-change controls."""
    from itertools import combinations
    required={frozenset(o.notation for o in pair) for pair in combinations(triple.mutations,2)}
    available={frozenset(o.notation for o in c.mutations):c for c in resolved_controls if c.order==2}
    if len(required)!=3 or not required<=available.keys():return None
    return [available[key] for key in sorted(required,key=lambda key:tuple(sorted(key)))]


def site_evidence(candidate, site_changes, area):
    changes=np.asarray(site_changes,bool);area=np.asarray(area,float)
    ordered=sorted(candidate.mutations,key=lambda o:(o.residue.chain_id,int(o.residue.auth_seq_id),(o.residue.insertion_code or '').strip()))
    positions=[o.residue.uniprot_pos for o in ordered]
    if len(positions)!=candidate.order or len(set(positions))!=len(positions):
        return {'complete':False,'reason':'MUTATION_SITES_NOT_DISTINCT','sites':[]}
    # direct_site_masks sorts native author-chain/residue/insertion identities.
    # Bind the actual masks using the exported mutation order supplied by caller.
    if changes.shape!=(candidate.order,len(area)):
        return {'complete':False,'reason':'DIRECT_SITE_EVIDENCE_MISSING','sites':[]}
    rows=[{'uniprot_pos':position,'notation':option.notation,
           'verified_direct_area_a2':float(area[mask].sum()),'repeat_rule':'same four of five'}
          for position,option,mask in zip(positions,ordered,changes)]
    return {'complete':all(r['verified_direct_area_a2']>0 for r in rows),
            'reason':None if all(r['verified_direct_area_a2']>0 for r in rows) else 'MUTATION_SITE_NO_VERIFIED_CHANGE',
            'sites':rows,'interpretation':'site-supported changes within the actual combination; not isolated causal effects'}


def primary_count(curve, available, initial, maximum=12):
    cap=min(int(available),int(maximum),12)
    start=min(cap,max(4,int(initial)))
    for count in range(start,cap+1):
        value=curve[count-1].get('C_hard_all')
        if value is not None and value>=.8:return count
    return cap


def goal_state(selected, coverage, evidence, *, audit_complete, fatal=False):
    if fatal:return 'INCOMPLETE_EVIDENCE'
    if not 4<=len(selected)<=12:return 'INSUFFICIENT_VERIFIED_CONSTRUCTS'
    if any(c.order not in (2,3) or not evidence.get(c.candidate_id,{}).get('complete') for c in selected):
        return 'INCOMPLETE_SITE_EVIDENCE'
    if coverage is None or coverage<.8:return 'COVERAGE_BELOW_TARGET'
    return 'ACHIEVED' if audit_complete else 'FINAL_AUDIT_PENDING'


def audit_target_panel(folder, summary):
    """Recompute the target from current artifacts instead of trusting a label."""
    import json
    from pathlib import Path
    from .paired import validate_pair
    folder=Path(folder);ids=summary.get('selected_candidate_ids',[])
    if not 4<=len(ids)<=12 or len(ids)!=len(set(ids)):raise ValueError('TARGET_PANEL_COUNT')
    candidates={c['candidate_id']:c for c in json.loads((folder/'panel_candidates.json').read_text(encoding='utf-8'))}
    pairs=json.loads((folder/'surface_perturbation/paired_manifest.json').read_text(encoding='utf-8'))
    with np.load(folder/'surface_perturbation/coverage_matrix.npz',allow_pickle=False) as data:
        names=data['candidate_ids'].tolist();area=data['area'];hard=data['hard']
        if hard.shape!=(len(names),len(area)) or not np.isfinite(area).all() or (area<0).any() or area.sum()<=0:
            raise ValueError('TARGET_PANEL_MATRIX')
        selected=[names.index(cid) for cid in ids]
        actual=float(area[np.any(hard[selected],axis=0)].sum()/area.sum())
    if actual<.8:raise ValueError('TARGET_COVERAGE_BELOW_80')
    for cid in ids:
        c=candidates[cid];positions=[m['residue']['uniprot_pos'] for m in c['mutations']]
        if c['order'] not in (2,3) or len(positions)!=c['order'] or len(set(positions))!=c['order'] or not c['formal_eligible']:
            raise ValueError('TARGET_MUTATION_ORDER_OR_QUALIFICATION')
        if len(pairs[cid])<5:raise ValueError('TARGET_MODEL_REPEATS')
        for pair in pairs[cid][:5]:validate_pair(pair)
        with np.load(folder/f'surface_perturbation/{cid}_direct_site_evidence.npz',allow_pickle=False) as data:
            masks=data['site_changes']
            if masks.shape!=(c['order'],len(area)) or not np.all(np.any(masks,axis=1)):
                raise ValueError('TARGET_SITE_NO_VERIFIED_CHANGE')
    if summary.get('final_apbs_evidence',{}).get('status')!='APBS_EXECUTION_COMPLETE':raise ValueError('TARGET_APBS_AUDIT_MISSING')
    return {'status':'PASS','C_hard_all':actual,'constructs':len(ids)}
