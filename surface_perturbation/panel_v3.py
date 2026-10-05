"""Three-state geometric signatures; no antibody-response or IID assumptions."""
import numpy as np

from .coverage import largest_blind, nested_panel, independent_bounds


def metrics(lower, upper, area, graph=None, direct=None):
    lower, upper, area = np.asarray(lower,bool), np.asarray(upper,bool), np.asarray(area,float)
    if lower.ndim != 2 or lower.shape != upper.shape or lower.shape[1:] != area.shape:
        raise ValueError('panel matrix/reference mismatch')
    if not np.isfinite(area).all() or np.any(area<0) or area.sum()<=0 or np.any(lower & ~upper):
        raise ValueError('invalid area or lower outside proven upper')
    count = lower.sum(axis=0); total=float(area.sum())
    covered=count>0
    worst=min((float(area[np.any(np.delete(lower,i,axis=0),axis=0)].sum()/total)
               for i in range(len(lower))),default=0.)
    codes=np.where(lower,1,np.where(upper,-1,0)).T.astype(np.int8)
    groups,inverse=np.unique(codes,axis=0,return_inverse=True)
    weights=np.bincount(inverse,weights=area,minlength=len(groups))
    # Pair denominator excludes a quadrature sample paired with itself. Same
    # code groups contribute zero distinction, but stay in the denominator.
    denominator=(total*total-float(np.dot(area,area)))/2
    numer=np.zeros(3)
    block=max(1,min(128,16*1024**2//max(1,len(groups)*max(1,lower.shape[0])*4)))
    for start in range(0,len(groups),block):
        left=groups[start:start+block]
        opposite=((left[:,None,:]==1)&(groups[None,:,:]==0)
                  |(left[:,None,:]==0)&(groups[None,:,:]==1)).sum(axis=2)
        w=weights[start:start+block,None]*weights[None,:]
        valid=np.arange(start,start+len(left))[:,None]<np.arange(len(groups))[None,:]
        for i in range(3):numer[i]+=np.sum(w*valid*(opposite>=i+1))
    result={'C_screen_lower_bound':float(area[covered].sum()/total),
        'double_covered_fraction':float(area[count>=2].sum()/total),
        'worst_single_construct_loss_coverage':worst,
        'geometric_pair_distinction':{str(i+1):float(numer[i]/denominator) if denominator>0 else 0. for i in range(3)},
        'signature_groups':len(groups),'constructs':len(lower),
        'interpretation':'1 verified scan opportunity; 0 outside contract upper; ? unknown; not measured binding codes'}
    if graph is not None:result['largest_blind_patch_a2']=largest_blind(graph,covered,area)
    if direct is not None:
        direct=np.asarray(direct,bool)
        if direct.shape!=lower.shape:raise ValueError('direct matrix mismatch')
        result['C_physical_lower_bound']=float(area[np.any(direct,axis=0)].sum()/total)
    return result


def score(indices, lower, upper, area, graph, risk, ids):
    m=metrics(lower[indices],upper[indices],area,graph)
    return (-m['C_screen_lower_bound'],m['largest_blind_patch_a2'],
            -m['geometric_pair_distinction']['1'], -m['worst_single_construct_loss_coverage'],
            float(np.sum(np.asarray(risk)[indices])),tuple(sorted(ids[i] for i in indices)))


def select(lower,upper,area,graph,risk,ids,maximum=12,solver_seconds=5.):
    lower,upper=np.asarray(lower,bool),np.asarray(upper,bool)
    metrics(lower,upper,area)
    greedy=nested_panel(lower,area,risk,graph,ids,maximum)
    # Four constructs are an explicit experimental contract, even if three
    # already cover the target or another adds only redundancy.
    greedy += sorted((i for i in range(len(ids)) if i not in greedy),key=lambda i:ids[i])[:max(0,4-len(greedy))]
    best=list(greedy);bounds=[]
    for k in range(4,min(maximum,len(ids))+1):
        result=independent_bounds(lower,area,k,seconds=solver_seconds);bounds.append(result)
        proposal=result['selected_indices']
        proposal+=sorted((i for i in range(len(ids)) if i not in proposal),key=lambda i:ids[i])[:max(0,4-len(proposal))]
        if len(proposal)<=k and metrics(lower[proposal],upper[proposal],area)['C_screen_lower_bound']>=.8:
            best=proposal;break
        if k==min(maximum,len(ids)) and score(proposal,lower,upper,area,graph,risk,ids)<score(best,lower,upper,area,graph,risk,ids):
            best=proposal
    for _ in range(20):
        incumbent=score(best,lower,upper,area,graph,risk,ids);winner=None
        for old in sorted(best):
            for new in sorted(set(range(len(ids)))-set(best),key=lambda i:ids[i]):
                proposal=[new if i==old else i for i in best]
                candidate=score(proposal,lower,upper,area,graph,risk,ids)
                if candidate<incumbent:incumbent,winner=candidate,proposal
        if winner is None:break
        best=winner
    return best,{'status':'BEST_FOUND_FEASIBLE_PANEL','coverage_bounds':bounds,
                 'minimum_constructs_proven':bool(bounds and metrics(lower[best],upper[best],area)['C_screen_lower_bound']>=.8 and bounds[-1]['status']==0 and all(
                     r.get('upper_area_a2') is not None and r['upper_area_a2'] < .8*np.sum(area)
                     for r in bounds[:-1])), 'swap_round_limit':20}


def compression_information(before,after,lower,upper,area):
    a=metrics(np.array([lower[c] for c in before]),np.array([upper[c] for c in before]),area)
    b=metrics(np.array([lower[c] for c in after]),np.array([upper[c] for c in after]),area)
    reasons=[]
    for distance in ('1','2','3'):
        if b['geometric_pair_distinction'][distance]+1e-12<a['geometric_pair_distinction'][distance]:
            reasons.append('GEOMETRIC_DISTINCTION_LOSS' if distance=='1' else f'GEOMETRIC_DISTINCTION_{distance}_LOSS')
    if b['worst_single_construct_loss_coverage']+1e-12<a['worst_single_construct_loss_coverage']:
        reasons.append('SINGLE_CONSTRUCT_FAILURE_ROBUSTNESS_LOSS')
    return {'accepted':not reasons,'reasons':reasons,'before':a,'after':b}


def feasible_target(lower,upper,area,graph,risk,ids):
    """A verified feasible witness permits stopping; proxies never do."""
    lower=np.asarray(lower,bool)
    if len(ids)<4 or float(np.asarray(area)[lower.any(axis=0)].sum()/np.sum(area))<.8:return None
    indices,_=select(lower,upper,area,graph,risk,ids,solver_seconds=2.)
    return indices if 4<=len(indices)<=12 and metrics(lower[indices],np.asarray(upper)[indices],area)['C_screen_lower_bound']>=.8 else None


def frozen_compression_baseline(ids,regular,lower,area):
    """Use the declared achieved double panel, not a newly optimized pool."""
    lookup={c.candidate_id:c for c in regular}
    if (not 4<=len(ids)<=12 or len(set(ids))!=len(ids)
            or any(cid not in lookup or lookup[cid].order!=2 for cid in ids)):
        raise ValueError('FROZEN_COMPRESSION_BASELINE_NOT_A_QUALIFIED_DOUBLE_PANEL')
    if float(np.asarray(area)[np.any([lower[cid] for cid in ids],axis=0)].sum()/np.sum(area))<.8:
        raise ValueError('FROZEN_COMPRESSION_BASELINE_BELOW_TARGET')
    return [lookup[cid] for cid in ids]
