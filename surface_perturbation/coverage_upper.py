"""Conservative upper bound under the actual five-model hard scanning rule."""
import numpy as np
from .coverage import scan_support


def hard_scan_upper_bound(run_support, graph):
    support=np.asarray(run_support,bool)
    if support.ndim!=2 or len(support)<5:
        raise ValueError('hard upper bound requires five actual paired-model support masks')
    if graph.shape!=(support.shape[1],support.shape[1]) or np.any(graph.data<0):
        raise ValueError('invalid canonical graph for hard upper bound')
    # A verified per-residue seed needs direct support in the same four runs.
    # Aggregating residues can only enlarge that seed set, never shrink it.
    possible_seed=np.sum(support[:5],axis=0)>=4
    return scan_support(graph,possible_seed,np.zeros(support.shape[1],int))
