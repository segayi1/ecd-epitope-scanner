from types import SimpleNamespace
import numpy as np
from ecd_epitope_scanner.surface_perturbation.calibration import calibrate_channel
from ecd_epitope_scanner.surface_perturbation.null_store import FIELD_CHANNELS
from ecd_epitope_scanner.surface_perturbation.calibration_scope import local_scope
from ecd_epitope_scanner.surface_perturbation.calibration_schedule import null_append_needed,calibration_buckets
from ecd_epitope_scanner.scripts.continue_surface_validation import expansion_has_work


def partial_null(runs=20,scope=None):
    required=np.ones(3,bool) if scope is None else np.asarray(scope,bool)
    error=np.where(required,.1,np.nan)
    evidence={name:calibrate_channel(np.zeros((runs-1,3)),error,'matched',name,
               required_mask=required,area=np.ones(3)) for name in FIELD_CHANNELS}
    evidence['shape_sdf_a']['valid'][0]=False
    return evidence


def test_exact_capped_null_is_deferred_without_granting_scientific_completeness():
    mask=np.ones(3,bool);cal=partial_null()
    assert not local_scope(cal,mask,np.ones(3))['complete']
    assert not null_append_needed(cal,mask,np.ones(3))
    assert null_append_needed(partial_null(10),mask,np.ones(3))
    assert null_append_needed(partial_null(scope=[True,False,False]),mask,np.ones(3))
    assert null_append_needed({},mask,np.ones(3))


def test_unattempted_tasks_precede_retries_and_capped_cohorts_remain_unresolved():
    candidates=[SimpleNamespace(candidate_id=name) for name in ('capped','failed','fresh')]
    evidence={'capped':partial_null(),'failed':{},'fresh':{}}
    masks={c.candidate_id:np.ones(3,bool) for c in candidates}
    fresh,retry,deferred=calibration_buckets(candidates,evidence,masks,np.ones(3),{'failed'})
    assert [c.candidate_id for c in fresh]==['fresh']
    assert [c.candidate_id for c in retry]==['failed']
    assert deferred==['capped']
    assert not local_scope(evidence['capped'],masks['capped'],np.ones(3))['complete']


def test_unchanged_blindspot_requests_continue_when_batch_left_unattempted_work():
    requests=('candidate1','candidate2');previous=[requests]
    assert expansion_has_work(requests,previous,{'calibration_remaining_appendable_candidate_ids':['candidate2']})
    assert not expansion_has_work(requests,previous,{'calibration_remaining_appendable_candidate_ids':[]})
    assert expansion_has_work(('new',),previous,{})
    # A fresh coordinator must not retry an exhausted, exactly matched pool.
    exhausted={'calibration_requested_candidate_ids':list(requests),
               'calibration_remaining_appendable_candidate_ids':[]}
    assert not expansion_has_work(requests,[],exhausted)
    assert expansion_has_work(('new',),[],exhausted)
