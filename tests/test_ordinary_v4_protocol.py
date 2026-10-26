import pytest
from copy import deepcopy
from ecd_epitope_scanner.joint_design_v4.engine import _annotate_stability
from ecd_epitope_scanner.joint_design_v4.models import PredictionResult, V4Candidate, V4MutationOption, ResidueKey, V4Sector
from ecd_epitope_scanner.joint_design_v4.sequential_foldx import merge_repeats
from ecd_epitope_scanner.config import display_sectoring_parameters


def result(slots, value=.8):
    return PredictionResult(status='success', records=[{'mutations':['A1D','A2V'],'ddg':value,
        'paired_models':[{'sampling_slot':i,'ddg_kcal_mol':value,'batch_path':str(i),'group':0,'run':i} for i in slots]}])


def candidate():
    options=[V4MutationOption(ResidueKey('0','A',i,uniprot_pos=i,one_letter='A'),aa,f'A{i}{aa}', 'test',0.1,1.0,True) for i,aa in [(1,'D'),(2,'V')]]
    return V4Candidate('c','r',options,order=2,structural_risk=0)


@pytest.mark.parametrize('spurs,foldx,status',[(1,2,'CONSENSUS_PASS'),(.8,-8,'CONSENSUS_PASS'),(-8,1.8,'CONSENSUS_PASS'),(1.01,2,'CONSENSUS_REVIEW'),(1,2.01,'CONSENSUS_REVIEW'),(2.01,1,'COMBINATION_REJECT')])
def test_independent_gates_no_sign_or_gap_filter(spurs,foldx,status):
    c=candidate()
    multi=PredictionResult('success',[{'mutations':['A1D','A2V'],'ddg':spurs}])
    _annotate_stability([c],[],multi,result(range(5),foldx),{'stability_threshold_policy':'ordinary_v4_directional_v1','acceptable_max_ddg':1,'review_max_ddg':2},{'combination_gate_kcal_mol':2,'review_max_kcal_mol':3})
    assert c.stability.consensus_status == status
    assert c.impact_evidence['stability_repeat_audit']['foldx_evidence_complete']
    assert c.impact_evidence['stability_repeat_audit']['foldx_duplicate_content_count'] is None


def test_supplement_keeps_initial_pairs_and_rejects_duplicate_slots():
    first=result(range(3)); frozen=deepcopy(first.records)
    merged=merge_repeats(first,result([3,4],1.2))
    assert first.records == frozen
    assert len(merged.records[0]['paired_models']) == 5
    assert merged.records[0]['ddg'] == .8
    with pytest.raises(ValueError,match='unique'):
        merge_repeats(first,result([2,3]))


@pytest.mark.parametrize('length,target,minimum,maximum',[(119,7,6,8),(120,8,6,10),(150,8,6,10),(151,12,10,15)])
def test_length_boundaries(length,target,minimum,maximum):
    p=display_sectoring_parameters(length)
    assert (p['target_sector_residues'],p['preferred_minimum_sector_residues'],p['preferred_maximum_sector_residues']) == (target,minimum,maximum)

def test_sequential_budget_caps_and_failed_attempts_are_not_retried():
    from types import SimpleNamespace as NS
    from ecd_epitope_scanner.joint_design_v4.sequential_foldx import refine_foldx
    pool=[]; rows=[]
    for i in range(80):
        c=candidate(); c.candidate_id=f'c{i:03d}'
        c.mutations[0].notation=f'A{i+1}D'
        c.formal_eligible=True
        c.impact_evidence={'selection_rank':i+1}
        c.stability=NS(foldx_ddg_kcal_mol=.8)
        pool.append(c)
        rows.append({'mutations':[o.notation for o in c.mutations],'ddg':.8,
            'paired_models':[{'sampling_slot':j,'ddg_kcal_mol':.8} for j in range(3)]})
    calls=[]
    def predict(batch,index):
        calls.extend(c.candidate_id for c in batch)
        return PredictionResult('failed', errors=[{'code':'TEST_FAILURE'}])
    _, audit=refine_foldx(pool,[],PredictionResult('success',rows),
        choose_panel=lambda candidates:NS(selected_candidate_ids=[c.candidate_id for c in candidates[:8]]),
        predict=predict,recheck=lambda _:None,display_candidate_ids=['c079', 'unscored'])
    assert audit['total_slots'] == 300
    assert audit['supplemental_slots'] == 60
    assert len(calls) == len(set(calls)) == 30
    assert 'c079' in calls
    assert audit['display_candidates_not_refined'] == ['unscored']
    assert audit['display_candidate_budget_status'] == 'INCOMPLETE_FIRST_PASS_SELECTION'


def test_first_pass_recycles_global_shortage_without_duplicate_constructs():
    from ecd_epitope_scanner.joint_design_v4.candidate_search import select_v4_foldx_first_pass
    from types import SimpleNamespace as NS
    pool=[]
    for i in range(90):
        c=candidate(); c.candidate_id=f'c{i:03d}'
        c.mutations[0].notation=f'A{i+1}D'
        c.stability=NS(spurs_multi_ddg_kcal_mol=.2)
        pool.append(c)
    selected,audit=select_v4_foldx_first_pass(pool+[pool[0]],[],priority_candidate_ids=['c089'])
    assert len(selected) == 80
    assert len({tuple(sorted(o.notation for o in c.mutations)) for c in selected}) == 80
    assert 'c089' in {c.candidate_id for c in selected}
    assert audit['selected']['total'] == 80


def test_first_pass_priority_constructs_leave_room_for_local_exploration():
    from dataclasses import replace
    from types import SimpleNamespace as NS
    from ecd_epitope_scanner.joint_design_v4.candidate_search import select_v4_foldx_first_pass

    sector = V4Sector('local', list(range(1, 130)), 0, 1.0, 1.0, 0.0)
    pool = []
    for i in range(81):
        c = candidate()
        c.candidate_id = f'c{i:03d}'
        positions = (2 * i + 1, 2 * i + 2) if i < 61 else (300 + 2 * i, 301 + 2 * i)
        for option, position, mutant in zip(c.mutations, positions, ('D', 'V')):
            option.residue = replace(option.residue, uniprot_pos=position)
            option.notation = f'A{position}{mutant}'
        c.stability = NS(spurs_multi_ddg_kcal_mol=.8 if 45 <= i < 61 else .2)
        c.impact_evidence = {'selection_rank': i + 1}
        pool.append(c)

    selected, audit = select_v4_foldx_first_pass(
        pool, [sector], mechanism_percentile=1.0,
        priority_candidate_ids=[c.candidate_id for c in pool[:45]],
    )

    assert len(selected) == audit['selected']['total'] == 80
    assert audit['selected']['local'] == 60
    assert audit['selected']['global'] == 20
    assert audit['selected']['local_strong_exploration'] <= 16


def test_display_priority_covers_every_region_before_alternatives():
    from ecd_epitope_scanner.joint_design_v4.sequential_foldx import display_priority
    rows = [{'candidate_id': f'p{i}', 'alternatives': [
        {'candidate_id': f'a{i}'}, {'candidate_id': 'shared'}]} for i in range(12)]
    priority = display_priority(rows)
    assert priority[:12] == [f'p{i}' for i in range(12)]
    assert priority[12:24] == [f'a{i}' for i in range(12)]
    assert priority[24:] == ['shared']


def test_only_five_unique_finite_slots_enter_display_pool():
    from ecd_epitope_scanner.joint_design_v4.sequential_foldx import complete_display_pool
    pool = []
    for slots, values in [(list(range(5)), [1.] * 5), (list(range(3)), [1.] * 3),
                          ([0, 1, 2, 3, 3], [1.] * 5), (list(range(5)), [1.] * 4 + [float('nan')]),
                          (list(range(5)), [1.] * 4 + [None])]:
        c = candidate()
        c.impact_evidence = {'stability_repeat_audit': {
            'foldx_repeat_slots': slots, 'foldx_repeat_values_kcal_mol': values}}
        pool.append(c)
    assert complete_display_pool(pool) == [pool[0]]
