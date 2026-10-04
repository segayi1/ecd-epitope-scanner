from types import SimpleNamespace as NS
import numpy as np
from ecd_epitope_scanner.surface_perturbation.combo_contract import site_evidence,goal_state,primary_count,waiting_confirmation


def test_installed_distribution_includes_physical_combination_module():
    from pathlib import Path
    import tomllib
    config=tomllib.loads((Path(__file__).parents[1]/'pyproject.toml').read_text(encoding='utf-8'))
    assert 'ecd_epitope_scanner.surface_perturbation' in config['tool']['setuptools']['packages']


def candidate(order=2,positions=(20,10),cid='c'):
    options=[NS(notation=f'A{p}V',mutant_aa='V',residue=NS(uniprot_pos=p,chain_id='A',auth_seq_id=p,insertion_code='')) for p in positions]
    return NS(order=order,mutations=options,candidate_id=cid)


def test_site_contribution_uses_native_identity_not_option_order():
    item=candidate()
    report=site_evidence(item,[[1,0],[0,1]],[2.,7.])
    assert report['complete']
    assert [(r['uniprot_pos'],r['verified_direct_area_a2']) for r in report['sites']]==[(10,2.),(20,7.)]


def test_silent_site_and_duplicate_position_are_not_valid_combinations():
    assert not site_evidence(candidate(),[[1,0],[0,0]],[1.,1.])['complete']
    assert not site_evidence(candidate(positions=(10,10)),[[1,0],[0,1]],[1.,1.])['complete']
    assert not site_evidence(candidate(),np.empty((0,2),bool),[1.,1.])['complete']


def test_target_expands_from_initial_budget_without_exceeding_twelve():
    curve=[{'C_hard_all':i/10} for i in range(1,13)]
    assert primary_count(curve,12,5)==8
    assert primary_count(curve,7,5)==7
    assert primary_count([{'C_hard_all':.7}]*12,12,5)==12


def test_physical_target_can_be_achieved_without_search_closure():
    items=[candidate(cid=str(i)) for i in range(8)]
    evidence={c.candidate_id:{'complete':True} for c in items}
    assert goal_state(items,.8,evidence,audit_complete=True)=='ACHIEVED'
    assert goal_state(items,.79,evidence,audit_complete=True)=='COVERAGE_BELOW_TARGET'
    assert goal_state(items,.8,evidence,audit_complete=False)=='FINAL_AUDIT_PENDING'
    assert goal_state(items,.8,evidence,audit_complete=True,fatal=True)=='INCOMPLETE_EVIDENCE'


def test_singles_or_missing_site_evidence_cannot_meet_panel_target():
    items=[candidate(cid=str(i)) for i in range(8)]
    evidence={c.candidate_id:{'complete':True} for c in items}
    items[0].order=1
    assert goal_state(items,.99,evidence,audit_complete=True)=='INCOMPLETE_SITE_EVIDENCE'
    items[0].order=2;evidence['0']['complete']=False
    assert goal_state(items,.99,evidence,audit_complete=True)=='INCOMPLETE_SITE_EVIDENCE'
    assert goal_state(items[:3],.99,evidence,audit_complete=True)=='INSUFFICIENT_VERIFIED_CONSTRUCTS'


def test_unmodeled_combinations_wait_but_real_failures_do_not():
    item=candidate();item.formal_exclusion_reasons=['SPURS_ONLY_REVIEW_REQUIRED','HIGHER_ORDER_REQUIRES_FOLDX']
    item.stability=NS(consensus_status='SPURS_ONLY_REVIEW_REQUIRED')
    assert waiting_confirmation(item)
    item.formal_exclusion_reasons+=['AUTO_EXCLUDED_MUTANT']
    assert not waiting_confirmation(item)
    item.formal_exclusion_reasons=['COMBINATION_REJECT'];item.stability.consensus_status='COMBINATION_REJECT'
    assert not waiting_confirmation(item)


def test_confirmation_proxy_diversity_does_not_depend_on_display_window_ids():
    from ecd_epitope_scanner.joint_design_v4.engine import _physical_confirmation_priority
    def probe(cid,positions,risk):
        return NS(candidate_id=cid,covered_positions_by_scenario={'geodesic_12':positions},structural_risk=risk,epitope_score=0.)
    items=[probe('a',[1,2],0.),probe('b',[1,2],.1),probe('c',[3],.2)]
    assert [c.candidate_id for c in _physical_confirmation_priority(items,2,{1:2.,2:2.,3:3.})]==['a','c']
    assert [c.candidate_id for c in _physical_confirmation_priority(items,1,{1:2.,2:2.,3:3.},{1,2})]==['c']


def test_triples_have_three_actual_distinct_pair_controls_without_duplicate_ids():
    from ecd_epitope_scanner.tests.test_joint_v4 import _probe
    from ecd_epitope_scanner.joint_design_v4.candidate_search import add_pair_controls
    triple=_probe('triple',3,(1,2,3),{})
    graph={1:{2:1.},2:{1:1.,3:1.},3:{2:1.}}
    result=add_pair_controls([triple],graph)
    pairs=[c for c in result if c.order==2]
    assert len(pairs)==3 and len({c.candidate_id for c in pairs})==3
    assert {tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in pairs}=={(1,2),(1,3),(2,3)}
    assert all(c.stability is None for c in pairs)
    assert len(add_pair_controls(result,graph))==4


def test_triple_cannot_use_incomplete_counterfactuals_or_single_union():
    from ecd_epitope_scanner.surface_perturbation.combo_contract import actual_pair_controls
    triple=candidate(3,(1,2,3),'triple')
    controls=[candidate(2,p,str(p)) for p in [(1,2),(1,3),(2,3)]]
    assert actual_pair_controls(triple,controls[:2]) is None
    assert actual_pair_controls(triple,[candidate(1,(p,),str(p)) for p in [1,2,3]]) is None
    # A numerically resolved zero-change pair is a legitimate counterfactual,
    # even when it cannot qualify as a final construct with two active sites.
    assert len(actual_pair_controls(triple,controls))==3


def test_global_span_does_not_cut_a_path_at_display_window_boundary():
    from ecd_epitope_scanner.tests.test_joint_v4 import _probe
    from ecd_epitope_scanner.joint_design_v4.models import V4Sector
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
    options=_probe('pair',2,(1,3),{}).mutations
    sector=V4Sector('S',[1,3],0,2.,2.,0.)
    graph={1:{2:1.},2:{1:1.,3:1.},3:{2:1.}}
    assert not generate_candidates([sector],options,graph,{'orders':[2]})
    pairs=generate_candidates([sector],options,graph,{'orders':[2],'global_span_paths':True})
    assert len(pairs)==1 and pairs[0].connected


def test_cached_qualification_is_recomputed_from_raw_null(monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import engine
    raw={'shape_sdf_a':{'valid':np.ones(1,bool),'threshold':np.zeros(1),
         'null_blocks':np.array([[0.],[0.],[0.],[10.]]),'numerical_error':np.zeros(1),
         'required_area':12.}}
    monkeypatch.setattr(engine,'usable_calibration',lambda *args:raw)
    calibrated=engine.load_calibration('ignored','p')['shape_sdf_a']
    assert not calibrated['valid'][0]
    assert calibrated['threshold'][0]>0 and calibrated['required_area']==12.
