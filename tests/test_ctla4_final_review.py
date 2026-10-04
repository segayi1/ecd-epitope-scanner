import hashlib
import numpy as np
import pytest

from ecd_epitope_scanner.scripts import run_ctla4_coverage_v2 as runner
from ecd_epitope_scanner.scripts.audit_ctla4_case import verify_hard_curve
from ecd_epitope_scanner.util import write_json


def test_paired_execution_identity_preserves_equal_bytes_and_rejects_duplicates():
    pair={'batch_path':'batch1','group':1,'run':0,'wt_sha256':'wt','mut_sha256':'mut'}
    second={**pair,'batch_path':'batch2'}
    one=runner.paired_execution_identities({'a':[pair]})
    two=runner.paired_execution_identities({'a':[pair,second]})
    assert len(two['a'])==2 and one!=two
    assert two==runner.paired_execution_identities({'a':[second,pair]})
    for duplicate in (pair,{**pair,'mut_sha256':'changed'}):
        with pytest.raises(ValueError,match='Duplicate'):
            runner.paired_execution_identities({'a':[pair,duplicate]})
    with pytest.raises(ValueError,match='Missing'):
        runner.paired_execution_identities({'a':[{**pair,'batch_path':''}]})


@pytest.mark.parametrize('new_batch',[False,True])
def test_final_apbs_zero_model_claim_checks_execution_batch(tmp_path,monkeypatch,new_batch):
    path=tmp_path/'results/surface_perturbation/paired_manifest.json'
    pair={'batch_path':'batch1','group':1,'run':0,'wt_sha256':'wt','mut_sha256':'mut'}
    write_json(path,{'a':[pair]})
    write_json(tmp_path/'results/surface_perturbation_preview.json',{'selected_candidate_ids':['a']})
    reports=[{'engineering_status':'COMPLETED','configuration':{'phase':'production'}}]
    monkeypatch.setattr(runner,'case_remaining',lambda c:1000)
    monkeypatch.setattr(runner,'stage_attempt_number',lambda *a:1)
    def launch(*args):
        assert args[3]['schedule']==[5]
        assert args[3]['calibration_extension_candidate_ids']==[]
        write_json(path,{'a':[{**pair,'batch_path':'batch2' if new_batch else 'batch1'}]})
        return {'stage':'final','engineering_status':'COMPLETED'}
    monkeypatch.setattr(runner,'launch',launch)
    audit=tmp_path/'final_apbs_model_reuse_audit.json'
    if new_batch:
        with pytest.raises(ValueError,match='changed paired model identities'):
            runner.final_panel_audit(tmp_path,tmp_path,reports)
        assert not audit.exists()
    else:
        runner.final_panel_audit(tmp_path,tmp_path,reports)
        assert runner.read(audit)['new_paired_models']==0


def test_incomplete_goal_still_requires_mask_summary_consistency():
    summary={'goal_status':'INSUFFICIENT_VERIFIED_CONSTRUCTS','selected_candidate_ids':['a'],
             'coverage_curve':[{'K':1,'actual_constructs':1,'C_hard_all':.5}]}
    area=np.ones(2)
    verify_hard_curve(summary,{'a':np.array([True,False])},area)
    with pytest.raises(ValueError,match='masks differ'):
        verify_hard_curve(summary,{'a':np.array([False,False])},area)
    summary['coverage_curve'][0]['actual_constructs']=0
    with pytest.raises(ValueError,match='masks differ'):
        verify_hard_curve(summary,{'a':np.array([True,False])},area)


def test_completed_triples_require_production_and_keep_roles(tmp_path):
    write_json(tmp_path/'screened_candidate_history.json',{
        'ab':{'candidate':{'mutations':[{'notation':'S1A'},{'notation':'S2W'}]}},
        'ac':{'candidate':{'mutations':[{'notation':'S1A'},{'notation':'R3E'}]}}})
    config={'phase':'production','triple_role':'TRIPLE_EXPANSION','triple_parent_ids':['ab','ac']}
    write_json(tmp_path/'stages/old_report.json',{'engineering_status':'COMPLETED','configuration':config})
    cid='CV2_'+hashlib.sha256('R3E|S1A|S2W'.encode()).hexdigest()[:16]
    write_json(tmp_path/'stages/cal_report.json',{'engineering_status':'COMPLETED',
        'configuration':{**config,'phase':'calibration','triple_candidate_id':'cal_only'}})
    write_json(tmp_path/'stages/failed_report.json',{'engineering_status':'ENGINEERING_FAILURE',
        'configuration':{**config,'triple_candidate_id':'failed'}})
    assert runner.completed_triple_roles(tmp_path)=={(cid,'TRIPLE_EXPANSION')}
    write_json(tmp_path/'stages/compression_report.json',{'engineering_status':'COMPLETED',
        'configuration':{**config,'triple_role':'TRIPLE_COMPRESSION','triple_candidate_id':cid}})
    assert runner.completed_triple_roles(tmp_path)=={(cid,'TRIPLE_EXPANSION'),(cid,'TRIPLE_COMPRESSION')}


def test_final_review_local_before_remote_and_does_not_repeat_failed_engineering(tmp_path,monkeypatch):
    write_json(tmp_path/'results/surface_perturbation_preview.json',{'goal_status':'INSUFFICIENT_VERIFIED_CONSTRUCTS'})
    calls=[]
    def triple(*args,remote=False):
        calls.append(remote)
        if calls==[False]:return [{'engineering_status':'COMPLETED'}]
        if remote:return [{'engineering_status':'ENGINEERING_FAILURE'}]
        return []
    monkeypatch.setattr(runner,'triple_test',triple)
    reports=runner.review_final_triples(tmp_path,tmp_path,{},[],100)
    assert calls==[False,False,True]
    assert reports[-1]['engineering_status']=='ENGINEERING_FAILURE'


def test_compression_stops_after_one_acceptance_and_at_minimum(tmp_path,monkeypatch):
    summary={'goal_status':'ACHIEVED','selected_candidate_ids':list('abcde')}
    path=tmp_path/'results/surface_perturbation_preview.json';write_json(path,summary)
    calls=[]
    def triple(*args,**kwargs):
        calls.append(1);write_json(path,{**summary,'compression_audits':[{'accepted':True}]})
        return [{'engineering_status':'COMPLETED'}]
    monkeypatch.setattr(runner,'triple_test',triple)
    assert len(runner.review_final_triples(tmp_path,tmp_path,{},[],100))==1
    assert calls==[1]
    write_json(path,{**summary,'selected_candidate_ids':list('abcd')})
    assert runner.review_final_triples(tmp_path,tmp_path,{},[],100)==[]
    assert calls==[1]


def test_triple_proposal_exact_controls_completed_exclusion_and_remote_cap(tmp_path,monkeypatch):
    def mutation(p,wt,aa):
        return {'notation':f'{wt}{p}{aa}','mutant_aa':aa,
                'residue':{'uniprot_pos':p,'one_letter':wt}}
    parents=[{'candidate_id':cid,'order':2,'mutations':mutations,
              'stability':{'consensus_status':'CONSENSUS_PASS'}}
             for cid,mutations in [('ab',[mutation(1,'S','A'),mutation(2,'S','W')]),
                                   ('ac',[mutation(1,'S','A'),mutation(3,'R','E')])]]
    output=tmp_path/'results'
    graph={'edges':[{'source':1,'target':2,'distance_a':10},
                    {'source':2,'target':3,'distance_a':10}]}
    write_json(output/'surface_graph_v4.json',graph)
    write_json(output/'uniprot.json',{'sequence':{'value':'SSR'}})
    write_json(output/'panel_candidates.json',parents)
    write_json(output/'candidate_pool_v2.json',{'candidates':[]})
    write_json(output/'surface_perturbation/paired_manifest.json',{'ab':[{}]*5,'ac':[{}]*5})
    write_json(output/'surface_perturbation_preview.json',{'goal_status':'INSUFFICIENT_VERIFIED_CONSTRUCTS'})
    monkeypatch.setattr(runner,'read_budget',lambda p:{'elapsed_seconds':4000})
    monkeypatch.setattr(runner,'preserve_pair_index',lambda c:None)
    monkeypatch.setattr(runner,'record_completed_snapshot',lambda c,r:None)
    calls=[]
    def launch(c,a,stage,config,seconds):
        calls.append(config)
        report={'stage':stage,'engineering_status':'COMPLETED','configuration':config}
        write_json(c/'stages'/(stage+'_report.json'),report)
        return report
    monkeypatch.setattr(runner,'launch',launch)
    reports=runner.triple_test(tmp_path,tmp_path,{},['ab','ac'],100)
    assert len(reports)==2
    controls=[list(x) for x in [('S1A','S2W'),('R3E','S1A'),('R3E','S2W')]]
    assert sorted(map(sorted,calls[0]['mandatory_candidate_signatures']))==sorted(map(sorted,controls))
    assert calls[0]['triple_role']=='TRIPLE_EXPANSION'
    assert runner.triple_test(tmp_path,tmp_path,{},['ab','ac'],100)==[]
    assert len(calls)==2
    graph['edges'][1]['distance_a']=20
    write_json(output/'surface_graph_v4.json',graph)
    # A distinct role may be audited; a remote-only pool never disguises local layout.
    for p in (tmp_path/'stages').glob('*_report.json'):
        r=runner.read(p);r['configuration']['triple_role']='TRIPLE_COMPRESSION';write_json(p,r)
    write_json(tmp_path/'remote_triple_confirmation_ids.json',['old1','old2'])
    assert runner.triple_test(tmp_path,tmp_path,{},['ab','ac'],100,remote=True)==[]
    assert len(calls)==2
