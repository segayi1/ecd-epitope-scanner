from types import SimpleNamespace as NS
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from ecd_epitope_scanner.pipeline import PipelineRequest
from ecd_epitope_scanner.joint_design_v4.stability_v3 import decision
from ecd_epitope_scanner.surface_perturbation.calibration import resolved_changes
from ecd_epitope_scanner.surface_perturbation.panel_v3 import metrics, compression_information, select
from ecd_epitope_scanner.surface_perturbation.sampling_v3 import validate_roles
from ecd_epitope_scanner.surface_perturbation.engine_v3 import coupling


@pytest.mark.parametrize('s,f,status',[(.3,-3,'CONSENSUS_PASS'),(.9,.4,'CONSENSUS_PASS'),
    (1.2,0,'CONSENSUS_REVIEW'),(2.1,-3,'COMBINATION_REJECT'),(.3,6,'STRUCTURAL_ENERGY_CONFLICT'),
    (.3,None,'SPURS_ONLY_REVIEW_REQUIRED'),(np.nan,0,'NOT_SCORED')])
def test_stability_primary(s,f,status):assert decision(s,f)==status


def test_profile_threshold_and_budget():
    r=PipelineRequest(input_type='uniprot',identifier='P16410',optimization_profile='coverage_scan_v3_preview')
    assert r.surface_strategy=='surface_perturbation_v1_preview'
    assert r.metadata['surface_rsasa_threshold']==.2 and r.maximum_constructs==12
    with pytest.raises(ValueError,match='Conflicting'):
        PipelineRequest(input_type='uniprot',optimization_profile='coverage_scan_v3_preview',metadata={
            'surface_rsasa_threshold':.2,'joint_v4_ml':{'surface':{'rsasa_threshold':.25}}})


def test_four_positive_one_unknown_is_lower_bound():
    c={'status':'PARTIAL','protocol_id':'p','valid':np.ones(3,bool),'threshold':np.ones(3),
       'positive_lower_bound':True}
    values=np.array([[2,2,0],[2,2,0],[2,2,0],[2,np.nan,0],[np.nan,np.nan,np.nan]])
    hit,unknown=resolved_changes(values,c,'p')
    assert hit.tolist()==[True,False,False] and unknown.tolist()==[False,True,False]
    hit,unknown=resolved_changes(values,{**c,'positive_lower_bound':False},'p')
    assert not hit.any() and unknown.all()


def test_unknown_code_never_becomes_zero():
    lower=np.array([[1,0,0]],bool);upper=np.array([[1,1,0]],bool)
    m=metrics(lower,upper,np.ones(3))
    assert m['geometric_pair_distinction']['1']==pytest.approx(1/3)
    assert m['worst_single_construct_loss_coverage']==0
    with pytest.raises(ValueError,match='lower outside'):
        metrics(lower,np.zeros_like(lower),np.ones(3))


def test_redundancy_and_loss_are_separate():
    lower=np.array([[1,1,0],[0,1,1]],bool)
    m=metrics(lower,lower,np.ones(3))
    assert m['double_covered_fraction']==pytest.approx(1/3)
    assert m['worst_single_construct_loss_coverage']==pytest.approx(2/3)


def test_compression_can_preserve_area_and_lose_information():
    lo={'A':np.array([1,0,0],bool),'B':np.array([0,1,0],bool),
        'C':np.array([1,1,0],bool)}
    audit=compression_information(['A','B'],['C'],lo,lo,np.ones(3))
    assert not audit['accepted'] and 'GEOMETRIC_DISTINCTION_LOSS' in audit['reasons']


def test_fixed_sampling_roles_not_content_uniqueness():
    rows=[dict(batch_path='b',group=1,run=i,sampling_protocol='time_rotabase_v3',
               sampling_role='production',sampling_cohort='c',sampling_slot=i,wt_sha256='same') for i in range(5)]
    assert len(validate_roles(rows,'production',5,'c'))==5
    assert validate_roles(rows,'null_train',10,'c')==[]
    with pytest.raises(ValueError,match='duplicate'):
        validate_roles(rows+[rows[0]],'production',5,'c')
    with pytest.raises(ValueError,match='fixed prefix'):
        validate_roles(rows[1:],'production',5,'c')


def test_select_four_and_deterministic():
    lo=np.eye(5,dtype=bool);graph=csr_matrix(np.ones((5,5))-np.eye(5))
    a,report=select(lo,lo,np.ones(5),graph,np.zeros(5),list('abcde'),solver_seconds=1)
    b,_=select(lo,lo,np.ones(5),graph,np.zeros(5),list('abcde'),solver_seconds=1)
    assert a==b and len(a)==4
    assert metrics(lo[a],lo[a],np.ones(5))['C_screen_lower_bound']==.8


def test_coupling_is_proxy_not_distance_gate():
    assert coupling({'minimum_native_heavy_atom_distance_a':4,'proxy_support_overlap_fraction':0})=='CONTACT_COUPLED_PROXY'
    assert coupling({'minimum_native_heavy_atom_distance_a':8,'proxy_support_overlap_fraction':0})=='SPATIALLY_SEPARATED_PROXY'
    assert coupling({})=='UNKNOWN'


def test_disjoint_null_and_untuned_holdout(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.null_v3 import build
    from ecd_epitope_scanner.surface_perturbation.null_store import FIELD_CHANNELS
    ref={'points':np.zeros((2,3)),'normals':np.tile([0,0,1],(2,1)),
         'area':np.ones(2),'owner':np.array([1,2]),'domain':np.array(['E','E']),'graph':csr_matrix((2,2))}
    train=[{name:np.full(2,.2) for name in FIELD_CHANNELS} for _ in range(5)]
    held=[{name:np.full(2,.1) for name in FIELD_CHANNELS} for _ in range(5)]
    for r in held:r['shape_sdf_a'][0]=2
    numerical={name:np.full(2,.1) for name in FIELD_CHANNELS}
    result=build(tmp_path,'p',train,held,numerical,np.ones(2,bool),ref,[])
    value=result['shape_sdf_a']
    assert value['threshold'].tolist()==[.2,.2]
    assert value['valid'].tolist()==[False,True]
    assert result['occupancy_gain_a']['valid'].all()
    # Corruption is rejected, not silently reused as new evidence.
    (tmp_path/'p.npz').write_bytes(b'corrupt')
    from ecd_epitope_scanner.surface_perturbation.null_store import load_calibration
    with pytest.raises(ValueError,match='fingerprint'):load_calibration(tmp_path,'p')


def test_surface_graph_does_not_jump_between_nearby_components():
    from ecd_epitope_scanner.surface_perturbation.geometry import Mesh
    from ecd_epitope_scanner.surface_perturbation.reference import canonical_reference
    from ecd_epitope_scanner.surface_perturbation.coverage import scan_support
    xyz=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
    faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
    mesh=Mesh(np.vstack([xyz,xyz+[0,0,5]]),np.vstack([faces,faces+4]))
    atom_map={('A',i,'','CA'):(np.array([.25,.25,.25+z]),'ALA','C') for i,z in ((1,0),(2,5))}
    ref=canonical_reference(mesh,atom_map,[NS(uniprot_pos=i,domain='E',chain_id='A',auth_seq_id=i,insertion_code='') for i in (1,2)],{1,2},{1,2})
    covered=scan_support(ref['graph'],ref['owner']==1,ref['owner'])
    assert covered[ref['owner']==1].all() and not covered[ref['owner']==2].any()
    assert np.linalg.norm(ref['points'][0]-ref['points'][-1])<12


def test_folded_graph_is_not_an_euclidean_ball():
    from ecd_epitope_scanner.surface_perturbation.coverage import scan_support
    graph=csr_matrix(np.array([[0,6,0,0],[6,0,6,0],[0,6,0,6],[0,0,6,0]],float))
    assert scan_support(graph,np.array([1,0,0,0],bool),np.arange(4)).tolist()==[True,True,True,False]


def test_euler_difference_keeps_valid_local_channel(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import fields
    wt={('A',1,'','CB'):(np.array([0.,0.,0.]),'ALA','C')}
    mut={('A',1,'','CB'):(np.array([.1,0.,0.]),'ALA','C')}
    ref={'points':np.array([[0.,0.,0.]]),'normals':np.array([[0.,0.,1.]]),
         'area':np.ones(1),'owner':np.ones(1,int),'domain':np.array(['E']),'graph':csr_matrix((1,1))}
    monkeypatch.setattr(fields,'aligned_atoms',lambda *a:(wt,mut,{'core_rmsd_a':0.,'rotation':np.eye(3),'translation':[0,0,0]}))
    base={'sdf':np.zeros(1),'shell_sdf':np.zeros((1,7)),'shell_valid':np.ones((1,7),bool),
          'shell_owned':np.ones((1,7),bool),'chemistry':np.zeros((1,20)),
          'projected_normals':np.array([[0,0,1.]]),'local_ses_area':np.ones(1),
          'atomic_inside':np.zeros((1,7)),'ses_total_area_a2':1.}
    monkeypatch.setattr(fields,'spatial_field',lambda *a,**k:base)
    monkeypatch.setattr(fields,'difference',lambda *a:{'shape_sdf_a':np.ones(1)})
    pair={'wt_sha256':'w','mut_sha256':'m'}
    tool=NS(config={'optimization_profile':'coverage_scan_v3_preview'},nano='nano',identity=lambda p:'hash',
            ses=lambda a:NS(euler_characteristic=2 if a is wt else 4))
    run=fields.compute_pair(pair,{},[],ref,tool,tmp_path)
    assert run['topology_changed'] and run['shape_sdf_a'][0]==1
    tool.config={}
    assert np.isnan(fields.compute_pair(pair,{},[],ref,tool,tmp_path)['shape_sdf_a'][0])


def test_site_witness_requires_the_same_four(monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import witness_v3
    monkeypatch.setattr(witness_v3,'FIELD_CHANNELS',('scalar',))
    monkeypatch.setattr(witness_v3,'VECTOR_GROUPS',{})
    runs=[{'scalar':np.array([2.,2.])} for _ in range(4)]+[{'scalar':np.array([np.nan,0.])}]
    support=np.ones((5,1,2),bool)
    support[3,0,1]=False
    calibration={'scalar':{'valid':np.ones(2,bool),'threshold':np.ones(2)}}
    result=witness_v3.witnesses(runs,support,calibration)
    assert result['site0__scalar'].tolist()==[15,0]


def test_backup_stage_preserves_original_case_and_controls(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import run_ctla4_coverage_v3 as runner
    from ecd_epitope_scanner.util import write_json
    case=tmp_path/'case';route=case/'routes/coverage_v3'
    ab='CV2_f1408e6b440bb101';ac='CV2_ad118b953b660407';bc='CV2_9d617576509297ec';abc='CV2_efd18628d844df64'
    write_json(route/'results/surface_perturbation_preview.json',{
        'candidate_local_calibration':{ab:{'positive_lower_bound_qualified':True}}})
    config={'candidate_ids':[ab,ac,bc,abc],'triple_candidate_id':abc,'triple_parent_ids':[ab,ac],
        'triple_role':'TRIPLE_EXPANSION','mandatory_candidate_signatures':[
            ['A61W','M134R'],['A61W','T67R'],['M134R','T67R']]}
    monkeypatch.setattr(runner.old,'triple_test',lambda *a,**k:dict(config))
    calls=[]
    def stage(case_arg,assets,route_arg,label,ids,**extra):
        calls.append((case_arg,route_arg,ids,extra))
        return {'engineering_status':'COMPLETED','stage':'T5','physical':{'hard_target_status':'BELOW_TARGET'}}
    monkeypatch.setattr(runner,'run_stage',stage)
    assert runner.triple_review(case,tmp_path,route,'T5')
    assert calls[0][0]==case and set(calls[0][2])=={ab,ac,bc,abc}
    assert calls[0][3]['mandatory_candidate_signatures']==config['mandatory_candidate_signatures']


def test_compression_protects_higher_code_distances():
    masks={str(i):np.array([1,0,0],bool) if i<3 else np.array([0,1,0],bool) for i in range(5)}
    result=compression_information(list(masks),['1','2','3','4'],masks,masks,np.ones(3))
    assert result['before']['geometric_pair_distinction']['1']==result['after']['geometric_pair_distinction']['1']
    assert result['before']['worst_single_construct_loss_coverage']==result['after']['worst_single_construct_loss_coverage']
    assert not result['accepted'] and 'GEOMETRIC_DISTINCTION_3_LOSS' in result['reasons']


def test_v3_resource_interface_does_not_reinstate_old_case_cap():
    from ecd_epitope_scanner.api import JobCreate
    r=JobCreate(input_type='uniprot',identifier='P16410',optimization_profile='coverage_scan_v3_preview',
                case_time_budget_seconds=100000,case_byte_budget=40*1024**3)
    assert PipelineRequest(**r.model_dump()).case_time_budget_seconds==100000
    with pytest.raises(ValueError):JobCreate(input_type='uniprot',identifier='P16410',case_time_budget_seconds=100000)


def test_risk_annotation_preserves_coupling_measurements():
    from ecd_epitope_scanner.joint_design_v4.risk import annotate_candidate_risks
    option=NS(spurs_ddg_kcal_mol=.1,glycan_risk=False,residue=NS(uniprot_pos=1))
    candidate=NS(mutations=[option],stability=None,risk_components={
        'geodesic_span_a':18.,'proxy_support_overlap_fraction':0.,'minimum_native_heavy_atom_distance_a':8.})
    annotate_candidate_risks([candidate],[],[])
    assert candidate.risk_components['geodesic_span_a']==18.
    assert coupling(candidate.risk_components)=='SPATIALLY_SEPARATED_PROXY'


def test_archive_rejects_rehashed_missing_execution_identity(tmp_path):
    import hashlib,json
    from zipfile import ZipFile
    from ecd_epitope_scanner.scripts.audit_ctla4_coverage_v3 import archive
    payload=json.dumps({'candidate':[{'sampling_role':'production'}]}).encode()
    path=tmp_path/'bad.zip'
    with ZipFile(path,'w') as z:
        z.writestr('portable_pairs.json',payload)
        z.writestr('evidence_manifest.json',json.dumps({'files':{'portable_pairs.json':hashlib.sha256(payload).hexdigest()}}))
    with pytest.raises(ValueError,match='execution identity'):archive(path,'evidence_manifest.json')


def test_early_stop_needs_a_formal_four_construct_witness():
    from ecd_epitope_scanner.surface_perturbation.panel_v3 import feasible_target
    graph=csr_matrix(np.ones((5,5))-np.eye(5))
    lower=np.eye(5,dtype=bool)
    assert feasible_target(lower[:3],lower[:3],np.ones(5),graph,np.zeros(3),list('abc')) is None
    assert len(feasible_target(lower,lower,np.ones(5),graph,np.zeros(5),list('abcde')))==4


def test_next_batch_size_uses_worst_measured_group(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v3 import next_batch_size
    from ecd_epitope_scanner.util import write_json
    write_json(tmp_path/'stages/V3_test_report.json',{'physical':{'cascade':{'a':{'physical_seconds':1500},'b':{'physical_seconds':30}}}})
    assert next_batch_size(tmp_path)==4


def test_resume_triple_keeps_its_exact_controls_and_role(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import run_ctla4_coverage_v3 as runner
    from ecd_epitope_scanner.util import write_json
    route=tmp_path/'routes/coverage_v3'
    write_json(route/'results/panel_candidates.json',[{'candidate_id':'ABC','order':3,
        'mutations':[{'notation':s} for s in ('A1R','E2W','Q3Y')],
        'parent_candidate_ids':['AB','AC'],'triple_role':'TRIPLE_EXPANSION'}])
    calls=[]
    def launch(case,assets,stage,config,seconds):
        calls.append(config);return {'engineering_status':'ENGINEERING_FAILURE'}
    monkeypatch.setattr(runner.old,'launch',launch)
    runner.run_stage(tmp_path,tmp_path,route,'resume',['ABC'])
    assert {'ABC','AB','AC'}<=set(calls[0]['candidate_ids'])
    assert len(calls[0]['mandatory_candidate_signatures'])==3
    assert calls[0]['triple_roles_by_candidate_id']=={'ABC':'TRIPLE_EXPANSION'}
def test_compression_does_not_resume_unrelated_deferred_doubles():
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v3 import triple_confirmation_ids
    summary={'selected_candidate_ids':['panel'], 'candidate_local_calibration':{
        'unselected_positive':{'positive_lower_bound_qualified':True},
        'deferred':{'positive_lower_bound_qualified':False}}}
    config={'triple_role':'TRIPLE_COMPRESSION','triple_candidate_id':'triple',
            'triple_parent_ids':['parent'], 'mandatory_candidate_signatures':[['A1W','M2R']],
            'candidate_ids':['panel','unrelated_deferred']}
    ids=triple_confirmation_ids(summary,config)
    assert {'panel','parent','triple'}<=set(ids)
    assert len(ids)==4
    assert 'deferred' not in ids and 'unselected_positive' not in ids
    config['triple_role']='TRIPLE_EXPANSION'
    assert 'unselected_positive' in triple_confirmation_ids(summary,config)


def test_completed_compression_receipt_recovers_once_without_scientific_calls(tmp_path,monkeypatch):
    import json
    from ecd_epitope_scanner.scripts import run_ctla4_coverage_v3 as runner
    route=tmp_path/'route';route.mkdir();stages=tmp_path/'stages';stages.mkdir()
    report={'stage':'V3_T6_COMPRESSION_1','engineering_status':'COMPLETED','case_elapsed_seconds':123.,
            'configuration':{'triple_candidate_id':'triple','triple_role':'TRIPLE_COMPRESSION'}}
    (stages/'V3_T6_COMPRESSION_1_report.json').write_text(json.dumps(report))
    (stages/'V3_T6_FAILED_report.json').write_text(json.dumps({**report,'engineering_status':'ENGINEERING_FAILURE','case_elapsed_seconds':999.}))
    ledger=tmp_path/'case_budget.json';ledger.write_text(json.dumps({'elapsed_seconds':123.,'active':False}))
    before=ledger.read_bytes();calls=[]
    monkeypatch.setattr(runner,'audit_current',lambda r:calls.append(('audit',r)))
    monkeypatch.setattr(runner.old,'record_completed_snapshot',lambda r,x:calls.append(('record',x['stage'])))
    runner.recover_completed_search_stage(tmp_path,route)
    runner.recover_completed_search_stage(tmp_path,route)
    assert ledger.read_bytes()==before and len(calls)==2
    assert json.loads((route/'v3_triple_attempts.json').read_text())==[
        {'candidate_id':'triple','role':'TRIPLE_COMPRESSION','remote':False,'stage':'V3_T6_COMPRESSION_1'}]


def test_compression_uses_frozen_achieved_double_panel():
    from types import SimpleNamespace as NS
    from ecd_epitope_scanner.surface_perturbation.panel_v3 import frozen_compression_baseline
    regular=[NS(candidate_id=str(i),order=2) for i in range(5)]
    lower={str(i):np.eye(5,dtype=bool)[i] for i in range(5)}
    ids=['0','1','2','3']
    assert [c.candidate_id for c in frozen_compression_baseline(ids,regular,lower,np.ones(5))]==ids
    with pytest.raises(ValueError,match='QUALIFIED_DOUBLE_PANEL'):
        frozen_compression_baseline(['0','1','2','missing'],regular,lower,np.ones(5))
    with pytest.raises(ValueError,match='BELOW_TARGET'):
        frozen_compression_baseline(ids,regular,{**lower,'3':np.zeros(5,bool)},np.ones(5))
    with pytest.raises(ValueError,match='QUALIFIED_DOUBLE_PANEL'):
        frozen_compression_baseline(['0','1','2','4'],regular[:-1]+[NS(candidate_id='4',order=3)],lower,np.ones(5))


def test_old_compression_reaudit_reuses_original_proposal_and_baseline(tmp_path,monkeypatch):
    import json
    from ecd_epitope_scanner.scripts import run_ctla4_coverage_v3 as runner
    from ecd_epitope_scanner.surface_perturbation.engine_v3 import COMPRESSION_CONTRACT
    route=tmp_path/'route';route.mkdir();(tmp_path/'stages').mkdir()
    config={'triple_candidate_id':'ABC','triple_role':'TRIPLE_COMPRESSION',
            'triple_parent_ids':['A','B'],'compression_baseline_ids':['A','B','C','D'],
            'mandatory_candidate_signatures':[], 'stage':'V3_T6_COMPRESSION_1'}
    report={'stage':config['stage'],'engineering_status':'COMPLETED','case_elapsed_seconds':1.,
            'configuration':config,'physical':{}}
    (tmp_path/'stages/V3_T6_COMPRESSION_1_report.json').write_text(json.dumps(report))
    calls=[]
    def stage(case,assets,r,label,ids,**extra):
        calls.append((ids,extra));return {'engineering_status':'COMPLETED','stage':'corrected'}
    monkeypatch.setattr(runner,'run_stage',stage)
    assert runner.compression_review(tmp_path,tmp_path,route,{'selected_candidate_ids':['unrelated_new']})
    assert set(calls[0][0])=={'A','B','C','D','ABC'}
    assert calls[0][1]['compression_baseline_ids']==['A','B','C','D']
    assert json.loads((route/'v3_triple_attempts.json').read_text())[0]['compression_contract']==COMPRESSION_CONTRACT


def test_final_confirmation_uses_only_fully_paid_positive_doubles():
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v3 import final_confirmation_ids
    from ecd_epitope_scanner.surface_perturbation.sampling_v3 import PROTOCOL
    models=[{'sampling_protocol':PROTOCOL,'sampling_role':role,'sampling_cohort':cohort,
             'sampling_slot':slot,'batch_path':role,'group':1,'run':slot}
            for role,cohort,count in [('production','production_v3',5),
                                     ('null_train','matched',10),('null_holdout','matched',10)]
            for slot in range(count)]
    history={cid:{'candidate':{'order':order},'candidate_local_calibration':{
        'positive_lower_bound_qualified':positive}} for cid,order,positive in
        [('paid',2,True),('cold',2,True),('negative',2,False),('triple',3,True)]}
    pairs={'paid':models,'cold':models[:1],'negative':models,'triple':models}
    assert final_confirmation_ids({'selected_candidate_ids':['selected']},history,pairs)==['paid','selected']

