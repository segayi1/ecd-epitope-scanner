import copy
import json
import numpy as np
import pytest
from ecd_epitope_scanner.surface_perturbation import local_constructs as local


def evidence(monkeypatch):
    checked=[]
    monkeypatch.setattr(local,'validate_pair',lambda p:checked.append(p['run']))
    candidate={'candidate_id':'c','order':1,'formal_eligible':True,'warnings':[],
        'mutations':[{'eligible':True,'hard_exclusions':[],'mutant_aa':'W','notation':'T82W',
                      'residue':{'one_letter':'T','chain_id':'A','auth_seq_id':45,
                                 'insertion_code':'','uniprot_pos':82}}],
        'stability':{'consensus_status':'CONSENSUS_PASS','foldx_ddg_kcal_mol':.1,
                     'spurs_multi_ddg_kcal_mol':.05}}
    pairs=[{'batch_path':'batch','group':1,'run':i,'mutations':['T46W'],
            'input_sha256':'input','normalized_input_sha256':'normalized',
            'binary_sha256':'binary','ddg_kcal_mol':.1} for i in range(5)]
    scope={'complete':True,'measurement_local_complete':True,'hard_scan_bound_resolved':False}
    return candidate,pairs,scope,checked


def test_local_qualified_construct_does_not_promote_global_panel(tmp_path,monkeypatch):
    candidate,pairs,scope,checked=evidence(monkeypatch)
    record=local.local_construct(candidate,pairs,scope,[2,3,4],[1,1,0],[1,0,0])
    assert checked==list(range(5))
    assert record['mutation_notations']==['T82W']
    assert record['normalized_model_mutations']==['T46W']
    assert record['residue_mapping'][0]['auth_seq_id']==45
    assert record['actual_changed_area_a2']==2
    assert record['hard_scan_center_area_a2']==5
    local.export_local_constructs(tmp_path,[record],'CALIBRATION_EXPANSION_REQUIRED')
    result=json.loads((tmp_path/'surface_perturbation_local_constructs.json').read_text())
    assert result['panel_status']=='CALIBRATION_EXPANSION_REQUIRED'
    assert result['constructs'][0]['status']=='VERIFIED_LOCAL_CONSTRUCT'
    assert 'T82W' in (tmp_path/'surface_perturbation_local_constructs.csv').read_text(encoding='utf-8-sig')


@pytest.mark.parametrize('failure',['four_models','duplicate_identity','partial_scope',
                                   'missing_measurement','stale_ddg','cysteine','input_mismatch'])
def test_missing_or_inconsistent_local_evidence_is_never_exported(monkeypatch,failure):
    candidate,pairs,scope,_=evidence(monkeypatch)
    if failure=='four_models':pairs=pairs[:4]
    if failure=='duplicate_identity':pairs[4]=copy.deepcopy(pairs[0])
    if failure=='partial_scope':scope['complete']=False
    if failure=='missing_measurement':scope['measurement_local_complete']=False
    if failure=='stale_ddg':candidate['stability']['foldx_ddg_kcal_mol']=.2
    if failure=='cysteine':candidate['mutations'][0]['mutant_aa']='C'
    if failure=='input_mismatch':pairs[4]['input_sha256']='different'
    with pytest.raises(ValueError):local.local_construct(candidate,pairs,scope,[2,3],[1,1],[1,0])


def test_upper_bound_can_close_local_scope_without_filling_unknown(monkeypatch):
    candidate,pairs,scope,_=evidence(monkeypatch)
    scope.update(complete=False,measurement_local_complete=False,hard_scan_bound_resolved=True)
    record=local.local_construct(candidate,pairs,scope,[2,3],[1,1],[1,0])
    assert record['hard_scan_bound_resolved']
    assert not scope['complete']


def test_compression_export_requires_own_role_marker_without_claiming_gain(monkeypatch):
    candidate,pairs,scope,_=evidence(monkeypatch)
    candidate.update(order=3,triple_role='TRIPLE_COMPRESSION',parent_candidate_ids=['ab','ac'])
    base=candidate['mutations'][0]
    candidate['mutations']=[{**base,'notation':f'T{p}W','residue':{**base['residue'],
        'uniprot_pos':p,'auth_seq_id':p}} for p in (82,83,84)]
    for pair in pairs:pair['mutations']=['T82W','T83W','T84W']
    with pytest.raises(ValueError,match='role qualification'):
        local.local_construct(candidate,pairs,scope,[2,3],[1,1],[1,0])
    candidate['warnings']=['REAL_TRIPLE_EXTRA_GAIN_VERIFIED']
    with pytest.raises(ValueError,match='role qualification'):
        local.local_construct(candidate,pairs,scope,[2,3],[1,1],[1,0])
    candidate['warnings']=['REAL_TRIPLE_COMPRESSION_VERIFIED']
    record=local.local_construct(candidate,pairs,scope,[2,3],[1,1],[1,0])
    assert record['triple_role']=='TRIPLE_COMPRESSION' and record['parent_candidate_ids']==['ab','ac']
    assert 'REAL_TRIPLE_EXTRA_GAIN_VERIFIED' not in record['warnings']
