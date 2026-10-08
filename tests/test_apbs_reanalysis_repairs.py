import json
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from ecd_epitope_scanner.scripts.reanalysis_apbs_preparation_fix import surface_distances, replay_saved
from ecd_epitope_scanner.scripts.reconcile_case_ledger import reconcile
from ecd_epitope_scanner.scripts.prepare_cross_protein_validation import protein_chains
from ecd_epitope_scanner.surface_perturbation import electrostatics as es
from ecd_epitope_scanner.surface_perturbation.paired import atoms


@pytest.mark.parametrize('missing_anchor', [False, True])
def test_terminal_rotation_insertion_keeps_conditional_qc(tmp_path, monkeypatch, missing_anchor):
    native={('A',56,'B','CA'):(np.array([-2.,0.,0.]),'ASN','C'),
            ('A',56,'B','CG'):(np.zeros(3),'ASN','C'),
            ('A',56,'B','OD1'):(np.array([1.,1.,0.]),'ASN','O'),
            ('A',56,'B','ND2'):(np.array([1.,-1.,0.]),'ASN','N')}
    monkeypatch.setattr(es,'aligned_atoms',lambda *args:(native,native,{'mutated_author_residues':[]}))
    class Tools:
        def prepare_pqr(self, path):
            values={(k[0],k[1],k[3]):{'xyz':v[0].tolist(),'residue':v[1],'charge':0.,'radius':1.7}
                    for k,v in atoms(path).items()}
            values['A',1,'OD1']['xyz']=[1.,-1.,0.]
            values['A',1,'ND2']['xyz']=[1.,1.,0.]
            if missing_anchor:values.pop(('A',1,'CA'))
            return path,values
    _, audits, _=es.prepare_candidate('flip',[{}]*5,native,[],Tools(),tmp_path)
    for row in audits:
        assert row['heavy_atom_preparation_changes']
        diagnostics=row['sidechain_rotation_diagnostics']
        if missing_anchor:assert diagnostics==[]
        else:
            assert len(diagnostics)==2
            assert diagnostics[0]['residue']==['A',56,'B']
            assert diagnostics[0]['electrostatic_equivalence']=='NOT_ESTABLISHED'


def test_surface_distances_are_per_point_and_disconnected():
    graph=csr_matrix([[0,2,0,0],[2,0,3,0],[0,3,0,0],[0,0,0,0]])
    owner=np.array([10,20,30,40])
    assert np.allclose(surface_distances(graph,owner,[10]),[0,2,5,np.inf])
    assert np.allclose(surface_distances(graph,owner,[10,30]),[0,2,0,np.inf])
    assert np.isinf(surface_distances(graph,owner,[99])).all()
    with pytest.raises(ValueError):surface_distances(graph,np.array([10]),[10])


def test_saved_replay_uses_error_and_solvent_not_finiteness():
    data=dict(raw_preparation_conditional_offset_deltas=np.full((5,3,3),2.),
              numerical_error=np.array([.1,np.nan,.1]),solvent_valid=np.ones((5,3),bool),
              null_blocks=np.zeros((4,3)),medium_evaluated=np.array(True))
    data['solvent_valid'][:,2]=False
    report={'preparation_gate_rejected':False,'evidence_states':{'APBS_PREPARATION_VALID':True}}
    _,_,changed,unknown,reasons=replay_saved(data,report,'p')
    assert changed.tolist()==[True,False,False] and unknown.tolist()==[False,True,True]
    assert reasons['MISSING_OR_INVALID_NUMERICAL_ERROR']==1
    report['preparation_gate_rejected']=True
    assert not replay_saved(data,report,'p')[2].any()
    assert replay_saved(data,report,'p')[3].all()


def test_ledger_reconciliation_preserves_unknown_allocation(tmp_path):
    (tmp_path/'case_budget.json').write_text(json.dumps({'active':False,'elapsed_seconds':15.}))
    (tmp_path/'stages').mkdir()
    (tmp_path/'stages/a_report.json').write_text(json.dumps(
        {'stage':'a','case_elapsed_seconds':12.,'charged_seconds':10.}))
    result=reconcile(tmp_path)
    assert result['receipt_seconds']==10 and result['separately_reported_residual_seconds']==5
    assert result['accounting_error_seconds']==0 and result['allocation_status']=='PARTIAL'
    assert result['finalization_tail_seconds']==3
    (tmp_path/'case_budget.json').write_text(json.dumps({'active':True,'elapsed_seconds':15.}))
    with pytest.raises(ValueError,match='closed'):reconcile(tmp_path)


def test_polymer_chain_context_excludes_glycans():
    assert protein_chains({'_entity_poly.type':['polypeptide(L)','polysaccharide(D)'],
                           '_entity_poly.pdbx_strand_id':['A, Z','B,C']})==['A','Z']
