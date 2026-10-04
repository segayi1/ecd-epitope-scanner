from types import SimpleNamespace as NS
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from ecd_epitope_scanner.surface_perturbation import engine,export,benchmark
from ecd_epitope_scanner.surface_perturbation.calibration import calibrate_channel
from ecd_epitope_scanner.surface_perturbation.null_store import FIELD_CHANNELS
from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
from ecd_epitope_scanner.tests.test_surface_perturbation import cube


@pytest.mark.parametrize('target',[5,20,'capped_expansion','combo_capped_initial'])
def test_calibration_can_recheck_same_repeat_count_without_new_models_or_fields(tmp_path,monkeypatch,target):
    combo=target=='combo_capped_initial'
    expanding=target=='capped_expansion'
    capped=expanding or combo
    target=20 if capped else target
    pdb=tmp_path/'native.pdb';write_atoms(pdb,{('A',1,'','CA'):(np.array([.5,.5,.5]),'ALA','C')})
    class Tools:
        nano='fixture';calls={'nanoshaper':0,'apbs':0,'pdb2pqr':0,'cache_hits':0}
        def __init__(self,*args):pass
        def ses(self,*args):return cube()
        def identity(self,*args):return 'fixture'
    reference={'points':np.array([[0.,0,0],[1.,0,0]]),'normals':np.tile([0.,0,1.],(2,1)),
        'area':np.ones(2),'owner':np.ones(2,int),'designable':np.ones(2,bool),
        'domain':np.array(['ECD']*2),'graph':csr_matrix((2,2)),
        'excluded_area_a2':0.,'internal_cavity_area_a2':0.}
    monkeypatch.setattr(engine,'PhysicalTools',Tools)
    monkeypatch.setattr(engine,'canonical_reference',lambda *args,**kwargs:reference)
    monkeypatch.setattr(engine,'validate_pair',lambda *args:None)
    calls=[]
    monkeypatch.setattr(engine,'compute_pair',lambda *args:calls.append(1) or {n:np.ones(2) for n in FIELD_CHANNELS})
    monkeypatch.setattr(engine,'complete_scope',lambda pair,raw,*args:{**raw,
        'query_scope_sha256':np.array('fixture'),'supplemented_query_points':np.array(0)})
    monkeypatch.setattr(engine,'protocol_identity',lambda *args:'fixture')
    monkeypatch.setattr(engine,'direct_support',lambda *args:np.ones(2,bool))
    monkeypatch.setattr(engine,'direct_site_masks',lambda *args:[np.ones(2,bool)])
    calibration={n:calibrate_channel(np.zeros((19 if capped else 4,2)),np.full(2,.1),'fixture',n) for n in FIELD_CHANNELS}
    if capped:calibration['shape_sdf_a']['valid'][0]=False
    monkeypatch.setattr(engine,'load_calibration',lambda *args:calibration)
    monkeypatch.setattr(engine,'cached_numerical_bounds',lambda *args,**kwargs:{})
    counts=[]
    monkeypatch.setattr(engine,'save_calibration',lambda folder,pairs,runs,*args,**kwargs:
        counts.append(len(runs)) or (calibration if len(runs)>=target else {}))
    monkeypatch.setattr(benchmark,'controls',lambda *args:{})
    monkeypatch.setattr(export,'evidence_archive',lambda *args:None)
    record=NS(chain_id='A',auth_seq_id=1,insertion_code='',uniprot_pos=1,coordinate_eligible=True,review_reasons=[],domain='ECD')
    option=NS(mutant_aa='V',notation='A1V',residue=NS(one_letter='A',uniprot_pos=1))
    candidate=NS(candidate_id='repeated',order=1,formal_eligible=True,structural_risk=0.,mutations=[option])
    if combo:
        candidate.order=2
        option.residue.chain_id='A';option.residue.auth_seq_id=1;option.residue.insertion_code=''
        candidate.mutations.append(NS(mutant_aa='V',notation='A2V',residue=NS(one_letter='A',uniprot_pos=2,chain_id='A',auth_seq_id=2,insertion_code='')))
    notations=[o.notation for o in candidate.mutations]
    pairs=[{'batch_path':str(tmp_path),'group':1,'run':i,'wt_sha256':'wt','mut_sha256':'mut','chain':'A'} for i in range(target)]
    checks=[];registered=[]
    output=tmp_path/'calibration' if expanding else tmp_path
    output.mkdir(exist_ok=True)
    if expanding:
        from ecd_epitope_scanner.util import write_json
        write_json(tmp_path/'production/surface_perturbation/blindspot_calibration_requests.json',{'candidate_ids':['repeated']})
    _,report=engine.run_physical_preview(structure=pdb,records=[record],options=[],candidates=[candidate],
        paired_records=[{'mutations':notations,'paired_models':pairs}],notation=lambda c:notations,
        supplement=lambda *args:(_ for _ in ()).throw(AssertionError('existing five repeats must be reused')),
        recheck=lambda c,p:checks.append(len(p)) or True,output=output,
        register_pairs=lambda p:registered.append(len(p)),settings={'phase':'calibration','calibrate_expansion_requests':expanding,
            'panel_contract':'combo_2_3_hard80_v1' if combo else None})
    if capped:
        assert report['status']=='CALIBRATION_INSUFFICIENT' and len(calls)==1 and counts==[]
        assert report['calibration_deferred_candidate_ids']==['repeated']
        assert report['calibration_selected_candidate_ids']==[]
        assert report['calibration_remaining_appendable_candidate_ids']==[]
        assert report['calibration_expansion_scope_converged'] is False
        assert checks==[1] and registered==[20]
        return
    assert report['status']=='LOCAL_CALIBRATION_COMPLETE'
    assert counts==([5] if target==5 else [5,10,20]) and len(calls)==1
    assert report['trajectories'][0]['calibration_runs']==target
    # Null extension retains the five-model stability/perturbation evidence.
    assert checks==([1,5,5] if target==5 else [1,5,5,5,5])
