from __future__ import annotations

import json
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from ecd_epitope_scanner.surface_perturbation.geometry import Mesh, triangle_distance, align_pair
from ecd_epitope_scanner.surface_perturbation.calibration import calibrate_channel, resolved_changes
from ecd_epitope_scanner.surface_perturbation.coverage import budget_curve, nested_panel, independent_bounds, scan_support
from ecd_epitope_scanner.surface_perturbation.runtime import RunBudget


def cube():
    vertices=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    faces=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                    [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return Mesh(vertices,faces)


def test_triangle_distance_not_nearest_vertex():
    triangle=np.array([[[0,0,0],[10,0,0],[0,10,0]]],float)
    assert triangle_distance(np.array([2,2,1]),triangle)[0] == pytest.approx(1)


def test_mesh_deadline_preserves_exact_distances_and_interrupts():
    from ecd_epitope_scanner.surface_perturbation.deadlines import guard_mesh
    class Budget:
        elapsed=0
        limit=1
        def remaining(self):
            if self.elapsed>=self.limit:raise TimeoutError('exhausted')
            return self.limit-self.elapsed
    points=np.array([[.4,.3,.2],[2,.3,.2],[0,.3,.2]])
    expected=cube().signed_distance(points)
    budget=Budget(); guarded=guard_mesh(cube(),budget)
    actual=guarded.signed_distance(points)
    for a,b in zip(actual,expected):np.testing.assert_array_equal(a,b)
    budget.elapsed=2
    with pytest.raises(TimeoutError):guarded.signed_distance(points)


def test_parallel_potentials_deduplicate_and_preserve_pair_order(tmp_path,monkeypatch):
    import threading
    from ecd_epitope_scanner.surface_perturbation import execution
    monkeypatch.setattr(execution,'available_memory',lambda:32*1024**3)
    a=tmp_path/'a.pqr';b=tmp_path/'b.pqr';c=tmp_path/'same.pqr'
    a.write_text('a');b.write_text('b');c.write_text('a')
    class Tools:
        config={'apbs_workers':2}
        activity=[]
        calls=[]
        barrier=threading.Barrier(2)
        def potential(self,path,grid):
            self.calls.append(path)
            if path==b:self.barrier.wait(timeout=5)
            elif len(self.calls)<=2:self.barrier.wait(timeout=5)
            return path.read_text()
    tools=Tools();grid={'counts':[33,33,33]}
    assert execution.potential_pair(tools,a,b,grid)==('a','b')
    tools.calls=[]
    # No second worker/job for identical content.
    tools.potential=lambda path,grid:tools.calls.append(path) or path.read_text()
    assert execution.potential_pair(tools,a,c,grid)==('a','a')
    assert len(tools.calls)==1 and tools.activity[-1]['unique_jobs']==1


def test_exact_solvent_query_cache_invalidates_input_and_detects_corruption(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import execution
    class Tools:
        cache=tmp_path
        nano='nano'
        calls={'cache_hits':0}
        activity=[]
        def identity(self,*args):return 'a'*64
    measured=[]
    monkeypatch.setattr(execution,'pqr_surface',lambda *args:measured.append(1) or cube())
    tools=Tools();atoms={('A',1,'CA'):{'xyz':[0,0,0],'radius':1.7}}
    points=np.array([[.4,.3,.2],[2,.3,.2]])
    expected=execution.solvent_distances(tools,atoms,points)
    actual=execution.solvent_distances(tools,atoms,points)
    for a,b in zip(actual,expected):np.testing.assert_array_equal(a,b)
    assert len(measured)==1
    execution.solvent_distances(tools,atoms,points+.1)
    assert len(measured)==2
    data=next(tmp_path.glob('solvent_queries_*/distances.npz'));data.write_bytes(b'corrupt')
    with pytest.raises(ValueError):execution.solvent_distances(tools,atoms,points)


def test_nonclosed_pqr_surface_is_cached_as_unknown_with_reason(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import execution
    class Tools:
        cache=tmp_path
        nano='nano'
        calls={'cache_hits':0}
        activity=[]
        def identity(self,*args):return 'a'*64
    measured=[]
    def nonclosed(*_args):
        measured.append(1)
        raise RuntimeError('PQR SES is not closed')
    monkeypatch.setattr(execution,'pqr_surface',nonclosed)
    tools=Tools();atoms={('A',1,'CA'):{'xyz':[0,0,0],'radius':1.7}}
    points=np.array([[.4,.3,.2],[2,.3,.2]])
    first=execution.solvent_distances(tools,atoms,points,return_audit=True)
    second=execution.solvent_distances(tools,atoms,points,return_audit=True)
    for left,right in zip(first[:2],second[:2]):np.testing.assert_array_equal(left,right)
    assert first[2]==second[2]=={'status':'UNRESOLVED','reason':'PQR SES is not closed'}
    assert np.isnan(first[0]).all() and not first[1].any()
    assert measured==[1] and tools.calls['cache_hits']==1


def test_whole_triangle_cache_preserves_projection_unknown_and_mesh_identity(tmp_path):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation.query_cache import memoize_mesh
    tools=SimpleNamespace(cache=tmp_path,calls={'cache_hits':0},activity=[])
    points=np.array([[.4,.3,.2],[2,.3,.2],[0,.3,.2],[1.,1.,1.]])
    expected=cube().signed_distance(points,return_projection=True)
    measured=memoize_mesh(cube(),tools).signed_distance(points,return_projection=True)
    cached=memoize_mesh(cube(),tools).signed_distance(points,return_projection=True)
    for left,right in zip(expected,cached):np.testing.assert_array_equal(left,right)
    for left,right in zip(measured,cached):np.testing.assert_array_equal(left,right)
    assert tools.calls['cache_hits']==1
    moved=Mesh(cube().vertices+[3,0,0],cube().faces)
    actual=memoize_mesh(moved,tools).signed_distance(points,return_projection=True)
    for left,right in zip(actual,moved.__class__(moved.vertices,moved.faces).signed_distance(points,return_projection=True)):
        np.testing.assert_array_equal(left,right)
    assert tools.calls['cache_hits']==1
    partial=Mesh(cube().vertices,cube().faces[:-1])
    result=memoize_mesh(partial,tools).signed_distance(points,return_projection=True)
    again=memoize_mesh(partial,tools).signed_distance(points,return_projection=True)
    assert not again[1].any() and np.isnan(again[0]).all()
    for left,right in zip(result,again):np.testing.assert_array_equal(left,right)


def test_exact_point_state_cache_reuses_subset_but_computes_whole_partial_miss(tmp_path):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation.query_cache import memoize_mesh
    tools=SimpleNamespace(cache=tmp_path,calls={'cache_hits':0},activity=[])
    mesh=cube();original=mesh.signed_distance;calls=[]
    def measured(points,**kwargs):
        calls.append(len(points));return original(points,**kwargs)
    mesh.signed_distance=measured;mesh=memoize_mesh(mesh,tools)
    points=np.random.default_rng(17).uniform(-.2,1.2,(64,3))
    mesh.signed_distance(points,return_projection=True)
    subset=points[8:24][::-1]
    # BLAS may round a reversed subset differently; exact-state caches must miss then.
    same_state=((points@mesh.ray_basis)[8:24][::-1].tobytes()==(subset@mesh.ray_basis).tobytes()
        and mesh.tree.query(points)[1][8:24][::-1].tobytes()==mesh.tree.query(subset)[1].tobytes())
    expected=cube().signed_distance(subset,return_projection=True)
    actual=mesh.signed_distance(subset,return_projection=True)
    for left,right in zip(actual,expected):np.testing.assert_array_equal(left,right)
    expected_calls=[64] if same_state else [64,16]
    assert calls==expected_calls and tools.calls['cache_hits']==int(same_state)
    mixed=np.array([points[0],[2.1,.13,.27]])
    actual=mesh.signed_distance(mixed,return_projection=True)
    for left,right in zip(actual,cube().signed_distance(mixed,return_projection=True)):np.testing.assert_array_equal(left,right)
    assert calls==expected_calls+[2] # The original batch path, not a one-point miss query.
    bank=next((tmp_path/'triangle_queries').glob('points_*.npz'));bank.write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='point cache changed'):mesh.signed_distance(subset,return_projection=True)


def test_current_policy_invalidates_historical_verified_panel(tmp_path):
    import json
    from ecd_epitope_scanner.surface_perturbation.case_status import case_status
    summary = {'status': 'VERIFIED_PANEL', 'nested_verified_candidate_ids': ['old']}
    (tmp_path/'surface_perturbation').mkdir()
    (tmp_path/'surface_perturbation_preview.json').write_text(json.dumps(summary))
    manifest = tmp_path/'surface_perturbation/paired_manifest.json'
    manifest.write_text(json.dumps({'old': [{'mutations': ['S38C']}]}))
    assert case_status(tmp_path, summary) == 'INCOMPLETE_MUTATION_POLICY'
    manifest.write_text(json.dumps({'old': [{'mutations': ['S38A']}]}))
    assert case_status(tmp_path, summary) == 'VERIFIED_PANEL'


def test_evidence_audit_rejects_corruption_and_preserves_case_status(tmp_path):
    from zipfile import ZipFile
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    from ecd_epitope_scanner.surface_perturbation.archive_audit import validate_archive,safe_name
    folder=tmp_path/'surface_perturbation';folder.mkdir()
    (folder/'reference_contract.json').write_text('{"reference":"native"}')
    (tmp_path/'surface_perturbation_preview.json').write_text('{"status":"INCOMPLETE"}')
    archive=evidence_archive(tmp_path,{},tmp_path/'null')
    report=validate_archive(archive)
    assert report['status']=='PASS' and report['case_status']=='INCOMPLETE'
    damaged=tmp_path/'bad.zip'
    with ZipFile(archive) as source,ZipFile(damaged,'w') as target:
        for name in source.namelist():target.writestr(name,b'bad' if name=='reference_contract.json' else source.read(name))
    with pytest.raises(ValueError,match='hash changed'):validate_archive(damaged)
    for path in ('../x','C:/x','/x','a\\x'):
        with pytest.raises(ValueError):safe_name(path)


def test_legacy_pair_archive_includes_only_matching_original_ddg(tmp_path,monkeypatch):
    from zipfile import ZipFile
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    from ecd_epitope_scanner import stability
    batch=tmp_path/'batch';batch.mkdir()
    for name in ('wt.pdb','input_1.pdb','PdbList.fxout','Dif_input.fxout','Dif_other.fxout'):
        (batch/name).write_text(name)
    monkeypatch.setattr(stability,'_parse_dif_rows',lambda p:[('input_1.pdb',-.2)] if p.name=='Dif_input.fxout' else [('other.pdb',1.)])
    row={'batch_path':str(batch),'wt_path':str(batch/'wt.pdb'),'mut_path':str(batch/'input_1.pdb'),
        'pdb_list_path':str(batch/'PdbList.fxout'),'ddg_kcal_mol':-.2}
    with ZipFile(evidence_archive(tmp_path,{'c':[row]},tmp_path/'null')) as archive:
        pair=json.loads(archive.read('portable_pairs.json'))['c'][0]
        assert archive.read(pair['ddg_file_path'])==b'Dif_input.fxout'
        assert pair['ddg_file_sha256']
        assert not any(name.endswith('Dif_other.fxout') for name in archive.namelist())


def test_final_evidence_requires_all_audits_but_never_converts_unknown_to_coverage():
    from ecd_epitope_scanner.surface_perturbation.final_evidence import require_apbs_audits
    with pytest.raises(ValueError,match='INCOMPLETE'):require_apbs_audits(['a','b'],['a'])
    with pytest.raises(ValueError,match='NOT_REQUESTED'):require_apbs_audits(['a'],[],requested=False)
    report=require_apbs_audits(['a','b'],['b','a'])
    assert report['status']=='APBS_EXECUTION_COMPLETE' and 'unknown electrostatics' in report['interpretation']
    assert report['evidence_states']['APBS_EXECUTION_COMPLETE']
    assert not report['evidence_states']['EXTENDED_COVERAGE_RESOLVED']
    assert not any('coverage'==key for key in report)


def test_rejected_preparation_retains_five_solves_and_raw_values_without_medium_claim(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import electrostatics as es
    from ecd_epitope_scanner.surface_perturbation.paired import atoms
    source={('A',1,'','CA'):(np.array([.5,.3,.3]),'ALA','C')}
    mutant={('A',1,'','CA'):(np.array([1.5,.3,.3]),'VAL','C')}
    monkeypatch.setattr(es,'aligned_atoms',lambda *args:(source,mutant,{'mutated_author_residues':[('A',1,'')]}))
    monkeypatch.setattr(es,'pqr_surface',lambda *args:pytest.fail('rejected preparation must not run useless SES queries'))
    class Tools:
        config={'apbs_workers':1}
        solves=0
        def prepare_pqr(self,path):
            values={(k[0],k[1],k[3]):{'xyz':(v[0]+[3,0,0]).tolist(),'residue':v[1],'charge':0.,'radius':1.7}
                for k,v in atoms(path).items()}
            return path,values
        def potential(self,path,grid):
            self.solves+=1
            return lambda points:np.full(len(points),100. if 'Mut' in path.name else 0.)
    tools=Tools();reference={'points':np.array([[.5,.3,.3]]),'normals':np.array([[1.,0,0]])}
    changed,unknown,report=es.audit_candidate('rejected',[{}]*5,source,[],reference,tools,tmp_path)
    assert tools.solves==10 and unknown.all() and not changed.any()
    assert report['preparation_gate_rejected'] and report['skipped_unusable_SES_query_batches']==10
    with np.load(tmp_path/'rejected_electrostatics.npz') as field:
        assert not field['medium_evaluated'] and np.isnan(field['delta_phi']).all()
        assert (field['raw_preparation_conditional_offset_deltas']==100.).all()


def test_interrupted_checkpoint_preserves_last_readable_cumulative_ledger(tmp_path,monkeypatch):
    from pathlib import Path
    path=tmp_path/'budget.json';budget=RunBudget(path,7200)
    original=json.loads(path.read_text())
    with monkeypatch.context() as scoped:
        scoped.setattr(Path,'replace',lambda *args:(_ for _ in ()).throw(OSError('interrupted atomic checkpoint')))
        with pytest.raises(OSError):budget.save()
    assert json.loads(path.read_text())==original
    budget.save(False)
    assert not json.loads(path.read_text())['active']


def test_transient_windows_budget_reader_does_not_fail_the_case(tmp_path,monkeypatch):
    from pathlib import Path
    path=tmp_path/'budget.json';budget=RunBudget(path,7200);actual_replace=Path.replace;calls=[]
    def transient(source,target):
        calls.append(1)
        if len(calls)<3:raise PermissionError('transient reader')
        return actual_replace(source,target)
    with monkeypatch.context() as scoped:
        scoped.setattr(Path,'replace',transient);budget.save()
    assert len(calls)==3 and json.loads(path.read_text())['active']
    budget.save(False)


def test_whole_pipeline_export_time_cannot_pass_cumulative_deadline(tmp_path,monkeypatch):
    from ecd_epitope_scanner import pipeline
    from ecd_epitope_scanner.surface_perturbation.runtime import ACTIVE_BUDGET
    def completed_physics(*args):
        budget=ACTIVE_BUDGET.get();budget.carried=budget.limit+1
        return {'fixture':'physical stage complete, export exceeded budget'}
    monkeypatch.setattr(pipeline,'_run_pipeline',completed_physics)
    request=pipeline.PipelineRequest(input_type='structure',surface_strategy='surface_perturbation_v1_preview',output_dir=str(tmp_path))
    with pytest.raises(TimeoutError):pipeline.run_pipeline(request)
    verdict=json.loads((tmp_path/'case_acceptance.json').read_text())
    assert verdict['status']=='INCOMPLETE_RUNTIME_EXHAUSTED'
    assert not json.loads((tmp_path/'surface_production_cumulative_budget.json').read_text())['active']


def test_final_bundle_preserves_hashes_without_recompressing_nested_zip(tmp_path):
    import hashlib
    from zipfile import ZipFile,ZIP_STORED
    from ecd_epitope_scanner.exporting import finalize_scientific_bundle
    binary=b'PK engineering fixture';(tmp_path/'nested.zip').write_bytes(binary)
    (tmp_path/'analysis.json').write_text('{"status":"fixture"}')
    (tmp_path/'run_manifest.json').write_text(json.dumps({'files':{'analysis.json':'old','nested.zip':'old'}}))
    finalize_scientific_bundle(tmp_path)
    with ZipFile(tmp_path/'ecd_epitope_design.zip') as archive:
        manifest=json.loads(archive.read('run_manifest.json'))
        for name,digest in manifest['files'].items():assert hashlib.sha256(archive.read(name)).hexdigest()==digest
        assert archive.getinfo('nested.zip').compress_type==ZIP_STORED
        assert archive.read('nested.zip')==binary


@pytest.mark.parametrize('channel',['shape_sdf_a','aromatic_x'])
def test_contact_site_cannot_intersect_different_four_repeat_sets(monkeypatch,channel):
    from ecd_epitope_scanner.surface_perturbation import engine
    from ecd_epitope_scanner.surface_perturbation.classification import classify_fields
    from ecd_epitope_scanner.surface_perturbation.null_store import FIELD_CHANNELS
    pairs=[{'run':i} for i in range(5)]
    runs=[{name:np.array([.5 if name==channel and i in (0,1,2,4) else 0.]) for name in FIELD_CHANNELS} for i in range(5)]
    calibration={name:calibrate_channel(np.zeros((4,1)),np.array([.1]),'matched',name) for name in FIELD_CHANNELS}
    monkeypatch.setattr(engine,'direct_support',lambda *args:np.array([True]))
    monkeypatch.setattr(engine,'direct_site_masks',lambda pair,*args:[np.array([pair['run'] in (0,1,2,3)]),np.array([pair['run'] in (1,2,3,4)])])
    ref={'area':np.array([1.]),'graph':csr_matrix((1,1)),'owner':np.array([1])}
    result=classify_fields(runs,pairs,{},[],ref,calibration,'matched')
    assert result['direct'].all() and not result['hard'].any() and not result['two_direct_sites'].any()


def test_duplicate_coordinates_share_buffers_without_deduplicating_real_runs():
    from ecd_epitope_scanner.surface_perturbation.buffers import repeated_field
    calls=[];buffer={};evidence=[]
    def compute():calls.append(1);return {'shape_sdf_a':np.array([.2])}
    for run in range(5):
        pair={'run':run,'wt_sha256':'wt','mut_sha256':'mut','chain':'A'}
        value,reused=repeated_field(pair,buffer,compute);evidence.append(value)
        assert reused==(run>0)
    assert len(evidence)==5 and len(calls)==1 and all(value is evidence[0] for value in evidence)
    repeated_field({'wt_sha256':'wt','mut_sha256':'changed','chain':'A'},buffer,compute)
    repeated_field({'wt_sha256':'wt','mut_sha256':'mut','chain':'B'},buffer,compute)
    assert len(calls)==3


def test_sdf_inside_outside_and_open_unknown():
    mesh=cube()
    distance,valid=mesh.signed_distance([[.4,.3,.2],[2,.3,.2],[0,.3,.2]])
    assert valid.all()
    np.testing.assert_allclose(distance,[-.2,1,0],atol=1e-8)
    missing=Mesh(mesh.vertices,mesh.faces[:-1])
    distance,valid=missing.signed_distance([[.4,.3,.2]])
    assert np.isnan(distance).all() and not valid.any()


def test_pair_transform_preserves_mutation_delta():
    points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1]],float)
    reference={i:p for i,p in enumerate(points)}
    wt={i:p+[3,4,5] for i,p in enumerate(points)}
    mutant={i:p+[3,4,5]+([0,0,.3] if i==1 else [0,0,0]) for i,p in enumerate(points)}
    rotation,translation,error=align_pair(reference,wt,mutant,reference.keys())
    assert error < 1e-12
    np.testing.assert_allclose(mutant[1]@rotation+translation-reference[1],[0,0,.3],atol=1e-12)


def test_calibration_mismatch_missing_and_four_of_five():
    calibration=calibrate_channel(np.zeros((4,3)),np.array([.1,.1,np.nan]),"p","shape")
    delta=np.array([[.2,.2,.2],[.2,.2,.2],[.2,.2,.2],[.2,-.2,.2],[0,-.2,.2]])
    changed,unknown=resolved_changes(delta,calibration,"p")
    assert changed.tolist()==[True,False,False]
    assert unknown.tolist()==[False,False,True]
    changed,unknown=resolved_changes(delta,calibration,"another")
    assert not changed.any() and unknown.all()
    changed,unknown=resolved_changes(delta[:1],calibration,"p")
    assert not changed.any() and unknown.all()


def test_calibration_requires_actual_requested_local_area():
    from ecd_epitope_scanner.surface_perturbation.calibration_scope import local_scope
    blocks=np.array([[0.,np.nan,np.nan]]*4)
    value=calibrate_channel(blocks,np.zeros(3),'p','shape')
    assert not value['requested_scope_converged']
    scoped=calibrate_channel(blocks,np.zeros(3),'p','shape',
        required_mask=np.array([True,False,False]),area=np.array([2.,3.,5.]))
    assert scoped['requested_scope_converged']
    assert scoped['calibrated_local_area_fraction']==1.
    outside=local_scope({'shape':scoped},np.array([False,True,True]),np.array([2.,3.,5.]))
    assert not outside['complete']
    assert outside['calibrated_local_area_fraction']==0.


def test_identity_stability_cannot_waive_uncalibrated_blindspot_candidates():
    from ecd_epitope_scanner.surface_perturbation.convergence import convergence_state
    assert convergence_state(['verified'],['verified'],['unknown-positive-proxy'])=='CALIBRATION_EXPANSION_REQUIRED'
    assert convergence_state(['verified'],['verified'],['unknown-zero-proxy'])=='CALIBRATION_EXPANSION_REQUIRED'
    assert convergence_state([],None,['unknown'])=='CALIBRATION_EXPANSION_REQUIRED'
    assert convergence_state(['verified'],['verified'],[])=='VERIFIED_PANEL'


def test_hard_upper_requires_actual_five_models_and_dominates_every_allowed_seed():
    from ecd_epitope_scanner.surface_perturbation.coverage_upper import hard_scan_upper_bound
    from ecd_epitope_scanner.surface_perturbation.coverage import scan_support
    graph=csr_matrix([[0.,13.,0.],[13.,0.,13.],[0.,13.,0.]])
    support=np.array([[1,1,0],[1,1,0],[1,1,1],[1,0,1],[0,0,1]],bool)
    upper=hard_scan_upper_bound(support,graph)
    assert upper.tolist()==[True,False,False]
    allowed_seed=np.sum(support,axis=0)>=4
    lower=scan_support(graph,allowed_seed,np.zeros(3,int))
    assert not np.any(lower&~upper)
    with pytest.raises(ValueError,match='five actual'):hard_scan_upper_bound(support[:1],graph)


def test_legacy_raw_reuse_requires_identical_physical_source_and_keeps_history(tmp_path,monkeypatch):
    import hashlib,shutil
    from pathlib import Path
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import engine,legacy_fields
    from ecd_epitope_scanner.surface_perturbation.tools import fingerprint
    from ecd_epitope_scanner.util import sha256_file
    source_code=Path(engine.__file__).parent
    snapshot=tmp_path/'snapshot';current=tmp_path/'current'
    snapshot.mkdir();current.mkdir()
    names=('fields.py','geometry.py','chemistry.py','reference.py','calibration.py','tools.py','null_store.py')
    for name in names:
        shutil.copy2(source_code/name,snapshot/name);shutil.copy2(source_code/name,current/name)
    legacy=hashlib.sha256(b''.join((snapshot/n).read_bytes() for n in names)).hexdigest()
    (snapshot/'manifest.json').write_text(json.dumps({'implementation_id':legacy,
        'source_sha256':{n:sha256_file(snapshot/n) for n in names}}))
    with (current/'calibration.py').open('a') as output:output.write('\n# Revised acceptance metadata\n')
    monkeypatch.setattr(legacy_fields,'__file__',str(current/'legacy_fields.py'))
    monkeypatch.setattr(engine,'IMPLEMENTATION_ID','new-implementation')
    monkeypatch.setattr(engine,'reference_id',lambda reference:'same-reference')
    tools=SimpleNamespace(config={'legacy_numerical_source_snapshot':str(snapshot)},nano='nano',identity=lambda nano:'same-nano')
    fields=tmp_path/'fields';fields.mkdir()
    contract={'reference':'same-reference','implementation_sha256':legacy,'nano':'same-nano'}
    path=fields/('field_'+fingerprint(contract)+'.npz');np.savez(path,shape_sdf_a=np.array([.2,np.nan]))
    path.with_suffix('.json').write_text(json.dumps({'contract':contract,'sha256':sha256_file(path)}))
    original=path.read_bytes()
    migrated=engine.seed_raw_field_cache(fields,fields,{},tools)
    assert len(migrated)==1 and migrated[0]['raw_value_migration']
    target=fields/migrated[0]['field']
    assert target.read_bytes()==original and path.read_bytes()==original
    metadata=json.loads(target.with_suffix('.json').read_text(encoding='utf-8'))
    assert metadata['contract']['implementation_sha256']=='new-implementation'
    with (current/'geometry.py').open('a') as output:output.write('\n# Physical source changed\n')
    assert legacy_fields.verified_legacy_id(tools.config) is None
    assert engine.seed_raw_field_cache(fields,tmp_path/'different-physical-code',{},tools)==[
        {'field':target.name,'sha256':sha256_file(target),'source':str(target)}]


def test_real_engine_does_not_finalize_calibrated_subset_with_unresolved_gain(tmp_path,monkeypatch):
    from types import SimpleNamespace as NS
    from ecd_epitope_scanner.surface_perturbation import engine,export,benchmark
    from ecd_epitope_scanner.surface_perturbation.null_store import FIELD_CHANNELS
    from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
    pdb=tmp_path/'native.pdb'
    write_atoms(pdb,{('A',1,'','CA'):(np.array([.5,.5,.5]),'ALA','C')})
    class Tools:
        nano='fixture';calls={'nanoshaper':0,'apbs':0,'pdb2pqr':0,'cache_hits':0}
        def __init__(self,*args):pass
        def ses(self,*args):return cube()
        def identity(self,*args):return 'fixture'
    reference={'points':np.array([[0.,0,0],[1.,0,0],[2.,0,0]]),'normals':np.tile([0.,0,1.],(3,1)),
        'area':np.array([2.,3.,5.]),'owner':np.array([1,2,3]),'designable':np.ones(3,bool),
        'domain':np.array(['ECD']*3),'graph':csr_matrix((3,3)),
        'excluded_area_a2':0.,'internal_cavity_area_a2':0.}
    monkeypatch.setattr(engine,'PhysicalTools',Tools)
    monkeypatch.setattr(engine,'canonical_reference',lambda *args,**kwargs:reference)
    monkeypatch.setattr(engine,'validate_pair',lambda pair:None)
    monkeypatch.setattr(engine,'compute_pair',lambda *args:{n:np.ones(3) for n in FIELD_CHANNELS})
    # This fixture tests the verification state machine. Physical query
    # supplementation has separate geometry/cache regressions and saved-input replay.
    monkeypatch.setattr(engine,'complete_scope',lambda pair,raw,*args:{**raw,
        'query_scope_sha256':np.array('fixture'),'supplemented_query_points':np.array(0)})
    monkeypatch.setattr(engine,'protocol_identity',lambda pairs,*args:pairs[0]['candidate'])
    support=lambda pair,*args:np.array([True,True,False]) if pair['candidate']=='known' else np.array([False,False,True])
    monkeypatch.setattr(engine,'direct_support',support)
    monkeypatch.setattr(engine,'direct_site_masks',lambda pair,*args:[support(pair)])
    calibration={n:calibrate_channel(np.zeros((4,3)),np.full(3,.1),'known',n) for n in FIELD_CHANNELS}
    monkeypatch.setattr(engine,'load_calibration',lambda folder,protocol:calibration if protocol=='known' else {})
    monkeypatch.setattr(benchmark,'controls',lambda *args:{})
    monkeypatch.setattr(export,'evidence_archive',lambda *args:None)
    # This state-machine fixture has synthetic models and upstream evidence.
    # Individual native-evidence validation/export has separate strict regressions.
    from ecd_epitope_scanner.surface_perturbation import local_constructs
    monkeypatch.setattr(local_constructs,'local_construct',lambda candidate,*args:
                        {'candidate_id':candidate.candidate_id})
    monkeypatch.setattr(local_constructs,'export_local_constructs',lambda *args:None)
    record=NS(chain_id='A',auth_seq_id=1,insertion_code='',uniprot_pos=1,
        coordinate_eligible=True,review_reasons=[],domain='ECD')
    candidates=[];rows=[]
    for cid,mutant in [('known','V'),('unknown','L')]:
        option=NS(mutant_aa=mutant,notation='A1'+mutant,residue=NS(one_letter='A',uniprot_pos=1))
        candidates.append(NS(candidate_id=cid,order=1,formal_eligible=True,structural_risk=0.,mutations=[option]))
        rows.append({'mutations':[option.notation],'paired_models':[
            {'candidate':cid,'batch_path':str(tmp_path/cid),'group':1,'run':r,'wt_sha256':cid+'wt','mut_sha256':cid+'mut'} for r in range(5)]})
    final,report=engine.run_physical_preview(structure=pdb,records=[record],options=[],candidates=candidates,
        paired_records=rows,notation=lambda c:[c.mutations[0].notation],supplement=lambda *args:pytest.fail('reuse five real identities'),
        recheck=lambda *args:True,output=tmp_path)
    assert final==[] and report['status']=='CALIBRATION_EXPANSION_REQUIRED'
    assert report['verified_within_calibrated_pool_candidate_ids']==['known']
    assert report['calibration_expansion_candidate_ids']==['unknown']
    requests=json.loads((tmp_path/'surface_perturbation/blindspot_calibration_requests.json').read_text(encoding='utf-8'))
    assert requests['priorities']['unknown']['new_proxy_area_a2']==5.
    assert report['tool_calls']['apbs']==0


def test_unknown_electrostatics_not_zero_and_protected_denominator():
    graph=csr_matrix([[0,1,0],[1,0,1],[0,1,0]],dtype=float)
    hard=np.array([[1,0,0],[0,0,1]],bool)
    area=np.array([1.,8.,1.]); designable=np.array([1,0,1],bool)
    rows=budget_curve(hard,np.full((2,3),np.nan),[0,1],area,designable,graph)
    assert rows[1]["C_hard_all"]==pytest.approx(.2)
    assert rows[1]["C_hard_designable"]==1
    assert rows[1]["C_extended_all"] is None
    assert rows[1]["extended_unknown_area_a2"]==8


def test_distinct_contact_residues_not_duplicate_points():
    graph=csr_matrix([[0,1,0],[1,0,1],[0,1,0]],dtype=float)
    changed=np.array([1,1,0],bool); owner=np.array([2,2,3])
    assert not scan_support(graph,changed,owner,radius=3,minimum_sites=2).any()


def test_maximum_coverage_bounds_and_nested_complementarity():
    graph=csr_matrix([[0,1,0],[1,0,1],[0,1,0]],dtype=float)
    hard=np.array([[1,1,0],[0,0,1]],bool); area=np.array([2.,1.,4.])
    order=nested_panel(hard,area,[0,0],graph,["a","b"])
    assert order==[1,0]
    result=independent_bounds(hard,area,2)
    assert result["lower_area_a2"]==pytest.approx(7)
    assert result["upper_area_a2"]==pytest.approx(7)


def test_nested_tie_prefers_distinct_measured_footprint_before_risk():
    graph=csr_matrix(np.ones((4,4))-np.eye(4))
    # A wins first. B/C add the same cell and leave the same blind patch;
    # C is less confounded with A even though its risk tie-break is worse.
    masks=np.array([[1,1,1,0],[1,1,0,1],[0,0,0,1]],bool)
    assert nested_panel(masks,np.ones(4),[0.,0.,1.],graph,['a','b','c'])==[0,2]


def test_budget_cannot_reset_on_resume(tmp_path):
    path=tmp_path/"budget.json"
    path.write_text(json.dumps({"elapsed_seconds":7201,"limit_seconds":7200,"active":False}))
    budget=RunBudget(path)
    with pytest.raises(TimeoutError): budget.remaining()


def test_exact_duplicate_mesh_import_preserves_area_and_does_not_fill_holes(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.tools import read_off
    mesh=cube();vertices=np.vstack([mesh.vertices,mesh.vertices[0]])
    faces=mesh.faces.copy();faces[faces==0]=8
    faces=np.vstack([faces,[0,8,1]])
    path=tmp_path/'mesh.off'
    path.write_text('OFF+N\n9 13 0\n'+'\n'.join(' '.join(map(str,v)) for v in vertices)+'\n'+
                    '\n'.join('3 '+' '.join(map(str,f)) for f in faces))
    imported=read_off(path)
    assert imported.closed and imported.area.sum()==pytest.approx(6.)
    assert imported.import_audit['zero_area_faces_removed']==1
    assert not Mesh(imported.vertices,imported.faces[:-1]).closed


def test_canonical_area_and_protected_denominator_conserve_native_ses():
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation.reference import canonical_reference
    record=SimpleNamespace(chain_id='A',auth_seq_id=7,insertion_code='B',uniprot_pos=42,domain='ECD')
    ref=canonical_reference(cube(),{('A',7,'B','CA'):(np.array([.5,.5,.5]),'ALA','C')},[record],{42},set(),spacing=.2)
    assert ref['area'].sum()==pytest.approx(6.)
    assert not ref['designable'].any()
    assert (ref['owner']==42).all()


def test_actual_combination_difference_is_not_single_field_union():
    from ecd_epitope_scanner.surface_perturbation.fields import difference
    def field(value):
        return {'sdf':np.array([value]),'shell_sdf':np.array([[value]]),
                'shell_valid':np.array([[True]]),'chemistry':np.zeros((1,14))}
    # A and B individually grow; AB may compensate. The real AB stays unchanged.
    a=difference(field(0),field(-1));b=difference(field(0),field(-2))
    ab=difference(field(0),field(0))
    assert a['shape_sdf_a'][0] and b['shape_sdf_a'][0]
    assert ab['shape_sdf_a'][0]==0 and ab['occupancy_gain_a'][0]==0


def test_preview_and_ordinary_api_accept_twelve():
    from ecd_epitope_scanner.api import JobCreate
    assert JobCreate(identifier='P16410',surface_strategy='surface_perturbation_v1_preview',maximum_constructs=12).maximum_constructs==12
    assert JobCreate(identifier='P16410',maximum_constructs=12).maximum_constructs==12


def test_evidence_zip_manifest_matches_every_exported_file(tmp_path):
    from zipfile import ZipFile
    import hashlib
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    folder=tmp_path/'surface_perturbation';folder.mkdir()
    (folder/'reference_contract.json').write_text('{"basis":"SES"}')
    path=evidence_archive(tmp_path,{},tmp_path/'calibration')
    with ZipFile(path) as archive:
        manifest=json.loads(archive.read('evidence_manifest.json'))
        for name,digest in manifest['files'].items():
            assert hashlib.sha256(archive.read(name)).hexdigest()==digest


def test_pair_grid_exact_spacing_and_insertion_mapping(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
    from ecd_epitope_scanner.surface_perturbation.tools import make_grid
    from ecd_epitope_scanner.surface_perturbation.paired import atoms
    source={('A',1,'A','CA'):(np.array([1.,2.,3.]),'ALA','C'),
            ('A',1,'','CA'):(np.array([5.,6.,7.]),'MET','C')}
    path=tmp_path/'aligned.pdb';write_atoms(path,source,{('A',1,'A'):1,('A',1,''):2})
    assert set(atoms(path))=={('A',1,'','CA'),('A',2,'','CA')}
    grid=make_grid(np.array([[1,2,3],[5,6,7]]),spacing=.4)
    np.testing.assert_allclose(np.array(grid['lengths'])/(np.array(grid['counts'])-1),.4)


def test_single_run_foldx_pairing_uses_explicit_group_and_ddg(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
    from ecd_epitope_scanner.surface_perturbation.paired import scan_pairs,validate_pair
    folder=tmp_path/'engine_runs'/'foldx_batch_test';folder.mkdir(parents=True)
    (folder/'individual_list.txt').write_text('AA1V;\n')
    (folder/'PdbList_input.fxout').write_text('input_1.pdb\nWT_input_1.pdb\n')
    (folder/'Dif_input.fxout').write_text('Pdb\ttotal energy\ninput_1.pdb\t0.35\n')
    source={('A',1,'','CA'):(np.array([0.,0.,0.]),'ALA','C')}
    write_atoms(folder/'WT_input_1.pdb',source);write_atoms(folder/'input.pdb',source)
    write_atoms(folder/'input_1.pdb',{('A',1,'','CA'):(np.array([0.,0.,0.]),'VAL','C')})
    pairs,errors=scan_pairs(tmp_path)
    assert not errors and pairs[0]['group']==1 and pairs[0]['run']==0
    assert pairs[0]['ddg_kcal_mol']==pytest.approx(.35)
    validate_pair(pairs[0])
    with pytest.raises(ValueError,match='metadata'):validate_pair({**pairs[0],'mutations':['A1N']})
    from zipfile import ZipFile
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    with ZipFile(evidence_archive(tmp_path,{'single':pairs},tmp_path/'null')) as archive:
        review=tmp_path/'review';archive.extractall(review)
    portable=json.loads((review/'portable_pairs.json').read_text(encoding='utf-8'))['single'][0]
    for key in ('wt_path','mut_path','pdb_list_path','ddg_file_path','batch_path','input_path','mutant_file_path'):
        portable[key]=str(review/portable[key])
    validate_pair(portable)
    original_input=(folder/'input.pdb').read_bytes()
    (folder/'input.pdb').write_text('changed')
    with pytest.raises(ValueError,match='input fingerprint'): validate_pair(pairs[0])
    (folder/'input.pdb').write_bytes(original_input)
    (folder/'input_1.pdb').write_text('changed')
    with pytest.raises(ValueError): validate_pair(pairs[0])


def test_sdf_rigid_rotation_and_concave_notch():
    mesh=cube();angle=.63
    rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    queries=np.array([[.3,.4,.5],[1.2,.5,.6]])
    moved=Mesh(mesh.vertices@rotation+[2,4,6],mesh.faces)
    original,valid=mesh.signed_distance(queries)
    transformed,valid_moved=moved.signed_distance(queries@rotation+[2,4,6])
    assert valid.all() and valid_moved.all()
    np.testing.assert_allclose(original,transformed,atol=1e-8)
    # Two closed cubes leave a real groove: no convex hull fills the gap.
    groove=Mesh(np.vstack([mesh.vertices,mesh.vertices+[2,0,0]]),np.vstack([mesh.faces,mesh.faces+8]))
    distance,valid=groove.signed_distance([[1.5,.4,.4]])
    assert valid[0] and distance[0]==pytest.approx(.5)


def test_missing_null_cannot_create_final_construct_or_call_apbs(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import engine
    from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
    source={('A',1,'','CA'):(np.array([.5,.5,.5]),'ALA','C')}
    pdb=tmp_path/'native.pdb';write_atoms(pdb,source)
    class Tools:
        calls={'nanoshaper':0,'apbs':0,'pdb2pqr':0,'cache_hits':0}
        def __init__(self,*args): pass
        def ses(self,*args): return cube()
    monkeypatch.setattr(engine,'PhysicalTools',Tools)
    record=SimpleNamespace(chain_id='A',auth_seq_id=1,insertion_code='',uniprot_pos=1,
                           coordinate_eligible=True,review_reasons=[],domain='ECD')
    final,report=engine.run_physical_preview(structure=pdb,records=[record],options=[],candidates=[],
        paired_records=[],notation=lambda c:[],supplement=lambda c,n:pytest.fail('no eligible pool'),
        recheck=lambda c,p:False,output=tmp_path)
    assert final==[] and report['status']=='NO_VERIFIED_PANEL'
    assert report['tool_calls']['apbs']==0


def test_concurrent_retry_cannot_double_use_same_budget(tmp_path):
    path=tmp_path/'budget.json';budget=RunBudget(path)
    with pytest.raises(RuntimeError,match='CASE_ALREADY_RUNNING'): RunBudget(path)
    budget.save(False)
    resumed=RunBudget(path);assert resumed.carried>=0
    resumed.save(False)


@pytest.mark.parametrize('profile,expected',[(None,[0]),('coverage_scan_v2_preview',[0,1,2,3,4])])
def test_one_run_shortlist_does_not_calculate_twenty_cached_fields(tmp_path, monkeypatch,profile,expected):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import engine, export
    from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
    source={('A',1,'','CA'):(np.array([.5,.5,.5]),'ALA','C')}
    pdb=tmp_path/'native.pdb';write_atoms(pdb,source)
    class Tools:
        calls={'nanoshaper':0,'apbs':0,'pdb2pqr':0,'cache_hits':0}
        def __init__(self,*args): pass
        def ses(self,*args): return cube()
    monkeypatch.setattr(engine,'PhysicalTools',Tools)
    monkeypatch.setattr(engine,'validate_pair',lambda row:None)
    calculated=[]
    monkeypatch.setattr(engine,'compute_pair',lambda row,*args:calculated.append(row['run']) or {})
    monkeypatch.setattr(export,'evidence_archive',lambda *args:None)
    record=SimpleNamespace(chain_id='A',auth_seq_id=1,insertion_code='',uniprot_pos=1,
                           coordinate_eligible=True,review_reasons=[],domain='ECD')
    option=SimpleNamespace(mutant_aa='V',residue=SimpleNamespace(one_letter='A'))
    candidate=SimpleNamespace(candidate_id='cached',order=2,formal_eligible=True,structural_risk=0.,mutations=[option])
    stale=SimpleNamespace(candidate_id='forbidden-cached',order=1,formal_eligible=True,
        mutations=[SimpleNamespace(mutant_aa='C',residue=SimpleNamespace(one_letter='A'))])
    rows=[{'batch_path':str(tmp_path),'group':1,'run':i} for i in range(20)]
    final,report=engine.run_physical_preview(structure=pdb,records=[record],options=[],candidates=[candidate,stale],
        paired_records=[{'mutations':['A1V'],'paired_models':rows}],notation=lambda c:['A1V'],
        supplement=lambda c,n:pytest.fail('cached repeats must be reused'),recheck=lambda c,p:False,output=tmp_path,
        settings={'optimization_profile':profile})
    assert calculated==expected
    assert final==[] and report['tool_calls']['apbs']==0
    assert len(json.loads((tmp_path/'surface_perturbation'/'paired_manifest.json').read_text(encoding='utf-8'))['cached'])==20
    assert 'forbidden-cached' not in json.loads((tmp_path/'surface_perturbation'/'paired_manifest.json').read_text(encoding='utf-8'))


def test_portable_archive_retains_input_and_mutation_records(tmp_path):
    from zipfile import ZipFile
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    batch=tmp_path/'batch';batch.mkdir()
    for name in ('wt.pdb','mut.pdb','PdbList.fxout','Dif.fxout','input.pdb','individual_list.txt'):
        (batch/name).write_text(name,encoding='utf-8')
    pair={key:str(batch/name) for key,name in [('wt_path','wt.pdb'),('mut_path','mut.pdb'),
        ('pdb_list_path','PdbList.fxout'),('ddg_file_path','Dif.fxout')]}
    pair['batch_path']=str(batch)
    (tmp_path/'structure_for_energy.pdb').write_bytes(b'native author mapping')
    (tmp_path/'source_confidence.json').write_bytes(b'{"original":true}')
    failed=tmp_path/'null/task_failures/protocol.json'
    failed.parent.mkdir(parents=True)
    failed.write_text('{"status":"CALIBRATION_TASK_FAILED"}')
    with ZipFile(evidence_archive(tmp_path,{'candidate':[pair]},tmp_path/'null')) as archive:
        portable=json.loads(archive.read('portable_pairs.json'))['candidate'][0]
        assert archive.read(portable['input_path'])==b'input.pdb'
        assert archive.read(portable['mutant_file_path'])==b'individual_list.txt'
        assert portable['batch_path']=='models/candidate/000'
        assert archive.read('source/structure_for_energy.pdb')==b'native author mapping'
        assert archive.read('source/source_confidence.json')==b'{"original":true}'
        assert archive.read('calibration/task_failures/protocol.json')==failed.read_bytes()


def test_spurs_auto_cuda_cache_reuses_only_exact_weights(tmp_path,monkeypatch):
    import hashlib
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import upstream_cache as cache
    from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
    repo=tmp_path/'repo'; module=repo/'surface_perturbation'/'upstream_cache.py'
    module.parent.mkdir(parents=True);module.write_text('')
    monkeypatch.setattr(cache,'__file__',str(module))
    weights={}
    for relative in ('models/spurs/spurs/checkpoints/best.ckpt','models/spurs/esm2_t33_650M_UR50D.pt'):
        path=repo/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(b'weights')
        weights[relative]=hashlib.sha256(b'weights').hexdigest()
    for name in ('spurs_gpu_worker.py','gpu_runtime.py'):
        path=repo/'integrations'/'workers'/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('worker')
    monkeypatch.setattr(cache.subprocess,'run',lambda argv,**kwargs:SimpleNamespace(
        stdout=b'{"python":"test-python","torch":"test-torch"}' if '-c' in argv else b'test-driver, test-gpu'))
    model=tmp_path/'native.pdb';model.write_text('model')
    request=SimpleNamespace(structure_path=str(model),chain_id='A',sequence='A',structure_type='pdb',
                            mutations=[],parameters={},device='auto')
    calls=[]
    class Adapter:
        cache_dir=tmp_path/'cache'
        def predict(self,r):
            calls.append(r)
            return PredictionResult(status='success',records=[{'position':1}],provenance={
                'actual_device':'cuda','model_files_sha256':weights})
        def checked_result(self,r,v): return v
    adapter=Adapter();cache.predict_cached(adapter,request)
    assert cache.predict_cached(adapter,request).provenance['cache_hit']
    assert len(calls)==1
    frozen=tmp_path/'frozen_assets.json';frozen.write_text(json.dumps({'files':weights}))
    request.parameters['frozen_model_assets_path']=str(frozen)
    assert cache.predict_cached(adapter,request).provenance['cache_hit']
    assert len(calls)==1 # Adding a qualification assertion does not change raw-score identity.
    (repo/'models/spurs/spurs/checkpoints/best.ckpt').write_bytes(b'changed')
    with pytest.raises(ValueError,match='identity changed'):cache.predict_cached(adapter,request)
    assert len(calls)==1
    request.parameters.pop('frozen_model_assets_path')
    cache.predict_cached(adapter,request)
    assert len(calls)==2


def test_electrostatic_medium_switch_is_unknown_at_every_offset(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import electrostatics as es
    from ecd_epitope_scanner.surface_perturbation.paired import atoms
    wt={('A',1,'','CA'):(np.array([.5,.3,.3]),'ALA','C')}
    mut={('A',1,'','CA'):(np.array([1.5,.3,.3]),'VAL','C')}
    monkeypatch.setattr(es,'aligned_atoms',lambda *args:(wt,mut,{'mutated_author_residues':[('A',1,'')]}))
    class Tools:
        def prepare_pqr(self,path):
            values={(key[0],key[1],key[3]):{'xyz':v[0].tolist(),'residue':v[1],'charge':0.,'radius':1.7}
                    for key,v in atoms(path).items()}
            return path,values
        def potential(self,path,grid): return lambda points:np.ones(len(points))*('Mut' in path.name)
        def ses(self,source):
            mesh=cube()
            return mesh if source is wt else Mesh(mesh.vertices+[1,0,0],mesh.faces)
    reference={'points':np.array([[.5,.3,.3]]),'normals':np.array([[1.,0.,0.]])}
    monkeypatch.setattr(es,'pqr_surface',lambda tools,values:cube() if next(iter(values.values()))['xyz'][0]<1 else Mesh(cube().vertices+[1,0,0],cube().faces))
    changed,unknown,report=es.audit_candidate('switch',[{}]*5,wt,[],reference,Tools(),tmp_path)
    assert not changed.any() and unknown.all()
    fields=np.load(tmp_path/'switch_electrostatics.npz')
    assert np.isnan(fields['offset_deltas'][:,:2]).all()
    assert not fields['offset_solvent_valid'][:,:2].any()
    assert fields['offset_solvent_valid'][:,2].all()
    assert report['unknown_common_medium_queries']==5
    assert fields['null_blocks'].shape==(4,1)


def test_independent_route_search_preserves_cross_boundary_candidates():
    from ecd_epitope_scanner.surface_perturbation.routes import search_controls
    from ecd_epitope_scanner.joint_design_v4.models import ResidueKey,V4MutationOption,V4Sector,V4SurfaceHypothesis
    options=[V4MutationOption(residue=ResidueKey('0','A',p,uniprot_pos=p,one_letter='D'),
        mutant_aa='A',notation=f'D{p}A',chemistry_class='charge_removal',spurs_ddg_kcal_mol=0.,
        score=1.,eligible=True) for p in (1,2)]
    sectors=[V4Sector('left',[1],0,1.,1.,0.),V4Sector('right',[2],0,1.,1.,0.)]
    graph={1:{2:3.},2:{1:3.}}
    report=search_controls(sectors,options,graph,{'orders':[2],'exploratory_orders':[]},
        [V4SurfaceHypothesis('overlap',{1:1.,2:1.})])
    assert report['hard_partition']['unique_candidates']==0
    assert report['unpartitioned_local']['unique_candidates']==1
    assert report['overlapping_windows']['unique_candidates']==1
    assert report['unpartitioned_local']['unmeasured_candidates']==1


def test_posthoc_contact_labels_cannot_fill_unknown_reference_or_change_panel(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.sanity import contact_sanity
    from ecd_epitope_scanner.util import sha256_file
    path=tmp_path/'labels.json'
    path.write_text(json.dumps({'numbering':'UniProt precursor numbering','accession':'test',
        'contacts':{'hidden':{'uniprot_positions':[1,3,99]}}}),encoding='utf-8')
    hard=np.array([[1,0,0]],bool);before=hard.copy()
    result=contact_sanity(path,sha256_file(path),{'area':np.array([1.,3.,0.]),'owner':np.array([1,2,3])},hard,['probe'])
    assert result['contacts']['hidden']['curve'][0]['contact_center_overlap']==1.
    assert result['contacts']['hidden']['positions_outside_reliable_target']==[3,99]
    np.testing.assert_array_equal(hard,before)
    with pytest.raises(ValueError,match='fingerprint'):contact_sanity(path,'wrong',{},hard,[])


def test_independent_dense_reference_audit_conserves_area_and_qualifications():
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation.reference_audit import sampling_audit
    record=SimpleNamespace(chain_id='A',auth_seq_id=7,insertion_code='B',uniprot_pos=42,domain='ECD')
    source={('A',7,'B','CA'):(np.array([.5,.5,.5]),'ALA','C')}
    result=sampling_audit(cube(),cube(),source,[record],{42},set())
    assert all(r['target_area_a2']==pytest.approx(6.) for r in result['settings'])
    assert all(r['designable_area_a2']==0. for r in result['settings'])
    assert result['relative_target_area_resampling_difference']==0.


def test_apbs_medium_surface_cache_tracks_actual_pqr_radii(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.pqr_surface import pqr_surface
    class Tools:
        cache=tmp_path;nano='fixture-nano';calls={'nanoshaper':0,'cache_hits':0}
        def identity(self,path):return 'a'*64
        def run(self,argv,folder):
            mesh=cube()
            (folder/'triangulatedSurf.off').write_text('OFF\n8 12 0\n'+
                '\n'.join(' '.join(map(str,p)) for p in mesh.vertices)+'\n'+
                '\n'.join('3 '+' '.join(map(str,f)) for f in mesh.faces),encoding='ascii')
            return 'fixture mesher; test cache identity only'
    params={('A',1,'C'):{'xyz':[0.,0.,0.],'radius':1.7},
            ('A',1,'H'):{'xyz':[.1,.1,.1],'radius':0.}}
    tools=Tools();assert pqr_surface(tools,params).closed
    assert pqr_surface(tools,params).closed
    assert tools.calls=={'nanoshaper':1,'cache_hits':1}
    params[('A',1,'C')]['radius']=1.8
    pqr_surface(tools,params)
    assert tools.calls['nanoshaper']==2
    params[('A',1,'C')]['radius']=-1.
    with pytest.raises(ValueError,match='radius'):pqr_surface(tools,params)


def test_residue_name_change_does_not_make_unchanged_carbon_direct_support(monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import engine
    wt={('A',1,'','CB'):(np.array([0.,0.,0.]),'PHE','C')}
    mut={('A',1,'','CB'):(np.array([0.,0.,0.]),'TYR','C'),
         ('A',1,'','OH'):(np.array([10.,0.,0.]),'TYR','O')}
    monkeypatch.setattr(engine,'aligned_atoms',lambda *args:(wt,mut,{'mutated_author_residues':[('A',1,'')],'core_rmsd_a':0.}))
    ref={'points':np.array([[.1,0,0],[10.1,0,0]]),'area':np.ones(2)}
    assert engine.direct_support({}, {}, [],ref).tolist()==[False,True]


def test_vector_noise_uses_same_four_amplitudes_and_directions_and_rotates():
    from ecd_epitope_scanner.surface_perturbation.directions import resolve_vector
    names=('x','y','z')
    calibration={name:calibrate_channel(np.zeros((4,1)),np.zeros(1),'p',name) for name in names}
    vectors=np.array([[1.,10.,0.]]*3+[[1.,-10.,0.]]*2)[:,None,:]
    hit,_,_,_=resolve_vector(vectors,calibration,names,'p')
    assert not hit.any() # Four positive x changes do not establish one repeated vector direction.
    vectors[3]=vectors[0]
    hit,_,unknown,_=resolve_vector(vectors,calibration,names,'p')
    assert hit.all() and not unknown.any()
    rotation=np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]])
    rotated=resolve_vector(vectors@rotation,calibration,names,'p')[0]
    np.testing.assert_array_equal(hit,rotated)


def test_local_vector_qualification_is_rigid_rotation_invariant():
    from ecd_epitope_scanner.surface_perturbation.calibration_scope import local_scope
    from ecd_epitope_scanner.surface_perturbation.directions import resolve_vector
    names=('donor_x','donor_y','donor_z')
    null=np.array([[0.,1.,0.]]*3+[[1.,0.,0.]]*3)[:,None,:]
    rotation=np.array([[2**-.5,-2**-.5,0.],[2**-.5,2**-.5,0.],[0.,0.,1.]])
    for transform in (np.eye(3),rotation):
        transformed=null@transform
        calibration={name:calibrate_channel(transformed[:,:,i],np.zeros(1),'p',name)
                     for i,name in enumerate(names)}
        scope=local_scope(calibration,np.ones(1,bool),np.ones(1))
        assert scope['complete']
        assert set(scope['channels'])=={'donor_delta'}
        assert set(scope['component_diagnostics'])==set(names)
        hit,_,unknown,_=resolve_vector(np.tile([[[2.,0.,0.]]],(5,1,1))@transform,
                                      calibration,names,'p')
        assert hit.all() and not unknown.any()
    # The old separate x-component gate falsely rejected the first frame.
    assert not calibrate_channel(null[:,:,0],np.zeros(1),'p',names[0])['requested_scope_converged']


def test_joint_scope_preserves_unknown_and_unstable_vector_noise():
    from ecd_epitope_scanner.surface_perturbation.calibration_scope import local_scope
    names=('donor_x','donor_y','donor_z')
    # Each axis marginal is stable, but their joint amplitude distribution changes.
    null=np.array([[1.,0.,0.],[0.,1.,0.]]*2+
                  [[1.,1.,0.]]*2+[[0.,0.,0.]]*2)[:,None,:]
    calibration={name:calibrate_channel(null[:,:,i],np.zeros(1),'p',name)
                 for i,name in enumerate(names)}
    assert all(v['requested_scope_converged'] for v in calibration.values())
    scope=local_scope(calibration,np.ones(1,bool),np.ones(1))
    assert not scope['complete']
    assert scope['calibrated_local_area_fraction']==0.
    calibration.pop('donor_z')
    assert not local_scope(calibration,np.ones(1,bool),np.ones(1))['complete']


def test_direct_support_and_perturbation_must_share_four_repeats():
    from ecd_epitope_scanner.surface_perturbation.engine import supported_run_values
    calibration=calibrate_channel(np.zeros((4,1)),np.zeros(1),'p','shape')
    values=np.array([[2.],[2.],[2.],[2.],[0.]])
    support=np.array([[1],[1],[1],[0],[1]],bool)
    assert not resolved_changes(supported_run_values(values,support),calibration,'p')[0].any()
    support[3]=True
    assert resolved_changes(supported_run_values(values,support),calibration,'p')[0].all()
    np.testing.assert_array_equal(values[:,0],[2.,2.,2.,2.,0.])


def test_unconverged_apbs_grid_or_box_cannot_supply_a_finite_error_bound():
    from ecd_epitope_scanner.surface_perturbation.electrostatics import grid_resolution
    base=np.ones(4)
    variants=np.array([[1.02,1.2,1.01,1.04],[1.04,1.1,1.02,1.02],[1.01,1.01,1.3,1.01]])
    bound,valid,_=grid_resolution(base,variants)
    assert valid.tolist()==[True,False,False,False]
    assert bound[0]==pytest.approx(.04) and np.isnan(bound[1:]).all()


def test_numerical_sensitivity_cache_requires_the_exact_first_pair(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import engine
    calls=[]
    monkeypatch.setattr(engine,'protocol_identity',lambda *args:'matched')
    monkeypatch.setattr(engine,'numerical_bounds',lambda *args:calls.append(1) or
        {name:np.array([.1]) for name in engine.FIELD_CHANNELS})
    pair={'wt_sha256':'wt','mut_sha256':'mut'}
    args=(pair,{},[],{'area':np.ones(1)},None,{},tmp_path)
    engine.cached_numerical_bounds(*args);engine.cached_numerical_bounds(*args)
    assert len(calls)==1
    pair['mut_sha256']='different';engine.cached_numerical_bounds(*args)
    assert len(calls)==2


def test_canonical_queries_and_all_spatial_fields_follow_one_rigid_transform():
    from ecd_epitope_scanner.surface_perturbation.fields import spatial_field
    from ecd_epitope_scanner.surface_perturbation.chemistry import CHANNELS
    from ecd_epitope_scanner.surface_perturbation.directions import VECTOR_GROUPS
    source={('A',1,'',name):(np.array(xyz),'ALA',element) for name,xyz,element in
        [('N',[.2,.3,.4],'N'),('CA',[.5,.5,.5],'C'),('C',[.7,.6,.5],'C'),('O',[.8,.7,.6],'O')]}
    reference={'points':np.array([[1.,.23,.37],[.27,1.,.73]]),
               'normals':np.array([[1.,0,0],[0,1.,0]])}
    class Tools:
        def __init__(self,mesh):self.mesh=mesh
        def ses(self,*args,**kwargs):return self.mesh
    mesh=cube();original=spatial_field(source,reference,Tools(mesh))
    angle=.51;rotation=np.array([[np.cos(angle),-np.sin(angle),0],
                                 [np.sin(angle),np.cos(angle),0],[0,0,1.]])
    translation=np.array([3.,-7.,2.])
    moved_source={k:(v[0]@rotation+translation,v[1],v[2]) for k,v in source.items()}
    moved_reference={'points':reference['points']@rotation+translation,'normals':reference['normals']@rotation}
    moved=spatial_field(moved_source,moved_reference,Tools(Mesh(mesh.vertices@rotation+translation,mesh.faces)))
    for name in ('sdf','shell_sdf','shell_valid','shell_owned','atomic_inside','local_ses_area'):
        np.testing.assert_allclose(original[name],moved[name],atol=1e-8,equal_nan=True)
    np.testing.assert_allclose(original['projected_normals']@rotation,moved['projected_normals'],atol=1e-8)
    directional={n for names in VECTOR_GROUPS.values() for n in names}
    for i,name in enumerate(CHANNELS):
        if name not in directional:np.testing.assert_allclose(original['chemistry'][:,i],moved['chemistry'][:,i],atol=1e-8)
    for names in VECTOR_GROUPS.values():
        if not all(n in CHANNELS for n in names):continue
        indices=[CHANNELS.index(n) for n in names]
        np.testing.assert_allclose(original['chemistry'][:,indices]@rotation,moved['chemistry'][:,indices],atol=1e-8)


def test_hidden_contacts_use_label_alignment_and_exclude_only_declared_conflicts():
    from ecd_epitope_scanner.surface_perturbation.contact_labels import label_mapping
    data={'_struct_ref.id':['1'],'_struct_ref.db_name':['UNP'],'_struct_ref.pdbx_db_accession':['target'],
          '_struct_ref_seq.ref_id':['1'],'_struct_ref_seq.pdbx_strand_id':['A'],
          '_struct_ref_seq.seq_align_beg':['1'],'_struct_ref_seq.seq_align_end':['2'],
          '_struct_ref_seq.db_align_beg':['7'],'_struct_ref_seq.db_align_end':['8'],
          '_atom_site.auth_asym_id':['A','A'],'_atom_site.label_seq_id':['1','2'],
          '_atom_site.label_comp_id':['ALA','GLY']}
    assert label_mapping(data,'target','A','AAAAAAAG')=={1:7,2:8}
    data['_atom_site.label_comp_id'][0]='VAL'
    with pytest.raises(ValueError,match='mismatch'):label_mapping(data,'target','A','AAAAAAAG')
    data.update({'_struct_ref_seq_dif.pdbx_pdb_strand_id':['A'],'_struct_ref_seq_dif.seq_num':['1'],
        '_struct_ref_seq_dif.pdbx_seq_db_accession_code':['target'],'_struct_ref_seq_dif.mon_id':['VAL'],
        '_struct_ref_seq_dif.db_mon_id':['ALA'],'_struct_ref_seq_dif.pdbx_seq_db_seq_num':['7']})
    excluded=[]
    assert label_mapping(data,'target','A','AAAAAAAG',declared_conflicts=excluded)=={2:8}
    assert excluded[0]['uniprot_pos']==7


def test_physical_route_outputs_share_one_case_runtime_ledger(tmp_path,monkeypatch):
    from ecd_epitope_scanner import pipeline
    from ecd_epitope_scanner.surface_perturbation.runtime import ACTIVE_BUDGET
    ledger=tmp_path/'shared_case_budget.json'
    ledger.write_text(json.dumps({'elapsed_seconds':7100.,'limit_seconds':7200.,'active':False}))
    monkeypatch.setattr(pipeline,'_run_pipeline',lambda *args:{'path':str(ACTIVE_BUDGET.get().path)})
    elapsed=[]
    for route in ('hard_partition','unpartitioned_local'):
        request=pipeline.PipelineRequest(input_type='structure',surface_strategy='surface_perturbation_v1_preview',
            output_dir=str(tmp_path/route),metadata={'joint_v4_ml':{'surface_perturbation':{
                'candidate_route':route,'cumulative_budget_path':str(ledger)}}})
        assert pipeline.run_pipeline(request)['path']==str(ledger)
        saved=json.loads(ledger.read_text());elapsed.append(saved['elapsed_seconds'])
        assert saved['limit_seconds']==7200. and not saved['active']
    assert elapsed[1]>=elapsed[0]>=7100.


def test_cross_phase_raw_field_reuse_checks_content_and_implementation(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import engine
    from ecd_epitope_scanner.surface_perturbation.tools import fingerprint
    from ecd_epitope_scanner.util import sha256_file
    from types import SimpleNamespace
    monkeypatch.setattr(engine,'reference_id',lambda ref:'native-reference')
    source=tmp_path/'calibration';source.mkdir()
    contract={'reference':'native-reference','implementation_sha256':engine.IMPLEMENTATION_ID,'nano':'nano-sha'}
    path=source/('field_'+fingerprint(contract)+'.npz');np.savez_compressed(path,shape=np.array([1.]))
    path.with_suffix('.json').write_text(json.dumps({'contract':contract,'sha256':sha256_file(path)}))
    tools=SimpleNamespace(nano='NanoShaper',identity=lambda binary:'nano-sha')
    destination=tmp_path/'production'
    result=engine.seed_raw_field_cache(source,destination,{},tools)
    assert len(result)==1 and (destination/path.name).read_bytes()==path.read_bytes()
    monkeypatch.setattr(engine,'IMPLEMENTATION_ID','changed-scientific-code')
    assert engine.seed_raw_field_cache(source,tmp_path/'different',{},tools)==[]
    monkeypatch.setattr(engine,'IMPLEMENTATION_ID',contract['implementation_sha256'])
    path.write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='fingerprint'):engine.seed_raw_field_cache(source,tmp_path/'corrupt',{},tools)


def test_three_site_proxy_gain_cannot_replace_real_physical_gain_or_stability():
    from types import SimpleNamespace
    from ecd_epitope_scanner.joint_design_v4.engine import _formal_candidate_reasons
    candidate=SimpleNamespace(qualification='FORMAL_POOL',mutations=[],connected=True,order=3,
        coverage_by_scenario={},stability=SimpleNamespace(consensus_status='CONSENSUS_PASS',spurs_multi_ddg_kcal_mol=.1,foldx_ddg_kcal_mol=.2))
    kwargs={'allow_spurs_only':False,'minimum_coverage_gain':.05,'lower_order_coverage':0.}
    assert 'INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN' in _formal_candidate_reasons(candidate,**kwargs)
    assert _formal_candidate_reasons(candidate,physical_coverage_pending=True,**kwargs)==[]
    candidate.stability.consensus_status='CONSENSUS_FAIL'
    assert 'CONSENSUS_FAIL' in _formal_candidate_reasons(candidate,physical_coverage_pending=True,**kwargs)
    candidate.stability.foldx_ddg_kcal_mol=None
    assert 'HIGHER_ORDER_REQUIRES_FOLDX' in _formal_candidate_reasons(candidate,physical_coverage_pending=True,**kwargs)


def test_stable_enhancement_triple_passes_geometry_gate_but_still_requires_foldx():
    from types import SimpleNamespace
    from ecd_epitope_scanner.joint_design_v4.engine import _formal_candidate_reasons
    candidate=SimpleNamespace(qualification='FORMAL_POOL',mutations=[],connected=True,order=3,
        coverage_by_scenario={},risk_components={'triple_gate_passed':True,'triple_enhancement_delta':.2},
        impact_evidence={},
        stability=SimpleNamespace(consensus_status='CONSENSUS_PASS',spurs_multi_ddg_kcal_mol=.1,foldx_ddg_kcal_mol=.2))
    kwargs={'allow_spurs_only':False,'minimum_coverage_gain':.05,'lower_order_coverage':0.}
    assert _formal_candidate_reasons(candidate,**kwargs)==[]
    candidate.stability.foldx_ddg_kcal_mol=None
    assert 'HIGHER_ORDER_REQUIRES_FOLDX' in _formal_candidate_reasons(candidate,**kwargs)


def test_blindspot_expansion_targets_connected_area_without_discarding_zero_proxy():
    from ecd_epitope_scanner.surface_perturbation.expansion import blindspot_priority
    from types import SimpleNamespace
    candidates=[SimpleNamespace(candidate_id=name,structural_risk=0.) for name in ('small','large','zero')]
    reference={'area':np.array([10.,1.,1.]),'graph':csr_matrix([[0,0,0],[0,0,1],[0,1,0]])}
    proxies={'small':np.array([0,1,1],bool),'large':np.array([1,0,0],bool),'zero':np.zeros(3,bool)}
    order,audit=blindspot_priority(reference,np.zeros(3,bool),candidates,proxies)
    assert [c.candidate_id for c in order]==['large','small','zero']
    assert audit['large']['largest_blind_proxy_area_a2']==10.
    assert audit['small']['largest_blind_proxy_area_a2']==0.


def test_pqr_preserves_native_sequence_when_author_insertion_sorts_after_next_residue(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import electrostatics as es
    from ecd_epitope_scanner.surface_perturbation.paired import atoms
    native={('A',1,'A','CA'):(np.array([0.,0.,0.]),'ALA','C'),
            ('A',1,'','CA'):(np.array([3.,0.,0.]),'MET','C')}
    # Model dictionaries may have another insertion order; native chronology is authoritative.
    modeled=dict(reversed(list(native.items())))
    monkeypatch.setattr(es,'aligned_atoms',lambda *args:(modeled,modeled,{'mutated_author_residues':[]}))
    class Tools:
        def prepare_pqr(self,path):
            values={(k[0],k[1],k[3]):{'xyz':v[0].tolist(),'residue':v[1],'charge':0.,'radius':1.7}
                    for k,v in atoms(path).items()}
            return path,values
    prepared,audits,_=es.prepare_candidate('insert',[{}]*5,native,[],Tools(),tmp_path)
    assert [v['residue'] for v in prepared[0][2].values()]==['ALA','MET']
    assert audits[0]['residue_mapping']==[{'author':['A',1,'A'],'pqr_number':1},
                                        {'author':['A',1,''],'pqr_number':2}]


def test_case_acceptance_rejects_complete_physics_after_runtime_exhaustion(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.case_status import case_status
    summary={'status':'VERIFIED_PANEL','phase':'production'}
    ledger=tmp_path/'surface_production_cumulative_budget.json'
    ledger.write_text(json.dumps({'active':False,'elapsed_seconds':7201.,'limit_seconds':7200.}))
    assert case_status(tmp_path,summary)=='INCOMPLETE_RUNTIME_EXHAUSTED'
    ledger.write_text(json.dumps({'active':True,'elapsed_seconds':100.,'limit_seconds':7200.}))
    assert case_status(tmp_path,summary)=='RUNNING_NOT_ACCEPTED'
    ledger.write_text(json.dumps({'active':False,'elapsed_seconds':100.,'limit_seconds':7200.}))
    assert case_status(tmp_path,summary)=='VERIFIED_PANEL'
    (tmp_path/'case_acceptance.json').write_text(json.dumps({'status':'INCOMPLETE_WORKER_TIMEOUT'}))
    assert case_status(tmp_path,summary)=='INCOMPLETE_WORKER_TIMEOUT'


def test_outer_bundle_audit_detects_mutated_disk_file_and_manifest(tmp_path):
    import zipfile,hashlib
    from ecd_epitope_scanner.scripts.audit_run_bundle import audit_run_bundle
    contents=b'real saved evidence';(tmp_path/'analysis.json').write_bytes(contents)
    manifest={'files':{'analysis.json':hashlib.sha256(contents).hexdigest()}}
    (tmp_path/'run_manifest.json').write_text(json.dumps(manifest))
    with zipfile.ZipFile(tmp_path/'ecd_epitope_design.zip','w') as archive:
        for name in ('analysis.json','run_manifest.json'):archive.write(tmp_path/name,name)
    assert audit_run_bundle(tmp_path)['status']=='PASS'
    (tmp_path/'analysis.json').write_bytes(b'changed')
    with pytest.raises(ValueError,match='ZIP differs'):audit_run_bundle(tmp_path)
    (tmp_path/'analysis.json').write_bytes(contents)
    manifest['files']['analysis.json']='incorrect';(tmp_path/'run_manifest.json').write_text(json.dumps(manifest))
    with zipfile.ZipFile(tmp_path/'ecd_epitope_design.zip','w') as archive:
        for name in ('analysis.json','run_manifest.json'):archive.write(tmp_path/name,name)
    with pytest.raises(ValueError,match='manifest differs'):audit_run_bundle(tmp_path)


def test_outer_bundle_audit_accepts_safe_nested_offline_assets(tmp_path):
    import hashlib, json, zipfile
    from ecd_epitope_scanner.scripts.audit_run_bundle import audit_run_bundle

    assets = tmp_path / 'assets' / 'vendor'
    assets.mkdir(parents=True)
    content = b'offline dependency'
    (assets / 'viewer.js').write_bytes(content)
    manifest = {'files': {'assets/vendor/viewer.js': hashlib.sha256(content).hexdigest()}}
    (tmp_path / 'run_manifest.json').write_text(json.dumps(manifest), encoding='utf-8')
    with zipfile.ZipFile(tmp_path / 'ecd_epitope_design.zip', 'w') as archive:
        archive.write(assets / 'viewer.js', 'assets/vendor/viewer.js')
        archive.write(tmp_path / 'run_manifest.json', 'run_manifest.json')

    assert audit_run_bundle(tmp_path)['status'] == 'PASS'


@pytest.mark.parametrize('unsafe_name', ['../outside.json', 'C:/outside.json', 'assets\\viewer.js', 'assets/../outside.json'])
def test_outer_bundle_audit_rejects_unsafe_manifest_paths(tmp_path, unsafe_name):
    import json
    from ecd_epitope_scanner.scripts.audit_run_bundle import audit_run_bundle

    (tmp_path / 'run_manifest.json').write_text(json.dumps({'files': {unsafe_name: 'hash'}}), encoding='utf-8')
    with pytest.raises(ValueError, match='unsafe artifact name'):
        audit_run_bundle(tmp_path)


def test_benchmark_preview_migrates_legacy_engine_preserving_af3_inputs():
    from ecd_epitope_scanner.scripts.benchmark_surface_perturbation import preview_request
    saved={'design_engine':'joint_v3','input_type':'af3_upload','structure_path':'saved_model.cif',
           'confidence_path':'saved_confidence.json','chain':'A','metadata':{'source':'original'}}
    replay=preview_request(saved,238)
    assert saved['design_engine']=='joint_v3' and replay['design_engine']=='joint_v4_ml'
    assert all(replay[key]==saved[key] for key in ('input_type','structure_path','confidence_path','chain','metadata'))
    assert replay['surface_strategy']=='surface_perturbation_v1_preview'
    assert replay['target_constructs']==8 and replay['maximum_constructs']==12


@pytest.mark.parametrize('matching_input,isolated_cache,cached_sampling,requested_sampling',[(True,True,None,None),(False,True,None,None),(True,False,None,None),(True,True,'time_rotabase_v3',None),(True,True,'time_rotabase_v3','time_rotabase_v3')])
def test_foldx_worker_resumes_frozen_models_without_restarting_a_default_prefix(tmp_path,monkeypatch,capsys,matching_input,isolated_cache,cached_sampling,requested_sampling):
    import io,sys,hashlib
    from pathlib import Path
    from ecd_epitope_scanner.integrations.workers import foldx_worker as worker
    from ecd_epitope_scanner.surface_perturbation import paired
    from ecd_epitope_scanner import stability
    binary=tmp_path/'foldx.exe';binary.write_bytes(b'fixture binary')
    source=tmp_path/'source.pdb';source.write_text('fixture input')
    saved=[{'mutations':['A1V'],'batch_path':'saved_batch','group':1,'run':i,'ddg_kcal_mol':float(i),
        'execution_layout':'candidate_isolated_v1','buildmodel_request_sha256':'fixture-receipt'} for i in range(3)]
    if cached_sampling:
        for i,pair in enumerate(saved):pair.update(sampling_protocol=cached_sampling,sampling_role='production',sampling_cohort='test',sampling_slot=i)
    manifest={'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'input_structure_sha256':hashlib.sha256(source.read_bytes()).hexdigest() if matching_input else 'different',
        'records':[{'mutations':['A1V'],'paired_models':saved},
                   {'mutations':['A1L'],'paired_models':[{'fixture':'unrequested record must not be loaded'}]}]}
    root=worker._project_root()
    manifest['execution_identity']={
        'execution_layout':'candidate_isolated_v1' if isolated_cache else 'legacy_batch',
        'worker_sha256':hashlib.sha256(Path(worker.__file__).read_bytes()).hexdigest(),
        'repaired_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'rotabase_sha256':None,
        'parameters':{'sampling_protocol':cached_sampling},
        'support_sha256':{rel:hashlib.sha256((root/rel).read_bytes()).hexdigest() for rel in
            ('stability.py','surface_perturbation/paired.py','surface_perturbation/sampling_v3.py')}}
    cache=tmp_path/'cache';cache.mkdir();(cache/'paired_model_manifest.json').write_text(json.dumps(manifest))
    output=tmp_path/'output';output.mkdir();old=output/'paired_model_manifest.json';old.write_text(json.dumps(manifest))
    request={'mutations':[['A1V']],'sequence':'A','chain_id':'A','structure_path':str(source),'output_dir':str(output),
        'parameters':{'binary':str(binary),'emit_paired_models':True,'number_of_runs':5,'paired_cache_roots':[str(cache)],'sampling_protocol':requested_sampling,'sampling_cohort':'test'}}
    calls=[];validated=[]
    monkeypatch.setattr(sys,'stdin',io.StringIO(json.dumps(request)))
    monkeypatch.setattr(stability,'repair_with_foldx',lambda *args,**kw:source)
    def score(*args,**kw):
        calls.append((args[3],kw['runs']))
        args[3].mkdir(parents=True,exist_ok=True)
        (args[3]/'sampling_execution.json').write_text('{}')
        return [9.]
    monkeypatch.setattr(stability,'score_many_with_foldx',score)
    def scan(root):
        count=calls[-1][1]
        return ([{'mutations':['A1V'],'batch_path':str(root),'group':1,'run':i,'ddg_kcal_mol':9.,
            'execution_layout':'candidate_isolated_v1','buildmodel_request_sha256':'fixture-receipt'} for i in range(count)],[])
    monkeypatch.setattr(paired,'scan_pairs',scan)
    monkeypatch.setattr(paired,'validate_pair',lambda row:validated.append(row))
    assert worker.main()==0;result=json.loads(capsys.readouterr().out)
    assert result['status']=='success',result
    reused=matching_input and isolated_cache and cached_sampling==requested_sampling and requested_sampling is not None
    assert calls[0][1]==(2 if reused else 5)
    assert calls[0][0]!=output # Old model batches are not overwritten.
    rows=result['records'][0]['paired_models']
    assert len(rows)==5 and len({(r['batch_path'],r['group'],r['run']) for r in rows})==5
    assert result['records'][0]['ddg']==(2. if reused else 9.)
    if reused:assert rows[:3]==saved
    if matching_input and (not isolated_cache or cached_sampling != requested_sampling):
        assert result['provenance']['cache_rejections']['execution_protocol_mismatch']
    if matching_input and isolated_cache and cached_sampling is requested_sampling is None:
        assert result['provenance']['cache_rejections']['deterministic_prefix_restart']==3
    assert len(validated)==(8 if matching_input and isolated_cache and cached_sampling is requested_sampling is None else 5)
    history=list((output/'paired_manifest_history').glob('*/paired_model_manifest.json'))
    assert len(history)==1 and json.loads(history[0].read_text())==manifest
    # The completed own output is a lawful cache even when the caller only
    # configured a partial external cache. A repeated call builds no models.
    monkeypatch.setattr(sys,'stdin',io.StringIO(json.dumps(request)))
    assert worker.main()==0;again=json.loads(capsys.readouterr().out)
    assert again['status']=='success' and len(calls)==1
    assert again['provenance']['new_model_counts']==[0]
    assert again['records'][0]['paired_models']==rows


def test_failed_surface_null_job_does_not_abort_other_tasks_or_reuse_old_threshold(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import calibration_jobs as jobs
    monkeypatch.setattr(jobs,'load_calibration',lambda *args:{'shape':'previously saved threshold'})
    trajectories=[];errors=[];visited=[]
    for name in ('failed','valid'):
        with jobs.isolated_job(tmp_path,name,name,trajectories,errors,lambda:5):
            visited.append(name)
            if name=='failed':raise RuntimeError('SES mesh is not closed')
            trajectories.append({'candidate':name,'local_scope_converged':True})
    assert visited==['failed','valid'] and len(trajectories)==2 and len(errors)==1
    assert jobs.usable_calibration(tmp_path,'failed')=={}
    assert jobs.usable_calibration(tmp_path,'valid')=={'shape':'previously saved threshold'}
    jobs.clear_failure(tmp_path,'failed')
    assert jobs.usable_calibration(tmp_path,'failed')=={'shape':'previously saved threshold'}
    assert len(list((tmp_path/'task_failure_history').glob('*.json')))==1
    with pytest.raises(TimeoutError):
        with jobs.isolated_job(tmp_path,'deadline','deadline',trajectories,errors,lambda:10):raise TimeoutError('CUMULATIVE_RUNTIME_EXHAUSTED')
    assert not (tmp_path/'task_failures/deadline.json').exists() and len(errors)==1
