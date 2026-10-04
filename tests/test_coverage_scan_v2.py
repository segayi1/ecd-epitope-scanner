import math

import pytest

from ecd_epitope_scanner.pipeline import PipelineRequest
from ecd_epitope_scanner.joint_design_v4.coverage_profile import (
    PROFILE, descriptors, resolve_threshold, select_menu, sequence_reasons,
)
from ecd_epitope_scanner.joint_design_v4.evidence import _chemistry, stability_band
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, SingleMutationEvidence
from ecd_epitope_scanner.joint_design_v4.candidate_search import _coverage_from_distances


def test_threshold_is_resolved_before_analysis_and_conflicts_fail():
    request = PipelineRequest(input_type="structure", optimization_profile=PROFILE)
    assert request.surface_strategy == "surface_perturbation_v1_preview"
    assert request.metadata["surface_rsasa_threshold"] == .20
    assert resolve_threshold("legacy", {}) == .25
    with pytest.raises(ValueError, match="Conflicting"):
        PipelineRequest(input_type="structure", optimization_profile=PROFILE,
            metadata={"surface_rsasa_threshold": .25, "joint_v4_ml": {"surface": {"rsasa_threshold": .20}}})
    for threshold in (math.nan, math.inf, -.1, 1.1):
        with pytest.raises(ValueError):
            resolve_threshold(PROFILE, {"surface_rsasa_threshold": threshold})


def test_safe_singles_can_form_unsafe_complete_sequence():
    assert not sequence_reasons("SSR", [(1, "S", "N")])
    assert not sequence_reasons("SSR", [(3, "R", "S")])
    assert sequence_reasons("SSR", [(1, "S", "N"), (3, "R", "S")]) == ["glycosylation-sequon-created"]
    assert "glycosylation-sequon-destroyed" in sequence_reasons("NSS", [(1, "N", "A")])
    assert "duplicate-position" in sequence_reasons("SSR", [(1,"S","A"), (1,"S","W")])
    assert "WT_RESIDUE_MISMATCH" in sequence_reasons("SSR", [(1,"A","W")])
    assert "native-protected-residue" in sequence_reasons("C", [(1,"C","A")])


def test_nonfinite_scores_cannot_pass_stability():
    for score in (math.nan, math.inf, -math.inf):
        assert stability_band(score, {}) == "NOT_SCORED"


def test_menu_has_deterministic_mechanism_complement():
    key = ResidueKey("m", "A", 1, uniprot_pos=1, one_letter="E")
    def options():
        return [SingleMutationEvidence(key, aa, score, eligible=True,
            chemical_descriptors=descriptors("E", aa))
            for aa, score in [("D", -.5), ("Q", -.4), ("N", -.3), ("K", .1), ("W", .2)]]
    a, b = options(), options()
    first = select_menu(a, 3)
    second = select_menu(list(reversed(b)), 3)
    assert [x.mutant_aa for x in first] == [x.mutant_aa for x in second]
    assert first[0].mutant_aa == "D"
    assert {x.mutant_aa for x in first} != {"D", "Q", "N"}
    assert _chemistry("E", "K") == _chemistry("K", "E") == "charge_reversal"
    assert descriptors("H", "K")["charge_delta"] is None


def test_proxy_fraction_matches_its_denominator():
    coverage, support = _coverage_from_distances([{1: 0, 2: 3, 3: 5}], {1, 2})
    assert coverage["geodesic_12"] == 1
    assert support["geodesic_12"] == [1, 2, 3]


def test_stage_interfaces_validate_protocol():
    for schedule in ([], [20, 5], [5, 5], [6]):
        with pytest.raises(ValueError):
            PipelineRequest(input_type="structure", calibration_repeat_schedule=schedule)
    with pytest.raises(ValueError):
        PipelineRequest(input_type="structure", export_mode="invalid")


def test_compression_preserves_panel_without_five_percent_gain():
    import numpy as np
    from scipy.sparse import csr_matrix
    from ecd_epitope_scanner.surface_perturbation.compression import compression_audit
    area = np.ones(100)
    graph = csr_matrix(np.ones((100,100)) - np.eye(100))
    hard = {str(i): np.zeros(100, bool) for i in range(6)}
    for i in range(5):
        hard[str(i)][i*16:(i+1)*16] = True
    hard['5'] = hard['0'] | hard['1']
    baseline = np.any([hard[str(i)] for i in range(5)], axis=0)
    audit = compression_audit(list(map(str,range(5))), ['0','1'], '5', hard, area, graph, baseline)
    assert audit['accepted'] and audit['added_area_a2'] == 0
    assert audit['construct_count_change'] == -1
    hard['5'][0:8] = False
    hard['5'][80:85] = True
    assert not compression_audit(list(map(str,range(5))), ['0','1'], '5', hard, area, graph, baseline)['accepted']


def test_storage_object_sharing_preserves_run_paths(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.case_resources import EvidenceStore, owned_bytes, affordable_batch
    root=tmp_path/'task';root.mkdir();(root/'results').mkdir()
    for name in ('run1.pdb','run2.pdb'):
        (root/'results'/name).write_bytes(b'identical scientific coordinates')
    manifest=EvidenceStore(root).checkpoint('T3')
    assert manifest['results/run1.pdb']['object'] == manifest['results/run2.pdb']['object']
    assert (root/'results/run1.pdb').exists() and (root/'results/run2.pdb').exists()
    assert owned_bytes(root) < 2000
    (root/'results/run1.pdb').write_bytes(b'new mutable checkpoint')
    object_path=root/manifest['results/run1.pdb']['object']
    assert object_path.read_bytes()==b'identical scientific coordinates'
    assert affordable_batch(1800,10*1024**3,30,1000)[0] == 0


def test_deterministic_unknown_does_not_request_more_models():
    import numpy as np
    from ecd_epitope_scanner.surface_perturbation.unknowns import calibration_unknowns
    value={'null_blocks':np.zeros((4,2)), 'numerical_error':np.array([float('nan'),.1]),'valid':np.array([False,True])}
    assert not calibration_unknowns({'shape':value}, [True,True])['appendable']
    value['numerical_error'][0]=.1
    assert calibration_unknowns({'shape':value}, [True,True])['appendable']


def test_full_graph_pool_is_bounded_and_complete_sequence_safe():
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
    from ecd_epitope_scanner.joint_design_v4.models import V4MutationOption
    graph={1:{2:10},2:{1:10,3:10},3:{2:10}}
    def option(p,wt,aa,score=0):
        return V4MutationOption(ResidueKey('m','A',p,uniprot_pos=p,one_letter=wt),aa,f'{wt}{p}{aa}',
            'other',score,.1,True)
    options=[option(1,'S','N'),option(1,'S','W'),option(2,'S','A'),option(3,'R','S'),option(3,'R','E'),option(2,'S','L',math.nan)]
    diagnostics={}
    candidates=generate_candidates([],options,graph,{'coverage_scan_v2':True,
        'canonical_sequence':'SSR','global_shortlist':3},diagnostics=diagnostics)
    assert len(candidates)==3
    assert diagnostics['sequence_rejected']==1
    assert all(not sequence_reasons('SSR',[(m.residue.uniprot_pos,m.residue.one_letter,m.mutant_aa) for m in c.mutations]) for c in candidates)
    assert all(0<=v<=1 for c in candidates for v in c.coverage_by_scenario.values())
    # A shrinking scoring capacity cannot discard a validated requested node.
    requested=candidates[-1].candidate_id
    replay=generate_candidates([],options,graph,{'coverage_scan_v2':True,
        'canonical_sequence':'SSR','global_shortlist':1,'mandatory_candidate_ids':[requested]})
    assert [c.candidate_id for c in replay]==[requested]
    frontier_diagnostics={}
    frontier=generate_candidates([],options,graph,{'coverage_scan_v2':True,
        'canonical_sequence':'SSR','global_shortlist':2,
        'exclude_candidate_ids':[c.candidate_id for c in candidates]},diagnostics=frontier_diagnostics)
    assert len(frontier)==2
    assert not {c.candidate_id for c in frontier}&{c.candidate_id for c in candidates}
    assert frontier_diagnostics['full_enumerated_count']==diagnostics['full_enumerated_count']
    assert frontier_diagnostics['selected_count']==2


def test_real_ctla4_safe_singles_create_nss():
    from pathlib import Path
    import json
    source=Path(__file__).parent/'fixtures/ctla4_coverage_sequence.json'
    sequence=json.loads(source.read_text(encoding='utf-8'))['sequence']
    assert sequence[48:51]=='SSR'
    assert not sequence_reasons(sequence,[(49,'S','N')])
    assert not sequence_reasons(sequence,[(51,'R','S')])
    assert 'glycosylation-sequon-created' in sequence_reasons(sequence,[(49,'S','N'),(51,'R','S')])


def test_triple_expansion_requires_gain_over_actual_double_panel():
    import numpy as np
    from ecd_epitope_scanner.surface_perturbation.compression import expansion_audit
    hard={'ab':np.array([1,1,0,0],bool),'ac':np.array([0,0,1,0],bool),
        'abc':np.array([1,1,1,0],bool)}
    # Beats either parent separately but covers no new panel surface.
    assert not expansion_audit('abc',['ab','ac'],hard,np.ones(4))['accepted']
    hard['abc'][3]=True
    assert expansion_audit('abc',['ab','ac'],hard,np.ones(4))['actual_panel_added_area_a2']==1


def test_new_api_profile_does_not_inherit_fixed_legacy_assay_budget():
    from ecd_epitope_scanner.api import JobCreate
    new=JobCreate(identifier='P16410',optimization_profile='coverage_scan_v2_preview')
    assert new.assay_profile is None and new.maximum_constructs==12
    assert new.surface_strategy=='surface_perturbation_v1_preview'
    assert JobCreate(identifier='P16410').assay_profile is not None


def test_raw_reference_identity_ignores_design_menu_but_detects_surface_change():
    import numpy as np
    from scipy.sparse import csr_matrix
    from ecd_epitope_scanner.surface_perturbation.fields import reference_id
    ref={'points':np.zeros((2,3)),'normals':np.zeros((2,3)), 'area':np.ones(2),
        'owner':np.arange(2),'designable':np.array([False,False]),'graph':csr_matrix((2,2))}
    before=reference_id(ref)
    ref['designable'][:]=True
    assert reference_id(ref)==before
    ref['points'][0,0]=1
    assert reference_id(ref)!=before


def test_expanding_numerical_checkpoint_never_mutates_saved_content_object(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.case_resources import EvidenceStore
    folder=tmp_path/'results'/'fields';folder.mkdir(parents=True)
    source=folder/'numerical_cohort.npz';source.write_bytes(b'first scope')
    source.with_suffix('.json').write_text('{}')
    manifest=EvidenceStore(tmp_path).checkpoint('before')
    target=tmp_path/manifest['results/fields/numerical_cohort.npz']['object']
    source.write_bytes(b'expanded scope')
    assert target.read_bytes()==b'first scope'


def test_boundary_native_context_conflict_stops_before_combination_pool():
    from ecd_epitope_scanner.models import ResidueRecord
    from ecd_epitope_scanner.joint_design_v4.models import SurfaceEvidence
    from ecd_epitope_scanner.joint_design_v4.evidence import build_single_evidence
    record=ResidueRecord(1,'E','A',1,rsasa=.22,surface_accessible=True)
    key=ResidueKey('m','A',1,uniprot_pos=1,one_letter='E')
    surface=SurfaceEvidence(key,native_rsasa=.22,normalized_context_rsasa=.22,
        native_context_classification_changed=False)
    def table():return build_single_evidence([record],{1:key},[{'position':1,'mutant':'K','ddg':0}],{},
        sequence='E',coverage_profile=True,surface_by_position={1:surface})
    row=next(r for r in table() if r.mutant_aa=='K')
    assert row.eligible and 'RSASA_BOUNDARY_020_025' in row.soft_warnings
    surface.normalized_context_rsasa=.18;surface.native_context_classification_changed=True
    row=next(r for r in table() if r.mutant_aa=='K')
    assert not row.eligible and 'native-normalized-sasa-threshold-conflict' in row.hard_exclusions


def test_relaxed_designability_keeps_protected_wt_surface_in_denominator():
    import numpy as np
    from ecd_epitope_scanner.models import ResidueRecord
    from ecd_epitope_scanner.surface_perturbation.geometry import Mesh
    from ecd_epitope_scanner.surface_perturbation.reference import canonical_reference
    mesh=Mesh(np.array([[0,0,0],[3,0,0],[0,3,0],[0,0,3]]),
        np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]))
    records=[ResidueRecord(1,'C','A',1),ResidueRecord(2,'E','A',2,rsasa=.22)]
    atoms={('A',1,'','CA'):(np.array([0.,0.,0.]),'CYS','C'),
        ('A',2,'','CA'):(np.array([3.,0.,0.]),'GLU','C')}
    strict=canonical_reference(mesh,atoms,records,{1,2},set())
    relaxed=canonical_reference(mesh,atoms,records,{1,2},{2})
    assert 1 in strict['owner'] and strict['area'].sum()>0
    for name in ('points','normals','area','owner'):np.testing.assert_array_equal(strict[name],relaxed[name])
    assert not relaxed['designable'][relaxed['owner']==1].any()


def test_storage_scan_tolerates_owned_atomic_ledger_replacement(tmp_path,monkeypatch):
    import os
    from ecd_epitope_scanner.surface_perturbation.case_resources import owned_bytes
    stable=tmp_path/'model.pdb';stable.write_bytes(b'preserved')
    transient=tmp_path/'case_budget.checkpoint.tmp';transient.write_bytes(b'volatile')
    original=os.scandir
    class Entry:
        def __init__(self,entry):self.entry=entry;self.path=entry.path
        def is_dir(self,**kwargs):
            result=self.entry.is_dir(**kwargs)
            if self.path==str(transient):transient.unlink(missing_ok=True)
            return result
        def stat(self):
            if self.path==str(transient):raise FileNotFoundError('atomic ledger replacement')
            return self.entry.stat()
    class Scan:
        def __init__(self,path):self.scan=original(path)
        def __enter__(self):return (Entry(e) for e in self.scan)
        def __exit__(self,*args):self.scan.close()
    monkeypatch.setattr(os,'scandir',Scan)
    assert owned_bytes(tmp_path)==len(b'preserved')


def test_storage_accounting_counts_hardlink_bytes_once(tmp_path):
    import os
    from ecd_epitope_scanner.surface_perturbation.case_resources import owned_bytes
    original=tmp_path/'original';original.write_bytes(b'one physical object')
    nested=tmp_path/'nested';nested.mkdir();os.link(original,nested/'scientific_run_alias')
    assert owned_bytes(tmp_path)==len(b'one physical object')


def test_supervisor_inventory_cache_keeps_disk_checks_and_fresh_admissions(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import case_resources as resources
    clock=[0.];scans=[];free=[100*resources.GIB]
    monkeypatch.setattr(resources,'time',SimpleNamespace(monotonic=lambda:clock[0]))
    monkeypatch.setattr(resources,'owned_bytes',lambda root:scans.append(root) or 5*resources.GIB)
    monkeypatch.setattr(resources.shutil,'disk_usage',lambda root:SimpleNamespace(free=free[0]))
    guard=resources.StorageGuard(tmp_path,inventory_interval=30)
    guard.check();clock[0]=29;guard.check()
    assert len(scans)==1
    guard.check(reserved=1)
    assert len(scans)==2
    guard.check(force_inventory=True)
    assert len(scans)==3
    free[0]=resources.GIB
    with pytest.raises(TimeoutError,match='CASE_DISK_RESERVE_UNAVAILABLE'):
        guard.check()
    assert len(scans)==3
    free[0]=100*resources.GIB;clock[0]=59
    guard.check()
    assert len(scans)==4


def test_amended_inventory_cadence_applies_to_worker_guards_and_keeps_legacy_default(tmp_path):
    import json
    from ecd_epitope_scanner.surface_perturbation.case_resources import StorageGuard,GIB
    assert StorageGuard(tmp_path).inventory_interval==2
    (tmp_path/'resource_contract_amendment.json').write_text(json.dumps({
        'activation':'ACTIVE','storage_byte_budget':30*GIB,'minimum_free_disk_bytes':10*GIB}))
    assert StorageGuard(tmp_path).inventory_interval==30
    assert StorageGuard(tmp_path,inventory_interval=2).inventory_interval==2


def test_immutable_checkpoint_links_share_inventory_but_not_extra_bytes(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import case_resources as resources
    first=tmp_path/'first.pdb';second=tmp_path/'second.pdb'
    first.write_bytes(b'first immutable model');second.write_bytes(b'second immutable model')
    original_scan=resources.owned_bytes;scans=[]
    def scan(root):
        scans.append(root);return original_scan(root)
    monkeypatch.setattr(resources,'owned_bytes',scan)
    monkeypatch.setattr(resources,'time',SimpleNamespace(monotonic=lambda:0.))
    monkeypatch.setattr(resources.shutil,'disk_usage',lambda root:SimpleNamespace(free=100*resources.GIB))
    store=resources.EvidenceStore(tmp_path)
    for path in (first,second):
        item=store.ingest(path,immutable=True)
        assert path.stat().st_ino==(tmp_path/item['object']).stat().st_ino
        assert path.read_bytes()==(tmp_path/item['object']).read_bytes()
    assert len(scans)==1
    assert store.storage.check(force_inventory=True)==len(first.read_bytes())+len(second.read_bytes())
    assert len(scans)==2


def test_mutable_checkpoint_copies_force_separate_fresh_admissions(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from ecd_epitope_scanner.surface_perturbation import case_resources as resources
    first=tmp_path/'first.npz';second=tmp_path/'second.npz'
    first.write_bytes(b'first mutable result');second.write_bytes(b'second mutable result')
    original_scan=resources.owned_bytes;scans=[]
    def scan(root):
        scans.append(root);return original_scan(root)
    monkeypatch.setattr(resources,'owned_bytes',scan)
    monkeypatch.setattr(resources,'time',SimpleNamespace(monotonic=lambda:0.))
    monkeypatch.setattr(resources.shutil,'disk_usage',lambda root:SimpleNamespace(free=100*resources.GIB))
    store=resources.EvidenceStore(tmp_path)
    for path in (first,second):
        item=store.ingest(path,immutable=False)
        assert path.stat().st_ino!=(tmp_path/item['object']).stat().st_ino
        assert path.read_bytes()==(tmp_path/item['object']).read_bytes()
    assert len(scans)==2
    assert original_scan(tmp_path)==2*(len(first.read_bytes())+len(second.read_bytes()))


def test_remote_triple_controls_keep_exact_identity_and_independent_geometry():
    from ecd_epitope_scanner.joint_design_v4.models import V4MutationOption
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates,add_pair_controls
    graph={1:{2:10,3:20},2:{1:10},3:{1:20}}
    options=[V4MutationOption(ResidueKey('m','A',p,uniprot_pos=p,one_letter='S'),aa,
        f'S{p}{aa}','other',0,0,True) for p,aa in ((1,'W'),(2,'A'),(3,'R'))]
    base={'coverage_scan_v2':True,'canonical_sequence':'SSS','global_shortlist':6}
    parents=generate_candidates([],options,graph,base)
    ids=[c.candidate_id for c in parents]
    assert len(parents)==2
    assert not any(c.order==3 for c in generate_candidates([],options,graph,{**base,'triple_parent_ids':ids}))
    candidates=add_pair_controls(generate_candidates([],options,graph,{**base,
        'triple_parent_ids':ids,'allow_remote_triples':True}),graph)
    triple=next(c for c in candidates if c.order==3)
    controls=[c for c in candidates if triple.candidate_id in c.control_for_candidate_ids]
    assert len(controls)==3 and triple.search_role=='REMOTE_TRIPLE_BACKUP'
    remote=next(c for c in controls if c.risk_components['geodesic_span_a']>24)
    assert remote.order==2 and remote.search_role=='EXACT_TRIPLE_PAIR_CONTROL'
    remote.risk_components['geodesic_span_a']=31
    assert triple.risk_components['geodesic_span_a']==30


def test_failed_attempt_requires_later_exact_phase_recovery():
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import audit_attempt_resolution
    def row(name,status,phase,ids):
        return {'stage':name,'engineering_status':status,
            'configuration':{'phase':phase,'candidate_ids':ids}}
    failed=row('failed','FAILED','calibration',['a','retired'])
    wrong=row('wrong','COMPLETED','production',['a'])
    assert audit_attempt_resolution([failed,wrong],{'retired':{}})['status']=='INCOMPLETE'
    good=row('retry','COMPLETED','calibration',['a','b'])
    audit=audit_attempt_resolution([failed,wrong,good],{'retired':{}})
    assert audit['status']=='PASS' and audit['resolved_failures'][0]['recovered_by']=='retry'
    assert audit_attempt_resolution([good,failed],{'retired':{}})['status']=='INCOMPLETE'


def test_portable_audit_rejects_consistent_outer_hash_with_corrupt_inner_field(tmp_path):
    import json
    import numpy as np
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    from ecd_epitope_scanner.surface_perturbation.archive_audit import validate_archive
    folder=tmp_path/'surface_perturbation'/'fields';folder.mkdir(parents=True)
    np.savez_compressed(folder/'field_test.npz',valid=np.array([False,True]))
    (folder/'field_test.json').write_text(json.dumps({'sha256':'wrong-inner-identity'}))
    archive=evidence_archive(tmp_path,{},tmp_path/'null')
    owned=tmp_path/'audit';owned.mkdir()
    with pytest.raises(ValueError,match='inner identity mismatch'):
        validate_archive(archive,temporary_root=owned)
    assert not list(owned.iterdir())


def test_interrupted_bundle_export_keeps_previous_valid_archive(tmp_path,monkeypatch):
    import zipfile
    from ecd_epitope_scanner.exporting import _write_bundle_archive
    (tmp_path/'analysis.json').write_text('valid checkpoint')
    _write_bundle_archive(tmp_path,['analysis.json'])
    previous=(tmp_path/'ecd_epitope_design.zip').read_bytes()
    def fail(*args,**kwargs):raise PermissionError('simulated file in use')
    monkeypatch.setattr(zipfile.ZipFile,'write',fail)
    with pytest.raises(PermissionError):_write_bundle_archive(tmp_path,['analysis.json'])
    assert (tmp_path/'ecd_epitope_design.zip').read_bytes()==previous
    assert not list(tmp_path.glob('*.pending'))


def test_unknown_audit_reports_observed_causes_without_changing_appendability():
    import numpy as np
    from ecd_epitope_scanner.surface_perturbation.unknowns import calibration_unknowns
    required=np.ones(4,bool);chemistry=np.zeros((4,14))
    chemistry[0,:]=np.nan;chemistry[1,2:10]=np.nan;chemistry[2,4:10]=np.nan
    channel={'null_blocks':np.zeros((4,4)),'numerical_error':np.ones(4),'valid':np.zeros(4,bool)}
    baseline=calibration_unknowns({'channel':channel},required)
    audit=calibration_unknowns({'channel':channel},required,runs=[{'wt_chemistry':chemistry,'topology_changed':True}])
    modes=audit['observed_failure_modes']
    assert audit['appendable']==baseline['appendable']
    assert modes['topology_rejected_runs']==1 and modes['required_query_missing_points']==1
    assert modes['wt_histidine_protonation_unknown_points']==1 and modes['wt_direction_unknown_points']==1


def test_null_extensions_cannot_repair_two_topology_failures_in_fixed_five():
    import numpy as np
    from ecd_epitope_scanner.surface_perturbation.unknowns import calibration_unknowns
    channel={'null_blocks':np.zeros((4,1)),'numerical_error':np.ones(1),'valid':np.zeros(1,bool)}
    calibration={'channel':channel};required=np.ones(1,bool)
    runs=[{'topology_changed':i in (0,1)} for i in range(5)]
    audit=calibration_unknowns(calibration,required,runs=runs)
    assert not audit['appendable'] and audit['extension_stop_reason']=='FIRST_FIVE_TOPOLOGY_CANNOT_SUPPORT_FOUR_REPEATS'
    # Null-only later failures do not alter eligibility of the original five.
    runs=[{'topology_changed':False} for _ in range(5)]+[{'topology_changed':True} for _ in range(3)]
    assert calibration_unknowns(calibration,required,runs=runs)['appendable']
def test_checkpoint_proxy_pair_diagnostics_references_are_complete_and_verified(tmp_path):
    import copy,json
    from ecd_epitope_scanner.joint_design_v4.reporting import externalize_checkpoint_pair_diagnostics
    from ecd_epitope_scanner.util import write_json,sha256_file
    hierarchy={'raw':{'near_pairs':[['a','b',0.1]],'pairwise_distances':[['a','b',0.1]],
        'groups':[['a','b']]},'interpretation':'proxy only'}
    artifact=tmp_path/'hypothesis_hierarchy_v411.json'
    write_json(artifact,{'resolution_by_level':hierarchy})
    physical={'goal_status':'INSUFFICIENT_VERIFIED_CONSTRUCTS','actual_constructs':1}
    analysis={'partitioning':{'selected_panel':{'hierarchy_resolution':copy.deepcopy(hierarchy)}},
        'joint_v4_ml':{'selected_panel':{'hierarchy_resolution':copy.deepcopy(hierarchy)}},
        'design_summary':{'surface_perturbation':physical}}
    externalize_checkpoint_pair_diagnostics(analysis,tmp_path)
    assert analysis['design_summary']['surface_perturbation']==physical
    for key in ('partitioning','joint_v4_ml'):
        raw=analysis[key]['selected_panel']['hierarchy_resolution']['raw']
        assert raw['groups']==hierarchy['raw']['groups'] and 'near_pairs' not in raw
        for name in ('near_pairs','pairwise_distances'):
            reference=raw[name+'_reference'];assert reference['sha256']==sha256_file(artifact)
            restored=json.loads(artifact.read_text(encoding='utf-8'))
            for field in reference['json_path']:restored=restored[field]
            assert restored==hierarchy['raw'][name] and len(restored)==reference['count']
    bad={'partitioning':{'selected_panel':{'hierarchy_resolution':{'raw':{}}}}}
    with pytest.raises(ValueError,match='identity mismatch'):
        externalize_checkpoint_pair_diagnostics(bad,tmp_path)

def test_saved_frontier_replenishes_bounded_unmeasured_candidates(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import replenish_saved_frontier
    from ecd_epitope_scanner.util import write_json
    output=tmp_path/'results';output.mkdir()
    rows=[]
    for position,wt,aa in [(1,'S','W'),(2,'S','A'),(3,'R','E')]:
        rows.append({'residue':{'model_id':'m','chain_id':'A','auth_seq_id':position,
            'uniprot_pos':position,'one_letter':wt},'mutant_aa':aa,'chemistry_class':'other',
            'spurs_ddg_kcal_mol':0.,'eligible':True,'review_required':False,'menu_selected':True})
    write_json(output/'mutation_menu_v2.json',{'rows':rows})
    write_json(output/'uniprot.json',{'sequence':{'value':'SSR'}})
    write_json(output/'surface_graph_v4.json',{'nodes':[{'uniprot_pos':i} for i in (1,2,3,4)],
        'edges':[{'source':1,'target':2,'distance_a':10.},{'source':2,'target':3,'distance_a':10.}]})
    first=replenish_saved_frontier(tmp_path,set(),2)
    assert len(first)==2
    excluded={c['candidate_id'] for c in first}
    second=replenish_saved_frontier(tmp_path,excluded,2)
    assert len(second)==1 and second[0]['candidate_id'] not in excluded
    assert all(c['order']==2 for c in first+second)
    assert not replenish_saved_frontier(tmp_path,excluded|{second[0]['candidate_id']},2)


def test_completed_history_includes_legacy_pilots_without_whitelist():
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import completed_candidate_ids
    reports=[{'engineering_status':'COMPLETED','configuration':{'phase':'production'},
        'physical':{'candidate_site_evidence':{'old':{}}}},
        {'engineering_status':'COMPLETED','configuration':{'phase':'production','candidate_ids':['new']}},
        {'engineering_status':'ENGINEERING_FAILURE','configuration':{'phase':'production','candidate_ids':['failed']}},
        {'engineering_status':'COMPLETED','configuration':{'phase':'calibration','candidate_ids':['not_produced']}}]
    assert completed_candidate_ids(reports)=={'old','new'}


def test_final_partial_panel_apbs_scheduling_does_not_qualify_coverage():
    from ecd_epitope_scanner.surface_perturbation.final_evidence import apbs_panel_ready
    from ecd_epitope_scanner.surface_perturbation.combo_contract import goal_state
    assert not apbs_panel_ready(combo=True,construct_count=1,hard_coverage=.18)
    assert apbs_panel_ready(combo=True,construct_count=1,hard_coverage=.18,final_panel_audit=True)
    assert not apbs_panel_ready(combo=True,construct_count=0,hard_coverage=None,final_panel_audit=True)
    assert apbs_panel_ready(combo=True,construct_count=4,hard_coverage=.8)
    assert goal_state([object()],.18,{},audit_complete=True)=='INSUFFICIENT_VERIFIED_CONSTRUCTS'


def test_final_apbs_stage_reuses_all_pair_identities_without_extension(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import run_ctla4_coverage_v2 as runner
    from ecd_epitope_scanner.util import write_json
    physical=tmp_path/'results/surface_perturbation';physical.mkdir(parents=True)
    (tmp_path/'stages').mkdir()
    write_json(tmp_path/'case_budget.json',{'limit_seconds':14400,'elapsed_seconds':12000,'active':False})
    write_json(tmp_path/'results/surface_perturbation_preview.json',{'selected_candidate_ids':['a']})
    write_json(physical/'paired_manifest.json',{'a':[{'wt_sha256':'wt','mut_sha256':'mut','group':'g','run':1,'batch_path':'batch1'}]})
    calls=[]
    def launch(case,assets,stage,config,seconds):
        calls.append((config,seconds));return {'stage':stage,'engineering_status':'COMPLETED'}
    monkeypatch.setattr(runner,'launch',launch)
    reports=[{'engineering_status':'COMPLETED','configuration':{'phase':'production','schedule':[5,10,20],
        'calibration_extension_candidate_ids':['a']}}]
    assert len(runner.final_panel_audit(tmp_path,tmp_path,reports))==2
    config,seconds=calls[0]
    assert config['final_panel_audit'] and config['schedule']==[5]
    assert config['calibration_extension_candidate_ids']==[] and seconds==2100
    assert runner.read(tmp_path/'final_apbs_model_reuse_audit.json')['new_paired_models']==0


def test_unlimited_budget_is_explicit_and_preserves_charged_history(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.runtime import RunBudget,read_budget
    from ecd_epitope_scanner.surface_perturbation.case_resources import affordable_batch,GIB
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import case_limit,case_remaining
    from ecd_epitope_scanner.util import write_json
    assert case_limit(tmp_path)==14400
    write_json(tmp_path/'resource_contract_amendment.json',{'time_budget':'UNLIMITED'})
    write_json(tmp_path/'case_budget.json',{'limit_seconds':None,'elapsed_seconds':15000,'active':False})
    assert case_limit(tmp_path) is None and math.isinf(case_remaining(tmp_path))
    budget=RunBudget(tmp_path/'case_budget.json',None)
    assert budget.elapsed>=15000 and math.isinf(budget.remaining())
    budget.save(False)
    assert read_budget(tmp_path/'case_budget.json')['limit_seconds'] is None
    assert affordable_batch(math.inf,10*GIB,360,10*1024**2)==(8,256)


def test_explicit_time_amendment_requires_closed_boundary_and_preserves_old_ledger(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import activate_resource_amendment,read
    from ecd_epitope_scanner.util import write_json
    write_json(tmp_path/'resource_contract_amendment.json',{'time_budget':'UNLIMITED','authorization':'user'})
    old={'limit_seconds':14400,'elapsed_seconds':13000,'active':True}
    write_json(tmp_path/'case_budget.json',old)
    with pytest.raises(ValueError,match='closed stage'):activate_resource_amendment(tmp_path)
    old['active']=False;write_json(tmp_path/'case_budget.json',old)
    activate_resource_amendment(tmp_path)
    assert read(tmp_path/'case_budget_before_unlimited_amendment.json')==old
    current=read(tmp_path/'case_budget.json')
    assert current['limit_seconds'] is None and current['elapsed_seconds']>=13000
    assert read(tmp_path/'resource_contract_amendment.json')['activation']=='ACTIVE'


def test_full_replacement_frontier_includes_eligible_unselected_backups(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import replenish_saved_frontier
    from ecd_epitope_scanner.util import write_json
    output=tmp_path/'results';output.mkdir()
    rows=[{'residue':{'model_id':'m','chain_id':'A','auth_seq_id':i,'uniprot_pos':i,'one_letter':'S'},
        'mutant_aa':'A','chemistry_class':'other','spurs_ddg_kcal_mol':0.,
        'eligible':True,'review_required':False,'menu_selected':i==1} for i in (1,2)]
    write_json(output/'mutation_menu_v2.json',{'rows':rows})
    write_json(output/'uniprot.json',{'sequence':{'value':'SS'}})
    write_json(output/'surface_graph_v4.json',{'nodes':[{'uniprot_pos':i} for i in (1,2)],
        'edges':[{'source':1,'target':2,'distance_a':10.}]})
    assert not replenish_saved_frontier(tmp_path,set(),2)
    write_json(tmp_path/'resource_contract_amendment.json',{'exhaustion_scope':'ALL_ELIGIBLE_LOCAL_DOUBLES'})
    assert len(replenish_saved_frontier(tmp_path,set(),2))==1
    assert not rows[1]['menu_selected']


def test_site_pair_priority_selects_one_variant_and_excludes_completed_pairs():
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
    from ecd_epitope_scanner.joint_design_v4.models import V4MutationOption
    options=[V4MutationOption(ResidueKey('m','A',position,uniprot_pos=position,one_letter=wt),
        aa,f'{wt}{position}{aa}','other',0.,0.,True)
        for position,wt,letters in [(1,'S','AW'),(2,'S','AW'),(3,'R','EK')] for aa in letters]
    graph={1:{2:10},2:{1:10,3:10},3:{2:10}}
    params={'coverage_scan_v2':True,'canonical_sequence':'SSR','global_shortlist':32,'site_pair_priority':True}
    diagnostic={}
    first=generate_candidates([],options,graph,params,diagnostics=diagnostic)
    pairs=[tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in first]
    assert len(first)==3 and len(set(pairs))==3
    assert diagnostic['full_enumerated_count']==12 and diagnostic['excluded_screened_count']==0
    second=generate_candidates([],options,graph,{**params,'exclude_site_pairs':[[1,2]]})
    assert len(second)==2 and all(tuple(sorted(o.residue.uniprot_pos for o in c.mutations))!=(1,2) for c in second)


def test_amended_storage_preserves_disk_headroom_without_changing_legacy_defaults(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import case_resources as resources
    from ecd_epitope_scanner.util import write_json
    from types import SimpleNamespace
    legacy=resources.StorageGuard(tmp_path)
    assert legacy.maximum==20*resources.GIB and legacy.solve_limit==18*resources.GIB
    write_json(tmp_path/'resource_contract_amendment.json',{'activation':'ACTIVE',
        'storage_byte_budget':30*resources.GIB,'minimum_free_disk_bytes':10*resources.GIB})
    guard=resources.StorageGuard(tmp_path)
    assert guard.solve_limit==28*resources.GIB
    monkeypatch.setattr(resources,'owned_bytes',lambda _:25*resources.GIB)
    monkeypatch.setattr(resources.shutil,'disk_usage',lambda _:SimpleNamespace(free=9*resources.GIB))
    with pytest.raises(TimeoutError,match='DISK_RESERVE'):guard.check()
    monkeypatch.setattr(resources.shutil,'disk_usage',lambda _:SimpleNamespace(free=11*resources.GIB))
    assert guard.check()==25*resources.GIB


def test_cumulative_pair_index_retains_runs_and_rejects_identity_mutation(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import preserve_pair_index
    from ecd_epitope_scanner.util import write_json
    path=tmp_path/'results/surface_perturbation/paired_manifest.json'
    one={'batch_path':'batch','group':'g','run':1,'wt_sha256':'same','mut_sha256':'same'}
    two={**one,'run':2}
    write_json(path,{'old':[one,two]})
    assert len(preserve_pair_index(tmp_path)['old'])==2
    write_json(path,{'new':[{**one,'group':'new'}]})
    assert set(preserve_pair_index(tmp_path))=={'old','new'}
    write_json(path,{'old':[{**one,'mut_sha256':'changed'}]})
    with pytest.raises(ValueError,match='identity changed'):preserve_pair_index(tmp_path)


def test_bounded_active_pool_preserves_formal_positive_extendable_and_controls():
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import bounded_retained_ids
    pairs={str(i):[] for i in range(10)}
    summary={'independent_K_verified_candidate_ids':['1'],
        'shadow_positive_lower_bound':{'candidate_ids':['2']}}
    assert bounded_retained_ids(pairs,summary,['3'],['4'],maximum=4)=={'1','2','3','4'}
    with pytest.raises(ValueError,match='Necessary active'):
        bounded_retained_ids(pairs,summary,['3'],['4'],maximum=3)


def test_completed_snapshots_update_search_history_between_batches(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import record_completed_snapshot,read
    from ecd_epitope_scanner.util import write_json
    write_json(tmp_path/'results/surface_perturbation/paired_manifest.json',{})
    write_json(tmp_path/'results/panel_candidates.json',[{'candidate_id':'first','order':2,'mutations':[]}])
    report={'stage':'one','engineering_status':'COMPLETED','configuration':{'phase':'production'},
        'physical':{'candidate_site_evidence':{'first':{'complete':False,'reason':'UNKNOWN'}}}}
    record_completed_snapshot(tmp_path,report)
    write_json(tmp_path/'results/panel_candidates.json',[{'candidate_id':'second','order':2,'mutations':[]}])
    report={**report,'stage':'two','physical':{}}
    history=record_completed_snapshot(tmp_path,report)
    assert set(history)=={'first','second'}
    assert history['first']['candidate_site_evidence']['reason']=='UNKNOWN'
    assert read(tmp_path/'screened_candidate_history.json')['second']['last_completed_stage']=='two'
    with pytest.raises(ValueError,match='completed production'):
        record_completed_snapshot(tmp_path,{**report,'engineering_status':'ENGINEERING_FAILURE'})


def test_site_pair_scope_is_frozen_across_bounded_frontiers(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import replenish_saved_frontier,read
    from ecd_epitope_scanner.util import write_json
    output=tmp_path/'results';output.mkdir()
    rows=[{'residue':{'model_id':'m','chain_id':'A','auth_seq_id':i,'uniprot_pos':i,'one_letter':wt},
        'mutant_aa':aa,'chemistry_class':'other','spurs_ddg_kcal_mol':0.,'eligible':True,
        'review_required':False,'menu_selected':True} for i,wt,aa in [(1,'S','A'),(2,'S','W'),(3,'R','E')]]
    write_json(output/'mutation_menu_v2.json',{'rows':rows})
    write_json(output/'uniprot.json',{'sequence':{'value':'SSR'}})
    write_json(output/'surface_graph_v4.json',{'nodes':[{'uniprot_pos':i} for i in (1,2,3)],
        'edges':[{'source':1,'target':2,'distance_a':10.},{'source':2,'target':3,'distance_a':10.}]})
    write_json(tmp_path/'resource_contract_amendment.json',{'exhaustion_scope':'ELIGIBLE_LOCAL_SITE_PAIRS'})
    first=replenish_saved_frontier(tmp_path,set(),2)
    contract=read(tmp_path/'site_pair_search_contract.json')
    assert len(contract['site_pairs'])==3
    write_json(tmp_path/'screened_candidate_history.json',{c['candidate_id']:{'candidate':c} for c in first})
    second=replenish_saved_frontier(tmp_path,{c['candidate_id'] for c in first},2)
    assert len(second)==1 and read(tmp_path/'site_pair_search_contract.json')==contract


def test_incomplete_attempt_without_report_does_not_get_overwritten(tmp_path):
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import stage_attempt_number
    stages=tmp_path/'stages';stages.mkdir()
    assert stage_attempt_number(tmp_path,'T6_resume_72')==1
    (stages/'T6_resume_72_1_calibration_request.json').write_text('{}')
    assert stage_attempt_number(tmp_path,'T6_resume_72')==2
    (stages/'T6_resume_72_2_production.log').write_text('partial')
    assert stage_attempt_number(tmp_path,'T6_resume_72')==3


def test_stability_failure_schedules_one_nominal_stability_backup_not_all_variants(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import run_ctla4_coverage_v2 as runner
    from ecd_epitope_scanner.util import write_json
    from ecd_epitope_scanner.joint_design_v4.models import V4MutationOption
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
    from ecd_epitope_scanner.models import to_dict
    rows=[{'residue':{'model_id':'m','chain_id':'A','auth_seq_id':i,'uniprot_pos':i,'one_letter':'S'},
        'mutant_aa':aa,'chemistry_class':'other','spurs_ddg_kcal_mol':score,'eligible':True,
        'review_required':False,'menu_selected':True} for i in (1,2) for aa,score in [('W',-.1),('A',-1.)]]
    options=[V4MutationOption(ResidueKey(**row['residue']),row['mutant_aa'],f"S{row['residue']['uniprot_pos']}{row['mutant_aa']}",
        'other',row['spurs_ddg_kcal_mol'],0,True) for row in rows]
    pool=generate_candidates([],options,{1:{2:10},2:{1:10}},
        {'coverage_scan_v2':True,'canonical_sequence':'SS','global_shortlist':16})
    parent=next(to_dict(c) for c in pool if all(o.mutant_aa=='W' for o in c.mutations))
    parent['stability']={'consensus_status':'COMBINATION_REJECT'}
    output=tmp_path/'results';output.mkdir()
    write_json(output/'panel_candidates.json',[parent]);write_json(output/'mutation_menu_v2.json',{'rows':rows})
    write_json(output/'uniprot.json',{'sequence':{'value':'SS'}})
    write_json(output/'surface_graph_v4.json',{'nodes':[{'uniprot_pos':i} for i in (1,2)],
        'edges':[{'source':1,'target':2,'distance_a':10.}]})
    history={parent['candidate_id']:{'candidate':parent}}
    write_json(tmp_path/'screened_candidate_history.json',history)
    monkeypatch.setattr(runner,'rank_pending',lambda case,catalog,ids:ids)
    backup=runner.choose_failure_backup(tmp_path,set(history))
    assert backup['candidate_id']!=parent['candidate_id']
    assert all(m['mutant_aa']=='A' for m in backup['mutations'])
    assert backup['conditional_backup_parent_id']==parent['candidate_id']
    history[backup['candidate_id']]={'candidate':backup};write_json(tmp_path/'screened_candidate_history.json',history)
    assert runner.choose_failure_backup(tmp_path,set(history)) is None


def test_fatal_or_missing_physical_stage_cannot_mark_site_pairs_processed():
    from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import physical_stage_engineering_status
    for physical in (None,{}, {'status':'INCOMPLETE'}):
        assert physical_stage_engineering_status(physical)=='ENGINEERING_FAILURE'
    for status in ('CALIBRATION_EXPANSION_REQUIRED','NO_VERIFIED_PANEL','CALIBRATION_INSUFFICIENT'):
        assert physical_stage_engineering_status({'status':status})=='COMPLETED'

