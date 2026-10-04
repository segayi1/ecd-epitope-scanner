from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix
from ecd_epitope_scanner.surface_perturbation import scoped_fields as scope
from ecd_epitope_scanner.surface_perturbation.fields import difference
from ecd_epitope_scanner.surface_perturbation.chemistry import CHANNELS


def fixture(monkeypatch):
    def spatial(value,reference,tools,query_mask=None):
        active=np.ones(3,bool) if query_mask is None else query_mask
        def sampled(data):
            data=np.asarray(data,float);data[~active]=np.nan;return data
        chemistry=np.zeros((3,len(CHANNELS)));chemistry[:,0:2]=1
        chemistry[1,2:10]=np.nan # Genuine chemical uncertainty, not an unqueried cell.
        return {'sdf':sampled(np.full(3,value)),
            'shell_sdf':sampled(np.full((3,1),value)),
            'shell_valid':active[:,None].copy(),'shell_owned':np.ones((3,1),bool),
            'shell_step':1.,'chemistry':sampled(chemistry),
            'atomic_inside':sampled(np.full((3,1),value)),
            'projected_normals':sampled(np.tile([1.,0.,0.],(3,1))),
            'local_ses_area':sampled(np.full(3,value))}
    left=spatial(0,None,None,np.array([True,False,False]))
    right=spatial(1,None,None,np.array([True,False,False]))
    raw={**difference(left,right),**{target:left[source] for target,source in scope.WT_KEYS.items()},
         'topology_changed':np.array(False)}
    monkeypatch.setattr(scope,'aligned_atoms',lambda *args:(0,1,{'core_rmsd_a':0.}))
    monkeypatch.setattr(scope,'spatial_field',spatial)
    monkeypatch.setattr(scope,'reference_id',lambda reference:'reference')
    tools=SimpleNamespace(nano='nano',identity=lambda executable:'binary')
    return raw,tools


def test_missing_queries_are_really_computed_without_changing_raw_or_unknowns(tmp_path,monkeypatch):
    raw,tools=fixture(monkeypatch);before={k:v.copy() for k,v in raw.items()}
    pair={'wt_sha256':'wt','mut_sha256':'mut'};reference={'area':np.ones(3)}
    required=np.array([True,True,False])
    result=scope.complete_scope(pair,raw,{},[],reference,tools,required,tmp_path)
    assert result['shape_sdf_a'][1]==1 and result['supplemented_query_points']==1
    assert np.isnan(result['donor_support'][1]) # HIS-style abstention is retained.
    assert np.isnan(result['shape_sdf_a'][2]) # Outside scope remains unknown.
    for key in raw:np.testing.assert_array_equal(raw[key],before[key])
    monkeypatch.setattr(scope,'spatial_field',lambda *args,**kwargs:(_ for _ in ()).throw(AssertionError('valid cache reused')))
    reused=scope.complete_scope(pair,raw,{},[],reference,tools,required,tmp_path)
    np.testing.assert_array_equal(reused['shape_sdf_a'],result['shape_sdf_a'])


def test_topology_rejection_cannot_be_undone_by_more_queries(tmp_path,monkeypatch):
    raw,tools=fixture(monkeypatch);raw['topology_changed']=np.array(True)
    for name in scope.FIELD_CHANNELS:raw[name][:]=np.nan
    result=scope.complete_scope({'wt_sha256':'wt','mut_sha256':'mut'},raw,{},[],
        {'area':np.ones(3)},tools,np.array([True,True,False]),tmp_path)
    assert np.isnan(result['shape_sdf_a']).all()
    assert scope.query_present(result)[1] # WT background exists; mutation geometry remains rejected.


def test_expanded_numerical_scope_invalidates_cache_and_does_not_modify_raw(tmp_path,monkeypatch):
    from ecd_epitope_scanner.surface_perturbation import engine
    calls=[]
    def numerical(pair,native,records,reference,tools,run):
        calls.append(np.isfinite(run['wt_sdf']).copy())
        return {name:np.full(3,float(len(calls))) for name in engine.FIELD_CHANNELS}
    monkeypatch.setattr(engine,'numerical_bounds',numerical)
    monkeypatch.setattr(engine,'protocol_identity',lambda *args,**kwargs:'matched')
    pair={'wt_sha256':'wt','mut_sha256':'mut'};raw={'wt_sdf':np.zeros(3)}
    tools=SimpleNamespace(config={});reference={'area':np.ones(3)}
    first=np.array([True,False,False]);second=np.array([True,True,False])
    for required in (first,first,second):
        engine.cached_numerical_bounds(pair,{},[],reference,tools,raw,tmp_path,required_mask=required)
    assert len(calls)==2
    np.testing.assert_array_equal(calls,[first,second])
    assert np.isfinite(raw['wt_sdf']).all()
