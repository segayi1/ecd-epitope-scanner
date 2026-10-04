"""Actually sample missing candidate support; preserve source raw fields and unknowns."""
import hashlib,json
from pathlib import Path
import numpy as np
from ..util import sha256_file,write_json
from .fields import aligned_atoms,spatial_field,difference,reference_id,IMPLEMENTATION_ID,null_differences
from .null_store import FIELD_CHANNELS,protocol_identity as raw_protocol,load_calibration
from .calibration import calibrate_channel
from .tools import fingerprint

SCOPE_IMPLEMENTATION=sha256_file(Path(__file__))
WT_KEYS={'wt_sdf':'sdf','wt_shell_sdf':'shell_sdf','wt_shell_valid':'shell_valid',
    'wt_shell_owned':'shell_owned','wt_chemistry':'chemistry','wt_normals':'projected_normals',
    'wt_local_area':'local_ses_area','wt_atomic_inside':'atomic_inside'}


def query_present(run):
    # Nonpolar/polar channels are defined at queried rows, even when HIS
    # chemistry or donor/acceptor directions deliberately abstain.
    return np.isfinite(run['wt_chemistry'][:,1])


def protocol_identity(pairs,reference,tools,*,implementation_id=None):
    return fingerprint({'paired_protocol':raw_protocol(pairs,reference,tools,implementation_id=implementation_id),
        'query_contract':'complete_actual_candidate_support_1','scope_implementation':SCOPE_IMPLEMENTATION})


def complete_scope(pair,raw,native,records,reference,tools,required,cache):
    required=np.asarray(required,bool)
    if required.shape!=(len(reference['area']),):raise ValueError('required query scope/reference mismatch')
    scope=hashlib.sha256(required.tobytes()).hexdigest()
    digest=hashlib.sha256()
    for name,value in sorted(raw.items()):
        value=np.asarray(value);digest.update(name.encode());digest.update(str(value.dtype).encode())
        digest.update(str(value.shape).encode());digest.update(value.tobytes())
    contract={'reference':reference_id(reference),'wt':pair['wt_sha256'],'mut':pair['mut_sha256'],
        'implementation_sha256':IMPLEMENTATION_ID,'scope_implementation':SCOPE_IMPLEMENTATION,
        'query_scope_sha256':scope,'raw_array_sha256':digest.hexdigest(),'nano':tools.identity(tools.nano)}
    path=Path(cache)/('scoped_'+fingerprint(contract)+'.npz');metadata=path.with_suffix('.json')
    if path.exists() and metadata.exists():
        audit=json.loads(metadata.read_text(encoding='utf-8'))
        if audit['contract']!=contract or audit['sha256']!=sha256_file(path):raise ValueError('scoped field cache fingerprint mismatch')
        with np.load(path,allow_pickle=False) as data:return {k:data[k].copy() for k in data.files}
    missing=required&~query_present(raw)
    result=dict(raw)
    if missing.any():
        wt,mut,audit=aligned_atoms(pair,native,records)
        if audit['core_rmsd_a']>.5:raise ValueError('REFERENCE_UNRELIABLE: query supplement core RMSD')
        left=spatial_field(wt,reference,tools,query_mask=missing)
        right=spatial_field(mut,reference,tools,query_mask=missing)
        delta=difference(left,right)
        for target,source in WT_KEYS.items():
            result[target]=np.asarray(raw[target]).copy();result[target][missing]=left[source][missing]
        for name,value in delta.items():
            result[name]=np.asarray(raw[name]).copy()
            # A topology rejection is global; extra queries cannot undo it.
            result[name][missing]=np.nan if bool(raw['topology_changed']) and not bool(raw.get('local_channel_qc',False)) else value[missing]
    result['query_scope_sha256']=np.array(scope)
    result['supplemented_query_points']=np.array(int(missing.sum()))
    path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**result)
    write_json(metadata,{'contract':contract,'sha256':sha256_file(path),
        'supplemented_points':int(missing.sum()),'supplemented_area_a2':float(reference['area'][missing].sum()),
        'interpretation':'actual paired WT/Mut queries; original raw arrays unchanged; physical ambiguity remains unknown'})
    return result


def save_calibration(folder,pairs,runs,reference,tools,numerical,*,required_mask=None):
    protocol=protocol_identity(pairs,reference,tools);null=null_differences(runs)
    summaries={};arrays={}
    for channel in FIELD_CHANNELS:
        value=calibrate_channel(np.array([r[channel] for r in null]),numerical[channel],protocol,channel,
            required_mask=required_mask,area=reference['area'])
        summaries[channel]={k:v for k,v in value.items() if not isinstance(v,np.ndarray)}
        arrays.update({f'{channel}__{name}':data for name,data in value.items() if isinstance(data,np.ndarray)})
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    path=folder/(protocol+'.npz');np.savez_compressed(path,**arrays)
    write_json(path.with_suffix('.json'),{'protocol_id':protocol,'reference_id':reference_id(reference),
        'channels':summaries,'runs':len(runs),'array_sha256':sha256_file(path),
        'query_contract':'complete_actual_candidate_support_1',
        'query_scope_sha256':str(runs[0]['query_scope_sha256']),
        'pairs':[{'wt_sha256':p['wt_sha256'],'mut_sha256':p['mut_sha256'],'pdb_list_sha256':p['pdb_list_sha256']} for p in pairs[:len(runs)]],
        'dependencies':'consecutive WT differences share references; not independent biological samples'})
    return load_calibration(folder,protocol)
