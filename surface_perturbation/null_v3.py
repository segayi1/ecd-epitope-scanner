"""Candidate-matched disjoint WT differences with untouched holdout data."""
import json
from pathlib import Path

import numpy as np

from ..util import write_json, sha256_file
from .fields import aligned_atoms, spatial_field, difference, reference_id, IMPLEMENTATION_ID
from .scoped_fields import WT_KEYS
from .null_store import FIELD_CHANNELS, load_calibration
from .calibration import calibrate_channel, resolved_changes
from .tools import fingerprint
from .sampling_v3 import validate_roles


def identity(production,reference,tools,required):
    return fingerprint({'version':'disjoint_train10_holdout10_v3',
        'sampling':'time_rotabase_v3','raw_fields':IMPLEMENTATION_ID,
        'source':sha256_file(Path(__file__)),
        'classification':{name:sha256_file(Path(__file__).parent/name) for name in
                          ('calibration.py','directions.py','classification.py')},
        'reference':reference_id(reference),'query_scope':fingerprint(np.asarray(required,bool).tolist()),
        'production':[{k:p.get(k) for k in ('input_sha256','binary_sha256','batch_path','group','run','wt_sha256','mut_sha256','sampling_slot')} for p in production],
        'nano':tools.identity(tools.nano)})


def wt_field(pair,native,records,reference,tools,required,cache):
    contract={'wt':pair['wt_sha256'],'reference':reference_id(reference),
        'nano':tools.identity(tools.nano),'field_source':IMPLEMENTATION_ID,
        'scope':fingerprint(np.asarray(required,bool).tolist()),'version':'dedicated_WT_v3'}
    path=Path(cache)/('null_WT_'+fingerprint(contract)+'.npz')
    if path.exists():
        metadata=json.loads(path.with_suffix('.json').read_text(encoding='utf-8'))
        if metadata['contract']!=contract or sha256_file(path)!=metadata['sha256']:
            raise ValueError('WT calibration object corrupted')
        with np.load(path,allow_pickle=False) as data:return {k:data[k].copy() for k in data.files}
    wt,_,audit=aligned_atoms(pair,native,records)
    if audit['core_rmsd_a']>.5:raise ValueError('calibration WT core RMSD')
    field=spatial_field(wt,reference,tools,query_mask=required)
    arrays={name:field[key] for name,key in WT_KEYS.items()}
    path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**arrays)
    write_json(path.with_suffix('.json'),{'contract':contract,'sha256':sha256_file(path)})
    return arrays


def disjoint_blocks(rows):
    if len(rows)!=10:raise ValueError('ten WT models required per null role')
    return [difference({k:rows[i][name] for name,k in WT_KEYS.items()},
                       {k:rows[i+1][name] for name,k in WT_KEYS.items()}) for i in range(0,10,2)]


def build(folder,protocol,train,holdout,numerical,required,reference,pairs):
    if len(train)!=5 or len(holdout)!=5:raise ValueError('five disjoint blocks per split')
    arrays={};summaries={}
    for name in FIELD_CHANNELS:
        blocks=np.array([r[name] for r in train]);held=np.array([r[name] for r in holdout])
        value=calibrate_channel(blocks,numerical[name],protocol,name,required_mask=required,area=reference['area'])
        if 'valid' not in value:raise ValueError('null channel incomplete')
        fake_hit,_=resolved_changes(held,value,protocol)
        value['valid'] &= np.isfinite(held).all(axis=0) & ~fake_hit
        value.update(positive_lower_bound=True,holdout_blocks=held,holdout_false_positive=fake_hit,
                     holdout_rule='fixed four-of-five; held data never used to tune thresholds')
        # PARTIAL is sufficient for positive lower bounds; status reflects new mask.
        value['status']='COMPLETE' if value['valid'].all() else 'PARTIAL'
        summaries[name]={k:v for k,v in value.items() if not isinstance(v,np.ndarray)}
        arrays.update({name+'__'+k:v for k,v in value.items() if isinstance(v,np.ndarray)})
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    path=folder/(protocol+'.npz');np.savez_compressed(path,**arrays)
    write_json(path.with_suffix('.json'),{'protocol_id':protocol,'reference_id':reference_id(reference),
        'channels':summaries,'array_sha256':sha256_file(path),'runs':20,
        'pairs':pairs,'dependencies':'disjoint WT pairs; sidechain sampling, not independent biological samples',
        'interpretation':'empirical holdout diagnostic, no global error-control guarantee'})
    return load_calibration(folder,protocol)


def obtain(candidate,production,previous,supplement,native,records,reference,tools,required,numerical,folder,cache):
    protocol=identity(production,reference,tools,required)
    cohort=protocol
    existing=load_calibration(folder,protocol)
    roles=[]
    for role in ('null_train','null_holdout'):
        rows=validate_roles(previous,role,10,cohort)
        if len(rows)<10:
            rows+=supplement(candidate,10-len(rows),role=role,cohort=cohort,slot_start=len(rows))
            previous.extend(p for p in rows if (p['batch_path'],p['group'],p['run']) not in
                            {(q['batch_path'],q['group'],q['run']) for q in previous})
        if len(validate_roles(rows,role,10,cohort))!=10:raise ValueError('null role slots incomplete')
        roles.append(rows)
    production_ids={(p['batch_path'],p['group'],p['run']) for p in production}
    null_ids=[(p['batch_path'],p['group'],p['run']) for rows in roles for p in rows]
    if len(set(null_ids))!=20 or set(null_ids)&production_ids:raise ValueError('null/production role overlap')
    pairs=[{k:p.get(k) for k in ('wt_sha256','mut_sha256','batch_path','group','run','sampling_role','sampling_slot')} for rows in roles for p in rows]
    if existing:
        manifest=json.loads((Path(folder)/(protocol+'.json')).read_text(encoding='utf-8'))
        if manifest['pairs']!=pairs:raise ValueError('null cohort source identity changed')
        return protocol,existing
    field_roles=[[wt_field(p,native,records,reference,tools,required,cache) for p in rows] for rows in roles]
    return protocol,build(folder,protocol,*[disjoint_blocks(rows) for rows in field_roles],numerical,required,reference,pairs)
