"""Memoize exact whole triangle queries; preserve query batch and all unknowns."""
from __future__ import annotations
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import scipy
from ..util import sha256_file,write_json
from .tools import fingerprint


def _point_state(mesh,points):
    """Include all batch-computed inputs consumed by the original point loop.

    In particular preserve the exact rotated coordinates and nearest-face
    choice. A query with a different BLAS rounding or nearest tie is a miss.
    """
    if not getattr(mesh,'closed',False) or not hasattr(mesh,'ray_basis'):return None
    _,nearest=mesh.tree.query(points)
    rotated=points@mesh.ray_basis
    return np.concatenate([points,rotated],axis=1),np.asarray(nearest,dtype=np.int64)


def _point_keys(state):
    values,nearest=state
    return [row.tobytes()+index.tobytes() for row,index in zip(values,nearest)]


def _read_bank(path,contract):
    meta=path.with_suffix('.json')
    if not path.exists() or not meta.exists():return {}
    audit=json.loads(meta.read_text(encoding='utf-8'))
    if audit['contract']!=contract or audit['sha256']!=sha256_file(path):raise ValueError('triangle point cache changed')
    with np.load(path,allow_pickle=False) as data:
        keys=data['keys'];distance=data['distance'];valid=data['valid'];normals=data.get('normals')
        if keys.dtype!=np.dtype('V56') or distance.shape!=(len(keys),) or valid.shape!=distance.shape or valid.dtype!=np.dtype(bool):
            raise ValueError('triangle point cache identity or shape changed')
        if contract['projection'] and (normals is None or normals.shape!=(len(keys),3)):raise ValueError('triangle point normals shape changed')
        return {key.tobytes():(distance[i],valid[i],normals[i].copy() if normals is not None else None) for i,key in enumerate(keys)}


def _save_bank(path,contract,bank):
    rows=list(bank.values());values={'keys':np.array(list(bank),dtype='V56'),
        'distance':np.array([r[0] for r in rows]),'valid':np.array([r[1] for r in rows],bool)}
    if contract['projection']:values['normals']=np.array([r[2] for r in rows]).reshape(-1,3)
    temporary=path.with_name(path.stem+'.pending.npz');np.savez_compressed(temporary,**values);temporary.replace(path)
    write_json(path.with_suffix('.json'),{'contract':contract,'sha256':sha256_file(path)})


def memoize_mesh(mesh,tools):
    if getattr(mesh,'_exact_query_cache',False):return mesh
    digest=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()
    implementation=sha256_file(Path(__file__).parent/'geometry.py')
    original=mesh.signed_distance
    def query(points,*,return_projection=False):
        if getattr(tools,'budget',None):tools.budget.remaining()
        points=np.asarray(points,dtype=float).reshape(-1,3)
        contract={'mesh':digest,'geometry':implementation,'points':hashlib.sha256(points.tobytes()).hexdigest(),
            'shape':list(points.shape),'projection':bool(return_projection),'version':'whole_query_exact_1'}
        folder=tools.cache/'triangle_queries';folder.mkdir(exist_ok=True)
        path=folder/(fingerprint(contract)+'.npz');manifest=path.with_suffix('.json')
        point_contract={'mesh':digest,'geometry':implementation,'numpy':np.__version__,'scipy':scipy.__version__,
            'projection':bool(return_projection),'version':'exact_loop_state_1'}
        bank_path=folder/('points_'+fingerprint(point_contract)+'.npz')
        state=_point_state(mesh,points);keys=_point_keys(state) if state is not None else []
        bank=_read_bank(bank_path,point_contract) if keys else {}
        if path.exists() and manifest.exists():
            meta=json.loads(manifest.read_text(encoding='utf-8'))
            if meta['contract']!=contract or meta['sha256']!=sha256_file(path):raise ValueError('triangle query cache changed')
            with np.load(path,allow_pickle=False) as data:
                result=tuple(data[name].copy() for name in ('distance','valid','normals')[:3 if return_projection else 2])
            expected=[(len(points),),(len(points),),(len(points),3)][:len(result)]
            if [v.shape for v in result]!=expected or result[1].dtype!=np.dtype(bool):
                raise ValueError('triangle query cache shape or validity dtype changed')
            tools.calls['cache_hits']+=1
            if keys and any(key not in bank for key in keys):
                for i,key in enumerate(keys):bank[key]=(result[0][i],result[1][i],result[2][i].copy() if return_projection else None)
                _save_bank(bank_path,point_contract,bank)
            return result
        if keys and all(key in bank for key in keys):
            rows=[bank[key] for key in keys]
            result=(np.array([r[0] for r in rows]),np.array([r[1] for r in rows],bool))
            if return_projection:result=(*result,np.array([r[2] for r in rows]))
            tools.calls['cache_hits']+=1
            return result
        started=time.perf_counter();result=original(points,return_projection=return_projection)
        # A partial miss computes the original ENTIRE batch, preserving its
        # numerical path and unknowns. Never compute only the missing subset.
        if keys:
            for i,key in enumerate(keys):bank[key]=(result[0][i],result[1][i],result[2][i].copy() if return_projection else None)
            _save_bank(bank_path,point_contract,bank)
        np.savez_compressed(path,**dict(zip(('distance','valid','normals'),result)))
        write_json(manifest,{'contract':contract,'sha256':sha256_file(path)})
        if hasattr(tools,'activity'):tools.activity.append({'stage':'triangle_sdf','queries':len(points),
            'projection':return_projection,'runtime_seconds':time.perf_counter()-started})
        return result
    mesh.signed_distance=query;mesh._exact_query_cache=True
    return mesh
