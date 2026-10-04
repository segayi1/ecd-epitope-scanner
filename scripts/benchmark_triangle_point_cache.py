"""Exact query-state engineering test on saved native SES, without tool solves."""
import argparse,json,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from ..surface_perturbation.geometry import Mesh
from ..surface_perturbation.query_cache import memoize_mesh
from ..util import sha256_file,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case_output',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    folder=args.case_output/'surface_perturbation'
    with np.load(folder/'native_WT_SES_mesh.npz') as data:vertices=data['vertices'].copy();faces=data['faces'].copy()
    with np.load(folder/'canonical_reference.npz') as data:
        points=data['points'][::10];normals=data['normals'][::10]
    points=np.concatenate([points+1.5*normals,points-.5*normals,points[:100]])
    tools=SimpleNamespace(cache=args.output/'cache',calls={'cache_hits':0},activity=[]);tools.cache.mkdir(exist_ok=True)
    mesh=memoize_mesh(Mesh(vertices,faces),tools);mesh.signed_distance(points,return_projection=True)
    report={'scope':'one saved native mesh, exact subset query reuse; not full-pipeline acceleration or case acceptance',
        'mesh_file_sha256':sha256_file(folder/'native_WT_SES_mesh.npz'),'source_sha256':sha256_file(Path(__file__).resolve().parents[1]/'surface_perturbation/query_cache.py'),
        'no_scientific_worker_calls':True,'measurements':[]}
    for indices in (np.arange(0,len(points),2),np.arange(1,len(points),2),np.arange(len(points)-1,-1,-1)):
        subset=points[indices];original=Mesh(vertices,faces)
        started=time.perf_counter();expected=original.signed_distance(subset,return_projection=True);uncached=time.perf_counter()-started
        # A new mesh instance ensures the reuse is read from the cache on disk.
        cached=memoize_mesh(Mesh(vertices,faces),tools)
        started=time.perf_counter();actual=cached.signed_distance(subset,return_projection=True);elapsed=time.perf_counter()-started
        equal=[bool(np.array_equal(left,right,equal_nan=True)) for left,right in zip(actual,expected)]
        report['measurements'].append({'points':len(subset),'uncached_seconds':uncached,'cached_seconds':elapsed,
            'distance_valid_normals_exact':equal,'unknown_points':int((~actual[1]).sum()),'speedup':uncached/elapsed})
    report['status']='PASS' if all(all(row['distance_valid_normals_exact']) for row in report['measurements']) else 'FAIL'
    report['cache_hits']=tools.calls['cache_hits'];write_json(args.output/'report.json',report);print(json.dumps(report))
    if report['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
