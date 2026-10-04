"""Portable physical evidence archive, excluding intermediate tool caches."""
from __future__ import annotations
import json
import hashlib
import os
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from ..util import sha256_file


def evidence_archive(output, pairs, calibration_dir):
    output=Path(output); folder=output/'surface_perturbation'
    files={p.relative_to(folder).as_posix():p for p in folder.glob('*') if p.is_file() and p.suffix in {'.json','.npz','.pqr','.pdb'}}
    route_audit=output/'surface_perturbation_candidate_routes.json'
    if route_audit.is_file(): files[route_audit.name]=route_audit
    status=output/'surface_perturbation_preview.json'
    if status.is_file():files[status.name]=status
    for name in ('surface_perturbation_local_constructs.json','surface_perturbation_local_constructs.csv'):
        if (output/name).is_file():files[name]=output/name
    for name in ('frozen_input_manifest.json','launch_source_manifest.json','source_confidence.json','raw_null_import_manifest.json',
                 'native_input.cif','native_input.pdb','structure_for_energy.pdb','residue_mapping.csv'):
        if (output/name).is_file(): files['source/'+name]=output/name
    for p in (folder/'fields').glob('*.npz'): files['fields/'+p.name]=p
    for p in (folder/'fields').glob('*.json'): files['fields/'+p.name]=p
    calibration_dir=Path(calibration_dir)
    for p in calibration_dir.rglob('*'):
        if p.is_file() and p.suffix in {'.json','.npz'}:
            # Failed-task masks must travel with null arrays. Omitting them
            # would turn a blocked old partial calibration into reusable evidence.
            files['calibration/'+p.relative_to(calibration_dir).as_posix()]=p
    portable={}
    for candidate,rows in pairs.items():
        portable[candidate]=[]
        for index,row in enumerate(rows):
            item=dict(row)
            item['source_execution_identity']={key:row.get(key) for key in ('batch_path','group','run')}
            if row.get('ddg_kcal_mol') is not None and not item.get('ddg_file_path'):
                # Older lawful pairs retained the DDG value and batch directory,
                # but predate the explicit table-path manifest contract.
                from ..stability import _parse_dif_rows
                matching=[path for path in Path(row['batch_path']).glob('Dif_*.fxout')
                    if any(name==Path(row['mut_path']).name and abs(value-row['ddg_kcal_mol'])<=1e-12
                           for name,value in _parse_dif_rows(path))]
                if len(matching)!=1:raise ValueError('portable pair requires one matching original DDG table')
                item['ddg_file_path']=str(matching[0]);item['ddg_file_sha256']=sha256_file(matching[0])
            for key in ('wt_path','mut_path','pdb_list_path'):
                path=Path(row[key]); name=f'models/{candidate}/{index:03d}/{path.name}'
                files[name]=path;item[key]=name
            if item.get('ddg_file_path'):
                path=Path(item['ddg_file_path']);name=f'models/{candidate}/{index:03d}/{path.name}'
                files[name]=path;item['ddg_file_path']=name
            # A portable pair must retain the exact input and group mutation
            # records, rather than hashes pointing at an absent batch folder.
            batch=Path(row['batch_path']) if row.get('batch_path') else None
            for filename,key in (('input.pdb','input_path'),('individual_list.txt','mutant_file_path')):
                if batch is not None and (batch/filename).is_file():
                    name=f'models/{candidate}/{index:03d}/{filename}'
                    files[name]=batch/filename;item[key]=name
            item['batch_path']=f'models/{candidate}/{index:03d}'
            portable[candidate].append(item)
    manifest={name:sha256_file(path) for name,path in files.items()}
    portable_bytes=json.dumps(portable,sort_keys=True).encode('utf-8')
    manifest['portable_pairs.json']=hashlib.sha256(portable_bytes).hexdigest()
    destination=output/'surface_perturbation_evidence.zip'
    pending=destination.with_name(destination.name+f'.{os.getpid()}.pending')
    try:
        with ZipFile(pending,'w',ZIP_DEFLATED) as archive:
            for name,path in sorted(files.items()): archive.write(path,name)
            archive.writestr('portable_pairs.json',portable_bytes)
            archive.writestr('evidence_manifest.json',json.dumps({'schema':'1.0','files':manifest},sort_keys=True))
        pending.replace(destination)
    finally:
        pending.unlink(missing_ok=True)
    return destination
