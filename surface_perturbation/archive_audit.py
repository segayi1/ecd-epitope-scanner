"""Verify portable evidence as an independent extracted bundle."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path,PurePosixPath
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from .paired import validate_pair


def safe_name(name):
    path=PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or '\\' in name or ':' in name:
        raise ValueError('unsafe archive path')
    return path


def validate_archive(path,*,verify_pairs=True,temporary_root=None,
                     manifest_name='evidence_manifest.json',pairs_name='portable_pairs.json'):
    with ZipFile(path) as archive,TemporaryDirectory(prefix='surface_evidence_audit_',dir=temporary_root) as directory:
        names=archive.namelist()
        if len(names)!=len(set(names)):raise ValueError('duplicate archive member')
        for name in names:safe_name(name)
        safe_name(manifest_name);safe_name(pairs_name)
        manifest=json.loads(archive.read(manifest_name))
        expected=manifest['files']
        if set(names)!=set(expected)|{manifest_name}:raise ValueError('archive manifest membership differs')
        root=Path(directory);total=0
        for name,digest in expected.items():
            safe_name(name);destination=root/name;destination.parent.mkdir(parents=True,exist_ok=True)
            actual=hashlib.sha256()
            with archive.open(name) as source,destination.open('wb') as output:
                while chunk:=source.read(1024*1024):actual.update(chunk);output.write(chunk);total+=len(chunk)
            if actual.hexdigest()!=digest:raise ValueError('archive content hash changed: '+name)
        portable=json.loads((root/pairs_name).read_text(encoding='utf-8'))
        count=0
        if verify_pairs:
            for cid,rows in portable.items():
                for row in rows:
                    item=dict(row)
                    for key in ('batch_path','wt_path','mut_path','pdb_list_path','ddg_file_path','input_path','mutant_file_path'):
                        if item.get(key):safe_name(item[key]);item[key]=str(root/item[key])
                    try:validate_pair(item)
                    except ValueError as exc:raise ValueError(f'{cid} group={row.get("group")} run={row.get("run")}: {exc}') from exc
                    count+=1
        # Validate the inner field/null identities, not merely an outer manifest
        # that could consistently describe a corrupted cache generation.
        from ..util import sha256_file
        from .null_store import load_calibration
        import numpy as np
        null_count=0;field_count=0
        for metadata_path in (root/'fields').glob('*.json'):
            metadata=json.loads(metadata_path.read_text(encoding='utf-8'))
            digest=metadata.get('sha256')
            if digest is not None:
                arrays=metadata_path.with_suffix('.npz')
                if not arrays.is_file() or sha256_file(arrays)!=digest:
                    raise ValueError('archived field inner identity mismatch: '+metadata_path.name)
                with np.load(arrays,allow_pickle=False) as data:
                    for key in data.files:data[key]
                field_count+=1
        for metadata_path in (root/'calibration').rglob('*.json'):
            metadata=json.loads(metadata_path.read_text(encoding='utf-8'))
            if 'array_sha256' in metadata:
                channels=load_calibration(metadata_path.parent,metadata_path.stem)
                for channel in channels.values():
                    if 'null_blocks' in channel and len(channel['null_blocks'])!=channel['block_count']:
                        raise ValueError('archived null block count mismatch')
                    if 'required_mask' in channel and channel['required_mask'].dtype.kind!='b':
                        raise ValueError('archived null scope is not a boolean mask')
                null_count+=1
        summary=root/'surface_perturbation_preview.json'
        status=json.loads(summary.read_text(encoding='utf-8'))['status'] if summary.exists() else 'LEGACY_BUNDLE_STATUS_NOT_EMBEDDED'
        return {'status':'PASS','case_status':status,'members':len(names),'verified_pairs':count,
            'verified_field_objects':field_count,'verified_null_cohorts':null_count,
            'uncompressed_bytes':total,'scope':'portable hashes, actual paired models and DDG, field/null inner identities; not scientific coverage acceptance'}
