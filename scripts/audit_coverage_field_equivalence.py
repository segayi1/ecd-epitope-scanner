"""Read-only exact array audit across conservative source cache generations.

This does not migrate caches or certify differing calibration/classification
protocols. It tests the promised fixed-input numerical and missing-data invariant.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def same(left,right):
    return left.dtype==right.dtype and np.array_equal(left,right,equal_nan=left.dtype.kind in 'fc')


def audit(case):
    folder=Path(case)/'results/surface_perturbation/fields';groups={}
    for path in sorted(folder.glob('field_*.json')):
        metadata=json.loads(path.read_text(encoding='utf-8'));contract=metadata['contract']
        if hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest()!=metadata.get('sha256'):
            raise ValueError('Raw field content identity mismatch: '+str(path))
        key=tuple((k,json.dumps(v,sort_keys=True)) for k,v in sorted(contract.items())
            if k not in {'implementation_sha256','reference'})
        groups.setdefault(key,[]).append((path,metadata))
    # All references must have exactly the same physical arrays. A changed
    # designability mask is allowed; a changed surface is never compared.
    reference_paths=set()
    for path in Path(case).glob('*_object_manifest.json'):
        metadata=json.loads(path.read_text(encoding='utf-8'))
        row=metadata['files'].get('results/surface_perturbation/canonical_reference.npz')
        if row:reference_paths.add(Path(case)/row['object'])
    current=Path(case)/'results/surface_perturbation/canonical_reference.npz'
    with np.load(current,allow_pickle=False) as source:
        reference={k:source[k] for k in source.files if k!='designable'}
    for path in reference_paths:
        with np.load(path,allow_pickle=False) as other:
            for key,value in reference.items():
                if key not in other or not same(value,other[key]):
                    raise ValueError('Physical reference changed; comparison forbidden: '+key)
    checked=[]
    for rows in groups.values():
        if len(rows)<2:continue
        with np.load(rows[0][0].with_suffix('.npz'),allow_pickle=False) as data:
            baseline={k:data[k] for k in data.files}
        for path,metadata in rows[1:]:
            with np.load(path.with_suffix('.npz'),allow_pickle=False) as data:
                if set(data.files)!=set(baseline):raise ValueError('Field array membership changed')
                for key,value in baseline.items():
                    other=data[key]
                    if not same(value,other):
                        raise ValueError('Raw field/NaN/mask numerical difference: '+key)
            checked.append({'first':rows[0][0].name,'second':path.name,'arrays':len(baseline),
                'wt_sha256':metadata['contract']['wt'],'mut_sha256':metadata['contract']['mut']})
    return {'status':'PASS' if checked else 'NOT_YET_COMPARABLE',
        'reference_objects_checked':len(reference_paths),'comparisons':checked,
        'scope':'Exact raw arrays, dtype, NaN positions and valid masks. Changed calibration or repeat counts are not classification-equivalence claims.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    report=audit(args.case)
    args.output.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'status':report['status'],'comparisons':len(report['comparisons'])}))
if __name__=='__main__':main()
