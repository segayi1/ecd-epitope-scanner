"""Replay saved raw fields and nulls without running a scientific worker."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.sparse import load_npz
from ..surface_perturbation.fields import IMPLEMENTATION_ID,reference_id
from ..surface_perturbation.paired import atoms
from ..surface_perturbation.null_store import protocol_identity
from ..surface_perturbation.scoped_fields import protocol_identity as scoped_protocol,SCOPE_IMPLEMENTATION
from ..surface_perturbation.calibration_jobs import usable_calibration as load_calibration
from ..surface_perturbation.classification import classify_fields
from ..surface_perturbation.coverage import budget_curve
from ..util import sha256_file,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case_output',type=Path)
    p.add_argument('--calibration',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();folder=args.case_output/'surface_perturbation'
    summary=json.loads((args.case_output/'surface_perturbation_preview.json').read_text(encoding='utf-8'))
    identities=summary.get('apbs_audit_target_candidate_ids') or summary.get('selected_candidate_ids')
    if not identities:identities=summary['trajectories'][-1]['verified']
    records=[SimpleNamespace(**r) for r in json.loads((args.case_output/'analysis.json').read_text(encoding='utf-8'))['residues']]
    with np.load(folder/'canonical_reference.npz') as data:reference={k:data[k].copy() for k in data.files}
    reference['graph']=load_npz(folder/'reference_graph.npz')
    reference['context_identity']=json.loads((folder/'reference_contract.json').read_text(encoding='utf-8'))['context_identity']
    pairs=json.loads((folder/'paired_manifest.json').read_text(encoding='utf-8'));native=atoms(args.case_output/'structure_for_energy.pdb')
    catalog={}
    for path in list((folder/'fields').glob('field_*.json'))+list((folder/'fields').glob('scoped_*.json')):
        meta=json.loads(path.read_text(encoding='utf-8'));contract=meta['contract']
        if contract['implementation_sha256']!=IMPLEMENTATION_ID or contract['reference']!=reference_id(reference):continue
        if 'query_scope_sha256' in contract and contract.get('scope_implementation')!=SCOPE_IMPLEMENTATION:continue
        if meta['sha256']!=sha256_file(path.with_suffix('.npz')):raise ValueError('raw field content changed')
        catalog[(contract['wt'],contract['mut'],contract.get('query_scope_sha256'))]=(path.with_suffix('.npz'),contract['nano'])
    report={'scope':'saved-field classification equivalence and ablation only; not current formal qualification or a completed panel','candidates':{}}
    grouped={name:[] for name in ('atom_occupancy_shape','shape_chemistry','shape_chemistry_ses')}
    for cid in identities:
        rows=pairs[cid][:5];runs=[]
        audit=json.loads((folder/f'{cid}_field_audit.json').read_text(encoding='utf-8'))
        scope=audit.get('query_scope_sha256')
        for row in rows:
            path,nano=catalog[(row['wt_sha256'],row['mut_sha256'],scope)]
            with np.load(path) as values:runs.append({k:values[k].copy() for k in values.files})
        tools=SimpleNamespace(nano='saved_binary_identity',identity=lambda executable:nano)
        protocol=(scoped_protocol if scope is not None else protocol_identity)(rows,reference,tools)
        calibration=load_calibration(args.calibration,protocol)
        classified=classify_fields(runs,rows,native,records,reference,calibration,protocol)
        for name in grouped:grouped[name].append(classified['groups'][name])
        with np.load(folder/f'{cid}_coverage.npz') as previous:
            equal={name:bool(np.array_equal(classified[name],previous[name])) for name in ('direct','hard','two_direct_sites')}
            equal['unknown']=bool(np.array_equal(~classified['known'],previous['unknown']))
        report['candidates'][cid]={'arrays_equal':equal,'status':'PASS' if all(equal.values()) else 'DIFFERENT'}
    report['status']='PASS' if all(row['status']=='PASS' for row in report['candidates'].values()) else 'FAIL'
    report['nonfinal_saved_field_ablations']={name:budget_curve(np.asarray(masks),np.full(np.asarray(masks).shape,np.nan),
        list(range(len(identities))),reference['area'],reference['designable'],reference['graph']) for name,masks in grouped.items()}
    write_json(args.output,report);print(json.dumps(report))
    if report['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
