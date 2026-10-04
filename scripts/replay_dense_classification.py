"""Reclassify saved dense fields; no modeling, surface or electrostatic solve."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.sparse import load_npz
from ..surface_perturbation.fields import IMPLEMENTATION_ID,reference_id
from ..surface_perturbation.paired import atoms,validate_pair
from ..surface_perturbation.null_store import protocol_identity,load_calibration
from ..surface_perturbation.classification import classify_fields
from ..util import sha256_file,write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case_output',type=Path);parser.add_argument('dense_output',type=Path)
    parser.add_argument('--candidate',required=True);args=parser.parse_args();folder=args.dense_output
    with np.load(folder/'canonical_reference.npz') as data:reference={k:data[k].copy() for k in data.files}
    reference['graph']=load_npz(folder/'reference_graph.npz')
    reference['context_identity']=json.loads((folder/'reference_contract.json').read_text(encoding='utf-8'))['context_identity']
    pairs=json.loads((folder/'paired_manifest.json').read_text(encoding='utf-8'))[args.candidate]
    catalog={}
    for path in (folder/'fields').glob('field_*.json'):
        meta=json.loads(path.read_text(encoding='utf-8'));contract=meta['contract']
        if contract['implementation_sha256']!=IMPLEMENTATION_ID or contract['reference']!=reference_id(reference):continue
        if meta['sha256']!=sha256_file(path.with_suffix('.npz')):raise ValueError('saved raw field fingerprint differs')
        catalog[(contract['wt'],contract['mut'])]=(path.with_suffix('.npz'),contract['nano'])
    runs=[]
    for pair in pairs:
        validate_pair(pair);path,nano=catalog[(pair['wt_sha256'],pair['mut_sha256'])]
        with np.load(path) as data:runs.append({k:data[k].copy() for k in data.files})
    tools=SimpleNamespace(nano='saved_identity',identity=lambda executable:nano)
    protocol=protocol_identity(pairs,reference,tools);calibration=load_calibration(folder/'calibration',protocol)
    records=[SimpleNamespace(**r) for r in json.loads((args.case_output/'analysis.json').read_text(encoding='utf-8'))['residues']]
    classified=classify_fields(runs,pairs,atoms(args.case_output/'structure_for_energy.pdb'),records,reference,calibration,protocol)
    with np.load(folder/'dense_coverage.npz') as old:
        equal={name:bool(np.array_equal(classified[name],old[name])) for name in ('direct','hard','two_direct_sites')}
        equal['unknown']=bool(np.array_equal(~classified['known'],old['unknown']))
    area=reference['area'];coarse=args.case_output/'surface_perturbation'
    with np.load(coarse/'canonical_reference.npz') as ref,np.load(coarse/f'{args.candidate}_coverage.npz') as mask:
        coarse_fraction=float(ref['area'][mask['hard']].sum()/ref['area'].sum())
    dense_fraction=float(area[classified['hard']].sum()/area.sum())
    root=Path(__file__).resolve().parents[1]
    report={'status':'PASS_UNCHANGED' if all(equal.values()) else 'CLASSIFICATION_CHANGED',
        'arrays_equal':equal,'scope':'one saved candidate, dense repeat classification; not full-panel convergence',
        'candidate_id':args.candidate,'coarse_hard_center_fraction':coarse_fraction,
        'dense_hard_center_fraction':dense_fraction,'absolute_fraction_difference':abs(dense_fraction-coarse_fraction),
        'no_scientific_worker_calls':True,
        'loaded_source_sha256':{name:sha256_file(root/name) for name in
            ('surface_perturbation/classification.py','surface_perturbation/engine.py','surface_perturbation/directions.py','surface_perturbation/coverage.py')},
        'original_report_sha256':sha256_file(folder/'dense_audit.json')}
    write_json(folder/'classification_replay.json',report);print(json.dumps(report))


if __name__=='__main__':main()
