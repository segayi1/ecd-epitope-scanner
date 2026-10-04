"""Saved-input NMR and geodesic sampling audits; zero external solves."""
import argparse
from dataclasses import fields
import itertools
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

from ecd_epitope_scanner.models import ResidueRecord
from ecd_epitope_scanner.surface_perturbation.geometry import Mesh
from ecd_epitope_scanner.surface_perturbation.paired import atoms
from ecd_epitope_scanner.surface_perturbation.reference_audit import sampling_audit
from ecd_epitope_scanner.util import write_json,sha256_file
from ecd_epitope_scanner.scripts.run_ctla4_coverage_v2 import read


def audit(case,assets):
    path=assets/'workspace_data/8c068ad701a6/structure.cif'
    structure=MMCIFParser(QUIET=True).get_structure('1AH1',str(path))
    maps=[];model_ids=[]
    for model in structure:
        chain=model['A']
        maps.append({(r.id[1],r.id[2],r.resname):r['CA'].coord.astype(float)
                     for r in chain if r.id[0]==' ' and 'CA' in r})
        model_ids.append(model.serial_num)
    common=sorted(set.intersection(*(set(m) for m in maps)))
    assert len(maps)==20 and len(common)>20
    coordinates=[np.array([m[k] for k in common]) for m in maps]
    matrix=np.zeros((20,20))
    for i,j in itertools.combinations(range(20),2):
        x=coordinates[i]-coordinates[i].mean(axis=0);y=coordinates[j]-coordinates[j].mean(axis=0)
        u,_,vt=np.linalg.svd(x.T@y);correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
        rotation=u@correction@vt
        matrix[i,j]=matrix[j,i]=np.sqrt(np.mean(np.sum((x@rotation-y)**2,axis=1)))
    representative=min(itertools.combinations(range(20),3),key=lambda ids:(float(matrix[:,ids].min(axis=1).mean()),ids))
    output=case/'routes/coverage_v3';output.mkdir(parents=True,exist_ok=True)
    write_json(output/'T0_NMR_ensemble_audit.json',{'engineering_status':'PASS','input_sha256':sha256_file(path),
        'chain':'A','models':20,'common_CA_residues':len(common),'model_ids':model_ids,
        'aligned_pairwise_CA_RMSD_a':matrix.tolist(),
        'representative_model_ids':[model_ids[i] for i in representative],
        'new_external_solver_calls':0,'formal_ensemble_use':False,
        'interpretation':'submitted NMR conformers; not independent equilibrium samples or new production repeats'})
    folder=case/'results';names={f.name for f in fields(ResidueRecord)}
    records=[ResidueRecord(**{k:v for k,v in r.items() if k in names}) for r in read(folder/'analysis.json')['residues']]
    native=atoms(folder/'structure_for_energy.pdb')
    with np.load(folder/'surface_perturbation/native_WT_SES_mesh.npz',allow_pickle=False) as data:
        mesh=Mesh(data['vertices'],data['faces'])
    reliable={r.uniprot_pos for r in records if r.coordinate_eligible and not r.review_reasons}
    designable={r['residue']['uniprot_pos'] for r in read(folder/'mutation_menu_v2.json')['rows'] if r['eligible']}
    # Independent centers on the saved same mesh; no claim of finer SES.
    result=sampling_audit(mesh,mesh,native,records,reliable,designable,maximum_anchors=6)
    result['settings'][2]['setting']='same_mesh_control'
    result['new_external_solver_calls']=0
    assert abs(result['settings'][0]['target_area_a2']-6874.363277945836)<1e-8
    write_json(output/'T2_geodesic_sampling_audit.json',result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root',type=Path,required=True)
    parser.add_argument('--assets-root',type=Path,required=True)
    args=parser.parse_args();audit(args.case_root,args.assets_root)
