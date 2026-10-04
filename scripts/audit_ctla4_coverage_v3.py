"""Independent saved v3 checkpoint/package verification; no scientific solves."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from ..util import write_json,sha256_file
from ..surface_perturbation.paired import validate_pair
from ..surface_perturbation.sampling_v3 import validate_roles
from ..joint_design_v4.coverage_profile import sequence_reasons


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def archive(path,manifest_name):
    with ZipFile(path) as z:
        manifest=json.loads(z.read(manifest_name))['files']
        for name,digest in manifest.items():
            actual=hashlib.sha256()
            with z.open(name) as stream:
                while block:=stream.read(1024**2):actual.update(block)
            if actual.hexdigest()!=digest:raise ValueError('Archive hash mismatch: '+name)
        if 'portable_pairs.json' in z.namelist():
            identities=set()
            for rows in json.loads(z.read('portable_pairs.json')).values():
                for p in rows:
                    source=p.get('source_execution_identity',{})
                    identity=tuple(source.get(key) for key in ('batch_path','group','run'))
                    if not all(v is not None for v in identity) or identity in identities:
                        raise ValueError('Portable source execution identity missing or duplicated')
                    identities.add(identity)
                    for key,digest in (('wt_path','wt_sha256'),('mut_path','mut_sha256'),
                        ('pdb_list_path','pdb_list_sha256'),('ddg_file_path','ddg_file_sha256'),
                        ('input_path','input_sha256'),('mutant_file_path','mutant_file_sha256')):
                        if key in p and hashlib.sha256(z.read(p[key])).hexdigest()!=p[digest]:
                            raise ValueError('Portable model identity mismatch')
        return {'files':len(manifest),'bytes':path.stat().st_size,'sha256':sha256_file(path)}


def audit(route,packages=False):
    output=route/'results';folder=output/'surface_perturbation'
    summary=read(output/'surface_perturbation_preview.json')
    if summary['panel_contract']!='positive_lower_bound_information_v3':raise ValueError('Wrong evidence contract')
    pairs=read(folder/'paired_manifest.json');history=read(route/'screened_candidate_history.json')
    candidates={cid:r['candidate'] for cid,r in history.items()}
    candidates.update({c['candidate_id']:c for c in read(output/'panel_candidates.json')})
    sequence=read(output/'uniprot.json')['sequence']['value'];identities=set();counts={};diversity={}
    for cid,rows in pairs.items():
        c=candidates[cid]
        if sequence_reasons(sequence,[(m['residue']['uniprot_pos'],m['residue']['one_letter'],m['mutant_aa']) for m in c['mutations']]):
            raise ValueError('Full-sequence mutation protection failed')
        for p in rows:
            validate_pair(p)
            identity=(p['batch_path'],p['group'],p['run'])
            if identity in identities:raise ValueError('Execution shared between different genotypes/roles')
            identities.add(identity)
            counts[p['sampling_role']]=counts.get(p['sampling_role'],0)+1
        diversity[cid]={role:{'WT_content':len({p['wt_sha256'] for p in rows if p['sampling_role']==role}),
                        'mutant_content':len({p['mut_sha256'] for p in rows if p['sampling_role']==role}),
                        'logical_runs':sum(p['sampling_role']==role for p in rows)}
                        for role in ('production','null_train','null_holdout')}
    with np.load(folder/'coverage_matrix.npz',allow_pickle=False) as data:
        lower,upper,area,direct,ids=[data[k].copy() for k in ('hard','upper','area','direct','candidate_ids')]
    if np.any(lower&~upper) or not np.isclose(area.sum(),6874.363277945836,atol=1e-8,rtol=0):raise ValueError('Coverage contract mismatch')
    if ids.tolist()!=summary['selected_candidate_ids']:raise ValueError('Selected IDs/matrix differ')
    union=lower.any(axis=0);screen=float(area[union].sum()/area.sum());info=summary['panel_information']
    double=float(area[lower.sum(axis=0)>=2].sum()/area.sum())
    worst=min((float(area[lower[np.arange(len(lower))!=i].any(axis=0)].sum()/area.sum()) for i in range(len(lower))),default=0.)
    for value,key in ((screen,'C_screen_lower_bound'),(double,'double_covered_fraction'),(worst,'worst_single_construct_loss_coverage'),
                      (float(area[direct.any(axis=0)].sum()/area.sum()),'C_physical_lower_bound')):
        if not np.isclose(value,info[key],atol=1e-12,rtol=0):raise ValueError('Panel metric mismatch: '+key)
    # Independent integer-code enumeration; unknown contributes to neither bit.
    codes=np.stack((lower.T,~upper.T),axis=2)
    groups,inverse=np.unique(codes.reshape(len(area),-1),axis=0,return_inverse=True)
    weights=np.bincount(inverse,weights=area);numerator=np.zeros(3)
    for i,j in itertools.combinations(range(len(groups)),2):
        a=groups[i].reshape(-1,2);b=groups[j].reshape(-1,2)
        distance=np.count_nonzero((a[:,0]&b[:,1])|(a[:,1]&b[:,0]))
        numerator[:min(3,distance)]+=weights[i]*weights[j]
    denominator=(area.sum()**2-np.dot(area,area))/2
    if any(not np.isclose(numerator[k]/denominator,info['geometric_pair_distinction'][str(k+1)],atol=1e-12) for k in range(3)):
        raise ValueError('Ternary geometric code metric mismatch')
    for cid in ids:
        rows=pairs[cid];production=validate_roles(rows,'production',5,'production_v3')
        if len(production)!=5 or np.median([p['ddg_kcal_mol'] for p in production])>2:raise ValueError('Formal production stability failed')
        scope=summary['candidate_local_calibration'][cid]
        if not scope.get('positive_lower_bound_qualified') or not scope.get('holdout_validated'):raise ValueError('Formal lower-bound qualification absent')
        cohorts={p['sampling_cohort'] for p in rows if p['sampling_role']=='null_train'}
        if len(cohorts)!=1:raise ValueError('Null cohort inconsistent')
        cohort=next(iter(cohorts))
        if any(len(validate_roles(rows,role,10,cohort))!=10 for role in ('null_train','null_holdout')):raise ValueError('Null split missing')
        c=candidates[cid]
        if c['order']==3:
            required={frozenset(m['notation'] for m in options) for options in itertools.combinations(c['mutations'],2)}
            available={frozenset(m['notation'] for m in d['mutations']):key for key,d in candidates.items() if d['order']==2}
            if len(required)!=3 or not required<=available.keys():raise ValueError('Exact triple controls missing')
            gains=[]
            with np.load(folder/(cid+'_coverage.npz'),allow_pickle=False) as data:triple_lower=data['hard']
            for signature in required:
                control=available[signature]
                controls=validate_roles(pairs[control],'production',5,'production_v3')
                if len(controls)!=5 or np.median([p['ddg_kcal_mol'] for p in controls])>2:raise ValueError('Triple control stability missing')
                with np.load(folder/(control+'_coverage.npz'),allow_pickle=False) as data:
                    gains.append(float(area[triple_lower&~data['upper']].sum()))
            if c['triple_role']=='TRIPLE_EXPANSION' and min(gains)+1e-12<.05*area.sum():raise ValueError('Triple gain vs unknown upper not demonstrated')
    report={'engineering_status':'PASS','hard_target_status':'ACHIEVED' if 4<=len(ids)<=12 and screen>=.8 else 'BELOW_TARGET',
            'panel_information':info,'actual_pair_identities':len(identities),'role_counts':counts,'content_diversity':diversity,
            'new_scientific_calls':0,'scientific_summary_sha256':sha256_file(output/'surface_perturbation_preview.json')}
    if report['hard_target_status']!=summary['hard_target_status']:raise ValueError('Hard target label mismatch')
    if packages:report['packages']={name:archive(output/name,manifest) for name,manifest in
        (('surface_perturbation_evidence.zip','evidence_manifest.json'),('ecd_epitope_design.zip','run_manifest.json'))}
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('route',type=Path)
    parser.add_argument('--packages',action='store_true');parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();result=audit(args.route,args.packages);write_json(args.report,result)
    print(json.dumps({k:v for k,v in result.items() if k!='content_diversity'}))
