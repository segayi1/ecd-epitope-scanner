"""Independent current-result sequence, mapping, menu and coverage audit.

Candidate sequences are replay records, never qualified recommendations.
Run only against a completed checkpoint, not a concurrently changing output.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from ..joint_design_v4.coverage_profile import sequence_reasons
from ..surface_perturbation.paired import validate_pair
from ..util import sha256_file,write_json


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def verify_hard_curve(summary,masks,area):
    selected=summary['selected_candidate_ids'];curve=summary.get('coverage_curve',[])
    if not curve:raise ValueError('Missing official coverage curve')
    for row in curve:
        chosen=selected[:row['K']]
        covered=np.any([masks[cid] for cid in chosen],axis=0) if chosen else np.zeros(len(area),bool)
        fraction=float(area[covered].sum()/area.sum())
        if len(chosen)!=row['actual_constructs'] or not math.isclose(fraction,row['C_hard_all'],abs_tol=1e-12,rel_tol=0):
            raise ValueError('Official coverage masks differ from completed summary')


def audit(output):
    output=Path(output);sequence=read(output/'uniprot.json')['sequence']['value']
    def csv_rows(name):
        with (output/name).open(encoding='utf-8-sig',newline='') as stream:return list(csv.DictReader(stream))
    mapping={int(r['uniprot_pos']):r for r in csv_rows('residue_mapping.csv')}
    menu=read(output/'mutation_menu_v2.json')['rows'];table=csv_rows('mutation_menu_v2.csv')
    if len(menu)!=len(table):raise ValueError('JSON/CSV replacement membership differs')
    replacements=defaultdict(set)
    for row,text in zip(menu,table):
        residue=row['residue'];position=residue['uniprot_pos'];wt=residue['one_letter']
        if sequence[position-1]!=wt:raise ValueError('replacement table WT mismatch')
        aa=row['mutant_aa'];replacements[position].add(aa)
        if aa==wt:raise ValueError('replacement table includes WT as mutation')
        for key in ('residue','chemical_descriptors','hard_exclusions'):
            if json.loads(text[key])!=row[key]:raise ValueError('replacement CSV differs: '+key)
        for key in ('eligible','menu_selected'):
            if text[key]!=str(row[key]):raise ValueError('replacement CSV eligibility differs')
        value=row['spurs_ddg_kcal_mol']
        if row['eligible'] and (value is None or not math.isfinite(value)):
            raise ValueError('eligible substitution has missing/nonfinite score')
        if value is not None and float(text['spurs_ddg_kcal_mol'])!=value:
            raise ValueError('replacement CSV energy differs')
    if any(len(values)!=19 for values in replacements.values()):raise ValueError('incomplete nineteen-substitution table')
    candidates=read(output/'panel_candidates.json');by_id={c['candidate_id']:c for c in candidates}
    if len(by_id)!=len(candidates):raise ValueError('duplicate candidate identity')
    expected={};sequences=[]
    for candidate in candidates:
        mutations=candidate['mutations'];positions=[m['residue']['uniprot_pos'] for m in mutations]
        if len(positions)!=candidate['order'] or len(set(positions))!=len(positions):
            raise ValueError('candidate order/repeated position mismatch')
        changed=[];normalized=[];mutant=list(sequence)
        for mutation in mutations:
            residue=mutation['residue'];p=residue['uniprot_pos'];letter=mutation['mutant_aa'];source=mapping[p]
            if sequence[p-1]!=residue['one_letter']:raise ValueError('candidate WT mismatch')
            for key in ('chain_id','auth_seq_id','insertion_code'):
                if str(residue[key])!=source[key]:raise ValueError('candidate author mapping mismatch')
            changed.append((p,residue['one_letter'],letter));mutant[p-1]=letter
            normalized.append(f"{residue['one_letter']}{source['normalized_position']}{letter}")
        if sequence_reasons(sequence,changed):raise ValueError('candidate violates full-sequence protection')
        expected[candidate['candidate_id']]=sorted(normalized)
        sequences.append({'candidate_id':candidate['candidate_id'],
            'mutations':','.join(m['notation'] for m in mutations),'full_uniprot_sequence':''.join(mutant),
            'status':'SCREENED_REPLAY_ONLY'})
    physical=output/'surface_perturbation';pairs=read(physical/'paired_manifest.json');identities=set()
    for cid,rows in pairs.items():
        if cid not in expected:raise ValueError('paired candidate absent from current screened result')
        for pair in rows:
            validate_pair(pair)
            if sorted(pair['mutations'])!=expected[cid]:raise ValueError('normalized model differs from author mutation mapping')
            identity=(cid,pair['batch_path'],pair['group'],pair['run'])
            if identity in identities:raise ValueError('repeated actual model identity')
            identities.add(identity)
    with np.load(physical/'canonical_reference.npz',allow_pickle=False) as data:
        area=data['area'];denominator=float(area.sum())
    summary=read(output/'surface_perturbation_preview.json');selected=summary['selected_candidate_ids']
    covered=np.zeros(len(area),bool);hard_masks={}
    for cid in selected:
        if cid not in pairs or len(pairs[cid])<5:raise ValueError('official panel lacks five real model pairs')
        with np.load(physical/(cid+'_coverage.npz'),allow_pickle=False) as masks:
            if masks['hard'].dtype.kind!='b' or masks['hard'].shape!=area.shape:raise ValueError('invalid hard mask')
            covered|=masks['hard']
            hard_masks[cid]=masks['hard'].copy()
    fraction=float(area[covered].sum())/denominator
    verify_hard_curve(summary,hard_masks,area)
    achieved=4<=len(selected)<=12 and fraction>=.8
    if summary['goal_status']=='ACHIEVED':
        if not achieved:raise ValueError('official goal label differs from masks/construct count')
        from ..surface_perturbation.combo_contract import audit_target_panel
        audit_target_panel(output,summary)
    path=output/'candidate_sequence_replay.csv'
    with path.open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['candidate_id','mutations','full_uniprot_sequence','status'])
        writer.writeheader();writer.writerows(sequences)
    return {'engineering_status':'PASS','replacement_rows':len(menu),'positions':len(replacements),
        'screened_candidates':len(candidates),'verified_actual_pair_identities':len(identities),
        'reliable_wt_area_a2':denominator,'official_constructs':len(selected),'official_hard_fraction':fraction,
        'scientific_goal_status':summary['goal_status'],'search_status':summary.get('search_status'),
        'apbs_status':summary.get('final_apbs_evidence',{'status':'NOT_REPORTED'}),
        'replay_sequences_sha256':sha256_file(path),
        'scope':'Sequence, author-to-normalized model mapping, complete menu, actual pairs/DDG and official masks; no promotion of shadow evidence'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    parser.add_argument('--report',type=Path,required=True);args=parser.parse_args()
    result=audit(args.output);write_json(args.report,result);print(json.dumps(result))
if __name__=='__main__':main()
