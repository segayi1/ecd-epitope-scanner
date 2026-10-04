"""Same-eligibility, same-graph, same-cap menu search ablation; no new models."""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey,V4MutationOption
from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
from ecd_epitope_scanner.surface_perturbation.runtime import RunBudget
from ecd_epitope_scanner.util import sha256_file,write_json


def signature(candidate):
    return tuple(sorted(m.notation for m in candidate.mutations))


def audit(case):
    result=case/'results'
    rows=json.loads((result/'mutation_menu_v2.json').read_text(encoding='utf-8'))['rows']
    eligible=defaultdict(list)
    for row in rows:
        if row['eligible'] and not row['review_required']:eligible[row['residue']['uniprot_pos']].append(row)
    old=[row for group in eligible.values() for row in sorted(group,key=lambda r:(r['spurs_ddg_kcal_mol'],r['mutant_aa']))[:3]]
    new=[row for group in eligible.values() for row in group if row['menu_selected']]
    graph={}
    for e in json.loads((result/'surface_graph_v4.json').read_text())['edges']:
        graph.setdefault(e['source'],{})[e['target']]=e['distance_a']
        graph.setdefault(e['target'],{})[e['source']]=e['distance_a']
    sequence=json.loads((result/'uniprot.json').read_text(encoding='utf-8'))['sequence']['value']
    cap=json.loads((result/'candidate_pool_v2.json').read_text())['diagnostics']['hard_pool_capacity']
    variants={};diagnostics={}
    for name,menu in (('energy_three',old),('mechanism_menu_with_logged_supplements',new)):
        options=[V4MutationOption(ResidueKey(**row['residue']),row['mutant_aa'],
            row['residue']['one_letter']+str(row['residue']['uniprot_pos'])+row['mutant_aa'],
            row['chemistry_class'],row['spurs_ddg_kcal_mol'],0,True) for row in menu]
        diagnostic={}
        pool=generate_candidates([],options,graph,{'coverage_scan_v2':True,
            'canonical_sequence':sequence,'global_shortlist':cap},diagnostics=diagnostic)
        variants[name]={signature(c) for c in pool};diagnostics[name]=diagnostic
    candidates=json.loads((result/'panel_candidates.json').read_text())
    summary=json.loads((result/'surface_perturbation_preview.json').read_text())
    shadow=set(summary.get('shadow_positive_lower_bound',{}).get('candidate_ids',[]))
    measured=[]
    for c in candidates:
        key=tuple(sorted(m['notation'] for m in c['mutations']))
        measured.append({'candidate_id':c['candidate_id'],'mutation_notations':list(key),
            'stability_status':(c.get('stability') or {}).get('consensus_status'),
            'same_four_site_changes_complete':summary.get('candidate_site_evidence',{}).get(c['candidate_id'],{}).get('complete',False),
            'shadow_positive_lower_bound_candidate':c['candidate_id'] in shadow,
            'in_pools':{name:key in signatures for name,signatures in variants.items()}})
    return {'status':'PASS','new_foldx_calls':0,'new_spurs_calls':0,'pool_cap':cap,
        'menu_table_sha256':sha256_file(result/'mutation_menu_v2.json'),
        'graph_sha256':sha256_file(result/'surface_graph_v4.json'),'diagnostics':diagnostics,
        'measured_candidates':measured,
        'interpretation':'Isolates menu changes on identical eligibility/graph/cap. Includes only already measured outcomes; not an old-route rerun, physics speed benchmark, accuracy claim or scientific acceptance.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    from ecd_epitope_scanner.surface_perturbation.runtime import read_budget
    budget=RunBudget(args.case/'case_budget.json',read_budget(args.case/'case_budget.json').get('limit_seconds',14400))
    try:
        report=audit(args.case);write_json(args.output,report)
        print(json.dumps({'status':report['status'],'shadow_candidates':[r for r in report['measured_candidates'] if r['shadow_positive_lower_bound_candidate']]}))
    finally:budget.save(False)
if __name__=='__main__':main()
