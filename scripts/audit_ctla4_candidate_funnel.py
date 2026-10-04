"""Read-only diagnosis of the cumulative screened pool and current real evidence."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from ecd_epitope_scanner.util import sha256_file, write_json


def audit(case):
    results=case/'results'
    sources={'history':case/'screened_candidate_history.json',
        'candidates':results/'panel_candidates.json',
        'production':results/'surface_perturbation_preview.json',
        'pairs':results/'surface_perturbation/paired_manifest.json'}
    def read(path):return json.loads(path.read_text(encoding='utf-8'))
    history=read(sources['history']) if sources['history'].is_file() else {}
    candidates={cid:row['candidate'] for cid,row in history.items()}
    candidates.update({row['candidate_id']:row for row in read(sources['candidates'])})
    summary=read(sources['production']);pairs=read(sources['pairs'])
    formal=set(summary.get('independent_K_verified_candidate_ids',[]))
    shadow=set(summary.get('shadow_positive_lower_bound',{}).get('candidate_ids',[]))
    rows=[]
    for cid,candidate in sorted(candidates.items()):
        site=summary.get('candidate_site_evidence',{}).get(cid,history.get(cid,{}).get('candidate_site_evidence',{}))
        local=summary.get('candidate_local_calibration',{}).get(cid,history.get(cid,{}).get('candidate_local_calibration',{}))
        runs=len(pairs.get(cid,[]))
        rows.append({'candidate_id':cid,'mutations':';'.join(m['notation'] for m in candidate['mutations']),
            'candidate_source':'STABILITY_FAILURE_BACKUP' if history.get(cid,{}).get('conditional_backup_origin') else candidate.get('search_role'),
            'order':candidate['order'],'stability_status':candidate.get('stability',{}).get('consensus_status','UNKNOWN'),
            'upstream_formal_eligible':candidate.get('formal_eligible'),
            'upstream_exclusion_reasons':json.dumps(candidate.get('formal_exclusion_reasons',[])),
            'physical_pair_identities':runs,
            'evidence_tier':'FIVE_OR_MORE_PAIRED_MODELS' if runs>=5 else 'UPSTREAM_PRESCREEN_ONLY',
            'same_four_site_changes_complete':site.get('complete'),
            'site_reason':site.get('reason'),
            'full_local_calibration_complete':local.get('complete'),
            'measurement_local_complete':local.get('measurement_local_complete'),
            'hard_scan_bound_resolved':local.get('hard_scan_bound_resolved'),
            'unresolved_hard_upper_area_a2':local.get('hard_scan_unknown_upper_area_a2'),
            'formally_verified_pool':cid in formal,'shadow_positive':cid in shadow})
    destination=case/'candidate_failure_diagnosis.csv'
    with destination.open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else ['candidate_id'])
        writer.writeheader();writer.writerows(rows)
    report={'engineering_status':'PASS','scope':'Cumulative saved prescreens plus current physical evidence; '
        'no new solvers or scientific requalification. Missing evidence remains null.',
        'sources_sha256':{name:sha256_file(path) for name,path in sources.items() if path.is_file()},
        'screened_candidates':len(rows),'order_counts':dict(Counter(row['order'] for row in rows)),
        'upstream_stability_status_counts':dict(Counter(row['stability_status'] for row in rows)),
        'five_pair_groups':sum(row['physical_pair_identities']>=5 for row in rows),
        'current_physical_pair_identities':sum(map(len,pairs.values())),
        'same_four_site_positive':sum(row['same_four_site_changes_complete'] is True for row in rows),
        'formal_verified_pool':len(formal),'shadow_positive':len(shadow),
        'positive_but_unverified':[row for row in rows if row['shadow_positive'] and not row['formally_verified_pool']],
        'formal_goal_status':summary.get('goal_status'),'search_status':summary.get('search_status'),
        'interpretation':'Inactive prescreen rejection is not a five-model stability determination; '
        'shadow positives do not establish formal coverage. Pair totals are identities, not new FoldX call counts.'}
    index=case/'evidence_index/surface_perturbation/paired_manifest.json'
    if index.is_file():
        cumulative=read(index)
        report['cumulative_physical_pair_identities']=sum(map(len,cumulative.values()))
        report['cumulative_physical_candidate_groups']=len(cumulative)
        report['sources_sha256']['cumulative_pair_index']=sha256_file(index)
    write_json(case/'candidate_failure_diagnosis.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(audit(args.case_root),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
