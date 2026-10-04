"""Evaluate frozen hidden contacts after a verified panel; no design call."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.surface_perturbation.sanity import contact_sanity
from ecd_epitope_scanner.util import write_json
from ecd_epitope_scanner.surface_perturbation.case_status import case_status

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('case_output',type=Path);p.add_argument('labels',type=Path)
    p.add_argument('--label-sha256',required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    summary=json.loads((args.case_output/'surface_perturbation_preview.json').read_text(encoding='utf-8'))
    status=case_status(args.case_output,summary)
    if status not in {'VERIFIED_PANEL','TARGET_PANEL_VERIFIED_SEARCH_OPEN'}:
        write_json(args.output,{'status':'UNAVAILABLE','reason':'no accepted verified actual panel','case_status':status});return
    protein=json.loads((args.case_output/'analysis.json').read_text(encoding='utf-8'))['protein']
    labels=json.loads(args.labels.read_text(encoding='utf-8'))
    if labels['accession']!=protein['accession']:
        raise ValueError('contact label accession differs from frozen case')
    with np.load(args.case_output/'surface_perturbation'/'coverage_matrix.npz',allow_pickle=False) as matrix:
        reference={'area':matrix['area']}
        hard=matrix['hard'].copy();ids=matrix['candidate_ids'].tolist()
    with np.load(args.case_output/'surface_perturbation'/'canonical_reference.npz',allow_pickle=False) as canonical:
        reference['owner']=canonical['owner'].copy()
    write_json(args.output,contact_sanity(args.labels,args.label_sha256,reference,hard,ids))
if __name__=='__main__':main()
