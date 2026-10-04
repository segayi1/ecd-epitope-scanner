"""Audit saved canonical SES on dense centers and a finer SES sampling."""
from __future__ import annotations
import argparse,csv,json,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.surface_perturbation.reference_audit import sampling_audit
from ecd_epitope_scanner.surface_perturbation.paired import atoms
from ecd_epitope_scanner.surface_perturbation.tools import PhysicalTools
from ecd_epitope_scanner.util import write_json

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case_output',type=Path)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();started=time.perf_counter();folder=args.case_output/'surface_perturbation'
    with np.load(folder/'canonical_reference.npz',allow_pickle=False) as ref:
        reliable=set(map(int,ref['owner']));designable=set(map(int,ref['owner'][ref['designable']]))
        domains={int(pos):str(domain) for pos,domain in zip(ref['owner'],ref['domain'])}
    with (args.case_output/'residue_mapping.csv').open(encoding='utf-8-sig',newline='') as f:
        records=[SimpleNamespace(chain_id=r['chain_id'],auth_seq_id=int(r['auth_seq_id']),
            insertion_code=r['insertion_code'],uniprot_pos=int(r['uniprot_pos']),
            domain=domains.get(int(r['uniprot_pos']),'UNMAPPED')) for r in csv.DictReader(f)]
    native=atoms(args.case_output/'structure_for_energy.pdb');tools=PhysicalTools(args.cache)
    coarse=tools.ses(native,scale=2.);fine=tools.ses(native,scale=3.)
    result=sampling_audit(coarse,fine,native,records,reliable,designable)
    result.update(elapsed_seconds=time.perf_counter()-started,tool_calls=tools.calls)
    write_json(args.output,result);print(json.dumps(result))
if __name__=='__main__':main()
