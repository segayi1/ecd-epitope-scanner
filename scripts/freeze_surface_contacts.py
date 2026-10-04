"""Freeze explicit interface contacts before the case's physical panel design."""
from __future__ import annotations
import argparse,json,sys,datetime
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.surface_perturbation.contact_labels import freeze_contacts
from ecd_epitope_scanner.util import write_json,sha256_file

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('structure',type=Path)
    p.add_argument('--uniprot',type=Path,required=True);p.add_argument('--chain',required=True)
    p.add_argument('--partner-chain',action='append',required=True);p.add_argument('--name',required=True)
    p.add_argument('--kind',required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise ValueError('frozen labels already exist; do not silently replace')
    entry=json.loads(args.uniprot.read_text(encoding='utf-8'))
    result=freeze_contacts(args.structure,entry['primaryAccession'],args.chain,args.partner_chain,
        entry['sequence']['value'],args.name,args.kind)
    result['frozen_at_utc']=datetime.datetime.now(datetime.UTC).isoformat()
    write_json(args.output,result);print(json.dumps({'output':str(args.output),'sha256':sha256_file(args.output),
        'contact_residues':{name:len(value['uniprot_positions']) for name,value in result['contacts'].items()}}))
if __name__=='__main__':main()
