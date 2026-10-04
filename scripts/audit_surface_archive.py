"""Audit ZIP membership, hashes and pairs without changing saved evidence."""
import argparse,json
from pathlib import Path
from ..surface_perturbation.archive_audit import validate_archive
from ..util import write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('archive',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    try:report=validate_archive(args.archive)
    except (ValueError,OSError,KeyError) as exc:report={'status':'FAIL','reason':str(exc)}
    write_json(args.output,report);print(json.dumps(report))
    if report['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
