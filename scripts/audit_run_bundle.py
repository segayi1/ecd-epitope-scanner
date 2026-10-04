"""Verify the delivered outer manifest, ZIP and on-disk file identities."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from ..util import sha256_file,write_json
from ..surface_perturbation.case_status import case_status


def audit_run_bundle(folder):
    folder=Path(folder);manifest_path=folder/'run_manifest.json'
    manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    expected=set(manifest['files'])|{'run_manifest.json'}
    if any(Path(name).name!=name or '\\' in name for name in expected):raise ValueError('unsafe artifact name')
    with zipfile.ZipFile(folder/'ecd_epitope_design.zip') as archive:
        names=archive.namelist()
        if len(names)!=len(set(names)) or set(names)!=expected:raise ValueError('bundle membership differs')
        for name in expected:
            digest=sha256_file(folder/name)
            if hashlib.sha256(archive.read(name)).hexdigest()!=digest:raise ValueError('ZIP differs from disk: '+name)
            if name!='run_manifest.json' and manifest['files'][name]!=digest:raise ValueError('manifest differs: '+name)
    return {'status':'PASS','members':len(expected),'case_status':case_status(folder),
            'scope':'artifact consistency only; runtime and scientific acceptance remain separate'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('folder',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=audit_run_bundle(args.folder);write_json(args.output,result);print(json.dumps(result))


if __name__=='__main__':main()
