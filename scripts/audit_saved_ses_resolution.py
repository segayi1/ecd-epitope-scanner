"""Read-only topology sensitivity of already computed CTLA4 SES meshes."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.surface_perturbation.tools import read_off
from ecd_epitope_scanner.util import sha256_file,write_json


def audit(case):
    groups={}
    for manifest in (case/'physical_cache').glob('ses_*/complete.json'):
        row=json.loads(manifest.read_text(encoding='utf-8'))
        key=(row['xyzr_sha256'],row['binary_sha256'],row['probe_radius_a'],row['cavity_policy'])
        groups.setdefault(key,{})[float(row['scale'])]=(manifest,row)
    comparisons=[]
    for key,scales in groups.items():
        if not {2.,3.}<=scales.keys():continue
        values={}
        for scale in (2.,3.):
            manifest,row=scales[scale];path=manifest.parent/row['mesh']
            if sha256_file(path)!=row['mesh_sha256']:raise ValueError('Saved SES content changed')
            mesh=read_off(path)
            values[str(scale)]={'euler_characteristic':mesh.euler_characteristic,'closed':bool(mesh.closed),
                'area_a2':float(mesh.area.sum()),'mesh_sha256':row['mesh_sha256']}
        comparisons.append({'atomic_configuration_sha256':key[0],'scales':values,
            'topology_resolution_sensitive':values['2.0']['euler_characteristic']!=values['3.0']['euler_characteristic']})
    return {'engineering_status':'PASS','new_solver_calls':0,'comparisons':comparisons,
        'resolution_sensitive_count':sum(r['topology_resolution_sensitive'] for r in comparisons),
        'interpretation':'Observed mesh-resolution sensitivity only. Official topology masks and qualification are unchanged.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    report=audit(args.case);write_json(args.output,report)
    print(json.dumps({'comparisons':len(report['comparisons']), 'resolution_sensitive':report['resolution_sensitive_count'],'new_solver_calls':0}))
if __name__=='__main__':main()
