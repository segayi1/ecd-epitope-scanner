"""Separate unqueried WT backgrounds from genuine null/numerical uncertainty."""
import argparse,json
from pathlib import Path
import numpy as np
from ..util import sha256_file,write_json


def audit(case):
    folder=Path(case)/'calibration/surface_perturbation'
    candidates=[]
    for path in (folder/'calibration').glob('*.json'):
        meta=json.loads(path.read_text(encoding='utf-8'))
        if meta.get('runs')==20 and any('required_area' in v for v in meta.get('channels',{}).values()):
            candidates.append((path,meta))
    candidates.sort(key=lambda item:item[0].stat().st_mtime,reverse=True)
    catalog={}
    for path in (folder/'fields').glob('field_*.json'):
        meta=json.loads(path.read_text(encoding='utf-8'));contract=meta['contract']
        catalog[(contract['reference'],contract['wt'],contract['mut'])]=(path,meta)
    with np.load(folder/'canonical_reference.npz',allow_pickle=False) as data:area=data['area'].copy()
    reports=[]
    for path,meta in candidates:
        source=path.with_suffix('.npz')
        if sha256_file(source)!=meta['array_sha256']:raise ValueError('changed null arrays')
        with np.load(source,allow_pickle=False) as data:
            null={k:data[k].copy() for k in data.files}
        required=null['shape_sdf_a__required_mask'];denominator=float(area[required].sum())
        def fraction(mask):return float(area[required&mask].sum())/denominator if denominator else None
        masks=[];manifests=[];topology=[]
        for pair in meta['pairs']:
            field_path,field_meta=catalog[(meta['reference_id'],pair['wt_sha256'],pair['mut_sha256'])]
            values=field_path.with_suffix('.npz')
            if sha256_file(values)!=field_meta['sha256']:raise ValueError('changed raw WT arrays')
            with np.load(values,allow_pickle=False) as data:
                # The polar-support component is always numerically defined
                # at queried rows. Donor/direction components can deliberately
                # abstain near HIS or unavailable atom directions; those NaNs
                # must NOT be mislabeled as missing spatial queries.
                masks.append(np.isfinite(data['wt_chemistry'][:,1]))
                topology.append(bool(data['topology_changed']))
            manifests.append(str(field_path))
        queried=np.asarray(masks);all_queried=queried.all(axis=0)
        channels={}
        for channel in meta['channels']:
            blocks=null[channel+'__null_blocks'];error=null[channel+'__numerical_error']
            finite=np.isfinite(blocks).all(axis=0);numeric=np.isfinite(error)&(error>=0)
            valid=null[channel+'__valid']
            channels[channel]={'valid_local_area_fraction':fraction(valid),
                'missing_WT_query_area_fraction':fraction(~all_queried&~finite),
                'queried_but_null_invalid_area_fraction':fraction(all_queried&~finite),
                'numerical_error_unknown_area_fraction':fraction(~numeric),
                'finite_but_unconverged_area_fraction':fraction(finite&numeric&~valid)}
        reports.append({'protocol_id':meta['protocol_id'],'null_manifest':str(path),'runs':20,
            'required_area_a2':denominator,'all_WT_queries_present_local_area_fraction':fraction(all_queried),
            'first_WT_query_present_local_area_fraction':fraction(queried[0]),
            'topology_changed_repeats':topology,
            'raw_manifests':manifests,'channels':channels})
    return {'case':Path(case).name,'scope':'read-only saved input/query provenance diagnosis; no new scientific acceptance',
        'status':'DIAGNOSED' if reports else 'NO_20_RUN_LOCAL_CALIBRATION','cohorts':reports}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result={case:audit(args.root/case) for case in ('PD1','PDL1','B7-H6')}
    write_json(args.output,result)
    print(json.dumps({case:[{'WT_scope_fraction':r['all_WT_queries_present_local_area_fraction'],
        'shape':r['channels']['shape_sdf_a'],'donor':r['channels']['donor_support']} for r in value['cohorts']]
        for case,value in result.items()}))


if __name__=='__main__':main()
