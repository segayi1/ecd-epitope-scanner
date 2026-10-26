"""Submit exactly the six authorized structures as fresh jobs; never resume old IDs."""
import hashlib
import json
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT/'workspace_data/p0_p3_upgrade_20261016/rerun_manifest.json'
SOURCE = ROOT/'workspace_data/binder_input_preflight_20261016/article_batch_manifest.json'

def main():
    if TARGET.exists():
        raise SystemExit('Rerun manifest already exists; inspect it instead of submitting duplicates.')
    old=json.loads(SOURCE.read_text(encoding='utf8'))
    assert {r['job_id'] for r in old['jobs']} == {'2b522bedd0f5','66006a51bddb','f0143e2eeb28','5ec2cf77db08','eeb69c7eb656','d0606f9c5aac'}
    assert old['analysis_target'] == 'binder' and old['selected_chain'] == 'A'
    staged=[]
    for row in old['jobs']:
        request=json.loads((Path(row['workdir'])/'request.json').read_text(encoding='utf8'))
        path=Path(request['structure_path'])
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['source_sha256']:
            raise ValueError('Authorized source hash changed: '+str(path))
        staged.append((row,path))
    manifest={'protocol':'ordinary_v4_sequential_v1','status':'SUBMITTING',
        'source_manifest':str(SOURCE),'source_manifest_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'deferred_af3_files':old['deferred_af3_files'],'jobs':[],
        'budget_per_protein':{'unique_combinations':80,'initial_runs':3,'maximum_refined_combinations':30,'maximum_pair_slots':300}}
    def save():
        pending=TARGET.with_suffix('.pending')
        pending.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf8');pending.replace(TARGET)
    save()
    for row,path in staged:
        with path.open('rb') as handle:
            response=requests.post('http://127.0.0.1:8765/api/jobs/upload',
                files={'structure':(row['file'],handle,'chemical/x-pdb')},
                data={'analysis_target':'binder','chain':'A','ml_mode':'strict','surface_strategy':'legacy',
                      'biological_context_mode':'provided_structure','max_sectors':'7'},timeout=60)
        response.raise_for_status()
        job=response.json()
        manifest['jobs'].append({'file':row['file'],'source_job_id':row['job_id'],
            'source_sha256':row['source_sha256'],'job_id':job['id'],
            'workdir':str(ROOT/'workspace_data'/job['id']),
            'url':'http://127.0.0.1:8765/?job='+job['id']})
        save()
        print(row['file'],job['id'],flush=True)
    manifest['status']='SUBMITTED';save()

if __name__=='__main__':main()
