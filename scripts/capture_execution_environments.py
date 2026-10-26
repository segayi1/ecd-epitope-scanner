"""Capture installed package versions without recording package-source credentials."""
import datetime,hashlib,json,subprocess,sys
from pathlib import Path
TASK=Path(__file__).resolve().parents[1]/'workspace_data/p0_p3_upgrade_20261016'
SCRIPT="import importlib.metadata as m,json,platform,sys; print(json.dumps({'python':sys.version,'prefix':sys.prefix,'platform':platform.platform(),'packages':sorted([{'name':d.metadata.get('Name','UNKNOWN'),'version':d.version} for d in m.distributions()],key=lambda d:(d['name'].lower(),d['version']))},sort_keys=True))"
rows={}
for name,command in [('workbench',[sys.executable,'-c',SCRIPT]),('spurs_wsl',['wsl.exe','-d','Ubuntu-22.04','--exec','/home/segayi1/miniconda3/envs/spurs_v4/bin/python','-c',SCRIPT])]:
    run=subprocess.run(command,capture_output=True,text=True,timeout=45,check=True)
    value=json.loads(run.stdout);rows[name]={'environment':value,'sha256':hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()}
receipt={'captured_at_utc':datetime.datetime.now(datetime.UTC).isoformat(),'environments':rows,'interpretation':'Supplemental package-version lock captured during this batch, not a claim of capture before its first model call. No package source URLs or credentials are collected.'}
path=TASK/('environment_lock_'+datetime.datetime.now(datetime.UTC).strftime('%Y%m%dT%H%M%SZ')+'.json')
path.write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'path':str(path),'package_counts':{n:len(v['environment']['packages']) for n,v in rows.items()}},ensure_ascii=False))
