"""Monitor this manifest only; audit completed jobs without submitting or retrying models."""
import datetime,json,sqlite3,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TASK=ROOT/'workspace_data/p0_p3_upgrade_20261016'
manifest=json.loads((TASK/'rerun_manifest.json').read_text(encoding='utf8'))
assert len(manifest['jobs'])==6
previous=None;audited=set()
terminal={'COMPLETE','FAILED','WAITING_FOR_AF3_UPLOAD'}
while True:
    with sqlite3.connect((ROOT/'workspace_data/jobs.sqlite3').as_uri()+'?mode=ro',uri=True) as conn:
        statuses={row['job_id']:conn.execute('select status,stage,error from jobs where id=?',(row['job_id'],)).fetchone() for row in manifest['jobs']}
    assert all(statuses.values()),'Task identity missing from store'
    current={k:list(v) for k,v in statuses.items()}
    changed = current != previous
    new_audit = False
    if changed:
        print(json.dumps({'utc':datetime.datetime.now(datetime.UTC).isoformat(),'jobs':current},ensure_ascii=False),flush=True)
        previous=current
    for row in manifest['jobs']:
        identity=row['job_id'];status=statuses[identity][0]
        receipt=TASK/(identity+'_acceptance.json')
        if status=='COMPLETE' and identity not in audited:
            if not receipt.exists():
                result=subprocess.run([sys.executable,str(ROOT/'scripts/verify_p0_p3_binder_task.py'),identity],cwd=ROOT.parent,capture_output=True,text=True)
                (TASK/(identity+'_verification.log')).write_text(result.stdout+'\n'+result.stderr,encoding='utf8')
            audited.add(identity)
            new_audit = True
            if receipt.exists():
                verdict=json.loads(receipt.read_text(encoding='utf8'))
                print(json.dumps({'job_id':identity,'engineering_acceptance':verdict.get('status'),'errors':verdict.get('errors')},ensure_ascii=False),flush=True)
    status_receipt={'updated_at_utc':datetime.datetime.now(datetime.UTC).isoformat(),'jobs':current,'audited_job_ids':sorted(audited),'all_terminal':all(v[0] in terminal for v in statuses.values())}
    pending=TASK/'monitor_status.pending';pending.write_text(json.dumps(status_receipt,ensure_ascii=False,indent=2),encoding='utf8');pending.replace(TASK/'monitor_status.json')
    if changed or new_audit or status_receipt['all_terminal']:
        result=subprocess.run([sys.executable,str(ROOT/'scripts/write_p0_p3_upgrade_report.py')],cwd=ROOT.parent,capture_output=True,text=True)
        (TASK/'report_generation.log').write_text(result.stdout+'\n'+result.stderr,encoding='utf8')
    if status_receipt['all_terminal']:
        print('ALL_SIX_TERMINAL_NO_MORE_MODEL_SUBMISSIONS',flush=True);break
    time.sleep(30)
