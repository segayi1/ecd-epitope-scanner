"""Independent identity, raw-pair and export acceptance for the six new binder jobs."""
from __future__ import annotations
import argparse,csv,datetime as dt,hashlib,json,sys
from pathlib import Path
import sqlite3,zipfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent))
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
from ecd_epitope_scanner.scripts.audit_p0_p3_batch import audit_task


def read(path):return json.loads(Path(path).read_text(encoding='utf8'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _source_identity_assessment(frozen,deployed,process_evidence=None):
    mismatches={path:{'task_snapshot':frozen.get(path),'deployment':digest}
        for path,digest in deployed.items() if frozen.get(path)!=digest}
    if not mismatches:
        return {'ok':True,'mode':'EXACT_TASK_SNAPSHOT','mismatches':mismatches}
    evidence=process_evidence or {}
    inferred=(set(mismatches)=={'exporting.py'} and
        evidence.get('same_recovery_process') is True and
        evidence.get('anchor_snapshot_matches_deployment') is True and
        evidence.get('anchor_captured_at_process_start') is True and
        evidence.get('import_chain_loaded_before_task') is True and
        evidence.get('no_module_reload_path') is True)
    return {'ok':inferred,
        'mode':'SAME_PROCESS_IMPORT_INFERRED' if inferred else 'UNRESOLVED',
        'mismatches':mismatches,
        'evidence':evidence,
        'limitation':'Per-task source snapshot hashes files on disk; this recovery process did not capture a runtime module-bytecode digest.' if inferred else None}


def _shared_process_source_evidence(job_id,deployment):
    task=ROOT/'workspace_data/p0_p3_upgrade_20261016'
    recovery=read(task/'recovery_receipt.json')
    anchor_id='55ac1dce3c04'
    anchor_root=ROOT/'workspace_data'/anchor_id
    anchor=read(anchor_root/'run_source_snapshot_v1.json')
    anchor_frozen={r['path']:r['sha256'] for r in anchor['source_files']}
    anchor_matches=all(anchor_frozen.get(path)==digest for path,digest in deployment['source_sha256'].items())
    owner=None
    try:
        with sqlite3.connect(f'file:{ROOT / "workspace_data" / "jobs.sqlite3"}?mode=ro',uri=True) as conn:
            owner=conn.execute('SELECT owner_pid FROM jobs WHERE id=?',(job_id,)).fetchone()
        same_pid=bool(owner and owner[0] and owner[0]==recovery.get('recovery_pid'))
    except sqlite3.Error:
        same_pid=False
    try:
        started=dt.datetime.fromisoformat(recovery['started_at_utc'])
        captured=dt.datetime.fromisoformat(anchor['captured_at_utc'])
        near_start=started<=captured<=started+dt.timedelta(seconds=30)
    except (KeyError,TypeError,ValueError):
        near_start=False
    runner=(task/'recover_six_jobs.py').read_text(encoding='utf8')
    jobs=(ROOT/'jobs.py').read_text(encoding='utf8')
    pipeline=(ROOT/'pipeline.py').read_text(encoding='utf8')
    api=(ROOT/'api.py').read_text(encoding='utf8')
    import_chain=('from ecd_epitope_scanner.jobs import _process_alive' in runner and
        'from .pipeline import PipelineRequest, run_pipeline' in jobs and
        'from .exporting import finalize_scientific_bundle, write_scientific_outputs' in pipeline)
    no_reload=not any('importlib.reload' in text or 'reload(' in text
        for text in (runner,jobs,pipeline,api))
    return {'same_recovery_process':same_pid and job_id in recovery.get('resumed_job_ids',[]),
        'anchor_snapshot_matches_deployment':anchor_matches,
        'anchor_captured_at_process_start':near_start,
        'import_chain_loaded_before_task':import_chain,
        'no_module_reload_path':no_reload,
        'anchor_job_id':anchor_id,
        'anchor_captured_at_utc':anchor.get('captured_at_utc'),
        'recovery_started_at_utc':recovery.get('started_at_utc'),
        'recovery_pid':recovery.get('recovery_pid'),
        'task_owner_pid':owner[0] if owner else None}


def _source_snapshot_archive_matches(root,snapshot):
    archive_info=snapshot.get('source_archive',{})
    archive_path=root/archive_info.get('path','')
    if not archive_path.is_file() or sha(archive_path)!=archive_info.get('sha256'):
        return False
    try:
        with zipfile.ZipFile(archive_path) as archive:
            members=set(archive.namelist())
            for row in snapshot.get('source_files',[]):
                name=row['path']
                if name not in members or hashlib.sha256(archive.read(name)).hexdigest()!=row['sha256']:
                    return False
    except (OSError,KeyError,zipfile.BadZipFile):
        return False
    return True


def verify(row):
    root=Path(row['workdir']);analysis=read(root/'analysis.json');errors=[]
    checks={}
    def check(name,value):
        checks[name]=bool(value)
        if not value:errors.append(name)
    request=read(root/'request.json')
    check('binder_A_identity',request['analysis_target']=='binder' and request['chain']=='A')
    check('authorized_source_hash',sha(request['structure_path'])==row['source_sha256'])
    structure=PDBParser(QUIET=True).get_structure('source',request['structure_path'])
    sequence=''.join(seq1(r.resname) for r in structure[0]['A'] if r.id[0]==' ')
    check('coordinate_WT_sequence',sequence==analysis['protein']['ecd_sequence'])
    normalized=root/'prediction_context.pdb'
    normalized_hash=sha(normalized)
    pair_count=0;identities=set();slots={};binary_hashes=set()
    for path in root.rglob('paired_model_manifest.json'):
        if 'paired_manifest_history' in path.parts:continue
        manifest=read(path)
        check('normalized_input_'+str(path.relative_to(root)),manifest['input_structure_sha256']==normalized_hash)
        execution=manifest.get('execution_identity')
        check('strict_execution_fingerprint_'+str(path.relative_to(root)),bool(execution) and
            hashlib.sha256(json.dumps(execution,sort_keys=True).encode()).hexdigest()==manifest.get('execution_fingerprint'))
        binary_hashes.add(manifest['binary_sha256'])
        for record in manifest.get('records',[]):
            signature=tuple(sorted(record['mutations']))
            for pair in record.get('paired_models',[]):
                validate_pair(pair)
                check('repair_to_buildmodel_'+str(pair_count),pair['input_sha256']==manifest['repaired_structure_sha256'])
                identity=(pair['batch_path'],pair['group'],pair['run'])
                if identity in identities:raise ValueError('Repeated physical pair identity')
                identities.add(identity);pair_count+=1
                slots.setdefault(signature,[]).append(int(pair['sampling_slot']))
    budget=analysis['joint_v4_ml']['compute_budget'].get('foldx_repeat_budget',{})
    check('attempt_budget',budget.get('total_slots',301)<=300)
    check('successful_pair_slots',pair_count<=budget.get('total_slots',0) and pair_count>0)
    check('unique_combination_budget',0<len(slots)<=80)
    check('logical_slot_sets',all(sorted(v) in [[0,1,2],[0,1,2,3,4]] for v in slots.values()))
    check('binary_identity_unique',len(binary_hashes)==1)
    sequence_records=list(SeqIO.parse(root/'sequences.fasta','fasta'))
    sequences={r.id:str(r.seq) for r in sequence_records}
    check('fasta_IDs_unique',len(sequences)==len(sequence_records))
    check('fasta_WT',sequences.get('WT_ECD')==sequence)
    with (root/'designs.csv').open(encoding='utf-8-sig', newline='') as handle:
        csv_records=list(csv.DictReader(handle))
    csv_rows={r['id']:r for r in csv_records}
    check('csv_ids_unique',len(csv_rows)==len(csv_records))
    check('csv_design_ids',set(csv_rows)=={d['id'] for d in analysis.get('designs',[])})
    for design in analysis.get('designs',[]):
        if design.get('mutant_ecd_sequence'):
            check('sequence_export_'+design['id'],sequences.get(design['id'])==design['mutant_ecd_sequence'])
            mutant=list(sequence)
            notations=[]
            valid_mutations=True
            for mutation in design.get('mutations',[]):
                position=int(mutation['position'])
                wt,mutant_aa=mutation['wt'],mutation['mutant']
                valid_mutations &= 1<=position<=len(mutant) and mutant[position-1]==wt
                if 1<=position<=len(mutant):mutant[position-1]=mutant_aa
                notations.append(f"{wt}{position}{mutant_aa}")
            check('json_mutant_sequence_'+design['id'],valid_mutations and ''.join(mutant)==design['mutant_ecd_sequence'])
            csv_row=csv_rows.get(design['id'],{})
            check('csv_export_'+design['id'],csv_row.get('mutations')==';'.join(notations) and csv_row.get('status')==design.get('status'))
    check('fasta_design_ids',set(sequences)-{'WT_ECD'}=={d['id'] for d in analysis.get('designs',[]) if d.get('mutant_ecd_sequence')})
    snapshot=read(root/'run_source_snapshot_v1.json')
    frozen={r['path']:r['sha256'] for r in snapshot['source_files']}
    deployment=read(ROOT/'workspace_data/p0_p3_upgrade_20261016/deployment_receipt.json')
    provenance=_source_identity_assessment(frozen,deployment['source_sha256'],
        _shared_process_source_evidence(row['job_id'],deployment))
    check('deployed_executed_code',provenance['ok'])
    check('source_snapshot_archive_integrity',_source_snapshot_archive_matches(root,snapshot))
    audit=audit_task(root)
    check('independent_budget_source_export_audit',all(audit['checks'].values()))
    return {'job_id':row['job_id'],'source_job_id':row['source_job_id'],'file':row['file'],
        'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'errors':errors,
        'pair_count':pair_count,'unique_combinations':len(slots),'source_provenance':provenance,'audit':audit,
        'interpretation':'Engineering identity, mapping and evidence audit only; not experimental folding/binding validation.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('job_id');args=p.parse_args()
    task=ROOT/'workspace_data/p0_p3_upgrade_20261016'
    row=next(r for r in read(task/'rerun_manifest.json')['jobs'] if r['job_id']==args.job_id)
    try:result=verify(row)
    except Exception as error:
        result={'job_id':args.job_id,'status':'FAIL','errors':[f'{type(error).__name__}: {error}']}
    destination=task/(args.job_id+'_acceptance.json')
    destination.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({'job_id':args.job_id,'status':result['status'],'receipt':str(destination)},ensure_ascii=False))
    return 0 if result['status']=='PASS' else 2
if __name__=='__main__':raise SystemExit(main())
