from pathlib import Path
from ecd_epitope_scanner import jobs
from ecd_epitope_scanner.pipeline import PipelineRequest


def test_resumed_job_preserves_raw_models_and_records_queue(tmp_path,monkeypatch):
    store=jobs.JobStore(tmp_path)
    job=store.create(PipelineRequest(input_type='structure',identifier='TEST'))
    prior=Path(job['workdir'])/'engine_runs'
    prior.mkdir(); (prior/'evidence.txt').write_text('immutable')
    captured={}
    def run(request,progress):
        captured.update(request.metadata)
        return {'design_summary':{}}
    monkeypatch.setattr(jobs,'run_pipeline',run)
    runner=jobs.JobRunner(store)
    try:
        runner._run(job['id'])
    finally:
        runner.executor.shutdown(wait=True)
    assert store.get(job['id'])['status']=='COMPLETE'
    assert [p.read_text() for p in (Path(job['workdir'])/'prior_attempts').rglob('evidence.txt')]==['immutable']
    assert (prior/'evidence.txt').read_text()=='immutable'
    assert captured['execution_queue']['job_id']==job['id']

