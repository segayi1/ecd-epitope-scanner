"""Run the frozen 13-case manifest with the explicitly enabled site-first protocol."""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
import run_ordinary_adaptive_rerun as batch

OUT = ROOT / 'workspace_data/site_first_13_cases_20261029'
STATE = OUT / 'state.json'
RECEIPT = OUT / 'launcher_receipt.json'


def request_for(case, structure, confidence):
    request = original_request(case, structure, confidence)
    request.metadata['joint_v4_ml'] = {
        'execution_protocol': {'name': 'ordinary_site_first_v2'},
        'runtime': {'cache_dir': str(structure.parent.parent / 'site_first_score_cache')},
    }
    return request


def validate(workdir, analysis=None, expected_protocol='ordinary_site_first_v2'):
    return original_validate(workdir, analysis, expected_protocol=expected_protocol)


original_request = batch.request_for
original_validate = batch.validate_case_output


def main():
    manifest = json.loads(batch.MANIFEST.read_text(encoding='utf-8'))
    cases = batch.preflight(manifest)
    assert len(cases) == 13
    for case in cases:
        request = request_for(case, Path(case['structure_path']), Path(case['confidence_path']) if case.get('confidence_path') else None)
        assert request.metadata['joint_v4_ml']['execution_protocol']['name'] == 'ordinary_site_first_v2'
        if case['input_type'] == 'af3_upload':
            assert request.confidence_path
    batch.assert_no_active_jobs(ROOT / 'workspace_data/jobs.sqlite3')
    OUT.mkdir(parents=True, exist_ok=False)
    acceptance = json.loads((ROOT / 'workspace_data/site_first_acceptance_aib234_20261028_r1/report.json').read_text(encoding='utf-8'))
    if acceptance.get('status') != 'PASS' or acceptance.get('validation', {}).get('checks', {}).get('formal_display_eligibility') is not True:
        raise ValueError('AIB234 site-first acceptance evidence is not PASS')
    receipt = {'status': 'PASS', 'basis': 'frozen 13-input preflight and AIB234 site-first acceptance',
               'manifest_sha256': hashlib.sha256(batch.MANIFEST.read_bytes()).hexdigest(),
               'runner_sha256': batch.sha256(Path(batch.__file__).resolve()),
               'production_source_sha256': batch.production_source_sha256(),
               'launcher_sha256': batch.sha256(Path(__file__).resolve())}
    batch.atomic_json(RECEIPT, receipt)
    batch.request_for = request_for
    batch.validate_case_output = validate
    batch.OUT = OUT / 'report'
    batch.execute(manifest, RECEIPT, state_path=STATE)


if __name__ == '__main__':
    main()
