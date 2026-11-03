"""Continue the AIB234 comparison after independent re-export acceptance."""
import json
from pathlib import Path

import run_site_first_legacy_baseline as legacy
import run_site_first_single_acceptance as current


def run(acceptance: Path, execute=legacy.execute, baseline_dir: Path | None = None):
    evidence = json.loads((acceptance / 'derived_export_report.json').read_text(encoding='utf-8'))
    if evidence.get('status') != 'PASS' or evidence.get('model_calls') != 0:
        raise ValueError('Independent re-export has not passed without new model calls')
    baseline_dir = baseline_dir or acceptance
    if (baseline_dir / 'baseline_preflight_identity.json').exists():
        raise FileExistsError('Baseline continuation has already started')
    current.assert_no_active_jobs(current.DATA / 'jobs.sqlite3')
    plan = legacy.preflight(acceptance_dir=baseline_dir)
    identities = legacy._fingerprints(plan['model_identity_paths'])
    current.atomic_json(baseline_dir / 'baseline_preflight_identity.json', identities)
    return execute(acceptance_dir=baseline_dir)


if __name__ == '__main__':
    directory = current.DATA / 'site_first_acceptance_r1'
    result = run(directory, baseline_dir=current.DATA / 'sf_old1')
    print(json.dumps({'status': result.get('status'), 'errors': result.get('errors')}, ensure_ascii=False))
