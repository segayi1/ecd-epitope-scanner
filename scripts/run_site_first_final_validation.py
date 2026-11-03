"""Run a fresh AIB234 timing after the export fix and completed frozen baseline."""
import argparse
import json
from pathlib import Path

import run_site_first_single_acceptance as current


def run(baseline_dir: Path, output_dir: Path, execute=current.execute):
    baseline = json.loads((baseline_dir / 'baseline_report.json').read_text(encoding='utf-8'))
    if baseline.get('status') != 'COMPLETE':
        raise ValueError('Frozen baseline did not complete')
    evidence = json.loads((baseline_dir / 'derived_export_report.json').read_text(encoding='utf-8'))
    if evidence.get('status') != 'PASS' or evidence.get('model_calls') != 0:
        raise ValueError('Frozen baseline independent evidence did not pass')
    current.assert_no_active_jobs(current.DATA / 'jobs.sqlite3')
    return execute(acceptance_dir=output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path, default=current.DATA / 'sf_old1',
                        help='Frozen baseline directory (default: workspace_data/sf_old1)')
    parser.add_argument('--output-dir', type=Path, default=current.DATA / 'site_first_acceptance_r2')
    args = parser.parse_args()
    baseline = args.baseline_dir.resolve()
    output = args.output_dir.resolve()
    if output.parent != current.DATA.resolve() or not output.name.startswith('site_first_acceptance_'):
        raise SystemExit('Output must be a new site_first_acceptance_ directory under workspace_data')
    if not baseline.is_dir():
        raise SystemExit(f'Baseline directory does not exist: {baseline}')
    result = run(baseline, output)
    print(json.dumps({'status': result.get('status'), 'job_id': result.get('job_id'),
                      'errors': result.get('errors', [])}, ensure_ascii=False))
