"""Reconcile immutable stage receipts against a closed cumulative ledger."""
import argparse
import json
import math
from pathlib import Path

from ..util import write_json, sha256_file


def reconcile(case):
    ledger_path = case/'case_budget.json'
    ledger = json.loads(ledger_path.read_text(encoding='utf-8'))
    if ledger.get('active'):
        raise ValueError('requires a closed case ledger')
    total = float(ledger['elapsed_seconds'])
    rows = []
    for path in (case/'stages').glob('*_report.json'):
        r = json.loads(path.read_text(encoding='utf-8'))
        end, charge = float(r['case_elapsed_seconds']), float(r['charged_seconds'])
        if not all(math.isfinite(v) and v >= 0 for v in (end, charge)):
            raise ValueError('invalid receipt time: '+str(path))
        rows.append(dict(stage=r['stage'], end=end, charge=charge, sha256=sha256_file(path)))
    rows.sort(key=lambda r: (r['end'], r['stage']))
    if len({r['stage'] for r in rows}) != len(rows):
        raise ValueError('duplicate scientific stage receipts')
    intervals = []; previous = 0.
    for r in rows:
        difference = r['end']-previous-r['charge']
        intervals.append(dict(before_stage=r['stage'], previous_case_seconds=previous,
                              receipt_end_seconds=r['end'], receipt_charge_seconds=r['charge'],
                              difference_seconds=difference,
                              allocation='UNATTRIBUTED_BETWEEN_RECEIPTS' if difference >= 0 else 'RECEIPT_TIMING_OVERLAP'))
        previous = r['end']
    tail = total-previous
    if tail < -1e-6:
        raise ValueError('receipt exceeds closed ledger')
    receipts = math.fsum(r['charge'] for r in rows)
    residual = math.fsum(i['difference_seconds'] for i in intervals)+tail
    error = total-receipts-residual
    if abs(error) > 1e-6:
        raise ValueError('reconciliation arithmetic failed')
    return dict(schema='case_cost_reconciliation_v1', total_seconds=total,
                ledger_sha256=sha256_file(ledger_path), receipt_seconds=receipts,
                separately_reported_residual_seconds=residual, finalization_tail_seconds=tail,
                accounting_error_seconds=error, total_status='RECONCILED',
                allocation_status='PARTIAL' if abs(residual) > 1e-6 else 'COMPLETE',
                intervals=intervals, receipts=rows, new_scientific_calls=0,
                interpretation='Residual is not assigned to solvers or refunded; historical receipts unchanged')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError('preserve existing reconciliation report')
    result = reconcile(args.case)
    write_json(args.report, result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('intervals', 'receipts')}))
