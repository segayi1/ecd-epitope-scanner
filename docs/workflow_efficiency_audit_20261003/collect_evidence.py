"""Read-only metadata inventory. Does not launch workers or alter run artifacts."""
from __future__ import annotations

import collections
import datetime
import heapq
import json
import os
from pathlib import Path
import shutil
import time
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
CAMPAIGN = ROOT / 'workspace_data/surface_combo_four_case_20261003/batch_20261003_155853'


def inventory(root):
    totals = collections.defaultdict(lambda: [0, 0])
    extensions = collections.defaultdict(lambda: [0, 0])
    categories = collections.defaultdict(lambda: [0, 0])
    largest = []
    count = size = 0
    errors = []
    stack = [root]
    while stack:
        folder = stack.pop()
        try:
            with os.scandir(folder) as entries:
                for entry in entries:
                    try:
                        if entry.is_symlink():
                            continue
                        if entry.is_dir(follow_symlinks=False):
                            stack.append(Path(entry.path))
                            continue
                        if not entry.is_file(follow_symlinks=False):
                            continue
                        n = entry.stat(follow_symlinks=False).st_size
                        rel = Path(entry.path).relative_to(root)
                        count += 1
                        size += n
                        for depth in range(1, min(4, len(rel.parts) - 1) + 1):
                            key = '/'.join(rel.parts[:depth])
                            totals[key][0] += n
                            totals[key][1] += 1
                        suffix = Path(entry.name).suffix.lower()
                        extensions[suffix][0] += n
                        extensions[suffix][1] += 1
                        parts = set(rel.parts)
                        category = ('zip' if suffix == '.zip' else
                                    'physical_cache' if 'physical_cache' in parts else
                                    'fields' if 'fields' in parts else
                                    'foldx_short' if 'foldx_short' in parts else
                                    'other')
                        categories[category][0] += n
                        categories[category][1] += 1
                        row = (n, str(rel))
                        if len(largest) < 60:
                            heapq.heappush(largest, row)
                        elif row > largest[0]:
                            heapq.heapreplace(largest, row)
                    except OSError as exc:
                        errors.append({'path': entry.path, 'error': str(exc)})
        except OSError as exc:
            errors.append({'path': str(folder), 'error': str(exc)})
    return {'root': str(root), 'logical_bytes': size, 'file_count': count,
            'directories': sorted([{'path': k, 'bytes': v[0], 'files': v[1]}
                                    for k, v in totals.items() if v[1] > 1],
                                   key=lambda x: x['bytes'], reverse=True)[:220],
            'extensions': dict(extensions), 'categories': dict(categories),
            'largest_files': [{'bytes': n, 'path': p} for n, p in sorted(largest, reverse=True)],
            'errors': errors}


def main():
    started = time.perf_counter()
    roots = [ROOT]
    roots += sorted(p for p in ROOT.parent.iterdir()
                    if p.is_dir() and p.name.startswith(('surface_', '.codex_v411')))
    result = {'collected_local': datetime.datetime.now().astimezone().isoformat(),
              'scope': 'logical file sizes; symlinks skipped; no allocated-block or content dedup claim',
              'roots': []}
    for root in roots:
        row = inventory(root)
        result['roots'].append(row)
        print(f"{root.name}: {row['logical_bytes']/2**30:.3f} GiB, {row['file_count']} files", flush=True)
    result['disks'] = {drive: dict(zip(('total','used','free'), shutil.disk_usage(drive)))
                       for drive in ('C:/', 'E:/')}
    result['seconds'] = time.perf_counter() - started
    (OUT/'disk_inventory.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    summaries = {}
    for case in ('CTLA4', 'PD1', 'PDL1', 'B7-H6'):
        for phase in ('calibration', 'production'):
            folder = CAMPAIGN/case/phase
            row = {}
            for name in ('surface_perturbation_preview.json', 'phase_result.json', 'case_acceptance.json',
                         f'surface_{phase}_cumulative_budget.json'):
                path = folder/name
                if path.is_file():
                    row[name] = json.loads(path.read_text(encoding='utf-8'))
            summaries[f'{case}/{phase}'] = row
    for name in ('ctla4_hard_focus_result.json','final_combo_acceptance.json','ctla4_focus_request.json',
                 'initial_native_combination_audit.json','intermediate_native_site_diagnostic.json'):
        path = CAMPAIGN/name
        if path.is_file():
            summaries[name] = json.loads(path.read_text(encoding='utf-8'))
    (OUT/'campaign_snapshot.json').write_text(json.dumps(summaries, ensure_ascii=False, indent=2), encoding='utf-8')
    archives = {}
    for phase in ('calibration','production'):
        for name in ('surface_perturbation_evidence.zip','ecd_epitope_design.zip'):
            path = CAMPAIGN/'CTLA4'/phase/name
            if not path.is_file():
                continue
            try:
                with zipfile.ZipFile(path) as z:
                    items = z.infolist()
                    groups = collections.defaultdict(lambda: [0,0,0])
                    for info in items:
                        key = info.filename.split('/')[0]
                        groups[key][0] += info.file_size
                        groups[key][1] += info.compress_size
                        groups[key][2] += 1
                    archives[f'{phase}/{name}'] = {'bytes':path.stat().st_size,'members':len(items),
                        'uncompressed_bytes':sum(x.file_size for x in items),'groups':dict(groups),
                        'note':'central-directory metadata only; no CRC integrity assertion'}
            except Exception as exc:
                archives[f'{phase}/{name}'] = {'error':str(exc)}
    (OUT/'archive_inventory.json').write_text(json.dumps(archives, ensure_ascii=False, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
