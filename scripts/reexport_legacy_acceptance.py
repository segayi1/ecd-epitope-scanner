"""Independently recheck the frozen baseline's saved export without scoring models."""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from ecd_epitope_scanner.exporting import _write_bundle_archive
from ecd_epitope_scanner.joint_design_v4.impact_completeness import write_card_acceptance
from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
from ecd_epitope_scanner.util import sha256_file, write_json
from reexport_site_first_acceptance import extract_safe


def run():
    started = time.monotonic()
    root = ROOT / 'workspace_data/sf_old1'
    source = root / 'baseline/legacy_output'
    output = root / 'derived_export'
    output.mkdir(exist_ok=False)
    report = {'status': 'NEEDS_REVIEW', 'errors': [], 'model_calls': 0, 'source': str(source)}
    try:
        archive_hash = sha256_file(source / 'ecd_epitope_design.zip')
        extract_safe(source / 'ecd_epitope_design.zip', output)
        scientific = {name: sha256_file(source / name) for name in ('panel_candidates.json', 'sector_recommendations_v4.json')}
        if any(sha256_file(output / name) != value for name, value in scientific.items()):
            raise ValueError('Scientific evidence differs from source run')
        analysis_path = output / 'analysis.json'
        analysis_hash = sha256_file(analysis_path)
        analysis = json.loads(analysis_path.read_text(encoding='utf-8'))
        receipt = write_card_acceptance(output, analysis)
        if receipt['status'] != 'PASS' or not receipt['card_count']:
            raise ValueError('Recomputed legacy card evidence is incomplete')
        if sha256_file(analysis_path) != analysis_hash:
            raise ValueError('Recomputed receipt changed final legacy analysis')
        pairs = set()
        for path in [source / 'paired_model_manifest.json', *sorted((source / 'foldx_refinement').glob('batch_*/paired_model_manifest.json'))]:
            if not path.is_file():
                continue
            for row in json.loads(path.read_text(encoding='utf-8')).get('records', []):
                for pair in row.get('paired_models', []):
                    identity = (pair['batch_path'], pair['group'], pair['run'])
                    if identity not in pairs:
                        validate_pair(pair)
                        pairs.add(identity)
        if not pairs:
            raise ValueError('Legacy paired FoldX evidence is missing')
        manifest_path = output / 'run_manifest.json'
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        for name in ('analysis.json', 'card_evidence_completeness.json'):
            manifest['files'][name] = sha256_file(output / name)
        for name, digest in manifest['files'].items():
            if sha256_file(output / name) != digest:
                raise ValueError(f'Export hash mismatch: {name}')
        write_json(manifest_path, manifest)
        _write_bundle_archive(output, [*manifest['files'], 'run_manifest.json'])
        if sha256_file(source / 'ecd_epitope_design.zip') != archive_hash:
            raise ValueError('Original legacy archive changed during re-export')
        report.update(status='PASS', card_count=receipt['card_count'], validated_unique_pairs=len(pairs),
                      scientific_sha256=scientific, final_analysis_sha256=analysis_hash,
                      original_archive_sha256=archive_hash, output=str(output))
    except Exception as exc:
        report['errors'].append(f'{type(exc).__name__}: {exc}')
    report['elapsed_reexport_seconds'] = time.monotonic() - started
    write_json(root / 'derived_export_report.json', report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    run()
