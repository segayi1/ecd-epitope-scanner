"""Revalidate the completed AIB234 run in a new export, without model execution."""
import json
import stat
import sys
import time
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from ecd_epitope_scanner.exporting import _write_bundle_archive
from ecd_epitope_scanner.joint_design_v4.impact_completeness import write_card_acceptance
from ecd_epitope_scanner.util import sha256_file, write_json
from run_ordinary_adaptive_rerun import validate_case_output
from run_site_first_single_acceptance import _independent_geometry_evidence


def extract_safe(archive: Path, destination: Path):
    with zipfile.ZipFile(archive) as bundle:
        seen = set()
        for item in bundle.infolist():
            name = PurePosixPath(item.filename)
            if (not name.parts or name.is_absolute() or '..' in name.parts or ':' in name.parts[0]
                    or '\\' in item.filename or item.filename in seen
                    or stat.S_ISLNK(item.external_attr >> 16)):
                raise ValueError(f"Unsafe or duplicate archive entry: {item.filename}")
            seen.add(item.filename)
        bundle.extractall(destination)


def run():
    started = time.monotonic()
    source = ROOT / 'workspace_data/06710f1f861b'
    acceptance = ROOT / 'workspace_data/site_first_acceptance_r1'
    destination = acceptance / 'derived_export'
    destination.mkdir(exist_ok=False)
    report = {'status': 'NEEDS_REVIEW', 'source_job': str(source), 'errors': [], 'model_calls': 0}
    try:
        original_archive = source / 'ecd_epitope_design.zip'
        original_hash = sha256_file(original_archive)
        extract_safe(original_archive, destination)
        hashes = {name: sha256_file(source / name) for name in
                  ('panel_candidates.json', 'sector_recommendations_v4.json', 'ordinary_candidate_ledger.json')}
        if any(sha256_file(destination / name) != digest for name, digest in hashes.items()):
            raise ValueError('Scientific evidence differs from source run')
        analysis_path = destination / 'analysis.json'
        before_analysis = sha256_file(analysis_path)
        analysis = json.loads(analysis_path.read_text(encoding='utf-8'))
        receipt = write_card_acceptance(destination, analysis)
        if receipt['status'] != 'PASS' or not receipt['card_count']:
            raise ValueError('Recomputed recommendation card evidence is incomplete')
        if sha256_file(analysis_path) != before_analysis:
            raise ValueError('Recomputed card receipt changed the original final analysis')
        manifest_path = destination / 'run_manifest.json'
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        for name in ('analysis.json', 'card_evidence_completeness.json'):
            manifest['files'][name] = sha256_file(destination / name)
        write_json(manifest_path, manifest)
        _write_bundle_archive(destination, [*manifest['files'], 'run_manifest.json'])
        validation = validate_case_output(destination, analysis, expected_protocol='ordinary_site_first_v2')
        geometry = _independent_geometry_evidence(destination)
        if geometry['status'] != 'PASS':
            raise ValueError('Independent geometry did not pass')
        if validation['status'] != 'PASS':
            raise ValueError('Post-reexport validation failed: ' + str(validation['errors']))
        if sha256_file(original_archive) != original_hash:
            raise ValueError('Source archive changed during re-export')
        report.update(status='PASS', validation=validation, geometry=geometry, scientific_sha256=hashes,
                      original_archive_sha256=original_hash, final_analysis_sha256=before_analysis,
                      reexport_path=str(destination))
    except Exception as exc:
        report['errors'].append(f'{type(exc).__name__}: {exc}')
    report['elapsed_reexport_seconds'] = time.monotonic() - started
    write_json(acceptance / 'derived_export_report.json', report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if report['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    run()
