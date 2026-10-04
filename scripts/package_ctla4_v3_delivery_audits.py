"""Attach stage provenance to a closed v3 delivery without changing science."""
import argparse
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

from ..util import sha256_file, write_json
from ..exporting import finalize_scientific_bundle
from ..surface_perturbation.runtime import RunBudget, read_budget
from ..surface_perturbation.case_resources import StorageGuard
from .audit_ctla4_coverage_v3 import audit


def package(case):
    route = case / 'routes/coverage_v3'
    output = route / 'results'
    if read_budget(case / 'case_budget.json')['active']:
        raise ValueError('Delivery supplement requires a closed scientific boundary')
    if not (route / 'delivery_status.json').exists():
        raise ValueError('Core delivery must complete before supplementation')
    budget = RunBudget(case / 'case_budget.json', None)
    try:
        members = {}
        for folder, prefix in ((case / 'stages', 'stages'), (route, 'route')):
            paths = folder.glob('V3_*') if prefix == 'stages' else folder.glob('*.json')
            for path in paths:
                if path.is_file() and path.suffix in {'.json', '.log'}:
                    members[f'{prefix}/{path.name}'] = path
        for name in ('case_budget.json', 'site_pair_search_contract.json',
                     'resource_contract_amendment.json'):
            members['case/' + name] = case / name
        for path in (case / 'read_only_views').glob('*'):
            if path.is_file(): members['browser/' + path.name] = path
        # Pin every actual pilot/production source, not only the last snapshot.
        for path in sorted(case.glob('source_v3*/ecd_epitope_scanner/source_overlay_manifest.json')):
            pin = json.loads(path.read_text(encoding='utf-8'))
            prefix = path.parent.parent.name
            for name, digest in pin['files'].items():
                source = path.parent / name
                if sha256_file(source) != digest:
                    raise ValueError('Historical source snapshot changed: ' + str(source))
                members[f'{prefix}/{name}'] = source
            members[f'{prefix}/{path.name}'] = path
        members['supplement_packager.py'] = Path(__file__)
        traces = {}
        for path in output.rglob('sampling_execution.json'):
            row = json.loads(path.read_text(encoding='utf-8'))
            identity = (row['start_wall_time'], row['pid'])
            if identity in traces and traces[identity]['execution'] != row:
                raise ValueError('Sampling invocation identity conflict')
            traces[identity] = {'path': str(path.relative_to(case)), 'execution': row}
        execution_report = route / 'buildmodel_invocation_audit.json'
        write_json(execution_report, {
            'audited_time_seeded_buildmodel_invocations': len(traces),
            'invocations': list(traces.values()),
            'interpretation': 'Actual sampled BuildModel launches; excludes RepairPDB, '
                              'non-seeded tools and unlogged attempts. Multiple output '
                              'pairs from one launch are not extra launches.'})
        members['route/' + execution_report.name] = execution_report
        guard = StorageGuard(case)
        guard.check(sum(p.stat().st_size for p in members.values()) * 2, solving=False)
        manifest = {name: sha256_file(path) for name, path in members.items()}
        destination = output / 'v3_stage_provenance.zip'
        pending = destination.with_suffix('.zip.pending')
        try:
            with ZipFile(pending, 'w', ZIP_DEFLATED) as archive:
                for name, path in sorted(members.items()): archive.write(path, name)
                archive.writestr('provenance_manifest.json', json.dumps({'files': manifest}, sort_keys=True))
            with ZipFile(pending) as archive:
                import hashlib
                if any(hashlib.sha256(archive.read(name)).hexdigest() != digest
                       for name, digest in manifest.items()):
                    raise ValueError('Stage provenance package hash mismatch')
            pending.replace(destination)
        finally:
            pending.unlink(missing_ok=True)
        main = json.loads((output / 'run_manifest.json').read_text(encoding='utf-8'))
        main['files'][destination.name] = sha256_file(destination)
        write_json(output / 'run_manifest.json', main)
        finalize_scientific_bundle(output)
        write_json(route / 'independent_final_package_audit.json', audit(route, packages=True))
        guard.check(solving=False, force_inventory=True)
    finally:
        budget.save(False)
    status = json.loads((route / 'delivery_status.json').read_text(encoding='utf-8'))
    status['case_budget'] = read_budget(case / 'case_budget.json')
    status['stage_provenance_sha256'] = sha256_file(destination)
    write_json(route / 'delivery_status.json', status)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    package(parser.parse_args().case.resolve())
