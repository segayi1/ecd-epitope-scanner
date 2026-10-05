"""Attach stage provenance to a closed v3 delivery without changing science."""
import argparse
import json
import shutil
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
    closed_budget = route / 'delivery_budget_before_supplement.json'
    write_json(closed_budget, read_budget(case / 'case_budget.json'))
    budget = RunBudget(case / 'case_budget.json', None)
    try:
        members = {}
        for folder, prefix in ((case / 'stages', 'stages'), (route, 'route')):
            paths = folder.glob('V3_*') if prefix == 'stages' else folder.glob('*')
            for path in paths:
                if path.is_file() and path.suffix in {'.json', '.log', '.md'}:
                    members[f'{prefix}/{path.name}'] = path
        for name in ('site_pair_search_contract.json',
                     'resource_contract_amendment.json'):
            members['case/' + name] = case / name
        members['case/closed_budget_before_delivery_provenance.json'] = closed_budget
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
        members['independent_audit.py'] = Path(__file__).with_name('audit_ctla4_coverage_v3.py')
        members['report_writer.py'] = Path(__file__).with_name('write_ctla4_v3_run_report.py')
        members['core_delivery_packager.py'] = Path(__file__).with_name('run_ctla4_coverage_v3.py')
        members['delivery_tests/test_ctla4_delivery_genotype.py'] = Path(__file__).resolve().parents[1]/'tests/test_ctla4_delivery_genotype.py'
        documentation = Path(__file__).resolve().parents[1] / 'docs'
        for name in ('ctla4_coverage_v3_implementation_audit_20261004.md',
                     'ctla4_stability_and_physical_evidence_research_20261004.md'):
            members['documentation/' + name] = documentation / name
        traces = {}
        for path in (p for root in (output, case / 'fx') for p in root.rglob('sampling_execution.json')):
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
        invocation_copy=output/execution_report.name
        shutil.copy2(execution_report,invocation_copy)
        used=guard.check(solving=False,force_inventory=True)
        archive_reserve=(output/'ecd_epitope_design.zip').stat().st_size+invocation_copy.stat().st_size+65536
        guard.check(archive_reserve,solving=False)
        resource=output/'final_resource_audit.json'
        write_json(resource,{'engineering_status':'PASS',
            'measurement_phase':'provenance complete, before final resource receipt attachment',
            'case_owned_bytes_at_measurement':used,'free_disk_bytes_at_measurement':shutil.disk_usage(case).free,
            'minimum_free_disk_bytes':guard.minimum_free,
            'observed_stage_peak_owned_bytes':max(json.loads(p.read_text(encoding='utf-8'))['peak_owned_bytes']
                for p in (case/'stages').glob('V3_*_report.json')),
            'conservative_final_archive_write_peak_bytes':used+archive_reserve,
            'case_budget_before_resource_attachment':read_budget(case/'case_budget.json'),
            'component_archives':{p.name:{'bytes':p.stat().st_size,'sha256':sha256_file(p)}
                for p in (output/'surface_perturbation_evidence.zip',output/'implementation_source_v3.zip',destination)},
            'interpretation':'observed and conservative peak distinguished; closed final cost and exact final ZIP size in external delivery_status.json'})
        if resource.stat().st_size>65536:raise ValueError('Resource receipt exceeded archive reservation')
        for path in (resource,invocation_copy):main['files'][path.name]=sha256_file(path)
        for name in ('surface_perturbation_local_constructs.csv', 'CTLA4_v3_final_report.md'):
            if (output/name).is_file(): main['files'][name]=sha256_file(output/name)
        write_json(output / 'run_manifest.json', main)
        finalize_scientific_bundle(output)
        write_json(route / 'independent_final_package_audit.json', audit(route, packages=True))
        guard.check(solving=False, force_inventory=True)
    finally:
        budget.save(False)
    status = json.loads((route / 'delivery_status.json').read_text(encoding='utf-8'))
    status['case_budget'] = read_budget(case / 'case_budget.json')
    status['stage_provenance_sha256'] = sha256_file(destination)
    status['final_archive_bytes']=(output/'ecd_epitope_design.zip').stat().st_size
    status['final_archive_sha256']=sha256_file(output/'ecd_epitope_design.zip')
    status['final_resource_observation']={'case_owned_bytes':guard.check(solving=False,force_inventory=True),
        'free_disk_bytes':shutil.disk_usage(case).free,'minimum_free_disk_bytes':guard.minimum_free}
    write_json(route / 'delivery_status.json', status)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', type=Path)
    package(parser.parse_args().case.resolve())
