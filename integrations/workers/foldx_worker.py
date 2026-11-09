from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import shutil
import sys
import time
import traceback
from pathlib import Path
from uuid import uuid4


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _execution_root(root: Path, output: Path, parameters: dict, *, windows=None) -> Path:
    """Keep native FoldX below its measured Windows path-buffer limit."""
    windows = os.name == 'nt' if windows is None else windows
    # BuildModel appends fresh_models/engine_runs/private-batch directories.
    # The bundled binary crashes at a 239-character cwd (0xc0000409).
    if not windows or len(str(output.resolve())) + 90 <= 220:
        return output
    base = Path(parameters.get('engine_work_root') or root/'workspace_data/foldx_short')
    if not base.is_absolute():
        base = root/base
    if len(str(base.resolve())) + 110 > 220:
        raise RuntimeError('FoldX engine_work_root is too long; use a shorter persistent directory')
    job = base.resolve()/('job_'+uuid4().hex[:12])
    job.mkdir(parents=True, exist_ok=False)
    return job


def _binary(root: Path, request: dict) -> Path:
    configured = str(request.get("parameters", {}).get("binary") or "").strip()
    if configured:
        path = Path(configured)
        if not path.is_absolute():
            path = root / path
        if path.is_file():
            return path
        raise RuntimeError(f"Configured FoldX binary does not exist: {path}")
    matches = sorted((root / "vendor" / "foldx").glob("foldx_*.exe"))
    if not matches:
        raise RuntimeError("No FoldX executable was found under vendor/foldx")
    return matches[-1]


def _execution_binary(binary,execution_root,parameters,*,windows=None):
    """Use a verified short-named copy without changing the scientific binary."""
    windows=os.name=='nt' if windows is None else windows
    if not windows or not parameters.get('short_binary_alias') or binary.name=='foldx.exe':return binary
    digest=hashlib.sha256(binary.read_bytes()).hexdigest()
    base=Path(parameters.get('engine_work_root') or execution_root).resolve()
    target=base/('bin_'+digest[:12])/'foldx.exe'
    if len(str(target))>180:raise RuntimeError('FoldX binary execution alias is too long')
    target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists():shutil.copy2(binary,target)
    if hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
        raise RuntimeError('FoldX execution alias fingerprint changed')
    return target


def _parse_mutations(mutations: list[list[str]], chain_id: str, sequence: str):
    groups = []
    for group in mutations:
        parsed = []
        for notation in group:
            match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", notation.upper())
            if not match:
                raise RuntimeError(f"Invalid FoldX mutation notation: {notation}")
            wt, position_text, mutant = match.groups()
            position = int(position_text)
            if not 1 <= position <= len(sequence) or sequence[position - 1] != wt:
                raise RuntimeError(f"FoldX WT/position mismatch: {notation}")
            parsed.append((wt, position, chain_id, mutant))
        groups.append(parsed)
    return groups


def _site_first_repair(binary, source, workroot, binary_hash):
    from ecd_epitope_scanner.stability import repair_with_foldx
    from ecd_epitope_scanner.util import sha256_file, write_json
    input_hash = sha256_file(source)
    cache = Path(workroot).parent / 'foldx_repair_cache' / input_hash[:12] / binary_hash[:12]
    receipt = cache / 'repair_identity.json'
    repaired = cache / 'structure_Repair.pdb'
    reused = repaired.is_file()
    if reused:
        identity = json.loads(receipt.read_text(encoding='utf-8'))
        if identity != {'input_sha256': input_hash, 'binary_sha256': binary_hash,
                        'repaired_sha256': sha256_file(repaired)}:
            raise ValueError('Shared FoldX repair identity changed')
        return repaired, False
    repaired = repair_with_foldx(str(binary), source, Path(workroot), binary_sha256=binary_hash)
    if repaired is None:
        raise RuntimeError('FoldX RepairPDB failed')
    write_json(receipt, {'input_sha256': input_hash, 'binary_sha256': binary_hash,
                         'repaired_sha256': sha256_file(repaired)})
    return repaired, True


def _site_first_batch(request, binary, repaired, execution_root, output_dir,
                      invocation_dir, cache_identity, repair_seconds, repair_executed):
    """Collect isolated native candidate results without losing healthy peers."""
    import math
    import statistics
    from ecd_epitope_scanner.stability import score_many_with_foldx
    from ecd_epitope_scanner.integrations.paired import scan_pairs, validate_pair
    from ecd_epitope_scanner.util import sha256_file, write_json
    params = request['parameters']
    mutations = request['mutations']
    groups = _parse_mutations(mutations, str(request['chain_id']), str(request['sequence']))
    runs = int(params['number_of_runs'])
    start = int(params.get('sampling_slot_start', 0))
    if (runs, start) not in {(3, 0), (7, 3)} or params.get('sampling_protocol') != 'time_rotabase_v3':
        raise ValueError('Invalid site-first FoldX slot request')
    workers = min(8, max(1, int(params.get('max_workers', 8))))
    fresh_root = execution_root / ('fresh_' + uuid4().hex[:12])
    failures = []
    started = time.perf_counter()
    values = score_many_with_foldx(str(binary), repaired, groups, fresh_root,
        runs=runs, already_repaired=True, max_workers=workers,
        time_seed_rotabase=True, candidate_errors=failures)
    model_seconds = time.perf_counter() - started
    pairs, scan_errors = [], []
    listings = sorted(fresh_root.glob('engine_runs/*/PdbList_*.fxout'))
    for listing in listings or [None]:
        try:
            rows, errors = scan_pairs(fresh_root, listing_paths=[listing]) if listing else scan_pairs(fresh_root)
            pairs.extend(rows)
            scan_errors.extend(errors)
        except Exception as exc:
            scan_errors.append({'listing': str(listing), 'message': str(exc)})
    records, errors = [], []
    if scan_errors:
        errors.append({'code': 'PAIR_SCAN_ERRORS', 'message': str(scan_errors[:3])})
    for index, group in enumerate(mutations):
        selected = [pair for pair in pairs if tuple(sorted(pair['mutations'])) == tuple(sorted(group))]
        selected.sort(key=lambda pair: (str(pair['batch_path']), int(pair['group']), int(pair['run'])))
        row = {'mutations': group}
        try:
            if index >= len(values) or values[index] is None or not math.isfinite(float(values[index])):
                raise ValueError('FoldX candidate score missing or invalid')
            if len(selected) != runs:
                raise ValueError('FoldX candidate paired repeats incomplete')
            identities = {(str(pair['batch_path']), pair['group'], int(pair['run'])) for pair in selected}
            observed_runs = {int(pair['run']) for pair in selected}
            if (len({(pair['batch_path'], pair['group']) for pair in selected}) != 1
                    or len(identities) != runs
                    or observed_runs not in (set(range(runs)), set(range(1, runs + 1)))):
                raise ValueError('FoldX physical repeat identities are incomplete or duplicated')
            for slot, pair in enumerate(selected, start):
                pair['binary_sha256'] = cache_identity['binary_sha256']
                pair['normalized_input_sha256'] = cache_identity['structure_sha256']
                validate_pair(pair)
                sampling = Path(pair['batch_path']) / 'sampling_execution.json'
                pair.update(sampling_protocol='time_rotabase_v3', sampling_role='production',
                    sampling_cohort='ordinary_site_first_v2', sampling_slot=slot,
                    sampling_execution=json.loads(sampling.read_text(encoding='utf-8')))
            scores = [float(pair['ddg_kcal_mol']) for pair in selected]
            if not all(math.isfinite(value) for value in scores):
                raise ValueError('Non-finite paired energy')
            median = statistics.median(scores)
            # The shared FoldX API aggregates to four decimals; compare there, retain this raw median.
            if not math.isclose(round(median, 4), float(values[index]), abs_tol=1e-6, rel_tol=0):
                raise ValueError('FoldX aggregate differs from raw paired energy')
            row.update(status='success', ddg=median, paired_models=selected, ddg_runs_kcal_mol=scores)
        except Exception as exc:
            error = {'code': 'CANDIDATE_TECHNICAL_FAILURE', 'mutations': group, 'message': str(exc)}
            row.update(status='failed', errors=[error], paired_models=selected)
            errors.append(error)
        records.append(row)
    status = 'partial' if errors else 'success'
    manifest = {'schema_version': '1.0', 'execution_identity': cache_identity,
        'binary_sha256': cache_identity['binary_sha256'],
        'input_structure_sha256': cache_identity['structure_sha256'], 'number_of_runs': runs,
        'sampling_slot_contract': {'requested_slots': list(range(start, start + runs)),
                                   'raw_pair_records_are_append_only': True}, 'records': records}
    write_json(output_dir / 'paired_model_manifest.json', manifest)
    write_json(invocation_dir / 'execution_manifest.json', {
        **manifest, 'status': status, 'candidate_errors': failures, 'pair_scan_errors': scan_errors,
        'engine_receipts': {str(path): sha256_file(path) for path in fresh_root.rglob('buildmodel_request.json')}})
    return {'status': status, 'records': records, 'warnings': [], 'errors': errors,
        'provenance': {'worker': 'foldx_worker', 'execution_layout': 'candidate_isolated_v1',
            'sampling_cohort': 'ordinary_site_first_v2', 'attempted_slots': runs * len(groups),
            'repair_cache_path': str(repaired), 'repair_executed': repair_executed,
            'repair_seconds': repair_seconds, 'buildmodel_seconds': model_seconds,
            'max_workers': workers, 'invocation_directory': str(invocation_dir)},
        'runtime_seconds': repair_seconds + model_seconds}


def main() -> int:
    started = time.perf_counter()
    invocation_dir = None
    try:
        request = json.load(sys.stdin)
        root = _project_root()
        # Bind support imports to the same package that supplied this worker.
        spec = importlib.util.spec_from_file_location(
            'ecd_epitope_scanner', root / '__init__.py', submodule_search_locations=[str(root)]
        )
        package = importlib.util.module_from_spec(spec)
        sys.modules['ecd_epitope_scanner'] = package
        spec.loader.exec_module(package)
        from ecd_epitope_scanner.stability import repair_with_foldx, score_many_with_foldx

        mutations = request.get("mutations") or []
        if not mutations:
            raise RuntimeError("FoldX received no mutation combinations")
        parameters = request.get("parameters", {})
        if parameters.get('sampling_protocol') not in {None, 'time_rotabase_v3'}:
            raise ValueError('Unsupported FoldX sampling protocol')
        sampling = parameters.get('sampling_protocol') == 'time_rotabase_v3'
        role = parameters.get('sampling_role', 'production')
        cohort = parameters.get('sampling_cohort', 'production_v3')
        slot_start = int(parameters.get('sampling_slot_start', 0))
        if sampling and (role not in {'production', 'null_train', 'null_holdout'} or slot_start < 0):
            raise ValueError('Invalid v3 sampling role or slot')
        binary = _binary(root, request)
        binary_digest = hashlib.sha256()
        with binary.open("rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                binary_digest.update(chunk)
        output_dir = Path(request["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        invocation_dir = output_dir/'foldx_invocations'/uuid4().hex
        invocation_dir.mkdir(parents=True)
        (invocation_dir/'request.json').write_text(json.dumps(request, indent=2), encoding='utf-8')
        shutil.copy2(request['structure_path'], invocation_dir/'structure.pdb')
        for relative in ('integrations/workers/foldx_worker.py', 'stability.py',
                         'integrations/paired.py', 'integrations/sampling_v3.py'):
            target = invocation_dir/'source'/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root/relative, target)
        execution_root = _execution_root(root, output_dir, parameters)
        binary = _execution_binary(binary,execution_root,parameters)
        structure_path = Path(request["structure_path"])
        repair_started = time.perf_counter()
        site_first = cohort == 'ordinary_site_first_v2'
        repair_executed = False
        if site_first:
            repaired, repair_executed = _site_first_repair(binary, structure_path,
                parameters['repair_work_root'], binary_digest.hexdigest())
        else:
            repaired = repair_with_foldx(
                str(binary), structure_path, execution_root,
                binary_sha256=binary_digest.hexdigest(),
            )
        repair_seconds = time.perf_counter() - repair_started
        if repaired is None:
            raise RuntimeError("FoldX RepairPDB failed")
        groups = _parse_mutations(
            mutations, str(request["chain_id"]), str(request["sequence"])
        )
        runs = max(1, int(parameters.get("number_of_runs", 5)))
        workers = max(1, int(parameters.get("max_workers", 8)))
        cache_identity = {
            "execution_layout": "candidate_isolated_v1",
            "worker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "support_sha256": {str(rel): hashlib.sha256((root/rel).read_bytes()).hexdigest()
                for rel in ("stability.py", "integrations/paired.py", "integrations/sampling_v3.py") if (root/rel).is_file()},
            "binary_sha256": binary_digest.hexdigest(),
            "structure_sha256": hashlib.sha256(structure_path.read_bytes()).hexdigest(),
            "repaired_sha256": hashlib.sha256(repaired.read_bytes()).hexdigest(),
            "rotabase_sha256": (hashlib.sha256((repaired.parent/'rotabase.txt').read_bytes()).hexdigest()
                if (repaired.parent/'rotabase.txt').is_file() else None),
            "chain": request.get("chain_id"), "sequence": request.get("sequence"),
            "parameters": {k:v for k,v in parameters.items() if k not in
                {"paired_cache_roots", "paired_cache_exclude_ids", "max_workers", "timeout_seconds"}},
            "runtime": {"python": sys.version, "platform": sys.platform, "machine": __import__("platform").machine()},
        }
        cache_fingerprint = hashlib.sha256(json.dumps(cache_identity, sort_keys=True).encode()).hexdigest()
        if site_first:
            result = _site_first_batch(request, binary, repaired, execution_root, output_dir,
                invocation_dir, cache_identity, repair_seconds, repair_executed)
            print(json.dumps(result, sort_keys=True))
            return 0
        reused_pairs = {}
        requested_signatures={tuple(sorted(group)) for group in mutations}
        excluded = {tuple(x) for x in parameters.get("paired_cache_exclude_ids", [])}
        cache_rejections = {}
        if parameters.get("emit_paired_models"):
            from ecd_epitope_scanner.integrations.paired import validate_pair
            from ecd_epitope_scanner.util import sha256_file
            input_hash = sha256_file(structure_path)
            cache_roots=list(dict.fromkeys([str(output_dir),*parameters.get('paired_cache_roots',[])]))
            for cache_root in cache_roots:
                for manifest_path in sorted(Path(cache_root).rglob("paired_model_manifest.json")):
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    if manifest.get("binary_sha256") != binary_digest.hexdigest() or manifest.get("input_structure_sha256") != input_hash:
                        continue
                    identity = manifest.get('execution_identity', {})
                    if (identity.get('execution_layout') != 'candidate_isolated_v1'
                            or identity.get('worker_sha256') != cache_identity['worker_sha256']
                            or identity.get('repaired_sha256') != cache_identity['repaired_sha256']
                            or identity.get('rotabase_sha256') != cache_identity['rotabase_sha256']
                            or identity.get('parameters', {}).get('sampling_protocol') != parameters.get('sampling_protocol')
                            or identity.get('support_sha256') != cache_identity['support_sha256']):
                        cache_rejections['execution_protocol_mismatch'] = cache_rejections.get('execution_protocol_mismatch', 0) + 1
                        continue
                    if cohort == "ordinary_v4_sequential_v1" and manifest.get("execution_fingerprint") != cache_fingerprint:
                        continue
                    for row in manifest.get("records", []):
                        signature = tuple(sorted(row["mutations"]))
                        if signature not in requested_signatures:continue
                        for pair in row.get("paired_models", []):
                            if pair.get('execution_layout') != 'candidate_isolated_v1' or not pair.get('buildmodel_request_sha256'):
                                cache_rejections['execution_receipt_missing'] = cache_rejections.get('execution_receipt_missing', 0) + 1
                                continue
                            if sampling and (pair.get('sampling_protocol') != 'time_rotabase_v3'
                                    or pair.get('sampling_role') != role or pair.get('sampling_cohort') != cohort
                                    or not slot_start <= pair.get('sampling_slot', -1) < slot_start + runs):
                                continue
                            identity=(pair["batch_path"],pair["group"],pair["run"])
                            if identity in excluded: continue
                            if reused_pairs.get(signature,{}).get(identity)==pair:continue
                            validate_pair(pair)
                            if pair.get("ddg_kcal_mol") is not None:
                                reused_pairs.setdefault(signature, {})[identity] = pair
        cached = []
        for group in mutations:
            pairs = list(reused_pairs.get(tuple(sorted(group)), {}).values())
            if not sampling:
                # A default stream is one frozen invocation, never a mixture of cache prefixes.
                by_batch = {}
                for pair in pairs:by_batch.setdefault(pair['batch_path'], []).append(pair)
                complete = [sorted(rows, key=lambda p:int(p['run']))[:runs] for rows in by_batch.values()
                    if set(range(runs)).issubset({int(p['run']) for p in rows})]
                if pairs and not complete:
                    cache_rejections['deterministic_prefix_restart'] = cache_rejections.get('deterministic_prefix_restart', 0) + len(pairs)
                pairs = complete[0] if complete else []
            cached.append(pairs)
        if sampling:
            for pairs in cached:
                pairs.sort(key=lambda p: p['sampling_slot'])
                if len({p['sampling_slot'] for p in pairs}) != len(pairs):
                    raise ValueError('Duplicate sampling slots in cache')
        emit_pairs=bool(parameters.get('emit_paired_models'))
        missing_indices = [i for i,pairs in enumerate(cached) if len(pairs)<runs]
        values = [None]*len(mutations)
        fresh_roots=[]
        deficits={}
        for index in missing_indices:deficits.setdefault(runs-len(cached[index]),[]).append(index)
        for needed,indices in sorted(deficits.items()):
            # Immutable model directories keep prior run identities and hashes
            # valid when a later call appends evidence in the same output root.
            fresh_root=execution_root/f'fresh_models_{time.time_ns()}_{needed}' if emit_pairs else execution_root
            fresh_roots.append(fresh_root)
            fresh=score_many_with_foldx(str(binary),repaired,[groups[i] for i in indices],fresh_root,
                runs=needed,already_repaired=True,max_workers=workers,
                **({'time_seed_rotabase': True} if sampling else {}))
            if len(fresh)!=len(indices) or any(value is None for value in fresh):
                raise RuntimeError('FoldX fresh supplemental scores incomplete')
            for index,value in zip(indices,fresh):values[index]=value
        import statistics
        for index,pairs in enumerate(cached):
            if len(pairs)==runs: values[index]=statistics.median(p["ddg_kcal_mol"] for p in pairs)
        combined=None
        if emit_pairs:
            from ecd_epitope_scanner.integrations.paired import scan_pairs,validate_pair
            paired=[]
            for fresh_root in fresh_roots:
                rows,errors=scan_pairs(fresh_root)
                if errors:raise RuntimeError(f'FoldX paired structure contract failed: {errors[:3]}')
                paired.extend(rows)
            for pair in paired:
                pair['binary_sha256']=binary_digest.hexdigest()
                pair['normalized_input_sha256']=hashlib.sha256(structure_path.read_bytes()).hexdigest()
                validate_pair(pair)
            combined=[]
            for index,group in enumerate(mutations):
                signature=tuple(sorted(group))
                rows=cached[index]+[row for row in paired if tuple(sorted(row['mutations']))==signature]
                if sampling:
                    occupied={p['sampling_slot'] for p in cached[index]}
                    free=iter(s for s in range(slot_start, slot_start+runs) if s not in occupied)
                    for pair in rows[len(cached[index]):]:
                        execution=Path(pair['batch_path'])/'sampling_execution.json'
                        if not execution.is_file():raise ValueError('Sampling execution provenance missing')
                        pair.update(sampling_protocol='time_rotabase_v3', sampling_role=role,
                            sampling_cohort=cohort, sampling_slot=next(free),
                            sampling_execution=json.loads(execution.read_text(encoding='utf-8')))
                    rows.sort(key=lambda p:p['sampling_slot'])
                identities={(row['batch_path'],row['group'],row['run']) for row in rows}
                if len(rows)!=runs or len(identities)!=runs:raise RuntimeError('FoldX paired repeats incomplete or duplicated')
                if any(row.get('ddg_kcal_mol') is None for row in rows):raise RuntimeError('FoldX paired DDG incomplete')
                combined.append(rows);values[index]=statistics.median(row['ddg_kcal_mol'] for row in rows)
        if len(values) != len(mutations) or any(value is None for value in values):
            missing = sum(value is None for value in values)
            raise RuntimeError(
                f"FoldX returned incomplete scores: {missing} of {len(mutations)} missing"
            )
        records = [
            {"mutations": group, "ddg": float(values[index])}
            for index, group in enumerate(mutations)
        ]
        if parameters.get("emit_paired_models"):
            from ecd_epitope_scanner.util import write_json
            for index,record in enumerate(records):
                record["paired_models"] = combined[index]
                record["ddg_runs_kcal_mol"] = [row["ddg_kcal_mol"] for row in record["paired_models"]]
                if len(record["paired_models"]) < runs:
                    raise RuntimeError("FoldX paired structure output incomplete")
            manifest_path=output_dir/'paired_model_manifest.json'
            slot_contract = ({
                "protocol": parameters.get("sampling_protocol"),
                "cohort": cohort,
                "role": role,
                "requested_slot_start": slot_start,
                "requested_slot_end_exclusive": slot_start + runs,
                "requested_slots": list(range(slot_start, slot_start + runs)),
                "raw_pair_records_are_append_only": True,
            } if sampling else None)
            if manifest_path.exists():
                history=output_dir/'paired_manifest_history'/str(time.time_ns());history.mkdir(parents=True,exist_ok=True)
                shutil.copy2(manifest_path,history/'paired_model_manifest.json')
            write_json(manifest_path, {
                "schema_version": "1.0", "binary_sha256": binary_digest.hexdigest(),
                "execution_fingerprint": cache_fingerprint, "execution_identity": cache_identity,
                "input_structure_sha256": hashlib.sha256(structure_path.read_bytes()).hexdigest(),
                "repaired_structure_sha256": hashlib.sha256(repaired.read_bytes()).hexdigest(),
                "number_of_runs": runs, "sampling_slot_contract": slot_contract, "records": records,
            })
        (invocation_dir/'execution_manifest.json').write_text(json.dumps({
            'status': 'success', 'execution_identity': cache_identity,
            'execution_fingerprint': cache_fingerprint,
            'scores': [{'mutations': row['mutations'], 'ddg': row['ddg']} for row in records],
            'reused_model_counts': [len(rows) for rows in cached],
            'new_model_counts': [runs-len(rows) for rows in cached],
            'engine_receipts': {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                for fresh_root in fresh_roots for path in fresh_root.rglob('buildmodel_request.json')},
            'cache_rejections': cache_rejections,
        }, indent=2), encoding='utf-8')
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": [],
                    "errors": [],
                    "provenance": {
                        "worker": "foldx_worker",
                        "execution_layout": "candidate_isolated_v1",
                        "invocation_directory": str(invocation_dir),
                        "cache_rejections": cache_rejections,
                        "binary": str(binary),
                        "version": "5.1",
                        "weights_sha256": binary_digest.hexdigest(),
                        "device": "cpu",
                        "command": "RepairPDB + BuildModel",
                        "number_of_runs": runs,
                        "sampling_protocol": parameters.get("sampling_protocol") if sampling else None,
                        "sampling_role": role if sampling else None,
                        "sampling_cohort": cohort if sampling else None,
                        "sampling_slot_start": slot_start if sampling else None,
                        "sampling_slots": list(range(slot_start, slot_start + runs)) if sampling else None,
                        "max_workers": workers,
                        "reused_model_counts": [len(rows) for rows in cached],
                        "new_model_counts": [runs-len(rows) for rows in cached],
                        "repair_cache_path": str(repaired),
                        "execution_root": str(execution_root),
                        "output_dir": str(output_dir),
                    },
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0
    except Exception as exc:
        if invocation_dir is not None:
            (invocation_dir/'failure.json').write_text(json.dumps({'status':'failed', 'error':str(exc)}), encoding='utf-8')
        traceback.print_exc(file=sys.stderr)
        print(
            json.dumps(
                {
                    "status": "failed",
                    "records": [],
                    "warnings": [],
                    "errors": [{"code": "TOOL_CRASHED", "message": str(exc)}],
                    "provenance": {"worker": "foldx_worker"},
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
