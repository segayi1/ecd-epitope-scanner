from __future__ import annotations

import hashlib
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


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        root = _project_root()
        sys.path.insert(0, str(root.parent))
        from ecd_epitope_scanner.stability import repair_with_foldx, score_many_with_foldx

        mutations = request.get("mutations") or []
        if not mutations:
            raise RuntimeError("FoldX received no mutation combinations")
        parameters = request.get("parameters", {})
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
        execution_root = _execution_root(root, output_dir, parameters)
        binary = _execution_binary(binary,execution_root,parameters)
        structure_path = Path(request["structure_path"])
        repaired = repair_with_foldx(
            str(binary), structure_path, execution_root,
            binary_sha256=binary_digest.hexdigest(),
        )
        if repaired is None:
            raise RuntimeError("FoldX RepairPDB failed")
        groups = _parse_mutations(
            mutations, str(request["chain_id"]), str(request["sequence"])
        )
        runs = max(1, int(parameters.get("number_of_runs", 5)))
        workers = max(1, int(parameters.get("max_workers", 8)))
        reused_pairs = {}
        requested_signatures={tuple(sorted(group)) for group in mutations}
        excluded = {tuple(x) for x in parameters.get("paired_cache_exclude_ids", [])}
        if parameters.get("emit_paired_models"):
            from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
            from ecd_epitope_scanner.util import sha256_file
            input_hash = sha256_file(structure_path)
            cache_roots=list(dict.fromkeys([*parameters.get('paired_cache_roots',[]),str(output_dir)]))
            for cache_root in cache_roots:
                for manifest_path in sorted(Path(cache_root).rglob("paired_model_manifest.json")):
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    if manifest.get("binary_sha256") != binary_digest.hexdigest() or manifest.get("input_structure_sha256") != input_hash:
                        continue
                    for row in manifest.get("records", []):
                        signature = tuple(sorted(row["mutations"]))
                        if signature not in requested_signatures:continue
                        for pair in row.get("paired_models", []):
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
        cached = [list(reused_pairs.get(tuple(sorted(group)), {}).values())[:runs] for group in mutations]
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
            from ecd_epitope_scanner.surface_perturbation.paired import scan_pairs,validate_pair
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
            if manifest_path.exists():
                history=output_dir/'paired_manifest_history'/str(time.time_ns());history.mkdir(parents=True,exist_ok=True)
                shutil.copy2(manifest_path,history/'paired_model_manifest.json')
            write_json(manifest_path, {
                "schema_version": "1.0", "binary_sha256": binary_digest.hexdigest(),
                "input_structure_sha256": hashlib.sha256(structure_path.read_bytes()).hexdigest(),
                "repaired_structure_sha256": hashlib.sha256(repaired.read_bytes()).hexdigest(),
                "number_of_runs": runs, "records": records,
            })
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": [],
                    "errors": [],
                    "provenance": {
                        "worker": "foldx_worker",
                        "binary": str(binary),
                        "version": "5.1",
                        "weights_sha256": binary_digest.hexdigest(),
                        "device": "cpu",
                        "command": "RepairPDB + BuildModel",
                        "number_of_runs": runs,
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
