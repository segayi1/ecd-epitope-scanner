"""Portable, hash-checked FoldX source evidence for exported runs."""
from __future__ import annotations

import json
import shutil
from pathlib import Path, PurePosixPath, PureWindowsPath

from .util import sha256_file, write_json


def copy_foldx_source_evidence(source: Path, output: Path, source_manifest: dict) -> dict:
    """Copy and validate source FoldX batches from legacy and invocation layouts."""
    from ecd_epitope_scanner.surface_perturbation.paired import scan_pairs, validate_pair

    foldx = source_manifest.get("worker_runs", {}).get("foldx", {})
    provenance = foldx.get("provenance", {})
    if foldx.get("status") != "success":
        return {"status": "UNKNOWN", "reason": "source_foldx_worker_not_successful"}
    evidence_root = output / "source_evidence"
    existing_manifest = evidence_root / "foldx_source_evidence.json"
    if existing_manifest.is_file():
        try:
            payload = json.loads(existing_manifest.read_text(encoding="utf-8"))
            existing = {"manifest_path": existing_manifest.relative_to(output).as_posix(),
                        "manifest_sha256": sha256_file(existing_manifest),
                        "portable_pair_manifest_path": payload.get("portable_pair_manifest_path"),
                        "artifact_paths": [row["path"] for row in payload["raw_artifacts"]] +
                                          [existing_manifest.relative_to(output).as_posix()]}
        except (OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError("Existing portable FoldX evidence manifest is unreadable") from exc
        if not verify_portable_source_evidence(output, existing):
            raise ValueError("Existing portable FoldX source evidence failed hash verification")
        return existing
    copied, batches, destinations = [], set(), {}
    source_root = source.resolve()
    roots = [source_root]
    other_roots = set()
    for key in ("execution_root", "invocation_directory"):
        value = provenance.get(key)
        if value and Path(value).is_dir():
            resolved = Path(value).resolve()
            if resolved != source_root:
                other_roots.add(resolved)
    roots.extend(sorted(other_roots, key=str))

    def register_batch(batch, root):
        batch = Path(batch).resolve()
        if any(part in {"prior_attempts", "paired_manifest_history"} for part in batch.parts):
            return
        if batch in batches:
            return
        batches.add(batch)
        relative = batch.relative_to(root) if batch.is_relative_to(root) else None
        if relative is None:
            relative = Path("external_batches") / f"{batch.name}_{sha256_file(next(batch.glob('PdbList_*.fxout')))[:12]}"
        elif root != source_root:
            relative = Path("external_execution") / root.name / relative
        destinations[batch] = relative

    for root in roots:
        for listing in root.rglob("PdbList_*.fxout"):
            if "source_evidence" in listing.parts:
                continue
            register_batch(listing.parent, root)
        for invocation in sorted((root / "foldx_invocations").glob("*")):
            if invocation.is_dir():
                target = evidence_root / "foldx_invocations" / invocation.name
                for item in sorted(invocation.rglob("*")):
                    if item.is_file():
                        target_file = target / item.relative_to(invocation)
                        target_file.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(item, target_file)
                        digest = sha256_file(item)
                        if digest != sha256_file(target_file):
                            raise ValueError(f"FoldX invocation snapshot changed while copying: {item}")
                        copied.append({"path": target_file.relative_to(output).as_posix(),
                                       "sha256": digest, "source_sha256": digest,
                                       "size_bytes": target_file.stat().st_size})

    pair_manifest_path = source / "paired_model_manifest.json"
    pair_manifest = None
    manifest_pairs = {}
    if pair_manifest_path.is_file():
        try:
            pair_manifest = json.loads(pair_manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError("FoldX paired-model manifest is unreadable") from exc
        for record in pair_manifest.get("records", []):
            for pair in record.get("paired_models", []):
                validate_pair(pair)
                if pair.get("batch_path"):
                    batch = Path(pair["batch_path"]).resolve()
                    if not batch.is_dir():
                        raise ValueError(f"FoldX paired manifest references a missing batch: {batch}")
                    root = next((candidate for candidate in roots if batch.is_relative_to(candidate)), None)
                    register_batch(batch, root or batch.parent)
                    identity = (batch, int(pair["group"]), int(pair["run"]))
                    if identity in manifest_pairs:
                        raise ValueError(f"Duplicate source paired-model identity: {identity}")
                    manifest_pairs[identity] = pair

    if not batches:
        raise ValueError("FoldX is marked successful but its raw PdbList batches are missing")

    scan_roots = {}
    for batch in batches:
        engine = next((parent for parent in batch.parents if parent.name == "engine_runs"), None)
        if engine:
            scan_roots.setdefault(engine.parent.resolve(), set()).add(batch)
    pair_rows = {}
    for scan_root, expected_batches in scan_roots.items():
        rows, failures = scan_pairs(scan_root)
        if failures:
            raise ValueError(f"Invalid FoldX PdbList pair evidence: {failures[:3]}")
        for row in rows:
            batch = Path(row["batch_path"]).resolve()
            if batch in expected_batches:
                if row.get("ddg_kcal_mol") is None:
                    raise ValueError(f"FoldX pair has no matching Dif score: {batch}/{row['group']}/{row['run']}")
                validate_pair(row)
                pair_rows.setdefault(batch, []).append(row)

    buildmodel_requests, batch_metadata = [], []
    for batch in sorted(batches):
        if not (batch / "individual_list.txt").is_file():
            raise ValueError(f"FoldX batch lacks individual_list.txt: {batch}")
        request_paths = sorted(batch.glob("buildmodel_request.json"))
        requests = [json.loads(path.read_text(encoding="utf-8")) for path in request_paths]
        runs = requests[0].get("number_of_runs") if requests else provenance.get("number_of_runs")
        if not isinstance(runs, int) or runs < 1:
            raise ValueError(f"FoldX expected run count is unknown: {batch}")
        groups = (batch / "individual_list.txt").read_text(encoding="ascii").splitlines()
        rows = pair_rows.get(batch, [])
        found = {(int(row["group"]), int(row["run"])) for row in rows}
        expected = {(group, slot) for group in range(1, len(groups) + 1) for slot in range(runs)}
        if found != expected or len(rows) != len(expected):
            raise ValueError(f"Incomplete FoldX PdbList pairs in {batch}: missing {sorted(expected - found)[:4]}")
        for path, request in zip(request_paths, requests):
            if request.get("command"):
                buildmodel_requests.append({"package_path": (evidence_root / destinations[batch] /
                    path.name).relative_to(output).as_posix(), "command": request["command"],
                    "number_of_runs": request.get("number_of_runs"),
                    "time_seed_rotabase": request.get("time_seed_rotabase"),
                    "worker_sha256": request.get("worker_sha256")})
        target = evidence_root / destinations[batch]
        if target.exists():
            raise ValueError(f"Duplicate FoldX package path: {target}")
        shutil.copytree(batch, target)
        for item in sorted(batch.rglob("*")):
            if item.is_file():
                target_file = target / item.relative_to(batch)
                digest = sha256_file(item)
                if digest != sha256_file(target_file):
                    raise ValueError(f"FoldX evidence changed while copying: {item}")
                copied.append({"path": target_file.relative_to(output).as_posix(),
                               "sha256": digest, "source_sha256": digest,
                               "size_bytes": target_file.stat().st_size})
        batch_metadata.append({"source_batch": str(batch),
                               "package_path": target.relative_to(output).as_posix(),
                               "pair_count": len(rows), "number_of_runs_per_group": runs})

    # Keep the original manifest as an immutable source artifact, and expose a
    # second manifest whose evidence paths resolve entirely within this bundle.
    portable_records = []
    rows_by_group = {}
    for batch, rows in pair_rows.items():
        for row in rows:
            key = tuple(row["mutations"])
            rows_by_group.setdefault(key, []).append(row)

    def portable_pair(row, source_pair=None):
        batch = Path(row["batch_path"]).resolve()
        package_batch = evidence_root / destinations[batch]
        identity = dict((source_pair or {}).get("source_execution_identity") or {})
        identity.setdefault("batch_path", str(batch))
        identity.setdefault("group", int(row["group"]))
        identity.setdefault("run", int(row["run"]))
        result = dict(row)
        if source_pair:
            for key in ("sampling_protocol", "sampling_role", "sampling_cohort", "sampling_slot",
                        "sampling_execution"):
                if key in source_pair:
                    result[key] = source_pair[key]
        result["source_execution_identity"] = identity
        for key in ("wt_path", "mut_path", "pdb_list_path", "ddg_file_path",
                    "buildmodel_request_path"):
            value = row.get(key)
            if value:
                path = Path(value).resolve()
                if not path.is_relative_to(batch) or not path.is_file():
                    raise ValueError(f"FoldX pair path is outside its copied batch: {path}")
                result[key] = (package_batch / path.relative_to(batch)).relative_to(output).as_posix()
        result["batch_path"] = package_batch.relative_to(output).as_posix()
        result["input_path"] = (package_batch / "input.pdb").relative_to(output).as_posix()
        result["mutant_file_path"] = (package_batch / "individual_list.txt").relative_to(output).as_posix()
        return result

    if pair_manifest is not None:
        consumed = set()
        for record in pair_manifest.get("records", []):
            portable_record = dict(record)
            portable_rows = []
            for source_pair in record.get("paired_models", []):
                identity = (Path(source_pair["batch_path"]).resolve(),
                            int(source_pair["group"]), int(source_pair["run"]))
                row = next((candidate for candidate in pair_rows.get(identity[0], [])
                            if (identity[1], identity[2]) ==
                            (int(candidate["group"]), int(candidate["run"]))), None)
                if row is None:
                    raise ValueError(f"Source paired manifest identity is absent from validated PdbList: {identity}")
                consumed.add(identity)
                portable_rows.append(portable_pair(row, source_pair))
            portable_record["paired_models"] = portable_rows
            portable_records.append(portable_record)
        if consumed != set(manifest_pairs):
            raise ValueError("Source paired manifest identities were not fully represented in portable pairs")
    else:
        for mutations, rows in sorted(rows_by_group.items()):
            portable_records.append({
                "mutations": list(mutations),
                "paired_models": [portable_pair(row) for row in sorted(
                    rows, key=lambda item: (item["batch_path"], int(item["group"]), int(item["run"])))],
            })
    portable_pairs_payload = {
        "schema_version": "1.0",
        "source_manifest_status": "RECORDED" if pair_manifest is not None else "RECONSTRUCTED_FROM_VALIDATED_PDBLIST",
        "sampling_slot_contract": pair_manifest.get("sampling_slot_contract") if pair_manifest else None,
        "records": portable_records,
    }
    portable_pairs_path = evidence_root / "portable_pairs.json"
    write_json(portable_pairs_path, portable_pairs_payload)
    copied.append({"path": portable_pairs_path.relative_to(output).as_posix(),
                   "sha256": sha256_file(portable_pairs_path),
                   "source_sha256": sha256_file(portable_pairs_path),
                   "size_bytes": portable_pairs_path.stat().st_size})

    structure = source / "structure_for_energy.pdb"
    if not structure.is_file():
        raise ValueError("FoldX source structure dependency is missing")
    structure_target = evidence_root / "structure_for_energy.pdb"
    structure_target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(structure, structure_target)
    structure_hash = sha256_file(structure)
    if structure_hash != sha256_file(structure_target):
        raise ValueError("FoldX source structure changed while copying")
    copied.append({"path": structure_target.relative_to(output).as_posix(),
                   "sha256": structure_hash, "source_sha256": structure_hash,
                   "size_bytes": structure_target.stat().st_size})
    for name in ("prediction_context.pdb", "residue_mapping.csv"):
        original = source / name
        if original.is_file():
            target = evidence_root / "source_inputs" / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, target)
            digest = sha256_file(original)
            if digest != sha256_file(target):
                raise ValueError(f"FoldX source input changed while copying: {original}")
            copied.append({"path": target.relative_to(output).as_posix(), "sha256": digest,
                           "source_sha256": digest, "size_bytes": target.stat().st_size})

    repair_status = "UNKNOWN_NOT_AVAILABLE"
    engine_execution_status = "UNKNOWN_NOT_ARCHIVED"
    repair_path = provenance.get("repair_cache_path")
    repair_source = Path(repair_path) if repair_path else None
    if repair_source and repair_source.is_file():
        repair_dir = repair_source.parent
        repair_target = evidence_root / "foldx_repair_cache" / repair_dir.parent.name / repair_dir.name
        repair_status = "RECORDED"
        for name in ("structure_Repair.pdb", "structure.pdb", "rotabase.txt",
                     "engine_execution.json", "structure_Repair.fxout"):
            original = repair_dir / name
            if not original.is_file():
                if name in {"engine_execution.json", "structure_Repair.fxout"}:
                    repair_status = "PARTIAL"
                    continue
                raise ValueError(f"FoldX repair-cache evidence is missing: {name}")
            target_file = repair_target / name
            target_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, target_file)
            digest = sha256_file(original)
            if digest != sha256_file(target_file):
                raise ValueError(f"FoldX repair evidence changed while copying: {original}")
            copied.append({"path": target_file.relative_to(output).as_posix(),
                           "sha256": digest, "source_sha256": digest,
                           "size_bytes": target_file.stat().st_size})
            if name == "engine_execution.json":
                engine_execution_status = "RECORDED"

    source_v4 = source / "v4_manifest.json"
    if not source_v4.is_file():
        raise ValueError("FoldX source manifest is missing")
    source_v4_target = evidence_root / "source_v4_manifest.json"
    shutil.copy2(source_v4, source_v4_target)
    source_manifest_hash = sha256_file(source_v4)
    if source_manifest_hash != sha256_file(source_v4_target):
        raise ValueError("FoldX source manifest changed while copying")
    copied.append({"path": source_v4_target.relative_to(output).as_posix(),
                   "sha256": source_manifest_hash, "source_sha256": source_manifest_hash,
                   "size_bytes": source_v4_target.stat().st_size})
    if pair_manifest_path.is_file():
        pair_target = evidence_root / "paired_model_manifest.json"
        shutil.copy2(pair_manifest_path, pair_target)
        pair_hash = sha256_file(pair_manifest_path)
        if pair_hash != sha256_file(pair_target):
            raise ValueError("FoldX paired-model manifest changed while copying")
        copied.append({"path": pair_target.relative_to(output).as_posix(),
                       "sha256": pair_hash, "source_sha256": pair_hash,
                       "size_bytes": pair_target.stat().st_size})
    payload = {
        "schema_version": "1.0", "source_run": str(source),
        "source_v4_manifest_sha256": source_manifest_hash,
        "source_worker": provenance.get("worker"),
        "source_worker_script": {"status": "RECORDED" if any("foldx_invocations/" in row["path"] and "/source/" in row["path"] for row in copied)
                                 else "UNKNOWN_NOT_ARCHIVED_WITH_RUN",
                                 "recorded_path": source_manifest.get("tool_commands", {}).get("foldx", [])[-1:]},
        "recorded_worker_command": source_manifest.get("tool_commands", {}).get("foldx"),
        "actual_foldx_cli": {"status": "RECORDED" if buildmodel_requests else "UNKNOWN_NOT_RECORDED",
                              "buildmodel_requests": buildmodel_requests},
        "worker_provenance": provenance, "repair_cache_status": repair_status,
        "engine_execution_status": engine_execution_status,
        "worker_execution_manifest_status": "RECORDED" if any(
            Path(row["path"]).name == "execution_manifest.json" for row in copied)
            else "UNKNOWN_NOT_ARCHIVED",
        "foldx_parameters": source_manifest.get("parameters", {}).get("joint_v4_ml", {}).get("foldx_confirmation", {}),
        "batches": batch_metadata,
        "portable_pair_manifest_path": portable_pairs_path.relative_to(output).as_posix(),
        "historical_sampling": {"number_of_runs": provenance.get("number_of_runs"),
                                "time_seed_rotabase": [row["time_seed_rotabase"] for row in buildmodel_requests]
                                    if buildmodel_requests else "UNKNOWN_NOT_RECORDED",
                                "new_replay_model_calls": 0},
        "raw_artifacts": copied,
    }
    manifest_path = evidence_root / "foldx_source_evidence.json"
    write_json(manifest_path, payload)
    payload["manifest_path"] = manifest_path.relative_to(output).as_posix()
    payload["manifest_sha256"] = sha256_file(manifest_path)
    payload["artifact_paths"] = [row["path"] for row in copied] + [payload["manifest_path"]]
    return payload

def verify_portable_source_evidence(output: Path, evidence: dict) -> bool:
    """Verify every copied source artifact using only paths inside the replay package."""
    output = Path(output).resolve()

    def package_file(relative):
        value = str(relative)
        posix_path = PurePosixPath(value.replace("\\", "/"))
        windows_path = PureWindowsPath(value)
        if (posix_path.is_absolute() or windows_path.is_absolute() or windows_path.drive or
                ".." in posix_path.parts or not posix_path.parts):
            return None
        resolved = (output / Path(*posix_path.parts)).resolve()
        return resolved if resolved.is_relative_to(output) and resolved.is_file() else None

    try:
        manifest_path = package_file(evidence["manifest_path"])
        if manifest_path is None or sha256_file(manifest_path) != evidence["manifest_sha256"]:
            return False
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        for row in payload["raw_artifacts"]:
            path = package_file(row["path"])
            if path is None or sha256_file(path) != row["sha256"]:
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return False

