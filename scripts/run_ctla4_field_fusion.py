"""Freeze and run the independent CTLA4 field-fusion shadow study."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse import load_npz
from scipy.sparse.csgraph import dijkstra

from ..util import sha256_file, write_json
from ..surface_perturbation.chemistry import CHANNELS
from ..surface_perturbation.engine import direct_site_masks, direct_support, cached_numerical_bounds
from ..surface_perturbation.fields import aligned_atoms, difference, reference_id, spatial_field
from ..surface_perturbation.directions import VECTOR_GROUPS
from ..surface_perturbation.null_store import GEOMETRY_CHANNELS, FIELD_CHANNELS
from ..surface_perturbation.paired import atoms, scan_pairs, validate_pair
from ..surface_perturbation.tools import PhysicalTools, fingerprint, make_grid
from ..surface_perturbation.execution import potential_pair, solvent_distances
from ..surface_perturbation.electrostatics import grid_resolution
from ..surface_perturbation.research_electrostatics import (
    prepare_preserved_model, validate_pair_parameters, parameter_signature,
)
from ..surface_perturbation.field_fusion_research import (
    PROFILE, PREPARATION_PROTOCOL, fuse_scope, geometry_summary, validate_null_roles,
)

DEFAULT_SOURCE = Path(r"E:\user\New project\ctla4_cov2_20261004\routes\coverage_v3\results")
DEFAULT_STUDY = Path(r"E:\user\New project\ctla4_field_fusion_20261006")
PANEL_IDS = (
    "CV2_f1408e6b440bb101", "CV2_9d617576509297ec", "CV2_7d231189c520d563",
    "CV2_cf83338f020be23a", "CV2_02f92b7d78126be5",
)
LABELS = {
    "CV2_f1408e6b440bb101": "A61W/M134R",
    "CV2_9d617576509297ec": "T67R/M134R",
    "CV2_7d231189c520d563": "Q117R/E155W",
    "CV2_cf83338f020be23a": "S55F/S101R",
    "CV2_02f92b7d78126be5": "D78A/M122W",
    "CV2_5baddd4518a94929": "S55F/I102R (coverage extension)",
}
PILOT_IDS = ("CV2_f1408e6b440bb101", "CV2_02f92b7d78126be5")


def _copy(source: Path, target: Path) -> str:
    source, target = Path(source), Path(target)
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if sha256_file(source) != sha256_file(target):
            raise ValueError(f"refusing to replace frozen file with different content: {target}")
    else:
        shutil.copy2(source, target)
    return sha256_file(target)


def _source_snapshot(repo: Path, study: Path) -> dict:
    paths = [repo / "util.py", repo / "models.py", repo / "pyproject.toml", repo / "pytest.ini"]
    paths += sorted((repo / "surface_perturbation").glob("*.py"))
    paths += [repo / "integrations" / "wsl_bridge.py", repo / "scripts" / "run_ctla4_field_fusion.py",
              repo / "tests" / "test_field_fusion_research.py", repo / "tests" / "test_pdb2pqr_research_protocol.py"]
    expected = {source.relative_to(repo).as_posix(): sha256_file(source) for source in paths if source.is_file()}
    snapshot_id = hashlib.sha256(json.dumps(expected, sort_keys=True).encode()).hexdigest()
    target = study / "source_snapshots" / snapshot_id
    target.mkdir(parents=True, exist_ok=True)
    files = {}
    for source in paths:
        if not source.is_file():
            continue
        relative = source.relative_to(repo)
        digest = _copy(source, target / relative)
        files[relative.as_posix()] = digest
    if files != expected:
        raise ValueError("source snapshot file identity changed during freeze")
    return {"id": snapshot_id, "root": str(target), "files": files, "tree_sha256": snapshot_id}


def freeze(source: Path, study: Path) -> dict:
    source, study = Path(source).resolve(), Path(study).resolve()
    report_path = source / "CTLA4_v3_final_report.md"
    official = source / "surface_perturbation"
    pair_path = official / "paired_manifest.json"
    pair_map = json.loads(pair_path.read_text(encoding="utf-8"))
    selected = json.loads((source / "selected_panel_v4.json").read_text(encoding="utf-8"))["selected_candidate_ids"]
    if tuple(selected) != PANEL_IDS:
        raise ValueError("saved official panel identity/order differs from frozen CTLA4 reference")
    for cid in PANEL_IDS:
        if cid not in pair_map:
            raise ValueError(f"selected candidate has no saved paired models: {cid}")
    inputs = study / "inputs"
    for name in ("analysis.json", "structure_for_energy.pdb", "native_input.cif", "selected_panel_v4.json",
                 "CTLA4_v3_final_report.md", "run_manifest.json", "frozen_input_manifest.json"):
        _copy(source / name, inputs / name)
    for name in ("canonical_reference.npz", "reference_graph.npz", "native_WT_SES_mesh.npz",
                 "reference_contract.json", "paired_manifest.json"):
        _copy(official / name, inputs / name)
    _copy(official / "coverage_matrix.npz", inputs / "legacy_coverage_matrix.npz")
    for cid in PANEL_IDS:
        for suffix in ("coverage.npz", "field_audit.json", "electrostatics.npz", "electrostatics.json",
                       "grid_sensitivity.npz"):
            path = official / f"{cid}_{suffix}"
            if path.is_file():
                _copy(path, inputs / "legacy_evidence" / path.name)

    selected_pairs = {}; wt_content_diversity = {}
    batch_map = {}
    for cid in PANEL_IDS:
        rows = pair_map[cid]
        by_role = {}
        for role, expected in (("production", 5), ("null_train", 10), ("null_holdout", 10)):
            chosen = sorted((row for row in rows if row.get("sampling_role") == role),
                            key=lambda row: int(row["sampling_slot"]))
            if len(chosen) != expected or [int(r["sampling_slot"]) for r in chosen] != list(range(expected)):
                raise ValueError(f"{cid}: invalid fixed {role} slots ({len(chosen)})")
            identities = [(r["batch_path"], r["group"], r["run"]) for r in chosen]
            if len(set(identities)) != expected:
                raise ValueError(f"{cid}: repeated execution identity in {role}")
            by_role[role] = chosen
        train_hashes = {r["wt_sha256"] for r in by_role["null_train"]}
        held_hashes = {r["wt_sha256"] for r in by_role["null_holdout"]}
        wt_content_diversity[cid] = {
            "unique_WT_by_role": {role: len({r["wt_sha256"] for r in selected})
                                  for role, selected in by_role.items()},
            "train_holdout_shared_WT_contents": len(train_hashes & held_hashes),
            "production_train_shared_WT_contents": len({r["wt_sha256"] for r in by_role["production"]} & train_hashes),
            "production_holdout_shared_WT_contents": len({r["wt_sha256"] for r in by_role["production"]} & held_hashes),
            "independent_holdout_content": not bool(train_hashes & held_hashes),
            "interpretation": "logical roles remain separate; shared bytes reduce content independence and cannot be counted as new samples",
        }
        frozen_rows = []
        for role in ("production", "null_train", "null_holdout"):
            for original in by_role[role]:
                validate_pair(original)
                row = dict(original)
                original_batch = Path(original["batch_path"])
                batch_key = hashlib.sha256(str(original_batch.resolve()).encode()).hexdigest()[:20]
                target_batch = study / "inputs" / "foldx_batches" / batch_key
                if batch_key not in batch_map:
                    target_batch.mkdir(parents=True, exist_ok=True)
                    names = {"input.pdb", "individual_list.txt"}
                    names.update(p.name for p in original_batch.glob("PdbList_*.fxout"))
                    names.update(p.name for p in original_batch.glob("Dif_*.fxout"))
                    names.update(Path(r[key]).name for r in rows for key in ("wt_path", "mut_path")
                                 if Path(r["batch_path"]).resolve() == original_batch.resolve())
                    for name in sorted(names):
                        if (original_batch / name).is_file():
                            _copy(original_batch / name, target_batch / name)
                    batch_map[batch_key] = str(target_batch)
                row["batch_path"] = str(target_batch)
                for key in ("wt_path", "mut_path", "pdb_list_path", "ddg_file_path"):
                    if original.get(key):
                        row[key] = str(target_batch / Path(original[key]).name)
                validate_pair(row)
                row["sampling_role"] = role
                frozen_rows.append(row)
        selected_pairs[cid] = frozen_rows

    source_snapshot = _source_snapshot(Path(__file__).resolve().parents[1], study)
    manifest = {
        "research_profile": PROFILE, "study_id": "CTLA4-field-fusion-v1-20261006",
        "source_root": str(source), "source_report_sha256": sha256_file(report_path),
        "candidate_ids": list(PANEL_IDS), "candidate_labels": LABELS,
        "pilot_candidate_ids": list(PILOT_IDS), "selected_panel_sha256": sha256_file(source / "selected_panel_v4.json"),
        "paired_models": {cid: {role: sum(row["sampling_role"] == role for row in selected_pairs[cid])
                                 for role in ("production", "null_train", "null_holdout")} for cid in PANEL_IDS},
        "WT_content_diversity": wt_content_diversity,
        "pair_manifest": selected_pairs,
        "input_hashes": {p.relative_to(inputs).as_posix(): sha256_file(p)
                         for p in inputs.rglob("*") if p.is_file()},
        "source_snapshot": source_snapshot,
        "scientific_scope": {
            "official_12A_panel_preserved": True, "new_result_role": "research-only shadow control",
            "expression_binding_calibration": "unavailable; no binding probability or experimental hit-rate claims",
            "surface_denominator": "frozen reliable WT SES points and area weights",
            "direct_site_gate_A": 4.0, "full_surface_query": True,
            "apbs_preparation": PREPARATION_PROTOCOL,
            "apbs_conditions": {"force_field": "PARSE", "pH": 7.4, "temperature_K": 298.15,
                                "salt_molar": 0.15, "ion_valence": 1, "potential_only": True},
            "null_roles": {"production": 5, "estimate_WT": 10, "holdout_WT": 10,
                           "difference_blocks_each": 5},
        },
        "source_snapshot_sha256": source_snapshot["tree_sha256"],
        "status": "FROZEN_T0",
    }
    path = study / "study_manifest.json"
    if path.exists():
        prior = json.loads(path.read_text(encoding="utf-8"))
        if prior.get("candidate_ids") != list(PANEL_IDS) or prior.get("source_report_sha256") != manifest["source_report_sha256"]:
            raise ValueError("existing research manifest has a different frozen input identity")
        # Keep the original frozen manifest; never silently rebase an active study.
        manifest = prior
    else:
        write_json(path, manifest)
    return manifest


def _coordinate_sha256(path: Path) -> str:
    rows = []
    for key, value in sorted(atoms(path).items()):
        xyz = np.rint(np.asarray(value[0], float) * 1000).astype(np.int64).tolist()
        rows.append([*key, value[1], value[2], *xyz])
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def _foldx_short_alias(binary: Path) -> tuple[Path, str]:
    alias_root = Path(r"E:\cfnull_20261006\assets")
    alias_root.mkdir(parents=True, exist_ok=True)
    alias = alias_root / "foldx.exe"
    source_sha = sha256_file(binary)
    if alias.exists() and sha256_file(alias) != source_sha:
        raise ValueError("short FoldX execution alias has different binary content")
    if not alias.exists(): shutil.copy2(binary, alias)
    alias_sha = sha256_file(alias)
    if alias_sha != source_sha: raise ValueError("short FoldX alias verification failed")
    return alias, alias_sha


def _study_manifest(study: Path) -> dict:
    manifest = json.loads((study / "study_manifest.json").read_text(encoding="utf-8"))
    path = study / "evidence" / "null_models" / "independent_null_v1.json"
    if path.is_file():
        bank = json.loads(path.read_text(encoding="utf-8"))
        digest = hashlib.sha256(json.dumps(bank["candidates"], sort_keys=True,
            separators=(",", ":")).encode()).hexdigest()
        if digest != bank.get("candidates_sha256"):
            raise ValueError("independent null bank checksum mismatch")
        manifest["_independent_null_bank"] = bank
    return manifest


def _persist_null_bank(study: Path, checkpoint: dict) -> None:
    candidates = {}
    for cid, state in checkpoint.get("candidates", {}).items():
        roles = {role: value["selected_rows"] for role, value in state.get("roles", {}).items()
                 if value.get("status") == "COMPLETE"}
        if roles:
            candidates[cid] = roles
    digest = hashlib.sha256(json.dumps(candidates, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()
    write_json(study / "evidence" / "null_models" / "independent_null_v1.json", {
        "research_profile": PROFILE, "bank_version": "independent_null_v1",
        "sampling_protocol": "FoldX timeSeedRotabase; candidate-matched WT outputs",
        "split_contract": "10 distinct FoldX execution identities per role; five disjoint paired differences; production, estimate and holdout identities do not overlap; coordinate-content diversity is reported separately",
        "candidates": candidates, "candidates_sha256": digest})


def _collect_null_attempt(attempt_root: Path, mutations: list[str]) -> tuple[list[dict], list[dict]]:
    pairs, errors = scan_pairs(attempt_root)
    expected = tuple(sorted(mutations)); selected = []
    for row in pairs:
        if tuple(sorted(row["mutations"])) != expected:
            continue
        validate_pair(row)
        row["coordinate_sha256"] = _coordinate_sha256(Path(row["wt_path"]))
        selected.append(row)
    selected.sort(key=lambda row: int(row["run"]))
    return selected, errors


def _relocate_null_rows(rows: list[dict], source: Path, target: Path) -> None:
    source = source.resolve(); target = target.resolve()
    for row in rows:
        for key in ("batch_path", "wt_path", "mut_path", "pdb_list_path", "ddg_file_path"):
            value = row.get(key)
            if value:
                row[key] = str(target / Path(value).resolve().relative_to(source))


def _accept_null_rows(rows: list[dict], role: str, cid: str, selected: list[dict], seen: set) -> int:
    if len(selected) >= 10: return 0
    added = 0
    for row in rows:
        row.update(sampling_protocol="time_rotabase_v3", sampling_role=role,
                   sampling_cohort=f"field_fusion:{cid}:{role}")
        execution_path = Path(row["batch_path"]) / "sampling_execution.json"
        if not execution_path.is_file():
            raise ValueError("time-seed execution record missing")
        row["sampling_execution"] = json.loads(execution_path.read_text(encoding="utf-8"))
        if row["sampling_execution"].get("protocol") != "time_rotabase_v3":
            raise ValueError("time-seed protocol identity mismatch")
        identity = tuple(row.get("wt_execution_id", (row["batch_path"], row["group"], row["run"])))
        row["wt_execution_id"] = list(identity)
        if identity in seen: continue
        row["sampling_slot"] = len(selected); selected.append(row); seen.add(identity); added += 1
        if len(selected) == 10: break
    return added


def run_null_sampling(study: Path, candidate_ids: list[str]) -> dict:
    """Create fresh null cohorts; legacy duplicate WT slots remain audit-only."""
    from ..stability import find_engine_inventory, score_many_with_foldx

    study = Path(study).resolve(); manifest = _study_manifest(study)
    if not candidate_ids or set(candidate_ids) - set(manifest["candidate_ids"]):
        raise ValueError("null sampling candidates must be a nonempty subset of the frozen whitelist")
    inventory = find_engine_inventory()
    if not inventory.foldx_path:
        raise RuntimeError("FoldX is required to create independent null models")
    binary = Path(inventory.foldx_path).resolve(); binary_sha = sha256_file(binary)
    execution_binary, execution_binary_sha = _foldx_short_alias(binary)
    checkpoint_path = study / "evidence" / "null_models" / "sampling_checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8")) if checkpoint_path.is_file() else {
        "research_profile": PROFILE, "sampling_protocol": "time_rotabase_v3",
        "foldx_binary_sha256": binary_sha, "candidates": {}}
    if checkpoint.get("foldx_binary_sha256") != binary_sha:
        raise ValueError("FoldX binary changed during null sampling")
    ledger_path = study / "study_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    reports = {}
    for cid in candidate_ids:
        production = _role_rows(manifest, cid, "production")
        if len(production) != 5:
            raise ValueError(f"{cid}: fixed production slots incomplete")
        exemplar = production[0]
        input_path = Path(exemplar["batch_path"]) / "input.pdb"
        if not input_path.is_file() or sha256_file(input_path) != exemplar.get("input_sha256"):
            raise ValueError(f"{cid}: frozen FoldX input identity mismatch")
        mutation_group = []
        for notation in exemplar["mutations"]:
            match = re.fullmatch(r"([A-Z])(-?\d+)([A-Z])", notation)
            if not match: raise ValueError(f"{cid}: invalid mutation notation {notation}")
            wt, position, mutant = match.groups()
            mutation_group.append((wt, int(position), exemplar["chain"], mutant))
        state = checkpoint.setdefault("candidates", {}).setdefault(cid, {
            "mutation": exemplar["mutations"], "input_sha256": exemplar["input_sha256"], "roles": {}})
        if state.get("input_sha256") != exemplar["input_sha256"] or state.get("mutation") != exemplar["mutations"]:
            raise ValueError(f"{cid}: null sampling checkpoint input changed")
        production_ids = {(str(row["batch_path"]), int(row["group"]), int(row["run"])) for row in production}
        train_ids = set()
        for role in ("null_train", "null_holdout"):
            role_state = state.setdefault("roles", {}).setdefault(role, {
                "status": "PENDING", "attempts": [], "selected_rows": []})
            selected = role_state["selected_rows"]
            for row in selected:
                validate_pair(row)
                row["coordinate_sha256"] = row.get("coordinate_sha256") or _coordinate_sha256(Path(row["wt_path"]))
                row["wt_execution_id"] = row.get("wt_execution_id") or [row["batch_path"], row["group"], row["run"]]
            if role_state["status"] == "COMPLETE":
                if len(selected) != 10 or len({tuple(row["wt_execution_id"]) for row in selected}) != 10:
                    raise ValueError(f"{cid}: saved {role} cohort no longer meets unique-run contract")
                if role == "null_train": train_ids = {tuple(row["wt_execution_id"]) for row in selected}
                continue
            excluded = production_ids | (train_ids if role == "null_holdout" else set())
            seen = excluded | {tuple(row["wt_execution_id"]) for row in selected}
            # Recover complete prior attempts and restore all legal slots even when their
            # coordinate content repeats; diversity is reported, not used to erase runs.
            for prior in role_state["attempts"]:
                if prior.get("status") not in {"COMPLETE", "RECOVERED"}: continue
                prior_root = Path(prior["attempt_root"])
                if not prior_root.exists(): continue
                prior_rows, prior_errors = _collect_null_attempt(prior_root, exemplar["mutations"])
                if prior_errors: continue
                _accept_null_rows(prior_rows, role, cid, selected, seen)
            role_state["selected_rows"] = selected
            if len(selected) == 10:
                role_state["status"] = "COMPLETE"
                if role == "null_train": train_ids = {tuple(row["wt_execution_id"]) for row in selected}
                continue
            for prior in role_state["attempts"]:
                if prior.get("status") != "RUNNING": continue
                final_root = Path(prior["attempt_root"])
                staging_root = Path(prior["staging_root"])
                source_root = final_root if final_root.exists() else staging_root
                prior_rows = []; prior_errors = [{"reason": "interrupted FoldX attempt has no recoverable output"}]
                if source_root.exists(): prior_rows, prior_errors = _collect_null_attempt(source_root, exemplar["mutations"])
                if source_root == staging_root and staging_root.exists():
                    final_root.parent.mkdir(parents=True, exist_ok=True); staging_root.rename(final_root)
                    _relocate_null_rows(prior_rows, staging_root, final_root)
                if prior_rows:
                    _accept_null_rows(prior_rows, role, cid, selected, seen)
                prior.update(status="RECOVERED" if len(prior_rows) == 10 and not prior_errors else "FAILED_INTERRUPTED",
                             paired_outputs=len(prior_rows), recovery_errors=prior_errors)
                role_state["selected_rows"] = selected
                write_json(checkpoint_path, checkpoint)
            role_state["status"] = "RUNNING"
            while len(selected) < 10:
                attempt_number = len(role_state["attempts"]) + 1
                attempt_root = study / "evidence" / "null_models" / cid / role / f"attempt_{attempt_number:03d}"
                staging_root = Path(r"E:\cfnull_20261006") / cid[-8:] / role / f"attempt_{attempt_number:03d}"
                attempt = {"attempt_root": str(attempt_root), "staging_root": str(staging_root),
                           "status": "RUNNING", "started_unix": time.time()}
                role_state["attempts"].append(attempt); write_json(checkpoint_path, checkpoint)
                _check_disk(study); started = time.perf_counter(); scores = None
                source_root = attempt_root if attempt_root.exists() else staging_root
                if source_root.exists():
                    rows, errors = _collect_null_attempt(source_root, exemplar["mutations"])
                    if not rows and not errors: errors = [{"reason": "prior FoldX attempt has no paired output"}]
                else:
                    rows = []; errors = []
                    scores = score_many_with_foldx(str(execution_binary), input_path, [mutation_group], staging_root,
                        runs=10, already_repaired=False, max_workers=1, time_seed_rotabase=True)
                    rows, errors = _collect_null_attempt(staging_root, exemplar["mutations"])
                    if not scores or scores[0] is None or len(rows) != 10:
                        errors = errors or [{"reason": f"FoldX returned {len(rows)} of 10 paired WT outputs"}]
                if source_root == staging_root and staging_root.exists() and not attempt_root.exists():
                    attempt_root.parent.mkdir(parents=True, exist_ok=True); staging_root.rename(attempt_root)
                    _relocate_null_rows(rows, staging_root, attempt_root)
                added = 0
                try: added = _accept_null_rows(rows, role, cid, selected, seen)
                except (OSError, ValueError) as exc: errors.append({"reason": f"sampling provenance invalid: {exc}"})
                elapsed = time.perf_counter() - started
                attempt.update(status="COMPLETE" if len(rows) == 10 and not errors else "FAILED",
                    wall_seconds=elapsed, paired_outputs=len(rows), new_unique_WT_coordinates=added,
                    errors=errors, score_median=scores[0] if scores else None)
                role_state["selected_rows"] = selected
                role_state["status"] = "COMPLETE" if len(selected) == 10 else "INCOMPLETE_CONTENT_DIVERSITY"
                charged = ledger.setdefault("foldx_null_sampling", [])
                if not any(item.get("attempt_root") == str(attempt_root) for item in charged):
                    charged.append({"candidate_id": cid, "role": role, "attempt": attempt_number,
                        "attempt_root": str(attempt_root), "staging_root": str(staging_root),
                        "wall_seconds": elapsed, "status": attempt["status"],
                        "paired_outputs": len(rows), "new_unique_WT_coordinates": added,
                        "failed_or_duplicate_outputs_remain_charged": True})
                write_json(ledger_path, ledger); write_json(checkpoint_path, checkpoint)
                _persist_null_bank(study, checkpoint)
                if role_state["status"] == "COMPLETE" or errors or added == 0: break
            if role_state["status"] == "COMPLETE":
                selected[:] = selected[:10]
                if role == "null_train": train_ids = {tuple(row["wt_execution_id"]) for row in selected}
            if role_state["status"] != "COMPLETE": break
        if all(state.get("roles", {}).get(role, {}).get("status") == "COMPLETE"
               for role in ("null_train", "null_holdout")):
            try:
                state["validation"] = validate_null_roles(state["roles"]["null_train"]["selected_rows"],
                    state["roles"]["null_holdout"]["selected_rows"],
                    [{**row, "coordinate_sha256": _coordinate_sha256(Path(row["wt_path"]))} for row in production])
                state["status"] = "COMPLETE"
            except ValueError as exc:
                state["status"] = "INVALID_ROLE_OVERLAP"; state["validation_error"] = str(exc)
        reports[cid] = state
        write_json(checkpoint_path, checkpoint); _persist_null_bank(study, checkpoint)
    return {"research_profile": PROFILE,
        "status": "COMPLETE" if all(v.get("status") == "COMPLETE" for v in reports.values()) else "PARTIAL",
        "candidates": reports, "foldx_binary_sha256": binary_sha,
        "foldx_execution_alias": str(execution_binary), "foldx_execution_alias_sha256": execution_binary_sha}


def _reference(study: Path) -> tuple[dict, list, dict]:
    inputs = study / "inputs"
    records_json = json.loads((inputs / "analysis.json").read_text(encoding="utf-8"))
    records = [SimpleNamespace(**row) for row in records_json["residues"]]
    with np.load(inputs / "canonical_reference.npz", allow_pickle=False) as data:
        reference = {key: data[key].copy() for key in data.files}
    reference["graph"] = load_npz(inputs / "reference_graph.npz").tocsr()
    source_identity = sha256_file(inputs / "structure_for_energy.pdb")
    reference["context_identity"] = fingerprint({"research_profile": PROFILE, "native_sha256": source_identity,
        "official_reference_contract_sha256": sha256_file(inputs / "reference_contract.json")})
    if len(reference["area"]) != reference["graph"].shape[0]:
        raise ValueError("frozen reference and graph dimensions differ")
    return reference, records, atoms(inputs / "structure_for_energy.pdb")


def _atom_identity(atom_map: dict) -> str:
    digest = hashlib.sha256()
    for key in sorted(atom_map):
        xyz, residue, element = atom_map[key]
        digest.update(repr(key).encode()); digest.update(np.asarray(xyz, float).tobytes())
        digest.update(str(residue).encode()); digest.update(str(element).encode())
    return digest.hexdigest()


def _check_disk(path: Path) -> None:
    free = shutil.disk_usage(Path(path).anchor).free
    if free < 10 * 1024**3:
        raise RuntimeError("RESOURCE_STOP: less than 10 GiB free; no new physical solve started")


def _record_tool_delta(ledger: dict, tools: PhysicalTools, calls_before: dict, activity_start: int,
                       ledger_path: Path, node: str) -> None:
    activities = tools.activity[activity_start:]
    calls = {name: int(tools.calls.get(name, 0) - calls_before.get(name, 0))
             for name in tools.calls if name != "cache_hits"}
    seconds = sum(float(row.get("runtime_seconds", 0.0)) for row in activities if row.get("tool"))
    if "cumulative_solver_subprocess_seconds" not in ledger:
        ledger["cumulative_solver_subprocess_seconds"] = float(ledger.pop("cumulative_tool_seconds", 0.0))
    ledger["cumulative_solver_subprocess_seconds"] += seconds
    for name, count in calls.items():
        ledger.setdefault("calls", {})[name] = int(ledger["calls"].get(name, 0)) + count
    if activities or any(calls.values()):
        ledger.setdefault("tool_events", []).append({"node": node, "subprocess_seconds": seconds,
            "all_stage_activity_seconds": sum(float(row.get("runtime_seconds", 0.0)) for row in activities),
            "calls": calls, "activity": activities})
    write_json(ledger_path, ledger)


def _load_npz_verified(path: Path, manifest: Path, contract: dict) -> dict | None:
    if not path.exists() and not manifest.exists():
        return None
    if not path.is_file() or not manifest.is_file():
        raise ValueError(f"incomplete field cache pair: {path}")
    metadata = json.loads(manifest.read_text(encoding="utf-8"))
    if metadata.get("contract") != contract or metadata.get("sha256") != sha256_file(path):
        raise ValueError(f"field cache identity or checksum mismatch: {path}")
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key].copy() for key in data.files}


def _model_field(atom_map: dict, reference: dict, tools: PhysicalTools, cache: Path,
                 ledger: dict, label: str, chunk_size: int = 512,
                 ledger_path: Path | None = None) -> dict:
    _check_disk(cache)
    calls_before = dict(tools.calls); activity_start = len(tools.activity)
    ident = _atom_identity(atom_map)
    contract = {"profile": PROFILE, "atom_map_sha256": ident, "reference_sha256": reference_id(reference),
                "nanoshaper_sha256": tools.identity(tools.nano), "scale": 2.0, "shell_step_a": 1.0,
                "shell_offsets_a": [-2.0, 4.0], "query_points": len(reference["area"]),
                "query_scope": "all_frozen_reliable_surface_samples_v1", "chunk_size": chunk_size}
    path = cache / ("surface_" + fingerprint(contract) + ".npz")
    reused = _load_npz_verified(path, path.with_suffix(".json"), contract)
    if reused is not None:
        ledger.setdefault("nodes", {})[label] = {"status": "REUSED_VERIFIED", "array_sha256": sha256_file(path)}
        if ledger_path is not None: _record_tool_delta(ledger, tools, calls_before, activity_start, ledger_path, label)
        return reused
    started = time.perf_counter(); n = len(reference["area"]); arrays = {}; metadata = {}
    for start in range(0, n, chunk_size):
        ids = np.arange(start, min(n, start + chunk_size)); active = np.zeros(n, bool); active[ids] = True
        part = spatial_field(atom_map, reference, tools, scale=2.0, shell_step=1.0, query_mask=active)
        for key, value in part.items():
            if not isinstance(value, np.ndarray):
                metadata[key] = value
            elif value.ndim and value.shape[0] == n:
                if key not in arrays:
                    arrays[key] = np.empty(value.shape, dtype=value.dtype)
                arrays[key][ids] = value[ids]
        if ledger_path is not None and len(tools.activity) > activity_start:
            _record_tool_delta(ledger, tools, calls_before, activity_start, ledger_path, label + ":solver")
            calls_before = dict(tools.calls); activity_start = len(tools.activity)
    if set(arrays) < {"sdf", "shell_sdf", "shell_valid", "shell_owned", "chemistry", "atomic_inside",
                       "projected_normals", "local_ses_area"}:
        raise ValueError("full surface field returned an incomplete channel set")
    arrays["shell_step"] = np.asarray(1.0)
    _check_disk(path.parent)
    path.parent.mkdir(parents=True, exist_ok=True); np.savez_compressed(path, **arrays)
    digest = sha256_file(path); write_json(path.with_suffix(".json"), {"contract": contract, "sha256": digest, **metadata})
    ledger.setdefault("nodes", {})[label] = {"status": "COMPUTED", "array_sha256": digest,
        "wall_seconds": time.perf_counter() - started, "mesh_closed": metadata.get("mesh_closed"),
        "ses_total_area_a2": metadata.get("ses_total_area_a2")}
    ledger["cumulative_geometry_field_wall_seconds"] = float(
        ledger.get("cumulative_geometry_field_wall_seconds", 0.0)) + ledger["nodes"][label]["wall_seconds"]
    if ledger_path is not None: _record_tool_delta(ledger, tools, calls_before, activity_start, ledger_path, label + ":field_cache")
    return arrays


def _role_rows(manifest: dict, cid: str, role: str) -> list[dict]:
    if role in {"null_train", "null_holdout"} and manifest.get("_independent_null_bank"):
        candidate = manifest["_independent_null_bank"]["candidates"].get(cid)
        if candidate is None or role not in candidate:
            raise ValueError(f"{cid}: fresh independent {role} null cohort is not complete")
        return candidate[role]
    rows = [row for row in manifest["pair_manifest"][cid] if row["sampling_role"] == role]
    rows.sort(key=lambda row: int(row["sampling_slot"]))
    return rows


def _apbs_root(study: Path, manifest: dict) -> Path:
    return (study / "evidence" / "apbs" / "independent_null_v1"
            if manifest.get("_independent_null_bank") else study / "evidence" / "apbs")


def _grid_contains(grid: dict, coordinates: np.ndarray, tolerance: float = 1e-6) -> bool:
    xyz = np.asarray(coordinates, float)
    if xyz.size == 0 or xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        return False
    origin = np.asarray(grid["origin"], float)
    upper = origin + np.asarray(grid["lengths"], float)
    return bool(np.all(xyz >= origin - tolerance) and np.all(xyz <= upper + tolerance))


def _field_pair(pair: dict, native: dict, records: list, reference: dict, tools: PhysicalTools,
                cache: Path, ledger: dict, label: str, ledger_path: Path) -> tuple[dict, dict]:
    wt, mut, alignment = aligned_atoms(pair, native, records)
    if not np.isfinite(alignment["core_rmsd_a"]) or alignment["core_rmsd_a"] > .5:
        raise ValueError("PAIRED_CORE_RMSD_EXCEEDS_0.5A")
    left = _model_field(wt, reference, tools, cache, ledger, label + ":WT", ledger_path=ledger_path)
    right = _model_field(mut, reference, tools, cache, ledger, label + ":Mut", ledger_path=ledger_path)
    delta = difference(left, right)
    return delta, {**left, **delta, "wt_sdf": left["sdf"]}


def _difference_channels(rows: list[dict]) -> dict[str, np.ndarray]:
    names = set.intersection(*(set(row) for row in rows))
    return {name: np.stack([np.asarray(row[name], float) for row in rows]) for name in names}


def _aggregate_numerical_bounds(error_runs: list[dict]) -> dict[str, np.ndarray]:
    numerical = {}
    for channel in FIELD_CHANNELS:
        stack = np.stack([np.asarray(item[channel], float) for item in error_runs])
        count = np.isfinite(stack).sum(axis=0)
        bound = np.max(np.where(np.isfinite(stack), stack, -np.inf), axis=0)
        bound[count < 4] = np.nan
        numerical[channel] = bound
    return numerical


def _weighted_coverages(mask: np.ndarray, area: np.ndarray, graph) -> dict:
    ids = np.flatnonzero(mask)
    if not ids.size:
        distance = np.full(len(area), np.inf)
    else:
        distance = np.asarray(dijkstra(graph, indices=ids, min_only=True), float).reshape(-1)
    values = {}
    denominator = float(area.sum())
    for radius in (6., 9., 12., 15.):
        selected = distance <= radius
        values[str(int(radius))] = {"area_a2": float(area[selected].sum()),
            "fraction": float(area[selected].sum() / denominator) if denominator else None}
    return values


def run_geometry(study: Path, candidate_ids: list[str], *, chunk_size: int = 512) -> dict:
    study = Path(study).resolve(); manifest = _study_manifest(study)
    if manifest.get("research_profile") != PROFILE:
        raise ValueError("wrong or legacy study profile")
    if not candidate_ids or set(candidate_ids) - set(manifest["candidate_ids"]):
        raise ValueError("candidate IDs must be a nonempty subset of the frozen whitelist")
    source_snapshot = _source_snapshot(Path(__file__).resolve().parents[1], study)
    reference, records, native = _reference(study)
    area = np.asarray(reference["area"], float); graph = reference["graph"]
    mesh = PhysicalTools(study / "evidence" / "physical_cache", {
        "wsl_distro": "Ubuntu-22.04", "use_wsl": True,
        "optimization_profile": "coverage_scan_v3_preview", "research_profile": PROFILE,
        "mesh_object_cache_size": 2,
    })
    ledger_path = study / "study_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8")) if ledger_path.exists() else {
        "research_profile": PROFILE, "nodes": {}, "attempts": [], "cumulative_tool_seconds": 0.0,
        "calls": {"nanoshaper": 0, "apbs": 0, "pdb2pqr": 0},
    }
    summaries = {}
    versioned = bool(manifest.get("_independent_null_bank"))
    geometry_root = (study / "evidence" / "geometry" / "independent_null_v1"
                     if versioned else study / "evidence" / "geometry")
    for cid in candidate_ids:
        started = time.perf_counter(); out = geometry_root / cid; out.mkdir(parents=True, exist_ok=True)
        production = _role_rows(manifest, cid, "production")
        estimates = _role_rows(manifest, cid, "null_train")
        held = _role_rows(manifest, cid, "null_holdout")
        if (len(production), len(estimates), len(held)) != (5, 10, 10):
            raise ValueError(f"{cid}: frozen slot contract incomplete")
        null_identity_audit = None
        if versioned:
            production_identity_rows = [{**row, "coordinate_sha256": _coordinate_sha256(Path(row["wt_path"]))}
                                        for row in production]
            null_identity_audit = validate_null_roles(estimates, held, production_identity_rows)
        ledger.setdefault("attempts", []).append({"phase": "geometry", "candidate_id": cid, "status": "RUNNING"})
        write_json(ledger_path, ledger)
        prod_rows = []; prod_left = []; support = []; sites = []
        for row in production:
            validate_pair(row)
            delta, left = _field_pair(row, native, records, reference, mesh,
                                      study / "evidence" / "geometry" / "model_fields", ledger,
                                      f"{cid}:production:{row['sampling_slot']}", ledger_path)
            prod_rows.append(delta); prod_left.append(left)
            support.append(direct_support(row, native, records, reference))
            sites.append(direct_site_masks(row, native, records, reference))
        estimate_models = []; held_models = []
        for role, rows, target in (("null_train", estimates, estimate_models), ("null_holdout", held, held_models)):
            for row in rows:
                validate_pair(row)
                wt, _mut, alignment = aligned_atoms(row, native, records)
                if not np.isfinite(alignment["core_rmsd_a"]) or alignment["core_rmsd_a"] > .5:
                    raise ValueError(f"{cid}: {role} core alignment failed")
                role_label = f"independent_null_v1:{role}" if versioned else role
                target.append(_model_field(wt, reference, mesh, study / "evidence" / "geometry" / "model_fields",
                                          ledger, f"{cid}:{role_label}:{row['sampling_slot']}:WT", chunk_size, ledger_path))
        estimate_diffs = [difference(estimate_models[i], estimate_models[i + 1]) for i in range(0, 10, 2)]
        held_diffs = [difference(held_models[i], held_models[i + 1]) for i in range(0, 10, 2)]
        estimate_channels = _difference_channels(estimate_diffs); held_channels = _difference_channels(held_diffs)
        error_runs = []; numerical_grid_audit = []
        for index, row in enumerate(production):
            _check_disk(study)
            before_calls = dict(mesh.calls); activity_at = len(mesh.activity)
            try:
                values, audit = cached_numerical_bounds(row, native, records, reference, mesh, prod_left[index],
                    study / "evidence" / "geometry" / "numerical", required_mask=np.ones(len(area), bool),
                    return_audit=True)
                error_runs.append(values)
                numerical_grid_audit.append({"production_slot": int(row["sampling_slot"]), **audit})
            finally:
                _record_tool_delta(ledger, mesh, before_calls, activity_at, ledger_path,
                                   f"{cid}:numerical:{row['sampling_slot']}")
        numerical = _aggregate_numerical_bounds(error_runs)
        prod_channels = _difference_channels(prod_rows)
        support_runs = np.asarray(support, bool)
        direct = support_runs.sum(axis=0) >= 4
        site_runs = np.asarray(sites, bool)
        if site_runs.shape[:2] != (5, 2):
            raise ValueError(f"{cid}: direct site support identity mismatch")
        site_quorum = site_runs.sum(axis=0) >= 4
        if not site_quorum.any(axis=1).all():
            raise ValueError(f"{cid}: saved production slots no longer meet per-site direct contribution rule")

        scopes = {
            "geometry": (GEOMETRY_CHANNELS, {"surface_normal_delta": ("surface_normal_x", "surface_normal_y", "surface_normal_z")}),
            "chemistry": (CHANNELS, {k: v for k, v in {
                "donor_delta": ("donor_x", "donor_y", "donor_z"),
                "acceptor_delta": ("acceptor_x", "acceptor_y", "acceptor_z"),
                "aromatic_delta": ("aromatic_x", "aromatic_y", "aromatic_z"),
            }.items()}),
            "geometry_chemistry": ((*GEOMETRY_CHANNELS, *CHANNELS), None),
        }
        report_scopes = {}; output_arrays = {"area": area, "direct_support_4a": direct,
            "direct_support_runs": support_runs, "direct_site_support_runs": site_runs,
            "direct_site_quorum": site_quorum}
        raw_arrays = {f"production__{name}": value for name, value in prod_channels.items()}
        raw_arrays.update({f"estimate__{name}": value for name, value in estimate_channels.items()})
        raw_arrays.update({f"heldout__{name}": value for name, value in held_channels.items()})
        raw_arrays.update({f"numerical__{name}": value for name, value in numerical.items()})
        raw_path = out / "geometry_raw_fields.npz"
        _check_disk(out); np.savez_compressed(raw_path, **raw_arrays)
        raw_manifest = {"research_profile": PROFILE, "candidate_id": cid, "reference_id": reference_id(reference),
                        "production_slots": 5, "estimate_WT_slots": 10, "estimate_blocks": 5,
                        "heldout_WT_slots": 10, "heldout_blocks": 5,
                        "pairing_rule": "fixed ascending slots (0-1,2-3,4-5,6-7,8-9); independent null v1 validates ten distinct FoldX execution identities per role and reports coordinate-content diversity separately",
                        "null_identity_audit": null_identity_audit,
                        "numerical_grid_audit": numerical_grid_audit,
                        "sha256": sha256_file(raw_path), "channels": sorted(prod_channels)}
        write_json(out / "geometry_raw_fields.json", raw_manifest)
        for scope, (names, vector_groups) in scopes.items():
            production_scope = {k: prod_channels[k] for k in names}
            estimate_scope = {k: estimate_channels[k] for k in names}
            held_scope = {k: held_channels[k] for k in names}
            error_scope = {k: numerical[k] for k in names}
            vectors = VECTOR_GROUPS if vector_groups is None else vector_groups
            vector_components = {component for group in vectors.values() for component in group}
            scalar_names = tuple(k for k in names if k not in vector_components)
            fusion = fuse_scope(production_scope, estimate_scope, held_scope, error_scope, area,
                                scalar_names=scalar_names, vector_groups=vector_groups, protocol=PROFILE)
            mask = fusion["positive"]
            summary = geometry_summary(mask, area, graph, direct)
            summary.update(status=fusion["status"], known_area_a2=float(area[fusion["known"]].sum()),
                           unknown_area_a2=float(area[fusion["unknown"]].sum()),
                           nonlocal_positive_area_a2=float(area[mask & ~direct].sum()),
                           within_4a_positive_area_a2=float(area[mask & direct].sum()),
                           scan_opportunity_scenarios=_weighted_coverages(mask, area, graph),
                           channels=fusion["channels"], channel_audits=fusion["channel_audits"])
            summary["unknown_reason_area_a2"] = {
                reason: float(area[mask].sum()) for reason, mask in fusion["unknown_reason_masks"].items()}
            report_scopes[scope] = summary
            output_arrays[f"{scope}__positive"] = mask
            output_arrays[f"{scope}__unknown"] = fusion["unknown"]
            output_arrays[f"{scope}__known"] = fusion["known"]
            for name, value in fusion["channel_positive"].items(): output_arrays[f"{scope}__positive__{name}"] = value
            for name, value in fusion["channel_negative"].items(): output_arrays[f"{scope}__negative__{name}"] = value
            for name, value in fusion["channel_unknown"].items(): output_arrays[f"{scope}__unknown__{name}"] = value
            for name, reasons in fusion["channel_unknown_reasons"].items():
                for reason, value in reasons.items():
                    output_arrays[f"{scope}__unknown_reason__{name}__{reason}"] = value
            for reason, value in fusion["unknown_reason_masks"].items():
                output_arrays[f"{scope}__composite_unknown_reason__{reason}"] = value
        np.savez_compressed(out / "geometry_fusion_masks.npz", **output_arrays)
        result = {
            "research_profile": PROFILE, "candidate_id": cid, "label": LABELS[cid],
            "source_snapshot_sha256": source_snapshot["tree_sha256"],
            "status": "GEOMETRY_CHEMISTRY_COMPLETE", "samples": len(area),
            "reference_id": reference_id(reference), "fixed_denominator_area_a2": float(area.sum()),
            "production_pairs": 5, "estimate_WT_slots": 10, "estimate_difference_blocks": 5,
            "heldout_WT_slots": 10, "heldout_difference_blocks": 5,
            "unique_WT_content": {role: len({r["wt_sha256"] for r in rows})
                                   for role, rows in (("production", production), ("estimate", estimates), ("heldout", held))},
            "null_identity_audit": null_identity_audit,
            "null_bank_candidates_sha256": manifest.get("_independent_null_bank", {}).get("candidates_sha256"),
            "scope_results": report_scopes,
            "raw_field_bundle_sha256": raw_manifest["sha256"],
            "apbs": {"status": "NOT_RUN_IN_GEOMETRY_STAGE", "electrostatic_claims": "UNKNOWN"},
            "legacy_official_panel": {"coverage": 0.8602623530512513, "scan_kernel_radius_a": 12.0,
                                      "interpretation": "preserved reference, not overwritten"},
            "wall_seconds": time.perf_counter() - started,
            "model_cache_node_count": len(ledger.get("nodes", {})),
        }
        write_json(out / "geometry_fusion_report.json", result)
        summaries[cid] = result
        ledger["attempts"][-1].update(status="COMPLETE", wall_seconds=result["wall_seconds"])
        write_json(ledger_path, ledger)
    return summaries


def _residue_ids(atom_map: dict, native: dict) -> dict:
    residues = {key[:3] for key in atom_map}
    order = list(dict.fromkeys(key[:3] for key in native if key[:3] in residues))
    if set(order) != residues:
        raise ValueError("PQR author-residue chronology cannot be mapped uniquely")
    return {key: index for index, key in enumerate(order, 1)}


def _prepare_pair(row: dict, cid: str, role: str, slot: int, native: dict, records: list,
                  reference: dict, tools: PhysicalTools, study: Path, ledger: dict,
                  ledger_path: Path) -> tuple[dict, dict, dict, dict]:
    validate_pair(row)
    wt, mut, alignment = aligned_atoms(row, native, records)
    if not np.isfinite(alignment["core_rmsd_a"]) or alignment["core_rmsd_a"] > .5:
        raise ValueError("PAIRED_CORE_RMSD_EXCEEDS_0.5A")
    residue_ids = _residue_ids(wt, native)
    cache = _apbs_root(study, _study_manifest(study)) / "prepared_models"
    prepared = []
    for kind, atom_map in (("WT", wt), ("Mut", mut)):
        before_calls = dict(tools.calls); activity_start = len(tools.activity)
        item = prepare_preserved_model(atom_map, residue_ids, tools, cache)
        _record_tool_delta(ledger, tools, before_calls, activity_start, ledger_path,
                           f"{cid}:{role}:{slot}:PDB2PQR:{kind}")
        prepared.append(item)
    mutated_pqr = {(chain, int(residue_ids[(chain, position, insertion)]))
                   for chain, position, insertion in alignment["mutated_author_residues"]}
    parameter_audit = validate_pair_parameters(prepared[0]["atoms"], prepared[1]["atoms"], mutated_pqr)
    valid = (all(item["audit"]["status"] == "MATCHED" for item in prepared)
             and parameter_audit["status"] == "MATCHED")
    audit = {"candidate_id": cid, "role": role, "slot": slot,
             "pair_identity": {"wt_sha256": row["wt_sha256"], "mut_sha256": row["mut_sha256"],
                               "execution": [row["batch_path"], row["group"], row["run"]]},
             "heavy_atom_preparation": {kind: item["audit"] for kind, item in zip(("WT", "Mut"), prepared)},
             "nonmutated_parameter_audit": parameter_audit, "preparation_valid": bool(valid),
             "pqr_sha256": {kind: item["contract"]["prepared_pqr_sha256"]
                            for kind, item in zip(("WT", "Mut"), prepared)}}
    audit_dir = _apbs_root(study, _study_manifest(study)) / cid / "preparation"
    write_json(audit_dir / f"{role}_{slot:02d}.json", audit)
    return prepared[0], prepared[1], alignment, audit


def _prepare_null(row: dict, cid: str, role: str, slot: int, native: dict, records: list,
                  tools: PhysicalTools, study: Path, ledger: dict, ledger_path: Path) -> tuple[dict, dict]:
    validate_pair(row)
    wt, _mut, alignment = aligned_atoms(row, native, records)
    if not np.isfinite(alignment["core_rmsd_a"]) or alignment["core_rmsd_a"] > .5:
        raise ValueError(f"{role} core alignment exceeds 0.5 A")
    residue_ids = _residue_ids(wt, native)
    before_calls = dict(tools.calls); activity_start = len(tools.activity)
    prepared = prepare_preserved_model(wt, residue_ids, tools,
        _apbs_root(study, _study_manifest(study)) / "prepared_models")
    _record_tool_delta(ledger, tools, before_calls, activity_start, ledger_path,
                       f"{cid}:{role}:{slot}:PDB2PQR:WT")
    audit = {"candidate_id": cid, "role": role, "slot": slot, "wt_sha256": row["wt_sha256"],
             "execution": [row["batch_path"], row["group"], row["run"]],
             "heavy_atom_preparation": prepared["audit"],
             "parameter_signature": parameter_signature(prepared["atoms"]),
             "pqr_sha256": prepared["contract"]["prepared_pqr_sha256"]}
    return prepared, audit


def _fixed_apbs_grids(study: Path, manifest: dict, native: dict, records: list, reference: dict) -> tuple[dict, dict]:
    rows_digest = hashlib.sha256()
    xyz = [value[0] for value in native.values()]
    for cid in manifest["candidate_ids"]:
        production = _role_rows(manifest, cid, "production")
        for row in production:
            wt, mut, audit = aligned_atoms(row, native, records)
            if not np.isfinite(audit["core_rmsd_a"]) or audit["core_rmsd_a"] > .5:
                raise ValueError(f"{cid}: global APBS grid input alignment failed")
            xyz.extend(value[0] for value in wt.values()); xyz.extend(value[0] for value in mut.values())
            rows_digest.update(row["wt_sha256"].encode()); rows_digest.update(row["mut_sha256"].encode())
    coords = np.asarray(xyz, float)
    base = make_grid(coords, .5, 16.)
    variants = [make_grid(coords, spacing, margin) for spacing, margin in ((.4, 16.), (.6, 16.), (.5, 24.))]
    # ponytail: one frozen production-derived global grid; new null cohorts must fit its envelope.
    null_bank = manifest.get("_independent_null_bank")
    if null_bank:
        null_candidates = null_bank.get("candidates", {})
        roles_to_check = ((role, rows) for candidate in null_candidates.values()
                          for role, rows in candidate.items() if role in {"null_train", "null_holdout"})
    else:
        roles_to_check = ((role, _role_rows(manifest, cid, role)) for cid in manifest["candidate_ids"]
                          for role in ("null_train", "null_holdout"))
    for role, rows in roles_to_check:
        for row in rows:
            wt, _mut, audit = aligned_atoms(row, native, records)
            if not np.isfinite(audit["core_rmsd_a"]) or audit["core_rmsd_a"] > .5:
                raise ValueError(f"{role}: null APBS grid input alignment failed")
            null_xyz = np.asarray([value[0] for value in wt.values()], float)
            if not all(_grid_contains(grid, null_xyz) for grid in [base, *variants]):
                raise ValueError(f"{role}: null WT atoms extend beyond the frozen APBS sensitivity-grid envelope")
    contract = {"profile": PROFILE, "reference_sha256": reference_id(reference),
                "production_pair_inputs_sha256": rows_digest.hexdigest(),
                "candidate_ids": manifest["candidate_ids"], "global_grid": True,
                "grid_basis": "native plus all five frozen production WT/Mut pairs per panel candidate; null coordinates are containment-checked",
                "production_grid": {"spacing_a": .5, "margin_a": 16.},
                "sensitivity_grids": [{"spacing_a": .4, "margin_a": 16.},
                                       {"spacing_a": .6, "margin_a": 16.},
                                       {"spacing_a": .5, "margin_a": 24.}]}
    path = _apbs_root(study, manifest) / "apbs_grid_contract.json"
    if path.exists():
        saved = json.loads(path.read_text(encoding="utf-8"))
        if saved.get("contract") != contract:
            prior = saved.get("contract", {})
            prior_ids = prior.get("candidate_ids", [])
            appended_ids = contract["candidate_ids"][len(prior_ids):]
            stable_contract = {key: value for key, value in contract.items()
                               if key not in {"candidate_ids", "production_pair_inputs_sha256"}}
            prior_stable_contract = {key: value for key, value in prior.items()
                                     if key not in {"candidate_ids", "production_pair_inputs_sha256"}}
            same_grid = (saved.get("base_grid") == base
                         and saved.get("sensitivity_grids") == variants)
            is_append_only = (bool(appended_ids)
                              and contract["candidate_ids"][:len(prior_ids)] == prior_ids)
            if not (is_append_only and stable_contract == prior_stable_contract and same_grid):
                raise ValueError("global APBS grid contract changed after first solve; panel extension requires identical grids")
            prior_grid_identity = hashlib.sha256(json.dumps(
                {"base_grid": saved["base_grid"], "sensitivity_grids": saved["sensitivity_grids"]},
                sort_keys=True, separators=(",", ":")).encode()).hexdigest()
            saved.setdefault("contract_history", []).append({
                "contract": prior,
                "grid_geometry_sha256": prior_grid_identity,
                "appended_candidate_ids": appended_ids,
                "reason": "append-only panel extension; recomputed global grid geometry is byte-equivalent",
            })
            saved["contract"] = contract
            write_json(path, saved)
        grids = (saved["base_grid"], saved["sensitivity_grids"])
        if not all(_grid_contains(grid, np.asarray(xyz, float)) for grid in [grids[0], *grids[1]]):
            raise ValueError("saved APBS global grid no longer contains its frozen production inputs")
        return grids
    write_json(path, {"contract": contract, "base_grid": base, "sensitivity_grids": variants,
                      "reference_sample_count": len(reference["area"]),
                      "interpretation": "grid fixed from all panel production pairs; all null WT coordinates must pass containment"})
    return base, variants


def _interpolate(tools: PhysicalTools, pqr_path: Path, grid: dict, queries: np.ndarray,
                 ledger: dict, ledger_path: Path, node: str):
    _check_disk(ledger_path.parent)
    before_calls = dict(tools.calls); activity_start = len(tools.activity)
    try:
        field, _same = potential_pair(tools, pqr_path, pqr_path, grid)
        return np.asarray(field(queries), float)
    finally:
        _record_tool_delta(ledger, tools, before_calls, activity_start, ledger_path, node)


def _sample_solvent(tools: PhysicalTools, atom_map: dict, queries: np.ndarray,
                    ledger: dict, ledger_path: Path, node: str):
    before_calls = dict(tools.calls); activity_start = len(tools.activity)
    try:
        return solvent_distances(tools, atom_map, queries, return_audit=True)
    finally:
        _record_tool_delta(ledger, tools, before_calls, activity_start, ledger_path, node)


def run_apbs(study: Path, candidate_ids: list[str]) -> dict:
    study = Path(study).resolve(); manifest = _study_manifest(study)
    if not candidate_ids or set(candidate_ids) - set(manifest["candidate_ids"]):
        raise ValueError("APBS candidate IDs must be a nonempty subset of the frozen whitelist")
    for cid in candidate_ids:
        _role_rows(manifest, cid, "null_train"); _role_rows(manifest, cid, "null_holdout")
    apbs_root = _apbs_root(study, manifest)
    source_snapshot = _source_snapshot(Path(__file__).resolve().parents[1], study)
    reference, records, native = _reference(study); area = np.asarray(reference["area"], float)
    queries = np.asarray(reference["points"], float) + 1.5 * np.asarray(reference["normals"], float)
    grids = _fixed_apbs_grids(study, manifest, native, records, reference)
    tools = PhysicalTools(apbs_root / "physical_cache", {
        "wsl_distro": "Ubuntu-22.04", "use_wsl": True, "pdb2pqr_version3": False,
        "pdb2pqr_support_paths": ["/usr/share/pdb2pqr", "/usr/lib/pdb2pqr"],
        "pdb2pqr_support_python": "/usr/bin/python3", "apbs_workers": 2,
        "research_profile": PROFILE,
    })
    ledger_path = study / "study_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8")) if ledger_path.exists() else {
        "research_profile": PROFILE, "nodes": {}, "attempts": [], "calls": {"nanoshaper": 0, "apbs": 0, "pdb2pqr": 0}}
    apbs_sha = tools.identity(tools.apbs); pqr_sha = tools.identity(tools.pqr)
    report = {}
    for cid in candidate_ids:
        started = time.perf_counter(); out = apbs_root / cid; out.mkdir(parents=True, exist_ok=True)
        role_rows = {role: _role_rows(manifest, cid, role)
                     for role in ("production", "null_train", "null_holdout")}
        slot_contract = {role: [{"slot": int(row.get("sampling_slot", index)),
                                 "execution": [row["batch_path"], row["group"], row["run"]],
                                 "wt_sha256": row["wt_sha256"], "mut_sha256": row.get("mut_sha256"),
                                 "coordinate_sha256": row.get("coordinate_sha256"),
                                 "sampling_protocol": row.get("sampling_protocol")}
                                for index, row in enumerate(rows)]
                         for role, rows in role_rows.items()}
        slot_identity_sha256 = hashlib.sha256(json.dumps(slot_contract, sort_keys=True,
            separators=(",", ":")).encode()).hexdigest()
        fusion_protocol = hashlib.sha256(json.dumps({"profile": PROFILE, "candidate_id": cid,
            "reference": reference_id(reference), "grid": grids[0], "pdb2pqr": pqr_sha, "apbs": apbs_sha,
            "preparation": PREPARATION_PROTOCOL, "offset_a": 1.5,
            "evidence_slots_sha256": slot_identity_sha256}, sort_keys=True).encode()).hexdigest()
        ledger.setdefault("attempts", []).append({"phase": "apbs", "candidate_id": cid, "status": "RUNNING"})
        write_json(ledger_path, ledger)
        production = role_rows["production"]
        prepared_prod = []; production_audits = []
        for row in production:
            wt, mut, _alignment, audit = _prepare_pair(row, cid, "production", int(row["sampling_slot"]),
                native, records, reference, tools, study, ledger, ledger_path)
            prepared_prod.append((wt, mut)); production_audits.append(audit)

        raw_delta = np.full((5, len(area)), np.nan); solvent_valid = np.zeros((5, len(area)), bool)
        production_solvent_audits = []
        pair_preparation_valid = np.zeros(5, bool); grid_error = np.full((5, len(area)), np.nan)
        grid_valid = np.zeros((5, len(area)), bool); raw_grid_error = np.full((5, len(area), 3), np.nan)
        for index, (wt, mut) in enumerate(prepared_prod):
            _check_disk(study)
            before_calls = dict(tools.calls); activity_start = len(tools.activity)
            try:
                wt_field, mut_field = potential_pair(tools, wt["path"], mut["path"], grids[0])
                raw_delta[index] = mut_field(queries) - wt_field(queries)
            finally:
                _record_tool_delta(ledger, tools, before_calls, activity_start, ledger_path,
                                   f"{cid}:production:{index}:APBS:base")
            wdist, wvalid, waudit = _sample_solvent(tools, wt["atoms"], queries, ledger, ledger_path,
                                             f"{cid}:production:{index}:solvent:WT")
            mdist, mvalid, maudit = _sample_solvent(tools, mut["atoms"], queries, ledger, ledger_path,
                                             f"{cid}:production:{index}:solvent:Mut")
            solvent_valid[index] = wvalid & mvalid & (wdist > 0) & (mdist > 0)
            production_solvent_audits.append({"slot": index, "WT": waudit, "Mut": maudit,
                "valid_area_a2": float(area[solvent_valid[index]].sum()),
                "unknown_area_a2": float(area[~solvent_valid[index]].sum())})
            pair_preparation_valid[index] = bool(production_audits[index]["preparation_valid"])
            variant_delta = []
            for variant_index, grid in enumerate(grids[1]):
                before_calls = dict(tools.calls); activity_start = len(tools.activity)
                try:
                    wfield, mfield = potential_pair(tools, wt["path"], mut["path"], grid)
                    variant_delta.append(mfield(queries) - wfield(queries))
                finally:
                    _record_tool_delta(ledger, tools, before_calls, activity_start, ledger_path,
                        f"{cid}:production:{index}:APBS:grid_variant:{variant_index}")
            local_error, local_valid, local_errors = grid_resolution(raw_delta[index], np.asarray(variant_delta))
            grid_error[index] = local_error; grid_valid[index] = local_valid
            raw_grid_error[index] = local_errors.T
        production_valid = (solvent_valid & grid_valid & pair_preparation_valid[:, None]
                            & np.isfinite(raw_delta))
        production_delta = raw_delta.copy(); production_delta[~production_valid] = np.nan
        estimate_null = []; holdout_null = []; null_audits = {"null_train": [], "null_holdout": []}
        null_solvent_audits = {"null_train": [], "null_holdout": []}
        null_preparations = {}
        for role, target in (("null_train", estimate_null), ("null_holdout", holdout_null)):
            slots = []
            for row in _role_rows(manifest, cid, role):
                prepared, audit = _prepare_null(row, cid, role, int(row["sampling_slot"]),
                    native, records, tools, study, ledger, ledger_path)
                values = _interpolate(tools, prepared["path"], grids[0], queries, ledger, ledger_path,
                                      f"{cid}:{role}:{row['sampling_slot']}:APBS:base")
                distance, solvent_ok, solvent_audit = _sample_solvent(tools, prepared["atoms"], queries, ledger,
                    ledger_path, f"{cid}:{role}:{row['sampling_slot']}:solvent")
                valid = solvent_ok & (distance > 0) & np.isfinite(values)
                values[~valid] = np.nan
                slots.append((values, valid, prepared, audit))
                null_audits[role].append(audit)
                null_solvent_audits[role].append({"slot": int(row["sampling_slot"]), **solvent_audit,
                    "valid_area_a2": float(area[valid].sum()),
                    "unknown_area_a2": float(area[~valid].sum())})
            null_preparations[role] = slots
            for left_index in range(0, 10, 2):
                left, right = slots[left_index], slots[left_index + 1]
                diff = left[0] - right[0]
                signature_match = parameter_signature(left[2]["atoms"]) == parameter_signature(right[2]["atoms"])
                prep_match = (left[2]["audit"]["status"] == "MATCHED" and right[2]["audit"]["status"] == "MATCHED")
                valid = left[1] & right[1] & bool(signature_match and prep_match)
                diff[~valid] = np.nan
                target.append(diff)
        estimate_blocks = np.stack(estimate_null); heldout_blocks = np.stack(holdout_null)
        numerical_stack = grid_error.copy(); numerical_count = np.isfinite(numerical_stack).sum(axis=0)
        apbs_numerical = np.max(np.where(np.isfinite(numerical_stack), numerical_stack, -np.inf), axis=0)
        apbs_numerical[numerical_count < 4] = np.nan
        scalar_prod = {"delta_phi": production_delta}
        scalar_est = {"delta_phi": estimate_blocks}
        scalar_hold = {"delta_phi": heldout_blocks}
        scalar_error = {"delta_phi": apbs_numerical}
        electro_fusion = fuse_scope(scalar_prod, scalar_est, scalar_hold, scalar_error, area,
            scalar_names=("delta_phi",), vector_groups={}, protocol=fusion_protocol)

        geometry_root = (study / "evidence" / "geometry" / "independent_null_v1"
                         if manifest.get("_independent_null_bank") else study / "evidence" / "geometry")
        geometry_path = geometry_root / cid / "geometry_raw_fields.npz"
        geometry_meta = json.loads(geometry_path.with_suffix(".json").read_text(encoding="utf-8"))
        if (geometry_meta.get("candidate_id") != cid or geometry_meta.get("research_profile") != PROFILE
                or geometry_meta.get("reference_id") != reference_id(reference)
                or sha256_file(geometry_path) != geometry_meta["sha256"]):
            raise ValueError(f"{cid}: full-surface geometry raw bundle checksum mismatch")
        with np.load(geometry_path, allow_pickle=False) as data:
            geometry = {key: data[key].copy() for key in data.files}
        all_names = [name for name in FIELD_CHANNELS]
        production_fields = {name: geometry[f"production__{name}"] for name in all_names}
        estimate_fields = {name: geometry[f"estimate__{name}"] for name in all_names}
        heldout_fields = {name: geometry[f"heldout__{name}"] for name in all_names}
        error_fields = {name: geometry[f"numerical__{name}"] for name in all_names}
        production_fields["delta_phi"] = production_delta
        estimate_fields["delta_phi"] = estimate_blocks
        heldout_fields["delta_phi"] = heldout_blocks
        error_fields["delta_phi"] = apbs_numerical
        vector_components = {part for group in VECTOR_GROUPS.values() for part in group}
        three_channel = fuse_scope(production_fields, estimate_fields, heldout_fields, error_fields, area,
            scalar_names=tuple(name for name in (*FIELD_CHANNELS, "delta_phi") if name not in vector_components),
            vector_groups=VECTOR_GROUPS, protocol=fusion_protocol)
        with np.load(geometry_root / cid / "geometry_fusion_masks.npz",
                     allow_pickle=False) as mask_data:
            direct_seed = mask_data["direct_support_4a"].copy()
        summaries = {
            "electrostatic_only": geometry_summary(electro_fusion["positive"], area, reference["graph"]),
            "geometry_chemistry": json.loads((geometry_root / cid / "geometry_fusion_report.json").read_text(encoding="utf-8"))["scope_results"]["geometry_chemistry"],
            "three_mechanism_fusion": geometry_summary(three_channel["positive"], area, reference["graph"]),
        }
        for key, fusion in (("electrostatic_only", electro_fusion), ("three_mechanism_fusion", three_channel)):
            summaries[key].update(status=fusion["status"], known_area_a2=float(area[fusion["known"]].sum()),
                unknown_area_a2=float(area[fusion["unknown"]].sum()),
                scan_opportunity_scenarios=_weighted_coverages(fusion["positive"], area, reference["graph"]),
                channels=fusion["channels"], channel_audits=fusion["channel_audits"])
            summaries[key]["unknown_reason_area_a2"] = {
                reason: float(area[mask].sum()) for reason, mask in fusion["unknown_reason_masks"].items()}
        summaries["electrostatic_only"]["distance_interpretation"] = "electrostatic field change; not a geodesic propagation radius"
        summaries["three_mechanism_fusion"]["distance_interpretation"] = "mixed mechanisms; use geometry_chemistry for direct-seed geodesic summaries"
        result_arrays = {"area": area, "production_delta_phi": raw_delta,
            "production_qualified_delta_phi": production_delta, "production_solvent_valid": solvent_valid,
            "production_grid_converged": grid_valid, "production_grid_error": grid_error,
            "production_grid_error_variants": raw_grid_error, "production_preparation_valid": pair_preparation_valid,
            "estimate_null_delta_phi": estimate_blocks, "holdout_null_delta_phi": heldout_blocks,
            "apbs_numerical_error": apbs_numerical,
            "electrostatic_positive": electro_fusion["positive"], "electrostatic_unknown": electro_fusion["unknown"],
            "three_mechanism_positive": three_channel["positive"], "three_mechanism_unknown": three_channel["unknown"]}
        result_arrays["electrostatic_negative"] = electro_fusion["known"] & ~electro_fusion["positive"]
        result_arrays["three_mechanism_negative"] = three_channel["known"] & ~three_channel["positive"]
        for name, mask in electro_fusion["channel_positive"].items():
            result_arrays[f"electrostatic_positive__{name}"] = mask
            result_arrays[f"electrostatic_negative__{name}"] = electro_fusion["channel_negative"][name]
            result_arrays[f"electrostatic_unknown__{name}"] = electro_fusion["channel_unknown"][name]
            for reason, unknown_mask in electro_fusion["channel_unknown_reasons"][name].items():
                result_arrays[f"electrostatic_unknown_reason__{name}__{reason}"] = unknown_mask
        for reason, unknown_mask in electro_fusion["unknown_reason_masks"].items():
            result_arrays[f"electrostatic_composite_unknown_reason__{reason}"] = unknown_mask
        for name, mask in three_channel["channel_positive"].items():
            result_arrays[f"three_mechanism_positive__{name}"] = mask
            result_arrays[f"three_mechanism_negative__{name}"] = three_channel["channel_negative"][name]
            result_arrays[f"three_mechanism_unknown__{name}"] = three_channel["channel_unknown"][name]
            for reason, unknown_mask in three_channel["channel_unknown_reasons"][name].items():
                result_arrays[f"three_mechanism_unknown_reason__{name}__{reason}"] = unknown_mask
        for reason, unknown_mask in three_channel["unknown_reason_masks"].items():
            result_arrays[f"three_mechanism_composite_unknown_reason__{reason}"] = unknown_mask
        npz_path = out / "field_fusion_evidence.npz"; _check_disk(out)
        np.savez_compressed(npz_path, **result_arrays)
        result = {"research_profile": PROFILE, "candidate_id": cid, "label": LABELS[cid],
            "source_snapshot_sha256": source_snapshot["tree_sha256"], "fusion_protocol_sha256": fusion_protocol,
            "evidence_slot_identities_sha256": slot_identity_sha256,
            "null_bank_candidates_sha256": manifest.get("_independent_null_bank", {}).get("candidates_sha256"),
            "status": "APBS_COMPLETE_WITH_UNKNOWN_PRESERVED", "reference_id": reference_id(reference),
            "samples": len(area), "fixed_denominator_area_a2": float(area.sum()),
            "base_grid": grids[0], "sensitivity_grids": grids[1],
            "pdb2pqr_binary_sha256": pqr_sha, "apbs_binary_sha256": apbs_sha,
            "preparation_protocol": PREPARATION_PROTOCOL,
            "production_preparation": production_audits,
            "null_content_diversity": {role: {"logical_WT_slots": 10,
                "unique_WT_content": len({item[3]["wt_sha256"] for item in null_preparations[role]}),
                "unique_PQR_content": len({item[3]["pqr_sha256"] for item in null_preparations[role]})}
                for role in ("null_train", "null_holdout")},
            "null_sampling_interpretation": "FoldX execution slots are distinct; repeated coordinate or PQR content is reported as dependence, not independent structure evidence",
            "null_pair_audits": null_audits,
            "production_solvent_mesh_audits": production_solvent_audits,
            "null_solvent_mesh_audits": null_solvent_audits,
            "pair_convergence": [{"slot": index, "converged_area_a2": float(area[grid_valid[index]].sum()),
                                  "unknown_area_a2": float(area[~grid_valid[index]].sum()),
                                  "preparation_valid": bool(pair_preparation_valid[index])}
                                 for index in range(5)],
            "scope_results": summaries,
            "electrostatics_interpretation": "potential differences in kT/e, with matched medium and empirical WT null; no binding probability",
            "geometry_chemistry_not_cleared_by_electrostatic_unknown": True,
            "structure_contact_scale": "NOT_FIT_YET; 6/9/12/15 A are geodesic sensitivity scenarios only",
            "raw_evidence_sha256": sha256_file(npz_path), "raw_evidence_path": str(npz_path),
            "wall_seconds": time.perf_counter() - started,
            "legacy_official_12A_coverage": 0.8602623530512513,
            "legacy_panel_unchanged": True}
        write_json(out / "field_fusion_report.json", result)
        report[cid] = result
        ledger["attempts"][-1].update(status="COMPLETE", wall_seconds=result["wall_seconds"])
        write_json(ledger_path, ledger)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--study", type=Path, default=DEFAULT_STUDY)
    parser.add_argument("--phase", choices=("freeze", "null_sampling", "geometry", "apbs"), required=True)
    parser.add_argument("--candidate-id", action="append", dest="candidate_ids")
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()
    try:
        if args.phase == "freeze":
            result = freeze(args.source, args.study)
            result = {key: result[key] for key in ("research_profile", "study_id", "source_report_sha256",
                "candidate_ids", "pilot_candidate_ids", "WT_content_diversity", "source_snapshot_sha256", "status")}
        elif args.phase == "geometry":
            ids = args.candidate_ids or list(PILOT_IDS)
            result = run_geometry(args.study, ids, chunk_size=args.chunk_size)
        elif args.phase == "null_sampling":
            result = run_null_sampling(args.study, args.candidate_ids or list(PILOT_IDS))
        else:
            result = run_apbs(args.study, args.candidate_ids or list(PILOT_IDS))
        print(json.dumps(result, ensure_ascii=False, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else str(x)))
    except Exception as exc:
        ledger_path = args.study / "study_ledger.json"
        if ledger_path.exists():
            ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
            for attempt in reversed(ledger.get("attempts", [])):
                if attempt.get("status") == "RUNNING":
                    attempt.update(status="ENGINEERING_FAILED", reason=f"{type(exc).__name__}: {exc}")
                    break
            write_json(ledger_path, ledger)
        raise


if __name__ == "__main__":
    main()
