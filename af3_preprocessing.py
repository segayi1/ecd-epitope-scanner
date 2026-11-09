from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Any

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

from .models import ProteinContext
from .structure_analysis import (
    ChainMapping,
    ParsedStructure,
    best_chain_mapping,
    load_structure,
    standard_residues,
)
from .util import PipelineError

_MAX_MODELS = 100
_MAX_TOTAL_BYTES = 1024**3
_MAX_CONFIDENCE_BYTES = 64 * 1024**2
_MAX_TOKENS = 2048
_MAX_PAE_CELLS = _MAX_TOKENS**2
_EPS = 1e-6
_MODEL_RE = re.compile(r"(?:(?P<prefix>.*)_)?model_(?P<index>\d+)\.cif$", re.I)
_FULL_RE = re.compile(r"(?:(?P<prefix>.*)_)?full_data_(?P<index>\d+)\.json$", re.I)
_SUMMARY_RE = re.compile(r"(?:(?P<prefix>.*)_)?summary_confidences_(?P<index>\d+)\.json$", re.I)
_CONFIDENCES_RE = re.compile(r"(?:(?P<prefix>.*)_)?confidences_(?P<index>\d+)\.json$", re.I)
_SEED_SAMPLE_MODEL_RE = re.compile(r"(?:(?P<prefix>.*)_)?seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model\.cif$", re.I)
_SEED_SAMPLE_FULL_RE = re.compile(r"(?:(?P<prefix>.*)_)?seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_(?:full_data|confidences)\.json$", re.I)
_SEED_SAMPLE_SUMMARY_RE = re.compile(r"(?:(?P<prefix>.*)_)?seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_summary_confidences\.json$", re.I)
_SEED_RE = re.compile(r"seed[-_](\d+)", re.I)
_SAMPLE_RE = re.compile(r"sample[-_](\d+)", re.I)


@dataclass(slots=True)
class AF3Selection:
    parsed: ParsedStructure
    mapping: ChainMapping
    confidence: dict[str, Any]
    report: dict[str, Any]
    structure_path: Path
    confidence_path: Path
    report_path: Path
    files: list[str]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_files(root: Path) -> list[Path]:
    root = root.resolve(strict=True)
    files = []
    for path in root.rglob("*"):
        if path.is_symlink():
            raise PipelineError(f"AF3 input contains a symlink: {path.relative_to(root)}")
        if path.is_file():
            resolved = path.resolve(strict=True)
            if root not in resolved.parents:
                raise PipelineError("AF3 input path escapes the selected folder")
            files.append(path)
    return files


def _sidecars(files: list[Path]) -> tuple[dict[tuple[Path, str, int], Path], dict[tuple[Path, str, int], Path], dict[tuple[Path, str, int], Path]]:
    structures, confidence, summaries = {}, {}, {}
    for path in files:
        flat_model = _SEED_SAMPLE_MODEL_RE.fullmatch(path.name)
        flat_full = _SEED_SAMPLE_FULL_RE.fullmatch(path.name)
        flat_summary = _SEED_SAMPLE_SUMMARY_RE.fullmatch(path.name)
        if flat_model or flat_full or flat_summary:
            match = flat_model or flat_full or flat_summary
            key = (path.parent, match.group("prefix") or "", int(match.group("sample")))
            if flat_model:
                structures[key] = path
            elif flat_summary:
                summaries[key] = path
            else:
                if path.stat().st_size > _MAX_CONFIDENCE_BYTES:
                    raise PipelineError("AF3 confidence file exceeds 64 MiB")
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    target = confidence
                else:
                    target = confidence if isinstance(payload, dict) and "atom_plddts" in payload else summaries
                target[key] = path
            continue
        role, rx = ("structure", _MODEL_RE) if _MODEL_RE.fullmatch(path.name) else ("confidence", _FULL_RE) if _FULL_RE.fullmatch(path.name) else ("summary", _SUMMARY_RE) if _SUMMARY_RE.fullmatch(path.name) else ("confidences", _CONFIDENCES_RE) if _CONFIDENCES_RE.fullmatch(path.name) else (None, None)
        if not role or not rx:
            if path.name.lower() == "model.cif":
                key = (path.parent, "", 0)
                structures[key] = path
                continue
            if path.name.lower() not in {"full_data.json", "confidences.json", "summary_confidences.json"}:
                continue
            key = (path.parent, "", 0)
            if path.name.lower() == "full_data.json":
                target = confidence
            elif path.name.lower() == "summary_confidences.json":
                target = summaries
            else:
                if path.stat().st_size > _MAX_CONFIDENCE_BYTES:
                    raise PipelineError("AF3 confidence file exceeds 64 MiB")
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    target = confidence
                else:
                    target = confidence if isinstance(payload, dict) and "atom_plddts" in payload else summaries
            if key in target:
                raise PipelineError(f"Duplicate AF3 sidecar for {path.name}")
            target[key] = path
            continue
        match = rx.fullmatch(path.name)
        key = (path.parent, match.group("prefix") or "", int(match.group("index")))
        if role == "confidences":
            if path.stat().st_size > _MAX_CONFIDENCE_BYTES:
                raise PipelineError("AF3 confidence file exceeds 64 MiB")
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                target = confidence
            else:
                target = confidence if isinstance(payload, dict) and "atom_plddts" in payload else summaries
        else:
            target = {"structure": structures, "confidence": confidence, "summary": summaries}[role]
        if key in target:
            raise PipelineError(f"Duplicate AF3 {role} file for {path.name}")
        target[key] = path
    return structures, confidence, summaries


def _sample_identity(path: Path, root: Path) -> tuple[str, str, str]:
    encoded = _SEED_SAMPLE_MODEL_RE.fullmatch(path.name)
    if encoded:
        return encoded.group("prefix") or root.name, encoded.group("seed"), encoded.group("sample")
    parts = list(path.parent.relative_to(root).parts)
    rel = "/".join(parts)
    seed = _SEED_RE.search(rel)
    sample = _SAMPLE_RE.search(rel)
    parent_parts = [part for part in parts if not _SEED_RE.search(part) and not _SAMPLE_RE.search(part)]
    job = "/".join(parent_parts) or root.name
    return job, seed.group(1) if seed else "unknown", sample.group(1) if sample else path.parent.name


def _request_protein_sequences(metadata: dict[str, Any]) -> list[str]:
    found = []
    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if key in {"protein", "proteinChain"} and isinstance(child, dict) and isinstance(child.get("sequence"), str):
                    found.append(child["sequence"].upper())
                else:
                    visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)
    visit(metadata.get("sequences", metadata.get("sequences_and_bonds", [])))
    return sorted(found)


def _observed_system(path: Path) -> dict[str, Any]:
    parsed = load_structure(path)
    proteins = sorted("".join(seq1(residue.resname, custom_map={"UNK": "X"}) for residue in standard_residues(chain))
                      for chain in parsed.model if standard_residues(chain))
    ligands = sorted(f"{chain.id}:{residue.resname}" for chain in parsed.model for residue in chain
                     if residue.id[0] not in {" ", "W"} and residue.resname.upper() not in {"HOH", "WAT"})
    return {"protein_sequences": proteins, "ligands": ligands}


def _fingerprint(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _context_signature(metadata: dict[str, Any]) -> dict[str, Any]:
    aliases = {"model_params": ("model_parameters", "modelParams")}
    signature = {}
    for name in ("dialect", "version", "msa", "templates", "model_params", "input_scope", "protocol"):
        signature[name] = metadata.get(name)
        if signature[name] is None:
            signature[name] = next((metadata[key] for key in aliases.get(name, ()) if key in metadata), None)
    return signature


def _metadata(path: Path, files: list[Path], root: Path) -> dict[str, Any]:
    def metadata_name(name: str) -> bool:
        low = name.lower()
        return (low in {"job_request.json", "summary.json", "fold_input.json", "input.json", "job_data.json"}
                or low.endswith("_job_request.json")
                or (low.endswith("_data.json") and low != "full_data.json"))
    merged: dict[str, Any] = {}
    ancestors = [path.parent, *path.parent.parents]
    ancestors = [parent for parent in ancestors if parent == root or root in parent.parents]
    found = sorted((p for p in files if metadata_name(p.name) and p.parent in ancestors), key=lambda p: (len(p.parts), p.name))
    for item in found:
        try:
            data = json.loads(item.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise PipelineError(f"Invalid AF3 metadata file {item.name}: {exc}") from exc
        if isinstance(data, list):
            if len(data) != 1 or not isinstance(data[0], dict):
                raise PipelineError(f"AF3 metadata {item.name} contains multiple job records; submit them separately")
            data = data[0]
        if isinstance(data, dict):
            merged.update(data)
    return merged


def _metadata_sources(path: Path, files: list[Path], root: Path) -> list[Path]:
    def candidate(name: str) -> bool:
        low = name.lower()
        return low in {"job_request.json", "summary.json", "fold_input.json", "input.json", "job_data.json"} or low.endswith("_job_request.json") or (low.endswith("_data.json") and low != "full_data.json")
    ancestors = [parent for parent in (path.parent, *path.parent.parents) if parent == root or root in parent.parents]
    return sorted((p for p in files if candidate(p.name) and p.parent in ancestors), key=lambda p: (len(p.parts), p.name))


def _compatibility(rows: list[dict[str, Any]]) -> None:
    by_job: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_job[row["job_id"]].append(row)
    fields = ("sequences", "entities", "ptms", "ligands", "bondedAtomPairs", "dialect", "version", "msa", "templates", "model_params", "input_scope", "protocol")
    signatures = {}
    unknown = []
    for job, items in by_job.items():
        values = []
        metadata_values = [item["metadata"] for item in items]
        # These chemistry collections are defined as empty when omitted by an
        # input schema that provides sequences/entities; hidden protocol fields
        # remain unknown and cannot establish cross-job compatibility.
        schema_context = all(any(key in m for key in ("sequences", "entities", "sequences_and_bonds")) for m in metadata_values)
        for field in fields:
            def normalized(meta):
                aliases = {"sequences": ("sequences_and_bonds",), "model_params": ("model_parameters", "modelParams"), "bondedAtomPairs": ("bonded_atom_pairs",)}
                if field in meta:
                    return meta[field]
                for alias in aliases.get(field, ()):
                    if alias in meta:
                        return meta[alias]
                return [] if schema_context and field in {"ptms", "ligands", "bondedAtomPairs"} else None
            found = [normalized(item["metadata"]) for item in items]
            if found and any(value != found[0] for value in found[1:]):
                raise PipelineError(f"Conflicting AF3 {field} metadata within job {job}")
            values.append((field, found[0] if found else None))
        signatures[job] = tuple(values)
        if any(value is None for _, value in values):
            unknown.append(job)
    if len(by_job) < 2:
        return
    if unknown:
        raise PipelineError("AF3 jobs have unknown cross-job input context; submit them separately")
    first = next(iter(signatures.values()))
    for job, signature in signatures.items():
        for (field, left), (_, right) in zip(first, signature, strict=True):
            if left != right:
                raise PipelineError(f"AF3 jobs conflict on {field}; cannot combine {job}")


def import_af3_directory(directory: Path, output_dir: Path) -> Path:
    """Freeze AF3 model/confidence pairs and provenance into output_dir/af3_inputs."""
    root = Path(directory).resolve(strict=True)
    if not root.is_dir():
        raise PipelineError("AF3 input must be a directory")
    files = _safe_files(root)
    if sum(path.stat().st_size for path in files) > _MAX_TOTAL_BYTES:
        raise PipelineError("AF3 input exceeds 1 GiB")
    # Enforce the per-file JSON limit before sidecar discovery or metadata
    # parsing can read any JSON payload, including optional summaries.
    if any(path.suffix.lower() == ".json" and path.stat().st_size > _MAX_CONFIDENCE_BYTES for path in files):
        raise PipelineError("AF3 JSON file exceeds 64 MiB")
    model_count = sum(bool(_MODEL_RE.fullmatch(path.name) or _SEED_SAMPLE_MODEL_RE.fullmatch(path.name) or path.name.lower() == "model.cif") for path in files)
    if model_count > _MAX_MODELS:
        raise PipelineError(f"AF3 input exceeds {_MAX_MODELS} models")
    structures, confidences, summaries = _sidecars(files)
    if not structures:
        raise PipelineError("No AF3 *_model_i.cif files were found")
    keys = sorted(structures, key=lambda key: (str(key[0]), key[1], key[2]))
    if len(keys) > _MAX_MODELS:
        raise PipelineError(f"AF3 input exceeds {_MAX_MODELS} models")
    total = sum(p.stat().st_size for p in files)
    if total > _MAX_TOTAL_BYTES:
        raise PipelineError("AF3 input exceeds 1 GiB")
    for key in keys:
        if key in confidences and confidences[key].stat().st_size > _MAX_CONFIDENCE_BYTES:
            raise PipelineError("AF3 confidence file exceeds 64 MiB")
    dest = Path(output_dir).resolve() / "af3_inputs"
    if dest.exists():
        raise PipelineError("AF3 frozen destination already exists; choose a new output directory")
    dest.mkdir(parents=True, exist_ok=True)
    rows = []
    excluded_members = []
    aliases: dict[tuple[str, str, str], tuple[str, str, str]] = {}
    for key in keys:
        cif = structures[key]
        if key not in confidences:
            excluded_members.append({"source_name": cif.name, "reason": "missing confidence sidecar", "sha256": _sha(cif)})
            continue
        conf = confidences[key]
        summary = summaries.get(key)
        job_id, seed_id, sample_id = _sample_identity(cif, root)
        hashes = (_sha(cif), _sha(conf), _sha(summary) if summary else "")
        # Root-level convenience copies are aliases only when the byte-identical
        # model/confidence/summary triplet also exists in a seed/sample folder.
        if cif.parent == root:
            peer = next((candidate for candidate in keys if structures[candidate].parent != root and (_sha(structures[candidate]), _sha(confidences[candidate]), _sha(summaries[candidate]) if candidate in summaries else "") == hashes), None)
            if peer:
                aliases[(str(cif), *hashes)] = (str(structures[peer]), *hashes)
                continue
        rel_dir = Path("af3_inputs") / "models" / f"m{len(rows):03d}"
        out_cif = dest.parent / rel_dir / "model.cif"
        out_conf = dest.parent / rel_dir / "full_data.json"
        out_summary = dest.parent / rel_dir / "summary_confidences.json" if summary else None
        out_cif.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cif, out_cif)
        shutil.copyfile(conf, out_conf)
        if summary and out_summary:
            shutil.copyfile(summary, out_summary)
        if _sha(out_cif) != hashes[0] or _sha(out_conf) != hashes[1] or (out_summary and _sha(out_summary) != hashes[2]):
            raise PipelineError("AF3 source files changed while freezing the input")
        metadata = _metadata(cif, files, root)
        requested = _request_protein_sequences(metadata)
        try:
            observed = _observed_system(cif)
        except Exception:
            observed = None  # Structure parsing remains a per-model selection diagnostic.
        if requested and observed is not None and requested != observed["protein_sequences"]:
            for member in (out_cif, out_conf, out_summary):
                if member and member.exists():
                    member.unlink()
            excluded_members.append({"source_name": cif.name, "reason": "observed CIF protein chains disagree with AF3 request sequences", "sha256": hashes[0], "requested_system": requested, "observed_system": observed})
            continue
        seed_values = metadata.get("modelSeeds", metadata.get("model_seeds"))
        if seed_id == "unknown" and isinstance(seed_values, list) and len(seed_values) == 1:
            seed_id = str(seed_values[0])
        metadata_refs = []
        for metadata_source in _metadata_sources(cif, files, root):
            digest = _sha(metadata_source)
            metadata_rel = Path("af3_inputs") / "metadata" / f"{digest}.json"
            metadata_target = dest.parent / metadata_rel
            metadata_target.parent.mkdir(parents=True, exist_ok=True)
            if not metadata_target.exists():
                shutil.copyfile(metadata_source, metadata_target)
            if _sha(metadata_target) != digest:
                raise PipelineError("AF3 metadata changed while freezing the input")
            metadata_refs.append({"path": metadata_rel.as_posix(), "sha256": digest, "source_name": metadata_source.name})
        rows.append({
            "id": f"m{len(rows):03d}", "job_id": job_id, "seed_id": seed_id,
            "sample_id": sample_id, "model_index": key[2], "metadata": metadata,
            "structure": (rel_dir / "model.cif").as_posix(),
            "confidence": (rel_dir / "full_data.json").as_posix(),
            "summary": (rel_dir / "summary_confidences.json").as_posix() if out_summary else None,
            "sha256": {"structure": hashes[0], "confidence": hashes[1], "summary": hashes[2] or None},
            "metadata_files": metadata_refs,
            "source_name": cif.name,
            "source_relpath": cif.relative_to(root).as_posix(),
            "observed_system": observed or None,
            "requested_system": requested or None,
            "system_fingerprint": _fingerprint({"request": requested or None, "observed": observed}) if observed is not None else None,
            "context_fingerprint": _fingerprint(_context_signature(metadata)) if all(value is not None for value in _context_signature(metadata).values()) else None,
            "context_status": "KNOWN" if all(value is not None for value in _context_signature(metadata).values()) else "UNKNOWN",
        })
    systems_by_job: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["observed_system"] is not None:
            systems_by_job[row["job_id"]].add(_fingerprint(row["observed_system"]))
    if any(len(systems) > 1 for systems in systems_by_job.values()):
        raise PipelineError("AF3 models within a job have mixed observed protein/ligand systems")
    _compatibility(rows)
    if not rows:
        raise PipelineError("No AF3 models have a confidence sidecar")
    manifest = {"schema_version": 1, "source_directory_name": root.name, "models": rows, "excluded_members": excluded_members, "aliases": [{"source": k[0], "sha256": list(k[1:]), "canonical": v[0]} for k, v in aliases.items()]}
    path = dest / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _manifest(path: Path) -> tuple[dict[str, Any], Path]:
    path = Path(path).resolve(strict=True)
    if path.stat().st_size > _MAX_CONFIDENCE_BYTES:
        raise PipelineError("AF3 manifest exceeds 64 MiB")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PipelineError(f"Invalid AF3 manifest: {exc}") from exc
    if data.get("schema_version") != 1 or not isinstance(data.get("models"), list) or not data["models"]:
        raise PipelineError("Unsupported or empty AF3 manifest")
    if len(data["models"]) > _MAX_MODELS:
        raise PipelineError(f"AF3 manifest exceeds {_MAX_MODELS} models")
    base = path.parent.parent
    total_bytes = 0
    for row in data["models"]:
        frozen_metadata = {}
        for role in ("structure", "confidence", "summary"):
            rel = row.get(role)
            if not rel:
                continue
            lexical = base / rel
            if Path(rel).is_absolute() or ".." in Path(rel).parts or any(part.is_symlink() for part in (lexical, *lexical.parents) if part != base.parent):
                raise PipelineError("AF3 manifest contains a symlink or unsafe member path")
            candidate = lexical.resolve(strict=True)
            if base not in candidate.parents:
                raise PipelineError("AF3 manifest path escapes its frozen input directory")
            size = candidate.stat().st_size
            total_bytes += size
            if role == "confidence" and size > _MAX_CONFIDENCE_BYTES:
                raise PipelineError("AF3 confidence file exceeds 64 MiB")
            if total_bytes > _MAX_TOTAL_BYTES:
                raise PipelineError("AF3 manifest exceeds 1 GiB")
            expected = row.get("sha256", {}).get(role)
            if not expected:
                raise PipelineError(f"AF3 manifest is missing {role} provenance for {row.get('id')}")
            if expected and _sha(candidate) != expected:
                raise PipelineError(f"Frozen AF3 {role} file changed: {row.get('id')}")
        for item in row.get("metadata_files", []):
            rel = Path(item["path"])
            candidate = (base / rel).resolve(strict=True)
            if rel.is_absolute() or ".." in rel.parts or base not in candidate.parents or candidate.is_symlink():
                raise PipelineError(f"Frozen AF3 metadata changed or has unsafe path: {row.get('id')}")
            total_bytes += candidate.stat().st_size
            if total_bytes > _MAX_TOTAL_BYTES or candidate.stat().st_size > _MAX_CONFIDENCE_BYTES:
                raise PipelineError("AF3 metadata exceeds configured input limits")
            if _sha(candidate) != item.get("sha256"):
                raise PipelineError(f"Frozen AF3 metadata changed or has unsafe path: {row.get('id')}")
            try:
                metadata_item = json.loads(candidate.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise PipelineError(f"Invalid frozen AF3 metadata for {row.get('id')}: {exc}") from exc
            if isinstance(metadata_item, list):
                if len(metadata_item) != 1 or not isinstance(metadata_item[0], dict):
                    raise PipelineError(f"Frozen AF3 metadata contains multiple job records: {row.get('id')}")
                metadata_item = metadata_item[0]
            if not isinstance(metadata_item, dict):
                raise PipelineError(f"Frozen AF3 metadata is not an object: {row.get('id')}")
            frozen_metadata.update(metadata_item)
        if row.get("metadata_files") and frozen_metadata != row.get("metadata"):
            raise PipelineError(f"AF3 manifest metadata disagrees with frozen source metadata: {row.get('id')}")
        requested_system = _request_protein_sequences(frozen_metadata) or None
        if requested_system != row.get("requested_system"):
            raise PipelineError(f"AF3 manifest request system disagrees with frozen metadata: {row.get('id')}")
        if requested_system and row.get("observed_system") and requested_system != row["observed_system"].get("protein_sequences"):
            raise PipelineError(f"AF3 request and observed systems disagree: {row.get('id')}")
        expected_system_fp = _fingerprint({"request": requested_system, "observed": row.get("observed_system")}) if row.get("observed_system") else None
        context = _context_signature(frozen_metadata)
        context_known = all(value is not None for value in context.values())
        expected_context_fp = _fingerprint(context) if context_known else None
        if row.get("system_fingerprint") != expected_system_fp or row.get("context_fingerprint") != expected_context_fp or row.get("context_status") != ("KNOWN" if context_known else "UNKNOWN"):
            raise PipelineError(f"AF3 manifest fingerprints disagree with frozen source context: {row.get('id')}")
        if row.get("source_relpath"):
            source_root = Path("/af3_source") / data.get("source_directory_name", "source")
            source_path = source_root / row["source_relpath"]
            if Path(row["source_relpath"]).is_absolute() or ".." in Path(row["source_relpath"]).parts:
                raise PipelineError("AF3 manifest contains unsafe source identity")
            job_id, seed_id, sample_id = _sample_identity(source_path, source_root)
            seed_values = row.get("metadata", {}).get("modelSeeds", row.get("metadata", {}).get("model_seeds"))
            if seed_id == "unknown" and isinstance(seed_values, list) and len(seed_values) == 1:
                seed_id = str(seed_values[0])
            if (row.get("job_id"), row.get("seed_id"), row.get("sample_id")) != (job_id, seed_id, sample_id):
                raise PipelineError(f"AF3 manifest grouping disagrees with frozen source identity: {row.get('id')}")
    return data, base


def first_af3_structure(manifest_path: Path) -> Path:
    data, base = _manifest(manifest_path)
    rows = sorted(data["models"], key=lambda item: (item["job_id"], item["seed_id"], item["sample_id"], item["model_index"], item["sha256"]["structure"]))
    for row in rows:
        try:
            parsed = load_structure(base / row["structure"])
            confidence = json.loads((base / row["confidence"]).read_text(encoding="utf-8"))
            _raw_atom_map(parsed, confidence, row["sha256"]["confidence"])
            return parsed.path
        except Exception:
            continue
    raise PipelineError("No parseable AF3 structure with valid atom confidence was found")


def preserve_af3_manifest(manifest_path: Path, output_dir: Path) -> Path:
    """Copy only hash-verified frozen members while retaining manifest paths and identities."""
    data, source_base = _manifest(manifest_path)
    destination_base = Path(output_dir).resolve()
    expected_manifest = destination_base / "af3_inputs" / "manifest.json"
    if Path(manifest_path).resolve() == expected_manifest:
        return expected_manifest
    for row in data["models"]:
        members = [(role, row.get(role), row.get("sha256", {}).get(role)) for role in ("structure", "confidence", "summary")]
        members += [("metadata", item.get("path"), item.get("sha256")) for item in row.get("metadata_files", [])]
        for role, rel, digest in members:
            if not rel:
                continue
            source = source_base / rel
            target = destination_base / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.resolve() == target.resolve():
                continue
            shutil.copyfile(source, target)
            if _sha(target) != digest:
                raise PipelineError(f"Could not preserve frozen AF3 {role} file: {row.get('id')}")
    out_manifest = expected_manifest
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _manifest(out_manifest)
    return out_manifest


def _clean(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text in {".", "?"} else text


def _raw_atom_map(parsed: ParsedStructure, confidence: dict[str, Any], confidence_sha: str) -> tuple[dict[tuple[str, int, str], float], list[float]]:
    try:
        d = MMCIF2Dict(str(parsed.path))
    except Exception as exc:
        raise PipelineError(f"Could not read AF3 atom identities: {exc}") from exc
    cols = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id", "label_seq_id", "auth_atom_id", "auth_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "pdbx_PDB_model_num"]
    vectors = {name: d.get("_atom_site." + name) for name in cols}
    fallbacks = {"auth_atom_id": "label_atom_id", "auth_comp_id": "label_comp_id",
                 "auth_asym_id": "label_asym_id", "auth_seq_id": "label_seq_id"}
    for name, fallback in fallbacks.items():
        if vectors[name] is None:
            vectors[name] = d.get("_atom_site." + fallback)
    n = len(vectors["id"] or [])
    if not n or any(not isinstance(v, list) or len(v) != n for v in vectors.values()):
        raise PipelineError("AF3 CIF has incomplete _atom_site identity columns")
    vals = confidence.get("atom_plddts", confidence.get("atom_plddt"))
    if not isinstance(vals, list) or len(vals) != n:
        raise PipelineError("AF3 atom_plddts length does not match raw _atom_site rows")
    supplied_chains = confidence.get("atom_chain_ids")
    if supplied_chains is not None and (not isinstance(supplied_chains, list) or len(supplied_chains) != n):
        raise PipelineError("AF3 atom_chain_ids length does not match raw _atom_site rows")
    raw_by_id = {}
    raw_order_keys = []
    for i in range(n):
        try:
            chain, seq = _clean(vectors["auth_asym_id"][i]), int(float(vectors["auth_seq_id"][i]))
            resname = _clean(vectors["auth_comp_id"][i]).upper()
            atomname = _clean(vectors["auth_atom_id"][i]).upper()
            element = _clean(vectors["type_symbol"][i]).upper()
            model_num = int(float(vectors["pdbx_PDB_model_num"][i]))
            plddt = float(vals[i])
        except (TypeError, ValueError, IndexError) as exc:
            raise PipelineError("Invalid AF3 atom identity or pLDDT value") from exc
        if not math.isfinite(plddt) or not 0 <= plddt <= 100:
            raise PipelineError("AF3 atom pLDDT must be finite and between 0 and 100")
        if supplied_chains is not None and str(supplied_chains[i]).strip() not in {
            _clean(vectors["label_asym_id"][i]), chain,
        }:
            raise PipelineError("AF3 atom_chain_ids do not match raw CIF chain identities")
        key = (chain, seq, _clean(vectors["pdbx_PDB_ins_code"][i]), resname, atomname, _clean(vectors["label_alt_id"][i]), element, model_num)
        if key in raw_by_id:
            raise PipelineError("AF3 CIF contains duplicate atom identities")
        raw_by_id[key] = plddt
        raw_order_keys.append(key)
    parsed_atoms = list(parsed.model.get_atoms())
    parsed_by_id = {}
    parsed_keys = []
    for atom in parsed_atoms:
        residue, chain = atom.get_parent(), atom.get_parent().get_parent()
        key = (str(chain.id), int(residue.id[1]), _clean(residue.id[2]), residue.resname.upper(), atom.id.upper(), _clean(atom.altloc), atom.element.upper(), 1)
        if key in parsed_by_id:
            raise PipelineError("Parsed AF3 structure contains duplicate atom identities")
        parsed_by_id[key] = atom
        parsed_keys.append(key)
    if set(raw_by_id) != set(parsed_by_id):
        raise PipelineError("Raw CIF atom identities cannot be mapped exactly onto parsed atoms")
    mapped = [raw_by_id[key] for key in parsed_keys]
    for atom, plddt in zip(parsed_atoms, mapped, strict=True):
        if abs(float(atom.bfactor) - plddt) > 0.1 + _EPS:
            raise PipelineError("AF3 CIF B factor disagrees with mapped pLDDT by more than 0.1")
    grouped: dict[tuple[str, int, str], list[float]] = defaultdict(list)
    for key, value in raw_by_id.items():
        grouped[key[:3]].append(value)
    return {key: mean(items) for key, items in grouped.items()}, mapped


def _pae(confidence: dict[str, Any], parsed: ParsedStructure, token_limit: int) -> tuple[np.ndarray | None, dict[tuple[str, int], int], str]:
    chains, residues = confidence.get("token_chain_ids"), confidence.get("token_res_ids")
    if (isinstance(chains, list) and len(chains) > token_limit) or (isinstance(residues, list) and len(residues) > token_limit):
        raise PipelineError("AF3 token mapping exceeds 2048 tokens")
    raw = confidence.get("pae", confidence.get("predicted_aligned_error"))
    if raw is None or raw == []:
        return None, {}, "missing"
    try:
        matrix = np.asarray(raw, dtype=float)
    except (TypeError, ValueError):
        return None, {}, "invalid"
    if matrix.ndim == 2 and matrix.shape[0] > token_limit:
        raise PipelineError("AF3 PAE exceeds the 2048-token limit")
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.size > _MAX_PAE_CELLS or not np.isfinite(matrix).all() or (matrix < 0).any():
        return None, {}, "invalid"
    if not isinstance(chains, list) or len(chains) != matrix.shape[0]:
        return matrix, {}, "unmapped"
    if not isinstance(residues, list) or len(residues) != matrix.shape[0]:
        label_positions: dict[str, set[int]] = defaultdict(set)
        for model_chain in parsed.model:
            for residue in standard_residues(model_chain):
                auth_key = _mapping_residue_key(residue, str(model_chain.id))
                label_key = parsed.label_lookup.get(auth_key)
                if label_key and label_key[1] is not None:
                    label_positions[str(label_key[0])].add(int(label_key[1]))
        token_indices: dict[str, list[int]] = defaultdict(list)
        for index, chain in enumerate(chains):
            token_indices[str(chain)].append(index)
        if set(token_indices) != set(label_positions) or any(len(token_indices[chain]) != len(label_positions[chain]) for chain in token_indices):
            return matrix, {}, "unmapped"
        token_map = {}
        for chain, indices in token_indices.items():
            for label_seq, index in zip(sorted(label_positions[chain]), indices, strict=True):
                token_map[(chain, label_seq)] = index
        return matrix, token_map, "mapped"
    token_map = {}
    for i, (chain, resid) in enumerate(zip(chains, residues, strict=True)):
        try:
            token_map[(str(chain), int(resid))] = i
        except (TypeError, ValueError):
            return matrix, {}, "unmapped"
    if len(token_map) != len(chains):
        return matrix, {}, "unmapped"
    return matrix, token_map, "mapped"


def _atom_identity_to_token(parsed: ParsedStructure, token_map: dict[tuple[str, int], int]) -> dict[tuple[str, int, str], int]:
    try:
        d = MMCIF2Dict(str(parsed.path))
    except Exception:
        return {}
    chains, seqs, auth_chains, auth_seqs, ins = (d.get("_atom_site." + key) for key in ("label_asym_id", "label_seq_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"))
    auth_chains = auth_chains or chains
    auth_seqs = auth_seqs or seqs
    if not all(isinstance(x, list) for x in (chains, seqs, auth_chains, auth_seqs, ins)):
        return {}
    result = {}
    for label_chain, label_seq, chain, seq, icode in zip(chains, seqs, auth_chains, auth_seqs, ins, strict=True):
        try:
            token = token_map[(str(label_chain).strip(), int(float(label_seq)))]
            result[(str(chain).strip(), int(float(seq)), _clean(icode))] = token
        except (ValueError, KeyError, TypeError):
            continue
    return result


def _mapping_residue_key(residue: Any, chain_id: str) -> tuple[str, int, str]:
    return chain_id, int(residue.id[1]), _clean(residue.id[2])


def _metrics(parsed: ParsedStructure, mapping: ChainMapping, plddt: dict[tuple[str, int, str], float], protein: ProteinContext) -> tuple[dict[int, float], dict[int, np.ndarray], dict[int, int]]:
    scores, cas, tokens = {}, {}, {}
    for pos, residue in mapping.uniprot_to_residue.items():
        if not protein.ecd_start <= pos <= protein.ecd_end:
            continue
        key = _mapping_residue_key(residue, mapping.chain_id)
        if key in plddt:
            scores[pos] = plddt[key]
        if residue.has_id("CA"):
            cas[pos] = np.asarray(residue["CA"].coord, dtype=float)
        tokens[pos] = int(parsed.label_lookup.get((key[0], key[1], key[2]), ("", None))[1] or 0)
    return scores, cas, tokens


def _near_pairs(models: list[dict[str, Any]], positions: set[int]) -> set[tuple[int, int]]:
    pairs = set()
    for model in models:
        present = sorted(positions & set(model["cas"]))
        if len(present) < 2:
            continue
        coords = np.asarray([model["cas"][position] for position in present])
        pairs.update((present[i], present[j]) for i, j in cKDTree(coords).query_pairs(15 + _EPS) if present[j] - present[i] > 2)
    return pairs


def _is_clash(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value > 0
    return str(value).strip().lower() in {"true", "yes", "1"}


def _quality_order(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    metrics = ("ptm", "mean_plddt", "coverage", "summary_rank")
    comparable = [name for name in metrics if all(isinstance(row.get(name), (int, float)) and not isinstance(row.get(name), bool) and math.isfinite(float(row[name])) for row in rows)]
    return sorted(rows, key=lambda row: tuple(-float(row[name]) for name in comparable) + (row["row"]["sha256"]["structure"],))


def _summary_ptm(summary: dict[str, Any], parsed: ParsedStructure, chain_id: str) -> float | None:
    values = summary.get("chain_ptm")
    if isinstance(values, (int, float)) and not isinstance(values, bool) and math.isfinite(float(values)) and 0 <= values <= 1:
        return float(values)
    if isinstance(values, list):
        chain_ids = summary.get("chain_ids")
        unique_ids = list(dict.fromkeys(str(value) for value in chain_ids)) if isinstance(chain_ids, list) else []
        if chain_id in unique_ids and len(unique_ids) == len(values):
            value = values[unique_ids.index(chain_id)]
            return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and 0 <= value <= 1 else None
        protein_chains = [chain.id for chain in parsed.model if standard_residues(chain)]
        if len(values) == 1 and len(protein_chains) == 1 and protein_chains[0] == chain_id:
            value = values[0]
            return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and 0 <= value <= 1 else None
    value = summary.get("ptm")
    protein_chains = [chain.id for chain in parsed.model if standard_residues(chain)]
    if len(protein_chains) == 1 and protein_chains[0] == chain_id and isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and 0 <= value <= 1:
        return float(value)
    return None


def _quality(row: dict[str, Any], protein: ProteinContext, chain: str | None, admission: Callable[[ParsedStructure, ChainMapping], Any] | None) -> dict[str, Any]:
    structure, conf_path = row["_structure"], row["_confidence"]
    try:
        parsed = load_structure(structure)
    except PipelineError:
        raise
    except Exception as exc:
        raise PipelineError(f"Could not parse AF3 CIF {row['id']}: {exc}") from exc
    mapping = best_chain_mapping(parsed, protein, preferred_chains=[chain] if chain else None)
    summary_data = row.get("summary_data", {})
    score = _summary_ptm(summary_data, parsed, mapping.chain_id)
    summary_rank = summary_data.get("ranking_score")
    if not isinstance(summary_rank, (int, float)) or isinstance(summary_rank, bool) or not math.isfinite(float(summary_rank)) or not -100 <= summary_rank <= 1.5:
        summary_rank = None
    try:
        confidence = json.loads(conf_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"admitted": False, "selection_valid": False,
                "reason": f"Invalid AF3 full_data JSON for {row['id']}: {exc}",
                "row": row, "parsed": parsed, "mapping": mapping,
                "ptm": score, "summary_rank": summary_rank, "coverage": mapping.observed_ecd_coverage, "summary_data": summary_data}
    try:
        plddt, _ = _raw_atom_map(parsed, confidence, row["sha256"]["confidence"])
    except PipelineError as exc:
        return {"admitted": False, "selection_valid": False, "reason": str(exc), "row": row,
                "parsed": parsed, "mapping": mapping, "ptm": score, "summary_rank": summary_rank, "coverage": mapping.observed_ecd_coverage,
                "summary_data": summary_data}
    plddt = {key: value for key, value in plddt.items() if key[0] == mapping.chain_id}
    if admission:
        try:
            result = admission(parsed, mapping)
        except PipelineError as exc:
            return {"admitted": False, "selection_valid": False, "reason": str(exc), "row": row,
                    "parsed": parsed, "mapping": mapping, "ptm": score, "summary_rank": summary_rank, "coverage": mapping.observed_ecd_coverage, "summary_data": summary_data}
        if result is not None:
            return {"admitted": False, "selection_valid": False, "reason": str(result), "row": row,
                    "parsed": parsed, "mapping": mapping, "ptm": score, "summary_rank": summary_rank, "coverage": mapping.observed_ecd_coverage, "summary_data": summary_data}
    clashes = row.get("summary_data", {}).get("has_clash")
    if clashes is None:
        summary = row.get("_summary")
        if summary and summary.exists():
            try:
                clashes = json.loads(summary.read_text(encoding="utf-8")).get("has_clash")
            except (OSError, json.JSONDecodeError):
                pass
    scores, cas, _ = _metrics(parsed, mapping, plddt, protein)
    target = range(protein.ecd_start, protein.ecd_end + 1)
    f70 = sum(scores.get(pos, -1) >= 70 for pos in target) / protein.ecd_length
    target_scores = [scores[pos] for pos in range(protein.ecd_start, protein.ecd_end + 1) if pos in scores]
    mean_plddt = float(np.mean(target_scores)) if target_scores else 0.0
    pae, token_map, pae_status = _pae(confidence, parsed, _MAX_TOKENS)
    atom_token = _atom_identity_to_token(parsed, token_map) if token_map else {}
    pae_res = {}
    if atom_token:
        for pos, residue in mapping.uniprot_to_residue.items():
            key = _mapping_residue_key(residue, mapping.chain_id)
            label_key = parsed.label_lookup.get(key)
            if label_key and label_key[1] is not None:
                token = token_map.get((label_key[0], label_key[1]))
                if token is not None:
                    pae_res[pos] = token
    # summary ranking values are used only within a comparable token/template scope.
    summary_data = row.get("summary_data", {})
    return {"row": row, "parsed": parsed, "mapping": mapping, "selection_valid": True,
            "plddt": plddt, "scores": scores, "cas": cas, "f70": f70,
            "mean_plddt": mean_plddt, "has_clash": clashes,
            "coverage": mapping.observed_ecd_coverage,
            "pae_status": pae_status if pae_status in {"missing", "invalid"} or pae_res else ("unmapped" if pae is not None else "missing"),
            "pae_res": pae_res, "summary_data": summary_data, "ptm": score,
            "summary_rank": float(summary_rank) if summary_rank is not None else None}


def select_af3_model(manifest_path: Path, protein: ProteinContext, *, chain: str | None = None, allow_partial: bool = False, admission: Callable[[ParsedStructure, ChainMapping], Any] | None = None, output_dir: Path) -> AF3Selection:
    data, base = _manifest(manifest_path)
    candidates = []
    for row in data["models"]:
        candidate = dict(row)
        candidate["_structure"] = base / row["structure"]
        candidate["_confidence"] = base / row["confidence"]
        candidate["_summary"] = base / row["summary"] if row.get("summary") else None
        if candidate["_summary"]:
            try:
                candidate["summary_data"] = json.loads(candidate["_summary"].read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                candidate["summary_data"] = {}
        else:
            candidate["summary_data"] = {}
        try:
            candidates.append(_quality(candidate, protein, chain, admission))
        except PipelineError as exc:
            candidates.append({"admitted": False, "reason": str(exc), "row": candidate})
    admitted = [x for x in candidates if x.get("selection_valid")]
    if not admitted:
        reasons = "; ".join(f"{x['row']['id']}: {x.get('reason', 'not admitted')}" for x in candidates)
        raise PipelineError("No AF3 model passed structure, confidence, mapping, and admission checks" + (f" ({reasons})" if reasons else ""))
    admitted = [x for x in admitted if not _is_clash(x.get("has_clash"))]
    if not admitted:
        raise PipelineError("All AF3 models with valid evidence have has_clash=true")
    best_f70 = max(x["f70"] for x in admitted)
    q = [x for x in admitted if x["f70"] + _EPS >= best_f70 - 0.05]
    best_mean = max(x["mean_plddt"] for x in q)
    q = [x for x in q if x["mean_plddt"] + _EPS >= best_mean - 5.0]
    groups = sorted({(x["row"]["job_id"], x["row"]["seed_id"]) for x in q})
    samples_by_group: dict[tuple[str, str], set[str]] = defaultdict(set)
    for x in q:
        group = (x["row"]["job_id"], x["row"]["seed_id"])
        samples_by_group[group].add(x["row"]["sample_id"])
    weights = {}
    for x in q:
        group = (x["row"]["job_id"], x["row"]["seed_id"])
        weights[x["row"]["id"]] = 1 / len(groups) / len(samples_by_group[group]) / sum(1 for y in q if (y["row"]["job_id"], y["row"]["seed_id"], y["row"]["sample_id"]) == (*group, x["row"]["sample_id"]))
    rcore = set.intersection(*(set(x["cas"]) & {p for p, v in x["scores"].items() if v >= 70} for x in q)) if q else set()
    valid_models = [x for x in candidates if x.get("cas")]
    diag_counts: dict[int, int] = defaultdict(int)
    for x in valid_models:
        for position in x["cas"]:
            diag_counts[position] += 1
    rdiag = {position for position, count in diag_counts.items() if count >= 2}
    pair_ids = sorted(_near_pairs(q, rcore))
    edge_acc = {pair: [0.0, 0.0, 0.0, 0.0] for pair in pair_ids}
    # Process one full-data JSON and PAE matrix at a time; only edge summaries survive.
    for x in q:
        model_id = x["row"]["id"]
        model_weight = weights[model_id]
        for pair in pair_ids:
            a, b = pair
            if a not in x["cas"] or b not in x["cas"]:
                continue
            distance = float(np.linalg.norm(x["cas"][a] - x["cas"][b]))
            acc = edge_acc[pair]
            acc[0] += model_weight
            acc[1] += model_weight * distance
            if distance <= 15 + _EPS:
                acc[2] += model_weight
        if x["pae_status"] == "mapped":
            confidence_data = json.loads(x["row"]["_confidence"].read_text(encoding="utf-8"))
            pae_matrix, _, _ = _pae(confidence_data, x["parsed"], _MAX_TOKENS)
            for a, b in pair_ids:
                ti, tj = x["pae_res"].get(a), x["pae_res"].get(b)
                if pae_matrix is not None and ti is not None and tj is not None:
                    legal = max(float(pae_matrix[ti, tj]), float(pae_matrix[tj, ti])) <= 15 + _EPS
                    edge_acc[(a, b)][3] += model_weight if legal else 0.0
                    # The last field is reserved for the available PAE vote weight below.
                    edge_acc[(a, b)].append(model_weight)
            del confidence_data, pae_matrix
    edges = []
    for (pos_i, pos_j), (distance_weight, weighted_distance, distance_support, pae_support, *pae_weights) in edge_acc.items():
        if not distance_weight:
            continue
        record = {"residue_i": pos_i, "residue_j": pos_j,
                  "weighted_distance": weighted_distance / distance_weight,
                  "distance_support": distance_support, "pae_support": pae_support,
                  "pae_available_weight": sum(pae_weights)}
        if distance_support + _EPS >= 0.5 and pae_support + _EPS >= 0.5:
            edges.append(record)
    diagnostic_acc = {pair: [0, 0, []] for pair in sorted(_near_pairs(valid_models, rdiag))}
    for model in valid_models:
        for (a, b), acc in diagnostic_acc.items():
            if a not in model["cas"] or b not in model["cas"]:
                continue
            distance = float(np.linalg.norm(model["cas"][a] - model["cas"][b]))
            acc[0] += 1
            acc[1] += distance <= 15 + _EPS
            acc[2].append(model["row"]["id"])
    diag_edges = [{"residue_i": a, "residue_j": b,
                   "distance_support": acc[1] / acc[0], "observed_models": acc[0],
                   "near_models": acc[1]}
                  for (a, b), acc in diagnostic_acc.items() if acc[0] and acc[1]]
    ensemble_ids = [x["row"]["id"] for x in q]
    edge_pairs = [(edge["residue_i"], edge["residue_j"]) for edge in edges]
    pair_distances = {
        x["row"]["id"]: {
            (a, b): float(np.linalg.norm(x["cas"][a] - x["cas"][b]))
            for a, b in edge_pairs if a in x["cas"] and b in x["cas"]
        }
        for x in q
    }
    d_matrix = {left: {} for left in ensemble_ids}
    for left in q:
        for right in q:
            shared = set(pair_distances[left["row"]["id"]]) & set(pair_distances[right["row"]["id"]])
            d_matrix[left["row"]["id"]][right["row"]["id"]] = (
                float(np.mean([abs(pair_distances[left["row"]["id"]][edge] - pair_distances[right["row"]["id"]][edge]) for edge in shared]))
                if shared else 0.0
            )
    # Fixed Ecore and fixed sample weights; every model keeps its full vote.
    energies = {left: sum(d_matrix[left][right] * weights[right] for right in ensemble_ids) for left in ensemble_ids}
    core_ok = len(q) > 1 and len(rcore) >= 10 and len(edges) >= 20
    ranking = _quality_order(q)
    raw_medoid = min(q, key=lambda x: (energies[x["row"]["id"]], x["row"]["sha256"]["structure"]))
    near = [x for x in q if energies[x["row"]["id"]] <= energies[raw_medoid["row"]["id"]] + 0.05 + _EPS]
    selected = _quality_order(near)[0] if core_ok else ranking[0]
    if not core_ok and not allow_partial and (len(q) <= 1 or len(rcore) < 10 or len(edges) < 20):
        status = "INSUFFICIENT_ENSEMBLE_FALLBACK"
    else:
        status = "ENSEMBLE_CORE" if core_ok else "PARTIAL_ENSEMBLE_FALLBACK"
    # Global/core pair RMSD are direct CA RMSD after Kabsch alignment on shared core.
    def rmsd(left, right, positions):
        common = [p for p in positions if p in left["cas"] and p in right["cas"]]
        if len(common) < 3:
            return None
        a, b = np.array([left["cas"][p] for p in common]), np.array([right["cas"][p] for p in common])
        ac, bc = a.mean(0), b.mean(0)
        u, _, vt = np.linalg.svd((a-ac).T @ (b-bc))
        diag = np.eye(3)
        diag[-1, -1] = np.linalg.det(u @ vt)
        rot = u @ diag @ vt
        delta = (a-ac) @ rot - (b-bc)
        return float(np.sqrt(np.mean(np.sum(delta*delta, axis=1))))
    pair_rmsd = []
    for i, left in enumerate(q):
        for right in q[i+1:]:
            pair_rmsd.append({"left": left["row"]["id"], "right": right["row"]["id"], "core_rmsd": rmsd(left, right, rcore), "global_rmsd": rmsd(left, right, sorted(rdiag))})
    domains = [feature for feature in protein.features if feature.type.lower() in {"domain", "repeat"}]
    domain_geometry = {"domains": [], "inter_domain_distances": []}
    for feature in domains:
        positions = [p for p in range(max(protein.ecd_start, feature.start), min(protein.ecd_end, feature.end) + 1)]
        domain_geometry["domains"].append({
            "name": feature.description or f"{feature.type} {feature.start}-{feature.end}",
            "start": feature.start, "end": feature.end,
            "pairwise_rmsd": [{"left": left["row"]["id"], "right": right["row"]["id"],
                               "rmsd": rmsd(left, right, positions)}
                              for i, left in enumerate(valid_models) for right in valid_models[i+1:]],
        })
    for i, left_domain in enumerate(domains):
        left_positions = [p for p in range(max(protein.ecd_start, left_domain.start), min(protein.ecd_end, left_domain.end) + 1)]
        for right_domain in domains[i+1:]:
            right_positions = [p for p in range(max(protein.ecd_start, right_domain.start), min(protein.ecd_end, right_domain.end) + 1)]
            means = []
            for model in valid_models:
                distances = [float(np.linalg.norm(model["cas"][a] - model["cas"][b]))
                             for a in left_positions for b in right_positions if a in model["cas"] and b in model["cas"]]
                if distances:
                    means.append({"model_id": model["row"]["id"], "mean_distance": float(np.mean(distances))})
            domain_geometry["inter_domain_distances"].append({
                "left": left_domain.description or f"{left_domain.type} {left_domain.start}-{left_domain.end}",
                "right": right_domain.description or f"{right_domain.type} {right_domain.start}-{right_domain.end}",
                "by_model": means,
            })
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    selected_id = selected["row"]["id"]
    selected_dir = out / "af3_selected"
    selected_dir.mkdir(exist_ok=True)
    cif_out, conf_out = selected_dir / "selected.cif", selected_dir / "full_data.json"
    shutil.copyfile(selected["row"]["_structure"], cif_out)
    confidence = json.loads(selected["row"]["_confidence"].read_text(encoding="utf-8"))
    _, mapped_atoms = _raw_atom_map(selected["parsed"], confidence, selected["row"]["sha256"]["confidence"])
    confidence["atom_plddts"] = [float(x) for x in mapped_atoms]
    confidence["_validated_atom_plddts"] = {"|".join((key[0], str(key[1]), key[2])): value for key, value in selected["plddt"].items()}
    confidence["_source_sha256"] = selected["row"]["sha256"]["confidence"]
    confidence["_source_structure_sha256"] = selected["row"]["sha256"]["structure"]
    shutil.copyfile(selected["row"]["_confidence"], conf_out)
    ranking_evidence = [x for x in candidates if x.get("parsed") and x.get("summary_rank") is not None]
    original_rank1_score = max((x["summary_rank"] for x in ranking_evidence), default=None)
    original_rank1_ties = [x for x in ranking_evidence if original_rank1_score is not None and abs(x["summary_rank"] - original_rank1_score) <= _EPS]
    original_rank1 = min(original_rank1_ties, key=lambda x: (x["row"]["sha256"]["structure"], x["row"]["id"])) if original_rank1_ties else None
    ptm_evidence = [x for x in candidates if x.get("parsed") and x.get("ptm") is not None]
    ptm_best = min(ptm_evidence, key=lambda x: (-x["ptm"], x["row"]["sha256"]["structure"], x["row"]["id"]), default=None)
    plddt_best = max(admitted, key=lambda x: x["mean_plddt"])
    models_report = []
    q_ids = {x["row"]["id"] for x in q}
    original_rank1_exclusions = []
    for x in original_rank1_ties:
        reason = x.get("reason") or ("has_clash=true" if _is_clash(x.get("has_clash")) else "outside_Q" if x["row"]["id"] not in q_ids else None)
        if reason:
            original_rank1_exclusions.append({"model_id": x["row"]["id"], "reason": reason})
    for x in candidates:
        is_eligible = x in q
        models_report.append({
            "model_id": x["row"]["id"], "eligible": is_eligible,
            "exclusion_reasons": ([x.get("reason")] if x.get("reason") else
                                  ["has_clash=true"] if _is_clash(x.get("has_clash")) else
                                  ["outside_Q"] if x.get("parsed") and not is_eligible else []),
            "coverage": x.get("mapping").observed_ecd_coverage if x.get("mapping") else None,
            "f70": x.get("f70"), "mean_plddt": x.get("mean_plddt"),
            "ptm": x.get("ptm"), "ranking_score": x.get("summary_rank"),
            "consensus_score": energies.get(x["row"]["id"]),
        })
    per_residue_spread = {}
    for pos in sorted(rdiag):
        neighbors = [other for other in rdiag if other != pos]
        neighbor_spreads = []
        for other in neighbors:
            distances = [float(np.linalg.norm(x["cas"][pos] - x["cas"][other])) for x in valid_models if pos in x["cas"] and other in x["cas"]]
            if len(distances) >= 2:
                neighbor_spreads.append(float(np.std(distances)))
        covering = [x for x in valid_models if pos in x["cas"]]
        low_conf = [x for x in covering if x["scores"].get(pos, 101) < 70]
        missing_conf = [x for x in covering if pos not in x["scores"]]
        per_residue_spread[str(pos)] = {
            "n_valid": len(covering), "n_low_confidence": len(low_conf), "n_missing_confidence": len(missing_conf),
            "mean_internal_distance_spread": float(np.mean(neighbor_spreads)) if neighbor_spreads else None,
            "mapped_plddt": [x["scores"].get(pos) for x in valid_models if pos in x["cas"]],
            "low_confidence_models": [x["row"]["id"] for x in low_conf],
            "missing_confidence_models": [x["row"]["id"] for x in missing_conf],
        }
    domain_pae = {"global_mean": None, "domains": [], "inter_domain": []}
    if selected["pae_status"] == "mapped":
        selected_full_data = json.loads(selected["row"]["_confidence"].read_text(encoding="utf-8"))
        selected_pae, _, _ = _pae(selected_full_data, selected["parsed"], _MAX_TOKENS)
        target_tokens = [selected["pae_res"][pos] for pos in range(protein.ecd_start, protein.ecd_end + 1) if pos in selected["pae_res"]]
        if selected_pae is not None and target_tokens:
            values = selected_pae[np.ix_(target_tokens, target_tokens)]
            domain_pae["global_mean"] = float(np.mean(values))
            seen_domains = set()
            for feature in protein.features:
                if feature.type.lower() not in {"domain", "repeat"}:
                    continue
                positions = [pos for pos in range(max(protein.ecd_start, feature.start), min(protein.ecd_end, feature.end) + 1) if pos in selected["pae_res"]]
                token_ids = [selected["pae_res"][pos] for pos in positions]
                if len(token_ids) < 2:
                    continue
                label = feature.description or f"{feature.type} {feature.start}-{feature.end}"
                identity = (label, feature.start, feature.end)
                if identity in seen_domains:
                    continue
                seen_domains.add(identity)
                domain_values = selected_pae[np.ix_(token_ids, token_ids)]
                domain_pae["domains"].append({"name": label, "start": feature.start, "end": feature.end, "mean_pae": float(np.mean(domain_values)), "token_count": len(token_ids)})
            for index, left in enumerate(domain_pae["domains"]):
                left_tokens = [selected["pae_res"][pos] for pos in range(max(protein.ecd_start, left["start"]), min(protein.ecd_end, left["end"]) + 1) if pos in selected["pae_res"]]
                for right in domain_pae["domains"][index + 1:]:
                    right_tokens = [selected["pae_res"][pos] for pos in range(max(protein.ecd_start, right["start"]), min(protein.ecd_end, right["end"]) + 1) if pos in selected["pae_res"]]
                    forward = selected_pae[np.ix_(left_tokens, right_tokens)]
                    reverse = selected_pae[np.ix_(right_tokens, left_tokens)].T
                    domain_pae["inter_domain"].append({"domain_i": left["name"], "domain_j": right["name"],
                        "mean_pae_i_to_j": float(np.mean(forward)), "mean_pae_j_to_i": float(np.mean(reverse)),
                        "mean_symmetric_max_pae": float(np.mean(np.maximum(forward, reverse)))})
        del selected_full_data, selected_pae
    report = {
        "algorithm_version": "af3_preprocessing_v1",
        "status": status, "method": "AF3_ENSEMBLE_CORE" if core_ok else "AF3_QUALITY_FALLBACK",
        "protocol": {
            "admission": "current pipeline callback", "exclude_has_clash": True,
            "Q": {"f70_within_best": 0.05, "mean_target_plddt_within_best": 5.0},
            "Rcore": {"common_CA": True, "minimum_plddt": 70.0},
            "Rdiag": {"minimum_corresponding_CA_models": 2},
            "Ecore": {"sequence_gap_gt": 2, "distance_A_max": 15.0,
                      "distance_support_min": 0.5, "symmetric_max_pae_A": 15.0,
                      "pae_support_min": 0.5},
            "ensemble_weights": "equal job+seed groups, then equal samples, then equal models per sample",
            "medoid": {"distance": "mean absolute CA internal-distance difference on fixed Ecore",
                       "self_pairs_included": True, "near_medoid_delta": 0.05},
            "fallback": "quality ranking when Q/Rcore/Ecore is insufficient",
        },
        "parameters": {"max_models": _MAX_MODELS, "max_total_bytes": _MAX_TOTAL_BYTES,
                       "max_confidence_bytes": _MAX_CONFIDENCE_BYTES, "max_tokens": _MAX_TOKENS,
                       "max_pae_cells": _MAX_PAE_CELLS, "numerical_epsilon": _EPS},
        "selected_model_id": selected_id,
        "system_fingerprint": _fingerprint({"source": selected["row"].get("system_fingerprint"),
                                             "target": {"accession": protein.accession, "ecd_start": protein.ecd_start,
                                                        "ecd_end": protein.ecd_end, "ecd_sequence": protein.ecd_sequence,
                                                        "chain": selected["mapping"].chain_id}}),
        "context_fingerprint": selected["row"].get("context_fingerprint"),
        "context_status": selected["row"].get("context_status", "UNKNOWN"),
        "baselines": {"original_rank1": original_rank1["row"]["id"] if original_rank1 else None,
                      "original_rank1_ties": sorted(x["row"]["id"] for x in original_rank1_ties),
                      "original_rank1_exclusion_reason": original_rank1.get("reason") if original_rank1 else None,
                      "original_rank1_exclusions": original_rank1_exclusions,
                      "ptm_best": ptm_best["row"]["id"] if ptm_best else None,
                      "original_rank1_value": original_rank1_score,
                      "ptm_best_value": ptm_best["ptm"] if ptm_best else None,
                      "raw_medoid": raw_medoid["row"]["id"] if core_ok else None,
                      "plddt_best": plddt_best["row"]["id"], "plddt_best_value": plddt_best["mean_plddt"]},
        "warnings": (["No stable ensemble core; selected by model quality"] if not core_ok else []),
        "selected_sha256": selected["row"]["sha256"], "Q": [x["row"]["id"] for x in q], "weights": weights,
        "maxpTM": ptm_best["ptm"] if ptm_best else None,
        "Rcore": sorted(rcore), "Rdiag": sorted(rdiag), "Ecore": edges, "diagnostic_edges": diag_edges,
        "pair_rmsd": pair_rmsd, "domain_geometry": domain_geometry,
        "core_max_rmsd": max((r["core_rmsd"] or 0 for r in pair_rmsd), default=0),
        "core_max_rmsd_flag": max((r["core_rmsd"] or 0 for r in pair_rmsd), default=0) > 2.0,
        "models": models_report,
        "model_diagnostics": [{"model_id": x["row"]["id"], "has_clash": x.get("has_clash"), "weight": weights.get(x["row"]["id"]), "pae_status": x.get("pae_status"), "sha256": x["row"]["sha256"], "system_fingerprint": x["row"].get("system_fingerprint"), "context_fingerprint": x["row"].get("context_fingerprint"), "context_status": x["row"].get("context_status", "UNKNOWN")} for x in candidates],
        "low_confidence_positions": sorted(p for p, v in selected["scores"].items() if v < 70),
        "missing_confidence_positions": sorted(set(range(protein.ecd_start, protein.ecd_end+1)) - set(selected["scores"])),
        "per_residue_diagnostics": per_residue_spread,
        "domain_pae": domain_pae, "minority_contacts": [edge for edge in diag_edges if edge["distance_support"] < 0.5],
    }
    report_path = out / "af3_selection_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # Return parser loaded from the exact delivered CIF so callers can safely reuse it.
    delivered = load_structure(cif_out)
    result_conf = confidence
    return AF3Selection(delivered, best_chain_mapping(delivered, protein, preferred_chains=[chain] if chain else None), result_conf, report, cif_out, conf_out, report_path, [p.relative_to(out).as_posix() for p in (cif_out, conf_out, report_path)])
