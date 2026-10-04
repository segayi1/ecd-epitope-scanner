"""GRASP-JAX constraint formatting, FASTA mapping, and deterministic ZIP output."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import zipfile
from dataclasses import dataclass
from typing import Any

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
SOURCES = ("experimental_ir", "predicted_ir", "manual_rpr")


class ConstraintError(ValueError):
    code = "FASTA_CONSTRAINT_MISMATCH"

    def __init__(self, message: str, *, code: str | None = None):
        super().__init__(message)
        if code is not None:
            self.code = code


@dataclass(frozen=True)
class FastaChain:
    chain_index: int
    name: str
    sequence: str


def fasta_bytes(chains: list[FastaChain]) -> bytes:
    if [item.chain_index for item in chains] != [1, 2, 3]:
        raise ConstraintError("FASTA must contain chain order 1=ECD, 2=VH, 3=VL")
    lines = []
    for chain in chains:
        sequence = chain.sequence.upper()
        if not sequence or set(sequence) - AMINO_ACIDS or any(c in chain.name for c in "\r\n>"):
            raise ConstraintError(f"Invalid amino-acid sequence for chain {chain.chain_index}")
        lines.append(f">{chain.chain_index}|{chain.name}")
        lines.extend(sequence[start : start + 60] for start in range(0, len(sequence), 60))
    return ("\n".join(lines) + "\n").encode("ascii")


def format_ir(chain_index: int, residue_index: int, residue_type: str) -> str:
    aa = residue_type.upper()
    if chain_index not in {1, 2, 3} or residue_index < 1 or aa not in AMINO_ACIDS:
        raise ConstraintError("IR requires a valid chain index, 1-based position, and amino acid")
    return f"{chain_index}-{residue_index}-{aa}"


def format_rpr(
    chain_a: int,
    residue_a: int,
    residue_type_a: str,
    chain_b: int,
    residue_b: int,
    residue_type_b: str,
    distance_angstrom: float,
    fdr: float = 0.05,
) -> str:
    if chain_a == chain_b and residue_a == residue_b:
        raise ConstraintError("RPR cannot pair a residue with itself")
    if not math.isfinite(distance_angstrom) or not 4 <= distance_angstrom <= 32:
        raise ConstraintError("RPR distance cutoff must be finite and within GRASP bins [4, 32] Å")
    if not math.isfinite(fdr) or not 0 < fdr <= 1:
        raise ConstraintError("RPR FDR must be in (0, 1]")
    left = format_ir(chain_a, residue_a, residue_type_a)
    right = format_ir(chain_b, residue_b, residue_type_b)
    return f"{left}, {right}, {distance_angstrom:.3f}, {fdr:.5g}"


def parse_constraint(line: str) -> dict[str, Any]:
    value = line.strip()
    if not value or value.startswith("#"):
        raise ConstraintError("Constraint row is empty or a comment")
    if "," not in value:
        fields = value.split("-")
        if len(fields) != 3:
            raise ConstraintError(f"Invalid IR syntax: {value}")
        try:
            chain_index, residue_index = int(fields[0]), int(fields[1])
        except ValueError as exc:
            raise ConstraintError(f"Invalid IR position: {value}") from exc
        normalized = format_ir(chain_index, residue_index, fields[2])
        return {"kind": "IR", "constraint": normalized}
    fields = [item.strip() for item in value.split(",")]
    if len(fields) not in {3, 4}:
        raise ConstraintError(f"RPR requires 3 or 4 comma-separated fields: {value}")
    first = parse_constraint(fields[0])
    second = parse_constraint(fields[1])
    if first["kind"] != "IR" or second["kind"] != "IR":
        raise ConstraintError(f"RPR endpoints must use chain-position-amino-acid syntax: {value}")
    left_chain, left_position, left_aa = fields[0].split("-")
    right_chain, right_position, right_aa = fields[1].split("-")
    try:
        distance = float(fields[2])
        fdr = float(fields[3]) if len(fields) == 4 else 0.05
    except ValueError as exc:
        raise ConstraintError(f"RPR distance and FDR must be numbers: {value}") from exc
    normalized = format_rpr(
        int(left_chain), int(left_position), left_aa,
        int(right_chain), int(right_position), right_aa,
        distance, fdr,
    )
    return {"kind": "RPR", "constraint": normalized}


def _duplicate_key(parsed: dict[str, Any]) -> str:
    if parsed["kind"] == "IR":
        return parsed["constraint"]
    left, right, *_ = [item.strip() for item in parsed["constraint"].split(",")]
    return "RPR:" + ":".join(sorted((left, right)))


def validate_constraint_against_fasta(
    parsed: dict[str, Any], chain_sequences: dict[int, str]
) -> None:
    constraint = parsed["constraint"]
    endpoints = constraint.split(",") if parsed["kind"] == "RPR" else [constraint]
    for endpoint in endpoints[:2]:
        fields = endpoint.strip().split("-")
        chain_index, residue_index, residue_type = int(fields[0]), int(fields[1]), fields[2]
        sequence = chain_sequences.get(chain_index)
        if sequence is None:
            raise ConstraintError(f"Unknown FASTA chain {chain_index}")
        if residue_index > len(sequence):
            raise ConstraintError(
                f"Chain {chain_index} residue {residue_index} exceeds sequence length {len(sequence)}"
            )
        expected = sequence[residue_index - 1]
        if expected != residue_type:
            raise ConstraintError(
                f"FASTA mismatch at chain {chain_index} residue {residue_index}: "
                f"constraint has {residue_type}, FASTA has {expected}"
            )


def select_experimental_ir(
    effect_document: dict[str, Any],
    probe_id: str,
    mapping_rows: list[dict[str, Any]],
    ecd_sequence: str,
) -> tuple[str, dict[str, Any]]:
    """Map a reviewed single-site effect probe through the frozen ECD mapping."""
    projection = effect_document.get("surface_projection") or {}
    probe = next((item for item in projection.get("probes", []) if item.get("probe_id") == probe_id), None)
    if (
        probe is None or probe.get("attribution") != "SINGLE_SITE"
        or len(probe.get("mutated_positions", [])) != 1
        or probe.get("decline_pct") is None
        or probe.get("valid_decline_batch_count", 0) < 1
    ):
        raise ConstraintError(
            "EXPERIMENTAL_IR_SINGLE_SITE_REQUIRED: source probe must have a valid single-site decline",
            code="EXPERIMENTAL_IR_SINGLE_SITE_REQUIRED",
        )
    uniprot_position = int(probe["mutated_positions"][0])
    mapping = next(
        (item for item in mapping_rows if item.get("uniprot_pos") and int(item["uniprot_pos"]) == uniprot_position),
        None,
    )
    if mapping is None or not mapping.get("normalized_position"):
        raise ConstraintError(
            "EXPERIMENTAL_IR_POSITION_UNMAPPED: source position has no frozen ECD mapping",
            code="EXPERIMENTAL_IR_POSITION_UNMAPPED",
        )
    ecd_position = int(mapping["normalized_position"])
    sequence = ecd_sequence.upper()
    if ecd_position < 1 or ecd_position > len(sequence):
        raise ConstraintError(
            "EXPERIMENTAL_IR_POSITION_UNMAPPED: source position is outside the frozen ECD sequence",
            code="EXPERIMENTAL_IR_POSITION_UNMAPPED",
        )
    wt = sequence[ecd_position - 1]
    if mapping.get("one_letter") and mapping["one_letter"] != wt:
        raise ConstraintError(
            "Frozen mapping residue differs from the ECD reference sequence",
            code="EXPERIMENTAL_IR_FASTA_MAPPING_MISMATCH",
        )
    return format_ir(1, ecd_position, wt), {
        "probe_id": probe_id,
        "uniprot_position": uniprot_position,
        "ecd_position": ecd_position,
        "attribution": "SINGLE_SITE",
        "decline_pct": probe["decline_pct"],
        "valid_decline_batch_count": probe["valid_decline_batch_count"],
    }
def _tsv_bytes(records: list[dict[str, Any]]) -> bytes:
    columns = ["constraint", "evidence_source", "evidence_id", "reviewer", "rationale", "details_json"]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for item in records:
        writer.writerow(
            {
                "constraint": item["constraint"],
                "evidence_source": item["evidence_source"],
                "evidence_id": item.get("evidence_id", ""),
                "reviewer": item["reviewer"],
                "rationale": item["rationale"],
                "details_json": json.dumps(item.get("details", {}), ensure_ascii=False, sort_keys=True),
            }
        )
    return stream.getvalue().encode("utf-8")


def build_constraint_package(
    *,
    chains: list[FastaChain],
    experimental_ir: list[dict[str, Any]],
    predicted_ir: list[dict[str, Any]],
    manual_rpr: list[dict[str, Any]],
    merge_sources: list[str],
    manifest_context: dict[str, Any],
) -> tuple[bytes, dict[str, Any]]:
    fasta = fasta_bytes(chains)
    chain_sequences = {item.chain_index: item.sequence.upper() for item in chains}
    chain_to_source = {
        "experimental_ir": experimental_ir,
        "predicted_ir": predicted_ir,
        "manual_rpr": manual_rpr,
    }
    if len(set(merge_sources)) != len(merge_sources) or set(merge_sources) - set(SOURCES):
        raise ConstraintError("merge_sources contains duplicates or an unknown source")
    rows: dict[str, list[dict[str, Any]]] = {}
    for source in SOURCES:
        records = []
        seen = set()
        for index, raw in enumerate(chain_to_source[source], 1):
            item = dict(raw)
            if not item.get("reviewer", "").strip() or not item.get("rationale", "").strip():
                raise ConstraintError(f"{source} row {index}: reviewer and rationale are required")
            try:
                parsed = parse_constraint(item["constraint"])
            except (KeyError, TypeError) as exc:
                raise ConstraintError(f"{source} row {index}: constraint is required") from exc
            expected_kind = "RPR" if source == "manual_rpr" else "IR"
            if parsed["kind"] != expected_kind:
                raise ConstraintError(f"{source} row {index}: expected {expected_kind}")
            validate_constraint_against_fasta(parsed, chain_sequences)
            key = _duplicate_key(parsed)
            if key in seen:
                raise ConstraintError(f"{source} row {index}: duplicate constraint")
            seen.add(key)
            records.append(
                {
                    "constraint": parsed["constraint"],
                    "evidence_source": source,
                    "evidence_id": str(item.get("evidence_id", "")),
                    "reviewer": item["reviewer"].strip(),
                    "rationale": item["rationale"].strip(),
                    "details": item.get("details", {}),
                }
            )
        rows[source] = records
    merged = [record for source in merge_sources for record in rows[source]]
    merged_constraints = [record["constraint"] for record in merged]
    merged_keys = [_duplicate_key(parse_constraint(value)) for value in merged_constraints]
    if len(merged_keys) != len(set(merged_keys)):
        raise ConstraintError("Merged selections contain duplicate constraints across sources")
    parsed_rows = []
    for line_number, record in enumerate(merged, 1):
        parsed = parse_constraint(record["constraint"])
        validate_constraint_against_fasta(parsed, chain_sequences)
        parsed_rows.append({"line": line_number, **parsed, "evidence_source": record["evidence_source"]})
    file_map: dict[str, bytes] = {
        "constraints.fasta": fasta,
        "experimental_ir.tsv": _tsv_bytes(rows["experimental_ir"]),
        "predicted_ir.tsv": _tsv_bytes(rows["predicted_ir"]),
        "manual_rpr.tsv": _tsv_bytes(rows["manual_rpr"]),
        "grasp_constraints.txt": ("\n".join(merged_constraints) + ("\n" if merged else "")).encode("utf-8"),
        "evidence_records.json": json.dumps(
            {source: rows[source] for source in SOURCES},
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
        ).encode("utf-8"),
        "validation_report.json": json.dumps(
            {"status": "VALID", "row_count": len(parsed_rows), "parsed_rows": parsed_rows},
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
        ).encode("utf-8"),
    }
    files = {
        name: {"sha256": hashlib.sha256(data).hexdigest(), "size_bytes": len(data)}
        for name, data in sorted(file_map.items())
    }
    export_manifest = {
        "schema_version": "1.0",
        "status": "VALID" if manual_rpr else "IR_ONLY",
        "record_order": ["1=ECD", "2=VH", "3=VL"],
        "merge_sources": merge_sources,
        "fasta_sha256": hashlib.sha256(fasta).hexdigest(),
        "constraint_count": len(parsed_rows),
        "rpr_count": len(manual_rpr),
        "limitations": list(manifest_context.get("limitations", [])),
        "provenance": manifest_context.get("provenance", {}),
        "files": files,
    }
    file_map["export_manifest.json"] = json.dumps(
        export_manifest, ensure_ascii=False, sort_keys=True, indent=2
    ).encode("utf-8")
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for name, data in sorted(file_map.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100444 << 16
            bundle.writestr(info, data)
    return archive.getvalue(), export_manifest
