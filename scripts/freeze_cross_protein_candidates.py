"""Freeze research-only, stratified double-mutant candidates for the 2026 cohort."""

from __future__ import annotations

import csv
import hashlib
import http.client
import itertools
import json
import math
import time
import urllib.error
import urllib.request
from dataclasses import asdict, replace
from heapq import heappop, heappush
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from scipy.spatial import cKDTree

from .. import config
from ..integrations._adapter import request_for
from ..joint_design.surface_units import build_surface_units
from ..joint_design_v4.contexts import residue_keys, write_normalized_prediction_chain
from ..joint_design_v4.coverage_profile import mechanisms as descriptor_mechanisms
from ..joint_design_v4.coverage_profile import sequence_reasons
from ..joint_design_v4.engine import _adapter, _map_external_rows
from ..joint_design_v4.evidence import build_single_evidence, compressed_options
from ..models import FeatureRange, ProteinContext
from ..retrieval import _location_value
from ..structure_analysis import analyze_residues, load_structure, map_chain_to_uniprot
from ..surface_perturbation.engine_v3 import coupling
from ..util import sha256_file

DEFAULT_ROOT = Path(r"E:\user\New project\cross_protein_validation_20261006")
RSASA = 0.20
PAIR_SPAN_A = 24.0
SUPPORT_RADIUS_A = 12.0
TARGET_PER_STRATUM = 2
CTLA4_HISTORICAL_SOURCE_SHA256 = {
    "ctla4_candidate_pool_v2.json": "a0609b05765b81887de7a5ff95118f878ed60814c10e2ce0ddbc3295281cc21d",
    "ctla4_combination_stability.csv": "e1a7247f182c1888150da30894155ef90b5ae6f39967adb67427efcf24c66f31",
}
STRATA = ("CONTACT_COUPLED_PROXY", "INTERMEDIATE_PROXY", "SPATIALLY_SEPARATED_PROXY")
SPECS = {
    "CTLA4": {
        "pdb": "1AH1",
        "target": "CTLA4",
        "accession": "P16410",
        "chain": "A",
        "range": [37, 161],
        "group": "development",
        "excluded": [147],
    },
    "ubiquitin": {
        "pdb": "1UBQ",
        "target": "ubiquitin_mature",
        "accession": "P0CG48",
        "chain": "A",
        "range": [1, 76],
        "reference_range": [1, 76],
        "group": "development",
        "excluded": [],
    },
    "PD1": {
        "pdb": "3RRQ",
        "target": "PD1",
        "accession": "Q15116",
        "chain": "A",
        "range": [25, 170],
        "group": "confirmation",
        "excluded": [132],
    },
    "PDL1": {
        "pdb": "4Z18",
        "target": "PDL1",
        "accession": "Q9NZQ7",
        "chain": "A",
        "range": [19, 239],
        "group": "confirmation",
        "excluded": [],
    },
    "B7H6": {
        "pdb": "3PV7",
        "target": "B7H6",
        "accession": "Q68D85",
        "chain": "A",
        "range": [25, 262],
        "group": "confirmation",
        "excluded": [],
    },
    "hen_lysozyme": {
        "pdb": "1LYZ",
        "target": "hen_lysozyme",
        "accession": "P00698",
        "chain": "A",
        "range": [19, 147],
        "group": "confirmation",
        "excluded": [],
    },
    "MBP": {
        "pdb": "1OMP",
        "target": "MBP",
        "accession": "P0AEX9",
        "chain": "A",
        "range": [27, 396],
        "group": "held_out",
        "excluded": [],
    },
    "T4_lysozyme": {
        "pdb": "2LZM",
        "target": "T4_lysozyme",
        "accession": "P00720",
        "chain": "A",
        "range": [1, 164],
        "group": "held_out",
        "excluded": [12, 137],
    },
}
CTLA4_FIXED = (
    ("A61W", "M134R"),
    ("D78A", "M122W"),
    ("Q117R", "E155W"),
    ("L47W", "Q117R"),
    ("T67R", "M134R"),
    ("S55F", "S101R"),
)


def _sha(path: Path) -> str:
    return sha256_file(path)


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".part")
    temp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True, default=str),
        encoding="utf-8",
    )
    temp.replace(path)


def _fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"invalid frozen FASTA: {path}")
    return "".join(line.strip() for line in lines[1:]).upper()


def _uniprot_entry(
    accession: str, source: str, expected_sequence: str, out: Path
) -> tuple[dict, str]:
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(out.glob(f"{accession}_*.json"))
    for path in files:
        entry = json.loads(path.read_text(encoding="utf-8"))
        if entry.get("sequence", {}).get("value", "").upper() == expected_sequence:
            return entry, _sha(path)
    request = urllib.request.Request(
        f"https://rest.uniprot.org/uniprotkb/{accession}.json",
        headers={
            "User-Agent": "ecd-epitope-scanner-cross-protein-freeze/1.0",
            "Accept": "application/json",
        },
    )
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                raw = response.read()
            entry = json.loads(raw)
            break
        except (http.client.HTTPException, urllib.error.URLError, TimeoutError, OSError):
            if attempt == 3:
                raise
            time.sleep(2**attempt)
    if entry.get("sequence", {}).get("value", "").upper() != expected_sequence:
        raise ValueError(f"current UniProt sequence does not match frozen FASTA for {accession}")
    digest = hashlib.sha256(raw).hexdigest()
    path = out / f"{accession}_{digest[:16]}.json"
    temp = path.with_suffix(".json.part")
    temp.write_bytes(raw)
    temp.replace(path)
    return entry, digest


def _features(entry: dict) -> list[FeatureRange]:
    result = []
    for item in entry.get("features", []):
        location = item.get("location", {})
        start, end = _location_value(location, "start"), _location_value(location, "end")
        if start is None or end is None:
            continue
        evidence = ",".join(x.get("evidenceCode", "") for x in item.get("evidences", []) if x)
        result.append(
            FeatureRange(
                item.get("type", "Feature"), start, end, item.get("description", "") or "", evidence
            )
        )
    return result


def _dijkstra(graph: dict[int, dict[int, float]], start: int) -> dict[int, float]:
    distance = {start: 0.0}
    queue = [(0.0, start)]
    while queue:
        cost, node = heappop(queue)
        if cost != distance[node]:
            continue
        for neighbor, edge in graph[node].items():
            proposed = cost + float(edge)
            if proposed < distance.get(neighbor, math.inf):
                distance[neighbor] = proposed
                heappush(queue, (proposed, neighbor))
    return distance


def _rows_by_author_id(rows: list[dict], chain_id: str) -> dict[tuple[int, str], dict]:
    output = {}
    for row in rows:
        if row["author_chain"] != chain_id:
            continue
        try:
            author_number = int(row["author_residue"])
        except (TypeError, ValueError):
            continue
        output[(author_number, row["insertion_code"] or "")] = row
    return output


def select_stratified(pool: list[dict], target_per_stratum: int = TARGET_PER_STRATUM) -> dict:
    """Select only complete 2-per-class strata, never fill an unavailable class."""
    output = {}
    used_ids = set()
    for stratum in STRATA:
        eligible = [row for row in pool if row["coupling_class"] == stratum]
        site_pairs = {tuple(row["site_pair"]) for row in eligible}
        if len(site_pairs) < target_per_stratum:
            output[stratum] = {
                "status": "EMPTY",
                "eligible_genotypes": len(eligible),
                "distinct_site_pairs": len(site_pairs),
                "selected": [],
            }
            continue
        selected = []
        selected_pairs = set()
        seen_bands = set()
        used_mechanisms = set()
        while len(selected) < target_per_stratum:
            remaining = [
                row
                for row in eligible
                if row["candidate_id"] not in used_ids
                and tuple(row["site_pair"]) not in selected_pairs
            ]
            if not remaining:
                break
            winner = min(
                remaining,
                key=lambda row: (
                    0 if tuple(row["exposure_profile"]) not in seen_bands else 1,
                    0 if set(row["mechanisms"]) - used_mechanisms else 1,
                    -float(row["proxy_support_area_a2"]),
                    float(row["positive_single_spurs_risk"]),
                    row["candidate_id"],
                ),
            )
            selected.append(winner)
            selected_pairs.add(tuple(winner["site_pair"]))
            seen_bands.add(tuple(winner["exposure_profile"]))
            used_mechanisms.update(winner["mechanisms"])
            used_ids.add(winner["candidate_id"])
        output[stratum] = {
            "status": "FROZEN" if len(selected) == target_per_stratum else "EMPTY",
            "eligible_genotypes": len(eligible),
            "distinct_site_pairs": len(site_pairs),
            "selected": selected if len(selected) == target_per_stratum else [],
        }
    return output


def _mapping_for_candidate(spec: dict, input_root: Path, maps: dict, protein: ProteinContext):
    pdb = spec["pdb"]
    source = input_root / f"{pdb}.cif"
    frozen = maps["structures"][pdb]
    if _sha(source).casefold() != frozen["input_cif_sha256"].casefold():
        raise ValueError(f"frozen mmCIF hash mismatch: {pdb}")
    entity = next(
        (
            item
            for item in frozen["entity_mappings"]
            if item["target"] == spec["target"] and spec["chain"] in item["chains"]
        ),
        None,
    )
    if entity is None:
        raise ValueError(f"canonical map has no {pdb} chain {spec['chain']} for {spec['target']}")
    parsed = load_structure(source)
    mapped = map_chain_to_uniprot(parsed, spec["chain"], protein)
    by_key = _rows_by_author_id(entity["residues"], spec["chain"])
    query_to_uniprot = {}
    uniprot_to_residue = {}
    mismatches = []
    for index, residue in enumerate(mapped.residues):
        key = (int(residue.id[1]), str(residue.id[2]).strip())
        row = by_key.get(key)
        if row is None:
            raise ValueError(f"canonical map omitted observed residue {pdb}:{key}")
        expected = row["reference_position"]
        actual = mapped.query_to_uniprot.get(index)
        if expected is None:
            continue
        expected = int(expected)
        if actual != expected:
            raise ValueError(
                f"mapping disagreement {pdb}:{key}: frozen={expected}, runtime={actual}"
            )
        if row["pdb_aa"] != row["reference_aa"] and expected not in set(spec["excluded"]):
            raise ValueError(f"unreviewed WT mismatch {pdb}:{key} at canonical {expected}")
        query_to_uniprot[index] = expected
        uniprot_to_residue[expected] = residue
        if row["pdb_aa"] != row["reference_aa"]:
            mismatches.append(expected)
    mapped = replace(
        mapped, query_to_uniprot=query_to_uniprot, uniprot_to_residue=uniprot_to_residue
    )
    ignored = (
        set(spec["excluded"])
        | set(entity.get("candidate_freeze_excluded_positions", []))
        | set(mismatches)
    )
    return parsed, mapped, entity, sorted(ignored)


def _pair_pool(
    records, graph, areas, menus, protein, mapped_residues, structure_sha: str, source_identity: str
) -> tuple[list[dict], dict]:
    positions = sorted(p for p in menus if p in graph)
    distances = {position: _dijkstra(graph, position) for position in positions}
    supports = {p: {q for q, d in distances[p].items() if d <= SUPPORT_RADIUS_A} for p in positions}
    record_by_pos = {r.uniprot_pos: r for r in records}
    atom_coords = {}
    for position in positions:
        residue = mapped_residues[position]
        atom_coords[position] = np.asarray(
            [
                atom.coord
                for atom in residue
                if str(getattr(atom, "element", "")).upper() not in {"H", "D"}
            ],
            dtype=float,
        )
    rows = []
    rejects = {}
    sequence = protein.sequence
    for pos_a, pos_b in itertools.combinations(positions, 2):
        span = distances[pos_a].get(pos_b)
        if span is None or span > PAIR_SPAN_A:
            rejects["outside_local_surface_graph_span"] = (
                rejects.get("outside_local_surface_graph_span", 0) + 1
            )
            continue
        support_a, support_b = supports[pos_a], supports[pos_b]
        union = support_a | support_b
        union_area = sum(float(areas.get(p, 0.0)) for p in union)
        overlap_area = sum(float(areas.get(p, 0.0)) for p in support_a & support_b)
        overlap = overlap_area / union_area if union_area else None
        tree = cKDTree(atom_coords[pos_a])
        min_atom = float(tree.query(atom_coords[pos_b], k=1)[0].min())
        ca_a = np.asarray(record_by_pos[pos_a].ca, float)
        ca_b = np.asarray(record_by_pos[pos_b].ca, float)
        ca_distance = float(np.linalg.norm(ca_a - ca_b))
        components = {
            "minimum_native_heavy_atom_distance_a": min_atom,
            "proxy_support_overlap_fraction": overlap,
        }
        coupling_class = coupling(components)
        for left, right in itertools.product(menus[pos_a], menus[pos_b]):
            mutations = sorted((left, right), key=lambda item: item["position"])
            substitutions = [(item["position"], item["wt"], item["mutant"]) for item in mutations]
            reasons = sequence_reasons(sequence, substitutions)
            if reasons:
                for reason in reasons:
                    rejects[reason] = rejects.get(reason, 0) + 1
                continue
            notations = [item["notation"] for item in mutations]
            identifier = (
                "CPV1_"
                + hashlib.sha256(f"{source_identity}|{','.join(notations)}".encode()).hexdigest()[
                    :16
                ]
            )
            mechanisms = sorted({axis for item in mutations for axis in item["mechanisms"]})
            exposure = tuple(sorted(item["exposure_band"] for item in mutations))
            rows.append(
                {
                    "candidate_id": identifier,
                    "status": "RESEARCH_CANDIDATE_NOT_COMBINATION_QUALIFIED",
                    "mutations": notations,
                    "site_pair": [pos_a, pos_b],
                    "coupling_class": coupling_class,
                    "coupling_components": components,
                    "ca_distance_a": ca_distance,
                    "surface_graph_span_proxy_a": float(span),
                    "support_overlap_fraction": overlap,
                    "proxy_support_area_a2": float(sum(float(areas.get(p, 0.0)) for p in union)),
                    "proxy_support_positions": sorted(union),
                    "exposure_profile": list(exposure),
                    "mechanisms": mechanisms,
                    "positive_single_spurs_risk": sum(
                        max(0.0, float(item["spurs_ddg"])) for item in mutations
                    ),
                    "single_spurs": [item["spurs_ddg"] for item in mutations],
                    "stability_interpretation": "single-mutation menu prior only; no pair SPURS or FoldX computed",
                    "support_basis": "CA-residue surface graph proxy; scheduling only",
                    "native_heavy_atom_distance_a": min_atom,
                    "structure_sha256": structure_sha,
                }
            )
    return rows, {
        "rejected_combinations_by_reason": rejects,
        "position_count_with_menu": len(positions),
        "local_site_pair_count": len({tuple(row["site_pair"]) for row in rows}),
        "legal_genotype_count": len(rows),
    }


def _freeze_target(
    name: str, spec: dict, root: Path, maps: dict, fasta_root: Path, output_root: Path
) -> dict:
    input_root = root / "inputs"
    fasta_path = fasta_root / f"{spec['accession']}.fasta"
    full_sequence = _fasta(fasta_path)
    reference_range = spec.get("reference_range")
    sequence = (
        full_sequence[reference_range[0] - 1 : reference_range[1]]
        if reference_range
        else full_sequence
    )
    if (
        hashlib.sha256(sequence.encode()).hexdigest().casefold()
        != maps["references"][spec["target"]]
        .get("reference_sequence_sha256", hashlib.sha256(sequence.encode()).hexdigest())
        .casefold()
    ):
        raise ValueError(f"canonical reference sequence hash mismatch for {spec['target']}")
    entry, uniprot_hash = _uniprot_entry(
        spec["accession"], spec["target"], full_sequence, root / "mapping" / "uniprot_entries"
    )
    start, end = spec["range"]
    protein = ProteinContext(
        spec["accession"],
        name,
        "",
        sequence,
        start,
        end,
        end,
        "frozen-study-target-range",
        reviewed="reviewed" in entry.get("entryType", "").lower(),
        features=_features(entry),
        source_url=f"https://rest.uniprot.org/uniprotkb/{spec['accession']}.json",
    )
    parsed, mapping, entity, ignored = _mapping_for_candidate(spec, input_root, maps, protein)
    source_path = input_root / f"{spec['pdb']}.cif"
    input_hashes = {
        "coordinate_sha256": _sha(source_path),
        "canonical_maps_sha256": _sha(root / "mapping" / "canonical_residue_maps.json"),
        "map_manifest_sha256": _sha(root / "mapping" / "canonical_mapping_manifest.json"),
        "input_manifest_sha256": _sha(input_root / "input_manifest.json"),
        "fasta_sha256": _sha(fasta_path),
        "uniprot_json_sha256": uniprot_hash,
    }
    settings = dict(config.JOINT_V4_ML_PARAMETERS["spurs_single"])
    code_paths = [
        Path(__file__).resolve(),
        Path(__file__).parents[1] / "config.py",
        Path(__file__).parents[1] / "structure_analysis.py",
        Path(__file__).parents[1] / "joint_design_v4" / "coverage_profile.py",
        Path(__file__).parents[1] / "joint_design_v4" / "evidence.py",
        Path(__file__).parents[1] / "joint_design_v4" / "engine.py",
        Path(__file__).parents[1] / "joint_design_v4" / "contexts.py",
        Path(__file__).parents[1] / "joint_design" / "surface_units.py",
        Path(__file__).parents[1] / "integrations" / "_adapter.py",
        Path(__file__).parents[1] / "integrations" / "spurs_adapter.py",
        Path(__file__).parents[1] / "integrations" / "workers" / "spurs_gpu_worker.py",
    ]
    code_hashes = {str(path.name): _sha(path) for path in code_paths}
    source_identity = hashlib.sha256(
        json.dumps(
            {
                "input": input_hashes,
                "code": code_hashes,
                "spurs_settings": settings,
                "rsasa": RSASA,
                "pair_span_a": PAIR_SPAN_A,
                "support_radius_a": SUPPORT_RADIUS_A,
                "target": name,
                "chain": spec["chain"],
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()
    target_dir = output_root / name / source_identity[:16]
    records, analysis = analyze_residues(
        parsed,
        mapping,
        protein,
        target_dir / "work",
        rsasa_threshold=RSASA,
        ignored_positions=set(ignored),
    )
    if not records:
        raise ValueError(f"no mapped, non-excluded structure residues for {name}")
    normalized = write_normalized_prediction_chain(
        source_path, target_dir / "prediction_context.pdb", records, spec["chain"]
    )
    if normalized is None:
        raise ValueError(f"normalized prediction structure was not written for {name}")
    prediction_sequence = "".join(record.aa for record in records)
    adapter = _adapter("spurs_single", settings, target_dir / "cache", tool_name="spurs_single")
    attempt_dir = target_dir / "spurs_attempts" / str(time.time_ns())
    request = request_for(
        adapter,
        normalized,
        spec["chain"],
        prediction_sequence,
        attempt_dir,
        structure_type="experimental_normalized_pdb",
        parameters=settings,
    )
    prediction = adapter.predict(request)
    if not prediction.ok:
        raise RuntimeError(f"SPURS single failed for {name}: {prediction.errors}")
    weights_hash = prediction.provenance.get("weights_sha256")
    model_files_hash = prediction.provenance.get("model_files_sha256")
    if not isinstance(weights_hash, str) or len(weights_hash) != 64 or not model_files_hash:
        raise RuntimeError(f"SPURS model identity is incomplete for {name}")
    mapped_rows = _map_external_rows(prediction.records, records)
    score_hash = hashlib.sha256(
        json.dumps(mapped_rows, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    run_identity = hashlib.sha256(
        f"{source_identity}|{weights_hash}|{score_hash}".encode()
    ).hexdigest()
    result_dir = target_dir / "results" / run_identity[:16]
    result_path = result_dir / "candidate_freeze.json"
    keys = residue_keys(records, model_id="primary_model_0")
    singles = build_single_evidence(
        records, keys, mapped_rows, settings, sequence=sequence, coverage_profile=True
    )
    compressed_options(records, singles, 3, coverage_profile=True)
    by_position = {}
    for item in singles:
        by_position.setdefault(item.residue.uniprot_pos, []).append(item)
    singles_table = []
    menus = {}
    for record in records:
        options = by_position.get(record.uniprot_pos, [])
        menu = []
        for item in options:
            descriptors = item.chemical_descriptors
            row = {
                "position": record.uniprot_pos,
                "wt": record.aa,
                "rsasa": record.rsasa,
                "mutant": item.mutant_aa,
                "spurs_ddg": item.spurs_ddg_kcal_mol,
                "stability_band": item.stability_band,
                "eligible": item.eligible,
                "hard_exclusions": item.hard_exclusions,
                "warnings": item.soft_warnings,
                "chemistry_class": item.chemistry_class,
                "descriptors": descriptors,
                "pareto_layer": item.pareto_layer,
                "menu_selected": item.menu_selected,
                "selection_reason": item.selection_reason,
            }
            singles_table.append(row)
            if item.menu_selected:
                menu.append(
                    {
                        "position": record.uniprot_pos,
                        "wt": record.aa,
                        "mutant": item.mutant_aa,
                        "notation": f"{record.aa}{record.uniprot_pos}{item.mutant_aa}",
                        "spurs_ddg": float(item.spurs_ddg_kcal_mol),
                        "exposure_band": "BOUNDARY_020_025"
                        if 0.20 <= record.rsasa < 0.25
                        else "EXPOSED_025_PLUS",
                        "mechanisms": sorted(descriptor_mechanisms(descriptors)),
                        "descriptors": descriptors,
                    }
                )
        if menu:
            menus[record.uniprot_pos] = menu
    surface_units, graph, areas = build_surface_units(records, rsasa_threshold=RSASA)
    pool, pool_counts = _pair_pool(
        records,
        graph,
        areas,
        menus,
        protein,
        mapping.uniprot_to_residue,
        input_hashes["coordinate_sha256"],
        run_identity,
    )
    strata = select_stratified(pool)
    selected = [item for group in strata.values() for item in group["selected"]]
    selected_positions = {x for row in selected for x in row["site_pair"]}
    result = {
        "schema": "cross_protein_candidate_freeze_v1",
        "status": "FROZEN",
        "target": name,
        "target_group": spec["group"],
        "structure": {
            "pdb": spec["pdb"],
            "chain": spec["chain"],
            "model_policy": "first deposited coordinate model only",
            "sha256": input_hashes["coordinate_sha256"],
            "canonical_mapping_entity": entity["entity_id"],
            "excluded_positions": ignored,
            "residue_mapping_method": "runtime map_chain_to_uniprot verified position-by-position against frozen author-residue canonical map",
        },
        "reference": {
            "accession": spec["accession"],
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "sequence_length": len(sequence),
            "study_range": spec["range"],
            "rSASA_threshold": RSASA,
            "uniprot_json_sha256": uniprot_hash,
            "features_used": "current frozen UniProt feature annotations plus coordinate-observed disulfides and sequence sequons",
        },
        "input_hashes": input_hashes,
        "code_hashes": code_hashes,
        "mapping_qc": {
            "observed_records": len(records),
            "runtime_identity": mapping.identity,
            "excluded_mismatch_positions": ignored,
            "analysis": analysis,
            "surface_graph_nodes": len(graph),
            "surface_graph_method": "build_surface_units residue C-alpha surface graph; proxy for pair scheduling",
        },
        "spurs_single": {
            "status": prediction.status,
            "runtime_seconds": prediction.runtime_seconds,
            "provenance": prediction.provenance,
            "warnings": prediction.warnings,
            "tool_identity": asdict(adapter.inventory()),
            "settings": settings,
            "sequence_length": len(prediction_sequence),
            "mutation_matrix_rows": len(singles_table),
        },
        "spurs_score_sha256": score_hash,
        "pair_generation": {
            **pool_counts,
            "selected_count": len(selected),
            "strata": {
                key: {
                    "status": value["status"],
                    "eligible_genotypes": value["eligible_genotypes"],
                    "distinct_site_pairs": value["distinct_site_pairs"],
                    "selected_candidate_ids": [x["candidate_id"] for x in value["selected"]],
                }
                for key, value in strata.items()
            },
            "pair_span_proxy_cutoff_a": PAIR_SPAN_A,
            "support_kernel_radius_proxy_a": SUPPORT_RADIUS_A,
            "candidate_pool_csv": "candidate_pool.csv",
            "complete_single_substitution_table": "single_substitutions.json",
        },
        "selected_candidates": selected,
        "interpretation": "stratified research candidate freeze only; single-site SPURS and graph proxies do not establish pair stability or formal perturbation qualification",
        "base_source_identity_sha256": source_identity,
        "source_identity_sha256": run_identity,
    }
    if result_path.is_file():
        existing = json.loads(result_path.read_text(encoding="utf-8"))
        if (
            existing.get("source_identity_sha256") != run_identity
            or existing.get("spurs_score_sha256") != score_hash
            or [x["candidate_id"] for x in existing.get("selected_candidates", [])]
            != [x["candidate_id"] for x in selected]
        ):
            raise ValueError(f"existing result identity collision: {result_path}")
        for artifact_name, digest in existing.get("artifact_hashes", {}).items():
            if (
                not (result_dir / artifact_name).is_file()
                or _sha(result_dir / artifact_name).casefold() != digest.casefold()
            ):
                raise ValueError(
                    f"existing frozen candidate artifact hash mismatch: {artifact_name}"
                )
        return {
            "target": name,
            "status": "REUSED_REVERIFIED_SCORE_AND_ARTIFACTS",
            "path": str(result_path),
            "result_sha256": _sha(result_path),
            "selected_count": len(existing.get("selected_candidates", [])),
            "pool_count": pool_counts["legal_genotype_count"],
        }
    _write_json(
        result_dir / "single_substitutions.json",
        {"rows": singles_table, "selected_menu_by_position": {str(k): v for k, v in menus.items()}},
    )
    pool_path = result_dir / "candidate_pool.csv"
    pool_temp = pool_path.with_suffix(".csv.part")
    with pool_temp.open("w", newline="", encoding="utf-8") as stream:
        fields = sorted({key for row in pool for key in row}) or ["candidate_id", "status"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in pool:
            writer.writerow(
                {
                    key: json.dumps(value, ensure_ascii=False)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )
    pool_temp.replace(pool_path)
    result_dir.mkdir(parents=True, exist_ok=True)
    prediction_copy = result_dir / "prediction_context.pdb"
    prediction_temp = prediction_copy.with_suffix(".pdb.part")
    prediction_temp.write_bytes(normalized.read_bytes())
    prediction_temp.replace(prediction_copy)
    result["artifact_hashes"] = {
        path.name: _sha(path)
        for path in [
            result_dir / "single_substitutions.json",
            prediction_copy,
            pool_path,
        ]
    }
    _write_json(result_path, result)
    return {
        "target": name,
        "status": "FROZEN",
        "path": str(result_path),
        "result_sha256": _sha(result_path),
        "selected_count": len(selected),
        "pool_count": len(pool),
        "selected_site_positions": sorted(selected_positions),
    }


def _freeze_ctla4_fixed(root: Path, output_root: Path) -> dict:
    historical_root = root / "historical_sources"
    source = historical_root / "ctla4_candidate_pool_v2.json"
    stability_source = historical_root / "ctla4_combination_stability.csv"
    for path in (source, stability_source):
        if not path.is_file() or _sha(path).casefold() != CTLA4_HISTORICAL_SOURCE_SHA256[path.name]:
            raise ValueError(f"frozen CTLA4 historical source missing or changed: {path.name}")
    data = json.loads(source.read_text(encoding="utf-8"))
    with stability_source.open(newline="", encoding="utf-8") as stream:
        stability_rows = {row["candidate_id"]: row for row in csv.DictReader(stream)}
    by_mutations = {
        tuple(sorted(item["notation"] for item in row["mutations"])): row
        for row in data["candidates"]
    }
    selected = []
    for notation in CTLA4_FIXED:
        key = tuple(sorted(notation))
        row = by_mutations.get(key)
        if row is None:
            raise ValueError(f"fixed CTLA4 candidate missing from frozen historical pool: {key}")
        pair_row = stability_rows.get(row["candidate_id"])
        if pair_row is None:
            raise ValueError(f"historical CTLA4 pair stability record missing: {row['candidate_id']}")
        selected.append(
            {
                "candidate_id": row["candidate_id"],
                "mutations": list(key),
                "risk_components": row.get("risk_components", {}),
                "single_spurs": [item.get("spurs_ddg_kcal_mol") for item in row["mutations"]],
                "historical_pair_stability_screen": {
                    "spurs_multi_ddg_kcal_mol": float(pair_row["spurs_multi_ddg_kcal_mol"]),
                    "foldx_ddg_kcal_mol": float(pair_row["foldx_ddg_kcal_mol"]),
                    "order": int(pair_row["order"]),
                    "epistasis_estimate": float(pair_row["epistasis_estimate"]),
                    "extrapolation_level": pair_row["extrapolation_level"],
                    "historical_consensus_status": pair_row["consensus_status"],
                    "disagreement": float(pair_row["disagreement"]),
                    "qualification": "HISTORICAL_SCREEN_ONLY_NOT_REQUALIFIED",
                },
                "source_candidate_pool_sha256": _sha(source),
                "source_pair_stability_sha256": _sha(stability_source),
                "status": "FIXED_EXISTING_CANDIDATE_NOT_REQUALIFIED",
            }
        )
    payload = {
        "schema": "cross_protein_candidate_freeze_v1",
        "status": "FROZEN_FIXED_SET",
        "target": "CTLA4",
        "target_group": "development",
        "candidates": selected,
        "nmr_ensemble_sensitivity": _ctla4_nmr_sensitivity(root),
        "coupling_rule": "surface_perturbation.engine_v3.coupling applied to frozen candidate risk components",
        "source": str(source),
        "source_sha256": _sha(source),
        "pair_stability_source": str(stability_source),
        "pair_stability_source_sha256": _sha(stability_source),
        "historical_sources": {path.name: _sha(path) for path in (source, stability_source)},
        "interpretation": "fixed six-genotype comparison set; existing IDs and evidence retained without requalification",
    }
    path = output_root / "CTLA4" / "fixed_prior_set" / "candidate_freeze.json"
    _write_json(path, payload)
    return {
        "target": "CTLA4",
        "status": "FROZEN_FIXED_SET",
        "path": str(path),
        "selected_count": len(selected),
    }


def _ctla4_nmr_sensitivity(root: Path) -> dict:
    source = root / "inputs" / "1AH1.cif"
    maps = json.loads(
        (root / "mapping" / "canonical_residue_maps.json").read_text(encoding="utf-8")
    )
    entity = next(
        item
        for item in maps["structures"]["1AH1"]["entity_mappings"]
        if item["target"] == "CTLA4" and "A" in item["chains"]
    )
    rows = _rows_by_author_id(entity["residues"], "A")
    ignored = set(entity.get("candidate_freeze_excluded_positions", [])) | {147}
    structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(
        "1AH1", str(source)
    )
    models = list(structure.get_models())
    coordinate_sets = []
    for model in models:
        positions = {}
        if "A" not in model:
            continue
        for residue in model["A"]:
            row = rows.get((int(residue.id[1]), str(residue.id[2]).strip()))
            if row is None or row["reference_position"] is None:
                continue
            position = int(row["reference_position"])
            if 37 <= position <= 161 and position not in ignored and residue.has_id("CA"):
                positions[position] = np.asarray(residue["CA"].coord, dtype=float)
        coordinate_sets.append(positions)
    common = (
        set.intersection(*(set(item) for item in coordinate_sets)) if coordinate_sets else set()
    )
    common_positions = sorted(common)
    rmsds = []
    for left, right in combinations(coordinate_sets, 2):
        a = np.asarray([left[pos] for pos in common_positions])
        b = np.asarray([right[pos] for pos in common_positions])
        if len(a) < 3:
            continue
        ac = a - a.mean(axis=0)
        bc = b - b.mean(axis=0)
        u, _s, vt = np.linalg.svd(ac.T @ bc)
        correction = np.eye(3)
        correction[2, 2] = np.linalg.det(u @ vt)
        rotation = u @ correction @ vt
        rmsds.append(float(np.sqrt(np.mean(np.sum((ac @ rotation - bc) ** 2, axis=1)))))
    return {
        "coordinate_model_count": len(models),
        "pairwise_comparison_count": len(rmsds),
        "common_ca_residue_count": len(common_positions),
        "common_positions": common_positions,
        "pairwise_kabsch_ca_rmsd_a": {
            "minimum": min(rmsds) if rmsds else None,
            "median": float(np.median(rmsds)) if rmsds else None,
            "p90": float(np.quantile(rmsds, 0.9)) if rmsds else None,
            "maximum": max(rmsds) if rmsds else None,
        },
        "interpretation": "intra-entry conformer sensitivity only; the 20 models are not independent proteins or thermodynamic samples",
    }


def freeze(root: Path = DEFAULT_ROOT, targets: list[str] | None = None) -> dict:
    root = root.resolve()
    input_root = root / "inputs"
    output_root = root / "candidate_freeze"
    maps_path = root / "mapping" / "canonical_residue_maps.json"
    map_manifest_path = root / "mapping" / "canonical_mapping_manifest.json"
    maps = json.loads(maps_path.read_text(encoding="utf-8"))
    map_manifest = json.loads(map_manifest_path.read_text(encoding="utf-8"))
    input_manifest = json.loads((input_root / "input_manifest.json").read_text(encoding="utf-8"))
    if _sha(maps_path).casefold() != map_manifest["canonical_residue_maps_sha256"].casefold():
        raise ValueError("canonical mapping manifest hash mismatch")
    if (
        maps.get("input_manifest_sha256", "").casefold()
        != _sha(input_root / "input_manifest.json").casefold()
    ):
        raise ValueError("canonical map does not bind to the current frozen input manifest")
    mapped_digest = hashlib.sha256(
        json.dumps(maps["structures"], sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if mapped_digest.casefold() != maps.get("mapping_sha256", "").casefold():
        raise ValueError("canonical mapping content digest mismatch")
    for filename, expected in map_manifest.get("sequence_sources", {}).items():
        path = root / "mapping" / "sequence_sources" / filename
        if not path.is_file() or _sha(path).casefold() != expected.casefold():
            raise ValueError(f"frozen UniProt FASTA hash mismatch: {filename}")
    for name, expected in input_manifest["files"].items():
        if _sha(input_root / name).casefold() != expected.casefold():
            raise ValueError(f"frozen cross-protein input hash mismatch: {name}")
    output_root.mkdir(parents=True, exist_ok=True)
    results = [_freeze_ctla4_fixed(root, output_root)]
    payload = {
        "schema": "cross_protein_candidate_freeze_manifest_v1",
        "profile": "candidate_freeze_v1",
        "input_manifest_sha256": _sha(input_root / "input_manifest.json"),
        "canonical_maps_sha256": _sha(maps_path),
        "status": "PARTIAL",
        "results": results,
    }
    _write_json(output_root / "freeze_manifest.json", payload)
    names = targets or [name for name in SPECS if name != "CTLA4"]
    for name in names:
        if name == "CTLA4":
            continue
        if name not in SPECS:
            raise ValueError(f"unknown cross-protein target: {name}")
        try:
            result = _freeze_target(
                name, SPECS[name], root, maps, root / "mapping" / "sequence_sources", output_root
            )
            results.append(result)
            payload["results"] = results
            _write_json(output_root / "freeze_manifest.json", payload)
        except Exception as exc:
            attempt = {
                "target": name,
                "status": "FAILED",
                "error": f"{type(exc).__name__}: {exc}",
                "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            with (output_root / "attempts.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(attempt, ensure_ascii=False) + "\n")
            raise
    payload["status"] = "COMPLETE" if len(results) == 8 else "PARTIAL"
    _write_json(output_root / "freeze_manifest.json", payload)
    return payload


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--targets", nargs="*", choices=sorted(SPECS))
    args = parser.parse_args()
    result = freeze(args.root, args.targets)
    print(
        json.dumps({"status": result["status"], "results": result["results"]}, ensure_ascii=False)
    )


if __name__ == "__main__":
    main()
