"""Build a structure-only Fab/protein contact-scale reference for the CTLA4 study."""
from __future__ import annotations

import copy
import csv
import hashlib
import http.client
import json
import math
import shutil
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import Bio
import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, Model, ShrakeRupley, Structure
from Bio.PDB.Polypeptide import is_aa
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from ..util import sha256_file, write_json
from ..surface_perturbation.tools import PhysicalTools

DEFAULT_STUDY = Path(r"E:\user\New project\ctla4_field_fusion_20261006")
SEED = 731029
TARGET_CLUSTERS = 30
MIN_CLUSTERS = 20
IDENTITY_CUTOFF = 0.30
COVERAGE_CUTOFF = 0.80
FAB_TYPES = {"FAB"}
EXPERIMENTAL_METHODS = {"XRAY", "ELECTRON_MICROSCOPY", "SOLID_STATE_NMR", "SOLUTION_NMR", "ELECTRON_CRYSTALLOGRAPHY"}
CONTACT_CUTOFFS = (4.0, 4.5, 5.0)


def _digest_matches(path: Path, expected: str) -> bool:
    return sha256_file(path).casefold() == str(expected).casefold()
SURFACE_CONTACT_CUTOFF = 4.5
SURFACE_DISTANCE_SAMPLE = 128


def _sequence(chain) -> tuple[str, list]:
    residues = [r for r in chain if is_aa(r, standard=False)]
    return "".join(protein_letters_3to1_extended.get(r.resname.upper(), "X") for r in residues), residues


def _identity_coverage(a: str, b: str) -> tuple[float, float]:
    """Global sequence identity and shorter-chain coverage over aligned columns."""
    if not a or not b:
        return 0.0, 0.0
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner(); aligner.mode = "global"; aligner.match_score = 2.; aligner.mismatch_score = -1.
    aligner.open_gap_score = -3.; aligner.extend_gap_score = -.5
    alignment = aligner.align(a, b)[0]
    paired = valid = identical = 0
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        for x, y in zip(a[a0:a1], b[b0:b1]):
            paired += 1
            if x != "X" and y != "X":
                valid += 1; identical += x == y
    return (identical / valid if valid else 0.0, paired / min(len(a), len(b)))


def _homologous(sequences_a: tuple[str, ...], sequences_b: tuple[str, ...]) -> bool:
    for a in sequences_a:
        for b in sequences_b:
            identity, coverage = _identity_coverage(a, b)
            if identity >= IDENTITY_CUTOFF and coverage >= COVERAGE_CUTOFF:
                return True
    return False


def _heavy(atom) -> bool:
    element = str(atom.element or "").strip().upper()
    return element not in {"H", "D"} and not atom.get_name().strip().lstrip("0123456789").upper().startswith(("H", "D"))


def _chain_atoms(model, chain_ids: tuple[str, ...]) -> tuple[list, dict]:
    rows = []; atom_map = {}
    for chain_id in chain_ids:
        if chain_id not in model:
            raise ValueError(f"author chain absent from selected model: {chain_id}")
        chain = model[chain_id]
        for residue in chain:
            if not is_aa(residue, standard=False):
                continue
            reskey = (chain_id, int(residue.id[1]), str(residue.id[2]).strip())
            for atom in residue:
                if not _heavy(atom):
                    continue
                name = atom.get_name().strip(); element = str(atom.element or name[0]).strip().upper()
                coord = np.asarray(atom.coord, float)
                rows.append((coord, reskey, name, residue.resname, element))
                atom_map[(chain_id, reskey[1], reskey[2], name)] = (coord, residue.resname, element)
    if not rows:
        raise ValueError("no protein heavy atoms in selected chains")
    return rows, atom_map


def _selected_structure(model, chain_ids: tuple[str, ...]):
    structure = Structure.Structure("contact")
    model_copy = Model.Model(0); structure.add(model_copy)
    residues = {}
    for chain_id in chain_ids:
        if chain_id not in model:
            raise ValueError(f"chain absent for SASA: {chain_id}")
        chain = copy.deepcopy(model[chain_id])
        for residue in list(chain):
            if not is_aa(residue, standard=False):
                chain.detach_child(residue.id); continue
            for atom in list(residue):
                if not _heavy(atom): residue.detach_child(atom.id)
            residues[(chain_id, residue.id)] = residue
        model_copy.add(chain)
    return structure, residues


def _residue_sasa(model, chains: tuple[str, ...]) -> dict:
    structure, _ = _selected_structure(model, chains)
    ShrakeRupley(n_points=960).compute(structure, level="R")
    return {(chain.id, residue.id): float(residue.sasa)
            for chain in structure[0] for residue in chain if hasattr(residue, "sasa")}


def _interface_axes(coords: np.ndarray) -> dict:
    if len(coords) < 2:
        return {"long_axis_a": 0.0, "short_axis_a": 0.0, "contact_residue_ca_count": int(len(coords))}
    centered = coords - coords.mean(axis=0)
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    spans = np.ptp(centered @ axes.T, axis=0)
    ordered = np.sort(spans)[::-1]
    return {"long_axis_a": float(ordered[0]), "short_axis_a": float(ordered[1]),
            "contact_residue_ca_count": int(len(coords))}


def _sample_contact_vertices(contact: np.ndarray, vertex_area: np.ndarray, atom_map: dict) -> tuple[np.ndarray, bool]:
    contact = np.asarray(contact, int)
    weights = np.asarray(vertex_area, float)[contact]
    valid = np.isfinite(weights) & (weights > 0)
    contact, weights = contact[valid], weights[valid]
    if len(contact) <= SURFACE_DISTANCE_SAMPLE:
        return contact, True
    identity = hashlib.sha256()
    for key, (xyz, _residue, _element) in sorted(atom_map.items()):
        identity.update(repr(key).encode())
        identity.update(np.asarray(xyz, dtype="<f8").tobytes())
    seed = int.from_bytes(identity.digest()[:8], "big") ^ SEED
    sample = np.random.default_rng(seed).choice(contact, size=SURFACE_DISTANCE_SAMPLE,
                                                replace=False, p=weights / weights.sum())
    return sample, False


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    order = np.argsort(values[valid], kind="stable")
    sorted_values, sorted_weights = values[valid][order], weights[valid][order]
    index = min(int(np.searchsorted(np.cumsum(sorted_weights), quantile * sorted_weights.sum(), side="left")),
                len(sorted_values) - 1)
    return float(sorted_values[index])


def _surface_geodesics(tools, atom_map: dict, fab_xyz: np.ndarray) -> dict:
    try:
        mesh = tools.ses(atom_map, scale=2.0)
    except RuntimeError as exc:
        if str(exc) == "SES mesh is not closed":
            return {"status": "UNKNOWN_OPEN_SURFACE", "reason": str(exc)}
        raise
    if not mesh.closed:
        return {"status": "UNKNOWN_OPEN_SURFACE", "mesh_vertices": int(len(mesh.vertices))}
    faces = mesh.faces
    # Sort each edge's endpoints. Sorting axis=0 silently made self-edges,
    # leaving every mesh vertex disconnected while still producing "closed" SES.
    edges = np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1), axis=0)
    lengths = np.linalg.norm(mesh.vertices[edges[:, 0]] - mesh.vertices[edges[:, 1]], axis=1)
    if not len(edges) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        return {"status": "UNKNOWN_INVALID_SURFACE_GRAPH", "mesh_vertices": int(len(mesh.vertices)),
                "mesh_edges": int(len(edges)), "reason": "surface graph contains no edges or nonpositive/nonfinite edge lengths"}
    rows = np.r_[edges[:, 0], edges[:, 1]]; cols = np.r_[edges[:, 1], edges[:, 0]]
    graph = coo_matrix((np.r_[lengths, lengths], (rows, cols)), shape=(len(mesh.vertices),) * 2).tocsr()
    n_components, component = connected_components(graph, directed=False)
    distance, _ = cKDTree(fab_xyz).query(mesh.vertices, k=1)
    contact = np.flatnonzero(distance <= SURFACE_CONTACT_CUTOFF)
    vertex_area = np.zeros(len(mesh.vertices), float)
    np.add.at(vertex_area, faces.reshape(-1), np.repeat(np.asarray(mesh.area, float) / 3., 3))
    contact = contact[np.isfinite(vertex_area[contact]) & (vertex_area[contact] > 0)]
    if len(contact) < 2:
        return {"status": "UNKNOWN_TOO_FEW_CONTACT_SURFACE_POINTS", "mesh_vertices": int(len(mesh.vertices)),
                "contact_surface_vertices": int(len(contact)), "surface_components": int(n_components)}
    sample, complete_contact_mesh = _sample_contact_vertices(contact, vertex_area, atom_map)
    distances = dijkstra(graph, directed=False, indices=sample)
    pair_distances = distances[:, sample][np.triu_indices(len(sample), k=1)]
    finite = pair_distances[np.isfinite(pair_distances)]
    total_pairs = len(pair_distances)
    if complete_contact_mesh:
        pair_weights = np.outer(vertex_area[sample], vertex_area[sample])[np.triu_indices(len(sample), k=1)]
        connected_fraction = float(pair_weights[np.isfinite(pair_distances)].sum() / pair_weights.sum()) if pair_weights.sum() else 0.0
        quantiles = {str(q): _weighted_quantile(pair_distances, pair_weights, q) for q in (.5, .9, .95)}
        quantile_method = "exact vertex-area product weights over all contact vertices"
    else:
        connected_fraction = float(len(finite) / max(total_pairs, 1))
        quantiles = {str(q): float(np.quantile(finite, q)) for q in (.5, .9, .95)} if len(finite) else {}
        quantile_method = "128 deterministic area-proportional contact-vertex draws without replacement; unweighted pair quantiles"
    if not len(finite) or not np.isfinite(list(quantiles.values())).all() or any(value <= 0 for value in quantiles.values()):
        return {"status": "UNKNOWN_INVALID_GEODESIC_DISTANCES", "mesh_vertices": int(len(mesh.vertices)),
                "contact_surface_vertices": int(len(contact)), "sampled_contact_points": int(len(sample)),
                "connected_pair_fraction": connected_fraction, "surface_components": int(n_components),
                "reason": "resolved contact distances must be finite and strictly positive"}
    componentwise = connected_fraction < 1.0 - 1e-12
    return {"status": "RESOLVED_COMPONENTWISE" if componentwise else "RESOLVED",
            "distance_scope": "within connected antigen SES components; cross-component paths are undefined" if componentwise
            else "all sampled contact vertices are connected on the antigen SES mesh",
            "mesh_vertices": int(len(mesh.vertices)),
            "contact_surface_vertices": int(len(contact)), "sampled_contact_points": int(len(sample)),
            "contact_surface_area_a2": float(vertex_area[contact].sum()),
            "connected_pair_fraction": connected_fraction,
            "interface_surface_components": int(len(set(component[contact]))),
            "pairwise_surface_distance_quantiles_a": quantiles,
            "pairwise_surface_distance_max_a": float(finite.max()),
            "quantile_method": quantile_method,
            "interpretation": "geodesic along the isolated antigen SES mesh; disconnected components use within-component pairs only; structural scale prior, not experimental response calibration"}


def _resolved_geodesic(sample: dict) -> bool:
    return sample.get("geodesic", {}).get("status") in {"RESOLVED", "RESOLVED_COMPONENTWISE"}


def measure_complex(model, fab_chains: tuple[str, ...], antigen_chains: tuple[str, ...], tools) -> dict:
    fab_rows, _ = _chain_atoms(model, fab_chains)
    antigen_rows, atom_map = _chain_atoms(model, antigen_chains)
    fab_xyz = np.asarray([row[0] for row in fab_rows], float)
    antigen_xyz = np.asarray([row[0] for row in antigen_rows], float)
    neighbor_pairs = cKDTree(fab_xyz).query_ball_tree(cKDTree(antigen_xyz), max(CONTACT_CUTOFFS))
    pair_rows = []
    for fi, ids in enumerate(neighbor_pairs):
        for ai in ids:
            distance = float(np.linalg.norm(fab_xyz[fi] - antigen_xyz[ai]))
            f_res, a_res = fab_rows[fi][1], antigen_rows[ai][1]
            pair_rows.append((distance, f_res, a_res))
    contact_counts = {}; contact_residues = set(); fab_contact_residues = set()
    for cutoff in CONTACT_CUTOFFS:
        selected = [(d, f, a) for d, f, a in pair_rows if d <= cutoff]
        contact_counts[str(cutoff)] = {"heavy_atom_pairs": len(selected),
            "antigen_residues": len({a for _, _, a in selected}), "fab_residues": len({f for _, f, _ in selected})}
        if cutoff == 4.5:
            contact_residues = {a for _, _, a in selected}; fab_contact_residues = {f for _, f, _ in selected}
    if not contact_residues:
        raise ValueError("no Fab-protein heavy-atom contacts at 4.5 A")
    fab_free = _residue_sasa(model, fab_chains); antigen_free = _residue_sasa(model, antigen_chains)
    complex_free = _residue_sasa(model, tuple(dict.fromkeys((*fab_chains, *antigen_chains))))
    antigen_buried = sum(max(0.0, value - complex_free.get(key, value)) for key, value in antigen_free.items())
    fab_buried = sum(max(0.0, value - complex_free.get(key, value)) for key, value in fab_free.items())
    ca = []
    for chain_id, number, insertion in contact_residues:
        if chain_id in model and (" ", number, insertion or " ") in model[chain_id]:
            residue = model[chain_id][(" ", number, insertion or " ")]
            if "CA" in residue: ca.append(np.asarray(residue["CA"].coord, float))
    interface = _interface_axes(np.asarray(ca, float).reshape(-1, 3))
    geodesic = _surface_geodesics(tools, atom_map, fab_xyz)
    return {"fab_chain_ids": list(fab_chains), "protein_antigen_chain_ids": list(antigen_chains),
            "contact_cutoff_counts": contact_counts, "contact_residue_keys_4_5a": [list(x) for x in sorted(contact_residues)],
            "fab_contact_residue_count_4_5a": len(fab_contact_residues),
            "buried_surface_area_antigen_a2": float(antigen_buried),
            "buried_surface_area_fab_a2": float(fab_buried),
            "buried_surface_area_total_half_a2": float((antigen_buried + fab_buried) / 2),
            "sasa_method": "Bio.PDB Shrake-Rupley; heavy atoms only; 960 sphere points",
            "contact_residue_ca_spans": interface, "geodesic": geodesic,
            "axis_interpretation": "PCA spans of contacting antigen-residue CA atoms; not the contact-surface boundary"}


def _read_frozen_inputs(source_dir: Path) -> tuple[list[dict], dict]:
    with (source_dir / "sabdab_summary_all.csv").open(newline="", encoding="utf-8-sig") as handle:
        summary = list(csv.DictReader(handle))
    annotations = json.loads((source_dir / "sabdab_rcsb_annotations.json").read_text(encoding="utf-8"))
    by_key = defaultdict(list)
    for row in annotations:
        by_key[(str(row["PDB ID"]), int(row["model"]))].append(row)
    candidates = []
    for row in summary:
        if row.get("type") not in FAB_TYPES or row.get("antigen_type") != "PROTEIN": continue
        if row.get("method") not in EXPERIMENTAL_METHODS or not row.get("Hchain") or not row.get("Lchain"):
            continue
        antigen_labels = set(str(row.get("antigen_chain") or "").split("|"))
        if not antigen_labels - {"", "nan"}: continue
        for annotation in by_key.get((row["PDB"], int(float(row["model"]))), []):
            heavy = annotation.get("heavy chain") or {}; light = annotation.get("light chain") or {}
            if heavy.get("SAbDab auth_asym_id") != row["Hchain"] or light.get("SAbDab auth_asym_id") != row["Lchain"]:
                continue
            antigen = [x for x in annotation.get("antigen_instances", [])
                       if x.get("entity type") == "PROTEIN" and x.get("SAbDab auth_asym_id") in antigen_labels]
            if not antigen: continue
            candidates.append({"instance": row["INSTANCE"], "pdb_id": row["PDB"],
                "pdb_code": row["PDB"][-4:].lower(), "model": int(float(row["model"])),
                "fab_chains": tuple(dict.fromkeys([heavy["PDB auth_asym_id"], light["PDB auth_asym_id"]])),
                "antigen_chains": tuple(dict.fromkeys(x["PDB auth_asym_id"] for x in antigen)),
                "antigen_sequences": (), "antigen_name": row.get("antigen_name", ""),
                "method": row.get("method"), "resolution_a": _optional_float(row.get("resolution")),
                "row": row})
            break
    return candidates, by_key


def _write_incomplete_contact_scale_diagnostic(output_dir: Path, *, target_clusters: int, minimum_clusters: int,
                                               candidates: int, accepted: list[dict], excluded: list[dict],
                                               source_manifest_sha256: str, tools, started: float) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "contact_scale_attempt_diagnostic.json"
    diagnostic = {
        "profile": "contact_scale_reference_v1",
        "status": "INSUFFICIENT_SEQUENCE_CLUSTERS",
        "seed": SEED,
        "target_cluster_count": target_clusters,
        "minimum_required": minimum_clusters,
        "minimum_policy": {"default_minimum": MIN_CLUSTERS, "effective_minimum": minimum_clusters,
                           "override_used": minimum_clusters != MIN_CLUSTERS,
                           "scope": "exploratory structural prior only; not an assay-calibrated sample-size threshold"},
        "selected_structure_count": len(accepted),
        "candidate_population_mapping_valid": candidates,
        "source_manifest_sha256": source_manifest_sha256,
        "homology_rule": {
            "identity": IDENTITY_CUTOFF,
            "shorter_sequence_coverage": COVERAGE_CUTOFF,
            "method": "deterministic global alignment; any antigen-chain pair joins a cluster",
        },
        "selection_rule": "SHA256(seed|SAbDab INSTANCE) ascending; first eligible unique antigen sequence clusters; no contact-result ranking",
        "samples": accepted,
        "excluded_attempts": excluded,
        "software_identity": {
            "numpy_version": np.__version__,
            "biopython_version": Bio.__version__,
            "nanoshaper_sha256": tools.identity(tools.nano),
        },
        "field_costs": {
            "wall_seconds": float(time.perf_counter() - started),
            "calls": tools.calls,
            "activity_seconds": float(sum(row.get("runtime_seconds", 0.) for row in tools.activity)),
            "activity": tools.activity,
            "source_snapshot_sha256": sha256_file(Path(__file__).resolve()),
        },
        "interpretation": "diagnostic only; below the frozen minimum, so no contact-scale reference or radius scenarios are emitted",
    }
    write_json(path, diagnostic)
    return path


def _optional_float(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _download_cif(candidate: dict, folder: Path) -> Path:
    path = folder / f"{candidate['pdb_code']}.cif"
    if not path.exists():
        if shutil.disk_usage(folder.anchor).free < 10 * 1024**3:
            raise RuntimeError("RESOURCE_STOP: less than 10 GiB free before contact-structure download")
        request = urllib.request.Request(f"https://files.rcsb.org/download/{candidate['pdb_code'].upper()}.cif",
                                         headers={"User-Agent": "ctla4-field-fusion-shadow/1.0"})
        partial = path.with_suffix(".cif.part")
        for attempt in range(4):
            try:
                with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
                if partial.stat().st_size < 1000:
                    raise ValueError(f"invalid mmCIF download: {candidate['pdb_code']}")
                partial.replace(path)
                break
            except (http.client.HTTPException, urllib.error.URLError, TimeoutError, OSError):
                partial.unlink(missing_ok=True)
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
    if path.stat().st_size < 1000:
        raise ValueError(f"invalid mmCIF download: {candidate['pdb_code']}")
    return path


def _target_sequences(candidate: dict, model) -> tuple[str, ...]:
    sequences = []
    for chain_id in candidate["antigen_chains"]:
        seq, _ = _sequence(model[chain_id])
        if seq: sequences.append(seq)
    return tuple(dict.fromkeys(sequences))


def _ctla4_sequence(study: Path) -> tuple[str, list[int]]:
    data = json.loads((study / "inputs" / "analysis.json").read_text(encoding="utf-8"))
    residues = sorted(data["residues"], key=lambda row: int(row["uniprot_pos"]))
    return "".join(row["aa"] for row in residues), [int(row["uniprot_pos"]) for row in residues]


def _replay_frozen_cohort(selection_reference: Path, candidates: list[dict], cif_dir: Path,
                          target_sequence: str, source_manifest_sha256: str, tools) -> list[dict]:
    prior = json.loads(Path(selection_reference).read_text(encoding="utf-8"))
    if (prior.get("profile") != "contact_scale_reference_v1"
            or prior.get("source_manifest_sha256", "").casefold() != source_manifest_sha256.casefold()):
        raise ValueError("frozen cohort selection does not match this contact input manifest")
    samples = prior.get("samples", [])
    if len(samples) != int(prior.get("selected_structure_count", -1)):
        raise ValueError("frozen cohort selection has an inconsistent sample count")
    by_instance = {candidate["instance"]: candidate for candidate in candidates}
    accepted = []; cluster_sequences = []; current_structure = None; current_pdb = None
    for order, saved in enumerate(samples):
        instance = saved.get("instance")
        candidate = by_instance.get(instance)
        if candidate is None:
            raise ValueError(f"frozen cohort instance is absent from current input: {instance}")
        expected_rank = hashlib.sha256(f"{SEED}|{instance}".encode()).hexdigest()
        if (saved.get("rank_sha256") != expected_rank or saved.get("pdb_code") != candidate["pdb_code"]
                or int(saved.get("model", -1)) != candidate["model"]
                or tuple(saved.get("fab_chains", ())) != candidate["fab_chains"]
                or tuple(saved.get("antigen_chains", ())) != candidate["antigen_chains"]):
            raise ValueError(f"frozen cohort structure identity changed: {instance}")
        cif = _download_cif(candidate, cif_dir)
        if not _digest_matches(cif, saved.get("coordinate_file_sha256", "")):
            raise ValueError(f"frozen cohort coordinate hash changed: {instance}")
        if current_pdb != candidate["pdb_code"]:
            current_structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(
                candidate["pdb_code"], str(cif))
            current_pdb = candidate["pdb_code"]
        model = next((item for item in current_structure if int(item.id) == candidate["model"]), None)
        if model is None:
            raise ValueError(f"frozen cohort model is absent: {instance}")
        sequences = _target_sequences(candidate, model)
        expected_sequences = tuple(saved.get("protein_antigen_sequences", ()))
        if sequences != expected_sequences or _homologous(sequences, (target_sequence,)):
            raise ValueError(f"frozen cohort antigen sequence changed or maps to CTLA4: {instance}")
        if any(_homologous(sequences, prior_sequences) for prior_sequences in cluster_sequences):
            raise ValueError(f"frozen cohort is no longer sequence-distinct: {instance}")
        measured = measure_complex(model, candidate["fab_chains"], candidate["antigen_chains"], tools)
        sample = {key: value for key, value in candidate.items() if key != "row"}
        sample.update(measured)
        sample.update(coordinate_file_sha256=sha256_file(cif), protein_antigen_sequences=list(sequences),
                      selection_order=int(saved["selection_order"]), cluster_id=saved["cluster_id"])
        accepted.append(sample); cluster_sequences.append(sequences)
        print(f"contact-scale replayed frozen {order+1}/{len(samples)}: {instance}", flush=True)
    return accepted


def _ctla4_external_candidates(candidates: list[dict]) -> list[dict]:
    words = ("ctla-4", "cytotoxic t-lymphocyte protein 4", "cd152")
    seen = set(); output = []
    for candidate in candidates:
        text = " ".join(str(candidate["row"].get(key, "")) for key in ("antigen_name", "compound")).lower()
        identity = (candidate["pdb_code"], candidate["fab_chains"], candidate["antigen_chains"])
        if identity not in seen and any(word in text for word in words):
            seen.add(identity); output.append(candidate)
    return output


def _complex_residue_to_uniprot(candidate: dict, model, target_sequence: str, target_positions: list[int]) -> set[int]:
    from Bio.Align import PairwiseAligner
    from Bio.Data.PDBData import protein_letters_3to1_extended
    from Bio.PDB.Polypeptide import is_aa
    aligner = PairwiseAligner(); aligner.mode = "global"; aligner.match_score = 2.; aligner.mismatch_score = -1.
    aligner.open_gap_score = -3.; aligner.extend_gap_score = -.5
    target = set()
    for chain_id in candidate["antigen_chains"]:
        residues = [r for r in model[chain_id] if is_aa(r, standard=False)]
        seq = "".join(protein_letters_3to1_extended.get(r.resname.upper(), "X") for r in residues)
        alignment = aligner.align(target_sequence, seq)[0]
        mapping = {}
        for (t0, t1), (q0, q1) in zip(*alignment.aligned):
            for ti, qi in zip(range(t0, t1), range(q0, q1)):
                if target_sequence[ti] != "X" and seq[qi] != "X": mapping[qi] = target_positions[ti]
        for residue_index, residue in enumerate(residues):
            if (chain_id, int(residue.id[1]), str(residue.id[2]).strip()) in candidate["contact_residue_keys"] and residue_index in mapping:
                target.add(mapping[residue_index])
    return target


def _sample_mask_coverage(mask: np.ndarray, owner: np.ndarray, area: np.ndarray, positions: set[int]) -> dict:
    interface = np.isin(owner, list(positions))
    total = float(area[interface].sum())
    affected = float(area[interface & mask].sum())
    return {"mapped_residue_count": len(positions), "reference_interface_area_a2": total,
            "overlap_area_a2": affected, "overlap_fraction": affected / total if total else None,
            "interpretation": "sequence-mapped geometric intersection only; not binding-loss evidence"}


def _dilate(mask: np.ndarray, graph, radius: float) -> np.ndarray:
    seeds = np.flatnonzero(mask)
    if not len(seeds): return np.zeros(len(mask), bool)
    distances = np.asarray(dijkstra(graph, indices=seeds, min_only=True), float).reshape(-1)
    return distances <= radius


def _opportunity(mask: np.ndarray, unknown: np.ndarray, owner: np.ndarray, area: np.ndarray, graph,
                 radius: float) -> dict:
    lower = _dilate(mask, graph, radius)
    unknown_dilated = _dilate(unknown, graph, radius) & ~lower
    denominator = float(area.sum())
    lower_area = float(area[lower].sum()); unknown_area = float(area[unknown_dilated].sum())
    return {"radius_a": float(radius), "lower_area_a2": lower_area,
            "lower_fraction": lower_area / denominator if denominator else None,
            "unknown_opportunity_area_a2": unknown_area,
            "unknown_opportunity_fraction": unknown_area / denominator if denominator else None,
            "possible_upper_fraction": (lower_area + unknown_area) / denominator if denominator else None}


def _candidate_mask_overlaps(candidate_ids: list[str], positions: set[int], owner: np.ndarray,
                             area: np.ndarray, geometry_root: Path, apbs_root: Path) -> dict:
    overlaps = {}
    for cid in candidate_ids:
        geometry_path = geometry_root / cid / "geometry_fusion_masks.npz"
        if not geometry_path.is_file(): continue
        with np.load(geometry_path, allow_pickle=False) as masks:
            geometry_mask = masks["geometry_chemistry__positive"].copy()
            geometry_unknown = masks["geometry_chemistry__unknown"].copy()
        overlaps[cid] = {"geometry_chemistry": _sample_mask_coverage(geometry_mask, owner, area, positions)}
        overlaps[cid]["geometry_chemistry_unknown_opportunity"] = _sample_mask_coverage(
            geometry_unknown, owner, area, positions)
        fused_path = apbs_root / cid / "field_fusion_evidence.npz"
        if fused_path.is_file():
            with np.load(fused_path, allow_pickle=False) as masks:
                if "three_mechanism_positive" in masks.files:
                    overlaps[cid]["three_mechanism"] = _sample_mask_coverage(
                        masks["three_mechanism_positive"].copy(), owner, area, positions)
                if "three_mechanism_unknown" in masks.files:
                    overlaps[cid]["three_mechanism_unknown_opportunity"] = _sample_mask_coverage(
                        masks["three_mechanism_unknown"].copy(), owner, area, positions)
    return overlaps


def _apply_ctla4_contact_overlap(study: Path, accepted: list[dict], source_dir: Path,
                                 candidates: list[dict], tools,
                                 precomputed_external_results: list[dict] | None = None) -> dict:
    from scipy.sparse import load_npz
    manifest = json.loads((study / "study_manifest.json").read_text(encoding="utf-8"))
    reference = study / "inputs" / "canonical_reference.npz"
    with np.load(reference, allow_pickle=False) as data:
        area = data["area"].copy(); owner = data["owner"].copy()
    graph = load_npz(study / "inputs" / "reference_graph.npz").tocsr()
    target_sequence, target_positions = _ctla4_sequence(study)
    external = _ctla4_external_candidates(candidates)
    out = study / "evidence" / "contact_scale" / "results"
    versioned = (study / "evidence" / "null_models" / "independent_null_v1.json").is_file()
    geometry_root = study / "evidence" / "geometry" / "independent_null_v1" if versioned else study / "evidence" / "geometry"
    apbs_root = study / "evidence" / "apbs" / "independent_null_v1" if versioned else study / "evidence" / "apbs"
    source_structures = source_dir / "structures"
    external_results = []
    if precomputed_external_results is not None:
        for row in precomputed_external_results:
            item = dict(row)
            mapped = set(int(value) for value in item["mapped_CTLA4_uniprot_positions"])
            item["candidate_mask_overlaps"] = _candidate_mask_overlaps(
                manifest["candidate_ids"], mapped, owner, area, geometry_root, apbs_root)
            external_results.append(item)
    else:
        seen = set()
        for candidate in external:
            identity = (candidate["pdb_code"], candidate["fab_chains"], candidate["antigen_chains"])
            if identity in seen: continue
            seen.add(identity)
            cif = _download_cif(candidate, source_structures)
            structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(
                candidate["pdb_code"], str(cif))
            model = next((m for m in structure if int(m.id) == candidate["model"]), None)
            if model is None: continue
            measured = measure_complex(model, candidate["fab_chains"], candidate["antigen_chains"], tools)
            candidate["contact_residue_keys"] = {tuple(x) for x in measured["contact_residue_keys_4_5a"]}
            mapped = _complex_residue_to_uniprot(candidate, model, target_sequence, target_positions)
            overlaps = _candidate_mask_overlaps(
                manifest["candidate_ids"], mapped, owner, area, geometry_root, apbs_root) if mapped else {}
            external_results.append({"pdb_code": candidate["pdb_code"], "instance": candidate["instance"],
                "method": candidate["method"], "resolution_a": candidate["resolution_a"],
                "fab_chains": list(candidate["fab_chains"]), "antigen_chains": list(candidate["antigen_chains"]),
                "coordinate_file_sha256": sha256_file(cif), "mapped_CTLA4_uniprot_positions": sorted(mapped),
                "mapping_rule": "global sequence alignment of complex protein chain to frozen CTLA4 residue sequence",
                "contact_metrics": measured, "candidate_mask_overlaps": overlaps,
                "interpretation": "CTLA4 complex kept external to scale fitting; intersection does not imply binding loss"})

    resolved = [sample for sample in accepted if _resolved_geodesic(sample)]
    radius_values = [sample["geodesic"]["pairwise_surface_distance_quantiles_a"]["0.9"] / 2
                     for sample in resolved]
    radii = []
    for label, quantile in (("small", .25), ("medium", .5), ("large", .75)):
        if radius_values:
            radii.append((label, float(np.quantile(radius_values, quantile))))
    panel_ids = list(manifest["candidate_ids"])
    candidate_masks = {}; candidate_unknown = {}; impact_sources = {}
    for cid in panel_ids:
        geometry_path = geometry_root / cid / "geometry_fusion_masks.npz"
        if not geometry_path.is_file(): continue
        with np.load(geometry_path, allow_pickle=False) as masks:
            candidate_masks[cid] = masks["geometry_chemistry__positive"].copy()
            candidate_unknown[cid] = masks["geometry_chemistry__unknown"].copy()
        impact_sources[cid] = "geometry_chemistry_lower_bound"
        fused_path = apbs_root / cid / "field_fusion_evidence.npz"
        if fused_path.is_file():
            with np.load(fused_path, allow_pickle=False) as masks:
                if "three_mechanism_positive" in masks.files:
                    candidate_masks[cid] = masks["three_mechanism_positive"].copy()
                    candidate_unknown[cid] = masks["three_mechanism_unknown"].copy()
                    impact_sources[cid] = "geometry_chemistry_plus_APBS_delta_phi_lower_bound"
    scenario_results = {}
    for label, radius in [*radii, *((f"legacy_geodesic_{r:g}a", float(r)) for r in (6, 9, 12, 15))]:
        candidate_rows = {}; panel_lower = np.zeros(len(area), bool); panel_unknown = np.zeros(len(area), bool)
        for cid, mask in candidate_masks.items():
            low = _dilate(mask, graph, radius)
            unk = _dilate(candidate_unknown[cid], graph, radius) & ~low
            panel_lower |= low; panel_unknown |= unk
            candidate_rows[cid] = _opportunity(mask, candidate_unknown[cid], owner, area, graph, radius)
        panel_unknown &= ~panel_lower
        denominator = float(area.sum()); low_area = float(area[panel_lower].sum()); unk_area = float(area[panel_unknown].sum())
        scenario_results[label] = {"radius_a": radius, "candidate_results": candidate_rows,
            "panel_lower_area_a2": low_area, "panel_lower_fraction": low_area / denominator,
            "panel_unknown_opportunity_area_a2": unk_area, "panel_unknown_opportunity_fraction": unk_area / denominator,
            "panel_possible_upper_fraction": (low_area + unk_area) / denominator}
    result = {"CTLA4_external_complex_count": len(external_results), "CTLA4_external_complexes": external_results,
        "contact_scale_scan_scenarios": scenario_results, "candidate_impact_sources": impact_sources,
        "contact_scale_source_geodesic_sample_count": len(radius_values),
        "candidate_ids": panel_ids,
        "legacy_official_panel_12A_scan_coverage": 0.8602623530512513,
        "denominator_area_a2": float(area.sum()),
        "interpretation": "irregular verified masks are the physical range; geodesic scenarios expand them only to model assay-scale opportunity"}
    report_name = "ctla4_contact_overlap_and_scan_scenarios.json" if len(panel_ids) == 5 else (
        f"ctla4_contact_overlap_and_scan_scenarios_panel{len(panel_ids)}.json")
    result["report_file"] = report_name
    write_json(out / report_name, result)
    return result


def build(study: Path = DEFAULT_STUDY, *, target_clusters: int = TARGET_CLUSTERS,
          minimum_clusters: int = MIN_CLUSTERS, selection_reference: Path | None = None) -> dict:
    import time
    if minimum_clusters < 1 or target_clusters < minimum_clusters:
        raise ValueError("cluster counts must satisfy 1 <= minimum_clusters <= target_clusters")
    started = time.perf_counter()
    study = Path(study).resolve()
    contact_root = study / "evidence" / "contact_scale"
    source = contact_root / "source"
    source_manifest = json.loads((source / "source_manifest.json").read_text(encoding="utf-8"))
    if not _digest_matches(source / source_manifest["summary_file"], source_manifest["summary_sha256"]):
        raise ValueError("SAbDab summary hash does not match frozen source manifest")
    if not _digest_matches(source / source_manifest["annotations_file"], source_manifest["annotations_sha256"]):
        raise ValueError("SAbDab chain-annotation hash does not match frozen source manifest")
    candidates, _ = _read_frozen_inputs(source)
    if len(candidates) < minimum_clusters:
        raise RuntimeError(f"only {len(candidates)} mapping-valid Fab/protein complexes; need {minimum_clusters}")
    for candidate in candidates:
        candidate["rank_sha256"] = hashlib.sha256(f"{SEED}|{candidate['instance']}".encode()).hexdigest()
    candidates.sort(key=lambda x: (x["rank_sha256"], x["instance"]))
    target_sequence, target_positions = _ctla4_sequence(study)
    cif_dir = source / "structures"; cif_dir.mkdir(parents=True, exist_ok=True)
    tools = PhysicalTools(contact_root / "physical_cache", {"wsl_distro": "Ubuntu-22.04", "use_wsl": True,
                                                            "research_profile": "contact_scale_reference_v1"})
    accepted = []; excluded = []; cluster_sequences = []; seen_instances = set(); current_structure = None; current_pdb = None
    if selection_reference is not None:
        accepted = _replay_frozen_cohort(selection_reference, candidates, cif_dir, target_sequence,
                                         sha256_file(source / "source_manifest.json"), tools)
        candidates_for_scan = []
    else:
        candidates_for_scan = candidates
    for order, candidate in enumerate(candidates_for_scan):
        if candidate["instance"] in seen_instances: continue
        seen_instances.add(candidate["instance"])
        try:
            cif = _download_cif(candidate, cif_dir)
            if current_pdb != candidate["pdb_code"]:
                current_structure = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True).get_structure(
                    candidate["pdb_code"], str(cif))
                current_pdb = candidate["pdb_code"]
            model = next((m for m in current_structure if int(m.id) == candidate["model"]), None)
            if model is None: raise ValueError("requested experimental coordinate model absent")
            sequences = _target_sequences(candidate, model)
            if not sequences: raise ValueError("protein antigen chain sequence is empty")
            candidate["antigen_sequences"] = sequences
            if _homologous(sequences, (target_sequence,)):
                excluded.append({"instance": candidate["instance"], "pdb": candidate["pdb_code"],
                                 "reason": "CTLA4_TARGET_EXTERNAL_SET", "rank": order})
                continue
            match = next((idx for idx, seqs in enumerate(cluster_sequences) if _homologous(sequences, seqs)), None)
            if match is not None:
                excluded.append({"instance": candidate["instance"], "pdb": candidate["pdb_code"],
                                 "reason": f"HOMOLOGOUS_TO_CLUSTER_{match+1:02d}", "rank": order})
                continue
            measured = measure_complex(model, candidate["fab_chains"], candidate["antigen_chains"], tools)
            if not _resolved_geodesic(measured):
                excluded.append({"instance": candidate["instance"], "pdb": candidate["pdb_code"],
                                 "reason": measured["geodesic"]["status"], "rank": order})
                continue
            sample = {k: v for k, v in candidate.items() if k != "row"}
            sample.update(measured)
            sample["coordinate_file_sha256"] = sha256_file(cif)
            sample["protein_antigen_sequences"] = list(sequences)
            sample["selection_order"] = order
            sample["cluster_id"] = f"ANTIGEN_CLUSTER_{len(accepted)+1:02d}"
            accepted.append(sample); cluster_sequences.append(sequences)
            print(f"contact-scale accepted {len(accepted)}/{target_clusters}: {candidate['instance']} {candidate['pdb_code']}", flush=True)
            if len(accepted) >= target_clusters: break
        except (ValueError, KeyError, IndexError) as exc:
            excluded.append({"instance": candidate["instance"], "pdb": candidate["pdb_code"],
                             "reason": f"STRUCTURE_QC:{type(exc).__name__}:{exc}", "rank": order})

    if len(accepted) < minimum_clusters:
        diagnostic_path = _write_incomplete_contact_scale_diagnostic(
            contact_root / "results",
            target_clusters=target_clusters,
            minimum_clusters=minimum_clusters,
            candidates=len(candidates),
            accepted=accepted,
            excluded=excluded,
            source_manifest_sha256=sha256_file(source / "source_manifest.json"),
            tools=tools,
            started=started,
        )
        raise RuntimeError(
            f"only {len(accepted)} nonhomologous complexes have resolved contact-surface geodesics; "
            f"see diagnostic {diagnostic_path}"
        )
    resolved_samples = [sample for sample in accepted if _resolved_geodesic(sample)]
    radius_values = [sample["geodesic"]["pairwise_surface_distance_quantiles_a"]["0.9"] / 2
                     for sample in resolved_samples]
    radii = ({label: float(np.quantile(radius_values, q))
              for label, q in (("small", .25), ("medium", .50), ("large", .75))}
             if radius_values else {})
    methods = {name: sum(sample["method"] == name for sample in accepted) for name in sorted({s["method"] for s in accepted})}
    ctla4_overlap = _apply_ctla4_contact_overlap(study, accepted, source, candidates, tools)
    result = {"profile": "contact_scale_reference_v1", "seed": SEED,
        "source_dataset": "SAbDab2 API 2.1.4; RCSB PDB original mmCIF",
        "source_manifest_sha256": sha256_file(source / "source_manifest.json"),
        "candidate_population_mapping_valid": len(candidates), "target_cluster_count": target_clusters,
        "selected_structure_count": len(accepted), "geodesic_scale_sample_count": len(radius_values),
        "geodesic_reference_status": "RESOLVED" if len(radius_values) >= minimum_clusters else "PARTIAL_GEODESIC_REFERENCE",
        "minimum_required": minimum_clusters,
        "minimum_policy": {"default_minimum": MIN_CLUSTERS, "effective_minimum": minimum_clusters,
                           "override_used": minimum_clusters != MIN_CLUSTERS,
                           "scope": "exploratory structural prior only; not an assay-calibrated sample-size threshold"},
        "homology_rule": {"identity": IDENTITY_CUTOFF, "shorter_sequence_coverage": COVERAGE_CUTOFF,
                          "method": "deterministic global alignment; any antigen-chain pair joins a cluster"},
        "selection_rule": ("replay frozen sample identities/order from selection_reference; revalidate source hash, sequence, chains, model and coordinates"
                           if selection_reference is not None else
                           "SHA256(seed|SAbDab INSTANCE) ascending; first eligible unique antigen sequence clusters; no contact-result ranking"),
        "selection_reference_sha256": sha256_file(selection_reference) if selection_reference is not None else None,
        "experimental_methods": methods, "contact_cutoffs_a": list(CONTACT_CUTOFFS),
        "software_identity": {"numpy_version": np.__version__, "biopython_version": Bio.__version__,
                              "nanoshaper_sha256": tools.identity(tools.nano)},
        "sasa_method": "Bio.PDB Shrake-Rupley; heavy atoms only; 960 points",
        "surface_method": "NanoShaper SES scale 2.0; weighted Euclidean mesh edges; unique 4.5 A contact vertices sampled area-proportionally without replacement, capped at 128; disconnected components use within-component paths",
        "contact_scale_radius_scenarios_a": radii,
        "radius_definition": "quartiles of per-complex half the 90th percentile of finite within-component antigen-contact surface geodesic pairs; below the required resolved sample count the radius reference is marked partial",
        "median_contact_bsa_antigen_a2": float(np.median([x["buried_surface_area_antigen_a2"] for x in accepted])),
        "median_contact_bsa_total_half_a2": float(np.median([x["buried_surface_area_total_half_a2"] for x in accepted])),
        "excluded_attempts": excluded,
        "ctla4_external_validation": {"complex_count": ctla4_overlap["CTLA4_external_complex_count"],
                                      "report": ctla4_overlap["report_file"]},
        "field_costs": {"wall_seconds": float(time.perf_counter() - started), "calls": tools.calls,
                        "activity_seconds": float(sum(row.get("runtime_seconds", 0.) for row in tools.activity)),
                        "source_snapshot_sha256": sha256_file(Path(__file__).resolve())},
        "samples": accepted, "interpretation": "experimentally resolved complex geometry; no binding probability or assay response rate"}
    out = contact_root / "results"; out.mkdir(parents=True, exist_ok=True)
    suffix = "_frozen_cohort_replay" if selection_reference is not None else ""
    reference_name = f"contact_scale_reference{suffix}.json"
    samples_name = f"contact_scale_samples{suffix}.csv"
    exclusions_name = f"contact_scale_exclusions{suffix}.json"
    write_json(out / reference_name, result)
    with (out / samples_name).open("w", newline="", encoding="utf-8") as handle:
        columns = ("cluster_id", "instance", "pdb_code", "method", "resolution_a", "geodesic_status",
                   "connected_pair_fraction", "buried_surface_area_antigen_a2",
                   "buried_surface_area_total_half_a2", "contact_residue_ca_long_axis_a", "contact_residue_ca_short_axis_a", "pairwise_geodesic_p50_a",
                   "pairwise_geodesic_p90_a", "pairwise_geodesic_p95_a", "contact_pairs_4_0a", "contact_pairs_4_5a",
                   "contact_pairs_5_0a", "coordinate_file_sha256")
        writer = csv.DictWriter(handle, fieldnames=columns); writer.writeheader()
        for sample in accepted:
            geodesic = sample.get("geodesic", {})
            q = geodesic.get("pairwise_surface_distance_quantiles_a", {})
            row = {"cluster_id": sample["cluster_id"], "instance": sample["instance"], "pdb_code": sample["pdb_code"],
                   "method": sample["method"], "resolution_a": sample["resolution_a"],
                   "geodesic_status": geodesic.get("status"),
                   "connected_pair_fraction": geodesic.get("connected_pair_fraction"),
                   "buried_surface_area_antigen_a2": sample["buried_surface_area_antigen_a2"],
                   "buried_surface_area_total_half_a2": sample["buried_surface_area_total_half_a2"],
                   "contact_residue_ca_long_axis_a": sample["contact_residue_ca_spans"]["long_axis_a"],
                   "contact_residue_ca_short_axis_a": sample["contact_residue_ca_spans"]["short_axis_a"],
                   "pairwise_geodesic_p50_a": q.get("0.5"), "pairwise_geodesic_p90_a": q.get("0.9"),
                   "pairwise_geodesic_p95_a": q.get("0.95"), "coordinate_file_sha256": sample["coordinate_file_sha256"]}
            for cutoff, key in zip(CONTACT_CUTOFFS, ("contact_pairs_4_0a", "contact_pairs_4_5a", "contact_pairs_5_0a")):
                row[key] = sample["contact_cutoff_counts"][str(cutoff)]["heavy_atom_pairs"]
            writer.writerow(row)
    write_json(out / exclusions_name, {"count": len(excluded), "rows": excluded})
    return result


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=DEFAULT_STUDY)
    parser.add_argument("--target-clusters", type=int, default=TARGET_CLUSTERS)
    parser.add_argument("--minimum-clusters", type=int, default=MIN_CLUSTERS)
    parser.add_argument("--selection-reference", type=Path,
                        help="recompute a frozen cohort from an existing reference without rescanning candidates")
    args = parser.parse_args()
    result = build(args.study, target_clusters=args.target_clusters, minimum_clusters=args.minimum_clusters,
                   selection_reference=args.selection_reference)
    print(json.dumps({k: result[k] for k in ("profile", "selected_structure_count", "contact_scale_radius_scenarios_a",
                                               "median_contact_bsa_antigen_a2", "geodesic_scale_sample_count")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
