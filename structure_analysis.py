from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean
from typing import Any

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBIO, MMCIFParser, NeighborSearch, PDBParser, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

from . import config
from .models import ProteinContext, ResidueRecord
from .util import PipelineError

MAX_ASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
# Reference table MAX_ASA below is Tien et al. 2013 (PLoS ONE 8(11):e80635),
# generated with the DSSP radius set.  Any SASA path must use THE SAME radii,
# otherwise rsasa >= 0.25 means different things depending on which solver
# happened to be installed.  Bio.PDB's default ATOMIC_RADII are Bondi 1964
# (O 1.52 vs DSSP 1.40 -- an 8.6% difference on oxygen-rich side chains), so
# the element radii are pinned explicitly everywhere SASA is computed.
TIEN_DSSP_RADII = {"C": 1.80, "N": 1.60, "O": 1.40, "S": 1.85}
GLYCAN_RESNAMES = {
    "NAG",
    "NDG",
    "BMA",
    "MAN",
    "FUC",
    "FUL",
    "GAL",
    "GLC",
    "SIA",
    "NAN",
    "A2G",
    "BGC",
}
POSITIVE_ATOMS = {"ARG": {"NH1", "NH2", "NE"}, "LYS": {"NZ"}, "HIS": {"ND1", "NE2"}}
NEGATIVE_ATOMS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}


def _one_letter(resname: str) -> str:
    return protein_letters_3to1.get(resname.upper(), "M" if resname.upper() == "MSE" else "X")


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _clean_icode(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text in {".", "?"} else text


@dataclass(slots=True)
class ParsedStructure:
    path: Path
    format: str
    structure: Any
    model: Any
    raw_text: str
    label_lookup: dict[tuple[str, int, str], tuple[str, int | None]] = field(default_factory=dict)
    secondary_lookup: dict[tuple[str, int, str], str] = field(default_factory=dict)


@dataclass(slots=True)
class ChainMapping:
    chain_id: str
    residues: list[Any]
    query_to_uniprot: dict[int, int]
    uniprot_to_residue: dict[int, Any]
    identity: float
    observed_ecd_coverage: float
    max_internal_gap: int
    engineered_mismatches: int
    mismatch_positions: list[int] = field(default_factory=list)


def _parse_mmcif_metadata(path: Path) -> tuple[
    dict[tuple[str, int, str], tuple[str, int | None]], dict[tuple[str, int, str], str]
]:
    data = MMCIF2Dict(str(path))
    label_lookup: dict[tuple[str, int, str], tuple[str, int | None]] = {}
    keys = {
        "auth_chain": _as_list(data.get("_atom_site.auth_asym_id")),
        "auth_seq": _as_list(data.get("_atom_site.auth_seq_id")),
        "icode": _as_list(data.get("_atom_site.pdbx_PDB_ins_code")),
        "label_chain": _as_list(data.get("_atom_site.label_asym_id")),
        "label_seq": _as_list(data.get("_atom_site.label_seq_id")),
    }
    count = min((len(values) for values in keys.values() if values), default=0)
    for index in range(count):
        try:
            auth_seq = int(float(keys["auth_seq"][index]))
        except (TypeError, ValueError):
            continue
        try:
            label_seq = int(float(keys["label_seq"][index]))
        except (TypeError, ValueError):
            label_seq = None
        key = (
            str(keys["auth_chain"][index]).strip(),
            auth_seq,
            _clean_icode(keys["icode"][index]),
        )
        label_lookup.setdefault(
            key,
            (str(keys["label_chain"][index]).strip(), label_seq),
        )

    secondary: dict[tuple[str, int, str], str] = {}

    def add_ranges(prefix: str, code: str) -> None:
        chains = _as_list(data.get(f"{prefix}.beg_auth_asym_id"))
        starts = _as_list(data.get(f"{prefix}.beg_auth_seq_id"))
        ends = _as_list(data.get(f"{prefix}.end_auth_seq_id"))
        end_chains = _as_list(data.get(f"{prefix}.end_auth_asym_id"))
        conf_types = _as_list(data.get(f"{prefix}.conf_type_id"))
        n = min(len(chains), len(starts), len(ends), len(end_chains))
        for i in range(n):
            if str(chains[i]) != str(end_chains[i]):
                continue
            if code == "H" and conf_types and not str(conf_types[i]).upper().startswith("HELX"):
                continue
            try:
                start, end = int(float(starts[i])), int(float(ends[i]))
            except (TypeError, ValueError):
                continue
            for number in range(min(start, end), max(start, end) + 1):
                secondary[(str(chains[i]).strip(), number, "")] = code

    add_ranges("_struct_conf", "H")
    add_ranges("_struct_sheet_range", "E")
    return label_lookup, secondary


def _parse_pdb_secondary(text: str) -> dict[tuple[str, int, str], str]:
    secondary: dict[tuple[str, int, str], str] = {}
    for line in text.splitlines():
        try:
            if line.startswith("HELIX"):
                chain, start, end, code = line[19].strip(), int(line[21:25]), int(line[33:37]), "H"
            elif line.startswith("SHEET"):
                chain, start, end, code = line[21].strip(), int(line[22:26]), int(line[33:37]), "E"
            else:
                continue
        except (ValueError, IndexError):
            continue
        for number in range(min(start, end), max(start, end) + 1):
            secondary[(chain, number, "")] = code
    return secondary


def load_structure(path: Path) -> ParsedStructure:
    if not path.exists():
        raise PipelineError(f"Structure file does not exist: {path}")
    raw = path.read_text(encoding="utf-8", errors="replace")
    is_cif = path.suffix.lower() in {".cif", ".mmcif"} or raw.lstrip().startswith("data_")
    if is_cif:
        parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
        structure = parser.get_structure(path.stem, str(path))
        label_lookup, secondary = _parse_mmcif_metadata(path)
        fmt = "mmcif"
    else:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(path.stem, str(path))
        label_lookup, secondary = {}, _parse_pdb_secondary(raw)
        fmt = "pdb"
    models = list(structure.get_models())
    if not models:
        raise PipelineError("Structure contains no coordinate models")
    return ParsedStructure(path, fmt, structure, models[0], raw, label_lookup, secondary)


def standard_residues(chain: Any) -> list[Any]:
    return [
        residue
        for residue in chain
        if residue.id[0] == " " and is_aa(residue, standard=False) and _one_letter(residue.resname) != "X"
    ]


def chain_sequence(chain: Any) -> str:
    """Observed one-letter sequence for a coordinate chain.

    This is deliberately coordinate-derived.  It is suitable for a PDB-only
    ECD run, but must not be presented as a canonical full-length sequence.
    """
    return "".join(_one_letter(residue.resname) for residue in standard_residues(chain))


def _alignment_map(target: str, query: str) -> tuple[dict[int, int], float, int, list[int]]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.3
    aligner.end_insertion_score = 0.0
    aligner.end_deletion_score = 0.0
    alignment = aligner.align(target, query)[0]
    coordinates = np.asarray(alignment.coordinates)
    mapping: dict[int, int] = {}
    matches = 0
    mismatches = 0
    mismatch_positions: list[int] = []
    for column in range(coordinates.shape[1] - 1):
        t0, t1 = int(coordinates[0, column]), int(coordinates[0, column + 1])
        q0, q1 = int(coordinates[1, column]), int(coordinates[1, column + 1])
        if t1 - t0 <= 0 or q1 - q0 <= 0:
            continue
        length = min(t1 - t0, q1 - q0)
        for offset in range(length):
            mapping[q0 + offset] = t0 + offset + 1
            if target[t0 + offset] == query[q0 + offset]:
                matches += 1
            else:
                mismatches += 1
                mismatch_positions.append(t0 + offset + 1)
    identity = matches / max(1, matches + mismatches)
    return mapping, identity, mismatches, mismatch_positions


def map_chain_to_uniprot(
    parsed: ParsedStructure, chain_id: str, protein: ProteinContext
) -> ChainMapping:
    if chain_id not in parsed.model:
        raise PipelineError(f"Chain {chain_id!r} is not present in the structure")
    residues = standard_residues(parsed.model[chain_id])
    query = "".join(_one_letter(residue.resname) for residue in residues)
    if len(query) < 8:
        raise PipelineError(f"Chain {chain_id} contains too few mapped amino-acid residues")
    query_map, identity, mismatches, mismatch_positions = _alignment_map(protein.sequence, query)
    uniprot_to_residue = {query_map[i]: residue for i, residue in enumerate(residues) if i in query_map}
    observed = sorted(
        pos for pos in uniprot_to_residue if protein.ecd_start <= pos <= protein.ecd_end
    )
    coverage = len(observed) / protein.ecd_length
    if observed:
        observed_set = set(observed)
        longest = current = 0
        for position in range(observed[0], observed[-1] + 1):
            if position in observed_set:
                current = 0
            else:
                current += 1
                longest = max(longest, current)
    else:
        longest = protein.ecd_length
    return ChainMapping(
        chain_id=chain_id,
        residues=residues,
        query_to_uniprot=query_map,
        uniprot_to_residue=uniprot_to_residue,
        identity=identity,
        observed_ecd_coverage=coverage,
        max_internal_gap=longest,
        engineered_mismatches=mismatches,
        mismatch_positions=mismatch_positions,
    )


def best_chain_mapping(
    parsed: ParsedStructure,
    protein: ProteinContext,
    preferred_chains: Iterable[str] | None = None,
) -> ChainMapping:
    preferred = list(dict.fromkeys(preferred_chains or []))
    available = [chain.id for chain in parsed.model if standard_residues(chain)]
    order = preferred + [item for item in available if item not in preferred]
    mappings: list[ChainMapping] = []
    for chain_id in order:
        try:
            mappings.append(map_chain_to_uniprot(parsed, chain_id, protein))
        except PipelineError:
            continue
    if not mappings:
        raise PipelineError("No protein chain could be aligned to the UniProt sequence")
    return max(
        mappings,
        key=lambda item: (
            item.observed_ecd_coverage,
            item.identity,
            -item.max_internal_gap,
            -item.engineered_mismatches,
        ),
    )


class _FirstModelSelect(Select):
    def accept_model(self, model: Any) -> bool:
        return model.id == 0


class _ChainOnlySelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_model(self, model: Any) -> bool:
        return model.id == 0

    def accept_chain(self, chain: Any) -> bool:
        return chain.id == self.chain_id


class _MappedEcdSelect(Select):
    def __init__(self, chain_id: str, residue_ids: set[tuple[int, str]]):
        self.chain_id = chain_id
        self.residue_ids = residue_ids

    def accept_model(self, model: Any) -> bool:
        return model.id == 0

    def accept_chain(self, chain: Any) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue: Any) -> bool:
        return residue.id[0] == " " and (int(residue.id[1]), _clean_icode(residue.id[2])) in self.residue_ids

    def accept_atom(self, atom: Any) -> bool:
        return str(getattr(atom, "element", "")).upper() != "H"


def write_pdb(parsed: ParsedStructure, path: Path) -> Path:
    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(str(path), _FirstModelSelect())
    return path


def write_chain_pdb(parsed: ParsedStructure, chain_id: str, path: Path) -> Path:
    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(str(path), _ChainOnlySelect(chain_id))
    return path


def write_mapped_ecd_cif(
    parsed: ParsedStructure,
    mapping: ChainMapping,
    protein: ProteinContext,
    path: Path,
) -> Path:
    from Bio.PDB import MMCIFIO

    residue_ids = {
        (int(residue.id[1]), _clean_icode(residue.id[2]))
        for position, residue in mapping.uniprot_to_residue.items()
        if protein.ecd_start <= position <= protein.ecd_end
    }
    io = MMCIFIO()
    io.set_structure(parsed.structure)
    io.save(str(path), _MappedEcdSelect(mapping.chain_id, residue_ids))
    return path


def write_mapped_ecd_pdb(
    parsed: ParsedStructure,
    mapping: ChainMapping,
    protein: ProteinContext,
    path: Path,
) -> Path:
    residue_ids = {
        (int(residue.id[1]), _clean_icode(residue.id[2]))
        for position, residue in mapping.uniprot_to_residue.items()
        if protein.ecd_start <= position <= protein.ecd_end
    }
    io = PDBIO()
    io.set_structure(parsed.structure)
    io.save(str(path), _MappedEcdSelect(mapping.chain_id, residue_ids))
    return path


def calculate_sasa(parsed: ParsedStructure, workdir: Path) -> tuple[dict[tuple[str, int, str], tuple[float, float]], str]:
    """Production SASA: one deterministic, pinned-radius implementation."""
    calculator = ShrakeRupley(
        probe_radius=1.4,
        n_points=960,
        radii_dict=dict(TIEN_DSSP_RADII),   # pin to the reference-table radius set
    )
    calculator.compute(parsed.model, level="A")
    output: dict[tuple[str, int, str], tuple[float, float]] = {}
    for chain in parsed.model:
        for residue in chain:
            if residue.id[0] == "W":
                continue
            total = sum(float(getattr(atom, "sasa", 0.0)) for atom in residue)
            side = sum(
                float(getattr(atom, "sasa", 0.0))
                for atom in residue
                if atom.name not in BACKBONE_ATOMS
            )
            output[(chain.id, int(residue.id[1]), _clean_icode(residue.id[2]))] = (total, side)
    return output, (
        "Bio.PDB Shrake-Rupley fallback 1.4A/960 points, "
        "element radii pinned to the Tien-2013/DSSP set (C 1.80 N 1.60 O 1.40 S 1.85)"
    )


def _sequon_positions(sequence: str) -> set[int]:
    positions: set[int] = set()
    for index in range(len(sequence) - 2):
        if sequence[index] == "N" and sequence[index + 1] != "P" and sequence[index + 2] in "ST":
            positions.update({index + 1, index + 2, index + 3})
    return positions


def _feature_positions(protein: ProteinContext, patterns: Iterable[str]) -> set[int]:
    regexes = [re.compile(pattern, re.I) for pattern in patterns]
    positions: set[int] = set()
    for feature in protein.features:
        text = f"{feature.type} {feature.description}"
        if any(regex.search(text) for regex in regexes):
            if feature.type.lower() in {"disulfide bond", "cross-link"}:
                # UniProt encodes the two bonded residues as start/end, not a protected interval.
                positions.update({feature.start, feature.end})
            else:
                positions.update(range(feature.start, feature.end + 1))
    return positions


def _domain_for_position(protein: ProteinContext, position: int) -> str:
    domains = [
        item
        for item in protein.features
        if item.type.lower() in {"domain", "repeat"} and item.contains(position)
    ]
    if not domains:
        return "ECD"
    item = min(domains, key=lambda value: value.end - value.start)
    return item.description or f"{item.type} {item.start}-{item.end}"


def _mean_bfactor(residue: Any) -> float | None:
    values = [float(atom.bfactor) for atom in residue if atom.bfactor is not None]
    return mean(values) if values else None


def _af3_atom_plddt(parsed: ParsedStructure, confidence: dict[str, Any] | None) -> dict[int, float]:
    if not confidence:
        return {}
    values = confidence.get("atom_plddts") or confidence.get("atom_plddt")
    if not isinstance(values, list):
        return {}
    atoms = list(parsed.model.get_atoms())
    if len(values) != len(atoms):
        return {}
    grouped: dict[int, list[float]] = {}
    for atom, value in zip(atoms, values):
        grouped.setdefault(id(atom.get_parent()), []).append(float(value))
    return {key: mean(items) for key, items in grouped.items()}


def _charged_atoms(residue: Any, table: dict[str, set[str]]) -> list[Any]:
    wanted = table.get(residue.resname.upper(), set())
    return [atom for atom in residue if atom.name in wanted]


def _distance(a: tuple[float, float, float] | np.ndarray, b: tuple[float, float, float] | np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def _sidechain_direction(residue: Any) -> tuple[float, float, float] | None:
    """Return a local outward proxy without requiring a molecular surface mesh."""
    if not residue.has_id("CA"):
        return None
    sidechain_atoms = [
        np.asarray(atom.coord, dtype=float)
        for atom in residue
        if atom.name not in BACKBONE_ATOMS and atom.element != "H"
    ]
    if not sidechain_atoms:
        return None
    vector = np.mean(sidechain_atoms, axis=0) - np.asarray(residue["CA"].coord, dtype=float)
    length = float(np.linalg.norm(vector))
    if length < 1e-6:
        return None
    return tuple(float(value) for value in vector / length)


def _fibonacci_sphere(count: int) -> np.ndarray:
    """Deterministic unit-sphere points used only for local surface normals."""
    indices = np.arange(count, dtype=float)
    y = 1.0 - 2.0 * (indices + 0.5) / count
    radius = np.sqrt(np.maximum(0.0, 1.0 - y * y))
    theta = np.pi * (3.0 - np.sqrt(5.0)) * indices
    return np.column_stack((np.cos(theta) * radius, y, np.sin(theta) * radius))


def _atomic_radius(atom: Any) -> float:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    return float(TIEN_DSSP_RADII.get(element, 1.70))


def _solvent_surface_normal(residue: Any, atom_search: NeighborSearch) -> tuple[float, float, float] | None:
    """Estimate an outward normal from genuinely unoccluded atomic points.

    The same 1.4 A probe and pinned element radii as production SASA are used.
    This local construction is rotation invariant and does not assume that the
    domain is spherical.  A small point set is enough for direction; SASA
    itself continues to use the audited 960-point Shrake-Rupley calculation.
    """
    sphere = _fibonacci_sphere(config.COMPACT_SURFACE_NORMAL_POINTS)
    accumulated = np.zeros(3, dtype=float)
    exposed = 0
    for atom in residue:
        if str(getattr(atom, "element", "")).upper() == "H":
            continue
        atom_radius = _atomic_radius(atom) + 1.4
        center = np.asarray(atom.coord, dtype=float)
        for direction in sphere:
            point = center + atom_radius * direction
            occluded = False
            for other in atom_search.search(point, 3.7, level="A"):
                if other is atom or str(getattr(other, "element", "")).upper() == "H":
                    continue
                other_radius = _atomic_radius(other) + 1.4
                if float(np.linalg.norm(point - np.asarray(other.coord, dtype=float))) < other_radius - 1e-4:
                    occluded = True
                    break
            if not occluded:
                accumulated += direction
                exposed += 1
    length = float(np.linalg.norm(accumulated))
    if exposed == 0 or length < 1e-6:
        return None
    return tuple(float(value) for value in accumulated / length)


def analyze_residues(
    parsed: ParsedStructure,
    mapping: ChainMapping,
    protein: ProteinContext,
    workdir: Path,
    predicted: bool = False,
    confidence: dict[str, Any] | None = None,
    rsasa_threshold: float = config.RSASA_THRESHOLD,
    sidechain_sasa_threshold: float = config.SIDECHAIN_SASA_THRESHOLD_A2,
    plddt_threshold: float = 70.0,
) -> tuple[list[ResidueRecord], dict[str, Any]]:
    sasa, sasa_method = calculate_sasa(parsed, workdir)
    af3_plddt = _af3_atom_plddt(parsed, confidence)
    target_residue_ids = {id(residue) for residue in mapping.uniprot_to_residue.values()}

    partner_atoms: list[Any] = []
    glycan_atoms: list[Any] = []
    for chain in parsed.model:
        for residue in chain:
            if residue.id[0] == "W":
                continue
            if residue.resname.upper() in GLYCAN_RESNAMES:
                glycan_atoms.extend(list(residue.get_atoms()))
            if chain.id != mapping.chain_id or residue.id[0] != " ":
                partner_atoms.extend(list(residue.get_atoms()))
    partner_search = NeighborSearch(partner_atoms) if partner_atoms else None
    glycan_search = NeighborSearch(glycan_atoms) if glycan_atoms else None

    ecd_residues = [
        (position, residue)
        for position, residue in mapping.uniprot_to_residue.items()
        if protein.ecd_start <= position <= protein.ecd_end
    ]
    ecd_residues.sort(key=lambda pair: pair[0])
    ca_positions = [residue["CA"].coord for _, residue in ecd_residues if residue.has_id("CA")]
    ecd_heavy_atoms = [
        atom for _, residue in ecd_residues for atom in residue
        if str(getattr(atom, "element", "")).upper() != "H"
    ]
    ecd_atom_search = NeighborSearch(ecd_heavy_atoms)

    disulfide_positions = _feature_positions(protein, [r"disulfide bond"])
    glyco_positions = _feature_positions(protein, [r"glycosylation"])
    ptm_positions = _feature_positions(protein, [r"modified residue", r"lipidation", r"cross-link"])
    sequon_positions = _sequon_positions(protein.sequence)

    # Include coordinate-observed disulfides even when UniProt lacks annotation.
    cysteines = [(position, residue) for position, residue in ecd_residues if residue.resname == "CYS" and residue.has_id("SG")]
    for index, (pos_a, res_a) in enumerate(cysteines):
        for pos_b, res_b in cysteines[index + 1 :]:
            if _distance(res_a["SG"].coord, res_b["SG"].coord) <= 2.35:
                disulfide_positions.update({pos_a, pos_b})

    records: list[ResidueRecord] = []
    position_to_record: dict[int, ResidueRecord] = {}
    for query_index, residue in enumerate(mapping.residues):
        uniprot_pos = mapping.query_to_uniprot.get(query_index)
        if uniprot_pos is None or not (protein.ecd_start <= uniprot_pos <= protein.ecd_end):
            continue
        aa = _one_letter(residue.resname)
        auth_seq = int(residue.id[1])
        icode = _clean_icode(residue.id[2])
        label_chain, label_seq = parsed.label_lookup.get(
            (mapping.chain_id, auth_seq, icode), (mapping.chain_id, None)
        )
        total_sasa, side_sasa = sasa.get((mapping.chain_id, auth_seq, icode), (0.0, 0.0))
        ca = tuple(float(value) for value in residue["CA"].coord) if residue.has_id("CA") else None
        contacts = 0
        if ca is not None:
            contacts = sum(1 for other in ca_positions if 0.1 < _distance(ca, other) <= 10.0)

        interface_partners: set[str] = set()
        interface = False
        glycan_shielded = False
        for atom in residue:
            if partner_search:
                for partner in partner_search.search(atom.coord, 5.0, level="R"):
                    if id(partner) in target_residue_ids or partner.id[0] == "W":
                        continue
                    partner_chain = partner.get_parent().id
                    interface = True
                    interface_partners.add(
                        f"{partner_chain}:{partner.resname}{partner.id[1]}{_clean_icode(partner.id[2])}"
                    )
            if glycan_search and glycan_search.search(atom.coord, 8.0, level="A"):
                glycan_shielded = True

        positive = _charged_atoms(residue, POSITIVE_ATOMS)
        negative = _charged_atoms(residue, NEGATIVE_ATOMS)
        salt_bridge = False
        if positive or negative:
            for _, other in ecd_residues:
                if other is residue:
                    continue
                opposite = _charged_atoms(other, NEGATIVE_ATOMS if positive else POSITIVE_ATOMS)
                if any(_distance(a.coord, b.coord) <= 4.0 for a in (positive or negative) for b in opposite):
                    salt_bridge = True
                    break

        plddt = af3_plddt.get(id(residue)) if predicted else None
        if predicted and plddt is None:
            plddt = _mean_bfactor(residue)
        domain = _domain_for_position(protein, uniprot_pos)
        rsasa = total_sasa / MAX_ASA.get(aa, 200.0)
        # Exposure and residue size are separate: Gly has no side-chain heavy
        # atoms, and an absolute side-chain cutoff biases smaller residues.
        surface = rsasa >= rsasa_threshold
        hard: list[str] = []
        risks: list[str] = []
        review: list[str] = []
        if rsasa < rsasa_threshold:
            hard.append("insufficient-relative-SASA")
        if uniprot_pos in disulfide_positions:
            hard.append("disulfide-cysteine")
        if uniprot_pos in glyco_positions:
            hard.append("annotated-glycosylation")
        if uniprot_pos in ptm_positions:
            hard.append("annotated-PTM")
        if predicted and (plddt is None or plddt < plddt_threshold):
            review.append("low-pLDDT")
        if aa in "GP":
            risks.append("backbone-conformation-sensitive")
        secondary = parsed.secondary_lookup.get(
            (mapping.chain_id, auth_seq, icode),
            parsed.secondary_lookup.get((mapping.chain_id, auth_seq, ""), "U"),
        )
        if secondary == "E":
            risks.append("beta-strand")
        if salt_bridge:
            risks.append("salt-bridge")
        if interface:
            risks.append("native-interface")
        if glycan_shielded:
            risks.append("glycan-shielded")
        if uniprot_pos in sequon_positions:
            risks.append("N-glycosylation-sequon")

        record = ResidueRecord(
            uniprot_pos=uniprot_pos,
            aa=aa,
            chain_id=mapping.chain_id,
            auth_seq_id=auth_seq,
            insertion_code=icode,
            label_chain_id=label_chain,
            label_seq_id=label_seq,
            residue_name=residue.resname,
            ca=ca,
            sidechain_direction=_sidechain_direction(residue),
            surface_normal=_solvent_surface_normal(residue, ecd_atom_search),
            sasa=round(total_sasa, 3),
            sidechain_sasa=round(side_sasa, 3),
            rsasa=round(rsasa, 4),
            contacts=contacts,
            plddt=round(plddt, 2) if plddt is not None else None,
            secondary_structure=secondary,
            domain=domain,
            interface=interface,
            interface_partners=sorted(interface_partners),
            glycan_shielded=glycan_shielded,
            glycosylation_site=uniprot_pos in glyco_positions,
            glycosylation_sequon=uniprot_pos in sequon_positions,
            disulfide=uniprot_pos in disulfide_positions,
            ptm=uniprot_pos in ptm_positions,
            salt_bridge=salt_bridge,
            surface_accessible=surface,
            folding_screen_pass=surface and not hard and not review,
            mutation_qualification=("EXCLUDED" if hard else "REVIEW_REQUIRED" if review else "ELIGIBLE"),
            review_reasons=review,
            hard_exclusions=hard,
            risk_flags=risks,
            reasons=(
                [f"rSASA={rsasa:.3f}", f"side-chain-SASA={side_sasa:.1f}A2"]
                + ([f"pLDDT={plddt:.1f}"] if plddt is not None else [])
                + risks
                + hard
            ),
            token_index=query_index,
        )
        records.append(record)
        position_to_record[uniprot_pos] = record

    # Flag contacts that cross annotated domain boundaries.
    for index, record in enumerate(records):
        if record.ca is None:
            continue
        for other in records[index + 1 :]:
            if other.ca is None or other.domain == record.domain:
                continue
            if _distance(record.ca, other.ca) <= 8.0:
                record.domain_interface = True
                other.domain_interface = True
                if "domain-interface" not in record.risk_flags:
                    record.risk_flags.append("domain-interface")
                if "domain-interface" not in other.risk_flags:
                    other.risk_flags.append("domain-interface")

    flexible_ranges = flag_flexible_terminal_segments(records)
    metadata = {
        "sasa_method": sasa_method,
        "rsasa_threshold": rsasa_threshold,
        "sidechain_sasa_threshold_A2": sidechain_sasa_threshold,
        "sidechain_sasa_role": "descriptive_only_not_surface_gate",
        "plddt_threshold": plddt_threshold if predicted else None,
        "interface_definition": "any target heavy atom within 5.0 A of a non-target-chain or hetero residue",
        "glycan_shield_definition": "any target heavy atom within 8.0 A of a recognized glycan residue",
        "flexible_terminal_ranges": [
            {"start": start, "end": end} for start, end in flexible_ranges
        ],
        "mapped_ecd_residues": len(records),
    }
    return records, metadata


def flag_flexible_terminal_segments(records: list[ResidueRecord]) -> list[tuple[int, int]]:
    """Mark terminal runs with no observed secondary structure.

    Long unstructured N/C-terminal tails are the least trustworthy part of
    any structure — a crystal truncation artefact, one frame of an ensemble,
    or a low-pLDDT model segment — and they are poor mapping instruments: the
    tail conformation adapts, so a retained antibody binding cannot be read
    as "epitope not here".  Flagged residues are excluded from patch design
    (config.EXCLUDE_FLEXIBLE_TERMINAL) and reported separately as
    linear-epitope candidates.  The walk stops at the first helix/strand
    residue, which keeps the rule conservative.
    """
    from . import config as _config

    ordered = sorted(
        (record for record in records if record.ca is not None),
        key=lambda record: record.uniprot_pos,
    )
    if not ordered:
        return []
    # Guard: the terminal walk only means something when the structure actually
    # carries SSE annotation.  A PDB without HELIX/SHEET records (or an mmCIF
    # without parsed sec-struct) makes every residue look "C" -- flagging then
    # would mark the whole ECD as tail.  Require a minimal annotated fraction.
    structured = sum(1 for record in ordered if record.secondary_structure in {"H", "E"})
    if structured / len(ordered) < _config.FLEXIBLE_TERMINAL_MIN_SSE_FRACTION:
        return []
    ranges: list[tuple[int, int]] = []
    for from_start in (True, False):
        run: list[ResidueRecord] = []
        for record in (ordered if from_start else reversed(ordered)):
            if record.secondary_structure != "C":
                break
            run.append(record)
        if len(run) >= _config.FLEXIBLE_TERMINAL_MIN_RESIDUES:
            positions = sorted(record.uniprot_pos for record in run)
            ranges.append((positions[0], positions[-1]))
            for record in run:
                record.flexible_terminal = True
                if "flexible-terminal" not in record.risk_flags:
                    record.risk_flags.append("flexible-terminal")
    return ranges


def median_region_pae(region_positions: list[int], records: list[ResidueRecord], confidence: dict[str, Any] | None) -> float | None:
    if not confidence:
        return None
    pae = confidence.get("pae") or confidence.get("predicted_aligned_error")
    if not isinstance(pae, list):
        return None
    matrix = np.asarray(pae, dtype=float)
    by_position = {record.uniprot_pos: record for record in records}
    indices = [
        by_position[position].token_index
        for position in region_positions
        if position in by_position and by_position[position].token_index is not None
    ]
    indices = [int(index) for index in indices if 0 <= int(index) < matrix.shape[0]]
    if len(indices) < 2:
        return None
    values = matrix[np.ix_(indices, indices)]
    return round(float(np.median(values)), 3)
