"""Build a common-coordinate PDB and contact map from experimental complexes."""
import json
from io import StringIO
from pathlib import Path

import numpy as np
import requests
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web" / "pd1_antibodies.pdb"
CONTACTS_OUT = ROOT / "web" / "pd1_antibody_contacts.json"
SOURCES = [
    ("nivo", "5WT9", "G", "H", "L", "A", "B"),
    # 5GGS contains two complexes: Y contacts C/D; Z contacts A/B.
    ("pembro", "5GGS", "Y", "C", "D", "C", "D"),
    ("clone19", "9HK1", "A", "H", "L", "E", "F"),
]


def protein_residues(chain):
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def matching_ca(reference, moving):
    ref, mov = protein_residues(reference), protein_residues(moving)
    a = "".join(seq1(r.resname, custom_map={"UNK": "X"}) for r in ref)
    b = "".join(seq1(r.resname, custom_map={"UNK": "X"}) for r in mov)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(a, b)[0]
    pairs = [(ref[i]["CA"], mov[j]["CA"]) for i, j in zip(*alignment.indices)
             if i >= 0 and j >= 0 and a[i] == b[j]]
    return [p[0] for p in pairs], [p[1] for p in pairs], {
        mov[j].id: ref[i].id for i, j in zip(*alignment.indices)
        if i >= 0 and j >= 0 and a[i] == b[j]
    }


def heavy_atoms(chains):
    return np.array([atom.coord for chain in chains for residue in chain
                     if residue.id[0] == " " for atom in residue
                     if atom.element not in {"H", "D"}])


def interface(reference_chain, antibody_chains, numbering):
    tree = cKDTree(heavy_atoms(antibody_chains))
    contacts = []
    distances_by_reference = {}
    shortest = float("inf")
    for residue in protein_residues(reference_chain):
        distances, _ = tree.query(np.array([atom.coord for atom in residue
                                            if atom.element not in {"H", "D"}]))
        minimum = float(np.min(distances))
        shortest = min(shortest, minimum)
        if residue.id in numbering and minimum <= 8.0:
            reference_position = numbering[residue.id][1]
            distances_by_reference[reference_position] = min(
                minimum, distances_by_reference.get(reference_position, float("inf"))
            )
            if minimum <= 4.5:
                contacts.append(reference_position)
    if shortest > 5 or len(contacts) < 5:
        raise ValueError(f"Fab is not paired to PD-1: min distance {shortest:.1f} A; {len(contacts)} contacts")
    return sorted(set(contacts)), shortest, {
        str(position): round(distance, 2)
        for position, distance in sorted(distances_by_reference.items())
    }


parser = PDBParser(QUIET=True)
structures = []
for _, pdb, *_ in SOURCES:
    response = requests.get(f"https://files.rcsb.org/download/{pdb}.pdb", timeout=30)
    response.raise_for_status()
    structures.append(parser.get_structure(pdb, StringIO(response.text)))

reference = structures[0][0][SOURCES[0][2]]
numbering = [{r.id: r.id for r in protein_residues(reference)}]
for structure, (_, _, pd1, *_rest) in zip(structures[1:], SOURCES[1:]):
    fixed, moving, mapping = matching_ca(reference, structure[0][pd1])
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    if len(fixed) < 80 or sup.rms > 3:
        raise ValueError(f"Unreliable PD-1 alignment: {len(fixed)} residues, RMSD={sup.rms:.2f} Å")
    sup.apply(structure.get_atoms())
    numbering.append(mapping)
    print(f"{structure.id}: {len(fixed)} matched CA; RMSD {sup.rms:.2f} A")

contact_map = {}
for structure, (name, pdb, pd1, heavy, light, *_), mapping in zip(structures, SOURCES, numbering):
    residues, distance, residue_distances = interface(
        structure[0][pd1], [structure[0][heavy], structure[0][light]], mapping
    )
    contact_map[name] = {"pdb": pdb, "reference_residues": residues,
                         "min_distance_A": round(distance, 2),
                         "distance_by_residue_A": residue_distances}
    print(f"{name}: {len(residues)} PD-1 contact residues; min {distance:.2f} A")

combined = structures[0]
model = combined[0]
for chain in list(model):
    if chain.id not in {"G", "H", "L"}:
        model.detach_child(chain.id)
model["G"].id = "P"
model["H"].id = "A"
model["L"].id = "B"
for structure, (_, _, _, heavy, light, new_heavy, new_light) in zip(structures[1:], SOURCES[1:]):
    for old, new in ((heavy, new_heavy), (light, new_light)):
        chain = structure[0][old].copy()
        chain.id = new
        model.add(chain)

for chain in model:
    for residue in list(chain):
        if residue.id[0] != " ":
            chain.detach_child(residue.id)

io = PDBIO()
io.set_structure(combined)
io.save(str(OUT))
CONTACTS_OUT.write_text(json.dumps(contact_map, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"Wrote {OUT}")
