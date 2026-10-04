"""Collect PD-1 complex structures, align interfaces and export a bilingual report."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import requests
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, PDBParser, Superimposer
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "pd1_antibody_epitope_atlas"
STRUCTURES = OUT / "structures"
FIGURES = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
STRUCTURES.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

ANTIBODIES = [
    {"name": "Nivolumab", "zh": "纳武利尤单抗", "pdb": "5WT9", "source": "Bristol Myers Squibb (BMS)", "origin": "Fully human IgG4; BMS", "group": "N-terminal / loop face"},
    {"name": "Pembrolizumab", "zh": "帕博利珠单抗", "pdb": "5GGS", "source": "Merck Sharp & Dohme (MSD)", "origin": "Humanized IgG4; MSD", "group": "C′D / front β-sheet"},
    {"name": "Toripalimab", "zh": "特瑞普利单抗", "pdb": "6JBT", "source": "Shanghai Junshi Biosciences", "origin": "Humanized IgG4; Junshi", "group": "FG loop"},
    {"name": "Tislelizumab", "zh": "替雷利珠单抗", "pdb": "7CGW", "source": "BeiGene", "origin": "Humanized IgG4; BeiGene", "group": "CC′ / front β-sheet"},
    {"name": "Camrelizumab", "zh": "卡瑞利珠单抗", "pdb": "7CU5", "source": "Jiangsu Hengrui Medicine", "origin": "Humanized IgG4; Hengrui", "group": "BC / C′D / FG; N58 glycan"},
    {"name": "Cemiplimab", "zh": "西米普利单抗", "pdb": "8GY5", "source": "Regeneron / Sanofi", "origin": "Fully human IgG4; Regeneron", "group": "BC / C′D / FG; N58 glycan"},
    {"name": "Dostarlimab", "zh": "多塔利单抗", "pdb": "7WSL", "source": "AnaptysBio; licensed to GSK", "origin": "Humanized IgG4; AnaptysBio / GSK", "group": "BC / C′D / FG"},
    {"name": "Sasanlimab", "zh": "沙纳利单抗", "pdb": "6XKR", "source": "Pfizer", "origin": "Engineered human IgG4; Pfizer", "group": "CC′ / FG"},
    {"name": "MW11-h317", "zh": "MW11-h317", "pdb": "6JJP", "source": "Wang et al.; research antibody", "origin": "Research monoclonal antibody", "group": "BC / C′D / FG; glycan-adjacent"},
    {"name": "mAb059c", "zh": "mAb059c", "pdb": "6K0Y", "source": "Liu et al.; research antibody", "origin": "Research monoclonal antibody", "group": "BC / C′D / FG"},
    {"name": "Serplulimab", "zh": "斯鲁利单抗", "pdb": "7E9B", "source": "Shanghai Henlius Biotech", "origin": "Humanized anti-PD-1; Henlius", "group": "BC / C′D / FG"},
    {"name": "GY-5", "zh": "GY-5", "pdb": "6J15", "source": "Chen et al.; research antibody", "origin": "Research monoclonal antibody", "group": "N-terminal / multi-loop; FG hotspot"},
    {"name": "GY-14", "zh": "GY-14", "pdb": "6J14", "source": "Chen et al.; research antibody", "origin": "Research monoclonal antibody", "group": "FG hotspot"},
    {"name": "NBO1a", "zh": "NBO1a", "pdb": "6HIG", "source": "Fenwick et al.; research antibody", "origin": "Fully human research antibody", "group": "N-terminal / FG"},
    {"name": "D12", "zh": "D12", "pdb": "8AS0", "source": "AMG synthetic phage-display library", "origin": "Fully human research antibody", "group": "N-terminal conformation / multi-loop"},
    {"name": "Clone 19", "zh": "Clone 19", "pdb": "9HK1", "source": "Antibody agonism research study", "origin": "Research agonist antibody", "group": "Membrane-proximal backside"},
    {"name": "Prolgolimab", "zh": "普乐利单抗", "pdb": None, "source": "BIOCAD", "origin": "Humanized anti-PD-1; BIOCAD", "group": "FG / loop face", "evidence": "Modeled comparison; no PD-1 co-crystal located", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11452378/"},
    {"name": "UDIZ-007 / UDIZ-008", "zh": "UDIZ-007 / UDIZ-008", "pdb": None, "source": "Synthetic phage-display library", "origin": "Same V regions; distinct Fc formats", "group": "C′D / D-E strands; opposite PD-L1 face", "evidence": "Cross-linking MS epitope mapping; no co-crystal located", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12667665/"},
]

PARSER = PDBParser(QUIET=True)
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "PD1-epitope-atlas/1.0"})
CANONICAL_BY_REF_ID = {}


def residues(chain):
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def sequence(chain):
    return "".join(seq1(r.resname, custom_map={"UNK": "X"}) for r in residues(chain))


def atoms(chain):
    return [a for r in residues(chain) for a in r if a.element not in {"H", "D"}]


def pairwise_chain_distance(a, b):
    aa, bb = atoms(a), atoms(b)
    if not aa or not bb:
        return float("inf")
    tree = cKDTree(np.array([x.coord for x in aa]))
    d, _ = tree.query(np.array([x.coord for x in bb]), k=1)
    return float(np.min(d))


def sequence_pairs(reference, mobile):
    ref, mob = residues(reference), residues(mobile)
    a, b = sequence(reference), sequence(mobile)
    aligner = PairwiseAligner()
    aligner.mode, aligner.match_score, aligner.mismatch_score = "global", 2, -2
    aligner.open_gap_score, aligner.extend_gap_score = -8, -0.5
    alignment = aligner.align(a, b)[0]
    pairs = [(i, j) for i, j in zip(*alignment.indices)
             if i >= 0 and j >= 0 and a[i] == b[j]]
    return pairs, ref, mob


def extract_complex(ab: dict, reference_chain):
    pdb_id = ab["pdb"]
    target = STRUCTURES / f"{pdb_id}.pdb"
    if not target.exists():
        response = SESSION.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=45)
        response.raise_for_status()
        target.write_bytes(response.content)
    structure = PARSER.get_structure(pdb_id, str(target))
    model = structure[0]
    ref_seq = sequence(reference_chain)
    antigen_candidates = []
    for chain in model:
        mobile_seq = sequence(chain)
        if len(mobile_seq) < 60:
            continue
        pairs, _, _ = sequence_pairs(reference_chain, chain)
        identity = len(pairs) / max(1, min(len(ref_seq), len(mobile_seq)))
        if len(pairs) >= 75 and identity >= 0.65:
            antigen_candidates.append((len(pairs), chain))
    if not antigen_candidates:
        raise ValueError(f"No PD-1 chain found in {pdb_id}")
    _, antigen = max(antigen_candidates, key=lambda x: x[0])
    contacts = []
    for chain in model:
        if chain.id != antigen.id and len(residues(chain)) >= 60:
            d = pairwise_chain_distance(antigen, chain)
            if d < 6:
                contacts.append((d, chain))
    contacts.sort(key=lambda x: x[0])
    # Each entry is a single Fab/scFv complex. Keep the two closest contacting chains.
    binder_chains = [chain for _, chain in contacts[:2]]
    # 8AS0 contains several crystallographic copies. Chain A is paired with
    # heavy B and light C; light chain P belongs to another copy and is a
    # crystal-packing neighbor that also approaches A.
    if pdb_id == "8AS0" and antigen.id == "A":
        binder_chains = [model["B"], model["C"]]
    if not binder_chains:
        raise ValueError(f"No antibody chains contacting {pdb_id} chain {antigen.id}")
    ref_atoms, mobile_atoms = [], []
    pair_idx, ref_res, mob_res = sequence_pairs(reference_chain, antigen)
    for i, j in pair_idx:
        if "CA" in ref_res[i] and "CA" in mob_res[j]:
            ref_atoms.append(ref_res[i]["CA"])
            mobile_atoms.append(mob_res[j]["CA"])
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mobile_atoms)
    if len(ref_atoms) < 75 or sup.rms > 4:
        raise ValueError(f"Poor PD-1 fit {pdb_id}: {len(ref_atoms)} CA, RMSD {sup.rms:.2f}")
    sup.apply(structure.get_atoms())

    # Map PDB PD-1 residues to the canonical mature ectodomain numbering represented
    # by the human reference chain; keep canonical UniProt precursor positions.
    contacts_by_ref = {}
    binder_atoms = [a for c in binder_chains for a in atoms(c)]
    tree = cKDTree(np.array([a.coord for a in binder_atoms]))
    pair_idx, ref_res, mob_res = sequence_pairs(reference_chain, antigen)
    mobile_to_reference = {
        mob_res[j].id: CANONICAL_BY_REF_ID[ref_res[i].id]
        for i, j in pair_idx if ref_res[i].id in CANONICAL_BY_REF_ID
    }
    for residue in residues(antigen):
        ref_pos = mobile_to_reference.get(residue.id)
        if ref_pos is None:
            continue
        r_atoms = [a for a in residue if a.element not in {"H", "D"}]
        dist, _ = tree.query(np.array([a.coord for a in r_atoms]))
        minimum = float(np.min(dist))
        if minimum <= 8:
            contacts_by_ref[str(ref_pos)] = round(minimum, 2)
    if sum(d <= 4.5 for d in contacts_by_ref.values()) < 5:
        raise ValueError(f"Weak/inconsistent interface in {pdb_id}")

    # Save only PD-1 and directly paired binder chains for PyMOL.
    keep_ids = {antigen.id, *(c.id for c in binder_chains)}
    for chain in list(model):
        if chain.id not in keep_ids:
            model.detach_child(chain.id)
    for chain in model:
        for residue in list(chain):
            if residue.id[0] != " ":
                chain.detach_child(residue.id)
    # Standardize PD-1 numbering to the UniProt precursor numbering represented
    # by 5WT9 chain G; the same numbering is then usable for contact coloring.
    renumber = []
    for index, residue in enumerate(list(antigen)):
        ref_pos = mobile_to_reference.get(residue.id)
        if ref_pos is not None:
            renumber.append((residue, ref_pos, index))
    mapped_objects = {id(item[0]) for item in renumber}
    for residue in list(antigen):
        if id(residue) not in mapped_objects:
            antigen.detach_child(residue.id)
    for residue, _, index in renumber:
        residue.id = (" ", -1000 - index, " ")
    for residue, ref_pos, _ in renumber:
        residue.id = (" ", ref_pos, " ")
    io = PDBIO()
    io.set_structure(structure)
    aligned_path = STRUCTURES / f"{pdb_id}_aligned.pdb"
    io.save(str(aligned_path))
    return {
        "pdb": pdb_id,
        "antigen_chain": antigen.id,
        "binder_chains": [c.id for c in binder_chains],
        "rmsd_A": round(float(sup.rms), 2),
        "aligned_ca": len(ref_atoms),
        "contacts_4_5A": [int(k) for k, d in contacts_by_ref.items() if d <= 4.5],
        "contacts_8A": [int(k) for k in contacts_by_ref],
        "distance_by_residue_A": {int(k): v for k, v in contacts_by_ref.items()},
        "min_interface_distance_A": min(contacts_by_ref.values()),
        "aligned_path": str(aligned_path),
    }


def query_metadata(pdb_id):
    try:
        entry = SESSION.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=30).json()
        title = entry.get("struct", {}).get("title", "")
        citation = entry.get("rcsb_primary_citation", {})
        return {"title": title, "citation_title": citation.get("title", ""),
                "citation_year": citation.get("year", ""), "doi": citation.get("pdbx_database_id_DOI", ""),
                "resolution_A": entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]}
    except requests.RequestException:
        cache = OUT / "interface_data.json"
        if cache.exists():
            for item in json.loads(cache.read_text(encoding="utf-8")).get("antibodies", []):
                if item.get("pdb") == pdb_id:
                    return {k: item.get(k, "") for k in
                            ("title", "citation_title", "citation_year", "doi", "resolution_A")}
        return {"title": "", "citation_title": "", "citation_year": "", "doi": "", "resolution_A": None}


def main():
    global CANONICAL_BY_REF_ID
    canonical = SESSION.get("https://rest.uniprot.org/uniprotkb/Q15116.fasta", timeout=30).text
    canonical_seq = "".join(canonical.splitlines()[1:])
    # Build an artificial full precursor reference with the deposited extracellular
    # chain mapped back to UniProt sequence coordinates.
    pdb_ref = STRUCTURES / "5WT9.pdb"
    if not pdb_ref.exists():
        response = SESSION.get("https://files.rcsb.org/download/5WT9.pdb", timeout=45)
        response.raise_for_status()
        pdb_ref.write_bytes(response.content)
    ref_structure = PARSER.get_structure("reference", str(pdb_ref))[0]
    reference_chain = ref_structure["G"]
    aligner = PairwiseAligner()
    aligner.mode, aligner.match_score, aligner.mismatch_score = "local", 2, -2
    aligner.open_gap_score, aligner.extend_gap_score = -8, -0.5
    ref_seq = sequence(reference_chain)
    canonical_alignment = aligner.align(canonical_seq, ref_seq)[0]
    ref_residues = residues(reference_chain)
    for uni_idx, ref_idx in zip(*canonical_alignment.indices):
        if uni_idx >= 0 and ref_idx >= 0 and canonical_seq[uni_idx] == ref_seq[ref_idx]:
            CANONICAL_BY_REF_ID[ref_residues[ref_idx].id] = uni_idx + 1
    report = []
    for antibody in ANTIBODIES:
        row = dict(antibody)
        if antibody["pdb"]:
            row.update(query_metadata(antibody["pdb"]))
            if antibody["pdb"] == "8GY5" and not row.get("doi"):
                row["doi"] = "10.3390/biomedicines10123154"
                row["citation_title"] = "The High-Resolution Structure Reveals Remarkable Similarity in PD-1 Binding of Cemiplimab and Dostarlimab"
                row["citation_year"] = 2022
            row["structure"] = extract_complex(antibody, reference_chain)
        else:
            row["structure"] = None
        report.append(row)
        if row["structure"]:
            print(antibody["name"], antibody["pdb"], row["structure"]["antigen_chain"],
                  ",".join(row["structure"]["binder_chains"]), row["structure"]["rmsd_A"],
                  len(row["structure"]["contacts_4_5A"]))
    payload = {"reference": {"pdb": "5WT9", "chain": "G", "uniprot": "Q15116",
                             "sequence_length": len(canonical_seq)}, "antibodies": report}
    (OUT / "interface_data.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {OUT / 'interface_data.json'}")


if __name__ == "__main__":
    main()
