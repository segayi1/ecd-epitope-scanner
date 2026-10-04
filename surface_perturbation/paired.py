"""Explicit FoldX PdbList identities; no positional pairing of PDB files."""
from __future__ import annotations

import re
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.Data.PDBData import protein_letters_3to1
from ..util import sha256_file


def model_identity(name):
    match=re.search(r"_(\d+)_(\d+)\.pdb$",name)
    if match: return tuple(map(int,match.groups()))
    # numberOfRuns=1 omits the run suffix. The group is still explicit;
    # WT/Mut residue identities are checked against that exact batch group.
    match=re.search(r"_(\d+)\.pdb$",name)
    return (int(match.group(1)),0) if match else None


def atoms(path):
    model = PDBParser(QUIET=True).get_structure("pair", str(path))[0]
    output = {}
    for atom in model.get_atoms():
        residue = atom.get_parent()
        if residue.id[0] != " " or str(atom.element).upper() == "H":
            continue
        key = (residue.get_parent().id, residue.id[1], residue.id[2].strip(), atom.name.strip())
        if key in output:
            raise ValueError("ambiguous atom identity")
        output[key] = (np.asarray(atom.coord, dtype=float), residue.resname, str(atom.element).upper())
    return output


def scan_pairs(root: Path):
    """Bind actual PdbList rows to batch mutation records and structural WT."""
    pairs, failures = [], []
    for listing in sorted(root.glob("engine_runs/*/PdbList_*.fxout")):
        batch = listing.parent
        from ..stability import _parse_dif_rows
        ddg_by_identity = {}
        ddg_files = {}
        for score_file in batch.glob("Dif_*.fxout"):
            for name, ddg in _parse_dif_rows(score_file):
                identity = model_identity(name)
                if identity:
                    if identity in ddg_by_identity: raise ValueError("duplicate paired DDG identity")
                    ddg_by_identity[identity] = ddg
                    ddg_files[identity]=score_file
        groups = (batch / "individual_list.txt").read_text(encoding="ascii").splitlines()
        mapping = {}
        listed = re.findall(r"[^\s,;]+\.pdb", listing.read_text(encoding="utf-8", errors="replace"))
        # FoldX versions emit either two-column rows or separate WT/Mut lists.
        # Join by the explicit group/run identity present in BOTH listed names.
        by_identity = {}
        for name in listed:
            identity = model_identity(name)
            if not identity:
                failures.append({"batch": str(batch), "reason": "unmapped PdbList name"}); continue
            by_identity.setdefault(identity, []).append(name)
        for names in by_identity.values():
            if len(names) != 2:
                failures.append({"batch": str(batch), "reason": "incomplete/duplicate PdbList identity"}); continue
            wt_names = [n for n in names if Path(n).name.startswith("WT_")]
            if len(wt_names) != 1:
                failures.append({"batch": str(batch), "reason": "ambiguous PdbList row"}); continue
            wt_name = wt_names[0]
            mut_name = next(n for n in names if n != wt_name)
            identity = model_identity(mut_name)
            if not identity:
                failures.append({"batch": str(batch), "reason": "unmapped group/run"}); continue
            group, run = identity
            if not 1 <= group <= len(groups):
                failures.append({"batch": str(batch), "reason": "group out of range"}); continue
            identity = (group, run)
            if identity in mapping:
                raise ValueError("duplicate PdbList group/run")
            mapping[identity] = True
            # PdbList cannot escape its originating batch directory.
            if Path(wt_name).name != wt_name or Path(mut_name).name != mut_name:
                raise ValueError("PdbList path escapes batch")
            wt_path, mut_path = batch / wt_name, batch / mut_name
            if not wt_path.is_file() or not mut_path.is_file():
                failures.append({"batch": str(batch), "reason": "missing paired structure"}); continue
            mutations = []
            for notation in groups[group-1].strip().rstrip(";").split(","):
                m = re.fullmatch(r"([A-Z])(.)(-?\d+)([A-Z])", notation)
                if not m:
                    raise ValueError("invalid batch mutation identity")
                wt, chain, position, mutant = m.groups()
                mutations.append((chain, int(position), wt, mutant))
            w, m = atoms(wt_path), atoms(mut_path)
            for chain, position, wt, mutant in mutations:
                key = (chain, position, "", "CA")
                if key not in w or key not in m or protein_letters_3to1.get(w[key][1]) != wt or protein_letters_3to1.get(m[key][1]) != mutant:
                    raise ValueError("paired model mutation/WT identity mismatch")
            pairs.append({
                "mutations": [f"{wt}{pos}{mut}" for chain, pos, wt, mut in mutations],
                "chain": mutations[0][0], "group": group, "run": run,
                "wt_path": str(wt_path.resolve()), "mut_path": str(mut_path.resolve()),
                "wt_sha256": sha256_file(wt_path), "mut_sha256": sha256_file(mut_path),
                "pdb_list_path": str(listing.resolve()), "pdb_list_sha256": sha256_file(listing),
                "input_sha256": sha256_file(batch / "input.pdb"),
                "ddg_kcal_mol": ddg_by_identity.get(identity),
                "ddg_file_path":str(ddg_files[identity].resolve()) if identity in ddg_files else None,
                "ddg_file_sha256":sha256_file(ddg_files[identity]) if identity in ddg_files else None,
                "mutant_file_sha256":sha256_file(batch/'individual_list.txt'),
                "batch_path": str(batch.resolve()),
            })
    return pairs, failures


def validate_pair(row):
    if row.get('ddg_kcal_mol') is not None and not np.isfinite(row['ddg_kcal_mol']):
        raise ValueError('nonfinite paired DDG')
    if row.get('input_sha256') and sha256_file(Path(row['batch_path'])/'input.pdb')!=row['input_sha256']:
        raise ValueError('paired input fingerprint changed')
    for kind in ("wt", "mut"):
        if sha256_file(Path(row[f"{kind}_path"])) != row[f"{kind}_sha256"]:
            raise ValueError("paired structure fingerprint changed")
    # Independent replay checks coordinate residue identity as well as bytes;
    # a self-consistent but wrongly annotated manifest cannot certify a model.
    coordinate_models={kind:atoms(Path(row[f'{kind}_path'])) for kind in ('wt','mut')}
    for notation in row['mutations']:
        match=re.fullmatch(r'([A-Z])(-?\d+)([A-Z])',notation)
        if not match:raise ValueError('invalid paired mutation notation')
        wt,position,mutant=match.groups();key=(row['chain'],int(position),'','CA')
        for kind,letter in (('wt',wt),('mut',mutant)):
            value=coordinate_models[kind].get(key)
            if value is None or protein_letters_3to1.get(value[1])!=letter:
                raise ValueError('paired coordinate WT/mutation metadata identity mismatch')
    if sha256_file(Path(row["pdb_list_path"])) != row["pdb_list_sha256"]:
        raise ValueError("PdbList fingerprint changed")
    if row.get('ddg_file_sha256') and sha256_file(Path(row['ddg_file_path']))!=row['ddg_file_sha256']:
        raise ValueError('paired DDG table fingerprint changed')
    if row.get('ddg_kcal_mol') is not None:
        from ..stability import _parse_dif_rows
        batch=Path(row['batch_path'])
        matching=[v for p in batch.glob('Dif_*.fxout') for name,v in _parse_dif_rows(p) if name==Path(row['mut_path']).name]
        if len(matching)!=1 or abs(matching[0]-row['ddg_kcal_mol'])>1e-12:
            raise ValueError('paired DDG value/model identity mismatch')
    if row.get('mutant_file_sha256') and sha256_file(Path(row['batch_path'])/'individual_list.txt')!=row['mutant_file_sha256']:
        raise ValueError('mutation manifest changed')
    if row.get('mutant_file_sha256'):
        groups=(Path(row['batch_path'])/'individual_list.txt').read_text(encoding='ascii').splitlines()
        group=int(row['group']);identity=(group,int(row['run']))
        if not 1<=group<=len(groups):raise ValueError('paired group identity mismatch')
        expected=[];chains=set()
        for text in groups[group-1].strip().rstrip(';').split(','):
            match=re.fullmatch(r'([A-Z])(.)(-?\d+)([A-Z])',text)
            if not match:raise ValueError('invalid paired mutation record')
            wt,chain,pos,mut=match.groups();expected.append(f'{wt}{pos}{mut}');chains.add(chain)
        if sorted(expected)!=sorted(row['mutations']) or chains!={row['chain']}:
            raise ValueError('paired mutation metadata disagrees with batch')
        if any(model_identity(Path(row[key]).name)!=identity for key in ('wt_path','mut_path')):
            raise ValueError('paired filename/group/run identity mismatch')
