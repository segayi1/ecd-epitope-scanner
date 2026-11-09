from __future__ import annotations

from collections import defaultdict
from typing import Any

from Bio.PDB.Polypeptide import is_aa

from ..models import ProteinContext, ResidueRecord
from .models import ResidueId, StructureAudit


def audit_structure(
    records: list[ResidueRecord],
    protein: ProteinContext,
    chain_id: str,
    model_id: str = "0",
    chain_ids: list[str] | None = None,
    structure: Any | None = None,
) -> StructureAudit:
    """Create an auditable, residue-level structure summary.

    Parsing and sequence mapping remain the responsibility of the existing
    pipeline.  This stage deliberately reports uncertainty instead of trying
    to infer missing coordinates or a canonical chain from ambiguous input.
    """
    observed = sorted(item.uniprot_pos for item in records)
    missing_segments: list[dict[str, Any]] = []
    if observed:
        run_start = previous = observed[0]
        for position in observed[1:]:
            if position > previous + 1:
                missing_segments.append({"start": run_start, "end": previous})
                run_start = position
            previous = position
        missing_segments.append({"start": run_start, "end": previous})
    expected = set(range(protein.ecd_start, protein.ecd_end + 1))
    missing = sorted(expected - set(observed))
    if missing:
        missing_segments = []
        start = previous = missing[0]
        for position in missing[1:]:
            if position > previous + 1:
                missing_segments.append({"start": start, "end": previous})
                start = position
            previous = position
        missing_segments.append({"start": start, "end": previous})

    nonstandard: list[dict[str, Any]] = []
    duplicate_groups: list[list[str]] = []
    available_chains = list(chain_ids or [])
    model_count = 1
    alternate_location_count = 0
    partial_occupancy_atom_count = 0
    low_information_regions: list[dict[str, Any]] = []
    if structure is not None:
        models = list(structure.get_models())
        model_count = len(models)
        selected_model = models[0] if models else None
        chain_sequences: dict[str, str] = {}
        if selected_model is not None:
            available_chains = [str(chain.id) for chain in selected_model.get_chains()]
            ca_only_residues: list[dict[str, Any]] = []
            for chain in selected_model.get_chains():
                letters = []
                for residue in chain.get_residues():
                    if not is_aa(residue, standard=False):
                        continue
                    if not is_aa(residue, standard=True):
                        nonstandard.append({
                            "chain_id": str(chain.id),
                            "auth_seq_id": int(residue.id[1]),
                            "insertion_code": str(residue.id[2]).strip(),
                            "residue_name": str(residue.resname),
                        })
                    letters.append(str(residue.resname))
                    atoms = list(residue.get_atoms())
                    if len(atoms) == 1 and atoms[0].get_name() == "CA":
                        ca_only_residues.append({
                            "chain_id": str(chain.id),
                            "auth_seq_id": int(residue.id[1]),
                            "insertion_code": str(residue.id[2]).strip(),
                        })
                    alternate_location_count += sum(bool(str(atom.get_altloc()).strip()) for atom in atoms)
                    partial_occupancy_atom_count += sum(
                        atom.get_occupancy() is not None and float(atom.get_occupancy()) < 0.999
                        for atom in atoms
                    )
                chain_sequences[str(chain.id)] = "|".join(letters)
            if ca_only_residues:
                low_information_regions.append({
                    "type": "ca_only_residues", "residues": ca_only_residues
                })
        by_sequence: dict[str, list[str]] = defaultdict(list)
        for candidate_chain, sequence in chain_sequences.items():
            if sequence:
                by_sequence[sequence].append(candidate_chain)
        duplicate_groups = [sorted(values) for values in by_sequence.values() if len(values) > 1]
    else:
        nonstandard = [
            {"position": item.uniprot_pos, "residue_name": item.residue_name}
            for item in records
            if item.aa == "X"
        ]
    warnings: list[str] = []
    if missing:
        warnings.append("ECD contains unmapped or missing residues; no virtual residues were added.")
    if nonstandard:
        warnings.append("Non-standard residue names are retained as audit records and excluded from mutation menus.")
    warnings.append("B-factors are not interpreted as pLDDT without explicit structure provenance.")
    if model_count > 1:
        warnings.append(f"Input contains {model_count} models; model {model_id} was selected.")
    if alternate_location_count:
        warnings.append("Alternate atom locations are present and follow the parser-selected conformer.")
    if len(available_chains) > 1 and not duplicate_groups:
        warnings.append("Multiple non-equivalent chains are present; the upstream mapping-selected target chain is recorded explicitly.")
    return StructureAudit(
        selected_model=model_id,
        selected_chains=[chain_id],
        chain_decision=("explicit_or_sequence_mapping_selection" if chain_id else "pipeline_selected_chain"),
        residue_mapping_status="COMPLETE_ECD" if not missing else "PARTIAL_ECD",
        missing_segments=missing_segments,
        nonstandard_residues=nonstandard,
        duplicate_chain_groups=duplicate_groups,
        quality_warnings=warnings,
        fatal_errors=[],
        model_count=model_count,
        available_chains=available_chains,
        low_information_regions=low_information_regions,
        alternate_location_count=alternate_location_count,
        partial_occupancy_atom_count=partial_occupancy_atom_count,
    )


def residue_id(record: ResidueRecord, model_id: str = "0") -> ResidueId:
    return ResidueId(
        model_id=model_id,
        chain_id=record.chain_id,
        auth_seq_id=record.auth_seq_id,
        insertion_code=record.insertion_code,
        one_letter=record.aa,
    )
