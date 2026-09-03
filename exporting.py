from __future__ import annotations

import csv
import datetime as dt
import json
import textwrap
import zipfile
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .models import DesignRecord, ProteinContext, RegionRecord, ResidueRecord, to_dict
from .sequences import coding_fragment
from .util import sha256_file, write_json


def _join(value: Any) -> Any:
    if isinstance(value, list):
        return ";".join(str(item) for item in value)
    return value


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _join(row.get(field, "")) for field in fields})


def write_scientific_outputs(
    output_dir: Path,
    analysis: dict[str, Any],
    protein: ProteinContext,
    records: list[ResidueRecord],
    regions: list[RegionRecord],
    designs: list[DesignRecord],
    original_cds: str | None = None,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "analysis.json", analysis)
    write_json(output_dir / "regions.json", [to_dict(item) for item in regions])
    write_json(output_dir / "designs.json", [to_dict(item) for item in designs])

    residue_rows = [to_dict(item) for item in records]
    _write_csv(
        output_dir / "residues.csv",
        residue_rows,
        [
            "uniprot_pos",
            "aa",
            "chain_id",
            "auth_seq_id",
            "insertion_code",
            "label_chain_id",
            "label_seq_id",
            "sasa",
            "sidechain_sasa",
            "rsasa",
            "contacts",
            "plddt",
            "secondary_structure",
            "domain",
            "surface_accessible",
            "folding_screen_pass",
            "interface",
            "interface_partners",
            "glycan_shielded",
            "glycosylation_site",
            "glycosylation_sequon",
            "disulfide",
            "ptm",
            "salt_bridge",
            "domain_interface",
            "region_id",
            "hard_exclusions",
            "risk_flags",
            "reasons",
        ],
    )
    _write_csv(
        output_dir / "mapping.csv",
        residue_rows,
        [
            "uniprot_pos",
            "aa",
            "chain_id",
            "auth_seq_id",
            "insertion_code",
            "label_chain_id",
            "label_seq_id",
        ],
    )
    design_rows = []
    for design in designs:
        design_rows.append(
            {
                "id": design.id,
                "name": design.name,
                "region_id": design.region_id,
                "status": design.status,
                "mutations": ";".join(item.notation for item in design.mutations),
                "stability_engine": design.stability_engine,
                "combined_stability_score": design.combined_stability_score,
                "perturbation_score": design.perturbation_score,
                "rationale": design.rationale,
                "warnings": design.warnings,
            }
        )
    _write_csv(
        output_dir / "designs.csv",
        design_rows,
        [
            "id",
            "name",
            "region_id",
            "status",
            "mutations",
            "stability_engine",
            "combined_stability_score",
            "perturbation_score",
            "rationale",
            "warnings",
        ],
    )

    wt_dna, dna_source = coding_fragment(protein, protein.ecd_sequence, original_cds)
    with (output_dir / "sequences.fasta").open("w", encoding="ascii") as handle:
        handle.write(f">WT_ECD {protein.accession}:{protein.ecd_start}-{protein.ecd_end}\n")
        handle.write("\n".join(textwrap.wrap(protein.ecd_sequence, 80)) + "\n")
        for design in designs:
            mutations = ";".join(item.notation for item in design.mutations) or "none"
            handle.write(f">{design.id} {mutations} status={design.status}\n")
            handle.write("\n".join(textwrap.wrap(design.mutant_ecd_sequence, 80)) + "\n")
    with (output_dir / "coding_sequences.fasta").open("w", encoding="ascii") as handle:
        handle.write(f">WT_ECD_coding_fragment source={dna_source}\n")
        handle.write("\n".join(textwrap.wrap(wt_dna, 90)) + "\n")
        for design in designs:
            handle.write(f">{design.id}_coding_fragment\n")
            handle.write("\n".join(textwrap.wrap(design.mutant_coding_sequence, 90)) + "\n")

    genbank_records: list[SeqRecord] = []
    all_sequences = [("WT_ECD", wt_dna, [])] + [
        (design.id, design.mutant_coding_sequence, design.mutations) for design in designs
    ]
    for identifier, dna, mutations in all_sequences:
        record = SeqRecord(
            Seq(dna),
            id=identifier,
            name=identifier[:16],
            description=f"ECD coding fragment for {protein.accession}; no signal peptide/start/stop/tag",
        )
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "linear"
        record.annotations["data_file_division"] = "SYN"
        record.annotations["date"] = dt.datetime.now(dt.timezone.utc).strftime("%d-%b-%Y").upper()
        for mutation in mutations:
            relative = mutation.position - protein.ecd_start
            record.features.append(
                SeqFeature(
                    FeatureLocation(relative * 3, relative * 3 + 3),
                    type="variation",
                    qualifiers={"note": mutation.notation, "replace": mutation.mutant},
                )
            )
        genbank_records.append(record)
    SeqIO.write(genbank_records, output_dir / "sequences.gb", "genbank")

    methods = analysis.get("methods", {})
    partition = analysis.get("partitioning", {})
    selection = analysis.get("structure_selection", {})
    report = f"""# ECD epitope-region mutation design report

## Interpretation boundary

This workflow designs antigen mutants for experimental testing. It does not infer an antibody-specific epitope and it does not prove that a mutant is correctly folded. A construct supports epitope attribution only when antibody binding decreases while expression, folding controls and appropriate native-function controls remain acceptable.

## Target and structure

- UniProt: {protein.accession}
- Protein: {protein.protein_name}
- ECD: {protein.ecd_start}-{protein.ecd_end} ({protein.ecd_length} aa)
- Topology: {protein.topology}
- Structure: {selection.get('selected_identifier', 'user-provided')}
- Chain: {selection.get('selected_chain', '')}
- Observed ECD coverage: {selection.get('observed_ecd_coverage', 0):.1%}
- Maximum internal coordinate gap: {selection.get('max_internal_gap', 'n/a')}

## Surface and folding-risk screen

- SASA method: {methods.get('sasa_method', '')}
- Surface criterion: rSASA >= {methods.get('rsasa_threshold', 0.25)} and side-chain SASA >= {methods.get('sidechain_sasa_threshold_A2', 20)} A2
- Protected features: disulfide cysteines, annotated glycosylation/PTMs, missing mappings and low-confidence AF3 residues
- Native interfaces participate in design but retain a function-risk flag

## Partition and design

- Partition method: {partition.get('method', '')}
- Regions: {partition.get('actual_regions', len(regions))}
- Eligible-surface coverage: {partition.get('coverage', 0):.1%}
- In-silico substitutions evaluated: {analysis.get('design_summary', {}).get('in_silico_substitutions_evaluated', 0)}
- Final constructs: {len(designs)}
- Stability interpretation: {analysis.get('design_summary', {}).get('score_interpretation', '')}

## Required experimental gates

1. Confirm total and surface expression relative to WT.
2. Confirm a non-competing conformational antibody or another domain-specific folding control remains bound.
3. Measure monodispersity/thermal stability when using purified ECD.
4. Measure native ligand/function when functional-interface mutations are interpreted.
5. Treat loss of all antibodies or controls as a folding/expression failure, not an epitope.
"""
    (output_dir / "methods_report.md").write_text(report, encoding="utf-8")

    export_names = [
        "analysis.json",
        "mapping.csv",
        "residues.csv",
        "regions.json",
        "designs.json",
        "designs.csv",
        "sequences.fasta",
        "coding_sequences.fasta",
        "sequences.gb",
        "methods_report.md",
        "structure.cif",
        "ecd_structure.cif",
    ]
    manifest = {
        "schema_version": "1.0",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "software": analysis.get("software", {}),
        "request": analysis.get("request", {}),
        "source": analysis.get("source", {}),
        "parameters": analysis.get("parameters", {}),
        "files": {
            name: sha256_file(output_dir / name)
            for name in export_names
            if (output_dir / name).exists()
        },
    }
    write_json(output_dir / "run_manifest.json", manifest)
    export_names.append("run_manifest.json")
    zip_path = output_dir / "ecd_epitope_design.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in export_names:
            path = output_dir / name
            if path.exists():
                archive.write(path, arcname=name)
    return {name: name for name in export_names + [zip_path.name] if (output_dir / name).exists()}
