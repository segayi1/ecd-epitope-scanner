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
from .sequences import (
    HUMAN_PREFERRED_CODONS,
    apply_mutations,
    coding_fragment,
    normalize_cds,
    translate_cds,
)
from .util import sha256_file, write_json

ROOT = Path(__file__).resolve().parent
WEB_ROOT = ROOT / "web"
OFFLINE_TEMPLATE = WEB_ROOT / "offline_template.html"
VENDOR_3DMOL = WEB_ROOT / "vendor" / "3Dmol-min.js"


def _js_safe(chunk: str) -> str:
    """Make an arbitrary text blob safe to inline inside a <script> block."""
    return chunk.replace("</script", "<\\/script").replace("<!--", "<\\!--")


def _json_inline(value: Any) -> str:
    """JSON for inline embedding; '</' is escaped the JSON-legal way."""
    return json.dumps(value, ensure_ascii=False).replace("</", "<\\/")


def build_offline_html(
    analysis: dict[str, Any],
    workdir: Path | None = None,
    lang: str = "zh",
) -> str:
    """Render one self-contained offline report (no network, file:// openable).

    Embeds 3Dmol.js, the ECD PDB text, the full analysis payload and both i18n
    dictionaries.  The offline page offers beads + 3Dmol-cartoon rendering
    (PDBe Mol* stays online-only: embedding it would roughly sextuple the file).
    """
    template = OFFLINE_TEMPLATE.read_text(encoding="utf-8")
    js = _js_safe(VENDOR_3DMOL.read_text(encoding="utf-8"))
    i18n = {
        name: json.loads((WEB_ROOT / "i18n" / f"{name}.json").read_text(encoding="utf-8"))
        for name in ("zh", "en")
    }
    pdb_path = (workdir or ROOT) / "structure_for_energy.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(f"structure_for_energy.pdb not found under {pdb_path.parent}")
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
    protein = analysis.get("protein") or {}
    task_name = protein.get("protein_name") or protein.get("accession") or "ECD Epitope Workbench"
    accession = protein.get("accession")
    task_label = f"{task_name} ({accession})" if accession and accession not in task_name else task_name

    html = template
    html = html.replace("__DATA__", _json_inline(analysis))
    html = html.replace("__I18N__", _json_inline(i18n))
    html = html.replace("__PDB__", json.dumps(pdb_text))
    html = html.replace("__LANG__", "en" if lang == "en" else "zh")
    html = html.replace("__TASK_TITLE__", f"{task_label} | ECD Epitope Workbench")
    html = html.replace("__TASK_BRAND__", f"{task_label} · ECD Epitope Workbench")
    # 3Dmol last: it is the biggest blob and the others' placeholders cannot
    # plausibly occur inside it, but the reverse is not true.
    html = html.replace("/*__3DMOL__*/", _js_safe(js))
    return html


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
    if analysis.get("assay_profile"):
        write_json(
            output_dir / "fixed_assay_profile.json",
            {
                "assay_profile": analysis["assay_profile"],
                "assay_geometry": analysis.get("assay_geometry"),
            },
        )
    # A real offline deliverable must be written and zipped with the scientific
    # artefacts.  Previously only the template function existed, leaving users
    # with no file://-openable report in a completed job.
    (output_dir / "offline_report.html").write_text(
        build_offline_html(analysis, output_dir), encoding="utf-8"
    )

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
            "flexible_terminal",
            "engineered_neighborhood",
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
                "parent_design_id": design.parent_design_id,
                "design_role": design.design_role,
                "status": design.status,
                "mutations": ";".join(item.notation for item in design.mutations),
                "anchor": design.anchor_notation,
                "stability_engine": design.stability_engine,
                "combined_stability_score": design.combined_stability_score,
                "score_origin": design.score_origin,
                "disruption_score": design.disruption_score,
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
            "parent_design_id",
            "design_role",
            "status",
            "mutations",
            "anchor",
            "stability_engine",
            "combined_stability_score",
            "score_origin",
            "disruption_score",
            "perturbation_score",
            "rationale",
            "warnings",
        ],
    )
    plan = analysis.get("epitope_mapping_plan") or {}
    _write_csv(
        output_dir / "panel_coverage_matrix.csv",
        plan.get("coverage_matrix", []),
        ["design_id", "patch_id", "perturbed_fraction"],
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
        record.annotations["date"] = dt.datetime.now(dt.UTC).strftime("%d-%b-%Y").upper()
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

    # Full-length constructs for cell-surface binding assays: flow cytometry on
    # intact cells needs the mutations expressed in the full membrane protein,
    # not an ECD fragment.
    full_length_records: list[SeqRecord] = []
    full_length_protein_lines: list[str] = []
    for design in designs:
        if not design.mutations:
            continue
        try:
            full_protein = apply_mutations(protein.sequence, design.mutations, 1)
        except Exception:
            continue  # a mutation outside the canonical sequence must never abort export
        full_length_protein_lines.append(
            f">{design.id} full-length {protein.accession} mutations="
             + ";".join(item.notation for item in design.mutations)
        )
        full_length_protein_lines.extend(textwrap.wrap(full_protein, 80))
        if original_cds:
            normalized = normalize_cds(original_cds)
            if translate_cds(normalized) == protein.sequence:
                codons = [normalized[i : i + 3] for i in range(0, len(normalized), 3)]
                for mutation in design.mutations:
                    codons[mutation.position - 1] = HUMAN_PREFERRED_CODONS[mutation.mutant]
                full_dna = "".join(codons)
                full_length_records.append(
                    SeqRecord(
                        Seq(full_dna),
                        id=f"{design.id}_full",
                        name=design.id[:16],
                        description=(
                            f"Full-length coding sequence for {protein.accession}; only the mapped "
                            "ECD codons are changed; verify signal peptide/TM boundaries before ordering"
                        ),
                    )
                )
    # Rebuilding a job must not leave stale full-length artefacts when the
    # source CDS is unavailable in a later custom-design save.
    for stale_name in ("full_length_protein_sequences.fasta", "full_length_coding_sequences.gb"):
        stale = output_dir / stale_name
        if stale.exists():
            stale.unlink()
    if full_length_protein_lines:
        (output_dir / "full_length_protein_sequences.fasta").write_text(
            "\n".join(full_length_protein_lines) + "\n", encoding="ascii"
        )
    if full_length_records:
        SeqIO.write(full_length_records, output_dir / "full_length_coding_sequences.gb", "genbank")

    methods = analysis.get("methods", {})
    partition = analysis.get("partitioning", {})
    selection = analysis.get("structure_selection", {})
    plan = analysis.get("epitope_mapping_plan", {}) or {}
    coverage_plan = partition.get("coverage_plan", [])
    tier_lines = "\n".join(
        f"- Patch {row.get('patch_id')}: tier {row.get('tier')}, epitope prior {row.get('epitope_prior')}, "
        f"probes required {row.get('required_probes')}, predicted accessibility "
        f"{row.get('predicted_accessibility', 'unknown')} ({row.get('membrane_accessibility_score', 'n/a')})"
        for row in coverage_plan
    ) or "- (no patch plan available)"
    probe_lines = "\n".join(
        f"- {probe.get('design_id')} (patch {probe.get('patch_id')}): "
        f"{', '.join(probe.get('mutations', [])) or 'none'}; anchor {probe.get('anchor') or 'none'}"
        for probe in plan.get("stage1_probes", [])
    ) or "- (no probes)"
    blind_lines = "\n".join(
        f"- Patch {item.get('patch_id')}: {item.get('reason')}"
        for item in plan.get("blind_spots", [])
    ) or "- (none; every patch carries at least one probe)"
    preset_lines = "\n".join(
        f"- {item.get('design_id')} [{item.get('role')}] parent {item.get('parent_design_id')}: "
        f"{', '.join(item.get('mutations', []))}"
        for item in plan.get("stage2_presets", [])
    ) or "- (none)"
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
- Engineered mismatch positions: {', '.join(str(value) for value in methods.get('engineered_mismatch_positions', [])) or 'none'}

## Surface and folding-risk screen

- SASA method: {methods.get('sasa_method', '')}
- Surface criterion: rSASA >= {methods.get('rsasa_threshold', 0.25)} and side-chain SASA >= {methods.get('sidechain_sasa_threshold_A2', 20)} A2
- Protected features: disulfide cysteines, annotated glycosylation/PTMs, missing mappings and low-confidence AF3 residues
- Native interfaces participate in design but retain a function-risk flag

## Patch partition and probe design

- Partition method: {partition.get('method', '')}
- Partition mode: {partition.get('partition_mode', 'footprint')}
- Partition basis: {partition.get('partition_basis', 'residue surface graph')}
- Regions: {partition.get('actual_regions', len(regions))}; local Fab-footprint geodesic radius {partition.get('patch_radius_a', 13)} A
- Partition trigger: {'; '.join(partition.get('trigger_reasons', [])) or 'standard footprint workflow'}
- Assay format: {partition.get('assay_format', 'cell_surface')}
- Cell-surface accessibility model: {partition.get('membrane_accessibility_model', 'legacy binary model')}
- Eligible-surface coverage: {partition.get('coverage', 0):.1%}
- Flexible terminal residues excluded from patches: {partition.get('excluded_flexible_terminal_residues', 0)} (screened separately as linear alanine-scan candidates, see epitope_mapping_plan)
- In-silico substitutions evaluated: {analysis.get('design_summary', {}).get('in_silico_substitutions_evaluated', 0)}
- Final constructs: {len(designs)}
- Stability interpretation: {analysis.get('design_summary', {}).get('score_interpretation', '')}

### Patch tiers and epitope priors

{tier_lines}

### Stage-1 probe panel

Selection policy: {plan.get('selection_policy', 'Pareto front over (disruption, stability); primary = knee, alternate = aggressive end.')}

{probe_lines}

### Declared blind spots

{blind_lines}

### Stage-2 presets (order only after a Stage-1 hit)

{preset_lines}

## How to read the binding assay

{plan.get('assay_note', 'See analysis.json epitope_mapping_plan.')}

- panel_coverage_matrix.csv holds the construct-by-patch perturbation matrix; use it to weigh negative results.
- Identical loss patterns across antibodies are binning evidence, not proof of identical epitopes.
- A patch with no probe is UNTESTED, not negative.

## Required experimental gates

1. Confirm total and surface expression relative to WT.
2. Confirm a non-competing conformational antibody or another domain-specific folding control remains bound.
3. Measure monodispersity/thermal stability when using purified ECD.
4. Measure native ligand/function when functional-interface mutations are interpreted.
5. Treat loss of all antibodies or controls as a folding/expression failure, not an epitope.
"""
    if analysis.get('structure_reliability'):
        report += "\n## Structural reliability routing\n\nMulti-evidence heuristic states are recorded in residue_reliability.csv and structural_segments.csv; they are not calibrated disorder probabilities. See reliability_qc.json for thresholds, PAE mapping validity and imported sequence evidence provenance. Original SASA and stability methods are unchanged. In strict surface mode excluded ECD atoms do not occlude the partition surface; boundary residues retain occlusion but contribute no sampled vertices. IDR sequence windows are review-only, with SASA unavailable and no automatic FoldX evaluation. Missing coordinates are UNRESOLVED, not evidence of disorder.\n"
    if analysis.get("assay_profile"):
        geometry = analysis.get("assay_geometry") or {}
        report += (
            "\n## Fixed dual-cell assay topology\n\n"
            "Antibody H/K are displayed on Jurkat; WT or mutant antigen is displayed on K562. "
            "The antigen ECD is N-terminal/extracellular and C-terminally anchored through "
            "PDGFRB. The Fab is extracellular and its heavy chain is C-terminally anchored "
            "through the CD8A/Notch1 receptor. Reporter signal is read on Jurkat.\n\n"
            f"PDB orientation status: {geometry.get('status', 'unavailable')}. "
            "The C-terminal-anchor axis is a relative label only, not a physical membrane "
            "plane or a quantitative Fab-reach estimate. Membrane-proximal regions are not "
            "down-weighted. N-terminal Myc controls antigen display but does not prove native "
            "folding; require an independent conformation-sensitive control.\n"
        )
    (output_dir / "methods_report.md").write_text(report, encoding="utf-8")

    export_names = [
        "analysis.json",
        "mapping.csv",
        "residues.csv",
        "regions.json",
        "designs.json",
        "fixed_assay_profile.json",
        "offline_report.html",
        "designs.csv",
        "panel_coverage_matrix.csv",
        "sequences.fasta",
        "coding_sequences.fasta",
        "sequences.gb",
        "full_length_protein_sequences.fasta",
        "full_length_coding_sequences.gb",
        "methods_report.md",
        "structure.cif",
        "ecd_structure.cif",
        "surface_vertices.csv",
        "surface_edges.npz",
        "raw_footprints.json",
        "footprints.json",
        "footprints.csv",
        "footprint_membership.csv",
        "footprint_adjacency.csv",
        "surface_regions.json",
        "surface_region_membership.csv",
        "surface_v2_qc.json",
        "residue_reliability.csv",
        "structural_segments.csv",
        "surface_mask.csv",
        "idr_candidates.csv",
        "reliability_qc.json",
    ]
    manifest = {
        "schema_version": "1.0",
        "created_at": dt.datetime.now(dt.UTC).isoformat(),
        "software": analysis.get("software", {}),
        "request": analysis.get("request", {}),
        "source": analysis.get("source", {}),
        "assay_profile": analysis.get("assay_profile"),
        "assay_geometry": analysis.get("assay_geometry"),
        "parameters": analysis.get("parameters", {}),
        # Auditable solver/engine provenance (P1-3/P1-4): which SASA path and
        # radius set produced the numbers, and whether the FoldX engine was
        # actually used or silently degraded.
        "sasa_method": (analysis.get("methods") or {}).get("sasa_method", ""),
        "stability": {
            "engine": (analysis.get("stability") or {}).get("primary", ""),
            "foldx_available": bool((analysis.get("stability") or {}).get("foldx_available")),
            "evoef2_available": bool((analysis.get("stability") or {}).get("evoef2_available")),
            "foldx_diagnostics": (analysis.get("design_summary") or {}).get("foldx_diagnostics", {}),
        },
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
