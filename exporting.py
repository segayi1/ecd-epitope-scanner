from __future__ import annotations

import csv
import datetime as dt
import json
import re
import shutil
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


OFFLINE_SITE_ASSETS = (
    "workbench.css",
    "app.css",
    "renderer.js",
    "vendor/3Dmol-min.js",
    "vendor/pdbe-molstar-plugin-3.10.1.js",
    "vendor/pdbe-molstar-light-3.10.1.css",
    "vendor/PDBe-Molstar-LICENSE",
)


def _offline_classic_bundle() -> str:
    """Build the app's existing single-import frontend as a classic file:// bundle."""
    app = (WEB_ROOT / "app.js").read_text(encoding="utf-8")
    recommendations = (WEB_ROOT / "recommendations.js").read_text(encoding="utf-8")
    import_pattern = r'^import\s*\{\s*createRecommendationRenderer\s*\}\s*from\s*["\']\.\/recommendations\.js(?:\?[^"\']*)?["\'];\s*'
    app, import_count = re.subn(import_pattern, "", app, count=1, flags=re.MULTILINE)
    if import_count != 1 or re.search(r'^\s*import\s|^\s*export\s', app, re.MULTILINE):
        raise ValueError("Offline bundle assertion failed: app.js import/export contract changed")
    recommendations, export_count = re.subn(
        r'^export\s+function\s+createRecommendationRenderer\b',
        "function createRecommendationRenderer",
        recommendations,
        count=1,
        flags=re.MULTILINE,
    )
    if export_count != 1 or re.search(r'^\s*(?:import|export)\s', recommendations, re.MULTILINE):
        raise ValueError("Offline bundle assertion failed: recommendations.js module contract changed")
    return "/* local offline classic bundle; generated from the online app sources */\n" + recommendations + "\n" + app


def build_offline_site(output_dir: Path, analysis: dict[str, Any], lang: str = "zh") -> list[str]:
    """Build a self-contained copy of the same workbench for direct file:// use.

    Scientific data and structures are embedded in the root page so no local
    fetch, localhost API, or network access is needed to open a report.
    """
    from .joint_design_v4.reporting import project_display_analysis

    assets_dir = output_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    for rel_name in OFFLINE_SITE_ASSETS:
        source = WEB_ROOT / rel_name
        if not source.is_file():
            raise FileNotFoundError(f"Required offline frontend asset is missing: {source}")
        target = assets_dir / rel_name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    (assets_dir / "app.bundle.js").write_text(_offline_classic_bundle(), encoding="utf-8")

    structure_path = output_dir / "structure.cif"
    if not structure_path.is_file():
        structure_path = output_dir / "ecd_structure.cif"
    pdb_path = output_dir / "structure_for_energy.pdb"
    if not pdb_path.is_file():
        raise FileNotFoundError(f"No PDB structure is available for offline viewer in {output_dir}")
    if not structure_path.is_file():
        structure_path = pdb_path

    i18n = {
        name: json.loads((WEB_ROOT / "i18n" / f"{name}.json").read_text(encoding="utf-8"))
        for name in ("zh", "en")
    }
    page = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    page = re.sub(
        r'(?P<attr>href|src)="/assets/(?P<path>[^"]+)"',
        lambda match: f'{match.group("attr")}="./assets/{match.group("path").split("?", 1)[0]}"',
        page,
    )
    page = page.replace('href="/"', 'href="./index.html"')
    page = re.sub(r'href="/(?!assets/)[^"]*"', 'href="#offline-unavailable"', page)
    available_assets = set(OFFLINE_SITE_ASSETS) | {"app.bundle.js"}
    page = re.sub(
        r'href="\./assets/([^"?]+)(?:\?[^\"]*)?"',
        lambda match: match.group(0) if match.group(1) in available_assets else 'href="#offline-unavailable"',
        page,
    )
    page = re.sub(
        r'<script type="module" src="\./assets/app\.js[^\"]*"></script>',
        '<script src="./assets/app.bundle.js"></script>',
        page,
    )
    if 'type="module"' in page or re.search(
        r'(?:href|src)="(?:/|https?:|//)|/api/', page, flags=re.IGNORECASE
    ):
        raise ValueError("Offline site assertion failed: absolute/module network resource remains")
    page = page.replace(
        '<meta charset="utf-8">',
        '<meta charset="utf-8">\n<meta http-equiv="Content-Security-Policy" '
        'content="default-src \'self\' file: data: blob:; connect-src \'self\' file: blob: data:; '
        'script-src \'self\' \'unsafe-inline\' \'unsafe-eval\' file: blob:; '
        'style-src \'self\' \'unsafe-inline\' file:; img-src \'self\' file: data: blob:; '
        'font-src \'self\' file: data:; worker-src \'self\' file: blob:">',
    )
    snapshot = {
        "analysis": project_display_analysis(analysis),
        "i18n": i18n,
        "structure_data": structure_path.read_text(encoding="utf-8", errors="replace"),
        "structure_pdb": pdb_path.read_text(encoding="utf-8", errors="replace"),
        "language": "en" if lang == "en" else "zh",
    }
    embedded = "<script>window.APP_OFFLINE = " + _json_inline(snapshot) + ";</script>\n"
    page = page.replace('<body data-workspace-step="1">', '<body class="offline-mode" data-workspace-step="1">\n' + embedded)
    (output_dir / "index.html").write_text(page, encoding="utf-8")
    return ["index.html", *(f"assets/{name}" for name in OFFLINE_SITE_ASSETS), "assets/app.bundle.js"]


def build_offline_html(
    analysis: dict[str, Any],
    workdir: Path | None = None,
    lang: str = "zh",
) -> str:
    """Render one self-contained offline report (no network, file:// openable).

    Embeds 3Dmol.js, the ECD PDB text, the display analysis and both i18n
    dictionaries.  The offline page offers beads + 3Dmol-cartoon rendering
    (PDBe Mol* stays online-only: embedding it would roughly sextuple the file).
    """
    from .joint_design_v4.reporting import project_display_analysis

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
    html = html.replace("__DATA__", _json_inline(project_display_analysis(analysis)))
    html = html.replace("__I18N__", _json_inline(i18n))
    html = html.replace("__PDB__", json.dumps(pdb_text))
    html = html.replace("__LANG__", "en" if lang == "en" else "zh")
    html = html.replace("__TASK_TITLE__", f"{task_label} | ECD Epitope Workbench")
    html = html.replace("__TASK_BRAND__", f"{task_label} · ECD Epitope Workbench")
    # 3Dmol last: it is the biggest blob and the others' placeholders cannot
    # plausibly occur inside it, but the reverse is not true.
    html = html.replace("/*__3DMOL__*/", _js_safe(js))
    recommendation_js = (WEB_ROOT / "recommendations.js").read_text(encoding="utf-8")
    html = html.replace("/*__RECOMMENDATIONS__*/", _js_safe(recommendation_js.replace("export function", "function", 1)))
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


def write_analysis_artifact(output_dir: Path, analysis: dict[str, Any]) -> None:
    """Persist the display projection; canonical evidence remains in sidecars."""
    from .joint_design_v4.reporting import project_display_analysis

    write_json(
        output_dir / "analysis.json",
        project_display_analysis(analysis, compact_bulk_identity=False),
    )


def write_scientific_outputs(
    output_dir: Path,
    analysis: dict[str, Any],
    protein: ProteinContext,
    records: list[ResidueRecord],
    regions: list[RegionRecord],
    designs: list[DesignRecord],
    original_cds: str | None = None,
    *, defer_archive: bool = False,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_artifact(output_dir, analysis)
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
    offline_site_names = build_offline_site(output_dir, analysis)

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

## Structure admission and chain selection

- Admission: {selection.get('admission', {}).get('status', 'not recorded')}
- Chain selection: {selection.get('chain_selection', {}).get('method', 'not recorded')} ({selection.get('chain_selection', {}).get('confidence', 'unknown')} confidence)
- Admission warnings: {'; '.join(selection.get('admission', {}).get('warnings', [])) or 'none'}
- Chain-selection warnings: {'; '.join(selection.get('chain_selection', {}).get('warnings', [])) or 'none'}
- PDB/UniProt sequence differences: {', '.join(str(value) for value in selection.get('engineered_mismatch_positions', [])) or 'none'}; original coordinates retained and differing positions omitted from per-residue scoring, surface partitioning and mutation design

## Surface and folding-risk screen

- SASA method: {methods.get('sasa_method', '')}
- Surface criterion: native total rSASA >= {methods.get('rsasa_threshold', 0.25)}; side-chain SASA is descriptive only. Coordinate eligibility and confidence review are separate from mutation stability.
- Protected features: disulfide cysteines, annotated glycosylation/PTMs and missing mappings. Geometry-ineligible segments remain excluded; eligible internal low-confidence loops require review and are excluded from the strict formal panel.
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
    af3 = (analysis.get("methods") or {}).get("af3_preprocessing")
    if af3:
        baselines = af3.get("baselines", {})
        report += ("\n## AF3 multi-model preprocessing\n\n"
            f"Selected original sample: `{af3.get('selected_model_id')}`. "
            f"Selection status: `{af3.get('status')}`.\n\n"
            f"Original top rank: `{baselines.get('original_rank1')}`; "
            f"best chain pTM: `{baselines.get('ptm_best')}`; "
            f"raw weighted medoid: `{baselines.get('raw_medoid')}`.\n\n"
            "Original coordinates and the paired confidence file are retained. "
            "Core consensus, global disagreement and local diagnostics are recorded separately; "
            "these diagnostics do not add residue exclusion rules or change scoring budgets.\n")
        report += ("\n| Model | Coverage | f70 | Mean pLDDT | pTM | Original ranking | Consensus (Å) | Status |\n"
                   "|---|---:|---:|---:|---:|---:|---:|---|\n")
        for model in af3.get("models", []):
            cells = [model.get(key) for key in ("model_id", "coverage", "f70", "mean_plddt", "ptm", "ranking_score", "consensus_score")]
            cells.append("eligible" if model.get("eligible") else "; ".join(model.get("exclusion_reasons", [])))
            report += "| " + " | ".join("—" if value is None else str(value).replace("|", "\\|") for value in cells) + " |\n"
        report += "\nFull diagnostics and input identities: `af3_selection_report.json` and `af3_inputs/manifest.json`.\n"
    (output_dir / "methods_report.md").write_text(report, encoding="utf-8")

    export_names = [
        "analysis.json",
        "mapping.csv",
        "residues.csv",
        "regions.json",
        "designs.json",
        "fixed_assay_profile.json",
        "offline_report.html",
        *offline_site_names,
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
        # Joint V3 machine-readable audit files.  They are optional so legacy
        # runs keep the same output contract.
        "manifest.json",
        "structure_audit.json",
        "topology_context.json",
        "surface_residues.csv",
        "microfootprints.json",
        "probe_sectors.json",
        "candidate_probes.csv",
        "coverage_by_scenario.csv",
        "selected_panel.json",
        "exclusions.csv",
        "scientific_report.md",
        # Joint V4-ML evidence and compatibility projection files.
        "tool_inventory.json",
        "structure_contexts.json",
        "residue_mapping.csv",
        "prediction_context.pdb",
        "surface_evidence.csv",
        "sasa_comparison.csv",
        "spurs_single_matrix.csv.gz",
        "single_mutation_options.csv",
        "mutation_menu_v2.json",
        "mutation_menu_v2.csv",
        "candidate_pool_v2.json",
        "discotope_residues.csv",
        "epitope_islands.json",
        "probe_sectors_v4.json",
        "candidate_combinations.csv.gz",
        "combination_stability.csv",
        "panel_candidates.json",
        "ordinary_candidate_ledger.json",
        "selected_panel_v4.json",
        "sector_recommendations_v4.json",
        "sector_recommendations_v4.csv",
        "v4_manifest.json",
        "paired_model_manifest.json",
        "v4_scientific_report.md",
        "v4_manifest_v41.json",
        "surface_graph_v4.json",
        "surface_graph_v41.json",
        "surface_display_regions_v41.json",
        "partition_diagnostics_v4.json",
        "surface_hypotheses_v4.json",
        "hypothesis_probe_matrix_v4.json",
        "epitope_hypotheses_v41.json",
        "hypothesis_membership_v41.csv",
        "probe_hypothesis_matrix_v41.csv",
        "panel_identifiability_v41.json",
        "recommendation_diagnostics_v41.json",
        "qualification_audit_v411.json",
        "hypothesis_hierarchy_v411.json",
        "global_candidate_pool_v411.json",
        "review_candidates_v411.json",
        "panel_budget_decision_v411.json",
        "panel_budget_curve_v411.csv",
        "raw_hypotheses_v411.json",
        "raw_surface_footprints_v411.json",
        "footprint_families_v411.json",
        "inference_hypotheses_v411.json",
        "surface_complexity_v411.json",
        "retrospective_benchmark_contract_v411.json",
        "surface_perturbation_preview.json",
        "surface_perturbation_local_constructs.json",
        "raw_null_import_manifest.json",
        "surface_perturbation_local_constructs.csv",
        "surface_perturbation_evidence.zip",
        "frozen_input_manifest.json",
        "launch_source_manifest.json",
        "source_confidence.json",
        "source_input_identity.json",
        "native_input.cif",
        "native_input.pdb",
        "structure_for_energy.pdb",
    ]
    export_names.extend(path.relative_to(output_dir).as_posix()
        for path in sorted((output_dir / "foldx_refinement").glob("**/paired_model_manifest.json"))
        if not {"prior_attempts", "paired_manifest_history"}.intersection(path.relative_to(output_dir).parts))
    if af3:
        export_names.extend((analysis.get("source", {}).get("af3_preprocessing") or {}).get("files", []))
        export_names.extend(path.relative_to(output_dir).as_posix()
            for path in sorted((output_dir / "af3_inputs").rglob("*")) if path.is_file())
        export_names.extend(path.name for path in sorted(output_dir.glob("af3_preprocessing*")) if path.is_file())
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
    if not defer_archive:
        _write_bundle_archive(output_dir,export_names)
    return {name: name for name in export_names + [zip_path.name]
            if (output_dir / name).exists() or defer_archive and name==zip_path.name}


def _write_bundle_archive(output_dir: Path,names) -> None:
    import os
    destination=output_dir/'ecd_epitope_design.zip'
    pending=destination.with_name(destination.name+f'.{os.getpid()}.pending')
    try:
        with zipfile.ZipFile(pending,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1) as archive:
            for name in dict.fromkeys(names):
                path=output_dir/name
                if path.exists():
                    compression=zipfile.ZIP_STORED if path.suffix in {'.zip','.npz','.gz'} else zipfile.ZIP_DEFLATED
                    archive.write(path,arcname=name,compress_type=compression)
        pending.replace(destination)
    finally:
        pending.unlink(missing_ok=True)


def finalize_scientific_bundle(output_dir: Path, analysis: dict | None = None) -> dict | None:
    """Bind checksums and the ZIP to the final analysis written by the pipeline."""
    manifest_path = output_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    analysis_path = output_dir / "analysis.json"
    if analysis is None:
        analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    v4_manifest_path = output_dir / "v4_manifest.json"
    if v4_manifest_path.is_file():
        v4_manifest = json.loads(v4_manifest_path.read_text(encoding="utf-8"))
        if v4_manifest.get("worker_runs", {}).get("foldx", {}).get("status") == "success":
            from .foldx_source_evidence import copy_foldx_source_evidence
            try:
                evidence = copy_foldx_source_evidence(output_dir, output_dir, v4_manifest)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                evidence = {"status": "UNKNOWN", "reason": str(exc)}
            v4_manifest.setdefault("provenance", {})["portable_source_evidence"] = evidence
            analysis["portable_source_evidence"] = evidence
            if evidence.get("status") != "UNKNOWN":
                for name in evidence["artifact_paths"]:
                    analysis.setdefault("files", {})[name] = name
                    manifest.setdefault("files", {})[name] = ""
            write_json(v4_manifest_path, v4_manifest)
            write_analysis_artifact(output_dir, analysis)
    has_structure = (output_dir / "structure_for_energy.pdb").is_file()
    if has_structure:
        (output_dir / "offline_report.html").write_text(
            build_offline_html(analysis, output_dir), encoding="utf-8"
        )
        manifest.setdefault("files", {}).setdefault("offline_report.html", "")
    site_names = (build_offline_site(output_dir, analysis)
        if analysis.get("joint_v4_ml") or has_structure else [])
    for name in site_names:
        manifest.setdefault("files", {}).setdefault(name, "")
    card_receipt = None
    card_receipt_path = output_dir / "card_evidence_completeness.json"
    if card_receipt_path.is_file():
        from .joint_design_v4.impact_completeness import write_card_acceptance
        card_receipt = write_card_acceptance(output_dir, analysis)
        manifest.setdefault("files", {})[card_receipt_path.name] = ""
    manifest["files"].pop("run_manifest.json", None)
    manifest["files"].pop("ecd_epitope_design.zip", None)
    for name in manifest["files"]:
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Invalid bundle artifact name: {name}")
        artifact = output_dir / relative
        if not artifact.resolve().is_relative_to(output_dir.resolve()) or artifact.is_symlink():
            raise ValueError(f"Bundle artifact escapes output directory: {name}")
        if not artifact.is_file():
            raise FileNotFoundError(f"Bundle artifact missing at finalization: {name}")
        manifest["files"][name] = sha256_file(artifact)
    write_json(manifest_path, manifest)
    _write_bundle_archive(output_dir,[*manifest['files'],'run_manifest.json'])
    return card_receipt
