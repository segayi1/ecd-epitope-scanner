from __future__ import annotations

import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
from typing import Any

from .models import JointV3Result


def _jsonable(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return {key: _jsonable(item) for key, item in value.__dict__.items()}
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def write_joint_v3_outputs(output_dir: Path, result: JointV3Result, parameters: dict[str, Any]) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    payloads = {
        "structure_audit.json": result.audit,
        "topology_context.json": result.topology,
        "microfootprints.json": result.microfootprints,
        "probe_sectors.json": result.sectors,
        "selected_panel.json": result.panel,
    }
    for name, value in payloads.items():
        (output_dir / name).write_text(json.dumps(_jsonable(value), ensure_ascii=False, indent=2), encoding="utf-8")
    with (output_dir / "surface_residues.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["position", "chain_id", "auth_seq_id", "rasa", "mapped_area_a2", "designable", "component_id", "exclusion_reasons"])
        writer.writeheader()
        for unit in result.surface_units:
            writer.writerow({"position": unit.position, "chain_id": unit.residue.chain_id, "auth_seq_id": unit.residue.auth_seq_id, "rasa": unit.rasa, "mapped_area_a2": unit.mapped_area_a2, "designable": unit.designable, "component_id": unit.component_id, "exclusion_reasons": ";".join(unit.exclusion_reasons)})
    with (output_dir / "candidate_probes.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "sector_id", "mutations", "mutation_evidence", "connected", "structural_risk", "stability_value", "stability_unit", "stability_status", "stability_provenance", "evidence_level", "coverage_by_scenario", "warnings"])
        writer.writeheader()
        for candidate in result.candidates:
            writer.writerow({
                "candidate_id": candidate.candidate_id,
                "sector_id": candidate.sector_id,
                "mutations": ";".join(item.notation for item in candidate.mutations),
                "mutation_evidence": json.dumps([
                    {"notation": item.notation, "perturbation_class": item.perturbation_class, "rule_risk": item.rule_risk, "evidence_level": item.evidence_level}
                    for item in candidate.mutations
                ], sort_keys=True),
                "connected": candidate.connected,
                "structural_risk": candidate.structural_risk,
                "stability_value": candidate.stability_value,
                "stability_unit": candidate.stability_unit,
                "stability_status": candidate.stability_status,
                "stability_provenance": json.dumps(candidate.stability_provenance, sort_keys=True),
                "evidence_level": candidate.evidence_level,
                "coverage_by_scenario": json.dumps(candidate.coverage_by_scenario, sort_keys=True),
                "warnings": ";".join(candidate.warnings),
            })
    with (output_dir / "coverage_by_scenario.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["candidate_id", "scenario", "covered_positions", "covered_unit_count", "designable_area_a2", "fraction_of_designable_surface", "fraction_of_total_surface", "footprint_count"])
        writer.writeheader()
        for candidate in result.candidates:
            for scenario, values in sorted(candidate.coverage_by_scenario.items()):
                writer.writerow({
                    "candidate_id": candidate.candidate_id,
                    "scenario": scenario,
                    "covered_positions": ";".join(str(item) for item in candidate.coverage_units_by_scenario.get(scenario, [])),
                    "covered_unit_count": len(candidate.coverage_units_by_scenario.get(scenario, [])),
                    "designable_area_a2": values.get("designable_area_a2"),
                    "fraction_of_designable_surface": values.get("fraction_of_designable_surface"),
                    "fraction_of_total_surface": values.get("fraction_of_total_surface"),
                    "footprint_count": values.get("footprint_count"),
                })
    with (output_dir / "exclusions.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sector_id", "position", "reason"])
        writer.writeheader()
        for sector_id, entries in sorted(result.candidate_exclusions.items()):
            for entry in entries:
                position, separator, reason = entry.partition(":")
                writer.writerow({
                    "sector_id": sector_id,
                    "position": position if separator and position.isdigit() else "",
                    "reason": reason if separator else entry,
                })
    input_fingerprint_payload = {
        "audit": _jsonable(result.audit),
        "surface_residues": [
            [item.residue.model_id, item.residue.chain_id, item.residue.auth_seq_id, item.residue.insertion_code, item.residue.one_letter, item.position, item.mapped_area_a2]
            for item in result.surface_units
        ],
    }
    summary = {
        "schema_version": "3.0",
        "created_at": dt.datetime.now(dt.UTC).isoformat(),
        "engine": "joint_v3",
        "configuration": parameters,
        "configuration_sha256": hashlib.sha256(json.dumps(parameters, sort_keys=True).encode()).hexdigest(),
        "analysis_input_fingerprint_sha256": hashlib.sha256(json.dumps(input_fingerprint_payload, sort_keys=True).encode()).hexdigest(),
        "solver_status": result.panel.solver_status,
        "surface_graph": result.graph_provenance,
        "source_provenance": result.source_provenance,
        "stability": result.design_summary.get("foldx_diagnostics", {}),
        "runtime": {
            "python": platform.python_version(),
            "numpy": importlib.metadata.version("numpy"),
            "scipy": importlib.metadata.version("scipy"),
            "biopython": importlib.metadata.version("biopython"),
        },
        "warnings": ["V3 provides coarse localization evidence, not exact antibody epitopes or binding probabilities."],
    }
    (output_dir / "manifest.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (output_dir / "scientific_report.md").write_text(
        "# Joint V3 scientific report\n\n"
        "This report distinguishes measured structure features, geometry proxies, and unknowns. "
        "Probe sectors are coarse localization regions; they are not antibody footprints.\n\n"
        f"- Selected probes: {len(result.panel.selected_candidate_ids)}\n"
        f"- Solver status: `{result.panel.solver_status}`\n"
        f"- Robust coverage scenarios: {json.dumps(result.panel.achieved_coverage, ensure_ascii=False)}\n"
        f"- Surface graph: `{result.graph_provenance.get('basis', 'unknown')}`\n"
        f"- Total/designable surface area: {sum(item.mapped_area_a2 for item in result.surface_units):.1f} / {sum(item.mapped_area_a2 for item in result.surface_units if item.designable):.1f} A^2\n"
        f"- Protected surface residues: {sum(not item.designable for item in result.surface_units)}\n"
        f"- Candidate exclusions: {sum(len(value) for value in result.candidate_exclusions.values())}\n"
        f"- Alternative panels: {json.dumps(result.panel.alternative_panels, ensure_ascii=False)}\n\n"
        "Experimental interpretation requires WT, no-antigen/irrelevant-antigen, expression-normalized binding, and an independent conformation-sensitive control. A positive probe can reflect local folding, glycosylation, or allostery; a negative probe does not exclude epitope membership.\n",
        encoding="utf-8",
    )
    return {name: name for name in [*payloads, "surface_residues.csv", "candidate_probes.csv", "coverage_by_scenario.csv", "exclusions.csv", "manifest.json", "scientific_report.md"] if (output_dir / name).exists()}
