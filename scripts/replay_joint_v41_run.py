"""Replay a saved run's selected structure and canonical annotations with live workers.

Only remote input retrieval is replaced with the saved data. Surface geometry,
candidate generation, all evidence workers, gates, selection and exports run
normally. Historical run files are never changed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ecd_epitope_scanner.pipeline import PipelineRequest, run_pipeline
from ecd_epitope_scanner.retrieval import parse_uniprot_entry
from ecd_epitope_scanner.util import sha256_file, write_json


def replay_request_data(source: Path, output: Path, analysis: dict, request_data: dict) -> dict:
    """Preserve predicted input semantics and confidence during frozen replay."""
    data = dict(request_data)
    predicted_upload = data.get("input_type") == "af3_upload"
    structure_path = source / "structure.cif"
    if predicted_upload:
        structure_path = source / "af3_upload" / "official_af3_model.cif"
        confidence_path = source / "af3_upload" / "official_af3_confidences.json"
        if not structure_path.is_file() or not confidence_path.is_file():
            raise ValueError("AF3 replay requires the original model and validated confidence file")
        data["confidence_path"] = str(confidence_path)
    data.update(
        input_type="af3_upload" if predicted_upload else "structure",
        identifier=data.get("identifier") if predicted_upload else str(structure_path),
        structure_path=str(structure_path), accession=analysis["protein"]["accession"],
        chain=analysis["source"].get("selected_chain") or analysis["structure_selection"].get("selected_chain"),
        output_dir=str(output),
    )
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_run", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--profile-v3", action="store_true", help="Explicitly replay frozen inputs under the new adaptive V3 design policy")
    parser.add_argument("--surface-rsasa-threshold", type=float, choices=(0.15, 0.20, 0.25))
    parser.add_argument("--geometry-only", action="store_true", help="Recompute real geometry and exports with all prediction workers explicitly disabled; no mutation recommendation")
    parser.add_argument("--optimization-profile", choices=["legacy", "coverage_scan_v2_preview", "coverage_scan_v3_preview"], default="legacy")
    parser.add_argument("--candidate-ids", nargs="+")
    parser.add_argument("--calibration-repeat-schedule", nargs="+", type=int)
    parser.add_argument("--export-mode", choices=["checkpoint", "deliverable"], default="deliverable")
    args = parser.parse_args()
    source = args.source_run.resolve()
    output = args.output_dir.resolve()
    if source == output or source in output.parents or output in source.parents:
        parser.error("Replay output must be separate from the historical run")
    if (output / "analysis.json").exists():
        parser.error("Replay output already contains an analysis; choose a new directory")
    analysis = json.loads((source / "analysis.json").read_text(encoding="utf-8"))
    entry = json.loads((source / "uniprot.json").read_text(encoding="utf-8"))
    request_data = json.loads((source / "request.json").read_text(encoding="utf-8"))
    protein = analysis["protein"]
    request_data = replay_request_data(source, output, analysis, request_data)
    request_data.update(optimization_profile=args.optimization_profile, candidate_ids=args.candidate_ids,
        calibration_repeat_schedule=args.calibration_repeat_schedule, export_mode=args.export_mode)
    if args.profile_v3:
        request_data.update(assay_profile="jurkat_fab_k562_antigen_v3", max_primary_constructs=None,
                            minimum_constructs=4, target_constructs=7, maximum_constructs=8, adaptive_budget=True)
    if args.surface_rsasa_threshold is not None:
        request_data["metadata"] = {**request_data.get("metadata", {}), "surface_rsasa_threshold": args.surface_rsasa_threshold}
    if args.geometry_only:
        request_data["metadata"] = {
            **request_data.get("metadata", {}),
            "disable_external_engines": True,
            "validation_purpose": "geometry_only_replay_no_mutation_recommendation",
        }
    conservation = {
        item["uniprot_pos"]: item["conservation_score"]
        for item in analysis["residues"] if item.get("conservation_score") is not None
    }

    def frozen_protein(accession, client, output_dir, ecd_override=None):
        del client
        if accession != protein["accession"]:
            raise ValueError("Replay accession does not match frozen canonical annotations")
        write_json(output_dir / "uniprot.json", entry)
        return parse_uniprot_entry(entry, ecd_override)

    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "replay_input_provenance.json", {
        "source_run": str(source),
        "source_structure_sha256": sha256_file(source / "structure.cif"),
        "replayed_model_sha256": sha256_file(Path(request_data["structure_path"])),
        "source_uniprot_sha256": sha256_file(source / "uniprot.json"),
        "source_request_sha256": sha256_file(source / "request.json"),
        "source_confidence_sha256": (
            sha256_file(Path(request_data["confidence_path"]))
            if request_data["input_type"] == "af3_upload" else None
        ),
        "remote_inputs": "frozen canonical entry, selected structure and conservation",
        "scientific_computation": (
            "live geometry and exports; prediction workers disabled; no mutation recommendation"
            if args.geometry_only else "live geometry, workers, candidate gates, panel and exports"
        ),
    })
    with patch("ecd_epitope_scanner.pipeline.fetch_protein_context", frozen_protein), patch(
        "ecd_epitope_scanner.conservation.fetch_conservation", return_value=conservation or None
    ):
        result = run_pipeline(PipelineRequest(**request_data), lambda stage, value, message: print(
            f"{stage} {value:.0%}: {message}", flush=True
        ))
    print(json.dumps({
        "recommendation_status": result["design_summary"]["recommendation_status"],
        "diagnostics": result["design_summary"]["recommendation_diagnostics"],
        "output_dir": str(output),
    }, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
