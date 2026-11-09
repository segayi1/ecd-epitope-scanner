from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .assay_profiles import DEFAULT_ASSAY_PROFILE_ID, SUPPORTED_ASSAY_PROFILE_IDS
from .pipeline import PipelineRequest, run_pipeline
from .util import WaitingForAF3


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ecd_epitope_scanner")
    sub = parser.add_subparsers(dest="command", required=True)
    serve = sub.add_parser("serve", help="Start the local web workbench")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8765)

    cleanup = sub.add_parser("cleanup", help="Preview or remove successful FoldX diagnostics")
    cleanup.add_argument(
        "--execute", action="store_true", help="Actually remove the previewed diagnostics"
    )
    cleanup.add_argument("--include-test-output", action="store_true")

    analyze = sub.add_parser("analyze", help="Run a reproducible ECD regional design")
    source = analyze.add_mutually_exclusive_group(required=True)
    source.add_argument("--uniprot")
    source.add_argument("--pdb-id")
    source.add_argument("--structure")
    source.add_argument("--af3-directory", help="Extracted official AlphaFold 3 output folder")
    analyze.add_argument("--accession")
    analyze.add_argument("--chain")
    analyze.add_argument("--analysis-target", choices=["antigen", "binder"], default="antigen",
                         help="Analyze an antigen ECD or a binder chain (complexes require --chain)")
    analyze.add_argument("--sequence")
    analyze.add_argument("--ecd", help="Override ECD as inclusive UniProt start:end")
    analyze.add_argument("--confidence-json")
    analyze.add_argument(
        "--iupred-output", help="Full canonical-sequence IUPred output to validate and import"
    )
    analyze.add_argument("--disable-reliability", action="store_true")
    analyze.add_argument(
        "--reliability-surface-mode", choices=["strict", "occlusion_sensitivity"], default="strict"
    )
    analyze.add_argument("--original-cds")
    analyze.add_argument("--allow-partial", action="store_true")
    analyze.add_argument("--allow-afdb-fallback", action="store_true")
    analyze.add_argument("--min-regions", type=int, default=5)
    analyze.add_argument("--max-regions", type=int, default=8)
    analyze.add_argument(
        "--design-engine",
        choices=["joint_v4_ml"],
        default="joint_v4_ml",
        help="Evidence-aware Joint V4-ML design engine",
    )
    analyze.add_argument(
        "--ml-mode",
        choices=["strict", "preview"],
        default="strict",
        help="Joint V4 external-evidence policy",
    )
    analyze.add_argument("--max-sectors", type=int, default=7)
    analyze.add_argument("--optimization-profile", choices=["legacy"], default="legacy")
    analyze.add_argument("--surface-strategy", choices=["legacy"], default="legacy")
    analyze.add_argument("--candidate-ids", nargs="+")
    analyze.add_argument("--calibration-repeat-schedule", nargs="+", type=int)
    analyze.add_argument("--export-mode", choices=["checkpoint", "deliverable"], default="deliverable")
    analyze.add_argument("--model-assets-root")
    analyze.add_argument("--case-time-budget-seconds",type=float)
    analyze.add_argument("--case-byte-budget",type=int)
    analyze.add_argument("--max-primary-constructs", type=int, default=7)
    analyze.add_argument("--allow-spurs-only-final", action="store_true")
    analyze.add_argument("--biological-context-mode", default="provided_structure")
    analyze.add_argument("--procaliper-worker", nargs="+", help="Joint V4 ProCaliper worker argv")
    analyze.add_argument("--spurs-worker", nargs="+", help="Joint V4 SPURS worker argv")
    analyze.add_argument("--discotope-worker", nargs="+", help="Joint V4 DiscoTope worker argv")
    analyze.add_argument(
        "--foldx-worker", nargs="+", help="Joint V4 FoldX confirmation worker argv"
    )
    analyze.add_argument(
        "--assay-profile",
        choices=SUPPORTED_ASSAY_PROFILE_IDS,
        help="Fixed dual-cell topology (automatic for a PDB-only --structure run)",
    )
    analyze.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "serve":
        local_hosts = {"127.0.0.1", "localhost", "::1"}
        if args.host.lower() not in local_hosts:
            raise SystemExit(
                "This local workbench only permits a loopback host. Public deployment requires a dedicated public runtime profile."
            )
        import uvicorn

        uvicorn.run("ecd_epitope_scanner.api:app", host=args.host, port=args.port, reload=False)
        return 0
    if args.command == "cleanup":
        from .api import STORE

        preview = STORE.cleanup_preview(args.include_test_output)
        print(
            json.dumps(
                {"candidates": preview, "bytes": sum(item["bytes"] for item in preview)},
                ensure_ascii=False,
            )
        )
        if args.execute:
            removed = STORE.cleanup_successful_engine_runs()
            print(json.dumps({"removed_engine_run_bytes": removed}, ensure_ascii=False))
        return 0

    ecd_start = ecd_end = None
    if args.ecd:
        ecd_start, ecd_end = (int(value) for value in args.ecd.split(":", 1))
    if args.uniprot:
        input_type, identifier = "uniprot", args.uniprot
    elif args.pdb_id:
        input_type, identifier = "pdb_id", args.pdb_id
    elif args.af3_directory:
        if args.confidence_json:
            raise SystemExit("--confidence-json cannot be combined with --af3-directory")
        input_type, identifier = "af3_upload", args.af3_directory
    else:
        input_type, identifier = (
            ("af3_upload" if args.confidence_json else "structure"),
            args.structure,
        )
    original_cds = None
    if args.original_cds:
        path = Path(args.original_cds)
        original_cds = path.read_text(encoding="utf-8") if path.exists() else args.original_cds
    sequence = None
    if args.sequence:
        path = Path(args.sequence)
        if path.exists():
            sequence = "".join(
                line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
            )
        else:
            sequence = args.sequence
    request = PipelineRequest(
        optimization_profile=args.optimization_profile,
        surface_strategy=args.surface_strategy,
        candidate_ids=args.candidate_ids,
        calibration_repeat_schedule=args.calibration_repeat_schedule,
        export_mode=args.export_mode,
        model_assets_root=args.model_assets_root,
        case_time_budget_seconds=args.case_time_budget_seconds,
        case_byte_budget=args.case_byte_budget,
        input_type=input_type,
        analysis_target=args.analysis_target,
        identifier=identifier,
        accession=args.accession,
        chain=args.chain,
        structure_path=args.structure,
        confidence_path=args.confidence_json,
        sequence=sequence,
        ecd_start=ecd_start,
        ecd_end=ecd_end,
        original_cds=original_cds,
        strict_af3=not args.allow_afdb_fallback,
        allow_partial=args.allow_partial,
        minimum_regions=args.min_regions,
        maximum_regions=args.max_regions,
        assay_profile=(
            args.assay_profile
            or (
                DEFAULT_ASSAY_PROFILE_ID
                if args.structure and not args.accession and not sequence
                else None
            )
        ),
        output_dir=args.output,
        design_engine=args.design_engine,
        ml_mode=args.ml_mode,
        max_sectors=args.max_sectors,
        max_primary_constructs=args.max_primary_constructs,
        allow_spurs_only_final=args.allow_spurs_only_final,
        biological_context_mode=args.biological_context_mode,
        metadata={
            "iupred_output": args.iupred_output,
            "structure_reliability": {
                "enabled": not args.disable_reliability,
                "surface_mode": args.reliability_surface_mode,
            },
            "joint_v4_ml": {
                "procaliper": {"command": args.procaliper_worker or []},
                "spurs_single": {"command": args.spurs_worker or []},
                "spurs_multi": {"command": args.spurs_worker or []},
                "discotope": {"command": args.discotope_worker or []},
                "foldx_confirmation": {"command": args.foldx_worker or []},
            },
        },
    )

    if args.af3_directory:
        from .af3_preprocessing import import_af3_directory

        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        request.af3_manifest_path = str(import_af3_directory(Path(args.af3_directory), output_dir))

    def progress(stage: str, value: float, message: str) -> None:
        print(f"[{value:6.1%}] {stage}: {message}", file=sys.stderr)

    try:
        result = run_pipeline(request, progress)
    except WaitingForAF3 as exc:
        print(
            json.dumps(
                {"status": "WAITING_FOR_AF3_UPLOAD", "message": str(exc)}, ensure_ascii=False
            )
        )
        return 2
    print(
        json.dumps(
            {
                "status": "COMPLETE",
                "output": str(Path(args.output).resolve()),
                "summary": result["summary"],
            },
            ensure_ascii=False,
        )
    )
    return 0
