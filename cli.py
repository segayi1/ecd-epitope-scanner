from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .pipeline import PipelineRequest, run_pipeline
from .util import WaitingForAF3


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ecd_epitope_scanner")
    sub = parser.add_subparsers(dest="command", required=True)
    serve = sub.add_parser("serve", help="Start the local web workbench")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8765)

    analyze = sub.add_parser("analyze", help="Run a reproducible ECD regional design")
    source = analyze.add_mutually_exclusive_group(required=True)
    source.add_argument("--uniprot")
    source.add_argument("--pdb-id")
    source.add_argument("--structure")
    analyze.add_argument("--accession")
    analyze.add_argument("--chain")
    analyze.add_argument("--sequence")
    analyze.add_argument("--ecd", help="Override ECD as inclusive UniProt start:end")
    analyze.add_argument("--confidence-json")
    analyze.add_argument("--original-cds")
    analyze.add_argument("--allow-partial", action="store_true")
    analyze.add_argument("--allow-afdb-fallback", action="store_true")
    analyze.add_argument("--min-regions", type=int, default=5)
    analyze.add_argument("--max-regions", type=int, default=8)
    analyze.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "serve":
        import uvicorn

        uvicorn.run("ecd_epitope_scanner.api:app", host=args.host, port=args.port, reload=False)
        return 0

    ecd_start = ecd_end = None
    if args.ecd:
        ecd_start, ecd_end = (int(value) for value in args.ecd.split(":", 1))
    if args.uniprot:
        input_type, identifier = "uniprot", args.uniprot
    elif args.pdb_id:
        input_type, identifier = "pdb_id", args.pdb_id
    else:
        input_type, identifier = ("af3_upload" if args.confidence_json else "structure"), args.structure
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
        input_type=input_type,
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
        output_dir=args.output,
    )

    def progress(stage: str, value: float, message: str) -> None:
        print(f"[{value:6.1%}] {stage}: {message}", file=sys.stderr)

    try:
        result = run_pipeline(request, progress)
    except WaitingForAF3 as exc:
        print(json.dumps({"status": "WAITING_FOR_AF3_UPLOAD", "message": str(exc)}, ensure_ascii=False))
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
