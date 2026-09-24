from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import traceback
from pathlib import Path


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _binary(root: Path, request: dict) -> Path:
    configured = str(request.get("parameters", {}).get("binary") or "").strip()
    if configured:
        path = Path(configured)
        if not path.is_absolute():
            path = root / path
        if path.is_file():
            return path
        raise RuntimeError(f"Configured FoldX binary does not exist: {path}")
    matches = sorted((root / "vendor" / "foldx").glob("foldx_*.exe"))
    if not matches:
        raise RuntimeError("No FoldX executable was found under vendor/foldx")
    return matches[-1]


def _parse_mutations(mutations: list[list[str]], chain_id: str, sequence: str):
    groups = []
    for group in mutations:
        parsed = []
        for notation in group:
            match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", notation.upper())
            if not match:
                raise RuntimeError(f"Invalid FoldX mutation notation: {notation}")
            wt, position_text, mutant = match.groups()
            position = int(position_text)
            if not 1 <= position <= len(sequence) or sequence[position - 1] != wt:
                raise RuntimeError(f"FoldX WT/position mismatch: {notation}")
            parsed.append((wt, position, chain_id, mutant))
        groups.append(parsed)
    return groups


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        root = _project_root()
        sys.path.insert(0, str(root.parent))
        from ecd_epitope_scanner.stability import repair_with_foldx, score_many_with_foldx

        mutations = request.get("mutations") or []
        if not mutations:
            raise RuntimeError("FoldX received no mutation combinations")
        parameters = request.get("parameters", {})
        binary = _binary(root, request)
        binary_digest = hashlib.sha256()
        with binary.open("rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                binary_digest.update(chunk)
        output_dir = Path(request["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        structure_path = Path(request["structure_path"])
        repaired = repair_with_foldx(
            str(binary), structure_path, output_dir,
            binary_sha256=binary_digest.hexdigest(),
        )
        if repaired is None:
            raise RuntimeError("FoldX RepairPDB failed")
        groups = _parse_mutations(
            mutations, str(request["chain_id"]), str(request["sequence"])
        )
        runs = max(1, int(parameters.get("number_of_runs", 5)))
        workers = max(1, int(parameters.get("max_workers", 8)))
        values = score_many_with_foldx(
            str(binary),
            repaired,
            groups,
            output_dir,
            runs=runs,
            already_repaired=True,
            max_workers=workers,
        )
        if len(values) != len(mutations) or any(value is None for value in values):
            missing = sum(value is None for value in values)
            raise RuntimeError(
                f"FoldX returned incomplete scores: {missing} of {len(mutations)} missing"
            )
        records = [
            {"mutations": group, "ddg": float(values[index])}
            for index, group in enumerate(mutations)
        ]
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": [],
                    "errors": [],
                    "provenance": {
                        "worker": "foldx_worker",
                        "binary": str(binary),
                        "version": "5.1",
                        "weights_sha256": binary_digest.hexdigest(),
                        "device": "cpu",
                        "command": "RepairPDB + BuildModel",
                        "number_of_runs": runs,
                        "max_workers": workers,
                        "repair_cache_path": str(repaired),
                    },
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        print(
            json.dumps(
                {
                    "status": "failed",
                    "records": [],
                    "warnings": [],
                    "errors": [{"code": "TOOL_CRASHED", "message": str(exc)}],
                    "provenance": {"worker": "foldx_worker"},
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
