from __future__ import annotations

import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from typing import Any
from uuid import uuid4

from . import config
from .models import MutationChoice, ProteinContext, ResidueRecord


def _engine_run_dir(workdir: Path, prefix: str) -> Path:
    """Persistent per-invocation working directory for an external engine.

    FoldX drops ~25 model PDBs plus its own rotabase.txt/molecules copies per
    invocation.  Deleting that many files has proven fatal inside sandboxed
    hosts (bulk-delete guards kill the process AFTER the scores were computed),
    so run directories are intentionally KEPT: they live under the job folder,
    cost a few MB, and double as diagnostics for engine failures.  Never
    replace this with TemporaryDirectory cleanup.
    """
    root = Path(workdir) / "engine_runs"
    root.mkdir(parents=True, exist_ok=True)
    run = root / f"{prefix}_{uuid4().hex[:8]}"
    run.mkdir(parents=True, exist_ok=False)
    return run


def _run_engine(command: list[str], cwd: Path, timeout: int) -> subprocess.CompletedProcess | None:
    """Run an external engine step; a hang degrades to None instead of raising.

    The upstream binaries can hang on modern Windows hosts (see the EvoEF2
    opt-in note below).  A hung engine must cost a degraded score, never the
    whole job, so TimeoutExpired is caught here at every call site.
    """
    try:
        return subprocess.run(
            command,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return None


AA_PROPERTIES = {
    # hydropathy, side-chain volume A^3, charge, aromatic
    "A": (1.8, 67.0, 0, 0),
    "R": (-4.5, 148.0, 1, 0),
    "N": (-3.5, 96.0, 0, 0),
    "D": (-3.5, 91.0, -1, 0),
    "C": (2.5, 86.0, 0, 0),
    "Q": (-3.5, 114.0, 0, 0),
    "E": (-3.5, 109.0, -1, 0),
    "G": (-0.4, 48.0, 0, 0),
    "H": (-3.2, 118.0, 0.2, 1),
    "I": (4.5, 124.0, 0, 0),
    "L": (3.8, 124.0, 0, 0),
    "K": (-3.9, 135.0, 1, 0),
    "M": (1.9, 124.0, 0, 0),
    "F": (2.8, 135.0, 0, 1),
    "P": (-1.6, 90.0, 0, 0),
    "S": (-0.8, 73.0, 0, 0),
    "T": (-0.7, 93.0, 0, 0),
    "W": (-0.9, 163.0, 0, 1),
    "Y": (-1.3, 141.0, 0, 1),
    "V": (4.2, 105.0, 0, 0),
}
AMINO_ACIDS = tuple(AA_PROPERTIES)
AUTO_EXCLUDED_MUTANTS = {"C", "P", "G"}


def sequons(sequence: str) -> set[int]:
    return {
        index + 1
        for index in range(len(sequence) - 2)
        if sequence[index] == "N" and sequence[index + 1] != "P" and sequence[index + 2] in "ST"
    }


def changes_sequon(protein: ProteinContext, position: int, mutant: str) -> bool:
    before = sequons(protein.sequence)
    sequence = list(protein.sequence)
    sequence[position - 1] = mutant
    after = sequons("".join(sequence))
    relevant = range(max(1, position - 2), min(len(sequence) - 1, position) + 1)
    return any(index in before.symmetric_difference(after) for index in relevant)


def _menu_class(record: ResidueRecord, mutant: str) -> str:
    """Directed substitution class for the mapping assay.

    Bulky/hotspot WT residues are best probed by Ala removal (Bogan & Thorn
    1998 hot-spot enrichment); small WT residues need bulk insertion (Ala is
    nearly neutral there); charged WT residues are probed by full charge
    reversal, the most reliable interface disruptor.
    """
    wt = record.aa
    if wt in "KR" and mutant in "ED":
        return "charge_reversal"
    if wt in "ED" and mutant in "KR":
        return "charge_reversal"
    if wt in config.HOTSPOT_WT and mutant == "A":
        return "hotspot_to_ala"
    if wt in "WYF" and mutant == "A":
        return "bulky_to_ala"
    if wt in config.SMALL_WT and mutant == "W":
        return "trp_insert"
    if wt in config.BULKY_WT and mutant == "A":
        return "bulky_to_ala"
    return "other"


def disruption_score(record: ResidueRecord, mutant: str) -> float:
    """Instrument-power of one substitution for binding-disruption assays.

    Dimensionless by design: this is a probe-strength prior, not an energy.
    """
    score = config.D_W_EXPOSURE * min(1.2, record.rsasa * 2.0)
    volume_delta = AA_PROPERTIES[mutant][1] - AA_PROPERTIES[record.aa][1]
    score += config.D_W_VOLUME_GAIN * abs(volume_delta) / 85.0
    charge_delta = abs(AA_PROPERTIES[mutant][2] - AA_PROPERTIES[record.aa][2])
    if charge_delta >= 1.0:
        score += config.D_W_CHARGE_FLIP
    elif charge_delta > 0:
        score += config.D_W_PARTIAL_CHARGE
    if record.aa in config.HOTSPOT_WT and mutant == "A":
        score += config.D_W_HOTSPOT_REMOVAL
    if mutant in "WFY" and record.aa not in "WFY":
        score += config.D_W_AROMATIC_INSERT
    return round(score, 4)


def is_anchor(record: ResidueRecord, mutant: str) -> bool:
    """Anchors are substitutions with reliable interface-disruption chemistry:
    full charge reversal, hotspot removal (W/Y/F -> A), or charge removal
    from K/R (K/R -> A strips a full charge from the interface)."""
    menu = _menu_class(record, mutant)
    if menu in {"charge_reversal", "hotspot_to_ala"}:
        return True
    return menu == "bulky_to_ala" and record.aa in "KR"


def evaluate_substitution(
    record: ResidueRecord,
    mutant: str,
    protein: ProteinContext,
) -> MutationChoice:
    mutant = mutant.upper()
    wt = record.aa
    if mutant not in AA_PROPERTIES:
        raise ValueError(f"Unsupported amino acid: {mutant}")
    if mutant == wt:
        choice = MutationChoice(record.uniprot_pos, wt, mutant, 99.0, 0.0, 1.0, False, ["WT-to-WT"])
        choice.menu_class = "other"
        return choice
    wt_h, wt_v, wt_q, wt_ar = AA_PROPERTIES[wt]
    new_h, new_v, new_q, new_ar = AA_PROPERTIES[mutant]
    density = min(1.5, record.contacts / 14.0)
    exposure = min(1.2, max(0.0, record.rsasa))

    risk = max(0.0, new_v - wt_v) / 85.0 * density
    risk += max(0.0, new_h - 1.5) * exposure * 0.18
    if record.salt_bridge and abs(new_q - wt_q) >= 0.8:
        risk += 0.9
    if record.domain_interface:
        risk += abs(new_v - wt_v) / 180.0 + abs(new_q - wt_q) * 0.2
    if record.secondary_structure == "E" and mutant in {"P", "G"}:
        risk += 2.5
    if mutant in AUTO_EXCLUDED_MUTANTS:
        risk += 1.6
    if mutant == "W" and record.contacts >= 10:
        risk += 0.6

    perturbation = (
        abs(new_h - wt_h) / 4.5
        + abs(new_v - wt_v) / 90.0
        + abs(new_q - wt_q) * 1.2
        + abs(new_ar - wt_ar) * 0.45
    )
    menu = _menu_class(record, mutant)
    disruption = disruption_score(record, mutant)
    reasons: list[str] = []
    allowed = True
    if mutant in AUTO_EXCLUDED_MUTANTS:
        reasons.append(f"auto-design excludes introduced {mutant}")
        allowed = False
    if record.aa in {"P", "G"}:
        reasons.append("native backbone-sensitive residue")
        allowed = False
    if changes_sequon(protein, record.uniprot_pos, mutant):
        reasons.append("creates or disrupts an N-X-S/T sequon")
        risk += 2.0
        allowed = False
    if record.disulfide or record.glycosylation_site or record.ptm:
        reasons.append("protected structural/PTM site")
        allowed = False
    # Exposure-linked risk gate: a fully exposed side chain is structurally
    # cheap to change, so the same proxy number must not veto it as eagerly as
    # it vetoes a partially buried one.  Gates live in config (see the FoldX
    # RMSE citation in docs/改进方案.md).
    gate = config.risk_gate(record.rsasa)
    if risk > gate:
        reasons.append(f"high structural-risk proxy (>{gate:.2f})")
        allowed = False
    if record.conservation_score is not None and record.conservation_score >= 0.95:
        reasons.append("fully conserved residue; function/epitope prior reduced")
        risk += 0.3
        disruption = round(disruption * 0.8, 4)
    if record.interface:
        reasons.append("native interface; function-risk flag retained")
    if not reasons:
        reasons.append("passes deterministic chemistry and structural-risk constraints")
    choice = MutationChoice(
        position=record.uniprot_pos,
        wt=wt,
        mutant=mutant,
        stability_score=round(float(risk), 4),
        perturbation_score=round(float(perturbation), 4),
        stability_percentile=0.0,
        allowed=allowed,
        reasons=reasons,
    )
    choice.menu_class = menu
    choice.disruption_score = disruption
    return choice


def assign_percentiles(choices: list[MutationChoice]) -> None:
    ordered = sorted(choices, key=lambda item: (item.stability_score, item.mutant))
    denominator = max(1, len(ordered) - 1)
    for rank, choice in enumerate(ordered):
        choice.stability_percentile = round(rank / denominator, 4)


@dataclass(slots=True)
class EngineInventory:
    proxy: bool
    evoef2_path: str | None
    foldx_path: str | None

    @property
    def primary_name(self) -> str:
        if self.evoef2_path and self.foldx_path:
            return "EvoEF2 + FoldX"
        if self.foldx_path:
            return "FoldX 5.1"
        if self.evoef2_path:
            return "EvoEF2"
        return "builtin structural-risk proxy"


def find_engine_inventory() -> EngineInventory:
    # External engines are opt-in: the upstream 32-bit Windows binary can hang on
    # modern Windows hosts, so silently executing an unprobed bundled binary would
    # make the local job queue non-terminating. Set EVOEF2_PATH after validating it.
    evoef = os.environ.get("EVOEF2_PATH") or shutil.which("EvoEF2") or shutil.which("EvoEF2.exe")
    bundled_foldx = sorted(
        (Path(__file__).resolve().parent / "vendor" / "foldx").glob("foldx_*.exe")
    )
    foldx = (
        os.environ.get("FOLDX_PATH")
        or shutil.which("FoldX")
        or shutil.which("foldx")
        or (str(bundled_foldx[-1]) if bundled_foldx else None)
    )
    return EngineInventory(True, evoef, foldx)


def _parse_evoef_total(text: str) -> float | None:
    for line in reversed(text.splitlines()):
        if "total" not in line.lower() and "reference" not in line.lower():
            continue
        numbers = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?", line)
        if numbers:
            return float(numbers[-1])
    numbers = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", text)
    return float(numbers[-1]) if numbers else None


def score_with_evoef2(
    binary: str,
    pdb_path: Path,
    mutations: list[tuple[str, int, str, str]],
    workdir: Path,
) -> float | None:
    """Return EvoEF2 mutant-minus-WT energy units, or None when output is unparseable."""
    temp = _engine_run_dir(workdir, "evoef2")
    with nullcontext(temp):
        local_pdb = temp / "input.pdb"
        shutil.copy2(pdb_path, local_pdb)
        repair = _run_engine(
            [binary, "--command=RepairStructure", f"--pdb={local_pdb.name}", "--num_of_runs=3"],
            temp,
            300,
        )
        repaired = temp / "input_Repair.pdb"
        if repair is None or repair.returncode != 0 or not repaired.exists():
            return None
        mutation_text = ",".join(f"{wt}{chain}{number}{mutant}" for wt, number, chain, mutant in mutations) + ";\n"
        (temp / "individual_list.txt").write_text(mutation_text, encoding="ascii")
        build = _run_engine(
            [
                binary,
                "--command=BuildMutant",
                f"--pdb={repaired.name}",
                "--mutant_file=individual_list.txt",
                "--num_of_runs=10",
            ],
            temp,
            600,
        )
        mutant_files = sorted(temp.glob("*_Model_1.pdb"))
        if build is None or build.returncode != 0 or not mutant_files:
            return None

        def stability(path: Path) -> float | None:
            result = _run_engine(
                [binary, "--command=ComputeStability", f"--pdb={path.name}"],
                temp,
                300,
            )
            if result is None or result.returncode:
                return None
            return _parse_evoef_total(result.stdout + "\n" + result.stderr)

        wt_score = stability(repaired)
        mutant_score = stability(mutant_files[0])
        if wt_score is None or mutant_score is None:
            return None
        return round(mutant_score - wt_score, 4)


def _parse_dif_fxout(path: Path) -> list[float]:
    """Parse one FoldX BuildModel ``Dif_*.fxout`` into per-run delta-delta-G values.

    The file opens with a FoldX banner whose first line is ``FoldX 5.1 (2011)``.
    Splitting that line and reading ``fields[1]`` yields ``5.1``, which silently
    poisons the median of the first construct.  Real data rows are the ones that
    follow the ``Pdb<TAB>total energy...`` header and begin with a PDB filename,
    so both conditions are enforced here.
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = 0
    for index, line in enumerate(lines):
        if line.startswith("Pdb\t") or line.startswith("Pdb "):
            start = index + 1
            break
    values: list[float] = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "Pdb", "SD")):
            continue
        fields = re.split(r"\s+", stripped)
        if not fields or not fields[0].endswith(".pdb"):
            continue          # banner or metadata line, never a per-run result
        try:
            values.append(float(fields[1]))
        except (IndexError, ValueError):
            continue
    return values


def _buildmodel_timeout(constructs: int, runs: int) -> int:
    """Scale the BuildModel wall-clock guard with the batch size.

    Measured on the bundled FoldX: ~4 s per construct-run.  A fixed guard
    either kills large legitimate batches or leaves a hung engine unchecked;
    scale at 60 s per construct-run, floor 30 min, cap 4 h.
    """
    return int(min(max(1800, constructs * runs * 60), 4 * 3600))


def score_with_foldx(
    binary: str,
    pdb_path: Path,
    mutations: list[tuple[str, int, str, str]],
    workdir: Path,
    runs: int = 5,
    already_repaired: bool = False,
) -> float | None:
    """Return FoldX BuildModel median delta-delta-G in kcal/mol when parseable."""
    temp = _engine_run_dir(workdir, "foldx")
    with nullcontext(temp):
        local_pdb = temp / "input.pdb"
        shutil.copy2(pdb_path, local_pdb)
        if already_repaired:
            repaired = local_pdb
        else:
            repair = _run_engine(
                [binary, "--command=RepairPDB", f"--pdb={local_pdb.name}"],
                temp,
                1800,
            )
            repaired = temp / "input_Repair.pdb"
            if repair is None or repair.returncode or not repaired.exists():
                return None
        mutation_text = ",".join(
            f"{wt}{chain}{number}{mutant}" for wt, number, chain, mutant in mutations
        ) + ";\n"
        (temp / "individual_list.txt").write_text(mutation_text, encoding="ascii")
        build = _run_engine(
            [
                binary,
                "--command=BuildModel",
                f"--pdb={repaired.name}",
                "--mutant-file=individual_list.txt",
                f"--numberOfRuns={runs}",
            ],
            temp,
            _buildmodel_timeout(len(mutation_groups), runs) if mutation_groups else 1800,
        )
        if build is None or build.returncode:
            return None
        values: list[float] = []
        for path in temp.glob("Dif_*.fxout"):
            values.extend(_parse_dif_fxout(path))
        return round(float(median(values)), 4) if values else None


def score_many_with_foldx(
    binary: str,
    pdb_path: Path,
    mutation_groups: list[list[tuple[str, int, str, str]]],
    workdir: Path,
    runs: int = 5,
    already_repaired: bool = False,
    max_workers: int | None = None,
) -> list[float | None]:
    """Score all constructs, splitting the batch across parallel FoldX processes.

    FoldX is single-core per invocation; on a multi-core host the batch is cut
    into near-equal chunks that run in isolated engine directories and are
    re-aggregated in input order.  Group indices in Dif row names are per-chunk,
    so each chunk is parsed independently -- no cross-chunk coupling.
    """
    if not mutation_groups:
        return []
    workers = max_workers or config.FOLDX_MAX_WORKERS
    workers = max(1, min(int(workers), len(mutation_groups)))
    if workers <= 1:
        return _run_foldx_batch(binary, pdb_path, mutation_groups, workdir, runs, already_repaired)
    chunks = _split_groups(mutation_groups, workers)
    results: list[list[float | None] | None] = [None] * len(chunks)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            index: pool.submit(
                _run_foldx_batch, binary, pdb_path, chunk, workdir, runs, already_repaired
            )
            for index, chunk in enumerate(chunks)
        }
        for index, future in futures.items():
            results[index] = future.result()
    return [item for chunk_result in results for item in (chunk_result or [])]


def _split_groups(groups: list[Any], n_chunks: int) -> list[list[Any]]:
    size, remainder = divmod(len(groups), n_chunks)
    chunks: list[list[Any]] = []
    start = 0
    for index in range(n_chunks):
        end = start + size + (1 if index < remainder else 0)
        if end > start:
            chunks.append(groups[start:end])
        start = end
    return chunks


def _run_foldx_batch(
    binary: str,
    pdb_path: Path,
    mutation_groups: list[list[tuple[str, int, str, str]]],
    workdir: Path,
    runs: int,
    already_repaired: bool,
) -> list[float | None]:
    """One BuildModel invocation over `mutation_groups` in a private run dir."""
    temp = _engine_run_dir(workdir, "foldx_batch")
    with nullcontext(temp):
        local_pdb = temp / "input.pdb"
        shutil.copy2(pdb_path, local_pdb)
        if already_repaired:
            repaired = local_pdb
        else:
            repair = _run_engine(
                [binary, "--command=RepairPDB", f"--pdb={local_pdb.name}"],
                temp,
                1800,
            )
            repaired = temp / "input_Repair.pdb"
            if repair is None or repair.returncode or not repaired.exists():
                return [None] * len(mutation_groups)
        lines = [
            ",".join(f"{wt}{chain}{number}{mutant}" for wt, number, chain, mutant in group) + ";"
            for group in mutation_groups
        ]
        (temp / "individual_list.txt").write_text("\n".join(lines) + "\n", encoding="ascii")
        build = _run_engine(
            [
                binary,
                "--command=BuildModel",
                f"--pdb={repaired.name}",
                "--mutant-file=individual_list.txt",
                f"--numberOfRuns={runs}",
            ],
            temp,
            _buildmodel_timeout(len(mutation_groups), runs),
        )
        if build is None or build.returncode:
            return [None] * len(mutation_groups)
        rows: list[tuple[str, float]] = []
        for path in temp.glob("Dif_*.fxout"):
            rows.extend(_parse_dif_rows(path))
        return _median_per_group(rows, len(mutation_groups), runs)


def _parse_dif_rows(path: Path) -> list[tuple[str, float]]:
    """Parse one Dif_*.fxout into (pdb_row_name, ddG) pairs.

    Row names encode the construct group: ``structure_Repair_<group>_<run>.pdb``.
    Grouping must use that name -- positional chunking silently shifts every
    median by one construct when a single run fails to appear.
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = 0
    for index, line in enumerate(lines):
        if line.startswith("Pdb\t") or line.startswith("Pdb "):
            start = index + 1
            break
    rows: list[tuple[str, float]] = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "Pdb", "SD")):
            continue
        fields = re.split(r"\s+", stripped)
        if not fields or not fields[0].endswith(".pdb"):
            continue          # banner or metadata line, never a per-run result
        try:
            rows.append((fields[0], float(fields[1])))
        except (IndexError, ValueError):
            continue
    return rows


def _parse_dif_fxout(path: Path) -> list[float]:
    """Backwards-compatible view of :func:`_parse_dif_rows` (values only)."""
    return [value for _, value in _parse_dif_rows(path)]


def _median_per_group(
    rows: list[tuple[str, float]], n_groups: int, runs: int
) -> list[float | None]:
    """Reduce parsed BuildModel rows to one median per construct.

    Rows are grouped by the group index encoded in the FoldX row name
    (``..._<group>_<run>.pdb``).  Only when names carry no group index do we
    fall back to input order.  A group with fewer rows than ``runs`` yields
    None for that construct: a partial run count is a review flag, never a
    silently mis-chunked median.
    """
    if not rows:
        return [None] * n_groups
    groups: dict[int, list[float]] = {}
    positional = True
    for name, value in rows:
        stem = name[:-4] if name.endswith(".pdb") else name
        parts = stem.rsplit("_", 2)
        if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
            positional = False
            groups.setdefault(int(parts[1]), []).append(value)
    if positional:
        # No group indices in names: chunk positionally only when the row
        # count matches exactly, else refuse to guess.
        if len(rows) != n_groups * runs:
            return [None] * n_groups
        groups = {
            index: [value for _, value in rows[index * runs : (index + 1) * runs]]
            for index in range(n_groups)
        }
    elif set(groups) == set(range(1, n_groups + 1)):
        # FoldX names its constructs 1..N; the caller indexes 0..N-1.
        groups = {index - 1: values for index, values in groups.items()}
    elif set(groups) != set(range(n_groups)):
        # Group indices neither 0-based nor 1-based over the expected range:
        # refuse to guess which construct each median belongs to.
        return [None] * n_groups
    output: list[float | None] = []
    for index in range(n_groups):
        values = groups.get(index, [])
        if len(values) >= runs:
            output.append(round(float(median(values)), 4))
        else:
            # Missing or incomplete run count for this construct.
            output.append(None)
    return output


def repair_with_foldx(binary: str, pdb_path: Path, workdir: Path) -> Path | None:
    """Repair one structure and share the result across ALL jobs on this host.

    RepairPDB costs minutes, so the repaired PDB is cached under
    ``<jobs-root>/foldx_repair_cache/<sha256-12>/`` -- keyed by the input
    structure hash, independent of any single job.  Re-running a target with
    different design parameters skips the repair entirely.
    """
    from .util import sha256_file

    digest = sha256_file(pdb_path)[:12]
    cache_dir = Path(workdir).parent / "foldx_repair_cache" / digest
    repaired = cache_dir / "structure_Repair.pdb"
    if repaired.exists():
        return repaired
    cache_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(pdb_path, cache_dir / "structure.pdb")
    result = _run_engine(
        [binary, "--command=RepairPDB", "--pdb=structure.pdb"],
        cache_dir,
        1800,
    )
    if result is None or result.returncode or not repaired.exists():
        return None
    return repaired


def engine_metadata(inventory: EngineInventory) -> dict[str, Any]:
    return {
        "primary": inventory.primary_name,
        "evoef2_available": bool(inventory.evoef2_path),
        "foldx_available": bool(inventory.foldx_path),
        "proxy_units": "dimensionless within-target ranking score; not kcal/mol",
        "foldx_thresholds_if_enabled": {
            "single_kcal_mol_exposed": config.FOLDX_GATE_EXPOSED,
            "single_kcal_mol_edge": config.FOLDX_GATE_EDGE,
            "single_kcal_mol_default": config.FOLDX_GATE_DEFAULT,
            "combination_kcal_mol": config.FOLDX_COMBINATION_GATE,
            "gate_policy": "exposure-linked; literature FoldX RMSE 0.85-1.1 kcal/mol (docs/改进方案.md)",
        },
    }
