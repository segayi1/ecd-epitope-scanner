from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
from typing import Any

from .models import MutationChoice, ProteinContext, ResidueRecord


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
        return MutationChoice(record.uniprot_pos, wt, mutant, 99.0, 0.0, 1.0, False, ["WT-to-WT"])
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
    if risk > 1.45:
        reasons.append("high structural-risk proxy")
        allowed = False
    if record.interface:
        reasons.append("native interface; function-risk flag retained")
    if not reasons:
        reasons.append("passes deterministic chemistry and structural-risk constraints")
    return MutationChoice(
        position=record.uniprot_pos,
        wt=wt,
        mutant=mutant,
        stability_score=round(float(risk), 4),
        perturbation_score=round(float(perturbation), 4),
        stability_percentile=0.0,
        allowed=allowed,
        reasons=reasons,
    )


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
    with TemporaryDirectory(prefix="evoef2_", dir=workdir) as temp_name:
        temp = Path(temp_name)
        local_pdb = temp / "input.pdb"
        shutil.copy2(pdb_path, local_pdb)
        repair = subprocess.run(
            [binary, "--command=RepairStructure", f"--pdb={local_pdb.name}", "--num_of_runs=3"],
            cwd=temp,
            capture_output=True,
            text=True,
            timeout=300,
            check=False,
        )
        repaired = temp / "input_Repair.pdb"
        if repair.returncode != 0 or not repaired.exists():
            return None
        mutation_text = ",".join(f"{wt}{chain}{number}{mutant}" for wt, number, chain, mutant in mutations) + ";\n"
        (temp / "individual_list.txt").write_text(mutation_text, encoding="ascii")
        build = subprocess.run(
            [
                binary,
                "--command=BuildMutant",
                f"--pdb={repaired.name}",
                "--mutant_file=individual_list.txt",
                "--num_of_runs=10",
            ],
            cwd=temp,
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
        )
        mutant_files = sorted(temp.glob("*_Model_1.pdb"))
        if build.returncode != 0 or not mutant_files:
            return None

        def stability(path: Path) -> float | None:
            result = subprocess.run(
                [binary, "--command=ComputeStability", f"--pdb={path.name}"],
                cwd=temp,
                capture_output=True,
                text=True,
                timeout=300,
                check=False,
            )
            if result.returncode:
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


def score_with_foldx(
    binary: str,
    pdb_path: Path,
    mutations: list[tuple[str, int, str, str]],
    workdir: Path,
    runs: int = 5,
    already_repaired: bool = False,
) -> float | None:
    """Return FoldX BuildModel median delta-delta-G in kcal/mol when parseable."""
    with TemporaryDirectory(prefix="foldx_", dir=workdir) as temp_name:
        temp = Path(temp_name)
        local_pdb = temp / "input.pdb"
        shutil.copy2(pdb_path, local_pdb)
        if already_repaired:
            repaired = local_pdb
        else:
            repair = subprocess.run(
                [binary, "--command=RepairPDB", f"--pdb={local_pdb.name}"],
                cwd=temp,
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            repaired = temp / "input_Repair.pdb"
            if repair.returncode or not repaired.exists():
                return None
        mutation_text = ",".join(
            f"{wt}{chain}{number}{mutant}" for wt, number, chain, mutant in mutations
        ) + ";\n"
        (temp / "individual_list.txt").write_text(mutation_text, encoding="ascii")
        build = subprocess.run(
            [
                binary,
                "--command=BuildModel",
                f"--pdb={repaired.name}",
                "--mutant-file=individual_list.txt",
                f"--numberOfRuns={runs}",
            ],
            cwd=temp,
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
        )
        if build.returncode:
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
) -> list[float | None]:
    """Score all regional constructs in one BuildModel invocation.

    FoldX writes one row per run to Dif_*.fxout. Rows are kept in input order and
    reduced by median per construct, preserving the requested repeat count.
    """
    if not mutation_groups:
        return []
    with TemporaryDirectory(prefix="foldx_batch_", dir=workdir) as temp_name:
        temp = Path(temp_name)
        local_pdb = temp / "input.pdb"
        shutil.copy2(pdb_path, local_pdb)
        if already_repaired:
            repaired = local_pdb
        else:
            repair = subprocess.run(
                [binary, "--command=RepairPDB", f"--pdb={local_pdb.name}"],
                cwd=temp,
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            repaired = temp / "input_Repair.pdb"
            if repair.returncode or not repaired.exists():
                return [None] * len(mutation_groups)
        lines = [
            ",".join(f"{wt}{chain}{number}{mutant}" for wt, number, chain, mutant in group) + ";"
            for group in mutation_groups
        ]
        (temp / "individual_list.txt").write_text("\n".join(lines) + "\n", encoding="ascii")
        build = subprocess.run(
            [
                binary,
                "--command=BuildModel",
                f"--pdb={repaired.name}",
                "--mutant-file=individual_list.txt",
                f"--numberOfRuns={runs}",
            ],
            cwd=temp,
            capture_output=True,
            text=True,
            timeout=3600,
            check=False,
        )
        if build.returncode:
            return [None] * len(mutation_groups)
        values: list[float] = []
        for path in temp.glob("Dif_*.fxout"):
            values.extend(_parse_dif_fxout(path))
        if len(values) < len(mutation_groups):
            return [None] * len(mutation_groups)
        if len(values) >= len(mutation_groups) * runs:
            return [
                round(float(median(values[index * runs : (index + 1) * runs])), 4)
                for index in range(len(mutation_groups))
            ]
        return [round(float(value), 4) for value in values[: len(mutation_groups)]]


def repair_with_foldx(binary: str, pdb_path: Path, workdir: Path) -> Path | None:
    """Repair one structure once per job and cache the FoldX-ready PDB."""
    cache_dir = workdir / "foldx_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    repaired = cache_dir / "structure_Repair.pdb"
    if repaired.exists():
        return repaired
    shutil.copy2(pdb_path, cache_dir / "structure.pdb")
    result = subprocess.run(
        [binary, "--command=RepairPDB", "--pdb=structure.pdb"],
        cwd=cache_dir,
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )
    if result.returncode or not repaired.exists():
        return None
    return repaired


def engine_metadata(inventory: EngineInventory) -> dict[str, Any]:
    return {
        "primary": inventory.primary_name,
        "evoef2_available": bool(inventory.evoef2_path),
        "foldx_available": bool(inventory.foldx_path),
        "proxy_units": "dimensionless within-target ranking score; not kcal/mol",
        "foldx_thresholds_if_enabled": {"single_kcal_mol": 0.46, "combination_kcal_mol": 0.46},
    }
