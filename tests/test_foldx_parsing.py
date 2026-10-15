"""Regression tests for the FoldX Dif_*.fxout parser and median chunking.

The fixtures below are VERBATIM captures from the bundled FoldX binary
(``vendor/foldx/foldx_1_20270131.exe``, PD-L1/8AOM ECD repaired structure).
They reproduce the production bug where the banner line ``FoldX 5.1 (2011)``
was parsed as delta-delta-G = 5.1 and poisoned the first construct's median.
"""

from __future__ import annotations

from pathlib import Path

from ecd_epitope_scanner.stability import _median_per_group, _parse_dif_fxout


def test_buildmodel_guard_allows_slow_pdl1_batches_without_changing_runs():
    from ecd_epitope_scanner import config
    from ecd_epitope_scanner.stability import _buildmodel_timeout
    assert _buildmodel_timeout(5, 5) == 3600
    assert _buildmodel_timeout(6, 5) == 3600
    assert _buildmodel_timeout(8, 5) == 4800
    assert _buildmodel_timeout(1000, 5) == 14400
    settings = config.JOINT_V4_ML_PARAMETERS['foldx_confirmation']
    spurs = config.JOINT_V4_ML_PARAMETERS['spurs_multi']
    assert settings['timeout_seconds'] > _buildmodel_timeout(6, 5)
    assert settings['number_of_runs'] == 5
    assert settings['combination_gate_kcal_mol'] == 2.0
    assert settings['max_candidates'] == 60
    assert spurs['preferred_max_ddg'] == 0.5
    assert spurs['review_max_ddg'] == 1.0
    assert spurs['disagreement_gate_kcal_mol'] == 1.5


def test_native_timeout_is_recorded_and_never_returned_as_success(tmp_path, monkeypatch):
    import json
    import subprocess
    from ecd_epitope_scanner import stability
    def expire(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs['timeout'])
    monkeypatch.setattr(stability.subprocess, 'run', expire)
    assert stability._run_engine(['foldx.exe', '--command=BuildModel'], tmp_path, 3600) is None
    receipt = json.loads((tmp_path / 'engine_execution.json').read_text())
    assert receipt['status'] == 'TIMEOUT'
    assert receipt['timeout_seconds'] == 3600
    assert receipt['returncode'] is None

# Single construct, numberOfRuns=2 (captured 2026-09-03).
SINGLE_RUN_DIF = """FoldX 5.1 (2011)
by the FoldX Consortium
Jesper Borg, Frederic Rousseau, Joost Schymkowitz,
Luis Serrano and Francois Stricher
------------------------------------------------------

PDB file analysed: batch
Output type: BuildModel
Pdb\ttotal energy\tBackbone Hbond\tSidechain Hbond\tVan der Waals\tElectrostatics\tSolvation Polar\tSolvation Hydrophobic\tVan der Waals clashes\tentropy sidechain\tentropy mainchain\tsloop_entropy\tmloop_entropy\tcis_bond\ttorsional clash\tbackbone clash\thelix dipole\twater bridge\tdisulfide\telectrostatic kon\tpartial covalent bonds\tenergy Ionisation\tEntropy Complex
structure_Repair_1_0.pdb\t-0.597033\t-0.617447\t-0.617447\t0.155845\t-0.0174164\t-0.273664\t0.199517\t1.77636e-15\t0.489639\t0.0407704\t0\t0\t0\t0.0431711\t0\t0\t0\t0\t0\t0\t0\t0
structure_Repair_1_1.pdb\t-0.780479\t-0.690573\t-0.690573\t0.183328\t-0.0109187\t-0.300078\t0.255825\t0\t0.48852\t0.0304509\t0\t0\t0\t-0.046459\t0\t0\t0\t0\t0\t0\t0\t0
"""

# Five constructs, numberOfRuns=2 (captured 2026-09-03, values only; the
# intervening energy-term columns are trimmed to keep the fixture readable --
# parsing must not depend on them).
BATCH_FIRST_COLUMN = [
    "structure_Repair_1_0.pdb\t2.53489",
    "structure_Repair_1_1.pdb\t2.86757",
    "structure_Repair_2_0.pdb\t1.4499",
    "structure_Repair_2_1.pdb\t1.70996",
    "structure_Repair_3_0.pdb\t2.19008",
    "structure_Repair_3_1.pdb\t2.20384",
    "structure_Repair_4_0.pdb\t-0.545114",
    "structure_Repair_4_1.pdb\t-0.722922",
    "structure_Repair_5_0.pdb\t-0.165264",
    "structure_Repair_5_1.pdb\t-0.161045",
]


def _write_dif(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "Dif_structure.fxout"
    path.write_text(body, encoding="utf-8")
    return path


def _rows_from_fixture(rows: list[str]) -> list[tuple[str, float]]:
    return [(line.split("\t")[0], float(line.split("\t")[1])) for line in rows]


def test_parse_ignores_foldx_banner(tmp_path: Path) -> None:
    path = _write_dif(tmp_path, SINGLE_RUN_DIF)
    values = _parse_dif_fxout(path)
    assert values == [-0.597033, -0.780479]
    assert 5.1 not in values, "banner line 'FoldX 5.1 (2011)' leaked into scores"


def test_parse_expects_two_real_rows(tmp_path: Path) -> None:
    path = _write_dif(tmp_path, SINGLE_RUN_DIF)
    values = _parse_dif_fxout(path)
    assert len(values) == 2


def test_batch_values_in_input_order(tmp_path: Path) -> None:
    body = (
        "FoldX 5.1 (2011)\n"
        "Pdb\ttotal energy\tBackbone Hbond\n"
        + "\n".join(BATCH_FIRST_COLUMN)
        + "\n"
    )
    values = _parse_dif_fxout(_write_dif(tmp_path, body))
    assert values == [
        2.53489, 2.86757, 1.4499, 1.70996, 2.19008,
        2.20384, -0.545114, -0.722922, -0.165264, -0.161045,
    ]


def test_median_grouping_uses_row_names() -> None:
    rows = _rows_from_fixture(BATCH_FIRST_COLUMN)
    medians = _median_per_group(rows, n_groups=5, runs=2)
    assert medians == [
        round((2.53489 + 2.86757) / 2, 4),
        round((1.4499 + 1.70996) / 2, 4),
        round((2.19008 + 2.20384) / 2, 4),
        round((-0.545114 + -0.722922) / 2, 4),
        round((-0.165264 + -0.161045) / 2, 4),
    ]


def test_median_grouping_survives_a_missing_run() -> None:
    """One failed run must cost its own construct only, never shift medians.

    The old positional chunking moved every later median by one construct when
    a single run row was absent; name-based grouping isolates the damage.
    """
    rows = _rows_from_fixture([row for row in BATCH_FIRST_COLUMN if "structure_Repair_3_1" not in row])
    medians = _median_per_group(rows, n_groups=5, runs=2)
    assert medians[2] is None, "incomplete run count must flag the construct, not guess"
    assert medians[0] == round((2.53489 + 2.86757) / 2, 4)
    assert medians[1] == round((1.4499 + 1.70996) / 2, 4)
    assert medians[3] == round((-0.545114 + -0.722922) / 2, 4)
    assert medians[4] == round((-0.165264 + -0.161045) / 2, 4)


def test_median_chunking_returns_none_when_short() -> None:
    # Fewer rows than constructs and no group names -> no usable batch.
    assert _median_per_group([(f"a_{i}_0.pdb", 1.0) for i in range(1, 3)], n_groups=5, runs=2) == [None] * 5


def test_median_chunking_positional_requires_exact_count() -> None:
    # Names without a group index fall back to positional chunking, but only
    # when the row count matches exactly; otherwise refuse to guess.
    rows = [("model_0.pdb", 1.0), ("model_1.pdb", 2.0), ("model_2.pdb", 3.0)]
    assert _median_per_group(rows, n_groups=2, runs=2) == [None, None]
    exact = [("model_0.pdb", 1.0), ("model_1.pdb", 2.0)]
    assert _median_per_group(exact, n_groups=1, runs=2) == [1.5]
