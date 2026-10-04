from __future__ import annotations

import hashlib
import importlib.util
import io
import sys
import types
import zipfile
from pathlib import Path

import numpy as np
from ecd_epitope_scanner.experimental_grasp import FastaChain, build_constraint_package

FIXTURE = Path(__file__).parent / "fixtures" / "grasp" / "generate_restr.py"
EXPECTED_SHA256 = "dda53ab36cd88f639ec5ead73920dc71072ddc4c28e7b2bdbd77820208be0601"


def _load_official_parser():
    digest = hashlib.sha256(FIXTURE.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA256, "Pinned GRASP-JAX parser fixture hash changed"
    try:
        import absl.logging  # noqa: F401
    except ModuleNotFoundError:
        absl = types.ModuleType("absl")
        logging = types.ModuleType("absl.logging")
        logging.INFO = 20
        logging.info = lambda *args, **kwargs: None
        logging.set_verbosity = lambda *args, **kwargs: None
        absl.logging = logging
        sys.modules["absl"] = absl
        sys.modules["absl.logging"] = logging
    spec = importlib.util.spec_from_file_location("grasp_jax_generate_restr_fixture", FIXTURE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_serializer_restraints_round_trip_through_pinned_grasp_jax_parser(tmp_path):
    parser_module = _load_official_parser()
    fasta = tmp_path / "complex.fasta"
    package, _ = build_constraint_package(
        chains=[FastaChain(1, "ECD", "ACD"), FastaChain(2, "VH", "EFG"), FastaChain(3, "VL", "HIK")],
        experimental_ir=[{
            "constraint": "1-2-C", "evidence_id": "effect-1",
            "reviewer": "reviewer", "rationale": "single-site effect",
        }],
        predicted_ir=[{
            "constraint": "3-3-K", "evidence_id": "paraplume-1",
            "reviewer": "reviewer", "rationale": "reviewed prediction",
        }],
        manual_rpr=[{
            "constraint": "1-2-C, 2-3-G, 8.0, 0.05", "evidence_id": "complex-1",
            "reviewer": "reviewer", "rationale": "validated pair evidence",
        }],
        merge_sources=["experimental_ir", "predicted_ir", "manual_rpr"],
        manifest_context={"limitations": [], "provenance": {"fixture": "official-parser-roundtrip"}},
    )
    with zipfile.ZipFile(io.BytesIO(package)) as archive:
        fasta.write_bytes(archive.read("constraints.fasta"))
        restraints = archive.read("grasp_constraints.txt").decode("utf-8").splitlines()
    parser = parser_module.Restraints(str(fasta))

    initial = parser.textlist_to_arraydict(restraints)
    serialized = parser.arraydict_to_textlist(initial)
    round_tripped = parser.textlist_to_arraydict(serialized)

    assert {key for key, value in initial.items() if key.endswith("mask") and np.any(value)} == {
        "interface_mask",
        "sbr_mask",
    }
    for key in initial:
        np.testing.assert_allclose(round_tripped[key], initial[key])
    assert any(
        len(line.split(",")) in {3, 4}
        and line.replace(" ", "").startswith("1-2-C,2-3-G,")
        for line in serialized
    )
