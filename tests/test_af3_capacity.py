from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_2048_token_full_matrices_stay_below_one_gibibyte(tmp_path):
    child = r'''
import ctypes
import json
import os
import sys
from pathlib import Path

from ecd_epitope_scanner.af3_preprocessing import import_af3_directory, select_af3_model
from ecd_epitope_scanner.models import ProteinContext

root = Path(sys.argv[1])
source = root / "source" / "seed-1_sample-0"
source.mkdir(parents=True)
cif = source / "model.cif"
columns = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
    "auth_atom_id", "pdbx_PDB_model_num"]
with cif.open("w", encoding="ascii", newline="\n") as handle:
    handle.write("data_af3_capacity\n#\nloop_\n")
    for column in columns:
        handle.write(f"_atom_site.{column}\n")
    for pos in range(1, 2049):
        handle.write(f"ATOM {pos} C CA . ALA A 1 {pos} ? {pos * 3.8:.3f} 0.000 0.000 1.00 90.00 {pos} ALA A CA 1\n")
    handle.write("#\n")

def write_matrix(handle, size):
    handle.write("[")
    row = "[" + ",".join(["0"] * size) + "]"
    for index in range(size):
        if index:
            handle.write(",")
        handle.write(row)
    handle.write("]")

confidence = source / "full_data.json"
with confidence.open("w", encoding="ascii", newline="") as handle:
    handle.write('{"atom_plddts":[' + ",".join(["90"] * 2048) + '],')
    handle.write('"token_chain_ids":[' + ",".join(['"A"'] * 2048) + '],')
    handle.write('"token_res_ids":[' + ",".join(map(str, range(1, 2049))) + '],"pae":')
    write_matrix(handle, 2048)
    handle.write(',"contact_probs":')
    write_matrix(handle, 2048)
    handle.write("}")

manifest = import_af3_directory(root / "source", root / "frozen")
protein = ProteinContext("TEST", "capacity fixture", "test", "A" * 2048, 1, 2048, 2048, "unspecified")
selection = select_af3_model(manifest, protein, output_dir=root / "selected")
data = selection.confidence
cells = {}
for name in ("pae", "contact_probs"):
    matrix = data[name]
    assert len(matrix) == 2048 and all(len(row) == 2048 for row in matrix)
    cells[name] = sum(len(row) for row in matrix)
assert len(data["token_chain_ids"]) == len(data["token_res_ids"]) == 2048

if os.name == "nt":
    from ctypes import wintypes
    class PROCESS_MEMORY_COUNTERS_EX(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
            ("PrivateUsage", ctypes.c_size_t)]
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.GetCurrentProcess.restype = wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PROCESS_MEMORY_COUNTERS_EX), wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    counters = PROCESS_MEMORY_COUNTERS_EX()
    counters.cb = ctypes.sizeof(counters)
    if not psapi.GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    peak_rss = int(counters.PeakWorkingSetSize)
else:
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_rss = int(peak if sys.platform == "darwin" else peak * 1024)
print(json.dumps({"peak_rss_bytes": peak_rss, "token_count": 2048, "target_length": protein.ecd_length, "matrix_cells": cells,
    "model_count": len(json.loads(manifest.read_text())["models"]), "status": selection.report["status"]}))
'''
    completed = subprocess.run(
        [sys.executable, "-c", child, str(tmp_path)], cwd=ROOT.parent,
        capture_output=True, text=True, timeout=240,
    )
    assert completed.returncode == 0, completed.stderr
    result = json.loads(completed.stdout.strip().splitlines()[-1])
    assert result["token_count"] == 2048
    assert result["target_length"] == 2048
    assert result["matrix_cells"] == {"pae": 2048**2, "contact_probs": 2048**2}
    assert result["model_count"] == 1
    assert result["peak_rss_bytes"] < 1024**3
    receipt = ROOT / "test_output" / "af3_preprocessing_20261021" / "capacity_receipt.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(result | {"limit_bytes": 1024**3}, indent=2) + "\n", encoding="utf-8")


def test_token_and_confidence_size_gates_reject_before_large_work(tmp_path, monkeypatch):
    from ecd_epitope_scanner import af3_preprocessing as af3
    from ecd_epitope_scanner.util import PipelineError

    with pytest.raises(PipelineError, match="2048 tokens"):
        af3._pae({"token_chain_ids": ["A"] * 2049, "token_res_ids": list(range(2049))}, object(), 2048)

    source = tmp_path / "source" / "seed-1_sample-0"
    source.mkdir(parents=True)
    (source / "model.cif").write_text("placeholder", encoding="ascii")
    (source / "full_data.json").write_text('{"atom_plddts":[]}', encoding="ascii")
    monkeypatch.setattr(af3, "_MAX_CONFIDENCE_BYTES", 4)
    with pytest.raises(PipelineError, match="64 MiB"):
        af3.import_af3_directory(tmp_path / "source", tmp_path / "frozen")
    assert not (tmp_path / "frozen" / "af3_inputs").exists()
