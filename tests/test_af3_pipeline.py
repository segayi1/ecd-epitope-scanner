import json
from pathlib import Path

import pytest

from ecd_epitope_scanner.pipeline import PipelineRequest, _preserve_uploaded_inputs, run_pipeline
from ecd_epitope_scanner.structure_analysis import _af3_atom_plddt, load_structure
from ecd_epitope_scanner.util import PipelineError, sha256_file


ROOT = Path(__file__).resolve().parents[1]


def test_verified_atom_map_uses_residue_identity_and_refuses_changed_coordinates(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text("ATOM      1  CA  ALA A  78       0.000   0.000   0.000  1.00 20.00           C  \nEND\n")
    parsed = load_structure(source)
    confidence = {"_source_structure_sha256": sha256_file(source),
                  "_validated_atom_plddts": {("A", 78, ""): 93.0}, "atom_plddts": [1.0]}
    assert _af3_atom_plddt(parsed, confidence) == {id(parsed.model["A"][78]): 93.0}
    source.write_text(source.read_text().replace("0.000", "1.000", 1))
    with pytest.raises(PipelineError, match="different structure"):
        _af3_atom_plddt(parsed, confidence)


def test_manifest_request_requires_predicted_input_semantics():
    with pytest.raises(ValueError, match="input_type=af3_upload"):
        PipelineRequest(input_type="structure", af3_manifest_path="manifest.json")


def test_display_keeps_selection_and_references_full_diagnostics_without_mutating_evidence():
    from ecd_epitope_scanner.joint_design_v4.reporting import project_display_analysis
    result = {"joint_v4_ml": {"selected_panel": {}}, "methods": {"af3_preprocessing": {
        "selected_model_id": "seed-1/sample-2", "models": [{"model_id": "seed-1/sample-2"}],
        "per_residue_diagnostics": {"2": {"low_confidence_models": ["sample-1"]}},
        "Ecore": [[1, 5]]}}, "source": {"af3_preprocessing": {"report": "af3_selection_report.json"}}}
    original = json.dumps(result, sort_keys=True)
    compact = project_display_analysis(result)["methods"]["af3_preprocessing"]
    assert compact["selected_model_id"] == "seed-1/sample-2"
    assert compact["per_residue_diagnostics_reference"]["file"] == "af3_selection_report.json"
    assert "per_residue_diagnostics" not in compact
    assert json.dumps(result, sort_keys=True) == original


def test_archive_includes_frozen_manifest_identity(tmp_path):
    manifest = tmp_path / "af3_inputs" / "manifest.json"
    manifest.parent.mkdir()
    manifest.write_text('{"models": []}')
    request = PipelineRequest(input_type="af3_upload", af3_manifest_path=str(manifest))
    identity = _preserve_uploaded_inputs(request, tmp_path)
    assert identity["files"]["af3_manifest"]["artifact"] == "af3_inputs/manifest.json"
    assert identity["files"]["af3_manifest"]["sha256"] == sha256_file(manifest)


@pytest.mark.integration
def test_real_folder_and_selected_pair_have_identical_scientific_outputs(tmp_path, monkeypatch):
    from ecd_epitope_scanner.af3_preprocessing import import_af3_directory
    from ecd_epitope_scanner.scripts.audit_run_bundle import audit_run_bundle
    from ecd_epitope_scanner import sequences

    # Freeze the unrelated stochastic CDS step while testing both structural entry points.
    monkeypatch.setattr(sequences, "optimize_human_coding_fragment",
        lambda sequence: (sequences.reverse_translate(sequence), "fixed codons for AF3 equivalence test"))

    batch = ROOT / "workspace_data/af3_binder_batch_20261017/batch_manifest.json"
    if not batch.is_file():
        pytest.skip("Saved AF3 source manifest is not available")
    item = next(row for row in json.loads(batch.read_text(encoding="utf-8"))["preflight_items"]
                if row["binder"] == "AIB203")
    if not Path(item["source_directory"]).is_dir():
        pytest.skip("Original AF3 model folder is not available")
    folder = tmp_path / "folder"
    manifest = import_af3_directory(Path(item["source_directory"]), folder)
    shared = dict(input_type="af3_upload", analysis_target="binder", chain="A", ml_mode="preview",
                  metadata={"disable_external_engines": True})
    result = run_pipeline(PipelineRequest(**shared, af3_manifest_path=str(manifest), output_dir=str(folder)))
    direct = tmp_path / "direct"
    pair = run_pipeline(PipelineRequest(**shared, structure_path=result["request"]["structure_path"],
                        confidence_path=result["request"]["confidence_path"], output_dir=str(direct)))
    for name in ("structure_for_energy.pdb", "residue_mapping.csv", "residue_reliability.csv",
                 "surface_mask.csv", "surface_display_regions_v41.json"):
        assert (folder / name).read_bytes() == (direct / name).read_bytes(), name
    assert result["residues"] == pair["residues"]
    assert result["regions"] == pair["regions"]
    from Bio.Seq import Seq
    def scientific_designs(designs):
        normalized = []
        for design in designs:
            item = dict(design)
            coding = item.pop("mutant_coding_sequence", None)
            if coding:
                # Synonymous codon optimization is stochastic in the existing common flow.
                item["translated_coding_sequence"] = str(Seq(coding).translate())
            normalized.append(item)
        return normalized
    assert scientific_designs(result["designs"]) == scientific_designs(pair["designs"])
    assert result["request"]["input_type"] == "af3_upload"
    assert result["methods"]["af3_preprocessing"]["selected_model_id"]
    assert audit_run_bundle(folder)["status"] == "PASS"
    assert audit_run_bundle(direct)["status"] == "PASS"
    assert "af3_inputs/manifest.json" in result["files"]
