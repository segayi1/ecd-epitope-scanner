from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from ecd_epitope_scanner.af3_preprocessing import (
    first_af3_structure,
    import_af3_directory,
    preserve_af3_manifest,
    select_af3_model,
)
from ecd_epitope_scanner.models import FeatureRange, ProteinContext
from ecd_epitope_scanner.util import PipelineError


def _write_model(path: Path, bfactor: float = 90.0, scale: float = 1.0) -> list[float]:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
        "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
        "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
        "auth_atom_id", "pdbx_PDB_model_num",
    ]
    lines = ["data_af3", "#", "loop_", *[f"_atom_site.{field}" for field in fields]]
    scores = []
    serial = 1
    for residue in range(1, 13):
        x = ((residue - 1) % 4) * 3.8 * scale
        y = ((residue - 1) // 4) * 3.8 * scale
        for atom, element, dx, dy, dz in (("N", "N", -1.2, 0.0, 0.0), ("CA", "C", 0.0, 0.0, 0.0), ("C", "C", 1.3, 0.0, 0.0), ("O", "O", 1.9, 0.8, 0.0)):
            lines.append(f"ATOM {serial} {element} {atom} . ALA A 1 {residue} ? {x+dx*scale:.3f} {y+dy*scale:.3f} {dz*scale:.3f} 1.00 {bfactor:.2f} {residue} ALA A {atom} 1")
            serial += 1
            scores.append(bfactor)
    lines.append("#")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return scores


def _write_pair(folder: Path, index: int, *, plddt: float = 90.0, ptm: float = 0.5, pae: bool = True, pae_value: float = 5.0, scale: float = 1.0) -> tuple[Path, Path]:
    cif = folder / f"job_model_{index}.cif"
    atom_scores = _write_model(cif, plddt, scale)
    n = 12
    confidence = {"atom_plddts": atom_scores}
    if pae:
        confidence.update({"pae": [[0.0 if i == j else pae_value for j in range(n)] for i in range(n)], "token_chain_ids": ["A"] * n, "token_res_ids": list(range(1, n + 1))})
    conf = folder / f"job_full_data_{index}.json"
    conf.write_text(json.dumps(confidence), encoding="utf-8")
    summary = folder / f"job_summary_confidences_{index}.json"
    summary.write_text(json.dumps({"ptm": ptm, "chain_ptm": ptm, "ranking_score": ptm, "has_clash": False, "mean_plddt": plddt, "coverage": 1.0}), encoding="utf-8")
    return cif, conf


def _protein() -> ProteinContext:
    return ProteinContext("TEST", "test", "test", "A" * 12, 1, 12, 12, "type I")


def test_import_dedupes_root_alias_but_retains_independent_samples(tmp_path: Path):
    source = tmp_path / "source"
    root_cif, root_conf = _write_pair(source, 0, plddt=90)
    sample0_cif, sample0_conf = _write_pair(source / "seed-1_sample-0", 0, plddt=90)
    sample1_cif, sample1_conf = _write_pair(source / "seed-1_sample-1", 0, plddt=90)
    # A root alias duplicates sample 0 byte-for-byte; sample 1 remains an independent draw.
    shutil.copyfile(sample0_cif, root_cif)
    shutil.copyfile(sample0_conf, root_conf)
    manifest_path = import_af3_directory(source, tmp_path / "frozen")
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert len(data["models"]) == 2
    assert {row["sample_id"] for row in data["models"]} == {"0", "1"}
    assert all(row["system_fingerprint"] and row["context_status"] == "UNKNOWN"
               for row in data["models"])
    assert first_af3_structure(manifest_path).is_file()


def test_root_model_alias_is_deduped_even_when_sample_index_differs(tmp_path: Path):
    source = tmp_path / "source"
    root_cif, root_conf = _write_pair(source, 0)
    root_summary = source / "job_summary_confidences_0.json"
    nested = source / "seed-1_sample-3"
    nested.mkdir(parents=True)
    shutil.copyfile(root_cif, nested / "seed-1_sample-3_model.cif")
    shutil.copyfile(root_conf, nested / "seed-1_sample-3_full_data.json")
    shutil.copyfile(root_summary, nested / "seed-1_sample-3_summary_confidences.json")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert len(data["models"]) == 1
    assert data["models"][0]["sample_id"] == "3"
    assert len(data["aliases"]) == 1


def test_select_uses_atom_identity_validated_confidence_and_equal_sample_weights(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, plddt=90, ptm=0.55)
    _write_pair(source / "seed-1_sample-1", 0, plddt=91, ptm=0.95)
    manifest_path = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest_path, _protein(), output_dir=tmp_path / "selected")
    assert result.report["selected_model_id"] in {row["id"] for row in json.loads(manifest_path.read_text())["models"]}
    assert result.report["status"] == "ENSEMBLE_CORE"
    assert len(result.report["Rcore"]) == 12
    assert len(result.report["Ecore"]) >= 20
    assert sum(result.report["weights"].values()) == pytest.approx(1.0)
    assert sorted(result.report["weights"].values()) == pytest.approx([0.5, 0.5])
    assert result.confidence["_source_sha256"] == result.report["selected_sha256"]["confidence"]
    assert result.confidence["_source_structure_sha256"] == result.report["selected_sha256"]["structure"]
    assert len(result.confidence["_validated_atom_plddts"]) == 12
    # The delivered confidence JSON stays byte-identical to its frozen source.
    selected_row = next(row for row in json.loads(manifest_path.read_text())["models"] if row["id"] == result.report["selected_model_id"])
    assert result.confidence_path.read_bytes() == (manifest_path.parent.parent / selected_row["confidence"]).read_bytes()


def test_pae_maps_by_unambiguous_chain_order_when_residue_ids_are_absent(tmp_path: Path):
    source = tmp_path / "source"
    _, confidence_path = _write_pair(source / "seed-1_sample-0", 0)
    confidence = json.loads(confidence_path.read_text(encoding="utf-8"))
    confidence.pop("token_res_ids")
    confidence_path.write_text(json.dumps(confidence), encoding="utf-8")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert result.report["model_diagnostics"][0]["pae_status"] == "mapped"
    assert result.report["Ecore"]


def test_atom_confidence_must_agree_with_cif_bfactor(tmp_path: Path):
    source = tmp_path / "source"
    _, confidence = _write_pair(source / "seed-1_sample-0", 0, plddt=90)
    data = json.loads(confidence.read_text())
    data["atom_plddts"][0] = 89.0
    confidence.write_text(json.dumps(data), encoding="utf-8")
    manifest_path = import_af3_directory(source, tmp_path / "frozen")
    with pytest.raises(PipelineError, match="B factor"):
        select_af3_model(manifest_path, _protein(), output_dir=tmp_path / "selected")


def test_verified_residue_plddt_exactly_matches_direct_aggregation(tmp_path: Path):
    from ecd_epitope_scanner.af3_preprocessing import _raw_atom_map
    from ecd_epitope_scanner.structure_analysis import _af3_atom_plddt, load_structure
    cif, confidence_path = _write_pair(tmp_path / "source", 0)
    data = json.loads(confidence_path.read_text())
    data["atom_plddts"][:4] = [90.0, 90.01, 90.05, 90.08]
    parsed = load_structure(cif)
    verified, _ = _raw_atom_map(parsed, data, "fixture")
    direct = _af3_atom_plddt(parsed, data)
    for residue in parsed.model["A"]:
        assert verified[("A", residue.id[1], "")] == direct[id(residue)]


def test_bfactor_confidence_difference_of_point_one_is_within_tolerance(tmp_path: Path):
    source = tmp_path / "source"
    _, confidence_path = _write_pair(source / "seed-1_sample-0", 0, plddt=63.5)
    confidence = json.loads(confidence_path.read_text(encoding="utf-8"))
    confidence["atom_plddts"][0] = 63.6
    confidence_path.write_text(json.dumps(confidence), encoding="utf-8")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert result.report["models"][0]["eligible"] is True


def test_missing_pae_is_diagnostic_and_falls_back_without_core_gate(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, pae=False)
    manifest_path = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest_path, _protein(), output_dir=tmp_path / "selected")
    assert result.report["status"] == "INSUFFICIENT_ENSEMBLE_FALLBACK"
    assert result.report["Ecore"] == []
    assert result.report["models"][0]["eligible"] is True


def test_preserve_copies_only_manifest_members_and_keeps_group_identity(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-2_sample-4", 0)
    manifest_path = import_af3_directory(source, tmp_path / "frozen")
    preserved = preserve_af3_manifest(manifest_path, tmp_path / "replay")
    original = json.loads(manifest_path.read_text())
    copied = json.loads(preserved.read_text())
    assert copied == original
    assert copied["models"][0]["seed_id"] == "2"
    assert (tmp_path / "replay" / copied["models"][0]["structure"]).is_file()


def test_preserve_is_idempotent_when_manifest_already_in_destination(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    preserved = preserve_af3_manifest(manifest, manifest.parent.parent)
    assert preserved == manifest


def test_import_refuses_to_overwrite_frozen_inputs(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0)
    import_af3_directory(source, tmp_path / "frozen")
    with pytest.raises(PipelineError, match="already exists"):
        import_af3_directory(source, tmp_path / "frozen")


def test_import_accepts_open_source_bare_seed_sample_layout(tmp_path: Path):
    sample = tmp_path / "source" / "seed-4_sample-2"
    cif = sample / "model.cif"
    atom_scores = _write_model(cif)
    (sample / "full_data.json").write_text(json.dumps({"atom_plddts": atom_scores}), encoding="utf-8")
    (sample / "summary_confidences.json").write_text(json.dumps({"ptm": 0.8}), encoding="utf-8")
    manifest = import_af3_directory(tmp_path / "source", tmp_path / "frozen")
    row = json.loads(manifest.read_text(encoding="utf-8"))["models"][0]
    assert row["seed_id"] == "4"
    assert row["sample_id"] == "2"
    assert row["source_name"] == "model.cif"


def test_import_rejects_oversized_summary_before_reading_json(tmp_path: Path):
    source = tmp_path / "source"
    cif, _ = _write_pair(source, 0)
    summary = cif.with_name("job_summary_confidences_0.json")
    with summary.open("wb") as handle:
        handle.truncate(64 * 1024**2 + 1)

    with pytest.raises(PipelineError, match="JSON file exceeds 64 MiB"):
        import_af3_directory(source, tmp_path / "frozen")


def test_import_accepts_prefixed_seed_sample_model_filename(tmp_path: Path):
    sample = tmp_path / "source" / "seed-123_sample-0"
    cif = sample / "job_seed-123_sample-0_model.cif"
    atom_scores = _write_model(cif)
    (sample / "job_seed-123_sample-0_confidences.json").write_text(json.dumps({"atom_plddts": atom_scores}), encoding="utf-8")
    manifest = import_af3_directory(tmp_path / "source", tmp_path / "frozen")
    row = json.loads(manifest.read_text(encoding="utf-8"))["models"][0]
    assert (row["job_id"], row["seed_id"], row["sample_id"]) == ("job", "123", "0")


def test_import_freezes_suffixed_job_metadata_and_uses_single_model_seed(tmp_path: Path):
    source = tmp_path / "source"
    sample = source / "seed-unknown_sample-0"
    _write_pair(sample, 0)
    metadata = source / "fold_2026_job_request.json"
    metadata.write_text(json.dumps([{"sequences": [{"protein": {"sequence": "A" * 12}}], "modelSeeds": [123], "dialect": "alphafold3", "version": 1}]), encoding="utf-8")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    data = json.loads(manifest.read_text(encoding="utf-8"))
    row = data["models"][0]
    assert row["seed_id"] == "123"
    assert row["metadata_files"]
    metadata_member = manifest.parent.parent / row["metadata_files"][0]["path"]
    assert metadata_member.read_bytes() == metadata.read_bytes()
    copied = preserve_af3_manifest(manifest, tmp_path / "replay")
    assert (copied.parent.parent / row["metadata_files"][0]["path"]).read_bytes() == metadata.read_bytes()


def test_import_skips_model_missing_confidence_and_records_exclusion(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0)
    missing = source / "seed-1_sample-1" / "job_model_0.cif"
    _write_model(missing)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    data = json.loads(manifest.read_text(encoding="utf-8"))
    assert len(data["models"]) == 1
    assert data["excluded_members"] == [{"source_name": "job_model_0.cif", "reason": "missing confidence sidecar", "sha256": data["excluded_members"][0]["sha256"]}]


def test_first_structure_skips_unparseable_hashed_member(tmp_path: Path):
    source = tmp_path / "source"
    bad_cif, _ = _write_pair(source / "seed-1_sample-0", 0)
    _write_pair(source / "seed-1_sample-1", 0)
    bad_cif.write_text("not a CIF\n", encoding="utf-8")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    chosen = first_af3_structure(manifest)
    assert "m001" in chosen.parts


def test_weights_are_equal_by_seed_then_sample_and_contacts_need_half_support(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, scale=1.0)
    _write_pair(source / "seed-1_sample-1", 0, scale=1.0)
    _write_pair(source / "seed-2_sample-0", 0, scale=2.0, pae_value=20.0)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert sorted(result.report["weights"].values()) == pytest.approx([0.25, 0.25, 0.5])
    edge = next(edge for edge in result.report["diagnostic_edges"] if (edge["residue_i"], edge["residue_j"]) == (1, 4))
    assert edge["distance_support"] == pytest.approx(2 / 3)
    assert any((candidate["residue_i"], candidate["residue_j"]) == (1, 4)
               for candidate in result.report["Ecore"])


def test_per_residue_diagnostics_keep_low_confidence_models_outside_q(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, plddt=90)
    _write_pair(source / "seed-1_sample-1", 0, plddt=60)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    first = result.report["per_residue_diagnostics"]["1"]
    assert len(result.report["Rdiag"]) == 12
    assert first["n_valid"] == 2
    assert first["n_low_confidence"] == 1
    assert first["low_confidence_models"] == ["m001"]


def test_diagnostic_contacts_and_domain_geometry_include_low_confidence_models(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, plddt=90, scale=1.5)
    _write_pair(source / "seed-1_sample-1", 0, plddt=90, scale=1.5)
    _write_pair(source / "seed-2_sample-0", 0, plddt=60, scale=0.5)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    protein = _protein()
    protein.features = [FeatureRange("domain", 1, 6, "D1"), FeatureRange("domain", 7, 12, "D2")]
    result = select_af3_model(manifest, protein, output_dir=tmp_path / "selected")
    pair = next(edge for edge in result.report["diagnostic_edges"]
                if (edge["residue_i"], edge["residue_j"]) == (1, 12))
    assert pair["distance_support"] == pytest.approx(1 / 3)
    assert pair in result.report["minority_contacts"]
    assert len(result.report["domain_geometry"]["domains"][0]["pairwise_rmsd"]) == 3
    assert len(result.report["domain_geometry"]["inter_domain_distances"][0]["by_model"]) == 3
    assert result.report["domain_pae"]["inter_domain"][0]["mean_symmetric_max_pae"] == pytest.approx(5.0)


def test_ensemble_selection_stays_within_point_zero_five_of_raw_medoid(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, ptm=0.50, scale=1.0)
    _write_pair(source / "seed-2_sample-0", 0, ptm=0.60, scale=1.2)
    _write_pair(source / "seed-3_sample-0", 0, ptm=0.99, scale=1.4)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert result.report["status"] == "ENSEMBLE_CORE"
    assert result.report["baselines"]["raw_medoid"] == "m001"
    assert result.report["selected_model_id"] == "m001"
    assert result.report["baselines"]["ptm_best"] == "m002"


def test_baselines_include_original_rank_one_even_when_confidence_is_rejected(tmp_path: Path):
    source = tmp_path / "source"
    _, bad_conf = _write_pair(source / "seed-1_sample-0", 0, plddt=90, ptm=0.95)
    bad_data = json.loads(bad_conf.read_text(encoding="utf-8"))
    bad_data["atom_plddts"][0] = 80.0
    bad_conf.write_text(json.dumps(bad_data), encoding="utf-8")
    _write_pair(source / "seed-1_sample-0", 1, plddt=90, ptm=0.80)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    baseline = result.report["baselines"]
    assert result.report["selected_model_id"] == "m001"
    assert baseline["original_rank1"] == "m000"
    assert baseline["original_rank1_ties"] == ["m000"]
    assert baseline["original_rank1_exclusion_reason"] == "AF3 CIF B factor disagrees with mapped pLDDT by more than 0.1"
    assert baseline["original_rank1_exclusions"] == [{"model_id": "m000", "reason": baseline["original_rank1_exclusion_reason"]}]
    assert baseline["ptm_best"] == "m000"
    assert baseline["plddt_best"] == "m001"
    assert result.report["models"][0]["eligible"] is False


def test_original_rank_one_ties_are_reported_with_deterministic_member(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0, ptm=0.8)
    _write_pair(source / "seed-2_sample-0", 0, ptm=0.8)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert result.report["baselines"]["original_rank1_ties"] == ["m000", "m001"]
    assert result.report["baselines"]["original_rank1"] == "m000"


def test_bare_server_filenames_and_distinct_ptm_ranking_baselines(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source, 0, ptm=0.7)
    _write_pair(source, 1, ptm=0.95)
    summaries = sorted(source.glob("job_summary_confidences_*.json"))
    first, second = (json.loads(path.read_text()) for path in summaries)
    first["ranking_score"], second["ranking_score"] = 1.2, -50
    for path, value in zip(summaries, (first, second), strict=True):
        path.write_text(json.dumps(value), encoding="utf-8")
    for index, source_path in enumerate(sorted(source.glob("job_model_*.cif"))):
        source_path.rename(source / source_path.name.replace("job_", ""))
    for index, source_path in enumerate(sorted(source.glob("job_full_data_*.json"))):
        source_path.rename(source / source_path.name.replace("job_", ""))
    for source_path in sorted(source.glob("job_summary_confidences_*.json")):
        source_path.rename(source / source_path.name.replace("job_", ""))
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert result.report["selected_model_id"] == "m001"
    assert result.report["baselines"]["original_rank1"] == "m000"
    assert result.report["baselines"]["original_rank1_value"] == pytest.approx(1.2)
    assert result.report["baselines"]["ptm_best"] == "m001"
    assert result.report["models"][0]["ptm"] == pytest.approx(0.7)
    assert result.report["models"][0]["ranking_score"] == pytest.approx(1.2)
    assert result.report["models"][1]["ranking_score"] == pytest.approx(-50)
    assert result.report["algorithm_version"] == "af3_preprocessing_v1"


def test_ptm_outside_unit_interval_is_not_reported_as_legal(tmp_path: Path):
    source = tmp_path / "source"
    _, _ = _write_pair(source / "seed-1_sample-0", 0, ptm=1.2)
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    assert result.report["models"][0]["ptm"] is None
    assert result.report["baselines"]["ptm_best"] is None


def test_import_rejects_symlink_and_confidence_size_limits(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0)
    outside = tmp_path / "outside.json"
    outside.write_text("{}", encoding="utf-8")
    link = source / "external.json"
    try:
        link.symlink_to(outside)
    except OSError as exc:
        pytest.skip(f"Symlink creation is unavailable: {exc}")
    with pytest.raises(PipelineError, match="symlink"):
        import_af3_directory(source, tmp_path / "frozen")


def test_malformed_prefixed_confidence_is_reported_for_its_model(tmp_path: Path):
    source = tmp_path / "source"
    good = source / "seed-1_sample-0"
    atom_scores = _write_model(good / "job_seed-1_sample-0_model.cif")
    (good / "job_seed-1_sample-0_confidences.json").write_text(json.dumps({"atom_plddts": atom_scores}), encoding="utf-8")
    bad = source / "seed-1_sample-1"
    _write_model(bad / "job_seed-1_sample-1_model.cif")
    (bad / "job_seed-1_sample-1_confidences.json").write_text("{broken", encoding="utf-8")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    result = select_af3_model(manifest, _protein(), output_dir=tmp_path / "selected")
    broken = next(row for row in result.report["models"] if row["model_id"] == "m001")
    assert broken["eligible"] is False
    assert "Invalid AF3 full_data JSON" in broken["exclusion_reasons"][0]


def test_import_excludes_cif_whose_protein_system_disagrees_with_request(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0)
    request = source / "job_request.json"
    request.write_text(json.dumps({"sequences": [{"protein": {"sequence": "G" * 12}}]}), encoding="utf-8")
    with pytest.raises(PipelineError, match="No AF3 models"):
        import_af3_directory(source, tmp_path / "frozen")


def test_import_rejects_mixed_observed_systems_within_one_job(tmp_path: Path):
    source = tmp_path / "source"
    _write_pair(source / "seed-1_sample-0", 0)
    second, _ = _write_pair(source / "seed-1_sample-1", 0)
    second.write_text(second.read_text(encoding="utf-8").replace("ALA", "GLY"), encoding="utf-8")
    with pytest.raises(PipelineError, match="mixed observed"):
        import_af3_directory(source, tmp_path / "frozen")


def test_manifest_rederives_metadata_and_group_identity_from_frozen_sources(tmp_path: Path):
    source = tmp_path / "source"
    sample = source / "seed-2_sample-4"
    _write_pair(sample, 0)
    request = source / "job_request.json"
    request.write_text(json.dumps({"sequences": [{"protein": {"sequence": "A" * 12}}], "modelSeeds": [2]}), encoding="utf-8")
    manifest = import_af3_directory(source, tmp_path / "frozen")
    data = json.loads(manifest.read_text(encoding="utf-8"))
    data["models"][0]["seed_id"] = "99"
    manifest.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(PipelineError, match="grouping disagrees"):
        first_af3_structure(manifest)
    data["models"][0]["seed_id"] = "2"
    data["models"][0]["metadata"]["version"] = 99
    manifest.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(PipelineError, match="metadata disagrees"):
        first_af3_structure(manifest)
