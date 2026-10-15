import copy
import json

from ecd_epitope_scanner.joint_design_v4.reporting import project_display_analysis


def test_display_projection_preserves_cards_and_does_not_mutate_saved_evidence():
    model = {"version": "1.0.0", "weights_sha256": "weight-hash",
             "source_files_sha256": {"source.py": "source-hash"},
             "model_files_sha256": {"checkpoint.pt": "model-hash"},
             "git_tracked_modifications": [" M requirements.txt"]}
    card = {"candidate_id": "C1", "status": "PASSED", "gate_checks": [{"passed": True}],
            "site_details": [{"position": 8}], "stability_evidence": {"foldx": 0.1},
            "provenance": {"models": {"spurs_multi": {"provenance": model}}},
            "alternatives": [{"candidate_id": "C2", "status": "NEAR_THRESHOLD",
                              "provenance": {"source_files_sha256": {"analysis.json": "replay-source-hash"},
                                             "model_files_sha256_reference": {"file": "old-model.json", "count": 5}}}]}
    panel = {"sector_recommendations": [card], "selected_candidate_ids": ["C1"],
             "nearby_hypothesis_pairs": [["h1", "h2"]],
             "hierarchy_resolution": {"raw": {"near_pairs": [["h1", "h2"]],
                 "pairwise_distances": [0.1], "groups": [["h1", "h2"]]}}}
    analysis = {"joint_v4_ml": {"selected_panel": panel, "candidates": [{"candidate_id": "C1"}],
                               "surface_hypotheses": [{"hypothesis_id": "h1"}]},
                "partitioning": {"selected_panel": panel, "actual_regions": 1},
                "design_summary": {"sector_recommendations": [card], "epitope_mapping_plan": {"x": 1}},
                "epitope_mapping_plan": {"x": 1}, "residues": [{"uniprot_pos": 8}],
                "files": {"v4_manifest.json": "v4_manifest.json"}}
    saved = copy.deepcopy(analysis)
    display = project_display_analysis(analysis)
    assert analysis == saved
    assert project_display_analysis(display) == display
    v4 = display["joint_v4_ml"]
    assert "candidates" not in v4 and v4["candidates_reference"]["count"] == 1
    assert "selected_panel" not in display["partitioning"]
    shown = v4["selected_panel"]["sector_recommendations"][0]
    for field in ("gate_checks", "site_details", "stability_evidence"):
        assert shown[field] == card[field]
    assert shown["alternatives"][0]["candidate_id"] == "C2"
    assert shown["alternatives"][0]["provenance"]["model_files_sha256_reference"] == {
        "file": "old-model.json", "count": 5}
    shown_model = shown["provenance"]["models"]["spurs_multi"]["provenance"]
    assert shown_model["version"] == "1.0.0" and shown_model["weights_sha256"] == "weight-hash"
    ref = shown_model["source_files_sha256_reference"]
    assert ref["file"] == "analysis.json"
    original_value = analysis
    for field in ref["json_path"]:
        original_value = original_value[field]
    assert original_value == model["source_files_sha256"] and len(original_value) == ref["count"]
    references = []

    def check_identity_references(value):
        if isinstance(value, list):
            for item in value:
                check_identity_references(item)
        elif isinstance(value, dict):
            for key, item in value.items():
                if key in {"source_files_sha256_reference", "model_files_sha256_reference",
                           "git_tracked_modifications_reference"} and item.get("file") == "analysis.json":
                    original = analysis
                    for part in item["json_path"]:
                        original = original[part]
                    assert len(original) == item["count"]
                    references.append(original)
                else:
                    check_identity_references(item)

    check_identity_references(display)
    assert references == [model["source_files_sha256"], model["model_files_sha256"],
                          model["git_tracked_modifications"], {"analysis.json": "replay-source-hash"}]
    hierarchy = v4["selected_panel"]["hierarchy_resolution"]["raw"]
    assert hierarchy["groups"] == panel["hierarchy_resolution"]["raw"]["groups"]
    assert hierarchy["near_pairs_reference"]["json_path"] == ["resolution_by_level", "raw", "near_pairs"]
    json.dumps(display, allow_nan=False)


def test_legacy_display_projection_keeps_existing_payload():
    legacy = {"regions": [{"id": "R1"}], "designs": [], "residues": []}
    assert project_display_analysis(legacy) == legacy


def test_externalization_without_hierarchy_needs_no_hierarchy_artifact(tmp_path):
    from ecd_epitope_scanner.joint_design_v4.reporting import externalize_checkpoint_pair_diagnostics
    analysis = {"joint_v4_ml": {"selected_panel": {"hierarchy_resolution": {}}}}
    saved = copy.deepcopy(analysis)
    externalize_checkpoint_pair_diagnostics(analysis, tmp_path)
    assert analysis == saved
