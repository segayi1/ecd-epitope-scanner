import json

from scripts.audit_p0_p3_batch import _load_pairs


def test_pair_audit_ignores_historical_attempt_manifests(tmp_path):
    for relative in (
        "paired_model_manifest.json",
        "paired_manifest_history/old/paired_model_manifest.json",
        "prior_attempts/old/paired_model_manifest.json",
    ):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"records": []}), encoding="utf-8")

    _, errors, active_manifest_count, *_ = _load_pairs(tmp_path, {})

    assert errors == []
    assert active_manifest_count == 1
