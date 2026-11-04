import json
from dataclasses import asdict

from fastapi.testclient import TestClient

from ecd_epitope_scanner.jobs import JobStore
from ecd_epitope_scanner.pipeline import PipelineRequest


def test_legacy_job_summary_is_backfilled_once_without_full_list_reads(tmp_path, monkeypatch):
    store = JobStore(tmp_path)
    store._summary_migration.join(timeout=10)
    job = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    result = {"protein": {"accession": "Q9NZQ7"},
              "summary": {"regions": 2, "recommended_constructs": 4},
              "large_evidence": ["unchanged"]}
    store.update(job["id"], result_json=result)
    with store._connect() as connection:
        connection.execute("ALTER TABLE jobs DROP COLUMN summary_json")
        connection.commit()
    migrated = JobStore(tmp_path)
    migrated._summary_migration.join(timeout=10)
    assert not migrated._summary_migration.is_alive()
    assert migrated.get(job["id"])["result"] == result
    monkeypatch.setattr(migrated, "get", lambda *a, **k: (_ for _ in ()).throw(AssertionError("full read")))
    assert migrated.summaries()[0]["summary"] == {
        "protein": result["protein"], "summary": result["summary"]}
    assert migrated.summaries(owner_id="other") == []
    # A later restart must not overwrite an existing compact summary.
    with migrated._connect() as connection:
        connection.execute("UPDATE jobs SET result_json='invalid-but-unused'")
        connection.commit()
    assert JobStore(tmp_path).summaries()[0]["summary"]["protein"]["accession"] == "Q9NZQ7"


def test_custom_design_keeps_full_api_evidence_with_projected_disk(tmp_path, monkeypatch):
    import ecd_epitope_scanner.api as api
    import ecd_epitope_scanner.exporting as exporting
    from ecd_epitope_scanner.models import ProteinContext, ResidueRecord, MutationChoice

    store = JobStore(tmp_path)
    job = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    result = {
        "protein": asdict(ProteinContext("Q9NZQ7", "test", "human", "AA", 1, 2, 3, "test")),
        "residues": [asdict(ResidueRecord(1, "A", "A", 1))],
        "designs": [], "regions": [], "design_summary": {}, "summary": {"regions": 1},
        "joint_v4_ml": {"candidates": [{"candidate_id": "C1", "evidence": [1, 2]}],
                        "surface_hypotheses": [{"hypothesis_id": "H1"}]},
    }
    result["protein"]["ecd_sequence"] = "AA"
    store.update(job["id"], status="COMPLETE", stage="COMPLETE", result_json=result)
    exporting.write_analysis_artifact(tmp_path / job["id"], result)
    assert "candidates" not in json.loads((tmp_path / job["id"] / "analysis.json").read_text())["joint_v4_ml"]
    monkeypatch.setattr(api, "STORE", store)
    monkeypatch.setattr(api, "evaluate_substitution", lambda *a, **k: MutationChoice(1, "A", "L", 0, 1, None, True))
    monkeypatch.setattr(exporting, "write_scientific_outputs", lambda *a, **k: {})
    monkeypatch.setattr(exporting, "finalize_scientific_bundle", lambda *a, **k: None)
    client = TestClient(api.app)
    response = client.post(f"/api/jobs/{job['id']}/designs", json={"mutations": [{"position": 1, "mutant": "L"}]})
    assert response.status_code == 200, response.text
    full = client.get(f"/api/jobs/{job['id']}").json()["result"]
    assert full["joint_v4_ml"] == result["joint_v4_ml"]
    assert full["designs"][0]["mutations"][0]["mutant"] == "L"
    assert store.get(job["id"])["result"]["joint_v4_ml"] == result["joint_v4_ml"]


def test_real_queue_explores_hidden_residue_before_failure_retries_exhaust_budget():
    from ecd_epitope_scanner.tests.test_site_first_candidates import _inputs
    from ecd_epitope_scanner.joint_design_v4.site_first import build_site_first_queues
    from ecd_epitope_scanner.joint_design_v4.site_first_validation import SiteFirstValidation
    from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
    from ecd_epitope_scanner.tests.test_site_first_validation import pair_row

    queue = build_site_first_queues(*_inputs())
    queue._lanes = {key: value for key, value in queue._lanes.items() if key[1] == "epitope_support"}
    attempted = []
    spurs_calls = []

    def spurs(batch):
        spurs_calls.extend(batch)
        return PredictionResult("success", records=[
            {"mutations": [m.notation for m in row.mutations], "ddg": .2} for row in batch])

    def foldx(batch, *, slot_start, runs):
        attempted.extend({m.residue.uniprot_pos for m in row.mutations} for row in batch)
        return PredictionResult("success", records=[pair_row(row, [3.] * runs, slot_start) for row in batch])

    runner = SiteFirstValidation()
    runner.run_indexed(queue, spurs, foldx, spurs_request_limit=10, foldx_slot_limit=30, max_workers=1)
    assert attempted[0] == {1, 2}
    assert 3 in attempted[1]
    assert runner.foldx_ledger["attempted_slots"] <= 30
    assert len(spurs_calls) <= 10
    assert runner.foldx_ledger["explored_positions"] == [1, 2, 3]


def test_rotated_queue_resumes_pending_candidates_without_repeating_decisions():
    from ecd_epitope_scanner.tests.test_site_first_candidates import _inputs
    from ecd_epitope_scanner.joint_design_v4.site_first import build_site_first_queues

    queue = build_site_first_queues(*_inputs(3))
    lane = ("s", "epitope_support")
    queue._lanes = {lane: queue._lanes[lane]}
    first = queue.head(*lane, limit=99)
    suspended = {c.candidate_id for c in first[1:]}
    queue.record_outcome(*lane, first[0].candidate_id, "SPURS_SCIENTIFIC_FAILURE")
    queue.prioritize_uncovered({1, 2})
    assert 3 in {m.residue.uniprot_pos for m in queue.head(*lane)[0].mutations}
    seen = {first[0].candidate_id}
    for _ in range(100):
        head = queue.head(*lane, limit=1)
        if not head:
            break
        candidate = head[0]
        assert candidate.candidate_id not in seen
        seen.add(candidate.candidate_id)
        queue.record_outcome(*lane, candidate.candidate_id, "SPURS_SCIENTIFIC_FAILURE")
    else:
        raise AssertionError("queue did not exhaust")
    assert suspended.issubset(seen)


def test_spurs_failure_also_advances_to_hidden_positions_with_a_small_request_budget():
    from ecd_epitope_scanner.tests.test_site_first_candidates import _inputs
    from ecd_epitope_scanner.joint_design_v4.site_first import build_site_first_queues
    from ecd_epitope_scanner.joint_design_v4.site_first_validation import SiteFirstValidation
    from ecd_epitope_scanner.joint_design_v4.models import PredictionResult

    queue = build_site_first_queues(*_inputs())
    queue._lanes = {key: value for key, value in queue._lanes.items() if key[1] == "epitope_support"}
    positions = []

    def predict(batch):
        positions.extend({m.residue.uniprot_pos for m in c.mutations} for c in batch)
        return PredictionResult("success", records=[
            {"mutations": [m.notation for m in c.mutations], "ddg": 3.} for c in batch])

    def forbidden_foldx(*a, **k):
        raise AssertionError("SPURS scientific failures must not enter FoldX")

    runner = SiteFirstValidation()
    runner.run_indexed(queue, predict, forbidden_foldx, spurs_request_limit=2, foldx_slot_limit=30)
    assert positions[0] == {1, 2}
    assert 3 in positions[1]
    assert runner.spurs_ledger["unique_requests"] == 2
    assert runner.foldx_ledger["attempted_slots"] == 0
