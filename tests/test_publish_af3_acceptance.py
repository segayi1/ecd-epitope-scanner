import importlib.util
import json
from pathlib import Path
import sqlite3


def test_publication_is_idempotent_and_preserves_source_job(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[1] / "workspace_data/af3_binder_batch_20261017/finish_af3_acceptance.py"
    spec = importlib.util.spec_from_file_location("finish_af3_acceptance", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path / "repo"
    task = root / "workspace_data/batch"
    output = task / "completion_deliveries/source1"
    output.mkdir(parents=True)
    (output / "analysis.json").write_text(json.dumps({"request": {"metadata": {}}}))
    monkeypatch.setattr(module, "ROOT", root)
    monkeypatch.setattr(module, "TASK", task)
    with sqlite3.connect(root / "workspace_data/jobs.sqlite3") as conn:
        conn.execute("""CREATE TABLE jobs (id TEXT PRIMARY KEY, owner_id TEXT, owner_pid INTEGER,
            status TEXT,stage TEXT,progress REAL,message TEXT,request_json TEXT,result_json TEXT,error TEXT,
            workdir TEXT,created_at TEXT,updated_at TEXT)""")
        conn.execute("INSERT INTO jobs(id,owner_id,status,result_json,workdir) VALUES ('source1','local','COMPLETE','original','historical')")
    delivery = {"jobs": [{"source_job_id": "source1", "workdir": str(output)}]}
    first = module.publish(delivery)[0]["job_id"]
    second = module.publish(delivery)[0]["job_id"]
    assert first == second != "source1"
    with sqlite3.connect(root / "workspace_data/jobs.sqlite3") as conn:
        assert conn.execute("SELECT COUNT(*) FROM jobs").fetchone()[0] == 2
        assert conn.execute("SELECT result_json,workdir FROM jobs WHERE id='source1'").fetchone() == ("original", "historical")
        published = json.loads(conn.execute(
            "SELECT result_json FROM jobs WHERE id=?", (first,)).fetchone()[0])
    assert published == json.loads((output / "analysis.json").read_text())
