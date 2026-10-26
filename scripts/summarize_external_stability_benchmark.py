from __future__ import annotations

"""Independent score/identity audit and descriptive metrics for benchmark ledger."""
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from ecd_epitope_scanner.surface_perturbation.paired import validate_pair  # noqa: E402

BASE = ROOT / "workspace_data" / "p0_p3_upgrade_20261016" / "benchmark"


def read():
    candidates = json.loads((BASE / "candidate_manifest.json").read_text(encoding="utf-8"))["candidates"]
    rows = {str(r["source_variant_id"]): r for r in candidates}
    with (BASE / "candidates.csv").open(encoding="utf-8-sig", newline="") as f:
        audit = {str(r["source_variant_id"]): r for r in csv.DictReader(f)}
    for cid, row in rows.items():
        row["_csv_audit"] = audit.get(cid, {})
    events = [json.loads(x) for x in (BASE / "execution_ledger.jsonl").read_text(encoding="utf-8").splitlines() if x.strip()]
    return rows, events


def finite(values):
    return [float(x) for x in values if x is not None and math.isfinite(float(x))]


def stats(pred, obs):
    from scipy.stats import pearsonr, spearmanr
    pairs = [(float(p), float(o)) for p, o in zip(pred, obs) if p is not None and o is not None and math.isfinite(float(p)) and math.isfinite(float(o))]
    if not pairs:
        return {"n": 0}
    p, o = zip(*pairs)
    err = [a-b for a,b in pairs]
    result = {"n": len(pairs), "mae_kcal_mol": statistics.mean(abs(x) for x in err),
              "rmse_kcal_mol": math.sqrt(statistics.mean(x*x for x in err)),
              "mean_signed_error_kcal_mol": statistics.mean(err),
              "mean_predicted_kcal_mol": statistics.mean(p), "mean_experimental_kcal_mol": statistics.mean(o)}
    result["pearson_r"] = float(pearsonr(p,o).statistic) if len(pairs)>1 and len(set(p))>1 and len(set(o))>1 else None
    result["spearman_rho"] = float(spearmanr(p,o).statistic) if len(pairs)>1 and len(set(p))>1 and len(set(o))>1 else None
    return result


def candidate_pair_identity_matches(candidate, pair):
    expected = candidate["mutations_normalized"]
    if isinstance(expected, str):
        expected = json.loads(expected)
    # The normalized structures are independently verified single-chain A files.
    return (pair.get("chain") == candidate.get("normalized_chain_id", "A") and
            sorted(pair.get("mutations") or []) == sorted(expected) and
            pair.get("normalized_input_sha256") == candidate.get("normalized_pdb_sha256"))


def benchmark_status(all_started, complete, runner_status):
    if all_started:
        return "COMPLETE" if complete else "COMPLETE_WITH_MODEL_FAILURES"
    if runner_status in {"FAILED", "INTERRUPTED_NO_TERMINAL_RECEIPT"}:
        return "INTERRUPTED_INCOMPLETE"
    return "PILOT_OR_RUNNING"


def main():
    rows, events = read()
    runner_contract = json.loads((BASE / "runner_execution_manifest.json").read_text(encoding="utf-8"))
    approved_pair_slots = min(300, int(runner_contract.get("max_pair_slots", 300)),
                              int(runner_contract.get("approved_pair_slots", 300)))
    predictions = {"spurs": {}, "foldx": {}}
    stage_time = defaultdict(float)
    foldx_attempted = set()
    foldx_result_ids = set()
    foldx_failed_ids = set()
    spurs_attempted = set()
    spurs_failures = []
    foldx_pair_count = 0
    pair_identity_failures = []
    spurs_results = 0
    for e in events:
        if e["event"] == "spurs_result":
            stage_time["spurs_worker_elapsed_seconds"] += float(e.get("elapsed_seconds",0))
            if e.get("status") != "success":
                spurs_failures.append({
                    "candidate_ids": [str(cid) for cid in e.get("candidate_ids", [])],
                    "key": e.get("key"),
                    "errors": e.get("payload", {}).get("errors", []),
                    "elapsed_seconds": e.get("elapsed_seconds"),
                })
                continue
            spurs_results += 1
            payload=e["payload"]
            records=payload.get("records",[])
            ids=e.get("candidate_ids",[])
            if payload.get("provenance",{}).get("mode") == "single":
                # SPURS single worker returns a full L x 19 matrix, not one row
                # per selected variant. Map each blinded candidate onto that
                # matrix by normalized sequence position and mutant amino acid.
                by_position={int(r.get("position",0)):r for r in records}
                for cid in ids:
                    row=rows[str(cid)]; mut=row["mutations_normalized"][0]
                    wt,pos,aa=mut[0],int(mut[1:-1]),mut[-1]
                    rec=by_position.get(pos)
                    if not rec or rec.get("wt")!=wt or aa not in rec.get("scores",{}):
                        pair_identity_failures.append({"kind":"SPURS_SINGLE_MAPPING","candidate_id":cid,"position":pos})
                        continue
                    predictions["spurs"][str(cid)]=float(rec["scores"][aa])
                    row["spurs_worker_provenance"] = payload.get("provenance",{})
                continue
            if len(records)!=len(ids):
                pair_identity_failures.append({"kind":"SPURS_RECORD_COUNT","key":e.get("key"),"ids":len(ids),"records":len(records)})
                continue
            for cid, rec in zip(ids,records):
                row=rows[str(cid)]
                notation=row["mutations_normalized"]
                if int(row["order"])==1:
                    mut=notation[0]
                    wt,pos,aa=mut[0],int(mut[1:-1]),mut[-1]
                    if rec.get("position")!=pos or rec.get("wt")!=wt or aa not in rec.get("scores",{}):
                        pair_identity_failures.append({"kind":"SPURS_SINGLE_MAPPING","candidate_id":cid})
                        continue
                    value=rec["scores"][aa]
                else:
                    if sorted(rec.get("mutations",[]))!=sorted(notation):
                        pair_identity_failures.append({"kind":"SPURS_MULTI_MAPPING","candidate_id":cid,"record":rec.get("mutations"),"expected":notation})
                        continue
                    value=rec.get("ddg")
                predictions["spurs"][str(cid)]=float(value)
                row["spurs_worker_provenance"] = payload.get("provenance",{})
        elif e["event"] == "foldx_start":
            foldx_attempted.add(str(e["candidate_id"]))
        elif e["event"] == "spurs_start":
            spurs_attempted.update(str(x) for x in e.get("candidate_ids",[]))
        elif e["event"] == "foldx_result":
            cid=str(e["candidate_id"])
            foldx_result_ids.add(cid)
            stage_time["foldx_worker_elapsed_seconds"] += float(e.get("elapsed_seconds",0))
            if e.get("status")!="success":
                foldx_failed_ids.add(cid)
                continue
            rec_pairs=e.get("paired_models") or []
            run_values=e.get("ddg_runs_kcal_mol") or []
            if len(rec_pairs)!=5 or len(run_values)!=5:
                pair_identity_failures.append({"kind":"FOLDX_REPEAT_COUNT","candidate_id":cid,"pairs":len(rec_pairs),"scores":len(run_values)})
                foldx_failed_ids.add(cid)
                continue
            candidate = rows.get(cid)
            if candidate is None:
                pair_identity_failures.append({"kind":"FOLDX_UNKNOWN_CANDIDATE","candidate_id":cid})
                foldx_failed_ids.add(cid)
                continue
            runs_by_id = {}
            candidate_pairs_valid = True
            for pair in rec_pairs:
                try:
                    validate_pair(pair)
                    run_id = int(pair["run"])
                    if run_id in runs_by_id or not candidate_pair_identity_matches(candidate, pair):
                        raise ValueError("paired model is mapped to a different source candidate")
                    pair_score = pair.get("ddg_kcal_mol")
                    if pair_score is None or not math.isfinite(float(pair_score)):
                        raise ValueError("paired model has no finite mapped DDG")
                    runs_by_id[run_id] = float(pair_score)
                    foldx_pair_count+=1
                except Exception as exc:
                    pair_identity_failures.append({"kind":"FOLDX_PAIR_INVALID","candidate_id":cid,"reason":str(exc)})
                    candidate_pairs_valid = False
            provenance = e.get("provenance", {}) or {}
            try:
                valid_runs = (set(runs_by_id) == set(range(5)) and
                              all(abs(runs_by_id[i] - float(run_values[i])) <= 1e-8 for i in range(5)) and
                              abs(float(e.get("ddg_kcal_mol")) - statistics.median(float(x) for x in run_values)) <= 1e-8 and
                              int(provenance.get("number_of_runs", 0)) == 5 and
                              provenance.get("sampling_slots") == list(range(5)))
            except (TypeError, ValueError):
                valid_runs = False
            if not valid_runs:
                pair_identity_failures.append({"kind":"FOLDX_REPEAT_VALUE_MAPPING","candidate_id":cid})
                candidate_pairs_valid = False
            if candidate_pairs_valid:
                predictions["foldx"][cid]={"median":e.get("ddg_kcal_mol"),"runs":run_values,"pairs":len(rec_pairs),"provenance":provenance}
            else:
                foldx_failed_ids.add(cid)

    pair_slots_reserved = sum(int(e.get("reserved_pair_slots", 5)) for e in events if e.get("event") == "foldx_start")
    if pair_slots_reserved > approved_pair_slots:
        raise RuntimeError(f"FoldX slot ledger exceeds approved cap: {pair_slots_reserved}>{approved_pair_slots}")

    cand_summary=[]
    for cid,row in rows.items():
        sp=predictions["spurs"].get(cid)
        fx=predictions["foldx"].get(cid)
        cand_summary.append({"candidate_id":cid,"source":row["source"],"order":int(row["order"]),"protein":row.get("protein",row.get("cluster","unknown")),
            "chain_length_aa":len(row["sequence"]),"near_100aa_80_120":80<=len(row["sequence"])<=120,
            "experimental_ddg_kcal_mol":row["experimental_ddg_kcal_mol"],
            "experimental_ddg_source_value":row["experimental_ddg_source_value"],
            "ph":row.get("ph"),"temperature_c":row.get("temperature_c"),"method":row.get("method"),"measure":row.get("measure"),
            "reference_train_cluster_overlap":json.loads(row["_csv_audit"].get("train_cluster_overlap") or "{}"),
            "reference_train_single_component_overlap":json.loads(row["_csv_audit"].get("train_single_component_overlap") or "{}"),
            "actual_spurs_megascale_train_homology":"UNKNOWN_NOT_AUDITED",
            "spurs_kcal_mol":sp,"foldx_kcal_mol":fx.get("median") if fx else None,
            "foldx_repeat_values_kcal_mol":fx.get("runs") if fx else None,
            "foldx_repeat_range_kcal_mol":(max(fx["runs"])-min(fx["runs"])) if fx and fx.get("runs") else None,
            "spurs_approx_stable_gate_le_0_5":(sp is not None and sp<=0.5),
            "spurs_acceptance_gate_le_1":(sp is not None and sp<=1.0),
            "foldx_acceptance_gate_le_2":(fx is not None and fx.get("median") is not None and fx["median"]<=2.0),
            "foldx_review_band_2_to_3":(fx is not None and fx.get("median") is not None and 2.0<fx["median"]<=3.0),
            "foldx_pair_proof_count":fx.get("pairs",0) if fx else 0})

    candidate_manifest = json.loads((BASE / "candidate_manifest.json").read_text(encoding="utf-8"))
    out={"schema":"external_stability_benchmark_results_v1","status":"PILOT_OR_RUNNING",
        "runner_execution_status":runner_contract.get("status"),
        "runner_script_sha256":runner_contract.get("runner_script_sha256"),
        "preflight_prediction_status":candidate_manifest.get("summary",{}).get("prediction_status"),
        "score_blind_candidate_count":len(rows),"chain_length_aa":{"min":min(len(r['sequence']) for r in rows.values()),"max":max(len(r['sequence']) for r in rows.values()),
            "median":statistics.median(len(r['sequence']) for r in rows.values()),"n_80_to_120":sum(80<=len(r['sequence'])<=120 for r in rows.values()),
            "n_le_150":sum(len(r['sequence'])<=150 for r in rows.values())},
        "source_order_counts":{str(o):sum(int(r['order'])==o for r in rows.values()) for o in (1,2,3)},
        "actual_spurs_training_overlap":"UNKNOWN_NOT_AUDITED; official model trained on Megascale single and double training splits",
        "PTMUL_reference_train_clusters":"S2648/VB1432 membership is relative to the PTMUL benchmark reference splits, not proof of SPURS Megascale independence",
        "foldx_development_overlap":"FoldX is an empirical force-field method, not a learned ML split model; parameterization/calibration-development data overlap and this benchmark's domain transfer were not fully audited",
        "foldx":{"attempted_candidates":len(foldx_attempted),"result_event_count":len(foldx_result_ids),
            "missing_result_candidates":sorted(foldx_attempted-foldx_result_ids),"failed_result_candidates":sorted(foldx_failed_ids),
            "successful_candidates":len(predictions['foldx']),"pair_slots_reserved":pair_slots_reserved,
            "approved_pair_slots":approved_pair_slots,"hard_pair_slot_cap":300,
            "pair_proofs_validated":foldx_pair_count,
            "pair_validation_failures":sum(x['kind'] in ('FOLDX_PAIR_INVALID','FOLDX_REPEAT_VALUE_MAPPING','FOLDX_REPEAT_COUNT') for x in pair_identity_failures)},
        "spurs":{"attempted_candidate_count":len(spurs_attempted),"successful_worker_batches":spurs_results,"scored_candidate_count":len(predictions['spurs']),"missing_score_count":len(rows)-len(predictions['spurs']),
            "failed_events":spurs_failures,"failed_candidate_ids":sorted({cid for event in spurs_failures for cid in event["candidate_ids"]}),
            "failed_candidate_count":len({cid for event in spurs_failures for cid in event["candidate_ids"]})},
        "prior_attempt_recovery":runner_contract.get("prior_attempt_recovery"),
        "worker_elapsed_seconds":dict(stage_time),"identity_or_pair_failures":pair_identity_failures,
        "metrics_by_order":{},"metrics_by_length_group":{},"candidates":cand_summary,
        "interpretation":"Descriptive external comparison only; sample/source clustering and training overlap are incomplete. Not a universal threshold calibration or binder-specific guarantee."}
    for order in (1,2,3):
        sub=[r for r in cand_summary if r['order']==order]
        out['metrics_by_order'][str(order)]={model:stats([r[key] for r in sub],[r['experimental_ddg_kcal_mol'] for r in sub]) for model,key in [('SPURS','spurs_kcal_mol'),('FoldX','foldx_kcal_mol')]}
        s=[r for r in sub if r['spurs_kcal_mol'] is not None]
        fx=[r for r in sub if r['foldx_kcal_mol'] is not None]
        out['metrics_by_order'][str(order)]['prediction_gate_counts']={
            'SPURS':{'scored':len(s),'le_0_5':sum(r['spurs_kcal_mol']<=0.5 for r in s),'le_1_0':sum(r['spurs_kcal_mol']<=1.0 for r in s),'gt_1_0':sum(r['spurs_kcal_mol']>1.0 for r in s)},
            'FoldX':{'scored':len(fx),'le_2_0':sum(r['foldx_kcal_mol']<=2.0 for r in fx),'2_to_3_review':sum(2.0<r['foldx_kcal_mol']<=3.0 for r in fx),'gt_3_0':sum(r['foldx_kcal_mol']>3.0 for r in fx)}}
    for label, fn in [("80_120aa",lambda x:80<=x<=120),("le_150aa",lambda x:x<=150),("gt_150aa",lambda x:x>150)]:
        sub=[r for r in cand_summary if fn(r['chain_length_aa'])]
        out['metrics_by_length_group'][label]={"n":len(sub),**{model:stats([r[key] for r in sub],[r['experimental_ddg_kcal_mol'] for r in sub]) for model,key in [('SPURS','spurs_kcal_mol'),('FoldX','foldx_kcal_mol')]}}
    all_started = len(foldx_attempted) == len(rows) and len(spurs_attempted) == len(rows)
    complete = (len(foldx_attempted-foldx_result_ids) == 0 and not foldx_failed_ids and
                len(predictions['foldx']) == len(rows) and len(predictions['spurs']) == len(rows) and
                not pair_identity_failures)
    out['status'] = benchmark_status(all_started, complete, runner_contract.get("status"))
    out["prediction_status"] = out["status"]
    out["evidence_sha256"] = {
        "candidate_manifest": hashlib.sha256((BASE / "candidate_manifest.json").read_bytes()).hexdigest(),
        "execution_ledger": hashlib.sha256((BASE / "execution_ledger.jsonl").read_bytes()).hexdigest(),
        "runner_execution_manifest": hashlib.sha256((BASE / "runner_execution_manifest.json").read_bytes()).hexdigest(),
    }
    (BASE/'benchmark_results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
    with (BASE/'benchmark_results.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(cand_summary[0].keys()));writer.writeheader()
        for r in cand_summary:writer.writerow({k:(json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v) for k,v in r.items()})
    print(json.dumps({k:out[k] for k in ('status','chain_length_aa','foldx','spurs','metrics_by_order','metrics_by_length_group','identity_or_pair_failures')},indent=2,ensure_ascii=False))


if __name__=='__main__':
    main()
