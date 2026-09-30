"""Replay frozen inputs at three exposure thresholds with real workers.

Every output is separate, and recommendation gates stay enabled. A completed
reference run may supply the 0.25 result; it is checked against frozen input
hashes, V3 policy, and threshold before reuse.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.util import sha256_file, write_json


def summarize(run):
    analysis=json.loads((run/"analysis.json").read_text(encoding="utf-8"))
    residues=analysis["residues"]
    hierarchy=json.loads((run/"surface_complexity_v411.json").read_text(encoding="utf-8"))
    budget=json.loads((run/"panel_budget_decision_v411.json").read_text(encoding="utf-8"))
    candidates=json.loads((run/"panel_candidates.json").read_text(encoding="utf-8"))
    return dict(output_dir=str(run), threshold=analysis["methods"]["rsasa_threshold"],
                exposed_positions=[r["uniprot_pos"] for r in residues if r["surface_accessible"]],
                formal_candidate_count=sum(c["formal_eligible"] for c in candidates),
                review_candidate_count=sum(c.get("qualification")=="REVIEW_REQUIRED" for c in candidates),
                candidate_count=len(candidates),
                exclusion_counts=analysis["design_summary"]["recommendation_diagnostics"]["formal_exclusion_counts"],
                selected_candidate_ids=analysis["design_summary"]["selected_candidate_ids"],
                surface_complexity=hierarchy, budget_decision=budget,
                recommendation_status=analysis["design_summary"]["recommendation_status"],
                worker_failures=analysis["design_summary"]["recommendation_diagnostics"]["worker_failures"],
                model_sha256=sha256_file(run/"structure.cif"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_run",type=Path)
    parser.add_argument("output_root",type=Path)
    parser.add_argument("--reference-run",type=Path)
    args=parser.parse_args()
    source=args.source_run.resolve()
    root=args.output_root.resolve()
    if source==root or source in root.parents or root in source.parents:
        parser.error("Sensitivity outputs must be separate from historical inputs")
    root.mkdir(parents=True,exist_ok=True)
    rows=[]
    for threshold in (.15,.20,.25):
        output=(args.reference_run.resolve() if threshold==.25 and args.reference_run else root/f"rsasa_{threshold:.2f}")
        if not (output/"analysis.json").is_file():
            subprocess.run([sys.executable,str(Path(__file__).with_name("replay_joint_v41_run.py")),str(source),str(output),
                            "--profile-v3","--surface-rsasa-threshold",str(threshold)],check=True)
        analysis=json.loads((output/"analysis.json").read_text(encoding="utf-8"))
        provenance=json.loads((output/"replay_input_provenance.json").read_text(encoding="utf-8"))
        if (analysis["assay_profile"]["version"]!=3 or analysis["methods"]["rsasa_threshold"]!=threshold
            or provenance["source_request_sha256"]!=sha256_file(source/"request.json")
            or provenance["source_structure_sha256"]!=sha256_file(source/"structure.cif")):
            raise ValueError("Sensitivity run does not match frozen inputs, V3 policy or requested threshold")
        rows.append(summarize(output))
    baseline=set(rows[-1]["exposed_positions"])
    for row in rows:
        row["added_vs_025"]=sorted(set(row["exposed_positions"])-baseline)
    write_json(root/"surface_sensitivity_summary_v411.json", {
        "default_threshold_unchanged":.25,"runs":rows,
        "interpretation":"Observed candidate and panel changes at uncalibrated engineering thresholds; worker failures remain explicit.",
    })
    print(json.dumps([{k:row[k] for k in ("threshold","candidate_count","formal_candidate_count","review_candidate_count","selected_candidate_ids")} for row in rows]))


if __name__=="__main__": main()
