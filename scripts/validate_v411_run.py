"""Check real run geometry, evidence isolation, physical identity and bundle."""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.joint_design_v4.sectors import _distance
from ecd_epitope_scanner.util import write_json


def validate(run):
    load=lambda name: json.loads((run/name).read_text(encoding="utf-8"))
    analysis=load("analysis.json")
    snapshot=load("surface_graph_v4.json")
    # Snapshot edges are serialized as position-based undirected pairs.
    graph={int(n["uniprot_pos"]):{} for n in snapshot["nodes"]}
    for edge in snapshot["edges"]:
        a,b,d=int(edge["source"]),int(edge["target"]),float(edge["distance_a"])
        graph[a][b]=d
        graph[b][a]=d
    sectors=load("probe_sectors_v4.json")
    owners=[p for s in sectors for p in s["residue_positions"]]
    assert len(owners)==len(set(owners)) and set(owners)==set(graph)
    for sector in sectors:
        positions=set(sector["residue_positions"])
        assert all(set(_distance(graph,p,positions))==positions for p in positions)
        assert sector["geodesic_diameter_a"]<=26+1e-3
        assert sector["euclidean_diameter_a"]<=18+1e-3
    candidates=load("panel_candidates.json")
    signatures=[tuple(sorted((m["residue"]["uniprot_pos"],m["mutant_aa"]) for m in c["mutations"])) for c in candidates]
    assert len(signatures)==len(set(signatures))
    panel=load("selected_panel_v4.json")
    selected=set(panel["selected_candidate_ids"])
    for candidate in candidates:
        if candidate["candidate_id"] in selected:
            assert candidate["qualification"]!="REVIEW_REQUIRED"
            assert not any(m.get("review_required") for m in candidate["mutations"])
            if analysis["request"]["ml_mode"]=="strict": assert candidate["formal_eligible"]
    budget=load("panel_budget_decision_v411.json")
    assert len(selected)<=budget["maximum_constructs"]<=8
    rows=budget["budget_curve"]
    assert all(r["net_gain"]>0 for r in rows)
    assert all(a["objective"]<b["objective"] for a,b in zip(rows,rows[1:]))
    residues={r["uniprot_pos"]:r for r in analysis["residues"]}
    assert not any(not residues[p]["coordinate_eligible"] for p in graph)
    manifest=load("run_manifest.json")
    with zipfile.ZipFile(run/"ecd_epitope_design.zip") as archive:
        assert set(archive.namelist())==set(manifest["files"])|{"run_manifest.json"}
        for name,digest in manifest["files"].items():
            assert hashlib.sha256((run/name).read_bytes()).hexdigest()==digest
            assert hashlib.sha256(archive.read(name)).hexdigest()==digest
    return dict(output_dir=str(run), recommendation_status=analysis["design_summary"]["recommendation_status"],
                complexity=load("surface_complexity_v411.json"), selected_count=len(selected),
                budget_decision=budget, unique_candidates=len(signatures), bundle_files=len(manifest["files"]),
                invariants="PASS", omitted_geometry_count=sum(not r["coordinate_eligible"] for r in residues.values()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir",type=Path)
    args=parser.parse_args()
    result=validate(args.run_dir.resolve())
    write_json(args.run_dir/"upgrade_validation.json",result)
    print(json.dumps(result,ensure_ascii=False))


if __name__=="__main__":main()
