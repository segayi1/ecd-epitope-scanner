from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue

from ecd_epitope_scanner import structure_analysis
from ecd_epitope_scanner.assay_profiles import get_assay_profile, DEFAULT_ASSAY_PROFILE_ID
from ecd_epitope_scanner.joint_design_v4.budget import resolve_budget
from ecd_epitope_scanner.joint_design_v4.benchmark import validate_benchmark_case
from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates, hypothesis_balanced_shortlist
from ecd_epitope_scanner.joint_design_v4.evidence import build_single_evidence, compressed_options
from ecd_epitope_scanner.joint_design_v4.engine import _formal_candidate_reasons
from ecd_epitope_scanner.joint_design_v4.hierarchy import build_hierarchy
from ecd_epitope_scanner.joint_design_v4.identifiability import panel_metrics, exhaustive_comparison
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4MutationOption, V4Candidate, V4SurfaceHypothesis, V4Sector
from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel, _risk_penalty
from ecd_epitope_scanner.joint_design_v4.qualification import apply_coordinate_qualification
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord
from ecd_epitope_scanner.surface_patch.models import Footprint


def option(p, review=False):
    return V4MutationOption(ResidueKey("0","A",p,uniprot_pos=p,one_letter="D"), "A", f"D{p}A", "charge_removal", 0., 1., True,
                            review_required=review, accessibility_score=.8, chemical_disruption_score=.6)


def candidate(i, response, review=False):
    return V4Candidate(str(i), "DISPLAY", [option(i, review)], order=2,
                       hypothesis_response={f"H{j}": v for j,v in enumerate(response)},
                       qualification="REVIEW_REQUIRED" if review else "FORMAL_POOL")


def analyze_one(monkeypatch, tmp_path, name, confidence=90):
    model, chain = Model.Model(0), Chain.Chain("A")
    model.add(chain)
    residue = Residue.Residue((" ",1," "), name, " ")
    residue.add(Atom.Atom("CA",np.array([0.,0.,0.]),confidence,1.," "," CA ",1,element="C"))
    chain.add(residue)
    parsed = SimpleNamespace(model=model, label_lookup={}, secondary_lookup={})
    mapping = SimpleNamespace(chain_id="A", residues=[residue], uniprot_to_residue={1:residue}, query_to_uniprot={0:1})
    sequence = "G" if name == "GLY" else "D"
    protein = ProteinContext("TEST","test","test",sequence,1,1,1,"test")
    monkeypatch.setattr(structure_analysis, "calculate_sasa", lambda *a: ({("A",1,""):(100.,0.)}, "fixture_native"))
    return structure_analysis.analyze_residues(parsed, mapping, protein, tmp_path, predicted=True)[0][0]


def test_gly_exposure_not_sidechain_size(monkeypatch,tmp_path):
    record = analyze_one(monkeypatch,tmp_path,"GLY")
    assert record.surface_accessible and record.sidechain_sasa == 0
    assert "insufficient-sidechain-SASA" not in record.hard_exclusions
    assert "backbone-conformation-sensitive" in record.risk_flags


def test_low_confidence_is_review_not_fold_damage(monkeypatch,tmp_path):
    record = analyze_one(monkeypatch,tmp_path,"ASP",60)
    assert record.surface_accessible and not record.folding_screen_pass
    assert "low-pLDDT" not in record.hard_exclusions
    assert record.mutation_qualification == "REVIEW_REQUIRED"
    keys = {1: option(1).residue}
    single = build_single_evidence([record],keys,[{"position":1,"mutant":"A","ddg":0.}],{})
    options = compressed_options([record],single)
    assert options and all(o.review_required for o in options)
    c = V4Candidate("review","S",options,order=2)
    assert "STRUCTURAL_CONFIDENCE_REVIEW_REQUIRED" in _formal_candidate_reasons(c, allow_spurs_only=False, minimum_coverage_gain=.05, lower_order_coverage=0.)
    apply_coordinate_qualification([record], {"residues":[{"residue_id":1,"state":"UNCERTAIN","use_for_3d_footprint":False}]})
    assert not any(s.eligible for s in build_single_evidence([record],keys,[{"position":1,"mutant":"A","ddg":0.}],{}))


def hierarchy_fixture(supports=((1,2,3,4),(1,2,3,4))):
    records = [ResidueRecord(p,"D","A",p,ca=(float(p),0.,0.),rsasa=.8,sasa=20.,surface_accessible=True,coordinate_status="STRUCTURED") for p in range(1,7)]
    graph = {p: {q:1. for q in (p-1,p+1) if 1<=q<=6} for p in range(1,7)}
    footprints = []
    for center, support in enumerate(supports,1):
        for radius in (9,12,15):
            footprints.append(Footprint(f"F{center}_{radius}",center,center,"fixed_radius",radius,None,[],
                residues=list(support),residue_surface_area={p:10. for p in support},
                mean_normal=(1.,0.,0.),surface_component_id=1,total_surface_area=10*len(support)))
    return SimpleNamespace(raw_footprints=footprints,footprints=footprints),records,graph


def test_family_scales_and_deterministic_complete_link():
    surface,records,graph = hierarchy_fixture()
    raw,families,hypotheses = build_hierarchy(surface,records,[],graph)
    assert len(raw)==6 and len(families)==2 and len(hypotheses)==1
    assert len(hypotheses[0].source_footprint_ids)==6
    assert families[0].membership == {p:1. for p in range(1,5)}
    surface.raw_footprints.reverse()
    assert build_hierarchy(surface,records,[],graph)[1:] == (families,hypotheses)
    surface,records,graph = hierarchy_fixture(((1,2,3,4),(1,2,3,4,5),(2,3,4,5)))
    assert len(build_hierarchy(surface,records,[],graph,{"weighted_jaccard":.7})[2])==2


@pytest.mark.parametrize("difference",["component","domain","normal","missing_normal","reliability","diameter"])
def test_no_invalid_family_merge(difference):
    surface,records,graph = hierarchy_fixture(((1,2,3,4),(1,2,3,4,5)))
    params = {}
    if difference=="component":
        for f in surface.raw_footprints[3:]: f.surface_component_id=2
    elif difference=="domain": records[4].domain="OTHER"
    elif difference=="normal":
        for f in surface.raw_footprints[3:]: f.mean_normal=(-1.,0.,0.)
    elif difference=="missing_normal":
        for f in surface.raw_footprints[3:]: f.mean_normal=(0.,0.,0.)
    elif difference=="reliability": records[4].coordinate_status="UNCERTAIN"
    else: params["maximum_nominal_geodesic_a"]=3.
    assert len(build_hierarchy(surface,records,[],graph,params)[2])==2


def test_global_pool_independent_of_display_ownership_and_order():
    graph = {p:{q:1. for q in range(1,7) if p!=q} for p in range(1,7)}
    hypotheses = [V4SurfaceHypothesis("H",{p:10. for p in graph})]
    params = {"orders":[2,3],"beam_width_by_order":{2:8,3:8},"global_shortlist_per_order":3}
    one = [V4Sector("ONE",list(graph),1,60,60,0)]
    split = [V4Sector("LEFT",[1,2,3],1,30,30,0),V4Sector("RIGHT",[4,5,6],1,30,30,0)]
    a=generate_candidates(one,[option(p) for p in graph],graph,params,hypotheses)
    b=generate_candidates(split,[option(p) for p in graph],graph,params,hypotheses)
    assert [c.candidate_id for c in a]==[c.candidate_id for c in b]
    assert {c.order for c in a}=={2,3}
    params["orders"]=[3]
    triples=generate_candidates(one,[option(p) for p in graph],graph,params,hypotheses)
    assert [c.candidate_id for c in a if c.order==3]==[c.candidate_id for c in triples]


def test_global_capacity_expands_for_required_reservations():
    hypotheses=[V4SurfaceHypothesis(f"H{i}",{i:1}) for i in range(4)]
    pool=[]
    for i in range(4):
        for j in range(2):
            c=candidate(10*i+j,[0]*4)
            c.targeted_hypotheses=[f"H{i}"]
            pool.append(c)
    audit={}
    chosen=hypothesis_balanced_shortlist(pool,hypotheses,{"global_shortlist_per_order":2},audit)
    assert len(chosen)==8
    assert audit["global_order_pools"][0]["effective_capacity"]==8


def test_adaptive_no_forced_fill_or_review_leak():
    hypotheses=[V4SurfaceHypothesis("H0",{1:1}),V4SurfaceHypothesis("H1",{2:1})]
    pool=[candidate(1,[1,0]),candidate(2,[1,0]),candidate(3,[0,1],True)]
    budget=resolve_budget()
    panel,selected=select_panel(pool,hypotheses=hypotheses,panel_budget=budget)
    assert [c.candidate_id for c in selected]==["1"]
    assert not panel.budget_decision.minimum_target_met
    assert panel.budget_decision.stop_reason=="NO_FEASIBLE_NONREDUNDANT_CANDIDATES"
    assert panel.minimum_signature_distance==pytest.approx(1/np.sqrt(8),abs=1e-6)


def test_adaptive_reaches_eight_with_useful_gain_and_zero_minimum():
    # Two immutable duplicate signatures leave min distance at zero, while
    # independent probes still resolve the other eight hypotheses.
    hypotheses=[V4SurfaceHypothesis(f"H{i}",{i:1}) for i in range(10)]
    pool=[candidate(i+1,[1. if i==j else 0. for j in range(10)]) for i in range(8)]
    panel,selected=select_panel(pool,hypotheses=hypotheses,panel_budget=resolve_budget())
    assert len(selected)==8 and panel.minimum_signature_distance==0
    assert panel.budget_decision.stop_reason=="MAXIMUM_CONSTRUCT_BUDGET_REACHED"
    curve=panel.budget_decision.budget_curve
    assert all(row["net_gain"]>0 for row in curve)
    assert all(a["objective"]<b["objective"] for a,b in zip(curve,curve[1:]))
    assert any(set(group)=={"H8","H9"} for group in panel.unresolved_hypothesis_groups)


def test_exhaustive_small_matrix_and_pending_group_semantics():
    hypotheses=[V4SurfaceHypothesis(f"H{i}",{i:1}) for i in range(4)]
    pool=[candidate(i+1,row) for i,row in enumerate(([1,0,0,0],[0,1,0,0],[0,0,1,0],[.5,.5,.5,.5]))]
    result=exhaustive_comparison(pool,hypotheses,2,_risk_penalty)
    assert result["absolute_gap"]<1e-12 and len(result["optimal_candidate_ids"])==2
    # A-B and B-C can be close without A-C being close.
    metrics=panel_metrics([candidate(9,[0,.2,.4])],hypotheses[:3])
    assert metrics["groups"]==[["H0","H1","H2"]]
    assert len(metrics["near_pairs"])==2


def test_identifiability_diagnostics_separate_near_from_exact_pairs():
    hypotheses=[V4SurfaceHypothesis(f"H{i}",{i:1}) for i in range(8)]
    pool=[candidate(1,[0.]*8),candidate(2,[.05]+[0.]*7),candidate(3,[0.]*8)]
    diagnostics=panel_metrics(pool,hypotheses)["identifiability_diagnostics"]
    assert diagnostics["pair_count"]==28
    assert diagnostics["near_pair_count_rms_lt_0_1"]==28
    assert diagnostics["exact_response_vector_pair_count"]==21
    assert diagnostics["rms_distance_distribution"]["maximum"]==pytest.approx(.05/np.sqrt(8))


@pytest.mark.parametrize("values",[dict(minimum=5,target=4),dict(maximum=9),dict(minimum=0),dict(target=8,maximum=7)])
def test_budget_rejects_invalid_order(values):
    with pytest.raises(ValueError): resolve_budget(**values)


def test_versioned_profile_and_old_parameter_precedence():
    assert DEFAULT_ASSAY_PROFILE_ID.endswith("v3")
    assert get_assay_profile("jurkat_fab_k562_antigen_v2")["design_policy"]["stage1_maximum_constructs"]==6
    assert resolve_budget(legacy_max=6)["adaptive_budget"] is False
    assert resolve_budget(legacy_max=2,maximum=8,target=7,minimum=4,adaptive=True)["maximum_constructs"]==8
    from ecd_epitope_scanner.joint_design_v4.budget import optimizer_budget
    for version in (1,2):
        profile=get_assay_profile(f"jurkat_fab_k562_antigen_v{version}")
        budget=resolve_budget(policy=profile["design_policy"])
        assert budget["maximum_constructs"]==6 and optimizer_budget(profile,budget) is None
        with pytest.raises(ValueError): optimizer_budget(profile,resolve_budget())
    panel,selected=select_panel([candidate(1,[1,0]),candidate(2,[0,1],True)],max_primary=6)
    assert [c.candidate_id for c in selected]==["1"]


def test_benchmark_requires_separate_hidden_truth():
    from ecd_epitope_scanner.joint_design_v4.benchmark import benchmark_contract
    case={key:"test" for key in benchmark_contract()["required_fields"]}
    case["residue_mapping"]=[dict(chain_id="A",auth_seq_id=25,insertion_code="",uniprot_pos=25)]
    with pytest.raises(ValueError): validate_benchmark_case(case)
    case["hidden_truth_path"]="separate_truth.json"
    assert validate_benchmark_case(case)["antigen_chain"]=="test"


def test_api_budget_validation_and_defaults(monkeypatch,tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.jobs import JobStore
    from fastapi.testclient import TestClient
    monkeypatch.setattr(api,"STORE",JobStore(tmp_path/"jobs"))
    monkeypatch.setattr(api.RUNNER,"submit",lambda job_id:None)
    with TestClient(api.app) as client:
        response=client.post("/api/jobs",json={"identifier":"Q15116"})
        assert response.status_code==202
        assert response.json()["request"]["assay_profile"].endswith("v3")
        bad=client.post("/api/jobs",json={"identifier":"Q15116","minimum_constructs":5,"target_constructs":4})
        assert bad.status_code==422
        bad_upload=client.post("/api/jobs/upload",data={"minimum_constructs":"5","target_constructs":"4"},
                               files={"structure":("input.pdb",b"END\n","chemical/x-pdb")})
        assert bad_upload.status_code==422
        fixed=client.post("/api/jobs",json={"identifier":"Q15116","max_primary_constructs":6})
        assert fixed.status_code==202 and fixed.json()["request"]["adaptive_budget"] is None


def test_empty_order_pool_reports_hypothesis_quota_deficits():
    from types import SimpleNamespace
    from ecd_epitope_scanner.joint_design_v4.candidate_search import hypothesis_balanced_shortlist
    audit={}
    assert hypothesis_balanced_shortlist([], [SimpleNamespace(hypothesis_id="H0")], {"orders":[2,3]}, audit)==[]
    assert [(row["order"],row["quota_deficits"]["H0"]) for row in audit["global_order_pools"]]==[(2,2),(3,1)]


def test_missing_predicted_pae_is_not_low_uncertainty():
    surface,records,graph = hierarchy_fixture()
    _,families,hypotheses = build_hierarchy(surface,records,[],graph,predicted=True,confidence={"unmapped":True})
    assert all(f.structural_uncertainty["pae_band"]=="UNAVAILABLE" for f in families)
    assert len(hypotheses)==2
