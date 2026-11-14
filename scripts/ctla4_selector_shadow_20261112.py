import ast, csv, json, copy
from pathlib import Path
from ecd_epitope_scanner.joint_design_v4.models import V4Candidate,V4MutationOption,ResidueKey,CombinationStabilityEvidence,V4Sector,InferenceHypothesis
from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
root=Path(__file__).resolve().parents[1]
p=root/'workspace_data'/'1d812aa2f9be'
def load(name): return json.loads((p/name).read_text(encoding='utf-8'))
raw=load('panel_candidates.json'); candidates=[]
for row in raw:
 if not row['formal_eligible'] or not row['display_eligible']: continue
 row=copy.deepcopy(row)
 row['mutations']=[V4MutationOption(**{**o,'residue':ResidueKey(**o['residue'])}) for o in row['mutations']]
 row['stability']=CombinationStabilityEvidence(**row['stability']) if row['stability'] else None
 candidates.append(V4Candidate(**row))
sectors=[V4Sector(**r) for r in load('probe_sectors_v4.json')]
hyp=[]
for row in load('inference_hypotheses_v411.json'):
 for key in ('residue_weights','membership_fraction','membership'):
  if key in row: row[key]={int(k):v for k,v in row[key].items()}
 hyp.append(InferenceHypothesis(**row))
areas={r['uniprot_pos']:r['surface_area_a2'] for r in load('surface_graph_v41.json')['nodes']}
with (p/'discotope_residues.csv').open(encoding='utf-8', newline='') as file:
 epitope=list(csv.DictReader(file))
priors={ast.literal_eval(x['residue'])['uniprot_pos']:float(x['surface_percentile']) for x in epitope if x['status']=='PREDICTED' and x['surface_percentile']}
saved_budget=load('selected_panel_v4.json')['budget_decision']
budget={key:saved_budget[key] for key in ('minimum_constructs','target_constructs','maximum_constructs','adaptive_budget')}
optimizer=load('analysis.json')['parameters']['joint_v4_ml']['panel_optimizer']
for mode in ('legacy_selector','perturbation_selector'):
 pool=copy.deepcopy(candidates)
 if mode=='legacy_selector':
  for c in pool: c.impact_evidence.pop('selector_policy',None)
 panel,selected=select_panel(pool,sectors=sectors,surface_areas=areas,hypotheses=hyp,panel_budget=budget,higher_order_max=optimizer['higher_order_max_primary'],risk_penalty_weights=optimizer['risk_penalty_weights'],stop_parameters=optimizer,epitope_priors=priors)
 print(json.dumps({'scope':'selector_only_saved_development_candidates','selector':mode,'budget':budget,'scored_surface_priors':len(priors),'qualified_pool':len(pool),'selected':len(selected),'ids':[c.candidate_id for c in selected],'sectors':[c.sector_id for c in selected], 'stop':panel.budget_decision.stop_reason}))
