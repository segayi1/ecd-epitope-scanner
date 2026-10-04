"""Recalibrate saved real pairs on independent dense SES centers; calibration cost only."""
import argparse,json,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.sparse import save_npz
from ..surface_perturbation.geometry import Mesh
from ..surface_perturbation.reference import canonical_reference
from ..surface_perturbation.fields import compute_pair,reference_id
from ..surface_perturbation.paired import atoms,validate_pair
from ..surface_perturbation.tools import PhysicalTools
from ..surface_perturbation.deadlines import guard_mesh
from ..surface_perturbation.query_cache import memoize_mesh
from ..surface_perturbation.runtime import RunBudget
from ..surface_perturbation.engine import cached_numerical_bounds
from ..surface_perturbation.null_store import save_calibration,protocol_identity
from ..surface_perturbation.classification import classify_fields
from ..util import write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case_output',type=Path)
    parser.add_argument('--candidate',required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--budget-path',type=Path,required=True);parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--paired-source',type=Path,required=True,help='original calibration paired manifest; detector uses the saved production five')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    budget=RunBudget(args.budget_path,14400);report={'status':'INCOMPLETE','scope':'one saved actual construct, independently dense calibrated hard field; not full-panel validation'}
    started=time.perf_counter();tools=PhysicalTools(args.cache,{},budget)
    original=tools.ses;tools.ses=lambda *a,**kw:memoize_mesh(guard_mesh(original(*a,**kw),budget),tools)
    try:
        folder=args.case_output/'surface_perturbation'
        records=[SimpleNamespace(**r) for r in json.loads((args.case_output/'analysis.json').read_text(encoding='utf-8'))['residues']]
        contract=json.loads((folder/'reference_contract.json').read_text(encoding='utf-8'))
        native=atoms(args.case_output/'structure_for_energy.pdb')
        with np.load(folder/'native_WT_SES_mesh.npz') as mesh:surface=Mesh(mesh['vertices'],mesh['faces'])
        reliable={r.uniprot_pos for r in records if r.coordinate_eligible and not r.review_reasons}
        reference=canonical_reference(surface,native,records,reliable,set(contract['designable_positions']),spacing=.5)
        reference['context_identity']=contract['context_identity']
        if not np.isclose(reference['area'].sum(),contract['area_all_a2'],rtol=1e-12):raise ValueError('dense frozen target denominator changed')
        if not np.isclose(reference['area'][reference['designable']].sum(),contract['area_designable_a2'],rtol=1e-12):raise ValueError('dense designable denominator changed')
        np.savez_compressed(args.output/'canonical_reference.npz',**{k:v for k,v in reference.items() if isinstance(v,np.ndarray)})
        save_npz(args.output/'reference_graph.npz',reference['graph'])
        write_json(args.output/'reference_contract.json',{'reference_id':reference_id(reference),'context_identity':reference['context_identity'],
            'reference_spacing_a':.5,'frozen_coarse_contract':contract})
        detector_pairs=json.loads((folder/'paired_manifest.json').read_text(encoding='utf-8'))[args.candidate][:5]
        pairs=json.loads(args.paired_source.read_text(encoding='utf-8'))[args.candidate]
        if len(pairs)<20:raise ValueError('requires 20 existing lawful pairs; no automatic new models in this audit')
        pairs=pairs[:20]
        for row in pairs+detector_pairs:validate_pair(row)
        if len(detector_pairs)!=5:raise ValueError('five real detector pairs required')
        if protocol_identity(pairs,reference,tools)!=protocol_identity(detector_pairs,reference,tools):
            raise ValueError('dense calibration and detector model conditions differ')
        cohort={(row.get('input_sha256'),row.get('normalized_input_sha256'),row.get('binary_sha256'),tuple(sorted(row['mutations']))) for row in pairs+detector_pairs}
        if len(cohort)!=1:raise ValueError('dense repeated model cohort differs')
        write_json(args.output/'paired_manifest.json',{args.candidate:detector_pairs,args.candidate+'_null':pairs})
        runs=[];numerical=None
        for count in (5,10,20):
            for pair in pairs[len(runs):count]:
                budget.remaining();runs.append(compute_pair(pair,native,records,reference,tools,args.output/'fields'))
            if numerical is None:numerical=cached_numerical_bounds(pairs[0],native,records,reference,tools,runs[0],args.output/'fields')
            calibration=save_calibration(args.output/'calibration',pairs[:count],runs,reference,tools,numerical)
            print(f'DENSE {args.candidate}: null repeats={count}',flush=True)
            if all(value.get('requested_scope_converged',False) for value in calibration.values()):break
        detector_runs=[compute_pair(pair,native,records,reference,tools,args.output/'fields') for pair in detector_pairs]
        classified=classify_fields(detector_runs,detector_pairs,native,records,reference,calibration,protocol_identity(pairs,reference,tools))
        np.savez_compressed(args.output/'dense_coverage.npz',direct=classified['direct'],hard=classified['hard'],unknown=~classified['known'],
            two_direct_sites=classified['two_direct_sites'],area=reference['area'],designable=reference['designable'])
        area=reference['area'];fraction=float(area[classified['hard']].sum()/area.sum())
        with np.load(folder/f'{args.candidate}_coverage.npz') as old, np.load(folder/'canonical_reference.npz') as coarse:
            coarse_fraction=float(coarse['area'][old['hard']].sum()/coarse['area'].sum())
            coarse_changed=float(coarse['area'][old['direct']].sum())
        report.update(status='MEASURED_DENSE_SENSITIVITY',candidate_id=args.candidate,null_runs=len(runs),samples=len(area),
            target_area_a2=float(area.sum()),designable_area_a2=float(area[reference['designable']].sum()),
            coarse_hard_center_fraction=coarse_fraction,dense_hard_center_fraction=fraction,
            absolute_hard_center_fraction_difference=abs(fraction-coarse_fraction),
            coarse_actual_changed_area_a2=coarse_changed,dense_actual_changed_area_a2=float(area[classified['direct']].sum()),
            calibrated_channels={name:value['status'] for name,value in calibration.items()},
            no_new_foldx_models=True,interpretation='sensitivity measurement; no automatic convergence or full-panel acceptance')
    except (TimeoutError,ValueError,RuntimeError) as exc:report['reason']=str(exc)
    finally:
        budget.save(False);report.update(audit_wall_seconds=time.perf_counter()-started,calibration_cumulative_seconds=budget.elapsed,
            tool_calls=tools.calls,tool_activity=tools.activity)
        write_json(args.output/'dense_audit.json',report)
    print(json.dumps(report))

if __name__=='__main__':main()
