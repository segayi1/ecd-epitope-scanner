"""Expose verified individual constructs separately from panel convergence."""
from __future__ import annotations

import csv
from pathlib import Path
import numpy as np
from ..stability import AUTO_EXCLUDED_MUTANTS
from ..util import write_json
from .paired import validate_pair


def local_construct(candidate, pairs, scope, area, hard, direct):
    """Summarize the existing strict local gate; never promote a partial panel."""
    mutations=candidate['mutations'];stability=candidate.get('stability') or {}
    if (not candidate.get('formal_eligible') or not 1<=candidate['order']<=3
            or stability.get('consensus_status')!='CONSENSUS_PASS'):
        raise ValueError('local construct lacks formal upstream qualification')
    if any(not m.get('eligible') or m.get('hard_exclusions')
           or m['mutant_aa'] in AUTO_EXCLUDED_MUTANTS
           or m['residue']['one_letter'] in 'CGP' for m in mutations):
        raise ValueError('local construct violates substitution policy')
    if candidate['order']==3:
        required=('REAL_TRIPLE_COMPRESSION_VERIFIED'
            if candidate.get('triple_role')=='TRIPLE_COMPRESSION' else 'REAL_TRIPLE_EXTRA_GAIN_VERIFIED')
        if required not in candidate.get('warnings',[]):
            raise ValueError('real triple role qualification is not verified')
    v3 = scope.get('qualification_contract') == 'positive_lower_bound_information_v3'
    if not ((v3 and scope.get('positive_lower_bound_qualified') and scope.get('holdout_validated')) or scope.get('hard_scan_bound_resolved') or
            (scope.get('complete') and scope.get('measurement_local_complete'))):
        raise ValueError('local physical calibration/measurement is incomplete')
    if len(pairs)<5:
        raise ValueError('local construct needs five actual paired models')
    models=pairs[:5]
    identities={(p['batch_path'],p['group'],p['run']) for p in models}
    if len(identities)!=5:
        raise ValueError('duplicate model run identities')
    for pair in models:
        validate_pair(pair)
        if (sorted(pair['mutations'])!=sorted(models[0]['mutations'])
                or pair['input_sha256']!=models[0]['input_sha256']
                or pair.get('normalized_input_sha256')!=models[0].get('normalized_input_sha256')
                or pair.get('binary_sha256')!=models[0].get('binary_sha256')):
            raise ValueError('inconsistent paired model context')
    values=np.asarray([p.get('ddg_kcal_mol',np.nan) for p in models],float)
    median=float(np.median(values))
    if (not np.isfinite(values).all() or stability.get('foldx_ddg_kcal_mol') is None
            or not np.isclose(median,stability['foldx_ddg_kcal_mol'],rtol=0,atol=1e-9)):
        raise ValueError('stability evidence does not match the five field models')
    area=np.asarray(area,float);hard=np.asarray(hard,bool);direct=np.asarray(direct,bool)
    if (area.ndim!=1 or hard.shape!=area.shape or direct.shape!=area.shape
            or not np.isfinite(area).all() or (area<0).any()
            or not hard.any() or not direct.any() or np.any(direct&~hard)):
        raise ValueError('missing or inconsistent verified perturbation masks')
    mapping=[{'notation':m['notation'], 'chain_id':m['residue']['chain_id'],
              'auth_seq_id':m['residue']['auth_seq_id'],
              'insertion_code':m['residue']['insertion_code'],
              'uniprot_pos':m['residue']['uniprot_pos']} for m in mutations]
    return {'candidate_id':candidate['candidate_id'], 'status':'VERIFIED_LOCAL_CONSTRUCT',
            'mutation_notations':[m['notation'] for m in mutations], 'residue_mapping':mapping,
            'normalized_model_mutations':models[0]['mutations'], 'order':candidate['order'],
            'triple_role':candidate.get('triple_role'),
            'parent_candidate_ids':candidate.get('parent_candidate_ids',[]),
            'model_repeats':5, 'foldx_ddg_kcal_mol':median,
            'foldx_runs_kcal_mol':values.tolist(),
            'spurs_ddg_kcal_mol':stability['spurs_multi_ddg_kcal_mol'],
            'stability_status':stability['consensus_status'],
            'actual_changed_area_a2':float(area[direct].sum()),
            'hard_scan_center_area_a2':float(area[hard].sum()),
            'local_scope_contract':scope.get('scope_contract'),
            'hard_scan_bound_resolved':bool(scope.get('hard_scan_bound_resolved')),
            'qualification_contract':scope.get('qualification_contract','legacy'),
            'footprint_characterization_status':scope.get('footprint_characterization_status','NOT_EVALUATED'),
            'warnings':candidate.get('warnings',[]),
            'mutation_warnings':[{'notation':m['notation'],'warnings':m.get('warnings',[]),
                                  'glycan_risk':m.get('glycan_risk',[])} for m in mutations],
            'paired_model_fingerprints':[
                {k:p.get(k) for k in ('batch_path','group','run','wt_sha256','mut_sha256',
                                      'input_sha256','normalized_input_sha256','binary_sha256')}
                for p in models],
            'interpretation':'Individual computational qualification only; not global panel convergence, resolved electrostatics or experimental binding validation.'}


def export_local_constructs(output, records, panel_status):
    output=Path(output)
    write_json(output/'surface_perturbation_local_constructs.json',
               {'panel_status':panel_status,'constructs':records,
                'scope':'verified individual constructs; the panel status remains independent'})
    with (output/'surface_perturbation_local_constructs.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        names=['candidate_id','mutations_uniprot','author_residues','status','panel_status',
               'model_repeats','spurs_ddg_kcal_mol','foldx_ddg_kcal_mol',
               'actual_changed_area_a2','hard_scan_center_area_a2']
        writer=csv.DictWriter(stream,fieldnames=names);writer.writeheader()
        for item in records:
            writer.writerow({**{k:item[k] for k in names if k in item},'panel_status':panel_status,
                'mutations_uniprot':';'.join(item['mutation_notations']),
                'author_residues':';'.join(f"{m['chain_id']}:{m['auth_seq_id']}{m['insertion_code']}"
                                          for m in item['residue_mapping'])})
