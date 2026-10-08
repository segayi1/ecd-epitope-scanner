"""Final-panel electrostatics audit; never used by the hard optimizer."""
from __future__ import annotations

from pathlib import Path
import shutil
import numpy as np
from scipy.spatial import cKDTree
from ..util import write_json
from .tools import make_grid
from .fields import aligned_atoms
from .chemistry import RADII
from .calibration import calibrate_channel, resolved_changes
from .pqr_surface import pqr_surface


def write_atoms(path,atom_map,residue_ids=None):
    lines=[]
    sequence_order=residue_ids or {key:i for i,key in enumerate(dict.fromkeys(k[:3] for k in atom_map))}
    for i,key in enumerate(sorted(atom_map,key=lambda key:(sequence_order[key[:3]],key[3])),1):
        chain,res,icode,name=key; xyz,resname,element=atom_map[key]
        if residue_ids is not None: res=residue_ids[key[:3]];icode=''
        x,y,z=xyz
        lines.append(f"ATOM  {i:5d} {name:>4s} {resname:3s} {chain:1s}{res:4d}{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}")
    Path(path).write_text("\n".join(lines)+"\nTER\nEND\n",encoding="ascii")


def protonation_audit(wt,mut,mutated_author):
    w={};m={}
    for key,value in wt.items(): w.setdefault(key[:2],[]).append((key[2],value["residue"],value["charge"],value["radius"]))
    for key,value in mut.items(): m.setdefault(key[:2],[]).append((key[2],value["residue"],value["charge"],value["radius"]))
    changes=[{"residue":list(key),"wt":sorted(w.get(key,[])),"mut":sorted(m.get(key,[]))}
             for key in sorted(w.keys()|m.keys()) if sorted(w.get(key,[]))!=sorted(m.get(key,[])) and key not in mutated_author]
    return {"nonmutated_parameter_changes":changes,"status":"CONDITIONAL_PREPROCESSING_CHANGE" if changes else "MATCHED"}


def solvent_mask(atom_map,points):
    keys=sorted(atom_map); xyz=np.array([atom_map[k][0] for k in keys]); tree=cKDTree(xyz)
    radii=np.array([RADII.get(atom_map[k][2],1.7) for k in keys])
    outside=np.ones(len(points),bool)
    for i,p in enumerate(points):
        ids=tree.query_ball_point(p,radii.max())
        if ids: outside[i]=bool(np.all(np.linalg.norm(xyz[ids]-p,axis=1)>=radii[ids]))
    return outside


def grid_resolution(base,variations,*,absolute_tolerance=.05,relative_tolerance=.1):
    """Versioned engineering convergence screen, not statistical confidence.

    Variations are fine 0.4 A, coarse 0.6 A, and enlarged-box 0.5 A.
    Refinement must decrease the discrepancy and every tested discrepancy must
    fit the declared absolute/relative budget. Failure leaves a point unknown.
    """
    base=np.asarray(base,float);variations=np.asarray(variations,float)
    if variations.shape!=(3,*base.shape):raise ValueError('three matched grid/box variations required')
    if absolute_tolerance<0 or relative_tolerance<0:raise ValueError('negative grid tolerance')
    errors=np.abs(variations-base)
    bound=np.max(errors,axis=0)
    finite=np.isfinite(base)&np.isfinite(variations).all(axis=0)
    valid=finite&(errors[0]<=errors[1]+np.finfo(float).eps*64)&(
        bound<=absolute_tolerance+relative_tolerance*np.abs(base))
    return np.where(valid,bound,np.nan),valid,errors


def prepare_candidate(candidate_id,pairs,native_atoms,records,tools,folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    prepared=[]; audit=[]; all_xyz=[]
    for index,pair in enumerate(pairs[:5]):
        wt,mut,alignment=aligned_atoms(pair,native_atoms,records)
        # Author IDs need not sort in biological sequence order (CTLA4 has
        # author 1A before author 1). Native PDB record chronology is retained.
        pair_residues={k[:3] for k in wt}|{k[:3] for k in mut}
        native_order=list(dict.fromkeys(k[:3] for k in native_atoms if k[:3] in pair_residues))
        if set(native_order)!=pair_residues:
            raise ValueError('PQR residue sequence order unresolved in native reference')
        residue_ids={key:i for i,key in enumerate(native_order,1)}
        wt_path=folder/f"{candidate_id}_{index}_WT.pdb";mut_path=folder/f"{candidate_id}_{index}_Mut.pdb"
        write_atoms(wt_path,wt,residue_ids);write_atoms(mut_path,mut,residue_ids)
        wp,wq=tools.prepare_pqr(wt_path);mp,mq=tools.prepare_pqr(mut_path)
        # Preserve the actual atom charge/radius inputs outside the tool cache.
        portable_w=folder/f'{candidate_id}_{index}_WT.pqr';portable_m=folder/f'{candidate_id}_{index}_Mut.pqr'
        shutil.copyfile(wp,portable_w);shutil.copyfile(mp,portable_m)
        # PDB2PQR may move/add heavy atoms. Keep this as a distinct preparation audit.
        # A geometry-preserving terminal rotation remains conditional: chemical
        # orientation can alter electrostatics even with an unchanged backbone.
        drifts=[];rotations=[]
        backbone={'N','CA','C','O','OXT'}
        for kind,source,pqr in (("WT",wt,wq),("Mut",mut,mq)):
            source_keys={(k[0],residue_ids[k[:3]],k[3]) for k in source}
            displaced={}
            for key,value in source.items():
                pkey=(key[0],residue_ids[key[:3]],key[3])
                if pkey not in pqr: drifts.append({"kind":kind,"atom":list(key),"status":"MISSING"})
                else:
                    distance=float(np.linalg.norm(np.array(pqr[pkey]["xyz"])-value[0]))
                    if distance>.002:
                        entry={"kind":kind,"atom":list(key),"displacement_a":distance}
                        drifts.append(entry);displaced.setdefault(key[:3],[]).append(entry)
            for key in pqr:
                if key not in source_keys and not key[2].lstrip('0123456789').startswith('H'):
                    drifts.append({'kind':kind,'pqr_atom':list(key),'status':'ADDED_HEAVY_ATOM'})
            # Diagnose terminal rotations without waiving preparation QC.
            for residue,moves in displaced.items():
                names={m['atom'][3] for m in moves}
                residue_atoms={k[3] for k in source if k[:3]==residue}
                anchor=(backbone&residue_atoms)|{a for a in residue_atoms if a in {'CB','CG','CG1','CG2','CD','CD1','CD2','CE','CE1','CE2','CZ','NE','ND1','NE1','NE2','NZ','NH1','NH2','OD1','OG','OG1','OH','SD','SG'}}
                anchor_moved=bool(names&anchor-{'OD1','OG','OG1','OH'})
                terminal={n for n in names if n in {'OD1','ND2','OE1','OE2','OG','OG1','OH','NE','ND1','NE2','NZ','NH1','NH2'}}
                if not terminal or names-terminal or anchor_moved:continue
                fixed=[a for a in anchor if a not in names]
                if not fixed or any((residue[0],residue_ids[residue],a) not in pqr for a in fixed):continue
                preserved=all(
                    np.isclose(np.linalg.norm(np.array(source[(*residue,a)][0])-
                                                np.array(source[(*residue,b)][0])),
                               np.linalg.norm(np.array(pqr[(residue[0],residue_ids[residue],a)]["xyz"])-
                                              np.array(pqr[(residue[0],residue_ids[residue],b)]["xyz"])),
                               atol=1e-3,rtol=0)
                    for a in terminal for b in fixed)
                if not preserved:continue
                rotations.append({'kind':kind,'residue':[residue[0],residue[1],residue[2]],
                                  'atoms':[m['atom'] for m in moves],
                                  'max_displacement_a':max(m['displacement_a'] for m in moves),
                                  'status':'CONDITIONAL_TERMINAL_ROTATION',
                                  'backbone_unchanged':True,'anchor_distances_preserved':True,
                                  'electrostatic_equivalence':'NOT_ESTABLISHED'})
        parameter=protonation_audit(wq,mq,{(x[0],residue_ids[tuple(x)]) for x in alignment["mutated_author_residues"]})
        parameter['residue_mapping']=[{'author':list(k),'pqr_number':v} for k,v in residue_ids.items()]
        parameter['residue_order_policy']='native input record chronology; never lexicographic author ID order'
        audit.append({"run":index,**parameter,"heavy_atom_preparation_changes":drifts,
                      "sidechain_rotation_diagnostics":rotations})
        all_xyz.extend(v["xyz"] for v in wq.values());all_xyz.extend(v["xyz"] for v in mq.values())
        prepared.append((portable_w,portable_m,wq,mq))
    if len(prepared)<5: raise ValueError("APBS requires five verified paired models")
    return prepared,audit,all_xyz


def audit_candidate(candidate_id,pairs,native_atoms,records,reference,tools,folder,convergence=False,*,prepared_data=None,shared_grid=None):
    from .execution import potential_pair,solvent_distances
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    prepared,audit,all_xyz=prepared_data or prepare_candidate(candidate_id,pairs,native_atoms,records,tools,folder)
    # One box/grid covering all runs of this candidate. Never per-molecule auto-centering.
    grid=shared_grid or make_grid(all_xyz,.5)
    deltas=[];wt_maps=[];masks=[];offset_deltas=[];offset_masks=[];raw_offsets=[]
    preparation_rejected=any(row['heavy_atom_preparation_changes'] or row['nonmutated_parameter_changes'] for row in audit)
    queries=[reference['points']+d*reference['normals'] for d in (1.,1.5,2.)]
    for wp,mp,wt,mut in prepared:
        wfield,mfield=potential_pair(tools,wp,mp,grid)
        # Query common points; mediator crossing is unknown, never moved to fit Mut.
        # SES triangle occupancy, not atom-ball clearance masquerading as SES.
        if preparation_rejected:
            # This confirmed preparation gate already rejects the entire
            # electrostatic channel. Preserve all five paired solves, but do
            # not spend SES queries assigning a medium we cannot validate.
            same_solvent=np.zeros((3,len(reference['points'])),bool)
        elif hasattr(tools,'cache'):
            wd,wvalid=solvent_distances(tools,wt,np.vstack(queries))
            md,mvalid=solvent_distances(tools,mut,np.vstack(queries))
        else:
            wd,wvalid=pqr_surface(tools,wt).signed_distance(np.vstack(queries))
            md,mvalid=pqr_surface(tools,mut).signed_distance(np.vstack(queries))
        if not preparation_rejected:same_solvent=(wvalid&mvalid&(wd>0)&(md>0)).reshape(3,-1)
        known=same_solvent[1]
        values=mfield(queries[1])-wfield(queries[1]);values[~known]=np.nan
        deltas.append(values);wt_maps.append(wfield(queries[1]));masks.append(known)
        sensitivities=np.stack([mfield(q)-wfield(q) for q in queries])
        raw_offsets.append(sensitivities.copy())
        sensitivities[~same_solvent]=np.nan
        offset_deltas.append(sensitivities);offset_masks.append(same_solvent)
    numerical=np.full(len(reference["points"]),np.nan)
    convergence_audit={"status":"NOT_REQUESTED"}
    if convergence:
        wp,mp,wt,mut=prepared[0];points=reference["points"]+1.5*reference["normals"]
        variations=[];grids=[]
        for spacing,margin in ((.4,16.),(.6,16.),(.5,24.)):
            alternate=make_grid(all_xyz,spacing,margin);grids.append(alternate)
            wfield,mfield=potential_pair(tools,wp,mp,alternate)
            variations.append(mfield(points)-wfield(points))
        settings=getattr(tools,'config',{})
        absolute=float(settings.get('apbs_grid_absolute_tolerance_kt_e',.05))
        relative=float(settings.get('apbs_grid_relative_tolerance',.1))
        base=raw_offsets[0][1] if preparation_rejected else deltas[0]
        raw_numerical,raw_converged,grid_errors=grid_resolution(base,variations,
            absolute_tolerance=absolute,relative_tolerance=relative)
        numerical=np.full_like(raw_numerical,np.nan) if preparation_rejected else raw_numerical
        converged=np.zeros_like(raw_converged) if preparation_rejected else raw_converged
        np.savez_compressed(folder/f'{candidate_id}_grid_sensitivity.npz',
            base_delta=base,variant_deltas=variations,absolute_errors=grid_errors,
            converged=converged,numerical_error=numerical,raw_numerical_error=raw_numerical,raw_grid_converged=raw_converged)
        convergence_audit={"status":"NUMERICAL_ONLY_PREPARATION_REJECTED" if preparation_rejected else "COMPLETE" if converged.all() else "PARTIAL_OR_UNCONVERGED","grids":grids,
                           "converged_queries":int(converged.sum()),"unknown_queries":int((~converged).sum()),
                           "policy_version":"grid_box_resolution_preview_1",
                           "absolute_tolerance_kt_e":absolute,"relative_tolerance":relative,
                           "interpretation":"provisional numerical error budget; not biological confidence; freeze after benchmarks",
                           "max_tested_delta_difference_kt_e":float(np.nanmax(grid_errors))}
    null=np.diff(np.array(wt_maps),axis=0)
    protocol=candidate_id+":apbs:paired-grid-v1"
    calibration=calibrate_channel(null,numerical,protocol,"apbs")
    changed,unknown=resolved_changes(np.array(deltas),calibration,protocol)
    # Preprocessing changes can be a real protonation response, but do not silently
    # equate it with the fixed-protocol direct substitution channel.
    if preparation_rejected:
        unknown[:]=True;changed[:]=False
    np.savez_compressed(folder/f"{candidate_id}_electrostatics.npz",delta_phi=np.array(deltas),
                        wt_phi=np.array(wt_maps),solvent_valid=np.array(masks),offset_deltas=np.array(offset_deltas),
                        offset_solvent_valid=np.array(offset_masks),
                        raw_preparation_conditional_offset_deltas=np.array(raw_offsets),
                        medium_evaluated=np.array(not preparation_rejected),
                        numerical_error=numerical,resolved=changed,unknown=unknown,
                        null_blocks=null,q95=calibration['q95'],q99=calibration['q99'],
                        threshold=calibration['threshold'],calibration_valid=calibration['valid'])
    report={"candidate_id":candidate_id,"grid":grid,"unit":"kT/e","sampling_offset_a":1.5,
            "offset_sensitivity_a":[1.,2.],"runs":audit,"convergence":convergence_audit,
            "null_calibration_status":calibration['status'],
            "threshold_rule":"max(matched empirical Q95, measured numerical error); 4/5 same direction",
            "medium_classification":"NOT_EVALUATED_PREPARATION_CONTRACT_FAILED; false validity means unknown" if preparation_rejected else "NanoShaper SES with actual PQR radii; APBS smoothed interface remains discretization-dependent",
            "preparation_gate_rejected":preparation_rejected,
            "skipped_unusable_SES_query_batches":10 if preparation_rejected else 0,
            "unknown_common_medium_queries":int(np.sum(~np.array(masks))),
            "status":"COMPLETE" if not unknown.any() else "UNRESOLVED",
            "evidence_states":{
                "APBS_EXECUTION_COMPLETE":True,
                "APBS_PREPARATION_VALID":not preparation_rejected,
                "APBS_NUMERICALLY_RESOLVED":bool(not preparation_rejected and np.isfinite(numerical).all()),
                "EXTENDED_COVERAGE_RESOLVED":bool(not unknown.any())},
            "interpretation":"conditional electrostatic response, excluded from hard optimization"}
    write_json(folder/f"{candidate_id}_electrostatics.json",report)
    return changed,unknown,report
