"""Real combination fields on native WT samples, with matched WT subtraction."""
from __future__ import annotations

import hashlib
import inspect
import re
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .geometry import align_pair
from .paired import atoms, validate_pair
from .chemistry import chemical_field, CHANNELS, RADII
from .tools import fingerprint
from ..util import sha256_file, write_json

# Raw numerical fields do not depend on menu selection or null qualification.
# Conservative source fingerprints still invalidate on any numerical/tool code
# change; a classification-only change reuses the fields but re-audits nulls.
IMPLEMENTATION_ID=hashlib.sha256(b''.join((Path(__file__).parent/name).read_bytes()
    for name in ('fields.py','geometry.py','chemistry.py','tools.py'))+inspect.getsource(atoms).encode()).hexdigest()


def reference_id(reference):
    digest=hashlib.sha256()
    for key in ("points","normals","area","owner","domain"):
        if key not in reference: continue
        digest.update(np.asarray(reference[key]).tobytes())
    for value in (reference['graph'].data,reference['graph'].indices,reference['graph'].indptr):
        digest.update(value.tobytes())
    digest.update(reference.get('context_identity','').encode())
    return digest.hexdigest()


def aligned_atoms(pair, native_atoms, records):
    validate_pair(pair)
    wt,mut=atoms(Path(pair["wt_path"])),atoms(Path(pair["mut_path"]))
    lookup={i:(r.chain_id,int(r.auth_seq_id),(r.insertion_code or "").strip()) for i,r in enumerate(records,1)}
    def project(source):
        output={}
        for key,value in source.items():
            if key[1] not in lookup or key[0]!=pair["chain"]: continue
            output[(*lookup[key[1]],key[3])]=value
        return output
    wt,mut=project(wt),project(mut)
    mutated={int(re.fullmatch(r"[A-Z](\d+)[A-Z]",m).group(1)) for m in pair["mutations"]}
    excluded={lookup[p] for p in mutated}
    reliable={(r.chain_id,int(r.auth_seq_id),(r.insertion_code or "").strip()) for r in records if r.coordinate_eligible and not r.review_reasons}
    core=[k for k in native_atoms if k[3] in {"N","CA","C"} and k[:3] in reliable-excluded]
    rotation,translation,rmsd=align_pair({k:v[0] for k,v in native_atoms.items()},
                                       {k:v[0] for k,v in wt.items()},
                                       {k:v[0] for k,v in mut.items()},core)
    # Whole-pair transform. Retain unchanged native-context occluders omitted by ECD normalization.
    context={k:v for k,v in native_atoms.items() if k[:3] not in set(lookup.values())}
    wt=context|{k:(v[0]@rotation+translation,v[1],v[2]) for k,v in wt.items()}
    mut=context|{k:(v[0]@rotation+translation,v[1],v[2]) for k,v in mut.items()}
    return wt,mut,{"core_rmsd_a":rmsd,"core_atoms":len(core),"rotation":rotation.tolist(),
                   "translation":translation.tolist(),"mutated_author_residues":[list(k) for k in sorted(excluded)]}


def spatial_field(atom_map, reference, tools, scale=2., shell_step=1., query_mask=None):
    mesh=tools.ses(atom_map,scale=scale)
    all_points=reference['points']; total=len(all_points)
    active=np.arange(total) if query_mask is None else np.flatnonzero(query_mask)
    points,normals=all_points[active],reference["normals"][active]
    offsets=np.arange(-2.,4.01,shell_step)
    shell=(points[:,None,:]+offsets[None,:,None]*normals[:,None,:]).reshape(-1,3)
    # Fixed Voronoi ownership prevents two normal columns counting the same spatial query.
    _,assigned=cKDTree(all_points).query(shell)
    ownership=assigned==np.repeat(active,len(offsets))
    sdf,valid=mesh.signed_distance(shell)
    sdf=sdf.reshape(len(points),-1); valid=(valid&ownership).reshape(sdf.shape)
    middle=int(round(2/shell_step))
    central=sdf[:,middle].copy(); central[~valid[:,middle]]=np.nan
    chemistry=chemical_field(atom_map,points+normals*1.5,normals)
    atom_keys=sorted(atom_map);atom_xyz=np.array([atom_map[k][0] for k in atom_keys])
    radii=np.array([RADII.get(atom_map[k][2],1.7) for k in atom_keys])
    neighbors=cKDTree(atom_xyz).query_ball_point(shell,float(radii.max()))
    atomic_inside=np.array([bool(ids and np.any(np.linalg.norm(atom_xyz[ids]-point,axis=1)<radii[ids]))
                            for point,ids in zip(shell,neighbors)],float).reshape(sdf.shape)
    _, normal_valid, projected_normals=mesh.signed_distance(points,return_projection=True)
    projected_normals[~normal_valid]=np.nan
    face_distance, face_owner=cKDTree(all_points).query(mesh.centers)
    compatible=(face_distance<=4.) & (np.sum(mesh.normals*reference['normals'][face_owner],axis=1)>.5)
    local_area=np.bincount(face_owner[compatible],weights=mesh.area[compatible],minlength=total)[active]
    if getattr(tools,'config',{}).get('optimization_profile') == 'coverage_scan_v3_preview' and not mesh.closed:
        local_area[:]=np.nan
    def expand(value,fill=np.nan):
        output=np.full((total,*value.shape[1:]),fill,dtype=value.dtype);output[active]=value;return output
    return {"sdf":expand(central),"shell_sdf":expand(sdf),"shell_valid":expand(valid,False),"shell_owned":expand(ownership.reshape(sdf.shape),True),
            "shell_step":shell_step,"chemistry":expand(chemistry),
            "atomic_inside":expand(atomic_inside),
            "projected_normals":expand(projected_normals),"local_ses_area":expand(local_area),
            "ses_total_area_a2":float(mesh.area.sum()),"mesh_vertices":len(mesh.vertices),
            "mesh_faces":len(mesh.faces),"mesh_closed":mesh.closed}


def difference(wt,mut):
    valid=wt["shell_valid"]&mut["shell_valid"]
    inside_w=wt["shell_sdf"]<0; inside_m=mut["shell_sdf"]<0
    # Integral over sampled normal columns divided by area: equivalent thickness (A), not A^2.
    step=float(wt.get("shell_step",1.))
    gain=((inside_m&~inside_w)&valid).sum(axis=1).astype(float)*step
    loss=((inside_w&~inside_m)&valid).sum(axis=1).astype(float)*step
    known=(valid|~wt.get("shell_owned",np.ones_like(valid))).all(axis=1)
    gain[~known]=np.nan; loss[~known]=np.nan
    channels={"shape_sdf_a":mut["sdf"]-wt["sdf"],"occupancy_gain_a":gain,"occupancy_loss_a":loss}
    if 'atomic_inside' in wt:
        left,right=wt['atomic_inside'],mut['atomic_inside']
        known_atomic=np.isfinite(left).all(axis=1)&np.isfinite(right).all(axis=1)
        owned=wt.get('shell_owned',np.ones_like(valid))&mut.get('shell_owned',np.ones_like(valid))
        atomic_gain=((right>0)&(left==0)&owned).sum(axis=1).astype(float)*step
        atomic_loss=((left>0)&(right==0)&owned).sum(axis=1).astype(float)*step
        atomic_gain[~known_atomic]=np.nan;atomic_loss[~known_atomic]=np.nan
        channels.update(atom_occupancy_gain_a=atomic_gain,atom_occupancy_loss_a=atomic_loss)
    for i,name in enumerate(CHANNELS): channels[name]=mut["chemistry"][:,i]-wt["chemistry"][:,i]
    if 'projected_normals' in wt:
        for i,name in enumerate(('surface_normal_x','surface_normal_y','surface_normal_z')):
            channels[name]=mut['projected_normals'][:,i]-wt['projected_normals'][:,i]
        channels['local_ses_area_a2']=mut['local_ses_area']-wt['local_ses_area']
    return channels


def compute_pair(pair,native_atoms,records,reference,tools,cache):
    local_qc = getattr(tools, 'config', {}).get('optimization_profile') == 'coverage_scan_v3_preview'
    contract={"reference":reference_id(reference),"wt":pair["wt_sha256"],"mut":pair["mut_sha256"],
              "implementation_sha256":IMPLEMENTATION_ID,
              "nano":tools.identity(tools.nano),"field_version":5 if local_qc else 4,"scale":2.,"shell":[-2,4,1],"chemistry_resolution_a":1.,"local_query_a":7.}
    path=Path(cache)/("field_"+fingerprint(contract)+".npz")
    manifest_path=path.with_suffix('.json')
    if path.exists() and manifest_path.exists():
        manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest['sha256']!=sha256_file(path) or manifest['contract']!=contract:
            raise ValueError('field cache fingerprint mismatch')
        with np.load(path,allow_pickle=False) as data: return {k:data[k] for k in data.files}
    w,m,audit=aligned_atoms(pair,native_atoms,records)
    if audit["core_rmsd_a"]>.5:
        raise ValueError("REFERENCE_UNRELIABLE: core RMSD exceeds 0.5 A engineering guard")
    changed_xyz=[]
    for key in w.keys()|m.keys():
        left,right=w.get(key),m.get(key)
        if left is None or right is None or left[1]!=right[1] or np.linalg.norm(left[0]-right[0])>.001:
            if left is not None: changed_xyz.append(left[0])
            if right is not None: changed_xyz.append(right[0])
    active=np.zeros(len(reference['points']),bool)
    if changed_xyz:
        separation,_=cKDTree(np.asarray(changed_xyz)).query(reference['points'])
        active=separation<=7.
    # Outside the bounded local query is UNKNOWN, never declared unchanged.
    # Full SES meshes and their topology checks still use the entire structure.
    wf,mf=spatial_field(w,reference,tools,query_mask=active),spatial_field(m,reference,tools,query_mask=active)
    channels=difference(wf,mf)
    topology_changed=tools.ses(w).euler_characteristic!=tools.ses(m).euler_characteristic
    if topology_changed and not local_qc:
        # Euler mismatch proves topology change; matching Euler alone is not
        # proof of complete topological correspondence.
        for value in channels.values(): value[:]=np.nan
    result={**channels,"wt_sdf":wf["sdf"],"wt_shell_sdf":wf["shell_sdf"],
            "wt_shell_valid":wf["shell_valid"],"wt_shell_owned":wf["shell_owned"],"wt_chemistry":wf["chemistry"],
            "wt_normals":wf['projected_normals'],"wt_local_area":wf['local_ses_area'],
            "wt_atomic_inside":wf['atomic_inside'],
            "core_rmsd_a":np.array(audit["core_rmsd_a"]),"rotation":np.array(audit["rotation"]),
            "translation":np.array(audit["translation"]),
            "topology_changed":np.array(topology_changed),
            "local_channel_qc":np.array(local_qc),
            "ses_area_delta_a2":np.array(mf["ses_total_area_a2"]-wf["ses_total_area_a2"])}
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**result)
    write_json(manifest_path,{'contract':contract,'sha256':sha256_file(path)})
    return result


def null_differences(runs):
    """Consecutive matched WT references; dependencies remain explicit."""
    output=[]
    for left,right in zip(runs,runs[1:]):
        def wt(row):
            return {"sdf":row["wt_sdf"],"shell_sdf":row["wt_shell_sdf"],
                    "shell_valid":row["wt_shell_valid"],"shell_owned":row["wt_shell_owned"],"chemistry":row["wt_chemistry"]}
        def geometry(row):
            value=wt(row)
            if 'wt_normals' in row: value.update(projected_normals=row['wt_normals'],local_ses_area=row['wt_local_area'])
            if 'wt_atomic_inside' in row: value['atomic_inside']=row['wt_atomic_inside']
            return value
        output.append(difference(geometry(left),geometry(right)))
    return output
