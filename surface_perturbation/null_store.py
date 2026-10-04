"""Protocol-bound calibration files; null cost is separate from production."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
import numpy as np
from ..util import write_json, sha256_file
from .fields import null_differences, aligned_atoms, spatial_field, reference_id, difference, IMPLEMENTATION_ID
from .chemistry import CHANNELS, chemical_field
from .calibration import calibrate_channel
from .tools import fingerprint


GEOMETRY_CHANNELS=('shape_sdf_a','occupancy_gain_a','occupancy_loss_a','surface_normal_x','surface_normal_y','surface_normal_z','local_ses_area_a2','atom_occupancy_gain_a','atom_occupancy_loss_a')
FIELD_CHANNELS=(*GEOMETRY_CHANNELS,*CHANNELS)


def protocol_identity(pairs,reference,tools, *, implementation_id=None):
    return fingerprint({"reference":reference_id(reference),"input":pairs[0]["input_sha256"],
                        "implementation_sha256":implementation_id or IMPLEMENTATION_ID,
                        "null_protocol_sha256":hashlib.sha256(b''.join((Path(__file__).parent/name).read_bytes()
                            for name in ('null_store.py','calibration.py'))).hexdigest(),
                        "mutations":sorted(pairs[0]["mutations"]),"nanoshaper":tools.identity(tools.nano),
                        "foldx_binary":pairs[0].get("binary_sha256"),
                        "normalized_input":pairs[0].get("normalized_input_sha256"),
                        "version":4,"repacking_scope":"candidate_matched_WT","scale":2.,"shell_step":1.,"local_query_a":7.})


def numerical_bounds(pair,native_atoms,records,reference,tools,run):
    wt,mut,_=aligned_atoms(pair,native_atoms,records)
    active=np.isfinite(run['wt_sdf'])
    fine=spatial_field(wt,reference,tools,scale=3.,shell_step=.5,query_mask=active)
    fine_mut=spatial_field(mut,reference,tools,scale=3.,shell_step=.5,query_mask=active)
    fine_delta=difference(fine,fine_mut)
    bounds={name:np.abs(fine_delta[name]-run[name]) for name in GEOMETRY_CHANNELS}
    rounded={k:(np.round(v[0],6),v[1],v[2]) for k,v in wt.items()}
    chem=chemical_field(rounded,reference["points"]+1.5*reference["normals"],reference["normals"])
    rounded_mut={k:(np.round(v[0],6),v[1],v[2]) for k,v in mut.items()}
    chem_mut=chemical_field(rounded_mut,reference["points"]+1.5*reference["normals"],reference["normals"])
    # A second WT-side projection uses normals from the finer matched WT SES.
    # Its discrepancy is a conservative mapping sensitivity envelope, not an
    # assertion that machine rounding alone measures chemical projection error.
    projected=fine['projected_normals']
    projection_valid=np.isfinite(projected).all(axis=1)
    alternate_queries=reference['points'][projection_valid]+1.5*projected[projection_valid]
    alternate_wt=chemical_field(wt,alternate_queries,projected[projection_valid])
    alternate_mut=chemical_field(mut,alternate_queries,projected[projection_valid])
    alternate=np.full((len(reference['points']),len(CHANNELS)),np.nan)
    alternate[projection_valid]=alternate_mut-alternate_wt
    for i,name in enumerate(CHANNELS):
        rounding=np.abs((chem_mut[:,i]-chem[:,i])-run[name])+np.finfo(float).eps*64
        bounds[name]=np.maximum(rounding,np.abs(alternate[:,i]-run[name]))
    return bounds


def save_calibration(folder,pairs,runs,reference,tools,numerical, *, required_mask=None):
    protocol=protocol_identity(pairs,reference,tools)
    null=null_differences(runs)
    summaries={}; arrays={}
    for channel in FIELD_CHANNELS:
        blocks=np.array([r[channel] for r in null])
        value=calibrate_channel(blocks,numerical[channel],protocol,channel,
                                required_mask=required_mask,area=reference['area'])
        summaries[channel]={k:v for k,v in value.items() if not isinstance(v,np.ndarray)}
        for name,data in value.items():
            if isinstance(data,np.ndarray): arrays[f"{channel}__{name}"]=data
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    path=folder/(protocol+".npz"); np.savez_compressed(path,**arrays)
    write_json(folder/(protocol+".json"),{"protocol_id":protocol,"reference_id":reference_id(reference),
               "channels":summaries,"runs":len(runs),"array_sha256":sha256_file(path),
               "pairs":[{"wt_sha256":p["wt_sha256"],"mut_sha256":p["mut_sha256"],"pdb_list_sha256":p["pdb_list_sha256"]} for p in pairs],
               "dependencies":"consecutive WT differences share references; not independent biological samples"})
    return load_calibration(folder,protocol)


def load_calibration(folder,protocol):
    path=Path(folder)/(protocol+".json")
    if not path.exists(): return {}
    manifest=json.loads(path.read_text(encoding="utf-8"))
    data_path=path.with_suffix(".npz")
    if manifest["protocol_id"]!=protocol or sha256_file(data_path)!=manifest["array_sha256"]:
        raise ValueError("calibration fingerprint mismatch")
    with np.load(data_path,allow_pickle=False) as data:
        return {channel:{**summary,**{k.split("__",1)[1]:data[k] for k in data.files if k.startswith(channel+"__")}}
                for channel,summary in manifest["channels"].items()}
