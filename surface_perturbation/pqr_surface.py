"""NanoShaper SES using the actual per-atom radii supplied to APBS."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import numpy as np
from .tools import fingerprint,read_off
from ..util import sha256_file,write_json


def pqr_surface(tools,parameters):
    rows=[]
    for key,value in sorted(parameters.items()):
        radius=float(value['radius'])
        if radius<0 or not np.isfinite(radius):raise ValueError('invalid PQR radius')
        # Zero-radius charge sites remain in APBS, but have no sphere boundary.
        if radius>0:rows.append([*value['xyz'],radius])
    if not rows:raise ValueError('no positive-radius PQR atoms')
    rows=np.array(rows,float)
    settings={'Grid_scale':2.,'Grid_perfil':80.,'XYZR_FileName':'atoms.xyzr',
        'Surface':'ses','Probe_Radius':1.4,'Build_epsilon_maps':'false','Build_status_map':'true',
        'Triangulation':'true','Accurate_Triangulation':'true','Smooth_Mesh':'false',
        'Compute_Vertex_Normals':'true','Number_thread':4,'Cavity_Detection_Filling':'false',
        'Save_Status_map':'false','Save_Mesh_MSMS_Format':'false',
        'Surface_File_Name':'triangulatedSurf.off','Patch_Based_Algorithm':'true'}
    contract={'xyzr_sha256':hashlib.sha256(rows.tobytes()).hexdigest(),
        'binary_sha256':tools.identity(tools.nano),'settings':settings,
        'radius_policy':'actual PQR per-atom radius; zero-radius sites have no sphere'}
    folder=tools.cache/('pqr_ses_'+fingerprint(contract));folder.mkdir(exist_ok=True)
    mesh_path=folder/'triangulatedSurf.off';manifest_path=folder/'complete.json'
    if manifest_path.is_file() and mesh_path.is_file():
        manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest['contract']!=contract or manifest['mesh_sha256']!=sha256_file(mesh_path):
            raise ValueError('PQR SES cache fingerprint changed')
        tools.calls['cache_hits']+=1
    else:
        np.savetxt(folder/'atoms.xyzr',rows,fmt='%.9f')
        (folder/'surface.prm').write_text('\n'.join(f'{k} = {v}' for k,v in settings.items())+'\n',encoding='ascii')
        tools.calls['nanoshaper']+=1
        log=tools.run([tools.nano,'surface.prm'],folder)
        (folder/'nanoshaper.log').write_text(log,encoding='utf-8')
        if not mesh_path.is_file():raise RuntimeError('PQR SES output missing')
        mesh=read_off(mesh_path)
        if not mesh.closed:raise RuntimeError('PQR SES is not closed')
        write_json(manifest_path,{'contract':contract,'mesh_sha256':sha256_file(mesh_path)})
    mesh=read_off(mesh_path)
    if not mesh.closed:raise RuntimeError('PQR SES is not closed')
    from .deadlines import guard_mesh
    return guard_mesh(mesh,getattr(tools,'budget',None))
