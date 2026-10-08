"""Local/WSL SES and electrostatics with content-addressed physical settings."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
import re
from collections import OrderedDict
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from ..util import write_json, sha256_file
from ..integrations.wsl_bridge import _wsl_path, _decode_output
from .geometry import Mesh
from .chemistry import RADII


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


class PhysicalTools:
    def __init__(self, cache, config=None, budget=None):
        self.cache = Path(cache); self.cache.mkdir(parents=True, exist_ok=True)
        self.config = config or {}; self.budget = budget
        self.distro = self.config.get("wsl_distro", "Ubuntu-22.04")
        self.use_wsl = self.config.get("use_wsl", os.name == "nt")
        self.nano = self.config.get("nanoshaper", os.environ.get("ECD_NANOSHAPER", "/home/segayi1/ecd_surface_tools/nanoshaper/build/NanoShaper"))
        self.apbs = self.config.get("apbs", "/usr/bin/apbs")
        self.pqr = self.config.get("pdb2pqr", "/usr/bin/pdb2pqr")
        self.calls = {"nanoshaper": 0, "apbs": 0, "pdb2pqr": 0, "cache_hits": 0}
        self.identities = {}
        self.activity=[]
        self.mesh_objects = OrderedDict()

    def path(self, value):
        return _wsl_path(str(Path(value).resolve()), self.distro) if self.use_wsl else str(Path(value).resolve())

    def run(self, argv, cwd, timeout=300):
        if self.budget: timeout = min(timeout, self.budget.remaining())
        command = ["wsl.exe", "-d", self.distro, "--cd", self.path(cwd), "--exec", *argv] if self.use_wsl else argv
        started=time.perf_counter()
        try:
            completed = subprocess.run(command, cwd=cwd, capture_output=True, timeout=timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            self.activity.append({'tool':argv[0],'runtime_seconds':time.perf_counter()-started,'status':'TIMEOUT'})
            raise TimeoutError('physical worker timeout: '+argv[0]) from exc
        text = _decode_output(completed.stdout or b"") + "\n" + _decode_output(completed.stderr or b"")
        peak=re.search(r'([\d.]+) MB high water',text)
        self.activity.append({'tool':argv[0],'runtime_seconds':time.perf_counter()-started,
                              'exit_code':completed.returncode,'peak_solver_mb':float(peak.group(1)) if peak else None})
        if completed.returncode:
            raise RuntimeError(f"physical worker exit {completed.returncode}: {text[-1800:]}")
        return text

    def identity(self, executable):
        if executable not in self.identities:
            if self.use_wsl:
                value = self.run(["sha256sum", executable], self.cache).split()[0]
            else: value = sha256_file(Path(executable))
            if len(value) != 64: raise RuntimeError("physical tool digest unavailable")
            self.identities[executable] = value
        return self.identities[executable]

    def ses(self, atom_map, *, scale=2.):
        from .case_resources import ACTIVE_STORAGE
        if ACTIVE_STORAGE.get() is not None:
            ACTIVE_STORAGE.get().check(reserved=128*1024**2)
        keys = sorted(atom_map)
        rows = np.array([[*atom_map[k][0], RADII.get(atom_map[k][2], 1.7)] for k in keys])
        contract = {"xyzr_sha256": hashlib.sha256(rows.tobytes()).hexdigest(),
                    "binary_sha256": self.identity(self.nano), "scale": scale,
                    "probe_radius_a": 1.4, "radius_policy": RADII, "smoothing": False,
                    "surface": "ses", "cavity_policy": "preserve_and_audit"}
        folder = self.cache / ("ses_" + fingerprint(contract))
        folder.mkdir(exist_ok=True)
        mesh_path = folder / "triangulatedSurf.off"
        if mesh_path.exists() and (folder/"complete.json").exists():
            self.calls["cache_hits"] += 1
        else:
            np.savetxt(folder / "atoms.xyzr", rows, fmt="%.9f")
            settings = {
                "Grid_scale": scale, "Grid_perfil": 80., "XYZR_FileName": "atoms.xyzr",
                "Surface": "ses", "Probe_Radius": 1.4, "Build_epsilon_maps": "false",
                "Build_status_map": "true", "Triangulation": "true", "Accurate_Triangulation": "true",
                "Smooth_Mesh": "false", "Compute_Vertex_Normals": "true", "Number_thread": 4,
                "Cavity_Detection_Filling": "false", "Save_Status_map": "false",
                "Save_Mesh_MSMS_Format": "false", "Surface_File_Name": "triangulatedSurf.off",
                "Patch_Based_Algorithm": "true"}
            (folder/"surface.prm").write_text("\n".join(f"{k} = {v}" for k,v in settings.items())+"\n", encoding="ascii")
            self.calls["nanoshaper"] += 1
            log = self.run([self.nano, "surface.prm"], folder)
            (folder/"nanoshaper.log").write_text(log, encoding="utf-8")
            if not mesh_path.exists():
                matches = sorted(folder.glob("*.off"))
                if len(matches) != 1: raise RuntimeError("NanoShaper did not produce an unambiguous OFF mesh")
                mesh_path = matches[0]
            mesh = read_off(mesh_path)
            if not mesh.closed: raise RuntimeError("SES mesh is not closed")
            write_json(folder/"complete.json", {**contract, "mesh_sha256": sha256_file(mesh_path), "mesh": mesh_path.name})
        manifest = json.loads((folder/"complete.json").read_text(encoding="utf-8"))
        mesh_path = folder / manifest["mesh"]
        if sha256_file(mesh_path) != manifest["mesh_sha256"]:
            raise RuntimeError("cached SES fingerprint changed")
        identity = manifest["mesh_sha256"]
        if identity in self.mesh_objects:
            self.mesh_objects.move_to_end(identity)
            return self.mesh_objects[identity]
        mesh = read_off(mesh_path)
        self.mesh_objects[identity] = mesh
        while len(self.mesh_objects) > max(1, int(self.config.get("mesh_object_cache_size", 4))):
            self.mesh_objects.popitem(last=False)
        return mesh

    def prepare_pqr(self, pdb_path, *, protocol="legacy_optimized"):
        from .case_resources import ACTIVE_STORAGE
        if ACTIVE_STORAGE.get() is not None:
            ACTIVE_STORAGE.get().check(reserved=max(2*1024**2,Path(pdb_path).stat().st_size*4))
        if protocol not in {"legacy_optimized", "heavy_atom_preserving_v1"}:
            raise ValueError(f"unknown PDB2PQR preparation protocol: {protocol}")
        support_paths=self.config.get('pdb2pqr_support_paths')
        if support_paths is None and self.pqr=='/usr/bin/pdb2pqr' and not self.config.get('pdb2pqr_version3'):
            support_paths=['/usr/share/pdb2pqr','/usr/lib/pdb2pqr']
        support_hash=None
        if support_paths:
            key='pqr_support:'+fingerprint(support_paths)
            if key not in self.identities:
                script=("import pathlib,hashlib,sys;h=hashlib.sha256(sys.version.encode());"
                        f"roots={support_paths!r};"
                        "files=sorted({p for r in roots for p in pathlib.Path(r).rglob('*') if p.is_file() and '__pycache__' not in p.parts});"
                        "[(h.update(str(p).encode()),h.update(p.read_bytes())) for p in files];print(h.hexdigest() if files else 'UNAVAILABLE')")
                support_hash=self.run([self.config.get('pdb2pqr_support_python','/usr/bin/python3'),'-c',script],self.cache).splitlines()[0].strip()
                self.identities[key]=support_hash if len(support_hash)==64 else None
            support_hash=self.identities[key]
        preserve_heavy = protocol == "heavy_atom_preserving_v1"
        contract = {"input_sha256": sha256_file(pdb_path), "binary_sha256": self.identity(self.pqr),
                    "support_files_sha256":support_hash,
                    "force_field": "PARSE", "pH": 7.4,
                    "optimize_hydrogens": not preserve_heavy,
                    "debump": not preserve_heavy, "protocol": protocol}
        folder = self.cache / ("pqr_"+fingerprint(contract)); folder.mkdir(exist_ok=True)
        output = folder / "model.pqr"
        if support_hash and output.exists() and (folder/"complete.json").exists():
            manifest = json.loads((folder/"complete.json").read_text())
            if sha256_file(output) != manifest["pqr_sha256"]: raise RuntimeError("PQR cache fingerprint changed")
            self.calls["cache_hits"] += 1
        else:
            # Ubuntu 22.04 packages PDB2PQR 2.x; callers can explicitly configure 3.x.
            version3 = bool(self.config.get("pdb2pqr_version3", False))
            flags = ["--ff=PARSE", "--with-ph=7.4", "--keep-chain" if version3 else "--chain"]
            if version3: flags += ["--titration-state-method=propka"]
            if preserve_heavy: flags += ["--nodebump", "--noopt"]
            self.calls["pdb2pqr"] += 1
            log = self.run([self.pqr, *flags, self.path(pdb_path), self.path(output)], folder)
            (folder/"pdb2pqr.log").write_text(log, encoding="utf-8")
            if not output.exists(): raise RuntimeError("PDB2PQR output missing")
            write_json(folder/"complete.json", {**contract, "pqr_sha256": sha256_file(output)})
        return output, read_pqr(output)

    def potential(self, pqr_path, grid):
        from .case_resources import ACTIVE_STORAGE
        if ACTIVE_STORAGE.get() is not None:
            ACTIVE_STORAGE.get().check(reserved=int(np.prod(grid['counts']))*24+Path(pqr_path).stat().st_size)
        contract = {"pqr_sha256": sha256_file(pqr_path), "binary_sha256": self.identity(self.apbs),
                    "grid": grid, "pdie": 2., "sdie": 78.54, "temperature": 298.15,
                    "ions": [[1,.15,2.],[-1,.15,2.]], "equation": "lpbe", "boundary": "sdh",
                    "surface": "smol", "srad": 1.4, "swin": .3, "charge_map": "spl2"}
        folder = self.cache / ("apbs_"+fingerprint(contract)); folder.mkdir(exist_ok=True)
        output = folder / "potential.dx"
        if (folder/"complete.json").exists():
            manifest = json.loads((folder/"complete.json").read_text())
            output=folder/manifest.get('dx_name','potential.dx')
            if sha256_file(output) != manifest["dx_sha256"]: raise RuntimeError("APBS cache fingerprint changed")
            self.calls["cache_hits"] += 1
        else:
            counts = " ".join(map(str, grid["counts"]))
            length = " ".join(f"{x:.8f}" for x in grid["lengths"])
            center = " ".join(f"{x:.8f}" for x in grid["center"])
            # Copy to a fixed basename: APBS input parsing must not interpret path spaces.
            import shutil
            shutil.copyfile(pqr_path, folder/"input.pqr")
            text = (f"read\n mol pqr input.pqr\nend\nelec\n mg-manual\n dime {counts}\n glen {length}\n gcent {center}\n"
                    " mol 1\n lpbe\n bcfl sdh\n pdie 2.0\n sdie 78.54\n chgm spl2\n srfm smol\n srad 1.4\n swin 0.3\n sdens 10.0\n"
                    " temp 298.15\n ion charge 1 conc 0.15 radius 2.0\n ion charge -1 conc 0.15 radius 2.0\n"
                    " calcenergy no\n calcforce no\n write pot dx potential\nend\nquit\n")
            (folder/"apbs.in").write_text(text, encoding="ascii")
            self.calls["apbs"] += 1
            log = self.run([self.apbs, "apbs.in"], folder, timeout=600)
            (folder/"apbs.log").write_text(log, encoding="utf-8")
            if "not converg" in log.lower() or "error:" in log.lower() or "fatal" in log.lower():
                raise RuntimeError("APBS convergence/error diagnostic")
            if not output.exists():
                maps=list(folder.glob('potential*.dx'))
                if len(maps)!=1: raise RuntimeError("APBS potential output missing or distributed maps ambiguous")
                output=maps[0]
            write_json(folder/"complete.json", {**contract, "dx_name":output.name,"dx_sha256": sha256_file(output)})
        field=read_dx(output)
        if list(map(len,field.grid))!=grid['counts']:
            raise ValueError('APBS actual grid counts differ from fixed contract')
        for axis,origin,spacing in zip(field.grid,grid['origin'],grid['spacing']):
            if not np.isclose(axis[0],origin,atol=1e-6) or not np.allclose(np.diff(axis),spacing,atol=1e-6):
                raise ValueError('APBS actual origin/spacing differ from fixed contract')
        return field


def read_off(path):
    lines = [line.strip() for line in Path(path).read_text().splitlines() if line.strip() and not line.startswith("#")]
    if lines[0].split()[0] not in {"OFF", "NOFF", "COFF", "OFF+N", "OFF+A", "OFF+N+A"}: raise ValueError("unsupported mesh header")
    nvertices, nfaces, _ = map(int, lines[1].split()[:3])
    vertices = np.array([[float(x) for x in line.split()[:3]] for line in lines[2:2+nvertices]])
    faces = []
    for line in lines[2+nvertices:2+nvertices+nfaces]:
        row = list(map(int, line.split()[:4]))
        if row[0] != 3: raise ValueError("non-triangle mesh")
        faces.append(row[1:])
    # NanoShaper may repeat EXACT coordinates and emit zero-area triangles at
    # grid coincidences. Merge exact duplicates only; never bridge an opening.
    unique, inverse = np.unique(vertices, axis=0, return_inverse=True)
    faces = inverse[np.asarray(faces, dtype=int)]
    triangles = unique[faces]
    area = np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],
                                  triangles[:,2]-triangles[:,0]), axis=1)/2
    keep = area > 0
    mesh = Mesh(unique, faces[keep])
    mesh.import_audit = {"exact_duplicate_vertices": len(vertices)-len(unique),
                         "zero_area_faces_removed": int((~keep).sum()),
                         "geometry_displacement_a": 0.}
    return mesh


def read_pqr(path):
    records = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] not in {"ATOM", "HETATM"}: continue
        if len(fields) == 11: chain, res = fields[4], fields[5]
        elif len(fields) == 10: chain, res = "", fields[4]
        else: raise ValueError("ambiguous PQR record")
        key = (chain, int(res), fields[2])
        if key in records: raise ValueError("duplicate PQR atom")
        records[key] = {"residue": fields[3], "xyz": list(map(float, fields[-5:-2])),
                        "charge": float(fields[-2]), "radius": float(fields[-1])}
    if not records: raise ValueError("empty PQR")
    return records


def make_grid(atom_xyz, spacing=.5, margin=16.):
    xyz = np.asarray(atom_xyz)
    extent = xyz.max(axis=0)-xyz.min(axis=0)+2*margin
    counts = (np.ceil(extent / spacing / 32).astype(int)*32+1)
    lengths = (counts-1)*spacing
    return {"counts": counts.tolist(), "lengths": lengths.tolist(),
            "center": ((xyz.max(axis=0)+xyz.min(axis=0))/2).tolist(), "spacing": [spacing]*3,
            "origin": (((xyz.max(axis=0)+xyz.min(axis=0))/2)-lengths/2).tolist()}


def read_dx(path):
    counts = origin = None; deltas = []; data = []; collecting = False; expected = None
    for line in Path(path).read_text().splitlines():
        if "gridpositions counts" in line: counts = tuple(map(int, line.split()[-3:]))
        elif line.startswith("origin"): origin = np.array(list(map(float, line.split()[1:4])))
        elif line.startswith("delta"): deltas.append(list(map(float, line.split()[1:4])))
        elif "data follows" in line:
            expected = int(line.split("items")[1].split()[0]); collecting = True
        elif collecting and len(data) < expected:
            data.extend(float(x) for x in line.split())
    if counts is None or origin is None or len(deltas) != 3 or len(data) != np.prod(counts):
        raise ValueError("invalid OpenDX map")
    matrix = np.array(deltas)
    if not np.allclose(matrix, np.diag(np.diag(matrix))): raise ValueError("non-axis-aligned OpenDX")
    axes = [origin[i]+np.arange(counts[i])*matrix[i,i] for i in range(3)]
    return RegularGridInterpolator(axes, np.array(data).reshape(counts), bounds_error=False, fill_value=np.nan)
