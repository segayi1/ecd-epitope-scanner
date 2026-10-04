"""Native SES quadrature with conserved areas and explicit residue identities."""
from __future__ import annotations

from collections import defaultdict
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def canonical_reference(mesh, atom_map, records, eligible_positions, designable_positions, spacing=2.):
    keys = sorted(atom_map)
    xyz = np.array([atom_map[k][0] for k in keys])
    tree = cKDTree(xyz)
    _, nearest = tree.query(mesh.centers)
    identities = {(r.chain_id, int(r.auth_seq_id), (r.insertion_code or "").strip()): r.uniprot_pos for r in records}
    owners = np.array([identities.get(keys[i][:3], -1) for i in nearest], dtype=int)
    # Mesh component identity and external cavity audit must precede aggregation.
    edge_faces = defaultdict(list)
    for face_id, face in enumerate(mesh.faces):
        for a, b in ((face[0],face[1]),(face[1],face[2]),(face[2],face[0])):
            edge_faces[tuple(sorted((a,b)))].append(face_id)
    neighbors = []
    for ids in edge_faces.values():
        if len(ids) == 2: neighbors.append(ids)
    edge = np.array(neighbors, dtype=int).reshape(-1,2)
    topology = coo_matrix((np.ones(2*len(edge)), (np.r_[edge[:,0],edge[:,1]],np.r_[edge[:,1],edge[:,0]])), shape=(len(owners),)*2).tocsr()
    count, component = connected_components(topology, directed=False)
    # A detached component enclosed by another closed component is an internal
    # cavity boundary. Determine containment geometrically, not by total area.
    external = np.ones(len(owners), bool)
    if count > 1:
        from .geometry import Mesh
        components = []
        for label in range(count):
            ids = np.flatnonzero(component == label)
            components.append((ids, Mesh(mesh.vertices, mesh.faces[ids])))
        for label, (ids, submesh) in enumerate(components):
            probe = submesh.centers[0]
            for other, (_, enclosing) in enumerate(components):
                if other == label: continue
                sdf, valid = enclosing.signed_distance([probe])
                if valid[0] and sdf[0] < 0:
                    external[ids] = False; break
    cell = np.floor(mesh.centers/spacing).astype(int)
    bins = np.where(mesh.normals >= 0, 1, 0)
    group_keys = [(int(owners[i]), int(component[i]), *cell[i], *bins[i]) for i in range(len(owners))]
    # Split disconnected fragments within each cell before area aggregation.
    parent = np.arange(len(owners))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for a,b in neighbors:
        if group_keys[a] == group_keys[b]:
            left,right = find(a),find(b)
            parent[max(left,right)] = min(left,right)
    groups = defaultdict(list)
    for i in range(len(owners)): groups[find(i)].append(i)
    face_to_sample = np.empty(len(owners), int)
    points=[]; normals=[]; areas=[]; owner=[]; reliable=[]; designed=[]; domain=[]
    domains = {r.uniprot_pos:r.domain for r in records}
    for ids in groups.values():
        index=len(points); face_to_sample[ids]=index
        weights=mesh.area[ids]; total=weights.sum()
        center=np.average(mesh.centers[ids], axis=0, weights=weights)
        representative=ids[np.argmin(np.linalg.norm(mesh.centers[ids]-center,axis=1))]
        normal=np.average(mesh.normals[ids],axis=0,weights=weights)
        normal /= max(np.linalg.norm(normal),1e-30)
        pos=int(owners[representative])
        points.append(mesh.centers[representative]); normals.append(normal); areas.append(total); owner.append(pos)
        reliable.append(pos in eligible_positions and bool(external[representative]))
        designed.append(pos in designable_positions); domain.append(domains.get(pos,"UNMAPPED"))
    points=np.array(points); normals=np.array(normals); area=np.array(areas)
    reliable=np.array(reliable); designed=np.array(designed)&reliable
    mapped=face_to_sample[edge]
    ids=np.flatnonzero(reliable); remap=np.full(len(points),-1); remap[ids]=np.arange(len(ids))
    keep=(mapped[:,0]!=mapped[:,1]) & reliable[mapped[:,0]] & reliable[mapped[:,1]]
    links=remap[mapped[keep]]
    lengths=np.linalg.norm(points[mapped[keep,0]]-points[mapped[keep,1]],axis=1)
    # Duplicate adjacency paths use minimum distance rather than adding weights.
    unique={}
    for (a,b),length in zip(links,lengths):
        key=tuple(sorted((int(a),int(b)))); unique[key]=min(unique.get(key,np.inf),float(length))
    links=np.array(list(unique),int).reshape(-1,2); lengths=np.array(list(unique.values()))
    graph=coo_matrix((np.r_[lengths,lengths],(np.r_[links[:,0],links[:,1]],np.r_[links[:,1],links[:,0]])),shape=(len(ids),)*2).tocsr()
    return {"points":points[ids],"normals":normals[ids],"area":area[ids],
            "owner":np.array(owner)[ids],"designable":designed[ids],"graph":graph,
            "domain":np.array(domain)[ids],"spacing":spacing,
            "excluded_area_a2":float(area[~reliable].sum()),
            "internal_cavity_area_a2":float(mesh.area[~external].sum()),
            "mesh_total_area_a2":float(mesh.area.sum()),"basis":"SES_NATIVE"}
