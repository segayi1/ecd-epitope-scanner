"""Spatial atom-type support and geometric directions, not contact energies."""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

CHANNELS = ("nonpolar_support", "polar_support", "donor_support", "acceptor_support",
            "donor_x", "donor_y", "donor_z", "acceptor_x", "acceptor_y", "acceptor_z",
            "aromatic_support", "aromatic_x", "aromatic_y", "aromatic_z")
DONORS = {"ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
          "HIS": {"ND1", "NE2"}, "LYS": {"NZ"}, "SER": {"OG"}, "THR": {"OG1"},
          "TRP": {"NE1"}, "TYR": {"OH"}, "CYS": {"SG"}}
ACCEPTORS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"},
             "GLN": {"OE1"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
             "MET": {"SD"}}
RINGS = {"PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
         "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
         "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
         "HIS": ("CG", "ND1", "CD2", "CE1", "NE2")}
RADII = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "P": 1.8}


def chemical_field(atom_map, points, normals, *, bandwidth=1.):
    keys = sorted(atom_map)
    xyz = np.array([atom_map[k][0] for k in keys])
    tree = cKDTree(xyz)
    types = np.zeros((len(keys), len(CHANNELS)))
    direction_valid = np.ones(len(keys), bool)
    radii = np.array([RADII.get(atom_map[k][2], 1.7) for k in keys])
    residues = {}
    for i, key in enumerate(keys):
        residues.setdefault(key[:3], []).append(i)
    for i, key in enumerate(keys):
        _, resname, element = atom_map[key]
        name = key[3]
        types[i, 0 if element == "C" else 1] = 1
        # Histidine donor/acceptor state cannot be asserted without protonation evidence.
        donor = name == "N" and resname != "PRO" or name in DONORS.get(resname, set()) and resname != "HIS"
        acceptor = name in {"O", "OXT"} or name in ACCEPTORS.get(resname, set())
        types[i, 2:4] = (donor, acceptor)
        bonded = [j for j in residues[key[:3]] if j != i and .5 < np.linalg.norm(xyz[j]-xyz[i]) < 1.95]
        if donor or acceptor:
            if bonded:
                vector = -(xyz[bonded] - xyz[i]) / np.linalg.norm(xyz[bonded]-xyz[i], axis=1)[:, None]
                vector = vector.sum(axis=0)
                norm = np.linalg.norm(vector)
                if norm > 1e-8:
                    vector /= norm
                    if donor: types[i, 4:7] = vector
                    if acceptor: types[i, 7:10] = vector
                else: direction_valid[i] = False
            else: direction_valid[i] = False
        if name in RINGS.get(resname, ()):
            ring = [j for j in residues[key[:3]] if keys[j][3] in RINGS[resname]]
            if len(ring) >= 3:
                centered = xyz[ring] - xyz[ring].mean(axis=0)
                _, _, vh = np.linalg.svd(centered)
                v = vh[-1]
                v *= 1 if v[np.argmax(np.abs(v))] >= 0 else -1
                types[i, 10] = 1
                # A plane has no signed face: project its normal onto the query-facing side below.
                types[i, 11:14] = v
    result = np.zeros((len(points), len(CHANNELS)))
    valid = np.ones((len(points), len(CHANNELS)), bool)
    for i, point in enumerate(points):
        candidates = tree.query_ball_point(point, float(radii.max()+3*bandwidth))
        for j in candidates:
            delta = point-xyz[j]; distance = np.linalg.norm(delta)
            if distance < 1e-8 or np.dot(delta, normals[i]) < 0:
                continue
            direction = delta / distance
            start = xyz[j] + direction * radii[j]
            segment = point-start
            length = np.linalg.norm(segment)
            if length > 0:
                neighbors = tree.query_ball_point((start+point)/2, length/2+radii.max())
                hidden = False
                for neighbor in neighbors:
                    if neighbor == j: continue
                    t = np.clip(np.dot(xyz[neighbor]-start, segment)/(length*length), 0, 1)
                    if np.linalg.norm(start+t*segment-xyz[neighbor]) < radii[neighbor]-1e-6:
                        hidden = True; break
                if hidden: continue
            weight = np.exp(-.5 * ((distance-radii[j])/bandwidth)**2)
            if atom_map[keys[j]][1]=='HIS' and keys[j][3] in {'ND1','NE2'} and weight>.01:
                valid[i,2:10]=False
            feature = types[j].copy()
            if np.dot(feature[11:14], normals[i]) < 0: feature[11:14] *= -1
            result[i] += weight*feature
            if not direction_valid[j] and weight > .01:
                valid[i, 4:10] = False
    result[~valid] = np.nan
    return result
