from __future__ import annotations

import numpy as np

from .models import Footprint, SurfaceRegion


def build_surface_regions(
    footprints: list[Footprint],
    edges=None,
    *,
    surface,
    graph,
    primary_candidates: list[Footprint] | None = None,
    parameters: dict | None = None,
) -> list[SurfaceRegion]:
    """Fixed-seed primary-scale Louvain with connected, bounded communities.

    Extent is measured between primary centers on the original surface graph.
    Display ownership is disjoint; scientific footprint memberships overlap.
    """
    from collections import defaultdict

    import networkx as nx
    from scipy.sparse.csgraph import dijkstra

    from .footprint_dedup import weighted_jaccard

    cfg = dict(
        resolution=0.65,
        maximum_center_geodesic=24.0,
        maximum_geodesic_diameter=60.0,
        domain_boundary_penalty=0.15,
        minimum_surface_area=700.0,
        minimum_shared_surface_area=40.0,
        minimum_weighted_jaccard=0.12,
    )
    cfg.update(parameters or {})
    if not footprints:
        return []
    by_center = defaultdict(list)
    for f in primary_candidates or footprints:
        by_center[f.center_vertex_id].append(f)

    def preference(f):
        valid = 8 <= len(f.residues) <= 30 and not any(
            x.startswith("ERROR") for x in f.quality_flags
        )
        scale = (
            0
            if f.generation_mode == "fixed_radius" and f.radius == 12
            else 1 if f.generation_mode == "adaptive_area" else 2
        )
        return (not valid, scale, abs(f.total_surface_area - 500), f.footprint_id)

    primary = [min(fs, key=preference) for _, fs in sorted(by_center.items())]
    distances = np.atleast_2d(
        dijkstra(graph, directed=False, indices=[f.center_vertex_id for f in primary])
    )
    cd = distances[:, [f.center_vertex_id for f in primary]]
    network = nx.Graph()
    network.add_nodes_from(range(len(primary)))
    for i, left in enumerate(primary):
        for j in range(i + 1, len(primary)):
            right = primary[j]
            if not np.isfinite(cd[i, j]) or cd[i, j] > cfg["maximum_center_geodesic"]:
                continue
            shared = float(
                surface.vertex_area[
                    sorted(set(left.vertex_ids) & set(right.vertex_ids))
                ].sum()
            )
            overlap = weighted_jaccard(left, right)
            if shared <= 0 or (
                shared < cfg["minimum_shared_surface_area"]
                and overlap < cfg["minimum_weighted_jaccard"]
            ):
                continue
            weight = (
                0.5 * overlap
                + 0.3 * np.exp(-cd[i, j] / 12)
                + 0.2
                * shared
                / max(1, min(left.total_surface_area, right.total_surface_area))
            )
            if not set(left.structural_domains) & set(right.structural_domains):
                weight *= 1 - cfg["domain_boundary_penalty"]
            network.add_edge(i, j, weight=float(weight))

    def extent(group):
        ids = sorted(group)
        return float(cd[np.ix_(ids, ids)].max())

    def split(group):
        ids = sorted(group)
        if len(ids) < 2 or extent(group) <= cfg["maximum_geodesic_diameter"]:
            return [set(group)]
        a, b = np.unravel_index(np.argmax(cd[np.ix_(ids, ids)]), (len(ids), len(ids)))
        left = {i for i in ids if cd[ids[a], i] <= cd[ids[b], i]}
        return [
            g
            for part in (left, set(ids) - left)
            for component in nx.connected_components(network.subgraph(part))
            for g in split(component)
        ]

    groups = []
    for component in nx.connected_components(network):
        sub = network.subgraph(sorted(component)).copy()
        communities = (
            nx.community.louvain_communities(
                sub, weight="weight", resolution=cfg["resolution"], seed=0
            )
            if sub.number_of_edges()
            else [component]
        )
        for community in communities:
            for connected in nx.connected_components(network.subgraph(community)):
                groups.extend(split(connected))

    def area(group):
        return float(
            surface.vertex_area[
                sorted({v for i in group for v in primary[i].vertex_ids})
            ].sum()
        )

    while True:
        options = []
        for i, group in enumerate(groups):
            if area(group) >= cfg["minimum_surface_area"] and len(group) > 1:
                continue
            for j, other in enumerate(groups):
                if i == j or extent(group | other) > cfg["maximum_geodesic_diameter"]:
                    continue
                boundary = sum(
                    network[a][b]["weight"]
                    for a in sorted(group)
                    for b in sorted(other)
                    if network.has_edge(a, b)
                )
                if boundary:
                    options.append((-boundary, min(group), min(other), i, j))
        if not options:
            break
        _, _, _, i, j = min(options)
        groups[i] |= groups[j]
        del groups[j]
    groups.sort(
        key=lambda g: (min(primary[i].center_residue_id or 0 for i in g), min(g))
    )
    assignment = {
        primary[i].center_vertex_id: k for k, g in enumerate(groups) for i in g
    }
    members = [{} for _ in groups]
    for f in footprints:
        k = assignment[f.center_vertex_id]
        members[k][f.footprint_id] = f
        f.region_id = f"SR{k+1:03d}"
    nearest = np.argmin(distances, axis=0)
    reachable = np.isfinite(np.min(distances, axis=0))
    contributions = defaultdict(lambda: np.zeros(len(groups)))
    for v in np.flatnonzero(reachable):
        k = assignment[primary[nearest[v]].center_vertex_id]
        contributions[int(surface.owner_residue[v])][k] += surface.vertex_area[v]
    display = [[] for _ in groups]
    for residue, values in sorted(contributions.items()):
        display[int(np.argmax(values))].append(residue)
    output = []
    for k, group in enumerate(groups):
        ids = sorted(group)

        def representative_key(i):
            f = primary[i]
            return (
                len(f.mutable_residues) < 5,
                not 8 <= len(f.residues) <= 30,
                float(cd[i, ids].mean()),
                abs(f.total_surface_area - 500),
                f.footprint_id,
            )

        rep = primary[min(ids, key=representative_key)]
        for i in ids:
            f = primary[i]
            members[k][f.footprint_id] = f
            f.region_id = f"SR{k+1:03d}"
        vertices = sorted({v for f in members[k].values() for v in f.vertex_ids})
        weights = surface.vertex_area[vertices]
        centroid = np.average(surface.vertices[vertices], axis=0, weights=weights)
        output.append(
            SurfaceRegion(
                region_id=f"SR{k+1:03d}",
                footprint_ids=sorted(members[k]),
                representative_footprint_id=rep.footprint_id,
                residues_union=sorted(
                    {r for f in members[k].values() for r in f.residues}
                ),
                surface_area_union=float(weights.sum()),
                centroid_xyz=tuple(map(float, centroid)),
                domain_ids=sorted(
                    {d for f in members[k].values() for d in f.structural_domains}
                ),
                quality_flags=["WARNING_SINGLETON_REGION"] if len(group) == 1 else [],
                display_residues=display[k],
                primary_footprint_ids=[primary[i].footprint_id for i in ids],
                geodesic_extent=extent(group),
            )
        )
    return output
