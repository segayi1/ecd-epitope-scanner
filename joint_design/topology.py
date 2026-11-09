from __future__ import annotations

import numpy as np

from ..models import ProteinContext, ResidueRecord
from .models import TopologyContext
from .structure_audit import residue_id


def build_topology_context(
    records: list[ResidueRecord],
    protein: ProteinContext,
    tm_anchor_xyz: tuple[float, float, float] | None = None,
    tm_axis_xyz: tuple[float, float, float] | None = None,
) -> TopologyContext:
    ordered = sorted((item for item in records if item.ca is not None), key=lambda item: item.uniprot_pos)
    if not ordered:
        return TopologyContext(
            mode="topology_only", anchor_residue=None, anchor_point_xyz=None,
            relative_axis_xyz=None, membrane_plane=None, uncertainty="not_estimable",
            assumptions=["No mapped C-alpha coordinates were available."],
        )
    tail = ordered[-3:]
    anchor = np.mean(np.asarray([item.ca for item in tail], dtype=float), axis=0)
    centroid = np.mean(np.asarray([item.ca for item in ordered], dtype=float), axis=0)
    vector = centroid - anchor
    norm = float(np.linalg.norm(vector))
    axis = tuple(float(value) for value in vector / norm) if norm > 1e-8 else None
    if tm_anchor_xyz is not None and tm_axis_xyz is not None:
        return TopologyContext(
            mode="tm_geometry", anchor_residue=residue_id(tail[-1]),
            anchor_point_xyz=tuple(float(value) for value in tm_anchor_xyz),
            relative_axis_xyz=tuple(float(value) for value in tm_axis_xyz),
            membrane_plane={"method": "observed_tm_coordinates", "status": "provided_by_pipeline"},
            uncertainty="medium",
            assumptions=["A membrane geometry mode is used only when TM coordinates were observed."],
        )
    return TopologyContext(
        mode="topology_only", anchor_residue=residue_id(tail[-1]),
        anchor_point_xyz=tuple(float(value) for value in anchor), relative_axis_xyz=axis,
        membrane_plane=None, uncertainty="high",
        assumptions=[
            "The C-terminal ECD anchor follows the construct topology.",
            "The isolated ECD does not define a unique membrane plane.",
            "Relative anchor-proximal/distal descriptions are not reach probabilities.",
        ],
    )
