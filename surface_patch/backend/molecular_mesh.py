"""Optional molecular-mesh backend boundary.

No mesher is bundled because the workbench must remain installable without a
platform-specific MSMS/EDTSurf binary.  This explicit backend boundary lets a
deployment provide one later; the audited SR point-cloud backend is the V2
fallback and is used by default.
"""
from __future__ import annotations


class MolecularMeshBackend:
    def build(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        raise RuntimeError(
            "molecular_mesh backend is not installed; use sr_pointcloud or configure a mesher"
        )
