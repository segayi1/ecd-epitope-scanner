"""Surface-geodesic footprint generation (V2).

The package deliberately keeps molecular surface geometry separate from
mutation eligibility.  Protected residues therefore remain part of the
surface and of footprints, while the legacy adapter exposes only mutable
members to the design engine.
"""

from .engine import run_surface_patch_v2
from .models import Footprint, SurfaceGeometry, SurfaceRegion

__all__ = ["Footprint", "SurfaceGeometry", "SurfaceRegion", "run_surface_patch_v2"]
