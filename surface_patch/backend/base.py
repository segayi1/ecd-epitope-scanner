from __future__ import annotations

from typing import Any, Protocol

from ..models import SurfaceGeometry


class SurfaceBackend(Protocol):
    def build(
        self,
        structure: Any,
        chain_id: str,
        ecd_residues: dict[int, Any],
        residue_confidence: dict[int, float | None],
        config: dict[str, Any],
    ) -> SurfaceGeometry: ...
