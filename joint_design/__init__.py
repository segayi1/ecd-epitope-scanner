"""Joint V3 surface-probe design engine.

The package is intentionally independent from the legacy ``RegionRecord``
semantics.  The pipeline adapts the final sectors and probes back to those
records so the existing web renderer can continue to display results.
"""

from .engine import run_joint_v3

__all__ = ["run_joint_v3"]
