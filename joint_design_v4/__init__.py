"""Opt-in Joint V4-ML design backend.

The package is dependency-light by design.  Heavy external predictors are
isolated behind JSON adapters in :mod:`ecd_epitope_scanner.integrations`.
"""

from .engine import run_joint_v4_ml
from .models import JointV4Result, ResidueKey

__all__ = ["JointV4Result", "ResidueKey", "run_joint_v4_ml"]
