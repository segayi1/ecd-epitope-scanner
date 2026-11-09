"""External predictor adapters used by the opt-in V4 backend."""

from .base import ExternalPredictor, ToolError, cache_key, load_cached_result, save_cached_result
from .discotope_adapter import DiscoTopeAdapter
from .foldx_adapter import FoldXAdapter
from .procaliper_adapter import ProCaliperAdapter
from .process_runner import JsonWorkerRunner
from .spurs_adapter import SPURSAdapter

__all__ = [
    "ExternalPredictor",
    "ToolError",
    "JsonWorkerRunner",
    "cache_key",
    "load_cached_result",
    "save_cached_result",
    "ProCaliperAdapter",
    "SPURSAdapter",
    "DiscoTopeAdapter",
    "FoldXAdapter",
]
