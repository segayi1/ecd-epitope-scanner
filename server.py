#!/usr/bin/env python3
"""Compatibility launcher for the FastAPI workbench."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ecd_epitope_scanner.cli import main

if __name__ == "__main__":
    raise SystemExit(main(["serve", *sys.argv[1:]]))
