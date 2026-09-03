#!/usr/bin/env python3
"""Compatibility launcher for command-line analysis."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ecd_epitope_scanner.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["analyze", *sys.argv[1:]]))
