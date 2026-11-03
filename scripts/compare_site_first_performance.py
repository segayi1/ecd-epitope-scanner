"""Compare completed matched runtime receipts; never starts model work."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.joint_design_v4.site_first_performance import compare_performance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    parser.add_argument('optimized', type=Path)
    args = parser.parse_args()
    baseline = json.loads(args.baseline.read_text(encoding='utf-8'))
    optimized = json.loads(args.optimized.read_text(encoding='utf-8'))
    print(json.dumps(compare_performance(baseline, optimized), indent=2))


if __name__ == '__main__':
    main()
