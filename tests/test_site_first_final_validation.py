import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_site_first_final_validation import run


def test_final_run_requires_completed_baseline(tmp_path):
    (tmp_path / 'baseline_report.json').write_text(json.dumps({'status': 'NEEDS_REVIEW'}))
    with pytest.raises(ValueError, match='did not complete'):
        run(tmp_path, tmp_path / 'output', execute=lambda **kwargs: pytest.fail('model must not start'))
