import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from continue_site_first_baseline import run


def test_rejects_failed_reexport_before_model_dispatch(tmp_path):
    (tmp_path / 'derived_export_report.json').write_text(json.dumps({'status': 'NEEDS_REVIEW', 'model_calls': 0}))
    with pytest.raises(ValueError, match='has not passed'):
        run(tmp_path, execute=lambda **kwargs: pytest.fail('model must not start'))
