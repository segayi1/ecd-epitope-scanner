import sys
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from reexport_site_first_acceptance import extract_safe


def test_extract_safe_checks_paths_before_extracting(tmp_path):
    archive = tmp_path / 'unsafe.zip'
    with zipfile.ZipFile(archive, 'w') as z:
        z.writestr('valid.txt', 'good')
        z.writestr('../escape.txt', 'bad')
    target = tmp_path / 'out'
    with pytest.raises(ValueError, match='Unsafe'):
        extract_safe(archive, target)
    assert not target.exists() and not (tmp_path / 'escape.txt').exists()
