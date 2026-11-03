import json
import pytest
from ecd_epitope_scanner.util import write_json


def test_stream_json_atomic_on_serialization_failure(tmp_path):
    path = tmp_path / 'value.json'
    value = {'unicode': '证据', 'rows': list(range(1000)), 'float': 0.123456789}
    write_json(path, value)
    before = path.read_bytes()
    assert json.loads(before) == value
    with pytest.raises(TypeError):
        write_json(path, {'partial': list(range(1000)), 'bad': object()})
    assert path.read_bytes() == before
    assert list(tmp_path.glob('*.tmp')) == []
