import json
import pytest
from ecd_epitope_scanner.scripts.continue_surface_validation import verify_frozen_source
from ecd_epitope_scanner.util import sha256_file


@pytest.mark.parametrize('changed',['source','checkout'])
def test_concurrent_lane_refuses_changed_code_without_overwriting_snapshot(tmp_path,changed):
    source=tmp_path/'source';checkout=tmp_path/'checkout'
    source.mkdir();checkout.mkdir()
    for directory in (source,checkout):(directory/'worker.py').write_text('evidence = 1\n')
    manifest={'source':str(source),'files':{'worker.py':sha256_file(source/'worker.py')}}
    (checkout/'source_overlay_manifest.json').write_text(json.dumps(manifest))
    verify_frozen_source(source,checkout)
    (source if changed=='source' else checkout).joinpath('worker.py').write_text('evidence = 2\n')
    original=(checkout/'worker.py').read_bytes()
    with pytest.raises(ValueError,match='FROZEN_SOURCE_CHANGED'):verify_frozen_source(source,checkout)
    assert (checkout/'worker.py').read_bytes()==original
