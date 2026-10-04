import pytest
from ecd_epitope_scanner.scripts.benchmark_surface_perturbation import frozen_worker_parameters
from ecd_epitope_scanner.util import sha256_file


def test_detached_replay_uses_source_binary_and_does_not_mutate_configuration(tmp_path):
    installed=tmp_path/'source/vendor/foldx';installed.mkdir(parents=True)
    binary=installed/'foldx_1.exe';binary.write_bytes(b'installed binary fixture')
    config={'upstream_parameters':{'foldx_confirmation':{'number_of_runs':5}}}
    parameters,audit=frozen_worker_parameters(tmp_path/'source',config)
    assert parameters['foldx_confirmation']['binary']==str(binary.resolve())
    assert audit['foldx_binary_sha256']==sha256_file(binary)
    assert config=={'upstream_parameters':{'foldx_confirmation':{'number_of_runs':5}}}
    binary.write_bytes(b'changed binary fixture')
    with pytest.raises(ValueError,match='FROZEN_DEPENDENCY_CHANGED'):
        frozen_worker_parameters(tmp_path/'source',{'foldx_binary_sha256':audit['foldx_binary_sha256']})


def test_missing_or_explicitly_misconfigured_dependency_fails_before_pipeline(tmp_path):
    with pytest.raises(ValueError,match='FROZEN_DEPENDENCY_MISSING'):
        frozen_worker_parameters(tmp_path,{})
    with pytest.raises(ValueError,match='FROZEN_DEPENDENCY_MISSING'):
        frozen_worker_parameters(tmp_path,{'upstream_parameters':{'foldx_confirmation':{'binary':'missing.exe'}}})
