from pathlib import Path
import pytest
from ecd_epitope_scanner.integrations.workers.foldx_worker import _execution_root


def test_long_native_execution_uses_persistent_short_root_without_moving_output(tmp_path):
    output=tmp_path/('original_'+('x'*100));output.mkdir()
    source=output/'prior_evidence.txt';source.write_text('retained')
    base=Path.cwd()/'foldx_path_test'
    # cwd can itself be long in CI: choose the filesystem temp root explicitly.
    import tempfile
    base=Path(tempfile.gettempdir())/'ecd_foldx_test'
    first=_execution_root(tmp_path,output,{'engine_work_root':str(base)},windows=True)
    second=_execution_root(tmp_path,output,{'engine_work_root':str(base)},windows=True)
    assert first!=second and first.is_dir() and second.is_dir()
    assert len(str(first))+90<=220
    assert source.read_text()=='retained' and output.is_dir()
    assert _execution_root(tmp_path,output,{},windows=False)==output
    with pytest.raises(RuntimeError,match='too long'):
        _execution_root(tmp_path,output,{'engine_work_root':str(output)},windows=True)
    # Empty unit-test directories only; physical run directories are retained.
    first.rmdir();second.rmdir()


def test_optin_binary_alias_is_byte_identical_persistent_and_rejects_corruption(tmp_path):
    from ecd_epitope_scanner.integrations.workers.foldx_worker import _execution_binary
    source=tmp_path/'foldx_long_version_name.exe';source.write_bytes(b'original scientific executable')
    parameters={'short_binary_alias':True,'engine_work_root':str(tmp_path)}
    alias=_execution_binary(source,tmp_path,parameters,windows=True)
    assert alias.name=='foldx.exe' and alias.read_bytes()==source.read_bytes()
    assert _execution_binary(source,tmp_path,parameters,windows=True)==alias
    assert _execution_binary(source,tmp_path,{},windows=True)==source
    alias.write_bytes(b'corrupt')
    with pytest.raises(RuntimeError,match='fingerprint'):
        _execution_binary(source,tmp_path,parameters,windows=True)
