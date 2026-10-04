import json
import os
from concurrent.futures import ThreadPoolExecutor
import pytest
from ecd_epitope_scanner import util


def test_concurrent_json_writers_publish_whole_payloads(tmp_path):
    path=tmp_path/'manifest.json'
    def publish(i):util.write_json(path,{'id':i,'payload':[i]*100})
    with ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(publish,range(40)))
    saved=json.loads(path.read_text(encoding='utf-8'))
    assert saved['payload']==[saved['id']]*100
    assert not list(tmp_path.glob('*.tmp'))


def test_transient_replace_failure_retains_previous_json_until_success(tmp_path,monkeypatch):
    path=tmp_path/'manifest.json';util.write_json(path,{'old':True})
    replace=util.os.replace;attempts=[]
    def occupied(source,destination):
        attempts.append(source)
        assert json.loads(path.read_text())=={'old':True}
        if len(attempts)<3:raise PermissionError('reader holds destination')
        replace(source,destination)
    monkeypatch.setattr(util.os,'replace',occupied)
    monkeypatch.setattr(util.time,'sleep',lambda duration:None)
    util.write_json(path,{'new':True})
    assert len(attempts)==3 and json.loads(path.read_text())=={'new':True}


def test_persistent_replace_failure_keeps_complete_old_and_attempted_payload(tmp_path,monkeypatch):
    path=tmp_path/'manifest.json';util.write_json(path,{'old':True})
    def denied(*args):raise PermissionError('persistent denial')
    monkeypatch.setattr(util.os,'replace',denied)
    monkeypatch.setattr(util.time,'sleep',lambda duration:None)
    with pytest.raises(PermissionError):util.write_json(path,{'new':True})
    assert json.loads(path.read_text())=={'old':True}
    assert json.loads(next(tmp_path.glob('*.tmp')).read_text())=={'new':True}


@pytest.mark.skipif(os.name!='nt',reason='Windows real non-delete-sharing handle')
def test_real_windows_reader_lock_is_retried_after_handle_closes(tmp_path):
    import ctypes,threading,time
    from ctypes import wintypes as w
    path=tmp_path/'paired_manifest.json';util.write_json(path,{'old':True})
    api=ctypes.WinDLL('kernel32',use_last_error=True)
    api.CreateFileW.argtypes=[w.LPCWSTR,w.DWORD,w.DWORD,w.LPVOID,w.DWORD,w.DWORD,w.HANDLE]
    api.CreateFileW.restype=w.HANDLE
    api.CloseHandle.argtypes=[w.HANDLE];api.CloseHandle.restype=w.BOOL
    handle=api.CreateFileW(str(path),0x80000000,3,None,3,0,None)
    if handle in (None,ctypes.c_void_p(-1).value):raise ctypes.WinError(ctypes.get_last_error())
    released=threading.Event()
    def release():
        time.sleep(.15);api.CloseHandle(handle);released.set()
    thread=threading.Thread(target=release);thread.start()
    try:util.write_json(path,{'new':True})
    finally:
        thread.join()
        if not released.is_set():api.CloseHandle(handle)
    assert json.loads(path.read_text())=={'new':True}
