"""Scoped Windows System/Execution power requests for active scientific work.

The display may turn off. Explicit user sleep and DC Modern Standby policy
can override these requests; they do not change the cumulative time contract.
"""
import ctypes
import os
from contextlib import contextmanager


class WindowsPowerRequest:
    def __init__(self):
        from ctypes import wintypes as w
        class Detailed(ctypes.Structure):
            _fields_=[('module',w.HMODULE),('reason_id',w.ULONG),('count',w.ULONG),
                      ('strings',ctypes.POINTER(w.LPWSTR))]
        class Reason(ctypes.Union):
            _fields_=[('detailed',Detailed),('simple',w.LPWSTR)]
        class Context(ctypes.Structure):
            _fields_=[('version',w.ULONG),('flags',w.DWORD),('reason',Reason)]
        self.api=ctypes.WinDLL('kernel32',use_last_error=True)
        self.api.PowerCreateRequest.argtypes=[ctypes.POINTER(Context)]
        self.api.PowerCreateRequest.restype=w.HANDLE
        for name in ('PowerSetRequest','PowerClearRequest'):
            function=getattr(self.api,name);function.argtypes=[w.HANDLE,w.DWORD];function.restype=w.BOOL
        self.api.CloseHandle.argtypes=[w.HANDLE];self.api.CloseHandle.restype=w.BOOL
        self.context=Context();self.context.version=0;self.context.flags=1
        self.reason='ECD surface perturbation: active model, calibration and evidence work'
        self.context.reason.simple=self.reason
        self.handle=self.api.PowerCreateRequest(ctypes.byref(self.context))
        if self.handle in (None,ctypes.c_void_p(-1).value):raise ctypes.WinError(ctypes.get_last_error())

    def set(self,kind):
        if not self.api.PowerSetRequest(self.handle,kind):raise ctypes.WinError(ctypes.get_last_error())

    def clear(self,kind):
        if not self.api.PowerClearRequest(self.handle,kind):raise ctypes.WinError(ctypes.get_last_error())

    def close(self):
        if not self.api.CloseHandle(self.handle):raise ctypes.WinError(ctypes.get_last_error())


@contextmanager
def prevent_sleep(*,backend=None):
    if backend is None and os.name!='nt':
        yield {'status':'NOT_APPLICABLE','platform':os.name}
        return
    request=backend if backend is not None else WindowsPowerRequest()
    active=[]
    try:
        for kind in (1,3): # PowerRequestSystemRequired, PowerRequestExecutionRequired.
            request.set(kind);active.append(kind)
        yield {'status':'POWER_REQUESTS_REGISTERED','requests':['SYSTEM','EXECUTION'],
               'display_required':False,'limits':'user-initiated sleep and DC Modern Standby policy can override'}
    finally:
        try:
            for kind in reversed(active):request.clear(kind)
        finally:request.close()
