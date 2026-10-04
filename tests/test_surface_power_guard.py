import os
import pytest
from ecd_epitope_scanner.surface_perturbation.power_guard import prevent_sleep


class FakeRequest:
    def __init__(self,fail=None):self.calls=[];self.fail=fail
    def set(self,kind):
        self.calls.append(('set',kind))
        if self.fail==kind:raise OSError('registration failed')
    def clear(self,kind):self.calls.append(('clear',kind))
    def close(self):self.calls.append(('close',))


def test_power_requests_released_after_worker_failure():
    api=FakeRequest()
    with pytest.raises(ValueError):
        with prevent_sleep(backend=api) as report:
            assert not report['display_required']
            raise ValueError('scientific worker failed')
    assert api.calls==[('set',1),('set',3),('clear',3),('clear',1),('close',)]


def test_partial_registration_failure_releases_first_request_and_handle():
    api=FakeRequest(fail=3)
    with pytest.raises(OSError):
        with prevent_sleep(backend=api):raise AssertionError('worker must not run')
    assert api.calls==[('set',1),('set',3),('clear',1),('close',)]


@pytest.mark.skipif(os.name!='nt',reason='Windows native API')
def test_native_power_request_registers_and_releases():
    with prevent_sleep() as report:
        assert report['status']=='POWER_REQUESTS_REGISTERED'
