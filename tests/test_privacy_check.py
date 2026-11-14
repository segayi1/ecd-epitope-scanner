import gzip
import io
import json
import os
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest


@pytest.mark.parametrize('root_parts', [('home',), ('tmp',), ('var', 'tmp')])
def test_scan_recurses_into_compressed_artifacts(tmp_path, root_parts):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    inner = io.BytesIO()
    with zipfile.ZipFile(inner, 'w') as archive:
        archive.writestr('report.json.gz', gzip.compress(
            json.dumps({'path': '/' + '/'.join((*root_parts, 'example_user', 'file'))}).encode()))
    outer = tmp_path / 'result.zip'
    with zipfile.ZipFile(outer, 'w') as archive:
        archive.writestr('nested.zip', inner.getvalue())
    result = subprocess.run([sys.executable, str(script), str(outer)], capture_output=True, text=True)
    assert result.returncode == 1
    assert 'Affected archive entries: 1' in result.stdout
    assert 'example_user' not in result.stdout


def test_scan_detects_escaped_windows_path_and_custom_marker(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    outer = tmp_path / 'result.zip'
    with zipfile.ZipFile(outer, 'w') as archive:
        archive.writestr('report.json', json.dumps({
            'path': 'C:' + chr(92) + 'Users' + chr(92) + 'sample' + chr(92) + 'data'}))
        archive.writestr('notes.txt', 'sample-private-marker')
    result = subprocess.run([sys.executable, str(script), str(outer)], capture_output=True,
                            text=True, env=dict(os.environ, ECD_PRIVACY_MARKERS='sample-private-marker'))
    assert result.returncode == 1
    assert 'Affected archive entries: 2' in result.stdout
    assert 'sample' not in result.stdout


def test_scan_handles_home_directory_without_trailing_slash(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    outer = tmp_path / 'result.zip'
    with zipfile.ZipFile(outer, 'w') as archive:
        archive.writestr('report.json', json.dumps({'home': '/' + 'home/' + 'example_user'}))
    result = subprocess.run([sys.executable, str(script), str(outer)], capture_output=True, text=True)
    assert result.returncode == 1
    assert 'example_user' not in result.stdout


def test_regex_source_is_not_reported_as_a_home_directory(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    outer = tmp_path / 'source.zip'
    with zipfile.ZipFile(outer, 'w') as archive:
        archive.writestr(
            'util.py',
            r'''USER_PATH = re.compile(r"/home/[^/\\s\\\"<>]+")'''
            '\n' + r'''OTHER = re.compile(r"/home/(?!<)[^/\\s\\\"<>]+")''',
        )
    result = subprocess.run([sys.executable, str(script), str(outer)], capture_output=True, text=True)
    assert result.returncode == 0
    assert 'Privacy scan: PASSED' in result.stdout


def test_vendor_code_path_examples_are_not_machine_data(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    outer = tmp_path / 'result.zip'
    with zipfile.ZipFile(outer, 'w') as archive:
        archive.writestr('assets/vendor/bundle.js', 'const example = ' + json.dumps(
            'C:' + chr(92) + 'Users' + chr(92) + 'sample' + chr(92) + 'file') + ';')
    result = subprocess.run([sys.executable, str(script), str(outer)], capture_output=True, text=True)
    assert result.returncode == 0


def test_custom_marker_is_checked_even_in_vendor_code(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    outer = tmp_path / 'result.zip'
    with zipfile.ZipFile(outer, 'w') as archive:
        archive.writestr('assets/vendor/bundle.js', 'private-marker')
    result = subprocess.run([sys.executable, str(script), str(outer)], capture_output=True, text=True,
                            env=dict(os.environ, ECD_PRIVACY_MARKERS='private-marker'))
    assert result.returncode == 1


def test_wheel_members_are_scanned(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    wheel = tmp_path / 'package.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('package/config.json', json.dumps({'path': '/' + 'home/' + 'example_user'}))
    result = subprocess.run([sys.executable, str(script), str(wheel)], capture_output=True, text=True)
    assert result.returncode == 1
    assert 'Affected archive entries: 1' in result.stdout
    assert 'example_user' not in result.stdout


def test_public_urls_and_app_routes_are_not_user_directories(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    archive_path = tmp_path / 'public.zip'
    with zipfile.ZipFile(archive_path, 'w') as archive:
        archive.writestr('links.json', json.dumps({
            'reference': 'https://example.org/tmp/reference',
            'home': 'https://example.org/home/reference',
            'structure': '/api/jobs/example/files/structure.cif',
            'shared_lock': '/tmp/ecd_v4_gpu.lock',
        }))
    result = subprocess.run([sys.executable, str(script), str(archive_path)], capture_output=True, text=True)
    assert result.returncode == 0
    assert 'Privacy scan: PASSED' in result.stdout


def test_tar_members_are_scanned_with_vendor_exemption(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'scripts' / 'check_export_privacy.py'
    archive_path = tmp_path / 'source.tar.gz'
    vendor = ('C:' + chr(92) + 'Users' + chr(92) + 'example' + chr(92) + 'file').encode()
    with tarfile.open(archive_path, 'w:gz') as archive:
        entry = tarfile.TarInfo('package/web/vendor/bundle.js')
        entry.size = len(vendor)
        archive.addfile(entry, io.BytesIO(vendor))
    result = subprocess.run([sys.executable, str(script), str(archive_path)], capture_output=True, text=True)
    assert result.returncode == 0
    with tarfile.open(archive_path, 'w:gz') as archive:
        entry = tarfile.TarInfo('package/config.json')
        entry.size = len(vendor)
        archive.addfile(entry, io.BytesIO(vendor))
    result = subprocess.run([sys.executable, str(script), str(archive_path)], capture_output=True, text=True)
    assert result.returncode == 1
