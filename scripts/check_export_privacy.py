"""Check archive text for user/temporary directories and caller-supplied private markers."""

import argparse
import gzip
import io
import os
import re
import zipfile
from pathlib import Path

USER_PATH = re.compile(
    r'(?i)(?<![A-Za-z0-9])(?:[a-z]:[\\/]+(?:users|userdata)[\\/]+[^\\/\s"<>]+|'
    r'/home/(?![<[(])[^/\s"<>]+(?:/|(?=[\s"<>]|$))|'
    r'/(?:tmp|var/tmp)/(?![<[(])[^/\s"<>]+/|'
    r'/mnt/[a-z]/(?:users|userdata)/[^/\s"<>]+(?:/|(?=[\s"<>]|$)))'
)


def findings(name, data, markers):
    if any(marker.casefold() in name.casefold() for marker in markers):
        yield name
    vendor = "/vendor/" in name.replace("\\", "/")
    if not vendor and USER_PATH.search(name):
        yield name
    if name.lower().endswith(('.zip', '.whl')):
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            for member in archive.infolist():
                if not member.is_dir():
                    yield from findings(name + '!' + member.filename, archive.read(member), markers)
    elif name.lower().endswith('.gz'):
        yield from findings(name[:-3], gzip.decompress(data), markers)
    else:
        try:
            text = data.decode('utf-8-sig')
        except UnicodeDecodeError:
            return
        text = text.replace('\\\\', '\\').replace('\\/', '/')
        marker_hit = any(marker.casefold() in text.casefold() for marker in markers)
        if marker_hit or (not vendor and USER_PATH.search(text)):
            yield name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    markers = [item for item in os.environ.get('ECD_PRIVACY_MARKERS', '').split(os.pathsep) if item]
    found = list(findings(args.archive.name, args.archive.read_bytes(), markers))
    if found:
        print('Affected archive entries:', len(found))
    print('Privacy scan:', 'FAILED' if found else 'PASSED')
    return 1 if found else 0


if __name__ == '__main__':
    raise SystemExit(main())
