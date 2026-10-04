"""Explicit raw-value reuse after an acceptance-only calibration revision."""
import ast
import hashlib
import json
from pathlib import Path


def verified_legacy_id(settings):
    location=settings.get('legacy_numerical_source_snapshot')
    if not location:return None
    snapshot=Path(location);root=Path(__file__).parent
    manifest=json.loads((snapshot/'manifest.json').read_text(encoding='utf-8'))
    names=('fields.py','geometry.py','chemistry.py','reference.py','calibration.py','tools.py','null_store.py')
    for name in names:
        if hashlib.sha256((snapshot/name).read_bytes()).hexdigest()!=manifest['source_sha256'][name]:
            raise ValueError('legacy numerical source snapshot changed')
    identity=hashlib.sha256(b''.join((snapshot/n).read_bytes() for n in names)).hexdigest()
    if identity!=manifest['implementation_id']:raise ValueError('legacy implementation identity mismatch')
    # Physical fields do not depend on calibration acceptance metadata.
    for name in ('fields.py','geometry.py','chemistry.py','reference.py','tools.py'):
        if (root/name).read_bytes()!=(snapshot/name).read_bytes():return None
    def numerical_ast(path):
        tree=ast.parse(path.read_text(encoding='utf-8'))
        return ast.dump(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='numerical_bounds'))
    if numerical_ast(root/'null_store.py')!=numerical_ast(snapshot/'null_store.py'):return None
    return identity
