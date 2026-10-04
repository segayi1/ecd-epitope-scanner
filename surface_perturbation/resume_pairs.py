"""Restore exact validated model identities, independently of field/calibration reuse."""
import json
import shutil
import time
from pathlib import Path
from .paired import validate_pair


def capture_manifests(folder, calibration_folder):
    folder=Path(folder);current=folder/'paired_manifest.json'
    if current.is_file():
        snapshot=folder/'paired_manifest_history'/str(time.time_ns())/'paired_manifest.json'
        snapshot.parent.mkdir(parents=True,exist_ok=False)
        shutil.copy2(current,snapshot)
    result=[]
    for root in dict.fromkeys((folder,Path(calibration_folder))):
        result.append(root/'paired_manifest.json')
        history=root/'paired_manifest_history'
        result.extend(sorted(history.glob('*/paired_manifest.json'),key=lambda p:p.stat().st_mtime,reverse=True))
    return result


def restore_pairs(candidate_id, signature, initial, manifests):
    initial=list(initial)
    keys=('input_sha256','normalized_input_sha256','binary_sha256')
    if not initial or any(not initial[0].get(key) for key in keys):
        return initial
    expected={key:initial[0][key] for key in keys}
    result=[];identities=set()
    for path in dict.fromkeys(map(Path,manifests)):
        if not path.is_file():continue
        for pair in json.loads(path.read_text(encoding='utf-8')).get(candidate_id,[]):
            if tuple(sorted(pair.get('mutations',[])))!=signature:continue
            if any(pair.get(key)!=value for key,value in expected.items()):continue
            validate_pair(pair)
            identity=(pair['batch_path'],pair['group'],pair['run'])
            if identity not in identities:result.append(pair);identities.add(identity)
    for pair in initial:
        identity=(pair['batch_path'],pair['group'],pair['run'])
        if identity not in identities:result.append(pair);identities.add(identity)
    return result
