import json
import pytest
from ecd_epitope_scanner.surface_perturbation import resume_pairs as resume


def pair(run,**changes):
    return {'batch_path':'saved','group':1,'run':run,'mutations':['A1V'],
            'input_sha256':'repair','normalized_input_sha256':'input','binary_sha256':'binary',**changes}


def test_model_resume_preserves_order_identity_and_rejects_different_inputs(tmp_path,monkeypatch):
    old=tmp_path/'paired_manifest.json';rows=[pair(3),pair(1),pair(2)]
    old.write_text(json.dumps({'candidate':rows}))
    validated=[];monkeypatch.setattr(resume,'validate_pair',lambda p:validated.append(p['run']))
    snapshots=resume.capture_manifests(tmp_path,tmp_path)
    old.write_text(json.dumps({'candidate':[pair(0)]}))
    result=resume.restore_pairs('candidate',('A1V',),[pair(0)],snapshots)
    assert [p['run'] for p in result]==[0,3,1,2]
    assert len({(p['batch_path'],p['group'],p['run']) for p in result})==4
    assert json.loads(snapshots[1].read_text())['candidate']==rows
    assert resume.restore_pairs('candidate',('A1V',),[pair(0,normalized_input_sha256='new')],snapshots)==[pair(0,normalized_input_sha256='new')]
    assert resume.restore_pairs('candidate',('A1L',),[pair(0,mutations=['A1L'])],snapshots)==[pair(0,mutations=['A1L'])]
    def corrupt(p):raise ValueError('cached model hash changed')
    monkeypatch.setattr(resume,'validate_pair',corrupt)
    with pytest.raises(ValueError,match='hash changed'):
        resume.restore_pairs('candidate',('A1V',),[pair(0)],snapshots)
