import json

from ecd_epitope_scanner.surface_perturbation.tools import PhysicalTools
from ecd_epitope_scanner.surface_perturbation.research_electrostatics import validate_preserved_pqr


def test_heavy_atom_preserving_protocol_is_separately_cached_and_uses_noopt_flags(tmp_path, monkeypatch):
    source = tmp_path / 'input.pdb'
    source.write_text('ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n')
    tools = PhysicalTools(tmp_path / 'cache', {'use_wsl': False, 'pdb2pqr': 'pdb2pqr-v2'})
    monkeypatch.setattr(tools, 'identity', lambda _path: 'a' * 64)
    calls = []

    def fake_run(argv, _cwd, timeout=300):
        calls.append(argv)
        output = argv[-1]
        with open(output, 'w', encoding='ascii') as handle:
            handle.write('ATOM      1   CA ALA A   1       1.000   2.000   3.000  0.0000  1.7000\nEND\n')
        return 'PDB2PQR test'

    monkeypatch.setattr(tools, 'run', fake_run)
    _path, pqr = tools.prepare_pqr(source, protocol='heavy_atom_preserving_v1')
    assert pqr[('A', 1, 'CA')]['xyz'] == [1., 2., 3.]
    assert '--nodebump' in calls[0] and '--noopt' in calls[0]
    contract = json.loads(next((tmp_path / 'cache').glob('pqr_*/complete.json')).read_text())
    assert contract['protocol'] == 'heavy_atom_preserving_v1'
    assert contract['optimize_hydrogens'] is False and contract['debump'] is False


def test_preserved_pqr_requires_exact_heavy_atom_coordinates_and_complete_parameters():
    source = {('A', 10, '', 'CA'): ([1., 2., 3.], 'ALA', 'C'),
              ('A', 10, '', 'CB'): ([2., 2., 3.], 'ALA', 'C')}
    residue_ids = {('A', 10, ''): 1}
    pqr = {('A', 1, 'CA'): {'xyz': [1., 2., 3.], 'residue': 'ALA', 'charge': 0., 'radius': 1.7},
           ('A', 1, 'CB'): {'xyz': [2., 2., 3.], 'residue': 'ALA', 'charge': 0., 'radius': 1.7},
           ('A', 1, 'H'): {'xyz': [1., 3., 3.], 'residue': 'ALA', 'charge': 0., 'radius': 1.0}}
    assert validate_preserved_pqr(source, residue_ids, pqr)['status'] == 'MATCHED'
    pqr[('A', 1, 'CB')]['xyz'] = [2.01, 2., 3.]
    assert validate_preserved_pqr(source, residue_ids, pqr)['status'] == 'UNKNOWN_PREPARATION_MISMATCH'
