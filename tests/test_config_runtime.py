from ecd_epitope_scanner import config
from ecd_epitope_scanner.cli import main
from ecd_epitope_scanner.structure_analysis import analyze_residues


def test_sasa_threshold_defaults_are_centralized():
    defaults = analyze_residues.__defaults__
    assert defaults is not None
    assert config.RSASA_THRESHOLD in defaults
    assert config.SIDECHAIN_SASA_THRESHOLD_A2 in defaults


def test_local_server_rejects_network_binding():
    try:
        main(["serve", "--host", "0.0.0.0"])
    except SystemExit as exc:
        assert "loopback" in str(exc)
    else:
        raise AssertionError("A local profile must reject public host binding")
