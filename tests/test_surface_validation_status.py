import json
import os
from ecd_epitope_scanner.scripts.summarize_surface_validation import phase_report


def test_historical_verified_summary_cannot_supply_current_coverage(tmp_path):
    folder=tmp_path/'production';folder.mkdir()
    summary=folder/'surface_perturbation_preview.json'
    summary.write_text(json.dumps({'status':'VERIFIED_PANEL','C_hard_all':.99}))
    launch=folder/'launch_source_manifest.json';launch.write_text('{}')
    os.utime(summary,(10,10));os.utime(launch,(20,20))
    result=phase_report(folder)
    assert result['status']=='INCOMPLETE_NO_CURRENT_SUMMARY'
    assert result['historical_summary_retained']
    assert result['coverage'] is None and not result['accepted_panel']


def test_generic_driver_error_does_not_erase_exhausted_budget(tmp_path):
    folder=tmp_path/'calibration';folder.mkdir()
    (folder/'phase_result.json').write_text(json.dumps({'status':'INCOMPLETE','reason':'timeout'}))
    (folder/'surface_calibration_cumulative_budget.json').write_text(json.dumps(
        {'active':False,'elapsed_seconds':14401,'limit_seconds':14400}))
    result=phase_report(folder)
    assert result['status']=='INCOMPLETE_RUNTIME_EXHAUSTED'
    assert result['driver_reason']=='timeout' and result['coverage'] is None


def test_saved_verified_status_is_rejected_when_blindspot_requests_remain(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.case_status import case_status
    (tmp_path/'surface_perturbation').mkdir()
    (tmp_path/'surface_perturbation/blindspot_calibration_requests.json').write_text(
        json.dumps({'candidate_ids':['unresolved']}))
    summary={'status':'VERIFIED_PANEL'}
    assert case_status(tmp_path,summary)=='CALIBRATION_EXPANSION_REQUIRED'
    summary['hard_blindspot_upper_bound_excluded_candidate_ids']=['unresolved']
    assert case_status(tmp_path,summary)=='VERIFIED_PANEL'


def test_accepted_primary_coverage_is_read_from_actual_nested_curve(tmp_path):
    folder=tmp_path/'production';folder.mkdir()
    (folder/'surface_perturbation_preview.json').write_text(json.dumps(
        {'status':'VERIFIED_PANEL','actual_constructs':5,'coverage_curve':[
            {'actual_constructs':5,'C_hard_all':.8,'C_hard_designable':.9,
             'C_extended_all':None,'C_extended_designable':None},
            {'actual_constructs':8,'C_hard_all':.95}]}))
    report=phase_report(folder)
    assert report['accepted_panel'] and report['coverage']['C_hard_all']==.8
    assert report['coverage']['C_extended_all'] is None


def test_benchmark_records_cumulative_timeout_as_specific_failure(tmp_path,monkeypatch):
    import sys
    from ecd_epitope_scanner.scripts import benchmark_surface_perturbation as driver
    def exhausted(*args):raise TimeoutError('CUMULATIVE_RUNTIME_EXHAUSTED')
    monkeypatch.setattr(driver,'run_case',exhausted)
    monkeypatch.setattr(sys,'argv',['benchmark','--case','PD1','--phase','production','--output',str(tmp_path)])
    driver.main()
    report=json.loads((tmp_path/'PD1/production/phase_result.json').read_text())
    assert report['status']=='INCOMPLETE_RUNTIME_EXHAUSTED'


def test_current_outer_failure_blocks_physical_panel_acceptance(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.case_status import case_status
    summary=tmp_path/'surface_perturbation_preview.json';summary.write_text('{"status":"VERIFIED_PANEL"}')
    launch=tmp_path/'launch_source_manifest.json';launch.write_text('{}')
    driver=tmp_path/'phase_result.json';driver.write_text('{"status":"INCOMPLETE","reason":"export failed"}')
    os.utime(launch,(20,20));os.utime(summary,(30,30));os.utime(driver,(40,40))
    assert case_status(tmp_path)=='INCOMPLETE'
    os.utime(driver,(10,10)) # A prior failure cannot reject a successful new attempt.
    assert case_status(tmp_path)=='VERIFIED_PANEL'
    os.utime(summary,(10,10))
    assert case_status(tmp_path)=='INCOMPLETE_NO_CURRENT_SUMMARY'


def test_coordinator_does_not_audit_stale_zip_as_a_current_delivery(tmp_path):
    from ecd_epitope_scanner.scripts.continue_surface_validation import current_bundle
    manifest=tmp_path/'run_manifest.json';manifest.write_text('{}')
    archive=tmp_path/'ecd_epitope_design.zip';archive.write_bytes(b'old bundle fixture')
    os.utime(manifest,(30,30));os.utime(archive,(10,10))
    assert not current_bundle(tmp_path,20)
    os.utime(archive,(30,30));assert current_bundle(tmp_path,20)


def test_production_continuation_requires_inactive_current_calibration_delivery(tmp_path):
    import pytest
    from ecd_epitope_scanner.scripts.continue_surface_validation import resumed_phases
    assert resumed_phases(tmp_path,'calibration')==['calibration','production']
    with pytest.raises(ValueError,match='COMPLETED_CURRENT'):
        resumed_phases(tmp_path,'production')
    folder=tmp_path/'calibration';folder.mkdir()
    ledger=folder/'surface_calibration_cumulative_budget.json'
    ledger.write_text(json.dumps({'active':True,'elapsed_seconds':12,'limit_seconds':14400}))
    for name in ('launch_source_manifest.json','surface_perturbation_preview.json','run_manifest.json','ecd_epitope_design.zip'):
        (folder/name).write_text('{}')
    with pytest.raises(ValueError,match='COMPLETED_CURRENT'):
        resumed_phases(tmp_path,'production')
    ledger.write_text(json.dumps({'active':False,'elapsed_seconds':12,'limit_seconds':14400}))
    # Partial scientific evidence never becomes complete via phase scheduling.
    (folder/'surface_perturbation_preview.json').write_text('{"status":"CALIBRATION_INSUFFICIENT"}')
    assert resumed_phases(tmp_path,'production')==['production']


def test_production_continuation_refuses_historical_summary_and_bundle(tmp_path):
    import pytest
    from ecd_epitope_scanner.scripts.continue_surface_validation import resumed_phases
    folder=tmp_path/'calibration';folder.mkdir()
    (folder/'surface_calibration_cumulative_budget.json').write_text('{"active":false}')
    for name in ('launch_source_manifest.json','surface_perturbation_preview.json','run_manifest.json','ecd_epitope_design.zip'):
        (folder/name).write_text('{}');os.utime(folder/name,(20,20))
    os.utime(folder/'surface_perturbation_preview.json',(10,10))
    with pytest.raises(ValueError,match='COMPLETED_CURRENT'):
        resumed_phases(tmp_path,'production')
    os.utime(folder/'surface_perturbation_preview.json',(30,30))
    os.utime(folder/'ecd_epitope_design.zip',(10,10))
    with pytest.raises(ValueError,match='COMPLETED_CURRENT'):
        resumed_phases(tmp_path,'production')
