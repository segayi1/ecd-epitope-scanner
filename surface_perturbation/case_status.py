"""Combine physical evidence with the cumulative runtime acceptance contract."""
import json
from pathlib import Path


def case_status(folder,summary=None):
    folder=Path(folder)
    if summary is None:
        path=folder/'surface_perturbation_preview.json'
        summary=json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    phase=summary.get('phase','production')
    acceptance=folder/'case_acceptance.json'
    if acceptance.exists():
        verdict=json.loads(acceptance.read_text(encoding='utf-8'))
        if verdict['status']!='WITHIN_CUMULATIVE_BUDGET':return verdict['status']
    ledger=folder/f'surface_{phase}_cumulative_budget.json'
    if ledger.exists():
        budget=json.loads(ledger.read_text(encoding='utf-8'))
        if budget.get('active'):return 'RUNNING_NOT_ACCEPTED'
        if budget['elapsed_seconds']>=budget['limit_seconds']:return 'INCOMPLETE_RUNTIME_EXHAUSTED'
    status = summary.get('status','INCOMPLETE_NO_CURRENT_SUMMARY')
    if summary.get('goal_status')=='ACHIEVED':
        from .combo_contract import audit_target_panel
        try:audit_target_panel(folder,summary)
        except (ValueError,OSError,KeyError):return 'INCOMPLETE_TARGET_PANEL_EVIDENCE'
    launch=folder/'launch_source_manifest.json'
    summary_path=folder/'surface_perturbation_preview.json'
    if launch.is_file() and (not summary_path.is_file() or summary_path.stat().st_mtime<launch.stat().st_mtime):
        return 'INCOMPLETE_NO_CURRENT_SUMMARY'
    driver=folder/'phase_result.json'
    if driver.is_file() and (not launch.is_file() or driver.stat().st_mtime>=launch.stat().st_mtime):
        driver_status=json.loads(driver.read_text(encoding='utf-8')).get('status','')
        if driver_status.startswith('INCOMPLETE'):return driver_status
    if status in {'VERIFIED_PANEL','TARGET_PANEL_VERIFIED_SEARCH_OPEN'} and (folder/'surface_perturbation/paired_manifest.json').exists():
        from ..scripts.audit_surface_mutation_policy import audit
        policy = audit(folder)
        if policy['status'] != 'CURRENT_POLICY_PASS':
            return policy['status']
    if status == 'VERIFIED_PANEL':
        requests=folder/'surface_perturbation/blindspot_calibration_requests.json'
        if requests.is_file():
            pending=set(json.loads(requests.read_text(encoding='utf-8')).get('candidate_ids',[]))
            closed=set(summary.get('hard_blindspot_upper_bound_excluded_candidate_ids',[]))
            if pending-closed:
                return 'CALIBRATION_EXPANSION_REQUIRED'
    return status
