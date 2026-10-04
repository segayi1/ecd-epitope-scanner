"""Final case status requires completed audits, including explicit unknowns."""


def apbs_panel_ready(*, combo, construct_count, hard_coverage, final_panel_audit=False):
    """Scheduling only: a final partial panel can be audited without qualifying it."""
    return bool(construct_count) and (final_panel_audit or not combo or
        (construct_count >= 4 and hard_coverage is not None and hard_coverage >= .8))


def require_apbs_audits(expected,completed,*,requested=True,reports=None):
    missing=sorted(set(expected)-set(completed))
    if not requested:raise ValueError('FINAL_APBS_AUDIT_NOT_REQUESTED: hard draft only')
    if missing:raise ValueError('FINAL_APBS_AUDIT_INCOMPLETE: '+','.join(missing))
    reports=reports or {}
    def all_resolved(flag):
        return bool(expected) and all(reports.get(cid,{}).get('evidence_states',{}).get(flag,False) for cid in expected)
    return {'status':'APBS_EXECUTION_COMPLETE','execution_status':'COMPLETE',
        'preparation_status':'VALID' if all_resolved('APBS_PREPARATION_VALID') else 'UNRESOLVED',
        'numerical_resolution_status':'RESOLVED' if all_resolved('APBS_NUMERICALLY_RESOLVED') else 'UNRESOLVED',
        'completed_candidate_ids':list(completed),
        'evidence_states':{'APBS_EXECUTION_COMPLETE':True,
            'APBS_PREPARATION_VALID':all_resolved('APBS_PREPARATION_VALID'),
            'APBS_NUMERICALLY_RESOLVED':all_resolved('APBS_NUMERICALLY_RESOLVED'),
            'EXTENDED_COVERAGE_RESOLVED':all_resolved('EXTENDED_COVERAGE_RESOLVED')},
        'interpretation':'audit completion may retain unknown electrostatics; not verified electrostatic coverage'}
