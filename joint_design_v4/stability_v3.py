"""Preview risk policy, not experimentally calibrated expression thresholds."""
import math

CONTRACT = 'spurs_primary_foldx_sanity_v3'


def decision(spurs, foldx, order=2):
    finite = lambda x: x is not None and math.isfinite(float(x))
    if not finite(spurs):
        return 'NOT_SCORED'
    if order == 3:
        # SPURS higher-order score is diagnostic. Exact double controls are
        # enforced downstream; complete triple FoldX remains necessary.
        if not finite(foldx):
            return 'SPURS_ONLY_REVIEW_REQUIRED'
        return 'CONSENSUS_PASS' if foldx <= 2 else 'STRUCTURAL_ENERGY_CONFLICT'
    if spurs > 2:
        return 'COMBINATION_REJECT'
    if spurs > 1:
        return 'CONSENSUS_REVIEW'
    if not finite(foldx):
        return 'SPURS_ONLY_REVIEW_REQUIRED'
    return 'CONSENSUS_PASS' if foldx <= 2 else 'STRUCTURAL_ENERGY_CONFLICT'


def risk_band(spurs):
    if spurs is None or not math.isfinite(spurs):
        return 'UNKNOWN'
    return ('PREFERRED' if spurs <= .5 else 'CONFIRMABLE' if spurs <= 1
            else 'REVIEW' if spurs <= 2 else 'HIGH_RISK')
