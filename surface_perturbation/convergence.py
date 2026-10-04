"""Verification convergence must account for the unresolved candidate pool."""


def convergence_state(current, previous, uncalibrated_candidate_ids):
    if uncalibrated_candidate_ids:
        return 'CALIBRATION_EXPANSION_REQUIRED'
    if current and current == previous:
        return 'VERIFIED_PANEL'
    return 'VERIFYING'
