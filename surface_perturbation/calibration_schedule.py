"""Append new matched null evidence without repeating exhausted cohorts.

This affects task priority only. A deferred unresolved field remains unknown
and still prevents global verification unless its real support bound closes.
"""
import numpy as np
from .null_store import FIELD_CHANNELS
from .calibration_scope import local_scope


def null_append_needed(calibration, required, area):
    required=np.asarray(required,bool)
    if not required.any() or any(name not in calibration for name in FIELD_CHANNELS):
        return True
    if local_scope(calibration,required,area)['complete']:
        return False
    for name in FIELD_CHANNELS:
        value=calibration[name]
        scope=np.asarray(value.get('required_mask',[]),bool)
        if scope.shape!=required.shape or np.any(required&~scope):
            return True
        if value.get('block_count',0)<19:
            return True
    # The exact matched 20-model cohort already reached the configured cap.
    # It is NOT promoted to calibrated/verified by this scheduling decision.
    return False


def calibration_buckets(candidates, evidence, supports, area, failed):
    fresh=[];retry=[];deferred=[]
    for candidate in candidates:
        cid=candidate.candidate_id
        if null_append_needed(evidence[cid],supports[cid],area):
            (retry if cid in failed else fresh).append(candidate)
        else:
            deferred.append(cid)
    return fresh,retry,deferred
