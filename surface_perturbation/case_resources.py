"""Task-owned storage admission and content identity, without deleting evidence."""
from __future__ import annotations

import os
import shutil
import stat
import time
import json
from contextvars import ContextVar
from pathlib import Path

from ..util import sha256_file, write_json

ACTIVE_STORAGE = ContextVar("coverage_case_storage", default=None)
GIB = 1024 ** 3


def owned_bytes(root):
    total, seen = 0, set()
    pending=[os.fspath(root)]
    while pending:
        directory=pending.pop()
        try:
            with os.scandir(directory) as entries:
                for entry in entries:
                    # Match Path.rglob: traverse owned directories, never a
                    # symlink directory into an external model asset root.
                    if entry.is_dir(follow_symlinks=False):
                        pending.append(entry.path);continue
                    # Windows DirEntry.stat reports zero file identifiers;
                    # os.stat retains the NTFS identity needed for hardlinks.
                    try:st=os.stat(entry.path) if os.name=='nt' else entry.stat()
                    except FileNotFoundError:continue # Owned atomic replacement.
                    if not stat.S_ISREG(st.st_mode):continue
                    identity=(st.st_dev,st.st_ino) if st.st_ino else str(Path(entry.path).resolve())
                    if identity not in seen:
                        total+=st.st_size;seen.add(identity)
        except FileNotFoundError:
            continue # An owned short work directory may have completed.
    return total


class StorageGuard:
    def __init__(self, root, maximum=None, solve_limit=None, *, inventory_interval=None):
        self.root = Path(root).resolve()
        amendment_path=self.root/'resource_contract_amendment.json'
        amendment=json.loads(amendment_path.read_text(encoding='utf-8')) if amendment_path.is_file() else {}
        amended=amendment.get('activation')=='ACTIVE' and amendment.get('storage_byte_budget') is not None
        if maximum is None:maximum=amendment['storage_byte_budget'] if amended else 20*GIB
        if solve_limit is None:solve_limit=max(0,maximum-2*GIB) if amended else 18*GIB
        self.maximum = int(maximum)
        self.solve_limit = min(int(solve_limit), self.maximum)
        self.minimum_free=int(amendment.get('minimum_free_disk_bytes',0)) if amended else 0
        if inventory_interval is None:inventory_interval=30 if amended else 2
        if not 0 < inventory_interval <= 30:
            raise ValueError('Storage inventory interval must be positive and at most 30 seconds')
        self.inventory_interval=float(inventory_interval)
        self.peak = 0
        self._checked_at = -float('inf')
        self._last_used = 0

    def check(self, reserved=0, *, solving=True, force_inventory=False):
        if not force_inventory and reserved == 0 and time.monotonic()-self._checked_at < self.inventory_interval:
            used = self._last_used
        else:
            used = owned_bytes(self.root)
            self._checked_at, self._last_used = time.monotonic(), used
        self.peak = max(self.peak, used)
        limit = self.solve_limit if solving else self.maximum
        if used + reserved > limit:
            raise TimeoutError("CASE_STORAGE_ADMISSION_EXHAUSTED")
        if shutil.disk_usage(self.root).free < max(reserved+self.minimum_free, self.maximum-used):
            raise TimeoutError("CASE_DISK_RESERVE_UNAVAILABLE")
        return used


class EvidenceStore:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.objects = self.root / "objects"
        self.objects.mkdir(parents=True, exist_ok=True)
        # Immutable ingestion adds a hardlink, not another copy of the bytes.
        # Reuse its zero-reservation inventory cache within this checkpoint;
        # mutable copies still demand a fresh scan before every admission.
        self.storage = StorageGuard(self.root)

    def ingest(self, path, *, immutable=False, prepaid=False):
        path = Path(path).resolve()
        if self.root not in path.parents:
            raise ValueError("Evidence store only owns this task's files")
        digest = sha256_file(path)
        target = self.objects / digest[:2] / digest
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if sha256_file(target) != digest:
                raise ValueError("Content-addressed evidence corrupted")
        else:
            if not prepaid:self.storage.check(reserved=0 if immutable else path.stat().st_size,solving=False)
            if immutable:
                os.link(path, target)
            else:
                shutil.copyfile(path, target)
        if immutable and path.stat().st_ino != target.stat().st_ino:
            temporary = path.with_name(path.name + ".dedup.pending")
            os.link(target, temporary)
            temporary.replace(path)
        return {"sha256": digest, "bytes": target.stat().st_size,
                "object": target.relative_to(self.root).as_posix()}

    def checkpoint(self, stage, *, source_roots=None):
        entries = {}
        roots=[self.root] if source_roots is None else [Path(p).resolve() for p in source_roots]
        if any(p!=self.root and self.root not in p.parents for p in roots):raise ValueError('Checkpoint root outside owned case')
        admitted=[]
        for path in (p for root in roots for p in root.rglob('*')):
            if path.is_file() and self.objects not in path.parents and path.suffix in {".npz", ".pdb", ".fxout", ".off", ".pqr"}:
                # Mutable output names (canonical_reference, prediction_context,
                # coverage masks) must never share writable inodes with objects.
                relative=path.relative_to(self.root)
                if relative.parts[0] not in {'results','results_v3','fx','foldx_repair_cache'} and not (relative.parts[:2] == ('routes','coverage_v3')):
                    continue
                immutable=relative.parts[0]=='fx' or ('fields' in relative.parts
                    and path.name.startswith(('field_','scoped_','null_WT_')) and path.with_suffix('.json').is_file())
                admitted.append((path,immutable))
        # Reserve the whole possible copy peak once. A fresh whole-case scan
        # per mutable file made checkpoint admission quadratic in file count.
        self.storage.check(reserved=sum(p.stat().st_size for p,immutable in admitted if not immutable),solving=False,force_inventory=True)
        for path,immutable in admitted:
            entries[path.relative_to(self.root).as_posix()] = self.ingest(path,immutable=immutable,prepaid=True)
        self.storage.check(solving=False,force_inventory=True)
        write_json(self.root / (stage + "_object_manifest.json"), {
            "schema_version": "case_objects_1", "files": entries,
            "identity_rule": "Byte sharing never replaces scientific batch/group/run identity"})
        return entries


def affordable_batch(remaining_seconds, remaining_bytes, measured_seconds, measured_bytes):
    if measured_seconds <= 0 or measured_bytes <= 0:
        raise ValueError("Positive measured costs are required")
    time_capacity=256 if remaining_seconds==float('inf') else int(max(0,remaining_seconds-1800)//measured_seconds)
    count = max(0, min(8, time_capacity,
                       int(max(0, remaining_bytes-2*GIB)//measured_bytes)))
    capacity = min(time_capacity,
                   int(max(0,remaining_bytes-2*GIB)//measured_bytes))
    return count, min(256, max(32, 4*capacity))
