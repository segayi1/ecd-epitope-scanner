from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import sqlite3
import tempfile
import time
from pathlib import Path, PurePosixPath, PureWindowsPath
from types import TracebackType
from typing import Any

import requests

USER_AGENT = "ecd-epitope-workbench/1.0 (+local scientific workflow)"


class ClosingSQLiteConnection(sqlite3.Connection):
    """Keep sqlite connection context managers from leaking file handles."""

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> bool | None:
        try:
            return super().__exit__(exc_type, exc_value, traceback)
        finally:
            self.close()


class PipelineError(RuntimeError):
    pass


class WaitingForAF3(PipelineError):
    def __init__(self, message: str, candidates: list[dict[str, Any]] | None = None):
        super().__init__(message)
        self.candidates = candidates or []


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return cleaned.strip("._") or "item"


def portable_command(command: list[str]) -> list[str]:
    known_names = {'python', 'python.exe', 'wsl_bridge.py', 'spurs_gpu_worker.py',
                   'discotope_gpu_worker.py', 'procaliper_worker.py', 'foldx_worker.py'}
    def redact(argument: str) -> str:
        prefix, separator, value = str(argument).partition('=')
        value = value if separator else prefix
        windows, posix = PureWindowsPath(value), PurePosixPath(value)
        if not windows.is_absolute() and not posix.is_absolute():
            return str(argument)
        name = windows.name if windows.is_absolute() else posix.name
        safe_tool = name.lower() in known_names or name.lower().endswith('.exe')
        replacement = f'<local>/{name}' if safe_tool else '<local-path>'
        return prefix + '=' + replacement if separator else replacement
    return [redact(argument) for argument in command]


_LOCAL_PATH = re.compile(
    r"(?i)(?<![A-Za-z0-9])(?:[A-Z]:[\\/][^\r\n\"'<>|,;]+|/(?:home|tmp|var/tmp)/[^\r\n\"'<>|,;]+|/mnt/[A-Za-z]/[^\r\n\"'<>|,;]+)"
)
_LOCAL_TOOL_NAMES = {
    "python", "python.exe", "wsl_bridge.py", "spurs_gpu_worker.py",
    "discotope_gpu_worker.py", "procaliper_worker.py", "foldx_worker.py",
}
_TEXT_ARTIFACT_SUFFIXES = {
    ".json", ".jsonl", ".md", ".csv", ".txt", ".log", ".xml", ".yaml", ".yml",
}


def portable_text(text: str) -> str:
    """Replace machine-local absolute paths while preserving useful basenames."""
    def replace(match: re.Match[str]) -> str:
        path = match.group(0).rstrip(".,;)")
        name = re.split(r"[\\/]", path.rstrip("\\/"))[-1].strip()
        safe_tool = name.lower() in _LOCAL_TOOL_NAMES or name.lower().endswith(".exe")
        replacement = f"<local>/{name}" if safe_tool else "<local-path>"
        return replacement + match.group(0)[len(path):]
    return _LOCAL_PATH.sub(replace, text)


def portable_projection(value: Any) -> Any:
    """Redact local paths in-place for JSON-safe data."""
    if isinstance(value, dict):
        for key, item in value.items():
            value[key] = portable_projection(item)
        return value
    if isinstance(value, list):
        for index, item in enumerate(value):
            value[index] = portable_projection(item)
        return value
    if isinstance(value, tuple):
        return [portable_projection(item) for item in value]
    if isinstance(value, str):
        return portable_text(value) if _LOCAL_PATH.search(value) else value
    return value


def sanitize_text_file(path: Path) -> bool:
    """Redact local paths in a text artifact in place and report whether it changed."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".gz":
        descriptor, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".sanitize", dir=path.parent)
        os.close(descriptor)
        tmp = Path(name)
        changed = False
        try:
            with (gzip.open(path, "rt", encoding="utf-8-sig", newline="") as source,
                  gzip.open(tmp, "wt", encoding="utf-8", newline="") as target):
                for line in source:
                    sanitized = portable_text(line)
                    changed = changed or sanitized != line
                    target.write(sanitized)
        except (OSError, UnicodeDecodeError):
            tmp.unlink(missing_ok=True)
            return False
        if not changed:
            tmp.unlink(missing_ok=True)
            return False
        os.replace(tmp, path)
        return True
    if suffix not in _TEXT_ARTIFACT_SUFFIXES:
        return False
    try:
        text = path.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeDecodeError):
        return False
    if suffix == ".json" and path.stat().st_size <= 32 * 1024 * 1024:
        try:
            value = json.loads(text)
        except json.JSONDecodeError:
            return False
        sanitized = json.dumps(portable_projection(value), ensure_ascii=False, indent=2) + "\n"
    elif suffix == ".jsonl":
        rows = []
        for line in text.splitlines():
            if line.strip():
                try:
                    rows.append(portable_projection(json.loads(line)))
                except json.JSONDecodeError:
                    return False
        sanitized = "".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in rows)
    else:
        # shortcut: very large JSON uses the string scanner to keep peak memory bounded.
        sanitized = portable_text(text)
    if sanitized == text:
        return False
    path.write_text(sanitized, encoding="utf-8")
    return True


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Independent writers must never share a temporary file; close it before
    # replacing on Windows, where readers can briefly deny the replacement.
    descriptor, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    tmp = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    for attempt in range(12):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if attempt == 11:
                # Preserve both the previous complete JSON and the attempted
                # payload for diagnosis; never truncate the destination.
                raise
            time.sleep(min(.005 * 2**attempt, .1))


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


class RemoteClient:
    def __init__(self, timeout: int = 45, retries: int = 3):
        self.timeout = timeout
        self.retries = retries
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": USER_AGENT, "Accept": "*/*"})

    def request(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        last: Exception | None = None
        for attempt in range(self.retries):
            try:
                response = self.session.request(method, url, timeout=self.timeout, **kwargs)
                response.raise_for_status()
                return response
            except requests.HTTPError as exc:
                # 4xx answers are deterministic; retrying them only burns time
                # and masks the real error.  Only 5xx/network faults retry.
                status = exc.response.status_code if exc.response is not None else None
                if status is not None and status < 500:
                    raise PipelineError(f"Remote request rejected: {url}: {status}") from exc
                last = exc
                if attempt + 1 < self.retries:
                    time.sleep(0.8 * (attempt + 1))
            except requests.RequestException as exc:
                last = exc
                if attempt + 1 < self.retries:
                    time.sleep(0.8 * (attempt + 1))
        raise PipelineError(f"Remote request failed: {url}: {last}") from last

    def json(self, url: str) -> Any:
        return self.request("GET", url, headers={"Accept": "application/json"}).json()

    def text(self, url: str) -> str:
        return self.request("GET", url).text

    def download(self, url: str, path: Path) -> Path:
        response = self.request("GET", url, stream=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        # Atomic: a previously interrupted half-written file must never be
        # mistaken for a complete structure on a later run.
        tmp = path.with_suffix(path.suffix + ".part")
        try:
            with tmp.open("wb") as handle:
                for chunk in response.iter_content(1024 * 256):
                    if chunk:
                        handle.write(chunk)
            os.replace(tmp, path)
        finally:
            tmp.unlink(missing_ok=True)
        return path
