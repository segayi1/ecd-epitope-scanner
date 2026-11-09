"""Single-flight JSON worker client with bounded request and exit waits."""
import json
import os
import queue
import signal
import subprocess
import tempfile
import threading
import uuid
from pathlib import Path

from .base import ToolError


class PersistentRunner:
    def __init__(self, command, timeout_seconds, log_path=None):
        self.timeout_seconds = timeout_seconds
        self.wsl_command = None
        self.control_path = None
        if '--persistent' in command and '--distro' in command and '--python' in command:
            self.wsl_command = ['wsl.exe', '-d', command[command.index('--distro') + 1],
                                '--exec', command[command.index('--python') + 1]]
            descriptor, control = tempfile.mkstemp(suffix='.spurs-process.json')
            os.close(descriptor)
            self.control_path = Path(control)
            command = [*command, '--control-file', control]
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
        self.stderr = (log_path.open('a', encoding='utf-8') if log_path is not None
                       else tempfile.TemporaryFile(mode='w+t', encoding='utf-8'))
        try:
            self.process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=self.stderr, text=True, encoding='utf-8', bufsize=1,
                start_new_session=os.name != 'nt')
        except BaseException:
            self.stderr.close()
            raise
        self.lines = queue.Queue()
        self.lock = threading.Lock()
        self.reader = threading.Thread(target=self._read, daemon=True)
        self.reader.start()

    def _read(self):
        try:
            for line in self.process.stdout:
                self.lines.put(line)
        finally:
            self.lines.put(None)

    def run(self, request):
        with self.lock:
            request_id = uuid.uuid4().hex
            try:
                self.process.stdin.write(json.dumps({'request_id': request_id, 'request': request}) + '\n')
                self.process.stdin.flush()
                line = self.lines.get(timeout=self.timeout_seconds)
                if line is None:
                    raise ValueError('worker exited without response')
                response = json.loads(line)
                if response.get('request_id') != request_id or not isinstance(response.get('result'), dict):
                    raise ValueError('worker response identity mismatch')
                if response['result'].get('status') in {'failed', 'partial'}:
                    self.close()
                return response['result']
            except BaseException as exc:
                self.close()
                if not isinstance(exc, (OSError, ValueError, queue.Empty)):
                    raise
                raise ToolError('TOOL_TIMEOUT' if isinstance(exc, queue.Empty) else 'TOOL_CRASHED', str(exc)) from exc

    def close(self):
        if self.process.poll() is None:
            try:
                self.process.stdin.close()
            except OSError:
                pass
            try:
                self.process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                if os.name == 'nt':
                    subprocess.run(['taskkill', '/PID', str(self.process.pid), '/T', '/F'],
                                   capture_output=True, check=False)
                else:
                    os.killpg(self.process.pid, signal.SIGKILL)
                self.process.wait(timeout=10)
        if self.control_path is not None:
            try:
                control = json.loads(self.control_path.read_text(encoding='utf-8'))
            except (OSError, ValueError):
                control = None
            if control:
                script = ("import os,signal,sys; p=int(sys.argv[1]); s=sys.argv[2]; "
                          "f='/proc/%d/stat'%p; "
                          "alive=os.path.exists(f) and open(f).read().split()[21]==s; "
                          "os.killpg(p,signal.SIGKILL) if alive else None")
                completed = subprocess.run([*self.wsl_command, '-c', script,
                    str(control['pid']), control['start_ticks']], capture_output=True, timeout=15)
                if completed.returncode:
                    raise ToolError('TOOL_CRASHED', 'WSL worker cleanup failed')
            self.control_path.unlink(missing_ok=True)
        self.process.stdout.close()
        self.stderr.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
