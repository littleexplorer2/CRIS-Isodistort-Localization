"""One process-local lock per project for raw and structured computations."""
from pathlib import Path
import threading

_guard = threading.Lock()
_locks = {}


def project_lock(root):
    key = str(Path(root).resolve())
    with _guard:
        return _locks.setdefault(key, threading.Lock())
