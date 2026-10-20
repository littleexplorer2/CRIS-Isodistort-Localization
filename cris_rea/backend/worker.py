"""Linux-only process supervisor, invoked as a script by the local backend.

This is an engineering supervisor, not a crystallographic implementation.
The parent writes its specification; scientific input is never executable code.
"""

import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time


def supervise(spec_path):
    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    result = {
        "status": "runtime_error",
        "exit_code": None,
        "platform": platform.platform(),
        "python": sys.version,
        "command": spec["command"],
        "cwd": spec["cwd"],
        "isodata": spec["isodata"],
    }
    started = time.monotonic()
    child = None
    try:
        if sys.platform != "linux":
            raise RuntimeError("The official executable requires Linux.")
        env = os.environ.copy()
        # A trailing slash is required by legacy Fortran path concatenation.
        env["ISODATA"] = spec["isodata"].rstrip("/") + "/"
        env["TMPDIR"] = spec["cwd"]
        env["LC_ALL"] = "C"
        result["environment_overrides"] = {
            key: env[key] for key in ("ISODATA", "TMPDIR", "LC_ALL")
        }
        with open(spec["input_path"], "rb") as source, \
                open(spec["stdout_path"], "wb") as out, \
                open(spec["stderr_path"], "wb") as err:
            child = subprocess.Popen(
                spec["command"], cwd=spec["cwd"], env=env,
                stdin=source if spec["input_mode"] == "stdin" else subprocess.DEVNULL,
                stdout=out, stderr=err, start_new_session=True, shell=False,
            )
            result["pid"] = child.pid
            while True:
                code = child.poll()
                output_bytes = sum(
                    Path(path).stat().st_size
                    for path in (spec["stdout_path"], spec["stderr_path"])
                )
                if output_bytes > spec["max_output_bytes"]:
                    result["status"] = "output_limit"
                    break
                if code is not None:
                    result["status"] = (
                        "completed_uninterpreted" if code == 0 else "process_failed"
                    )
                    break
                if time.monotonic() - started >= spec["timeout_seconds"]:
                    result["status"] = "timed_out"
                    break
                time.sleep(0.05)
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGKILL)
            result["exit_code"] = child.wait(timeout=10)
            # Terminate surviving subprocesses even when the group leader exited.
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait(timeout=10)
            except (OSError, subprocess.TimeoutExpired):
                pass
    finally:
        result["duration_seconds"] = round(time.monotonic() - started, 6)
        Path(spec["result_path"]).write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(supervise(sys.argv[1]))
