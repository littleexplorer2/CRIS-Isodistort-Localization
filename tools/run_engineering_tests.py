"""Run the project's engineering tests and retain exact logs with fingerprints."""

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = root / "evidence" / "validation" / (timestamp + "-" + sys.platform)
    target.mkdir(parents=True)
    command = [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"]
    started = datetime.now(timezone.utc).isoformat()
    test_temp = root / "runs" / "engineering-test-temp"
    test_temp.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({key: str(test_temp) for key in ("TEMP", "TMP", "TMPDIR")})
    result = subprocess.run(command, cwd=root, capture_output=True, check=False, env=env)
    (target / "stdout.log").write_bytes(result.stdout)
    (target / "stderr.log").write_bytes(result.stderr)
    source_files = sorted((root / "cris_rea").rglob("*.py"))
    source_files += sorted((root / "tests").rglob("*.py"))
    source_files += sorted((root / "tests").rglob("*.js"))
    record = {
        "scope": "New-project engineering tests; no scientific acceptance claim.",
        "reason": "Retain actual command output and evidence for final reviewed implementation.",
        "started_at_utc": started,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "cwd": str(root),
        "exit_code": result.returncode,
        "host": platform.platform(),
        "python": sys.version,
        "temporary_directory": str(test_temp),
        "stdout_sha256": hashlib.sha256(result.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(result.stderr).hexdigest(),
        "source_sha256": {
            path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in source_files
        },
    }
    (target / "record.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"record": str(target / "record.json"), "exit_code": result.returncode}))
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
