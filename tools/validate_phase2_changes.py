"""Targeted checks for phase-2 changes, not a complete project regression."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MODULES = [
    "test_modes", "test_mode_boundaries", "test_native_bomb",
    "test_local_contract.ExecutionEvidenceTests.test_zero_exit_syntax_error_is_diagnostic",
    "test_local_contract.ExecutionEvidenceTests.test_consecutive_iso_prompts_do_not_hide_syntax_error",
    "test_local_contract.CliTests", "test_local_contract.HttpTests",
    "test_local_contract.BrowserInputContractTests",
]


def main():
    dest = ROOT / "evidence" / "phase2" / "targeted" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    dest.mkdir(parents=True)
    tmp = ROOT / "runs" / "phase2-test-temp"; tmp.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update({name: str(tmp) for name in ("TEMP", "TMP", "TMPDIR")})
    env["PYTHONPATH"] = os.pathsep.join((str(ROOT), str(ROOT / "tests")))
    command = [sys.executable, "-m", "unittest", "-v", *MODULES]
    result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, check=False)
    (dest / "stdout.txt").write_bytes(result.stdout)
    (dest / "stderr.txt").write_bytes(result.stderr)
    record = {"scope": "Changed mode core, layer/lock/artifact boundaries, native fatal diagnostic, "
                       "and affected pre-existing CLI/HTTP/browser contracts. Not full regression.",
              "command": command, "exit_code": result.returncode,
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (ROOT / "cris_rea").rglob("*.py")}}
    (dest / "record.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print(result.stderr.decode(errors="replace"))
    print(dest / "record.json")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
