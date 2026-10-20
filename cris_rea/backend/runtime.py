"""Official binary execution, input validation and evidence capture."""

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import uuid

MAX_INPUT_BYTES = 1024 * 1024
MAX_TIMEOUT_SECONDS = 300
MAX_OUTPUT_BYTES = 20 * 1024 * 1024
REQUIRED_DATA = (
    "const.dat", "data_diperiodic.txt", "data_images.txt", "data_irreps.txt",
    "data_isotropy.txt", "data_little.txt", "data_magnetic.txt", "data_space.txt",
    "data_ssg.txt", "data_ssgmag.txt", "data_wyckoff.txt",
)
PROGRAMS = {
    "iso": {
        "input_mode": "stdin",
        "description": "ISOTROPY command language; raw queries, no structure workflow",
    },
    "findsym": {
        "input_mode": "file",
        "description": "Identify crystallographic or magnetic space-group symmetry",
    },
    "findsym_cifinput": {
        "input_mode": "file",
        "description": "Convert a CIF to FINDSYM input; does not run FINDSYM",
    },
    "smodes": {
        "input_mode": "stdin",
        "description": "Symmetry-adapted displacement modes using SMODES input",
    },
    "comsubs": {
        "input_mode": "stdin",
        "description": "Common-subgroup search using COMSUBS input and bounds",
    },
}


class InputError(ValueError):
    """Invalid local execution request."""


def _now():
    return datetime.now(timezone.utc).isoformat()


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _fingerprint(path):
    return {"bytes": path.stat().st_size, "sha256": _sha(path)}


def _write_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def validate_request(program, data, timeout):
    if not isinstance(program, str) or program not in PROGRAMS:
        raise InputError("Choose an allowed program: " + ", ".join(PROGRAMS))
    if not isinstance(data, bytes):
        raise InputError("Input must be bytes.")
    if not data or len(data) > MAX_INPUT_BYTES:
        raise InputError(f"Input must contain 1 to {MAX_INPUT_BYTES} bytes.")
    if b"\x00" in data:
        raise InputError("Binary input containing NUL bytes is not supported.")
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) \
            or not math.isfinite(timeout) or not 0 < timeout <= MAX_TIMEOUT_SECONDS:
        raise InputError(f"Timeout must be greater than 0 and at most {MAX_TIMEOUT_SECONDS}s.")


def _wsl_path(path):
    completed = subprocess.run(
        ["wsl.exe", "--exec", "wslpath", "-a", str(path.resolve())],
        capture_output=True, timeout=15, check=False,
    )
    if completed.returncode:
        detail = completed.stderr.decode("utf-8", errors="replace")
        raise RuntimeError("WSL path conversion failed: " + detail.strip())
    value = completed.stdout.decode("utf-8", errors="strict").strip()
    if not value.startswith("/"):
        raise RuntimeError("WSL did not return an absolute Linux path.")
    return value


def inspect_runtime(project_root, probe=True):
    root = Path(project_root).resolve()
    source = root / "vendor" / "isobyu"
    info = {
        "host_platform": platform.platform(),
        "host_python": sys.version,
        "project_root": str(root),
        "vendor": str(source),
        "programs": {
            name: {"present": (source / name).is_file(), **details}
            for name, details in PROGRAMS.items()
        },
        "data_files": len(list(source.glob("data_*.txt"))),
        "execution_supported": False,
        "scientific_validation": "not_performed",
    }
    missing = [
        name for name in (*PROGRAMS, *REQUIRED_DATA)
        if not (source / name).is_file() or (source / name).is_symlink()
    ]
    info["missing_vendor_files"] = missing
    info["vendor_directory_safe"] = source.resolve() == source.absolute() and not source.is_symlink()
    info["vendor_ready"] = not missing and info["vendor_directory_safe"]
    if sys.platform == "linux":
        command = [sys.executable, "-c", "import platform; print(platform.platform())"]
    elif sys.platform == "win32":
        command = ["wsl.exe", "--exec", "python3", "-c",
                   "import platform; print(platform.platform())"]
    else:
        info["error"] = "Only Linux or Windows with WSL is supported."
        return info
    info["runtime_probe_command"] = command
    if not probe:
        info["runtime_probe"] = "not_requested"
        return info
    try:
        completed = subprocess.run(command, capture_output=True, timeout=15, check=False)
        info["runtime_probe"] = {
            "exit_code": completed.returncode,
            "stdout": completed.stdout.decode("utf-8", errors="replace"),
            "stderr": completed.stderr.decode("utf-8", errors="replace"),
        }
        info["execution_supported"] = completed.returncode == 0 and info["vendor_ready"]
    except (OSError, subprocess.TimeoutExpired) as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
    return info


def _safe_source_files(source, program):
    if source.is_symlink() or not source.is_dir() or source.resolve() != source.absolute():
        raise RuntimeError("vendor/isobyu must be a real directory inside this project.")
    selected = [source / program, source / "const.dat"]
    selected.extend(sorted(source.glob("data_*.txt")))
    if len(selected) < 3:
        raise RuntimeError("No official data_*.txt database files found.")
    for path in selected:
        if path.is_symlink() or not path.is_file() or path.resolve().parent != source.resolve():
            raise RuntimeError("Missing or unsafe vendor file: " + str(path))
    if (source / program).read_bytes()[:4] != b"\x7fELF":
        raise RuntimeError("Selected program is not an ELF binary.")
    return selected


def _project_fingerprints(root):
    return {
        str(path.relative_to(root)).replace("\\", "/"): _fingerprint(path)
        for path in sorted((root / "cris_rea").rglob("*.py"))
        if path.is_file() and not path.is_symlink()
    }


def _collect_artifacts(run_dir):
    result = {}
    for path in sorted(run_dir.rglob("*")):
        relative = path.relative_to(run_dir).as_posix()
        if relative.startswith("isodata/") or relative == "evidence.json":
            continue
        if path.is_symlink():
            result[relative] = {"kind": "symlink", "hashed": False}
        elif path.is_file():
            if path.resolve().is_relative_to(run_dir.resolve()):
                result[relative] = _fingerprint(path)
    return result


def _output_diagnostics(run_dir):
    """Recognize explicit vendor/runtime errors, without interpreting science."""
    patterns = (
        ("syntax_error", re.compile(r"^\s*(?:\*\s*)*Syntax error:.*$", re.IGNORECASE)),
        ("fortran_runtime_error", re.compile(
            r"^\s*(?:Fortran runtime error:|forrtl: severe).*", re.IGNORECASE
        )),
        ("explicit_error", re.compile(r"^\s*(?:\*\s*)*Error(?:\s|:).*$", re.IGNORECASE)),
        ("program_bombed", re.compile(r"^\s*This program has bombed:.*$", re.IGNORECASE)),
    )
    diagnostics = []
    for name in ("stdout.txt", "stderr.txt"):
        path = run_dir / name
        if not path.is_file() or path.is_symlink():
            continue
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line in stream:
                for code, pattern in patterns:
                    if pattern.match(line):
                        diagnostics.append({
                            "code": code, "source": name, "excerpt": line.strip()[:500],
                        })
                        break
                if len(diagnostics) >= 20:
                    return diagnostics
    return diagnostics


def execute_program(project_root, program, data, timeout=30):
    """Run an allowlisted official binary, preserving byte input and raw output.

    completed_uninterpreted records process exit only. It does not assert that
    the vendor accepted input or produced a correct scientific result.
    """
    validate_request(program, data, timeout)
    root = Path(project_root).resolve()
    source = root / "vendor" / "isobyu"
    runs = root / "runs"
    if runs.is_symlink() or runs.resolve() != runs:
        raise InputError("runs must be a real directory inside this project.")
    runs.mkdir(exist_ok=True)
    run_id = str(uuid.uuid4())
    run_dir = runs / run_id
    run_dir.mkdir()
    work = run_dir / "work"
    isodata = run_dir / "isodata"
    work.mkdir()
    isodata.mkdir()
    input_path = run_dir / "input.dat"
    input_path.write_bytes(data)
    evidence = {
        "schema_version": 1,
        "run_id": run_id,
        "program": program,
        "created_at_utc": _now(),
        "status": "preparing",
        "scientific_validation": "not_performed",
        "input": {"path": "input.dat", **_fingerprint(input_path)},
        "input_transformation": "none; original bytes",
        "timeout_seconds": timeout,
        "max_output_bytes": MAX_OUTPUT_BYTES,
        "host": {"platform": platform.platform(), "python": sys.version},
        "source_code": _project_fingerprints(root),
        "resource_role": "official generic executable and databases; no website answers",
        "vendor_resources_before": {},
        "isolation": {
            "cwd": str(work), "isodata": str(isodata),
            "os_filesystem_sandbox": False,
        },
    }
    _write_json(run_dir / "evidence.json", evidence)
    started = time.monotonic()
    try:
        for path in _safe_source_files(source, program):
            fingerprint = _fingerprint(path)
            destination = isodata / path.name
            shutil.copyfile(path, destination)
            if _fingerprint(destination) != fingerprint:
                raise RuntimeError("Vendor file changed during copy: " + path.name)
            evidence["vendor_resources_before"][path.name] = {
                "source": str(path), **fingerprint,
            }
        executable = isodata / program
        if sys.platform == "linux":
            executable.chmod(executable.stat().st_mode | 0o700)
            convert = lambda path: str(path.resolve())
            launcher = [sys.executable]
        elif sys.platform == "win32":
            # One conversion per run; subpaths are known generated ASCII names.
            linux_run = _wsl_path(run_dir)
            convert = lambda path: linux_run + "/" + path.relative_to(run_dir).as_posix()
            launcher = ["wsl.exe", "--exec", "python3"]
            chmod = subprocess.run(
                ["wsl.exe", "--exec", "chmod", "u+x", convert(executable)],
                capture_output=True, timeout=15, check=False,
            )
            if chmod.returncode:
                raise RuntimeError("WSL could not mark the copied executable runnable.")
        else:
            raise RuntimeError("Only Linux or Windows with WSL is supported.")
        worker_path = run_dir / "worker.py"
        shutil.copyfile(Path(__file__).with_name("worker.py"), worker_path)
        # Old Fortran helpers truncate long arguments and ISODATA values.
        # Resolve paths relative to the isolated working directory instead.
        command = ["../isodata/" + program]
        if PROGRAMS[program]["input_mode"] == "file":
            command.append("../input.dat")
        spec = {
            "command": command,
            "cwd": convert(work),
            "isodata": "../isodata/",
            "input_mode": PROGRAMS[program]["input_mode"],
            "input_path": convert(input_path),
            "stdout_path": convert(run_dir / "stdout.txt"),
            "stderr_path": convert(run_dir / "stderr.txt"),
            "result_path": convert(run_dir / "process.json"),
            "timeout_seconds": timeout,
            "max_output_bytes": MAX_OUTPUT_BYTES,
        }
        _write_json(run_dir / "execution.json", spec)
        supervisor = launcher + [convert(worker_path), convert(run_dir / "execution.json")]
        evidence["supervisor_command"] = supervisor
        evidence["binary_command"] = command
        evidence["status"] = "running"
        _write_json(run_dir / "evidence.json", evidence)
        completed = subprocess.run(
            supervisor, capture_output=True, timeout=timeout + 30, check=False,
        )
        (run_dir / "supervisor.stdout.txt").write_bytes(completed.stdout)
        (run_dir / "supervisor.stderr.txt").write_bytes(completed.stderr)
        evidence["supervisor_exit_code"] = completed.returncode
        process_path = run_dir / "process.json"
        if process_path.is_file():
            process = json.loads(process_path.read_text(encoding="utf-8"))
            evidence["process"] = process
            evidence["status"] = process["status"]
            if completed.returncode != 0:
                evidence["status"] = "supervisor_failed"
        else:
            evidence["status"] = "supervisor_failed"
            evidence["error"] = "Supervisor returned without a process result."
    except subprocess.TimeoutExpired as exc:
        evidence["status"] = "supervisor_timeout"
        evidence["error"] = (
            "Supervisor did not return. Child termination could not be confirmed: " + str(exc)
        )
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        evidence["status"] = "runtime_error"
        evidence["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        evidence["finished_at_utc"] = _now()
        evidence["duration_seconds"] = round(time.monotonic() - started, 6)
        evidence["isodata_after"] = {
            path.name: _fingerprint(path) for path in sorted(isodata.iterdir())
            if path.is_file() and not path.is_symlink()
        }
        evidence["diagnostics"] = _output_diagnostics(run_dir)
        if evidence["status"] == "completed_uninterpreted" and evidence["diagnostics"]:
            evidence["status"] = "program_diagnostic"
        evidence["artifacts"] = _collect_artifacts(run_dir)
        _write_json(run_dir / "evidence.json", evidence)
    return evidence
