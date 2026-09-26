"""GD dependency/environment bootstrap script.

What it does:
1) Checks host Python version (>= 3.10) for running this script
2) Ensures the dedicated GD/.venv exists, and reuses it when it is already usable
3) Requires the GD/.venv interpreter to be Python 3.12.x
4) Installs only missing packages from requirements.txt
5) Optionally installs requirements-dev.txt with --dev

Version matching policy:
  只按 major.minor（3.12）匹配。`py` 启动器（Python Install Manager）只按
  major.minor 登记安装，`py -3.12.5` 会直接失败并报
  "No suitable Python runtime found"；因此本脚本永远只请求 `py -3.12`，
  具体补丁号由解释器上报后再校验，不把补丁号当版本标签使用。
  需要指定某个精确解释器时用 `--python <path>`。

Usage:
  python main_requirement.py
  python main_requirement.py --dev
  python main_requirement.py --recreate
  python main_requirement.py --python C:\\path\\to\\python.exe
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


# 只用 major.minor 匹配，避免把补丁号（3.12.5）当作启动器版本标签。
TARGET_PYTHON_SERIES = (3, 12)
TARGET_PYTHON_TEXT = "3.12"

# 依赖文件名：优先标准名，兼容旧名（早期脚本写成 requirement*.txt）。
REQUIREMENTS_RUNTIME_NAMES = ("requirements.txt", "requirement.txt")
REQUIREMENTS_DEV_NAMES = ("requirements-dev.txt", "requirement-dev.txt")


def _run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        check=check,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def _is_windows() -> bool:
    return sys.platform.startswith("win")


def _venv_python(venv_dir: Path) -> Path:
    if _is_windows():
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _requirement_lines(req_file: Path) -> list[str]:
    lines: list[str] = []
    for raw in req_file.read_text(encoding="utf-8").splitlines():
        text = raw.split("#", 1)[0].strip()
        if text:
            lines.append(text)
    return lines


def _distribution_name(req_line: str) -> str:
    name = re.split(r"[=<>!~\[]", req_line, maxsplit=1)[0].strip()
    return name.lower().replace("_", "-")


def _installed_distributions(python: Path) -> set[str]:
    cp = _run([str(python), "-m", "pip", "list", "--format=json"], check=False)
    if cp.returncode != 0 or not (cp.stdout or "").strip():
        return set()
    try:
        rows = json.loads(cp.stdout)
    except json.JSONDecodeError:
        return set()
    return {
        str(row.get("name", "")).lower().replace("_", "-")
        for row in rows
        if row.get("name")
    }


def _python_version(python: Path) -> tuple[int, int, int]:
    cp = _run(
        [
            str(python),
            "-c",
            "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}')",
        ]
    )
    text = (cp.stdout or "").strip()
    parts = text.split(".")
    if len(parts) != 3:
        raise RuntimeError(f"Unable to parse Python version from interpreter output: {text!r}")
    return int(parts[0]), int(parts[1]), int(parts[2])


def _python_series(python: Path) -> tuple[int, int]:
    major, minor, _ = _python_version(python)
    return major, minor


def _target_python_text() -> str:
    return TARGET_PYTHON_TEXT


def _pick_requirement_file(root: Path, names: tuple[str, ...]) -> Path:
    """返回第一个存在的依赖文件；都不存在时返回首选名，交由调用方提示缺失。"""
    for name in names:
        candidate = root / name
        if candidate.is_file():
            return candidate
    return root / names[0]


_PROBE = (
    "import sys; print(sys.executable); "
    "print(sys.version_info.major, sys.version_info.minor, sys.version_info.micro)"
)


def _candidate_launchers() -> list[list[str]]:
    """按“解释器名”列举候选命令，绝不使用补丁版本标签。

    `py -3.12` 合法；`py -3.12.5`/`py -3.14.6` 这类补丁号标签在 py 启动器上
    不合法（会报 "No suitable Python runtime found"）。
    """
    if _is_windows():
        return [["py", "-3.12"], ["py", "-3"], ["python3.12"], ["python"]]
    return [["python3.12"], ["python3"], ["python"]]


def _probe_interpreter(cmd: list[str]) -> tuple[Path, tuple[int, int, int]] | None:
    cp = _run([*cmd, "-c", _PROBE], check=False)
    if cp.returncode != 0 or not (cp.stdout or "").strip():
        return None
    lines = [line.strip() for line in (cp.stdout or "").splitlines() if line.strip()]
    if len(lines) < 2:
        return None
    tokens = lines[1].split()
    if len(tokens) != 3:
        return None
    try:
        version = (int(tokens[0]), int(tokens[1]), int(tokens[2]))
    except ValueError:
        return None
    return Path(lines[0]), version


def _resolve_target_python() -> Path:
    problems: list[str] = []
    for cmd in _candidate_launchers():
        probed = _probe_interpreter(cmd)
        label = " ".join(cmd)
        if probed is None:
            problems.append(f"{label}: not available")
            continue
        exe, version = probed
        if (version[0], version[1]) != TARGET_PYTHON_SERIES:
            problems.append(f"{label}: found Python {version[0]}.{version[1]}.{version[2]}")
            continue
        print(f"[python-probe] {label} -> {exe} (Python {version[0]}.{version[1]}.{version[2]})")
        return exe

    detail = "\n  ".join(problems) if problems else "(no candidates tried)"
    raise RuntimeError(
        f"No Python {_target_python_text()}.x interpreter found. Tried:\n  {detail}\n"
        f"Install Python {_target_python_text()} (any patch release) or pass an existing "
        "interpreter with --python <path>."
    )


def _requires_tensorflow(req_file: Path) -> bool:
    if not req_file.is_file():
        return False
    for line in _requirement_lines(req_file):
        if _distribution_name(line) == "tensorflow":
            return True
    return False


def _is_windows_file_lock_error(output: str) -> bool:
    text = output.lower()
    return ("winerror 32" in text) or ("being used by another process" in text)


def _install_one_package(python: Path, pkg: str, *, retries: int = 2) -> None:
    attempt = 0
    while True:
        attempt += 1
        try:
            _run([str(python), "-m", "pip", "install", pkg])
            return
        except subprocess.CalledProcessError as exc:
            out = (exc.stdout or "").strip()
            if _is_windows_file_lock_error(out) and attempt <= retries + 1:
                print(
                    f"[pip] File lock detected while installing {pkg} (attempt {attempt}/{retries + 1}). Retrying..."
                )
                time.sleep(3)
                continue
            raise RuntimeError(
                f"Failed to install package: {pkg}\n"
                f"pip output:\n{out}"
            ) from exc


def _check_python() -> None:
    if sys.version_info < (3, 10):
        raise RuntimeError(
            f"Python >= 3.10 required to run this script, got {sys.version_info.major}.{sys.version_info.minor}"
        )
    print(f"[python-host] {sys.version.split()[0]}")


def _ensure_venv(
    gd_root: Path,
    *,
    recreate: bool,
    python_exe: Path | None = None,
) -> tuple[Path, Path]:
    venv_dir = gd_root / ".venv"
    python = _venv_python(venv_dir)

    def usable() -> bool:
        if not python.exists():
            return False
        try:
            return _python_series(python) == TARGET_PYTHON_SERIES
        except Exception:
            return False

    if recreate and venv_dir.exists():
        print(f"[venv] Recreating: {venv_dir}")
        shutil.rmtree(venv_dir, ignore_errors=True)

    if venv_dir.exists() and not usable():
        print(f"[venv] Removing unusable .venv (missing or wrong Python series): {venv_dir}")
        shutil.rmtree(venv_dir, ignore_errors=True)

    if not venv_dir.exists():
        # 只有在必须新建时才去解析基础解释器：已有可用的 .venv 时不再依赖 `py`。
        if python_exe is not None:
            target_python = Path(python_exe)
            if not target_python.exists():
                raise RuntimeError(f"--python path does not exist: {target_python}")
            series = _python_series(target_python)
            if series != TARGET_PYTHON_SERIES:
                raise RuntimeError(
                    f"--python points at Python {series[0]}.{series[1]}, "
                    f"but GD requires {_target_python_text()}.x."
                )
            print(f"[python-target] {target_python} (from --python)")
        else:
            target_python = _resolve_target_python()
            print(f"[python-target] {target_python} ({_target_python_text()}.x)")
        print(f"[venv] Creating: {venv_dir}")
        _run([str(target_python), "-m", "venv", str(venv_dir)])
    else:
        print(f"[venv] Reusing: {venv_dir}")

    python = _venv_python(venv_dir)
    if not python.exists():
        raise RuntimeError(f"Virtualenv python not found: {python}")

    detected = _python_version(python)
    if (detected[0], detected[1]) != TARGET_PYTHON_SERIES:
        raise RuntimeError(
            "GD/.venv Python series mismatch: "
            f"found {detected[0]}.{detected[1]}, required {_target_python_text()}.x. "
            "Run with --recreate (optionally --python <path>) after installing the target Python."
        )
    print(f"[venv] python = {python} (Python {detected[0]}.{detected[1]}.{detected[2]})")
    return venv_dir, python


def _pip_install_missing(python: Path, req_file: Path) -> None:
    if not req_file.is_file():
        print(f"[pip] Skip missing file: {req_file}")
        return

    requested = _requirement_lines(req_file)
    if not requested:
        print(f"[pip] No packages in {req_file.name}")
        return

    installed = _installed_distributions(python)
    missing = [pkg for pkg in requested if _distribution_name(pkg) not in installed]

    if not missing:
        print(f"[pip] Already installed, skip: {req_file.name}")
        return

    print(f"[pip] Installing missing from {req_file.name}: {', '.join(missing)}")
    py_ver = _python_version(python)
    for pkg in missing:
        try:
            _install_one_package(python, pkg)
        except RuntimeError as exc:
            name = _distribution_name(pkg)
            message = str(exc)
            if name == "tensorflow" and py_ver >= (3, 13, 0):
                raise RuntimeError(
                    "TensorFlow is not available for this Python version "
                    f"({py_ver[0]}.{py_ver[1]}.{py_ver[2]}). "
                    f"Please use Python {_target_python_text()}.x for GD/.venv, then rerun with --recreate."
                ) from exc
            if "winerror 32" in message.lower() or "being used by another process" in message.lower():
                raise RuntimeError(
                    "Windows file lock prevented package installation. "
                    "Please close all Python/Jupyter processes using GD/.venv, pause OneDrive sync for this folder, "
                    "then rerun: python main_requirement.py --recreate"
                ) from exc
            raise


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare GD environment and dependencies")
    parser.add_argument("--dev", action="store_true", help="Also install requirements-dev.txt")
    parser.add_argument("--recreate", action="store_true", help="Recreate GD/.venv")
    parser.add_argument(
        "--python",
        dest="python_exe",
        default=None,
        help=f"Base interpreter used to create GD/.venv (default: auto-detect Python {_target_python_text()}.x)",
    )
    args = parser.parse_args()

    gd_root = Path(__file__).resolve().parent
    req_runtime = _pick_requirement_file(gd_root, REQUIREMENTS_RUNTIME_NAMES)
    req_dev = _pick_requirement_file(gd_root, REQUIREMENTS_DEV_NAMES)

    _check_python()
    _, venv_python = _ensure_venv(
        gd_root,
        recreate=args.recreate,
        python_exe=Path(args.python_exe) if args.python_exe else None,
    )

    # TensorFlow wheels are not published for Python >= 3.13 yet.
    if _requires_tensorflow(req_runtime):
        py_ver = _python_version(venv_python)
        if py_ver >= (3, 13, 0):
            raise RuntimeError(
                "Current GD/.venv uses Python "
                f"{py_ver[0]}.{py_ver[1]}.{py_ver[2]}, but TensorFlow in {req_runtime.name} "
                f"needs Python <= {_target_python_text()}. "
                f"Please recreate GD/.venv using Python {_target_python_text()}.x."
            )

    # Upgrade core packaging tools once per environment setup.
    _run([str(venv_python), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])

    _pip_install_missing(venv_python, req_runtime)
    if args.dev:
        _pip_install_missing(venv_python, req_dev)

    print("\n=== DONE ===")
    print(f"Use interpreter: {venv_python}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"[error] {exc}")
        raise SystemExit(1)
