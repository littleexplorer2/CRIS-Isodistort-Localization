"""
ISODISTORT / ISODISTORT_VALIDATE 统一依赖安装脚本

目标
----
1) 在 CRIS 仓库根目录只创建一份虚拟环境（默认：`CRIS/.venv`）
2) 安装并校验运行网页与终端交互所需的所有 Python 依赖
   （同时覆盖 ISODISTORT 与 ISODISTORT_VALIDATE 两部分）
3) 安装完成后给出最直接的运行方式提示

用法
----
从 CRIS 根目录执行：
  python ISODISTORT/main_requirement.py

可选：
  python ISODISTORT/main_requirement.py --recreate   # 强制重建 .venv
  python ISODISTORT/main_requirement.py --dev         # 额外安装 requirements-dev.txt（仅开发用途）
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _run(cmd: list[str], *, check: bool = True, cwd: str | None = None, env: dict[str, str] | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        check=check,
        cwd=cwd,
        env=env,
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


def _pip_install(python: Path, req_file: Path) -> None:
    print(f"\n[pip] Installing: {req_file}")
    _run([str(python), "-m", "pip", "install", "-r", str(req_file)])


def _pip_upgrade(python: Path) -> None:
    print("\n[pip] Upgrading pip/setuptools/wheel ...")
    _run([str(python), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])


def _ensure_venv(project_root: Path, *, recreate: bool) -> tuple[Path, Path]:
    venv_dir = project_root / ".venv"
    if recreate and venv_dir.exists():
        print(f"[venv] Recreating: {venv_dir}")
        shutil.rmtree(venv_dir, ignore_errors=True)

    if not venv_dir.exists():
        print(f"[venv] Creating: {venv_dir}")
        _run([sys.executable, "-m", "venv", str(venv_dir)])

    python = _venv_python(venv_dir)
    if not python.exists():
        raise RuntimeError(f"Virtualenv python not found: {python}")
    return venv_dir, python


def _check_wsl() -> None:
    # WSL 不是 Python 依赖，但本项目底层 iso/findsym/smodes 都在 WSL 运行。
    if not _is_windows():
        return
    if shutil.which("wsl.exe") is None:
        raise RuntimeError("Missing `wsl.exe` in PATH. Please install/initialize WSL first.")
    try:
        print("\n[wsl] Checking WSL status ...")
        cp = subprocess.run(
            ["wsl.exe", "--status"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if cp.returncode != 0:
            print("[wsl] WARNING: `wsl --status` failed:")
            print(cp.stdout[:8000])
        else:
            print("[wsl] OK")

        # 最基本的 WSL shell 可执行性（确保能跑 `wsl.exe -e sh -c ...`）
        print("[wsl] Checking `sh` inside WSL ...")
        sh_cp = subprocess.run(
            ["wsl.exe", "-e", "sh", "-c", "echo WSL_SH_OK"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
        )
        if sh_cp.returncode == 0 and "WSL_SH_OK" in sh_cp.stdout:
            print("[wsl] sh OK")
        else:
            raise RuntimeError("WSL shell test failed: `wsl.exe -e sh -c ...` did not run successfully.")
    except Exception as exc:
        raise RuntimeError(f"WSL check failed: {exc}") from exc


def _check_python_runtime() -> None:
    # 运行代码使用了 Python 3.10+ 的类型语法（如 list[str]、| 联合类型等）
    if sys.version_info < (3, 10):
        raise RuntimeError(f"Python >= 3.10 required, got {sys.version_info.major}.{sys.version_info.minor}")
    if shutil.which(sys.executable) is None:
        # 这条通常不会触发，但作为兜底提示
        print(f"[python] WARNING: unable to locate executable for: {sys.executable}")


def _check_project_paths(project_root: Path) -> None:
    # 必要配置文件存在性
    cfg = project_root / "ISODISTORT" / "config" / "settings.yaml"
    if not cfg.exists():
        raise RuntimeError(f"Missing config file: {cfg}")

    # 运行目录（wrapper 会尝试创建，但这里先尽量给出更早的错误提示）
    out_dir = project_root / "ISODISTORT" / "output"
    tmp_dir = out_dir / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)


def _check_isotropy_binaries(project_root: Path) -> None:
    """
    isobyu 目录下的 Linux ELF 二进制和数据库文件是“功能性必需项”。
    这里只做存在性检查；真正可运行性还受 WSL/权限/数据库完整性影响。
    """
    isobyu = project_root / "ISODISTORT" / "isobyu"
    if not isobyu.exists():
        print("\n[binary] WARNING: `ISODISTORT/isobyu/` not found. You need to deploy ISOTROPY suite there.")
        return

    # settings.yaml 里固定的二进制名
    required_bins = ["iso", "smodes"]
    optional_bins = ["findsym", "comsubs"]

    missing: list[str] = []
    for name in required_bins:
        if not (isobyu / name).exists():
            missing.append(name)

    # 数据库：至少应该有 data_*.txt（具体文件名随套件版本不同）
    data_files = list(isobyu.glob("data_*.txt")) + list(isobyu.glob("data*.txt"))
    if not data_files:
        raise RuntimeError("Missing ISOTROPY database files under `ISODISTORT/isobyu/` (expected `data_*.txt`).")

    for name in optional_bins:
        if not (isobyu / name).exists():
            print(f"[binary] INFO: optional binary missing `{name}` (not always required).")

    if missing:
        details = ", ".join(missing)
        raise RuntimeError(f"Missing required ISOTROPY binaries under `ISODISTORT/isobyu/`: {details}")
    else:
        print("\n[binary] OK: required binaries found (iso + smodes).")



def _post_import_smoke_check(project_root: Path, python: Path) -> None:
    # 只做 import 层面的快速校验（不触发 iso 二进制运行）。
    print("\n[check] Import smoke-test ...")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(project_root / "ISODISTORT") + os.pathsep + env.get("PYTHONPATH", "")
    _run(
        [
            str(python),
            "-c",
            "from isocore.api import IsoDistort; print('IsoDistort import OK')",
        ],
        env=env,
    )
    _run(
        [
            str(python),
            "-c",
            "import numpy, pymatgen, spglib; print('numpy/pymatgen/spglib OK')",
        ],
        env=env,
    )
    _run(
        [
            str(python),
            "-c",
            "import compare_cif; print('ISODISTORT_VALIDATE import OK')",
        ],
        cwd=str(project_root / "ISODISTORT_VALIDATE"),
        env=env,
    )

    _run(
        [
            str(python),
            "-c",
            "from web.server import main as web_main; print('web.server import OK')",
        ],
        cwd=str(project_root / "ISODISTORT"),
        env=env,
    )

    _run(
        [
            str(python),
            "-c",
            "import main_terminal; print('main_terminal import OK')",
        ],
        cwd=str(project_root / "ISODISTORT"),
        env=env,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recreate", action="store_true", help="Force recreate CRIS/.venv")
    ap.add_argument("--dev", action="store_true", help="Also install requirements-dev.txt (dev tests)")
    args = ap.parse_args()

    project_root = Path(__file__).resolve().parent.parent  # CRIS/
    print(f"[root] {project_root}")
    if not (project_root / "ISODISTORT").exists() or not (project_root / "ISODISTORT_VALIDATE").exists():
        raise RuntimeError("This script must be run inside the CRIS repo (missing ISODISTORT or ISODISTORT_VALIDATE).")

    _check_python_runtime()
    _check_project_paths(project_root)
    _check_wsl()
    _check_isotropy_binaries(project_root)

    venv_dir, python = _ensure_venv(project_root, recreate=args.recreate)
    _pip_upgrade(python)

    _pip_install(python, project_root / "ISODISTORT" / "requirements.txt")
    _pip_install(python, project_root / "ISODISTORT_VALIDATE" / "requirements.txt")

    if args.dev:
        _pip_install(python, project_root / "ISODISTORT" / "requirements-dev.txt")

    _post_import_smoke_check(project_root, python)

    print("\n=== DONE ===")
    if _is_windows():
        print(f"Use venv python: {python}")
        print("Run web:")
        print(f"  {python} ISODISTORT\\main_web.py")
        print("Run terminal:")
        print(f"  {python} ISODISTORT\\main_terminal.py")
        print("Run ISODISTORT_VALIDATE:")
        print(f"  {python} ISODISTORT_VALIDATE\\main.py")
    else:
        print(f"Use venv python: {python}")
        print("Run web:")
        print(f"  {python} ISODISTORT/main_web.py")
        print("Run terminal:")
        print(f"  {python} ISODISTORT/main_terminal.py")
        print("Run ISODISTORT_VALIDATE:")
        print(f"  {python} ISODISTORT_VALIDATE/main.py")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

