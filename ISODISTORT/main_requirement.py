"""
CRIS 统一依赖安装脚本（ISODISTORT / ISODISTORT_VALIDATE / ISOVIZ_INPUT）

目标
----
1) 在 CRIS 仓库根目录只创建一份虚拟环境（默认：`CRIS/.venv`）
2) 安装并校验运行网页、终端、IsoVIZ 导入所需的 Python 依赖
   （已安装的包不会重复下载）
3) 若缺少 ``ISODISTORT/output`` 则自动新建；若缺少 ``ISODISTORT/isobyu``
   则新建空目录，并提醒从 https://iso.byu.edu/isotropy.php 下载套件
4) 检查并配置 ``ISODATA``（ISOTROPY 数据库目录）
5) 安装完成后给出最直接的运行方式提示

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
import json
import os
import re
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


def _venv_usable(python: Path) -> bool:
    """True when the venv interpreter runs and matches this process's major.minor."""
    if not python.exists():
        return False
    cp = _run(
        [str(python), "-c", "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"],
        check=False,
    )
    if cp.returncode != 0:
        return False
    reported = (cp.stdout or "").strip()
    expected = f"{sys.version_info.major}.{sys.version_info.minor}"
    return reported == expected


def _requirement_lines(req_file: Path) -> list[str]:
    lines = []
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
    return {str(row.get("name", "")).lower().replace("_", "-") for row in rows if row.get("name")}


def _pip_install(python: Path, req_file: Path) -> None:
    if not req_file.is_file():
        print(f"\n[pip] Skip (missing): {req_file}")
        return
    pkgs = _requirement_lines(req_file)
    if not pkgs:
        print(f"\n[pip] Skip (no packages listed): {req_file}")
        return
    installed = _installed_distributions(python)
    missing = [line for line in pkgs if _distribution_name(line) not in installed]
    if not missing:
        print(f"\n[pip] Already installed, skip: {req_file}")
        return
    print(f"\n[pip] Installing missing from {req_file}: {', '.join(missing)}")
    _run([str(python), "-m", "pip", "install", *missing])


def _pip_upgrade(python: Path) -> None:
    print("\n[pip] Upgrading pip/setuptools/wheel ...")
    _run([str(python), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"])


def _ensure_venv(project_root: Path, *, recreate: bool) -> tuple[Path, Path]:
    venv_dir = project_root / ".venv"
    python = _venv_python(venv_dir)
    stale = venv_dir.exists() and not _venv_usable(python)
    if (recreate or stale) and venv_dir.exists():
        reason = "requested" if recreate else "broken or built with a different Python"
        print(f"[venv] Recreating ({reason}): {venv_dir}")
        shutil.rmtree(venv_dir, ignore_errors=True)
        if venv_dir.exists():
            raise RuntimeError(f"Could not remove old venv: {venv_dir}")

    if not venv_dir.exists():
        print(f"[venv] Creating: {venv_dir}")
        _run([sys.executable, "-m", "venv", str(venv_dir)])
    else:
        print(f"[venv] Reusing: {venv_dir}")

    python = _venv_python(venv_dir)
    if not _venv_usable(python):
        raise RuntimeError(f"Virtualenv python not usable: {python}")
    print(f"[venv] python = {python}")
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
    print(f"[python] {sys.version.split()[0]} ({sys.executable})")
    if shutil.which(sys.executable) is None:
        # 这条通常不会触发，但作为兜底提示
        print(f"[python] WARNING: unable to locate executable for: {sys.executable}")


def _check_project_paths(project_root: Path) -> None:
    """Ensure ISODISTORT/output exists; create it when missing."""
    cfg = project_root / "ISODISTORT" / "config" / "settings.yaml"
    if not cfg.exists():
        raise RuntimeError(f"Missing config file: {cfg}")

    out_dir = project_root / "ISODISTORT" / "output"
    tmp_dir = out_dir / "tmp"
    if not out_dir.exists():
        print(f"\n[paths] Creating missing folder: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    print(f"[paths] OK: output folder = {out_dir}")


_ISOTROPY_DOWNLOAD = "https://iso.byu.edu/isotropy.php"


def _remind_isobyu(isobyu: Path, *, created: bool) -> None:
    print("\n[binary] ISOTROPY Suite (iso / smodes / data_*.txt) is required.")
    if created:
        print(f"[binary] Created empty folder: {isobyu}")
    print(f"[binary] Download the Linux ISOTROPY Suite from {_ISOTROPY_DOWNLOAD}")
    print(f"[binary] Then move `iso`, `smodes`, and all `data_*.txt` into:\n         {isobyu}")


def _check_isotropy_binaries(project_root: Path) -> bool:
    """
    Check ISODISTORT/isobyu. If the folder is missing, create it and tell the
    user to download the suite from iso.byu.edu.

    Returns True when iso + smodes + data files are present.
    """
    isobyu = project_root / "ISODISTORT" / "isobyu"
    if not isobyu.exists():
        isobyu.mkdir(parents=True, exist_ok=True)
        _remind_isobyu(isobyu, created=True)
        return False

    required_bins = ["iso", "smodes"]
    optional_bins = ["findsym", "comsubs"]
    missing = [name for name in required_bins if not (isobyu / name).exists()]
    data_files = list(isobyu.glob("data_*.txt")) + list(isobyu.glob("data*.txt"))

    for name in optional_bins:
        if not (isobyu / name).exists():
            print(f"[binary] INFO: optional binary missing `{name}` (not always required).")

    if missing or not data_files:
        _remind_isobyu(isobyu, created=False)
        if missing:
            print(f"[binary] Missing binaries: {', '.join(missing)}")
        if not data_files:
            print("[binary] Missing database files (expected data_*.txt).")
        return False

    print("\n[binary] OK: required binaries found (iso + smodes).")
    return True


def _check_and_setup_isodata(project_root: Path, python: Path) -> None:
    """
    配置并校验 ``ISODATA``（iso / findsym / smodes 读取的数据库目录）。

    网页、终端和 Python API 在加载 ``config/settings.yaml`` 时会把
    ``isobyu.data_dir`` 写入进程环境变量 ``ISODATA``。本安装脚本：
    1) 确认数据库目录存在
    2) 在当前安装进程中设置 ``ISODATA``
    3) 用刚装好的 venv 解释器 import isocore，确认运行时也会设置该变量

    WSL 侧二进制还要求路径以 ``/`` 结尾；封装层会在首次运行时建立短路径
    符号链接，不必在 Windows 系统属性里永久 setx ``ISODATA``。
    """
    print("\n[env] Checking ISODATA ...")
    data_dir = project_root / "ISODISTORT" / "isobyu"
    if not data_dir.is_dir() or not (
        list(data_dir.glob("data_*.txt")) + list(data_dir.glob("data*.txt"))
    ):
        print("[env] SKIP: isobyu/data_*.txt not present yet; set ISODATA after you install the suite.")
        return
    os.environ["ISODATA"] = str(data_dir)
    print(f"[env] ISODATA (installer process) = {os.environ['ISODATA']}")

    env = os.environ.copy()
    env["PYTHONPATH"] = str(project_root / "ISODISTORT") + os.pathsep + env.get("PYTHONPATH", "")
    cp = _run(
        [
            str(python),
            "-c",
            (
                "from isocore.utils import get_config\n"
                "import os\n"
                "get_config()\n"
                "value = os.environ.get('ISODATA', '')\n"
                "assert value, 'ISODATA was not set by get_config()'\n"
                "print(value)\n"
            ),
        ],
        env=env,
    )
    imported = (cp.stdout or "").strip().splitlines()[-1] if cp.stdout else ""
    print(f"[env] ISODATA (after isocore import) = {imported}")
    print("[env] OK. Runtime sets ISODATA from config/settings.yaml; WSL uses a short symlink ending with /.")


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

    isoviz_root = project_root / "ISOVIZ_INPUT"
    if (isoviz_root / "isoviz_input").is_dir():
        env["PYTHONPATH"] = str(isoviz_root) + os.pathsep + env.get("PYTHONPATH", "")
        _run(
            [
                str(python),
                "-c",
                "from isoviz_input.amplitudes import read_amplitude_csv; print('ISOVIZ_INPUT import OK')",
            ],
            cwd=str(isoviz_root),
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
    has_isobyu = _check_isotropy_binaries(project_root)

    venv_dir, python = _ensure_venv(project_root, recreate=args.recreate)
    _pip_upgrade(python)

    _pip_install(python, project_root / "ISODISTORT" / "requirements.txt")
    _pip_install(python, project_root / "ISODISTORT_VALIDATE" / "requirements.txt")
    isoviz_req = project_root / "ISOVIZ_INPUT" / "requirements.txt"
    if isoviz_req.is_file():
        _pip_install(python, isoviz_req)

    if args.dev:
        _pip_install(python, project_root / "ISODISTORT" / "requirements-dev.txt")

    _post_import_smoke_check(project_root, python)
    if has_isobyu:
        _check_and_setup_isodata(project_root, python)
    else:
        print("\n[env] SKIP ISODATA until ISODISTORT/isobyu contains the ISOTROPY binaries.")

    print("\n=== DONE ===")
    if _is_windows():
        print(f"Use venv python: {python}")
        print("Run web:")
        print(f"  {python} ISODISTORT\\main_web.py")
        print("Run terminal:")
        print(f"  {python} ISODISTORT\\main_terminal.py")
        print("Run ISODISTORT_VALIDATE:")
        print(f"  {python} ISODISTORT_VALIDATE\\main.py")
        print("Run ISOVIZ_INPUT:")
        print(f"  {python} ISOVIZ_INPUT\\main.py")
    else:
        print(f"Use venv python: {python}")
        print("Run web:")
        print(f"  {python} ISODISTORT/main_web.py")
        print("Run terminal:")
        print(f"  {python} ISODISTORT/main_terminal.py")
        print("Run ISODISTORT_VALIDATE:")
        print(f"  {python} ISODISTORT_VALIDATE/main.py")
        print("Run ISOVIZ_INPUT:")
        print(f"  {python} ISOVIZ_INPUT/main.py")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

