"""
ISODISTORT_VALIDATE 依赖与环境准备脚本。

在 CRIS 仓库根目录使用同一份虚拟环境 ``CRIS/.venv``：
1) 检查 Python >= 3.10
2) 若 .venv 不存在则创建；已存在且可用则复用
3) 只 pip 安装 requirements.txt 里尚未安装的包（已下载的不重复下载）
4) 检查运行本工具不需要额外环境变量（无 ISODATA / WSL 要求）
5) 若缺少 compare/、compare/item、compare/true 则自动创建

启动入口只有 ``main.py``（交互菜单，或 ``main.py compare`` / ``main.py batch``）。
批量比较前须把 compare/true/ 中官网 CIF 改成与 compare/item/ 相同的相对路径和文件名。

用法（在 CRIS 根目录）：
  python ISODISTORT_VALIDATE/main_requirement.py
  python ISODISTORT_VALIDATE/main_requirement.py --dev
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


def _run(
    cmd: list[str],
    *,
    check: bool = True,
    cwd: str | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
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


def _pip_install_missing(python: Path, req_file: Path) -> None:
    if not req_file.is_file():
        print(f"[pip] Skip (missing): {req_file}")
        return
    pkgs = _requirement_lines(req_file)
    if not pkgs:
        print(f"[pip] Skip (no packages listed): {req_file}")
        return
    installed = _installed_distributions(python)
    missing = [line for line in pkgs if _distribution_name(line) not in installed]
    if not missing:
        print(f"[pip] Already present, skip download: {req_file.name}")
        return
    print(f"[pip] Installing missing: {', '.join(missing)}")
    _run([str(python), "-m", "pip", "install", *missing])


def _check_python() -> None:
    if sys.version_info < (3, 10):
        raise RuntimeError(
            f"Python >= 3.10 required, got {sys.version_info.major}.{sys.version_info.minor}"
        )
    print(f"[python] {sys.version.split()[0]} ({sys.executable})")


def _check_env_vars() -> None:
    """VALIDATE 不依赖 ISODATA / WSL；入口脚本会自行加入 sys.path。"""
    print("\n[env] Checking environment variables ...")
    extra = []
    for name in ("ISODATA", "PYTHONPATH"):
        value = os.environ.get(name)
        if value:
            extra.append(f"{name}={value}")
    if extra:
        print("[env] Present (not required for VALIDATE):")
        for line in extra:
            print(f"      {line}")
    else:
        print("[env] No extra variables required. PYTHONPATH is set by the entry scripts.")
    print("[env] OK")


def _ensure_compare_dirs(validate_root: Path) -> None:
    compare = validate_root / "compare"
    item = compare / "item"
    true = compare / "true"
    for folder in (compare, item, true):
        if not folder.exists():
            print(f"[paths] Creating missing folder: {folder}")
        folder.mkdir(parents=True, exist_ok=True)
    print(f"[paths] OK: compare/item = {item}")
    print(f"[paths] OK: compare/true = {true}")


def _smoke_import(validate_root: Path, python: Path) -> None:
    print("\n[check] Import smoke-test ...")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(validate_root) + os.pathsep + env.get("PYTHONPATH", "")
    cp = _run(
        [
            str(python),
            "-c",
            (
                "from isodistort_validate.compare_paths import ensure_compare_dirs\n"
                "from isodistort_validate.compare_cif import compare_cif\n"
                "from isodistort_validate.batch_compare import run_batch\n"
                "ensure_compare_dirs()\n"
                "print('ISODISTORT_VALIDATE import OK')\n"
            ),
        ],
        cwd=str(validate_root),
        env=env,
    )
    print(f"[check] {(cp.stdout or '').strip().splitlines()[-1] if cp.stdout else 'OK'}")
    _run(
        [
            str(python),
            "-c",
            "import numpy, pymatgen, spglib; print('numpy/pymatgen/spglib OK')",
        ],
        env=env,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recreate", action="store_true", help="Force recreate CRIS/.venv")
    parser.add_argument("--dev", action="store_true", help="Also install requirements-dev.txt")
    args = parser.parse_args()

    validate_root = Path(__file__).resolve().parent
    project_root = validate_root.parent
    if not (project_root / "ISODISTORT").is_dir():
        raise RuntimeError(
            "This script must be run from the CRIS repo "
            "(ISODISTORT_VALIDATE/ next to ISODISTORT/)."
        )

    print(f"[root] {project_root}")
    _check_python()
    _check_env_vars()
    _ensure_compare_dirs(validate_root)

    _venv_dir, python = _ensure_venv(project_root, recreate=args.recreate)
    _pip_install_missing(python, validate_root / "requirements.txt")
    if args.dev:
        _pip_install_missing(python, validate_root / "requirements-dev.txt")

    _smoke_import(validate_root, python)

    print("\n=== DONE ===")
    if _is_windows():
        print(f"Use venv python: {python}")
        print(f"  {python} ISODISTORT_VALIDATE\\main.py")
        print(f"  {python} ISODISTORT_VALIDATE\\main.py compare")
        print(f"  {python} ISODISTORT_VALIDATE\\main.py batch")
        print("  批量比较前请把 compare/true/ 中官网 CIF 改成与 compare/item/ 相同的相对路径和文件名。")
    else:
        print(f"Use venv python: {python}")
        print(f"  {python} ISODISTORT_VALIDATE/main.py")
        print(f"  {python} ISODISTORT_VALIDATE/main.py compare")
        print(f"  {python} ISODISTORT_VALIDATE/main.py batch")
        print("  批量比较前请把 compare/true/ 中官网 CIF 改成与 compare/item/ 相同的相对路径和文件名。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
