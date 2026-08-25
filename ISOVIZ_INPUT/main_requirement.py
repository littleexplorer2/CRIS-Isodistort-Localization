"""
ISOVIZ_INPUT 依赖与环境准备脚本。

在 CRIS 仓库根目录使用同一份虚拟环境 ``CRIS/.venv``：
1) 检查 Python >= 3.10，没有则提示
2) 若 .venv 不存在则创建
3) 只 pip 安装 requirements.txt 里尚未安装的包
4) 检查 Java（IsoVIZ 运行需要），以及 CRIS 根目录的 ISOViz.lnk
5) 若缺少 ISOVIZ_INPUT/output 则自动新建

用法（在 CRIS 根目录）：
  python ISOVIZ_INPUT/main_requirement.py
  python ISOVIZ_INPUT/main_requirement.py --dev
"""
from __future__ import annotations

import argparse
import json
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
) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        check=check,
        cwd=cwd,
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
        print(f"[pip] Skip (no extra packages in {req_file.name}; stdlib is enough).")
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
    print(f"[python] {sys.version.split()[0]}")


def _check_java() -> None:
    java = shutil.which("java") or shutil.which("javaw")
    if java is None:
        print("[java] WARNING: java/javaw not in PATH. IsoVIZ is a Java application.")
        print("        Install a JRE/JDK, or set ISOVIZ / ISOVIZ_JAR to the IsoVIZ launcher.")
        return
    cp = subprocess.run(
        [java, "-version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    first = (cp.stdout or "").strip().splitlines()[:1]
    print(f"[java] {java}" + (f" ({first[0]})" if first else ""))


def _check_isoviz_shortcut(project_root: Path) -> None:
    from isoviz_input.launcher import find_isoviz_launcher

    found = find_isoviz_launcher(project_root)
    if found is None:
        print("[isoviz] WARNING: no ISOViz.lnk / IsoViz.exe / IsoViz.jar in the CRIS root.")
        print("         Place the IsoVIZ shortcut there, or set environment variable ISOVIZ.")
        return
    print(f"[isoviz] launcher = {found}")


def _check_output(isoviz_root: Path) -> None:
    out_dir = isoviz_root / "output"
    if not out_dir.exists():
        print(f"[paths] Creating missing folder: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[paths] OK: output folder = {out_dir}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recreate", action="store_true", help="Force recreate CRIS/.venv")
    parser.add_argument("--dev", action="store_true", help="Also install requirements-dev.txt")
    args = parser.parse_args()

    isoviz_root = Path(__file__).resolve().parent
    project_root = isoviz_root.parent
    if not (project_root / "ISODISTORT").is_dir():
        raise RuntimeError("This script must be run from the CRIS repo (ISOVIZ_INPUT/ next to ISODISTORT/).")

    sys.path.insert(0, str(isoviz_root))
    print(f"[root] {project_root}")
    _check_python()
    _check_output(isoviz_root)
    _check_java()

    venv_dir, python = _ensure_venv(project_root, recreate=args.recreate)
    _pip_install_missing(python, isoviz_root / "requirements.txt")
    if args.dev:
        _pip_install_missing(python, isoviz_root / "requirements-dev.txt")

    cp = _run(
        [
            str(python),
            "-c",
            "from isoviz_input.amplitudes import read_amplitude_csv; print('ISOVIZ_INPUT import OK')",
        ],
        cwd=str(isoviz_root),
    )
    print(f"[check] {(cp.stdout or '').strip()}")
    _check_isoviz_shortcut(project_root)

    print("\n=== DONE ===")
    if _is_windows():
        print(f"Use venv python: {python}")
        print(f"  {python} ISOVIZ_INPUT\\main.py --data <csv> --structure <file.isoviz>")
    else:
        print(f"Use venv python: {python}")
        print(f"  {python} ISOVIZ_INPUT/main.py --data <csv> --structure <file.isoviz>")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
