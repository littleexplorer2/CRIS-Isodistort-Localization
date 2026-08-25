"""Locate and launch the Java IsoVIZ visualizer."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def cris_root() -> Path:
    return Path(__file__).resolve().parents[2]


def find_isoviz_launcher(root: Path | None = None) -> Path | None:
    """Return a shortcut, executable, or JAR if one is configured."""
    base = root or cris_root()
    env_path = os.environ.get("ISOVIZ") or os.environ.get("ISOVIZ_JAR")
    if env_path:
        candidate = Path(env_path).expanduser()
        if candidate.is_file():
            return candidate
    for name in ("ISOViz.lnk", "IsoVIZ.lnk", "ISOVIZ.lnk", "IsoViz.exe", "ISOViz.exe"):
        hit = base / name
        if hit.is_file():
            return hit
    for name in ("IsoViz.jar", "ISOViz.jar", "isoviz.jar"):
        hit = base / name
        if hit.is_file():
            return hit
    return None


def java_executable() -> str | None:
    return shutil.which("javaw") or shutil.which("java")


def open_isoviz(isoviz_file: Path, *, launcher: Path | None = None) -> None:
    """Open a ``.isoviz`` file in IsoVIZ."""
    target = isoviz_file.resolve()
    if not target.is_file():
        raise FileNotFoundError(f"Patched IsoVIZ file not found: {target}")
    app = launcher or find_isoviz_launcher()
    if sys.platform.startswith("win"):
        try:
            os.startfile(str(target))  # type: ignore[attr-defined]
            return
        except OSError:
            pass
        if app is not None and app.suffix.lower() == ".jar":
            _run_jar(app, target)
            return
        if app is not None:
            subprocess.Popen(["cmd", "/c", "start", "", str(app), str(target)], close_fds=True)
            return
        raise RuntimeError(
            "Windows could not open the .isoviz file. Associate it with IsoVIZ, "
            "or place ISOViz.lnk in the CRIS root, or set ISOVIZ / ISOVIZ_JAR."
        )
    if app is not None and app.suffix.lower() == ".jar":
        _run_jar(app, target)
        return
    if app is not None:
        subprocess.Popen([str(app), str(target)], close_fds=True)
        return
    raise RuntimeError(
        "IsoVIZ was not found. Install the Java IsoVIZ tool, or set ISOVIZ / ISOVIZ_JAR "
        "to the .exe / .jar path."
    )


def _run_jar(jar: Path, isoviz_file: Path) -> None:
    java = java_executable()
    if not java:
        raise RuntimeError("Java is required to launch IsoVIZ from a .jar file.")
    subprocess.Popen([java, "-jar", str(jar), str(isoviz_file)], close_fds=True)
