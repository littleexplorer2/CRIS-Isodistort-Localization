"""Guard ISODISTORT entry points against the OneDrive/WSL process restriction."""

from __future__ import annotations

import sys
from pathlib import Path


def require_cris_runner(script: str | Path) -> None:
    """Reject a direct Windows venv process that cannot create WSL children."""
    if not sys.platform.startswith("win"):
        return
    base_executable = str(getattr(sys, "_base_executable", "") or "").strip()
    if not base_executable:
        return
    try:
        current = Path(sys.executable).resolve()
        base = Path(base_executable).resolve()
        script_path = Path(script).resolve()
        repository_root = script_path.parent.parent
        venv_root = (repository_root / ".venv").resolve()
        current.relative_to(venv_root)
    except (OSError, ValueError):
        return
    if current == base:
        return
    relative_script = script_path.relative_to(repository_root)
    display_script = str(relative_script).replace("/", "\\")
    print(
        "This OneDrive checkout cannot run WSL from .venv\\Scripts\\python.exe.\n"
        "Start the same physical CRIS environment through the repository runner:\n"
        f"  .\\run_cris.ps1 {display_script}",
        file=sys.stderr,
        flush=True,
    )
    raise SystemExit(2)
