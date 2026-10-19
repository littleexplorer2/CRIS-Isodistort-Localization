from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

# CRIS subprojects also expose main_requirement.py. Bind the installer under
# test by its owned path, independent of pytest's cross-repository sys.path.
_installer_path = Path(__file__).resolve().parents[1] / "main_requirement.py"
_spec = importlib.util.spec_from_file_location("gd_environment_installer", _installer_path)
assert _spec is not None and _spec.loader is not None
main_requirement = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(main_requirement)


def test_run_strips_parent_python_environment(monkeypatch) -> None:
    for name in ("PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV", "__PYVENV_LAUNCHER__"):
        monkeypatch.setenv(name, f"sentinel-{name.lower()}")

    completed = main_requirement._run(
        [
            sys.executable,
            "-c",
            (
                "import os; "
                "names=('PYTHONHOME','PYTHONPATH','VIRTUAL_ENV',"
                "'__PYVENV_LAUNCHER__'); "
                "print('|'.join(os.environ.get(name, '<missing>') for name in names))"
            ),
        ]
    )

    assert completed.stdout.strip() == "<missing>|<missing>|<missing>|<missing>"
