from __future__ import annotations

import sys

import main_requirement


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
