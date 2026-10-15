"""Validate the installed console entrypoint independently of checkout paths."""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import pytest

PROJECT = Path(__file__).resolve().parents[1]


def test_wheel_runs_with_packaged_settings_outside_checkout(tmp_path):
    pytest.importorskip("setuptools", reason="Wheel build requires the declared development dependencies")
    pytest.importorskip("wheel", reason="Wheel build requires the declared development dependencies")
    source = tmp_path / "source"
    source.mkdir()
    for directory in ("isodistort_validate", "config"):
        shutil.copytree(
            PROJECT / directory, source / directory,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
    for filename in ("pyproject.toml", "setup.py", "MANIFEST.in", "README.md", "requirements.txt"):
        shutil.copyfile(PROJECT / filename, source / filename)
    environment = os.environ.copy()
    environment.update(TEMP=str(tmp_path), TMP=str(tmp_path), PYTHONDONTWRITEBYTECODE="1")
    wheels = tmp_path / "wheels"
    wheels.mkdir()
    build = subprocess.run(
        [sys.executable, "-c", "from setuptools.build_meta import build_wheel; "
         "import sys; build_wheel(sys.argv[1])", str(wheels)],
        cwd=source, env=environment, capture_output=True, text=True, check=False,
    )
    assert build.returncode == 0, build.stdout + build.stderr
    wheel, = wheels.glob("*.whl")
    installed = tmp_path / "installed"
    with zipfile.ZipFile(wheel) as archive:
        assert archive.read("isodistort_validate/_resources/settings.yaml") == (
            PROJECT / "config/settings.yaml"
        ).read_bytes()
        archive.extractall(installed)
    # A generic config directory from another installed package cannot override
    # this distribution's own bundled settings.
    (installed / "config").mkdir()
    (installed / "config/settings.yaml").write_text("unrelated: true\n", encoding="utf-8")
    script = """
import hashlib
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from importlib.metadata import distribution
from pymatgen.core import Lattice, Structure
from isodistort_validate.backend.config import CONFIG_PATH, get_config
from isodistort_validate import config_loader as legacy_config
from isodistort_validate.backend.config import service as canonical_config
assert legacy_config is canonical_config
assert CONFIG_PATH.parent.name == '_resources'
assert CONFIG_PATH.is_relative_to(Path(sys.argv[1]))
cfg = get_config()
assert cfg.compare_root == Path(sys.argv[1]) / 'compare'
entry = next(ep for ep in distribution('isodistort-validate').entry_points if ep.name == 'isodistort-validate')
main = entry.load()
assert main(['--help']) == 0
# Synthetic files only exercise installed orchestration, not external scientific acceptance.
structure = Structure(Lattice.from_parameters(4.1, 5.2, 6.3, 81, 93, 104),
                      ['Si', 'O'], [[.1, .2, .3], [.27, .36, .41]])
local = cfg.item_dir / 'packaging.cif'
reference = cfg.true_dir / 'packaging.cif'
structure.to(filename=str(local), fmt='cif')
reference.write_bytes(local.read_bytes())
before = [hashlib.sha256(path.read_bytes()).digest() for path in (local, reference)]
assert main(['compare', 'packaging.cif', '--json']) == 0
assert main(['batch', '--json']) == 0
assert before == [hashlib.sha256(path.read_bytes()).digest() for path in (local, reference)]
"""
    result = subprocess.run(
        [sys.executable, "-I", "-c", script, str(installed)],
        cwd=installed, env=environment, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert '"passed": true' in result.stdout
