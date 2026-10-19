"""Verify the three-layer package through an actual, isolated wheel install."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1] / "ISODISTORT"


def test_isodistort_wheel_contains_only_canonical_layers(tmp_path):
    pytest.importorskip("setuptools", reason="Wheel build requires setuptools and wheel")
    pytest.importorskip("wheel", reason="Wheel build requires setuptools and wheel")
    source = tmp_path / "source"
    source.mkdir()
    for name in ("pyproject.toml", "README.md", "requirements.txt"):
        shutil.copyfile(ROOT / name, source / name)
    for name in ("backend", "application", "frontend", "resources"):
        shutil.copytree(ROOT / name, source / name, ignore=shutil.ignore_patterns("__pycache__", "*.log"))
    wheels = tmp_path / "wheels"
    wheels.mkdir()
    build = subprocess.run(
        [sys.executable, "-B", "-c",
         "from setuptools.build_meta import build_wheel; import sys; build_wheel(sys.argv[1])",
         str(wheels)], cwd=source, capture_output=True, text=True, timeout=120,
    )
    assert build.returncode == 0, build.stdout + build.stderr
    installed = tmp_path / "installed"
    with zipfile.ZipFile(next(wheels.glob("*.whl"))) as archive:
        archive.extractall(installed)
    probe = """
import sys
from pathlib import Path
installed = Path(sys.argv[1])
sys.path[:0] = [str(installed), *sys.argv[2:]]
import application.session
from application.session import core_api
from backend.export.distortion_formats import cart_normalized_mode_matrix
from backend.utils.config_loader import CONFIG_PATH, get_config
from frontend.web.server import WEB_DIR
from application import workflows
assert Path(core_api.__file__).is_relative_to(installed)
assert callable(cart_normalized_mode_matrix)
for retired in ('features', 'frontend/i18n', 'backend/api', 'tests',
                'backend/method1/search_methods.py', 'backend/method1/affine_embeddings.py',
                'backend/method1/phase_path.py'):
    assert not (installed / retired).exists()
assert CONFIG_PATH == installed / 'resources/config/settings.yaml'
assert get_config().iso_bin.is_file()
assert (WEB_DIR / 'index.html').is_file()
assert (WEB_DIR / 'static/bootstrap.css').is_file()
assert workflows.messages()['load.done']
assert list(installed.glob('*.dist-info/entry_points.txt'))
print('installed canonical layers and runtime resources verified; retired packages absent')
"""
    dependencies = [path for path in sys.path if Path(path).name == "site-packages"]
    env = os.environ.copy()
    env["ISODISTORT_RUNTIME_ROOT"] = str(tmp_path / "runtime")
    result = subprocess.run(
        [sys.executable, "-I", "-S", "-B", "-c", probe, str(installed), *dependencies],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=45,
    )
    assert result.returncode == 0, result.stdout + result.stderr
