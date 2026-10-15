"""Exercise shipped entrypoints and wheel resources outside the source tree."""
from __future__ import annotations

import configparser
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


def _subprocess_environment(tmp_path: Path) -> dict[str, str]:
    environment = os.environ.copy()
    environment["TEMP"] = str(tmp_path)
    environment["TMP"] = str(tmp_path)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return environment


def test_repository_cli_prepares_copy_and_preserves_sources(tmp_path):
    sources = [FIXTURES / "sample.csv", FIXTURES / "sample.isoviz"]
    before = [hashlib.sha256(path.read_bytes()).digest() for path in sources]

    result = subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "main.py"), "--cli",
         "--data", str(sources[0]), "--structure", str(sources[1]), "--prepare-only"],
        cwd=tmp_path, env=_subprocess_environment(tmp_path),
        capture_output=True, text=True, check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "[matched]   2 mode(s)" in result.stdout
    generated = list(tmp_path.glob("*.isoviz"))
    assert len(generated) == 1
    assert "0.12345" in generated[0].read_text(encoding="utf-8")
    assert before == [hashlib.sha256(path.read_bytes()).digest() for path in sources]


def test_wheel_contains_module_packages_config_and_working_entrypoint(tmp_path):
    pytest.importorskip("setuptools", reason="Wheel build requires the declared development dependencies")
    pytest.importorskip("wheel", reason="Wheel build requires the declared development dependencies")
    source = tmp_path / "source"
    source.mkdir()
    for directory in ("isoviz_input", "config"):
        shutil.copytree(
            PROJECT_ROOT / directory, source / directory,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
    for filename in ("pyproject.toml", "setup.py", "MANIFEST.in", "README.md", "requirements.txt"):
        shutil.copyfile(PROJECT_ROOT / filename, source / filename)
    wheels = tmp_path / "wheels"
    wheels.mkdir()
    build = subprocess.run(
        [sys.executable, "-c", "from setuptools.build_meta import build_wheel; "
         "import sys; build_wheel(sys.argv[1])", str(wheels)],
        cwd=source, env=_subprocess_environment(tmp_path),
        capture_output=True, text=True, check=False,
    )
    assert build.returncode == 0, build.stdout + build.stderr
    wheel, = wheels.glob("*.whl")
    installed = tmp_path / "installed"
    with zipfile.ZipFile(wheel) as archive:
        assert archive.read("isoviz_input/_resources/settings.yaml") == (
            PROJECT_ROOT / "config/settings.yaml"
        ).read_bytes()
        names = set(archive.namelist())
        for layer, modules in {
            "backend": ("amplitudes", "inputs", "validation", "paths", "launcher", "config_loader"),
            "application": ("pipeline", "input_selection", "environment"),
            "frontend": ("cli", "gui", "entrypoint"),
        }.items():
            for module in modules:
                assert f"isoviz_input/{layer}/{module}/__init__.py" in names
        entrypoints = next(name for name in names if name.endswith(".dist-info/entry_points.txt"))
        config = configparser.ConfigParser()
        config.read_string(archive.read(entrypoints).decode("utf-8"))
        assert config["console_scripts"]["isoviz-input"] == "isoviz_input.frontend.entrypoint:main"
        archive.extractall(installed)

    # An unrelated package's generic config directory must not override the
    # resource bundled with this installed distribution.
    (installed / "config").mkdir()
    (installed / "config/settings.yaml").write_text("unrelated: true\n", encoding="utf-8")

    # Isolated Python plus an explicit wheel path prevents the checkout from
    # hiding a missing package, resource, import or console-entrypoint defect.
    script = """
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from importlib.metadata import distribution
from isoviz_input.backend.config_loader import CONFIG_PATH
from isoviz_input import amplitudes as legacy
from isoviz_input.backend import amplitudes
assert CONFIG_PATH.is_relative_to(Path(sys.argv[1]))
assert CONFIG_PATH.parent.name == '_resources'
assert legacy.apply_amplitudes is amplitudes.apply_amplitudes
entry = next(ep for ep in distribution('isoviz-input').entry_points if ep.name == 'isoviz-input')
assert entry.load()(['--help']) == 0
raise SystemExit(entry.load()(['--cli', '--data', sys.argv[2], '--structure', sys.argv[3], '--prepare-only']))
"""
    result = subprocess.run(
        [sys.executable, "-I", "-c", script, str(installed),
         str(FIXTURES / "sample.csv"), str(FIXTURES / "sample.isoviz")],
        cwd=installed, env=_subprocess_environment(tmp_path),
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "[matched]   2 mode(s)" in result.stdout
