"""Regression checks for the post-migration package and resource layout."""

from __future__ import annotations

import ast
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path, PurePosixPath

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
    import tomli as tomllib

import pytest

from backend.api import core_api
from backend.tables import irreps_cdml
from backend.utils.config_loader import CONFIG_PATH, PROJECT_ROOT, get_config

CRIS_ROOT = Path(__file__).resolve().parents[2]
ISODISTORT_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT_PATH = ISODISTORT_ROOT / "pyproject.toml"


def test_production_modules_do_not_silently_overwrite_top_level_definitions():
    """A duplicate helper definition can reverse an existing scientific API."""
    duplicates = []
    for package in ("backend", "features", "frontend"):
        for path in (ISODISTORT_ROOT / package).rglob("*.py"):
            tree = ast.parse(path.read_text(encoding="utf-8-sig"))
            counts = Counter(
                node.name for node in tree.body
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
            )
            duplicates.extend(
                (str(path.relative_to(ISODISTORT_ROOT)), name)
                for name, count in counts.items() if count > 1
            )
    assert not duplicates, duplicates

FRONTEND_PACKAGE_DATA = {
    "index.html",
    "static/bootstrap.css",
    "static/docs.css",
    "static/help.jpg",
}
RESOURCE_PACKAGE_DATA = {
    "config/settings.yaml",
    "isobyu/comsubs",
    "isobyu/comsubs.txt",
    "isobyu/comsubs_sample.in",
    "isobyu/const.dat",
    "isobyu/data_diperiodic.txt",
    "isobyu/data_images.txt",
    "isobyu/data_irreps.txt",
    "isobyu/data_isotropy.txt",
    "isobyu/data_little.txt",
    "isobyu/data_magnetic.txt",
    "isobyu/data_space.txt",
    "isobyu/data_ssg.txt",
    "isobyu/data_ssgmag.txt",
    "isobyu/data_wyckoff.txt",
    "isobyu/findsym",
    "isobyu/findsym.txt",
    "isobyu/findsym_cifinput",
    "isobyu/findsym_sample.in",
    "isobyu/iso",
    "isobyu/smodes",
    "isobyu/smodes.txt",
    "isobyu/smodes_sample.in",
}


def _load_pyproject() -> dict:
    with PYPROJECT_PATH.open("rb") as stream:
        return tomllib.load(stream)


def _matched_package_files(package_root: Path, patterns: list[str]) -> set[str]:
    return {
        path.relative_to(package_root).as_posix()
        for pattern in patterns
        for path in package_root.glob(pattern)
        if path.is_file()
    }


def _matches_any_package_pattern(relative_path: str, patterns: list[str]) -> bool:
    path = PurePosixPath(relative_path)
    return any(path.match(pattern) for pattern in patterns)


def test_reorganized_resource_roots_resolve_to_existing_project_paths() -> None:
    assert PROJECT_ROOT == ISODISTORT_ROOT
    assert CONFIG_PATH == ISODISTORT_ROOT / "resources" / "config" / "settings.yaml"
    assert irreps_cdml._ISO_DATA == (
        ISODISTORT_ROOT / "resources" / "isobyu" / "data_irreps.txt"
    )
    assert irreps_cdml._ISO_DATA.is_file()
    assert core_api._CRIS_ROOT == CRIS_ROOT
    assert (core_api._CRIS_ROOT / "experiment_data").is_dir()
    assert (core_api._CRIS_ROOT / "webpage_info").is_dir()


def test_output_paths_remain_outside_read_only_resources() -> None:
    cfg = get_config()
    assert cfg.resolve_path(cfg._cfg["runtime"]["output_dir"]) == (
        ISODISTORT_ROOT / "output"
    )
    assert cfg.resolve_path(cfg._cfg["runtime"]["temp_dir"]) == (
        ISODISTORT_ROOT / "output" / "tmp"
    )


@pytest.mark.parametrize(
    "module",
    [
        "features.method2.superspace",
        "features.method3.coupled_routes",
        "features.method3.inverse_landau_adapter",
    ],
)
def test_cross_feature_modules_import_first_in_a_clean_interpreter(module: str) -> None:
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(  # noqa: S603 - fixed interpreter and parametrized modules
        [sys.executable, "-B", "-c", f"import {module}"],
        cwd=ISODISTORT_ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr


def test_method1_public_api_does_not_reexport_other_feature_packages() -> None:
    import features.method1 as method1
    from features.method2 import DistortionMapper, ParametricModeResult
    from features.method3 import CoupledRouteWitness, FiniteGroup
    from features.method4 import DistortionEngine, HomogeneousStrainMode

    assert "DistortionMapper" not in method1.__all__
    assert "DistortionEngine" not in method1.__all__
    assert DistortionMapper.__module__ == "features.method2.distortion_mapper"
    assert ParametricModeResult.__module__ == "features.method2.superspace"
    assert CoupledRouteWitness.__module__ == "features.method3.coupled_routes"
    assert FiniteGroup.__module__ == "features.method3.inverse_landau"
    assert DistortionEngine.__module__ == "features.method4.distortion_engine"
    assert HomogeneousStrainMode.__module__ == "features.method4.strain_modes"


def test_setuptools_package_data_contract_is_complete_and_explicit() -> None:
    pyproject = _load_pyproject()
    setuptools = pyproject["tool"]["setuptools"]
    package_data = setuptools["package-data"]

    assert setuptools["include-package-data"] is False
    package_find = setuptools["packages"]["find"]
    assert "resources" in package_find["include"]
    assert package_find["namespaces"] is False
    assert (ISODISTORT_ROOT / "resources" / "__init__.py").is_file()

    frontend_patterns = package_data["frontend.web"]
    resource_patterns = package_data["resources"]
    assert _matched_package_files(
        ISODISTORT_ROOT / "frontend" / "web", frontend_patterns
    ) == FRONTEND_PACKAGE_DATA
    assert _matched_package_files(
        ISODISTORT_ROOT / "resources", resource_patterns
    ) == RESOURCE_PACKAGE_DATA

    for excluded in (
        "isobyu/comsubs_sample.log",
        "isobyu/findsym_sample.log",
        "isobyu/iso.log",
        "output/generated.cif",
    ):
        assert not _matches_any_package_pattern(excluded, resource_patterns)


def test_equivalent_installed_layout_resolves_packaged_runtime_assets(tmp_path: Path) -> None:
    """Exercise the package-data contract without requiring a network build backend."""
    pyproject = _load_pyproject()
    package_data = pyproject["tool"]["setuptools"]["package-data"]
    site_packages = tmp_path / "site-packages"

    for package_name in ("backend", "features", "frontend", "resources"):
        source_root = ISODISTORT_ROOT / package_name
        for source in source_root.rglob("*.py"):
            destination = site_packages / source.relative_to(ISODISTORT_ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)

    for package_name, patterns in package_data.items():
        package_root = ISODISTORT_ROOT.joinpath(*package_name.split("."))
        for pattern in patterns:
            for source in package_root.glob(pattern):
                if not source.is_file():
                    continue
                destination = site_packages / source.relative_to(ISODISTORT_ROOT)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)

    probe = """
import os
import sys
from pathlib import Path

site_packages = Path(sys.argv[1]).resolve()
for dependency_path in reversed(sys.argv[2:]):
    sys.path.insert(0, dependency_path)
sys.path.insert(0, str(site_packages))

import resources
from backend.tables import irreps_cdml
from backend.utils.config_loader import CONFIG_PATH, PROJECT_ROOT, get_config
from frontend.web.server import WEB_DIR

assert Path(resources.__file__).resolve().parent == site_packages / "resources"
assert PROJECT_ROOT == site_packages
assert CONFIG_PATH == site_packages / "resources" / "config" / "settings.yaml"
assert irreps_cdml._ISO_DATA == site_packages / "resources" / "isobyu" / "data_irreps.txt"
assert CONFIG_PATH.is_file()
assert irreps_cdml._ISO_DATA.is_file()
assert (WEB_DIR / "index.html").is_file()
assert (WEB_DIR / "static" / "bootstrap.css").is_file()
assert (WEB_DIR / "static" / "docs.css").is_file()
assert (WEB_DIR / "static" / "help.jpg").is_file()
assert get_config().iso_bin.parent == site_packages / "resources" / "isobyu"
runtime_root = Path(os.environ["ISODISTORT_RUNTIME_ROOT"])
assert get_config().output_dir == runtime_root / "output"
assert get_config().temp_dir == runtime_root / "output" / "tmp"
assert get_config().output_dir.is_dir()
assert get_config().temp_dir.is_dir()
assert not (site_packages / "output").exists()
assert not (site_packages / "resources" / "isobyu" / "iso.log").exists()
assert not (site_packages / "resources" / "output").exists()
"""
    dependency_paths = [path for path in sys.path if Path(path).name == "site-packages"]
    env = os.environ.copy()
    env["ISODISTORT_RUNTIME_ROOT"] = str(tmp_path / "wheel-runtime")
    completed = subprocess.run(  # noqa: S603 - fixed interpreter and local test paths
        [
            sys.executable,
            "-I",
            "-B",
            "-c",
            probe,
            str(site_packages),
            *dependency_paths,
        ],
        cwd=tmp_path,
        env=env,
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
