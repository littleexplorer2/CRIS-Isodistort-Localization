"""Keep folder boundaries and dependency direction executable."""

from __future__ import annotations

import ast
import importlib
import subprocess
import sys
from importlib.util import resolve_name
from pathlib import Path

import pytest

PROJECT = Path(__file__).resolve().parents[1]
PACKAGE = PROJECT / "isodistort_validate"
ALLOWED = {
    "backend": {"backend"},
    "application": {"application", "backend"},
    "frontend": {"frontend", "application"},
}


def _imports(path: Path):
    relative = path.relative_to(PROJECT).with_suffix("")
    parts = list(relative.parts)
    context = ".".join(parts[:-1])
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            base = "." * node.level + (node.module or "")
            base = resolve_name(base, context) if node.level else base
            # `from . import service` and `from ...backend import paths`
            # import modules, whereas imports from a module import its API.
            if len(base.split(".")) <= 2 or node.module is None:
                yield from (f"{base}.{alias.name}" for alias in node.names)
            else:
                yield base


def test_layers_use_module_folders_and_directed_public_imports() -> None:
    violations = []
    for layer, allowed in ALLOWED.items():
        layer_root = PACKAGE / layer
        assert (layer_root / "__init__.py").is_file()
        assert not [path for path in layer_root.glob("*.py") if path.name != "__init__.py"]
        modules = [path for path in layer_root.iterdir() if path.is_dir() and path.name != "__pycache__"]
        assert modules
        assert all((path / "__init__.py").is_file() for path in modules)
        for path in layer_root.rglob("*.py"):
            owner = path.relative_to(PROJECT).parts[:3]
            for dependency in _imports(path):
                parts = dependency.split(".")
                if parts[0] != "isodistort_validate":
                    continue
                if len(parts) < 2 or parts[1] not in allowed:
                    violations.append(f"{path.relative_to(PROJECT)} -> {dependency}")
                elif len(parts) > 3 and tuple(parts[:3]) != owner:
                    violations.append(f"private module import: {path.relative_to(PROJECT)} -> {dependency}")
    assert not violations, "\n".join(violations)


@pytest.mark.parametrize("name", ["compare_cif", "batch_compare", "compare_paths", "config_loader"])
def test_legacy_modules_only_forward(name: str) -> None:
    tree = ast.parse((PACKAGE / f"{name}.py").read_text(encoding="utf-8"))
    assert not any(isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) for node in ast.walk(tree))


def test_legacy_callers_share_implementations_and_mutable_paths(monkeypatch, tmp_path) -> None:
    old_paths = importlib.import_module("isodistort_validate.compare_paths")
    old_config = importlib.import_module("isodistort_validate.config_loader")
    old_single = importlib.import_module("isodistort_validate.compare_cif")
    old_batch = importlib.import_module("isodistort_validate.batch_compare")
    from isodistort_validate.application import workspace
    from isodistort_validate.application.batch import run_batch
    from isodistort_validate.application.single import compare_cif
    from isodistort_validate.backend.config import get_config
    from isodistort_validate.backend.paths import service

    assert old_paths is service
    assert old_config.get_config is get_config
    assert old_single.compare_cif is compare_cif
    assert old_batch.run_batch is run_batch
    monkeypatch.setattr(old_paths, "ITEM_DIR", tmp_path / "item")
    monkeypatch.setattr(old_paths, "TRUE_DIR", tmp_path / "true")
    assert workspace.resolve_pair("nested/a.cif") == (
        tmp_path / "item" / "nested" / "a.cif",
        tmp_path / "true" / "nested" / "a.cif",
    )


def test_importing_comparison_workflows_does_not_load_frontend() -> None:
    process = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; import isodistort_validate.application.single; "
                "import isodistort_validate.application.batch; "
                "assert not any(name.startswith('isodistort_validate.frontend') for name in sys.modules)"
            ),
        ],
        cwd=PROJECT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert process.returncode == 0, process.stderr


def test_main_is_a_thin_frontend_launcher() -> None:
    tree = ast.parse((PROJECT / "main.py").read_text(encoding="utf-8"))
    assert not any(isinstance(node, (ast.FunctionDef, ast.ClassDef)) for node in ast.walk(tree))
    assert "isodistort_validate.frontend.cli" in set(_imports(PROJECT / "main.py"))
