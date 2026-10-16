"""Executable module inventory and directed dependency checks for all CRIS tools."""
from __future__ import annotations

import ast
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = json.loads(Path(__file__).with_name("architecture_manifest.json").read_text(encoding="utf-8"))
PROJECTS = MANIFEST["projects"]
LAYERS = ("backend", "application", "frontend")
ALLOWED = {
    "backend": {"backend"},
    "application": {"application", "backend"},
    "frontend": {"frontend", "application"},
}


def _module_name(path: Path, project: dict) -> str:
    relative = path.relative_to(ROOT / project["root"]).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(filter(None, [project["namespace"], *parts]))


def _imports(tree: ast.AST, module: str, is_package: bool):
    package = module if is_package else module.rpartition(".")[0]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield node.lineno, alias.name
        elif isinstance(node, ast.ImportFrom):
            target = node.module or ""
            if node.level:
                target = importlib.util.resolve_name("." * node.level + target, package)
            yield node.lineno, target
            for alias in node.names:
                if alias.name != "*":
                    yield node.lineno, target + "." + alias.name
        elif isinstance(node, ast.Call) and node.args:
            name = (node.func.id if isinstance(node.func, ast.Name)
                    else node.func.attr if isinstance(node.func, ast.Attribute) else "")
            if name in {"import_module", "__import__"}:
                value = node.args[0]
                if isinstance(value, ast.Constant) and isinstance(value.value, str):
                    target = value.value
                    if target.startswith("."):
                        target = importlib.util.resolve_name(target, package)
                    yield node.lineno, target


def _violations(source: str, module: str, is_package: bool, layer: str, project: dict) -> list[str]:
    prefix = project["namespace"]
    compatibility = [".".join(filter(None, [prefix, path.removesuffix(".py").replace("/", ".")]))
                     for path in project["compatibility"]]
    compatibility += project.get("retired_imports", [])
    errors = []
    source_local = module.removeprefix(prefix + ".") if prefix else module
    source_module = ".".join(source_local.split(".")[:2])
    for line, target in _imports(ast.parse(source), module, is_package):
        if project.get("role") != "test_support" and (
            target == "tests" or target.startswith("tests.")
        ):
            errors.append(f"{module}:{line}: production imports validation/test code {target}")
        if any(target == legacy or target.startswith(legacy + ".") for legacy in compatibility):
            errors.append(f"{module}:{line}: production imports compatibility path {target}")
            continue
        local = target.removeprefix(prefix + ".") if prefix else target
        if prefix and local == target:
            continue
        target_layer = local.split(".", 1)[0]
        if target_layer in LAYERS and target_layer not in ALLOWED[layer]:
            errors.append(f"{module}:{line}: {layer} imports {target_layer}: {target}")
        parts = local.split(".")
        target_module = ".".join(parts[:2])
        if target_layer in LAYERS and source_module != target_module and len(parts) > 2:
            public_children = project.get("public_submodules", {}).get(target_module, [])
            child = ROOT / project["root"] / parts[0] / parts[1] / parts[2]
            # Python symbols are case-sensitive even on a Windows filesystem:
            # MESSAGES is an exported dict, not the messages.py submodule.
            is_submodule = any(item.name in {parts[2], parts[2] + ".py"}
                               for item in child.parent.iterdir())
            if is_submodule and parts[2] not in public_children:
                errors.append(f"{module}:{line}: imports private implementation {target}")
            if (parts[-1].startswith("_")
                    and local not in project.get("legacy_scientific_interfaces", [])):
                errors.append(f"{module}:{line}: imports private symbol {target}")
    return errors


@pytest.mark.parametrize("project", PROJECTS, ids=lambda item: item["root"])
def test_each_layer_module_has_its_own_registered_package(project):
    root = ROOT / project["root"]
    for layer, modules in project["modules"].items():
        layer_root = root / layer
        assert (layer_root / "__init__.py").is_file()
        actual = {path.name for path in layer_root.iterdir()
                  if path.is_dir() and path.name != "__pycache__"
                  and (path / "__init__.py").is_file()
                  and f"{layer}/{path.name}" not in project["compatibility"]}
        assert actual == set(modules), (layer_root, actual, modules)
        assert not [path for path in layer_root.glob("*.py") if path.name != "__init__.py"]
        for source in layer_root.rglob("*.py"):
            relative = source.relative_to(root).as_posix()
            if any(relative == legacy or relative.startswith(legacy + "/")
                   for legacy in project["compatibility"]):
                continue
            parts = source.relative_to(layer_root).parts
            assert parts == ("__init__.py",) or parts[0] in modules, source
        for module in modules:
            assert (layer_root / module / "__init__.py").is_file()


@pytest.mark.parametrize("project", PROJECTS, ids=lambda item: item["root"])
def test_production_dependencies_follow_the_three_layers(project):
    root = ROOT / project["root"]
    errors = []
    for layer, modules in project["modules"].items():
        paths = [root / layer / "__init__.py"]
        for module in modules:
            paths.extend((root / layer / module).rglob("*.py"))
        for path in paths:
            errors.extend(_violations(path.read_text(encoding="utf-8-sig"),
                                      _module_name(path, project), path.name == "__init__.py",
                                      layer, project))
    assert not errors, "\n".join(errors)


@pytest.mark.parametrize("project", PROJECTS, ids=lambda item: item["root"])
def test_compatibility_paths_contain_no_business_implementation(project):
    for relative in project["compatibility"]:
        path = ROOT / project["root"] / relative
        assert path.exists(), path
        paths = [path] if path.is_file() else path.rglob("*.py")
        for source in paths:
            tree = ast.parse(source.read_text(encoding="utf-8-sig"))
            for node in ast.walk(tree):
                assert not isinstance(node, (ast.ClassDef, ast.AsyncFunctionDef)), source
                if isinstance(node, ast.FunctionDef):
                    assert node.name in {"__getattr__", "__dir__"}, source
                    assert len(node.body) == 1 and isinstance(node.body[0], ast.Return), source
                assert not isinstance(node, (ast.For, ast.While, ast.Try, ast.With)), source


@pytest.mark.parametrize("source,layer,module", [
    ("from application.session import IsoDistort", "backend", "backend.utils.bad"),
    ("from frontend.web import server", "application", "application.session.bad"),
    ("from backend.wrappers import IsoWrapper", "frontend", "frontend.web.bad"),
    ("from ..api import IsoDistort", "backend", "backend.utils.bad"),
    ("import features.method1", "application", "application.session.bad"),
    ('importlib.import_module("backend.api")', "application", "application.session.bad"),
    ("from backend.utils.lattice import _unlisted_helper", "application", "application.session.bad"),
    ("from tests.validation.application.single import compare_cif", "backend", "backend.utils.bad"),
    ('importlib.import_module("tests.validation")', "application", "application.session.bad"),
])
def test_dependency_gate_rejects_reverse_and_legacy_imports(source, layer, module):
    assert _violations(source, module, False, layer, PROJECTS[0])


@pytest.mark.parametrize("source", [
    "from isoviz_input.backend.paths.resolver import resolve_isoviz_path",
    "from isoviz_input.backend.amplitudes import _header_match",
])
def test_dependency_gate_rejects_cross_module_private_imports(source):
    assert _violations(source, "isoviz_input.application.pipeline.service", False,
                       "application", PROJECTS[2])
