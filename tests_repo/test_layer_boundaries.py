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
        if project.get("role") == "test_support" and (
            target == "tests" or target.startswith("tests.")
        ) and not (target == prefix or target.startswith(prefix + ".")):
            errors.append(f"{module}:{line}: read-only comparison imports external test/helper code {target}")
        if any(target == forbidden or target.startswith(forbidden + ".")
               for forbidden in project.get("forbidden_imports", [])):
            errors.append(f"{module}:{line}: {project.get('role')} imports forbidden production role {target}")
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


def _private_session_accesses(source: str, factories: list[str]) -> list[str]:
    """Track registered session constructors, helper returns, and local aliases.

    This intentionally does not ban arbitrary underscore attributes. Unknown
    dynamic values remain unknown; the gate is a static regression check, not a
    proof about eval, reflection, runtime monkeypatching, or arbitrary data flow.
    """
    errors: list[str] = []
    function_results: dict[str, object] = {}
    factory_names = set(factories)

    class Visitor(ast.NodeVisitor):
        def __init__(self, names=None):
            self.names = dict(names or {})
            self.return_values = []

        def symbol(self, node):
            if isinstance(node, ast.Name):
                return self.names.get(node.id)
            if isinstance(node, ast.Attribute):
                base = self.symbol(node.value)
                return f"{base}.{node.attr}" if isinstance(base, str) else None
            if isinstance(node, (ast.Tuple, ast.List)):
                return tuple(self.symbol(element) for element in node.elts)
            if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant):
                value, index = self.symbol(node.value), node.slice.value
                if isinstance(value, tuple) and isinstance(index, int) and -len(value) <= index < len(value):
                    return value[index]
            if isinstance(node, ast.Call):
                target = self.symbol(node.func)
                if isinstance(target, str) and target in factory_names:
                    return "instance:" + target
                if isinstance(target, str):
                    return function_results.get(target)
            return None

        def bind(self, target, value):
            if isinstance(target, ast.Name):
                self.names[target.id] = value
            elif isinstance(target, (ast.Tuple, ast.List)):
                values = value if isinstance(value, tuple) and len(value) == len(target.elts) else [None] * len(target.elts)
                for child, item in zip(target.elts, values):
                    self.bind(child, item)

        def visit_ImportFrom(self, node):
            for alias in node.names:
                self.names[alias.asname or alias.name] = f"{node.module}.{alias.name}"

        def visit_Import(self, node):
            for alias in node.names:
                self.names[alias.asname or alias.name.split('.')[0]] = alias.name if alias.asname else alias.name.split('.')[0]

        def visit_Assign(self, node):
            self.visit(node.value)
            value = self.symbol(node.value)
            for target in node.targets:
                self.visit(target)
                self.bind(target, value)

        def visit_AnnAssign(self, node):
            if node.value is not None:
                self.visit(node.value)
            self.visit(node.target)
            self.bind(node.target, self.symbol(node.value))

        def visit_FunctionDef(self, node):
            child = Visitor(self.names)
            for arg in (*node.args.posonlyargs, *node.args.args, *node.args.kwonlyargs):
                annotation = child.symbol(arg.annotation)
                child.names[arg.arg] = "instance:" + annotation if annotation in factory_names else None
            for statement in node.body:
                child.visit(statement)
            key = f"function:{node.name}:{node.lineno}"
            values = child.return_values
            if values and all(value == values[0] for value in values):
                function_results[key] = values[0]
            self.names[node.name] = key

        visit_AsyncFunctionDef = visit_FunctionDef

        def visit_Return(self, node):
            self.visit(node.value) if node.value is not None else None
            self.return_values.append(self.symbol(node.value))

        def check_member(self, node, owner, member):
            if (isinstance(owner, str) and owner.startswith("instance:")
                    and isinstance(member, str) and member.startswith("_")
                    and not (member.startswith("__") and member.endswith("__"))):
                errors.append(f"line {node.lineno}: external consumer accesses private session member {member}")

        def visit_Attribute(self, node):
            self.check_member(node, self.symbol(node.value), node.attr)
            self.generic_visit(node)

        def visit_Call(self, node):
            if isinstance(node.func, ast.Name) and node.func.id in {"getattr", "setattr", "delattr", "hasattr"}:
                if len(node.args) >= 2 and isinstance(node.args[1], ast.Constant):
                    self.check_member(node, self.symbol(node.args[0]), node.args[1].value)
            self.generic_visit(node)

    Visitor().visit(ast.parse(source))
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


@pytest.mark.parametrize("source", [
    "from application.session import IsoDistort",
    "from backend.input_cif import read_cif",
    "import frontend.web.server",
    'importlib.import_module("application.session")',
    "from ISODISTORT.application.session import IsoDistort",
    "from tests.validate_method4_local import main",
    "from tests.support_backend_export import _export_data",
    "from tests import prepare_method4_inputs",
])
def test_read_only_validation_role_rejects_production_imports(source):
    assert _violations(source, "tests.validation.application.single.service", False, "application", PROJECTS[1])


@pytest.mark.parametrize("source", [
    "from tests.validation.backend.provenance import FileSnapshot",
    "from ...backend.results import VerificationPolicy",
    "from pymatgen.io.cif import CifParser",
    "import importlib.metadata",
])
def test_read_only_validation_role_allows_own_contracts_and_external_parsers(source):
    assert not _violations(source, "tests.validation.application.single.service", False, "application", PROJECTS[1])


def test_isodistort_checks_have_one_script_per_registered_module():
    """Keep check entry points aligned with production ownership."""
    checks = ROOT / "ISODISTORT" / "tests"
    for layer, modules in PROJECTS[0]["modules"].items():
        folder = checks / layer
        assert folder.is_dir(), folder
        assert {path.name for path in folder.glob("test_*.py")} == {
            f"test_{module}.py" for module in modules
        }
        assert not [path for path in folder.rglob("test_*.py") if path.parent != folder]
    assert not (checks / "manual").exists()


@pytest.mark.parametrize("source", [
    "from application.session import IsoDistort as Session\ns = Session()\ns._private()",
    "from application.session import IsoDistort\nFactory = IsoDistort\ns = Factory()\nalias = s\nalias._state",
    "import application.session as api\ns = api.IsoDistort()\ngetattr(s, '_state')",
    "def load():\n from application.session import IsoDistort\n return IsoDistort, None\nFactory, normalizer = load()\ns = Factory()\ns._private()",
    "from application.session import IsoDistort\ndef use(session: IsoDistort):\n alias = session\n alias._private()",
])
def test_external_consumer_gate_tracks_session_aliases_and_loader_returns(source):
    assert _private_session_accesses(source, ["application.session.IsoDistort"])


@pytest.mark.parametrize("source", [
    "from application.session import IsoDistort\ns = IsoDistort()\ns.export_mode_result()",
    "from application.session import IsoDistort\ns = IsoDistort()\nother = object()\nother._private()",
    "from application.session import IsoDistort\ns = IsoDistort()\ns = object()\ns._private()",
    "class Other:\n def method(self):\n  return self._private()",
])
def test_external_consumer_gate_preserves_public_contracts_and_unrelated_private_state(source):
    assert not _private_session_accesses(source, ["application.session.IsoDistort"])


@pytest.mark.parametrize("consumer", MANIFEST.get("external_consumers", []), ids=lambda item: item["path"])
def test_registered_external_consumers_use_public_session_contracts(consumer):
    path = (ROOT / consumer["path"]).resolve()
    if not path.is_file() and not consumer.get("required", True):
        pytest.skip(f"Optional sibling consumer is unavailable: {path}")
    assert path.is_file(), path
    if consumer.get("forbid_private_session_members"):
        errors = _private_session_accesses(path.read_text(encoding="utf-8-sig"), consumer["session_factories"])
        assert not errors, "\n".join(f"{path}: {error}" for error in errors)
