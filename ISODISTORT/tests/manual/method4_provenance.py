"""Reproducible provenance fingerprint for live Method 4 validation.

The local report is evidence produced by executable code, runtime dependencies,
ISOTROPY binaries/data, configuration, and frozen parent/daughter inputs. A
hash of the report alone cannot establish that provenance; both producer and
consumer fingerprint the complete immutable input surface independently.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
PROVENANCE_SCHEMA = 1
PROVENANCE_KIND = "method4-local-live-validation-v1"

RUNTIME_DISTRIBUTIONS = ("numpy", "scipy", "pymatgen", "spglib", "PyYAML")
ISO_DATA_BASENAMES = (
    "const.dat",
    "data_diperiodic.txt",
    "data_images.txt",
    "data_irreps.txt",
    "data_isotropy.txt",
    "data_little.txt",
    "data_magnetic.txt",
    "data_space.txt",
    "data_ssg.txt",
    "data_ssgmag.txt",
    "data_wyckoff.txt",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _runtime_state() -> dict[str, Any]:
    distributions: dict[str, str] = {}
    for name in RUNTIME_DISTRIBUTIONS:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
        except Exception as exc:  # noqa: BLE001 - broken metadata remains signable
            error_type = f"{type(exc).__module__}.{type(exc).__qualname__}"
            distributions[name] = f"<unavailable:{error_type}>"

    executable = Path(sys.executable).resolve()
    executable_record = {
        "path": str(executable),
        "exists": executable.is_file(),
        "sha256": sha256_file(executable) if executable.is_file() else None,
        "bytes": executable.stat().st_size if executable.is_file() else None,
    }
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": list(sys.version_info[:5]),
            "cache_tag": sys.implementation.cache_tag,
            "build": sys.version,
            "executable": executable_record,
            "prefix": sys.prefix,
            "base_prefix": sys.base_prefix,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "distributions": distributions,
    }


def _configured_path(config_path: Path, value: object, fallback: str) -> Path:
    raw = value if isinstance(value, str) and value.strip() else fallback
    path = Path(raw)
    if not path.is_absolute():
        path = config_path.parent / path
    return path.resolve()


def _backend_dependencies(package_root: Path) -> list[tuple[str, Path]]:
    """Return immutable ISO inputs while excluding generated logs and caches."""

    config_path = package_root / "resources" / "config" / "settings.yaml"
    isobyu: dict[str, Any] = {}
    try:
        loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict) and isinstance(loaded.get("isobyu"), dict):
            isobyu = loaded["isobyu"]
    except (OSError, UnicodeError, yaml.YAMLError):
        # The missing/malformed config and deterministic fallback paths are
        # both represented in the signed records below.
        pass

    bin_dir = _configured_path(config_path, isobyu.get("bin_dir"), "../isobyu")
    data_dir = _configured_path(config_path, isobyu.get("data_dir"), "../isobyu")
    executable_names: list[tuple[str, str]] = []
    for label, config_key, fallback in (
        ("iso", "iso_bin", "iso"),
        ("findsym", "findsym_bin", "findsym"),
        ("smodes", "smodes_bin", "smodes"),
    ):
        configured = isobyu.get(config_key)
        name = configured if isinstance(configured, str) and configured.strip() else fallback
        executable_names.append((label, name))

    dependencies = [
        (f"backend:{label}-binary", bin_dir / name)
        for label, name in executable_names
    ]
    data_names = set(ISO_DATA_BASENAMES)
    try:
        data_names.update(path.name for path in data_dir.glob("data_*.txt"))
    except OSError:
        pass
    dependencies.extend(
        (f"backend:iso-data:{name}", data_dir / name)
        for name in sorted(data_names, key=str.casefold)
    )
    return dependencies


def _input_paths(
    manifest: Mapping[str, Any],
    *,
    workspace: Path,
) -> list[tuple[str, Path]]:
    paths: list[tuple[str, Path]] = []
    parents = manifest.get("parents")
    if not isinstance(parents, Mapping):
        return [("input:<missing-parents>", workspace / "<missing-parents>")]
    for parent_name, parent in sorted(parents.items(), key=lambda item: str(item[0])):
        if not isinstance(parent, Mapping):
            paths.append((f"parent:{parent_name}:<malformed>", workspace / "<malformed>"))
            continue
        parent_cif = parent.get("parent_cif")
        paths.append(
            (
                f"parent:{parent_name}",
                workspace / str(parent_cif or "<missing-parent-cif>"),
            )
        )
        cases = parent.get("cases")
        if not isinstance(cases, Mapping):
            paths.append(
                (
                    f"daughter:{parent_name}:<missing-cases>",
                    workspace / "<missing-cases>",
                )
            )
            continue
        for case_id, case in sorted(cases.items(), key=lambda item: str(item[0])):
            daughter = case.get("daughter_cif") if isinstance(case, Mapping) else None
            paths.append(
                (
                    f"daughter:{parent_name}:{case_id}",
                    workspace / str(daughter or "<missing-daughter-cif>"),
                )
            )
    return paths


def _path_label(path: Path, workspace: Path) -> str:
    try:
        return path.resolve().relative_to(workspace.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def _file_record(label: str, path: Path, workspace: Path) -> dict[str, Any]:
    resolved = path.resolve()
    exists = resolved.is_file()
    return {
        "label": label,
        "path": _path_label(resolved, workspace),
        "exists": exists,
        "bytes": resolved.stat().st_size if exists else None,
        "sha256": sha256_file(resolved) if exists else None,
    }


def _production_python_inputs(package_root: Path) -> list[tuple[str, Path]]:
    inputs: list[tuple[str, Path]] = []
    for source_root in ("backend", "features"):
        root = package_root / source_root
        inputs.extend(
            (f"production:{path.relative_to(package_root).as_posix()}", path)
            for path in sorted(root.rglob("*.py"), key=lambda item: item.as_posix().casefold())
        )
    return inputs


def build_method4_provenance_state(
    manifest_path: Path,
    selection: Mapping[str, Any] | None,
    *,
    workspace: Path | None = None,
    package_root: Path | None = None,
) -> dict[str, Any]:
    """Fingerprint the complete immutable surface of one local Method 4 run."""

    workspace = (workspace or WORKSPACE).resolve()
    package_root = (package_root or workspace / "ISODISTORT").resolve()
    manifest_path = manifest_path.resolve()
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        manifest = {}

    normalized_selection = {
        "parent": (selection or {}).get("parent"),
        "context": (selection or {}).get("context"),
        "case_id": (selection or {}).get("case_id"),
    }
    inputs: list[tuple[str, Path]] = [
        ("validator", package_root / "tests" / "manual" / "validate_method4_local.py"),
        ("provenance-helper", package_root / "tests" / "manual" / "method4_provenance.py"),
        ("input-preparer", package_root / "tests" / "manual" / "prepare_method4_inputs.py"),
        ("input-manifest", manifest_path),
        ("runtime-config", package_root / "resources" / "config" / "settings.yaml"),
        ("declared-project", package_root / "pyproject.toml"),
        ("declared-requirements", package_root / "requirements.txt"),
    ]
    inputs.extend(_production_python_inputs(package_root))
    inputs.extend(_backend_dependencies(package_root))
    inputs.extend(_input_paths(manifest, workspace=workspace))

    records = [
        _file_record(label, path, workspace)
        for label, path in sorted(inputs, key=lambda item: item[0].casefold())
    ]
    runtime = _runtime_state()
    material = {
        "schema": PROVENANCE_SCHEMA,
        "kind": PROVENANCE_KIND,
        "selection": normalized_selection,
        "manifest_source_signature": manifest.get("signature"),
        "runtime": runtime,
        "inputs": records,
    }
    encoded = json.dumps(
        material,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        **material,
        "signature": hashlib.sha256(encoded, usedforsecurity=False).hexdigest(),
    }


def changed_provenance_surfaces(
    before: Mapping[str, Any],
    after: Mapping[str, Any],
) -> list[str]:
    """Return stable labels for every signed surface that changed."""

    changed: list[str] = []
    for key in ("schema", "kind", "selection", "manifest_source_signature", "runtime"):
        if before.get(key) != after.get(key):
            changed.append(key)
    before_records = {
        str(row.get("label")): row
        for row in before.get("inputs", [])
        if isinstance(row, Mapping)
    }
    after_records = {
        str(row.get("label")): row
        for row in after.get("inputs", [])
        if isinstance(row, Mapping)
    }
    for label in sorted(before_records.keys() | after_records.keys(), key=str.casefold):
        if before_records.get(label) != after_records.get(label):
            changed.append(label)
    return changed


__all__ = [
    "PROVENANCE_KIND",
    "PROVENANCE_SCHEMA",
    "build_method4_provenance_state",
    "changed_provenance_surfaces",
    "sha256_file",
]
