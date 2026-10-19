"""Fingerprint validation inputs and publish reports without executing production code.

Content hashes identify captured bytes, not scientific correctness. Every file
record derives its length and hash from the same read. Generated output/cache
trees are deliberately outside the automatically collected input surface.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

RUNTIME_DISTRIBUTIONS = ("numpy", "scipy", "pymatgen", "spglib", "PyYAML")
ISO_DATA_BASENAMES = (
    "const.dat", "data_diperiodic.txt", "data_images.txt", "data_irreps.txt",
    "data_isotropy.txt", "data_little.txt", "data_magnetic.txt", "data_space.txt",
    "data_ssg.txt", "data_ssgmag.txt", "data_wyckoff.txt",
)


@dataclass(frozen=True, slots=True)
class FileSnapshot:
    """One immutable content read shared by parsers and fingerprint consumers."""

    path: Path
    content: bytes = field(repr=False)

    @classmethod
    def capture(cls, path: str | Path) -> FileSnapshot:
        resolved = Path(path).expanduser().resolve()
        return cls(resolved, resolved.read_bytes())

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.content, usedforsecurity=False).hexdigest()

    @property
    def text(self) -> str:
        # Match Python text-mode newline decoding while retaining original bytes.
        return self.content.decode("utf-8-sig").replace("\r\n", "\n").replace("\r", "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validation_source_paths(package_root: Path) -> dict[str, Path]:
    """Canonical registry of executable comparison and evidence implementation."""
    root = package_root / "tests" / "validation"
    return {
        f"cif-validator:{path.relative_to(root).as_posix()}": path
        for layer in ("backend", "application")
        for path in sorted((root / layer).rglob("*.py"))
    }


def production_source_paths(package_root: Path) -> dict[str, Path]:
    """All executable production computation, session and workflow sources."""
    return {
        f"production:{path.relative_to(package_root).as_posix()}": path
        for layer in ("backend", "application")
        for path in sorted((package_root / layer).rglob("*.py"))
    }


def update_file_digest(digest: Any, path: Path, label: str) -> None:
    """Append a labelled file or explicit missing marker to an audit digest."""
    digest.update(label.encode("utf-8"))
    digest.update(b"\0")
    if path.is_file():
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
    else:
        digest.update(b"<missing>")
    digest.update(b"\0")


def write_json_atomic(path: Path, payload: Any, **json_options: Any) -> None:
    """Publish a complete audit JSON; failed serialization/replacement preserves the old file."""
    options = {"ensure_ascii": False, "indent": 2, **json_options}
    content = json.dumps(payload, **options) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", prefix=path.name + ".", suffix=".tmp",
            dir=path.parent, delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def runtime_state(names: Sequence[str] = RUNTIME_DISTRIBUTIONS) -> dict[str, Any]:
    distributions: dict[str, str] = {}
    for name in names:
        try:
            distributions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = "<missing>"
        except Exception as exc:  # noqa: BLE001 - broken metadata remains signable
            distributions[name] = f"<unavailable:{type(exc).__module__}.{type(exc).__qualname__}>"
    executable = Path(sys.executable).resolve()
    snapshot = FileSnapshot.capture(executable) if executable.is_file() else None
    return {
        "python": {
            "implementation": sys.implementation.name,
            "version": list(sys.version_info[:5]),
            "cache_tag": sys.implementation.cache_tag,
            "build": sys.version,
            "executable": {
                "path": str(executable), "exists": snapshot is not None,
                "sha256": snapshot.sha256 if snapshot else None,
                "bytes": len(snapshot.content) if snapshot else None,
            },
            "prefix": sys.prefix,
            "base_prefix": sys.base_prefix,
        },
        "platform": {
            "system": platform.system(), "release": platform.release(), "machine": platform.machine(),
        },
        "distributions": distributions,
    }


def _backend_dependencies(config_path: Path, config: FileSnapshot | None) -> list[tuple[str, Path]]:
    isobyu: dict[str, Any] = {}
    try:
        loaded = yaml.safe_load(config.text) if config else None
        if isinstance(loaded, dict) and isinstance(loaded.get("isobyu"), dict):
            isobyu = loaded["isobyu"]
    except (UnicodeError, yaml.YAMLError):
        pass  # The invalid configuration itself stays in the signed records.

    def configured_path(key: str, fallback: str) -> Path:
        value = isobyu.get(key)
        path = Path(value if isinstance(value, str) and value.strip() else fallback)
        return (path if path.is_absolute() else config_path.parent / path).resolve()

    bin_dir = configured_path("bin_dir", "../isobyu")
    data_dir = configured_path("data_dir", "../isobyu")
    dependencies = []
    for label in ("iso", "findsym", "smodes"):
        value = isobyu.get(f"{label}_bin")
        name = value if isinstance(value, str) and value.strip() else label
        dependencies.append((f"backend:{label}-binary", bin_dir / name))
    names = set(ISO_DATA_BASENAMES)
    names.update(path.name for path in data_dir.glob("data_*.txt"))
    dependencies.extend((f"backend:iso-data:{name}", data_dir / name) for name in sorted(names, key=str.casefold))
    return dependencies


def build_provenance_state(
    *,
    schema: int,
    kind: str,
    workspace: Path,
    package_root: Path,
    inputs: Sequence[tuple[str, Path]],
    metadata: Mapping[str, Any],
    captured_inputs: Sequence[FileSnapshot] = (),
) -> dict[str, Any]:
    """Sign common runtime/code/data plus an audit's independently supplied inputs.

    Missing files remain explicit records. Callers decide whether their task may
    execute with those records; a signature never constitutes an acceptance gate.
    Captured manifest/config inputs can be passed to keep parsing and signing on
    exactly the same bytes. No output directories are read automatically.
    """
    if {"schema", "kind", "inputs", "runtime", "signature"} & metadata.keys():
        raise ValueError("provenance metadata shadows reserved fields")
    workspace, package_root = workspace.resolve(), package_root.resolve()
    captured = {item.path.resolve(): item for item in captured_inputs}

    def capture(path: Path) -> FileSnapshot | None:
        path = path.resolve()
        if path not in captured:
            captured[path] = FileSnapshot.capture(path) if path.is_file() else None
        return captured[path]

    config_path = package_root / "resources" / "config" / "settings.yaml"
    collected = [
        *inputs,
        ("runtime-config", config_path),
        ("declared-project", package_root / "pyproject.toml"),
        ("declared-requirements", package_root / "requirements.txt"),
    ]
    collected.extend(production_source_paths(package_root).items())
    # The implementation of the signer is part of its own provenance surface.
    collected.extend(
        (f"provenance:{path.name}", path)
        for path in sorted(Path(__file__).parent.glob("*.py"))
    )
    collected.extend(_backend_dependencies(config_path, capture(config_path)))
    labels = [label for label, _ in collected]
    if len(labels) != len(set(labels)):
        raise ValueError("provenance input labels must be unique")
    records = []
    for label, path in sorted(collected, key=lambda item: item[0].casefold()):
        resolved = path.resolve()
        snapshot = capture(resolved)
        try:
            display = resolved.relative_to(workspace).as_posix()
        except ValueError:
            display = str(resolved)
        records.append({
            "label": label, "path": display, "exists": snapshot is not None,
            "bytes": len(snapshot.content) if snapshot else None,
            "sha256": snapshot.sha256 if snapshot else None,
        })
    material = {"schema": schema, "kind": kind, **metadata, "runtime": runtime_state(), "inputs": records}
    encoded = json.dumps(material, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return {**json.loads(encoded), "signature": hashlib.sha256(encoded, usedforsecurity=False).hexdigest()}


def changed_provenance_surfaces(before: Mapping[str, Any], after: Mapping[str, Any]) -> list[str]:
    """Identify changed metadata and file labels, ignoring the derived signature."""
    changed = [
        key for key in sorted((before.keys() | after.keys()) - {"inputs", "signature"})
        if before.get(key) != after.get(key)
    ]
    old = {str(row.get("label")): row for row in before.get("inputs", []) if isinstance(row, Mapping)}
    new = {str(row.get("label")): row for row in after.get("inputs", []) if isinstance(row, Mapping)}
    changed.extend(label for label in sorted(old.keys() | new.keys(), key=str.casefold) if old.get(label) != new.get(label))
    return changed
