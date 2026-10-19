"""Method 4 input/selection adapter for the shared validation provenance service.

The common service owns file capture, code/config/runtime/data fingerprinting,
and drift detection. This module only defines Method 4's independent inputs.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from tests.validation.backend.provenance import (
    FileSnapshot,
    build_provenance_state,
    changed_provenance_surfaces,
    sha256_file,
)

WORKSPACE = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
PROVENANCE_SCHEMA = 2
PROVENANCE_KIND = "method4-local-live-validation-v2"


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


def build_method4_provenance_state(
    manifest_path: Path | FileSnapshot,
    selection: Mapping[str, Any] | None,
    *,
    workspace: Path | None = None,
    package_root: Path | None = None,
) -> dict[str, Any]:
    """Describe Method 4 inputs and sign them through the common service."""
    workspace = (workspace or WORKSPACE).resolve()
    package_root = (package_root or workspace / "ISODISTORT").resolve()
    if isinstance(manifest_path, FileSnapshot):
        snapshot = manifest_path
        manifest_path = snapshot.path
    else:
        manifest_path = manifest_path.resolve()
        snapshot = FileSnapshot.capture(manifest_path) if manifest_path.is_file() else None
    try:
        manifest = json.loads(snapshot.text) if snapshot else {}
    except (UnicodeError, json.JSONDecodeError):
        manifest = {}
    if not isinstance(manifest, Mapping):
        manifest = {}
    inputs = [
        ("validator", package_root / "tests"  / "validate_method4_local.py"),
        ("provenance-helper", package_root / "tests"  / "method4_provenance.py"),
        ("input-preparer", package_root / "tests"  / "prepare_method4_inputs.py"),
        ("input-manifest", manifest_path),
        *_input_paths(manifest, workspace=workspace),
    ]
    return build_provenance_state(
        schema=PROVENANCE_SCHEMA,
        kind=PROVENANCE_KIND,
        workspace=workspace,
        package_root=package_root,
        inputs=inputs,
        metadata={
            "selection": {key: (selection or {}).get(key) for key in ("parent", "context", "case_id")},
            "manifest_source_signature": manifest.get("signature"),
        },
        captured_inputs=[snapshot] if snapshot else (),
    )


__all__ = [
    "PROVENANCE_KIND", "PROVENANCE_SCHEMA", "build_method4_provenance_state",
    "changed_provenance_surfaces", "sha256_file",
]
