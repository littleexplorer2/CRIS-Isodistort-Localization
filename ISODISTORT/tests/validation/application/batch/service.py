"""Pair and compare complete directory trees without CLI parsing or output."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Protocol

from ...backend import paths as cpaths
from ..single import VerificationPolicy, compare_cif


class BatchRequest(Protocol):
    """Input contract also accepted from legacy argparse.Namespace callers."""

    pattern: str
    lattice_tol: float
    coord_tol: float
    scalar_tol: float
    ignore_atom_order: bool
    hash_manifest: Path | None
    strict: bool


def _load_manifest(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in payload.items()):
        raise ValueError("hash manifest must be a JSON object mapping relative CIF paths to SHA-256 strings")
    return {key.replace("\\", "/"): value.lower() for key, value in payload.items()}


def run_batch(args: BatchRequest) -> tuple[list[dict[str, Any]], int]:
    local_dir = Path(getattr(args, "local_dir", None) or cpaths.ITEM_DIR).expanduser().resolve()
    reference_dir = Path(getattr(args, "reference_dir", None) or cpaths.TRUE_DIR).expanduser().resolve()
    if not local_dir.is_dir() or not reference_dir.is_dir():
        raise FileNotFoundError("both local_dir and reference_dir must be existing directories")
    manifest = _load_manifest(args.hash_manifest)
    local_files = sorted(path for path in local_dir.rglob(args.pattern) if path.is_file())
    reference_files = sorted(path for path in reference_dir.rglob(args.pattern) if path.is_file())
    if not local_files and not reference_files:
        raise ValueError("no CIF comparisons executed: the selected batch is empty")
    local_relative = {path.relative_to(local_dir).as_posix(): path for path in local_files}
    reference_relative = {path.relative_to(reference_dir).as_posix(): path for path in reference_files}
    results: list[dict[str, Any]] = []
    failed = 0
    policy = VerificationPolicy(strict=args.strict)

    for relative_path in sorted(set(local_relative) | set(reference_relative)):
        local_path = local_relative.get(relative_path)
        reference_path = reference_relative.get(relative_path)
        if local_path is None or reference_path is None:
            result = {"relative_path": relative_path, "passed": False, "issues": ["missing matching CIF"]}
            failed += 1
        else:
            expected_hash = manifest.get(relative_path)
            comparison = compare_cif(
                local_path,
                reference_path,
                args.lattice_tol,
                args.coord_tol,
                args.scalar_tol,
                args.ignore_atom_order,
                expected_hash,
            )
            decision = policy.evaluate(comparison)
            passed = decision.passed
            result = {
                "relative_path": relative_path,
                "passed": passed,
                "byte_exact": comparison.byte_exact,
                "structure_equal": comparison.structure_equal,
                "reference_hash_valid": comparison.reference_hash_valid,
                "issues": list(dict.fromkeys([*comparison.details["issues"], *decision.issues])),
                "verification": asdict(decision),
                "local_sha256": comparison.local_sha256,
                "reference_sha256": comparison.reference_sha256,
            }
            failed += not passed
        results.append(result)
    return results, failed
