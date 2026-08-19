"""Run compare_cif.py over a directory of local/reference CIF pairs."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

from compare_cif import compare_cif


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_manifest(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in payload.items()):
        raise ValueError("hash manifest must be a JSON object mapping relative CIF paths to SHA-256 strings")
    return {key.replace("\\", "/"): value.lower() for key, value in payload.items()}


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Batch-compare matching CIF files in two directories.")
    parser.add_argument("local_dir", help="directory containing local Isodistort_back CIF files")
    parser.add_argument("reference_dir", help="directory containing official/reference CIF files")
    parser.add_argument("--pattern", default="*.cif", help="recursive file pattern (default: *.cif)")
    parser.add_argument("--lattice-tol", type=float, default=1e-5)
    parser.add_argument("--coord-tol", type=float, default=1e-5)
    parser.add_argument("--scalar-tol", type=float, default=1e-5)
    parser.add_argument("--ignore-atom-order", action="store_true")
    parser.add_argument("--hash-manifest", type=Path, help="JSON mapping reference-relative paths to trusted SHA-256")
    parser.add_argument("--strict", action="store_true", help="also require byte-identical files")
    parser.add_argument("--json", action="store_true", help="print a machine-readable summary")
    return parser


def run_batch(args: argparse.Namespace) -> tuple[list[dict[str, Any]], int]:
    local_dir = Path(args.local_dir).expanduser().resolve()
    reference_dir = Path(args.reference_dir).expanduser().resolve()
    if not local_dir.is_dir() or not reference_dir.is_dir():
        raise FileNotFoundError("both local_dir and reference_dir must be existing directories")
    manifest = _load_manifest(args.hash_manifest)
    local_files = sorted(path for path in local_dir.rglob(args.pattern) if path.is_file())
    reference_files = sorted(path for path in reference_dir.rglob(args.pattern) if path.is_file())
    local_relative = {path.relative_to(local_dir).as_posix(): path for path in local_files}
    reference_relative = {path.relative_to(reference_dir).as_posix(): path for path in reference_files}
    results: list[dict[str, Any]] = []
    failed = 0

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
            passed = comparison.byte_exact and comparison.structure_equal if args.strict else comparison.passed
            result = {
                "relative_path": relative_path,
                "passed": passed,
                "byte_exact": comparison.byte_exact,
                "structure_equal": comparison.structure_equal,
                "reference_hash_valid": comparison.reference_hash_valid,
                "issues": comparison.details["issues"],
                "local_sha256": comparison.local_sha256,
                "reference_sha256": comparison.reference_sha256,
            }
            failed += not passed
        results.append(result)
    return results, failed


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        results, failed = run_batch(args)
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2
    if args.json:
        summary = {"total": len(results), "passed": len(results) - failed,
                   "failed": failed, "results": results}
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        for result in results:
            status = "PASS" if result["passed"] else "FAIL"
            print(f"{status:4s} {result['relative_path']}")
            for issue in result["issues"]:
                print(f"      - {issue}")
        print(f"Summary: total={len(results)}, passed={len(results) - failed}, failed={failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
