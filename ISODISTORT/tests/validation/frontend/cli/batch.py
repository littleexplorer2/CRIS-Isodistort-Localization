"""Command-line parsing and rendering for batch comparisons."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ...application import workspace as cpaths
from ...application.batch import run_batch
from ...application.workspace import get_config


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m tests.validation batch",
        description=(
            "Batch-compare CIF files in compare/item against matching files in compare/true (same relative paths)."
        ),
        epilog=cpaths.BATCH_PAIRING_HINT,
    )
    cfg = get_config()
    parser.add_argument("--local-dir", type=Path, help="frozen local CIF directory")
    parser.add_argument("--reference-dir", type=Path, help="read-only reference CIF directory")
    parser.add_argument(
        "--pattern",
        default=cfg.pattern,
        help=f"recursive file pattern (default: {cfg.pattern})",
    )
    parser.add_argument("--lattice-tol", type=float, default=cfg.lattice_tolerance)
    parser.add_argument("--coord-tol", type=float, default=cfg.coordinate_tolerance)
    parser.add_argument("--scalar-tol", type=float, default=cfg.scalar_tolerance)
    parser.add_argument(
        "--ignore-atom-order",
        action=argparse.BooleanOptionalAction,
        default=cfg.ignore_atom_order,
        help=(
            "ignore CIF atom row order and match by element plus fractional coordinates; "
            "default compares row-by-row. Does not ignore real coordinate or species differences"
        ),
    )
    parser.add_argument("--hash-manifest", type=Path, help="JSON mapping reference-relative paths to trusted SHA-256")
    parser.add_argument("--strict", action=argparse.BooleanOptionalAction, default=cfg.strict,
                        help="require semantic equality, trusted reference hashes, and byte-identical files")
    parser.add_argument("--json", action="store_true", help="print a machine-readable summary")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    print(cpaths.BATCH_PAIRING_HINT, file=sys.stderr)
    if bool(args.local_dir) != bool(args.reference_dir):
        parser.error("use --local-dir and --reference-dir together")
    if args.local_dir is None:
        cpaths.warn_unpaired_filenames()
        args.local_dir = cpaths.ITEM_DIR
        args.reference_dir = cpaths.TRUE_DIR
    try:
        results, failed = run_batch(args)
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2
    if args.json:
        summary = {"total": len(results), "passed": len(results) - failed, "failed": failed, "results": results}
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        for result in results:
            status = result.get("verification", {}).get("status", "pass" if result["passed"] else "fail").upper()
            print(f"{status:4s} {result['relative_path']}")
            for issue in result["issues"]:
                print(f"      - {issue}")
        print(f"Summary: total={len(results)}, passed={len(results) - failed}, failed={failed}")
    return 0 if failed == 0 else 1
