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
        prog="main.py batch",
        description=(
            "Batch-compare CIF files in compare/item against matching files in compare/true (same relative paths)."
        ),
        epilog=cpaths.BATCH_PAIRING_HINT,
    )
    cfg = get_config()
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
        action="store_true",
        help=(
            "ignore CIF atom row order and match by element plus fractional coordinates; "
            "default compares row-by-row. Does not ignore real coordinate or species differences"
        ),
    )
    parser.add_argument("--hash-manifest", type=Path, help="JSON mapping reference-relative paths to trusted SHA-256")
    parser.add_argument("--strict", action="store_true", help="also require byte-identical files")
    parser.add_argument("--json", action="store_true", help="print a machine-readable summary")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    print(cpaths.BATCH_PAIRING_HINT, file=sys.stderr)
    cpaths.warn_unpaired_filenames()
    cpaths.ensure_compare_dirs()
    args.local_dir = str(cpaths.ITEM_DIR)
    args.reference_dir = str(cpaths.TRUE_DIR)
    try:
        results, failed = run_batch(args)
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2
    if args.json:
        summary = {"total": len(results), "passed": len(results) - failed, "failed": failed, "results": results}
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        if not results:
            print(
                "compare/item 与 compare/true 中没有可配对的 CIF。"
                "请把官网参考 CIF 放入 compare/true/，把待验证 CIF 放入 compare/item/，"
                "并保持相对路径相同。"
            )
        for result in results:
            status = "PASS" if result["passed"] else "FAIL"
            print(f"{status:4s} {result['relative_path']}")
            for issue in result["issues"]:
                print(f"      - {issue}")
        print(f"Summary: total={len(results)}, passed={len(results) - failed}, failed={failed}")
    return 0 if failed == 0 else 1
