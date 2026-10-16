"""Command-line parsing and rendering for one CIF pair."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict

from ...application import workspace as cpaths
from ...application.single import ComparisonResult, compare_cif
from ...application.workspace import default_relative_path, get_config


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m tests.validation compare",
        description=(
            "Compare a local ISODISTORT CIF in compare/item with the matching "
            "official CIF in compare/true (same relative path)."
        ),
    )
    parser.add_argument(
        "relative_path",
        nargs="?",
        help="CIF path relative to compare/item and compare/true; omit if only one pair exists",
    )
    cfg = get_config()
    parser.add_argument("--local", help="explicit path to a frozen local CIF")
    parser.add_argument("--reference", help="explicit path to the read-only reference CIF")
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
    parser.add_argument("--reference-sha256", help="expected SHA-256 for the trusted reference CIF")
    parser.add_argument(
        "--strict",
        action=argparse.BooleanOptionalAction,
        default=cfg.strict,
        help="require byte-identical files in addition to semantic equality",
    )
    parser.add_argument(
        "--structure-only",
        action="store_true",
        help="deprecated alias; semantic equality is the default",
    )
    parser.add_argument("--json", action="store_true", help="print the complete result as JSON")
    return parser


def _print_result(result: ComparisonResult, strict: bool) -> None:
    print("=" * 72)
    print("CIF 比对结果")
    print(f"本地文件: {result.local_path}")
    print(f"官网文件: {result.reference_path}")
    print(f"字节完全一致: {'是' if result.byte_exact else '否'}")
    print(f"文本内容一致: {'是' if result.text_exact else '否'}")
    print(f"解析后结构与元数据一致: {'是' if result.structure_equal else '否'}")
    print(f"本地 SHA-256: {result.local_sha256}")
    print(f"官网 SHA-256: {result.reference_sha256}")
    hash_status = "未提供期望哈希" if result.reference_hash_valid is None else (
        "通过" if result.reference_hash_valid else "失败"
    )
    print(f"参考文件哈希校验: {hash_status}")
    details = result.details
    print("差异摘要:")
    print(f"  原子数: {details['local_atom_count']} / {details['reference_atom_count']}")
    print(f"  晶格最大绝对误差: {details['lattice_max_abs_difference']}")
    print(f"  分数坐标最大周期误差: {details['coordinate_max_periodic_difference']}")
    print(f"  占据率一致: {details['occupancies_equal']}")
    print(f"  磁矩一致: {details['magnetic_moments_equal']}")
    print(f"  空间群声明一致: {details['declared_space_group_equal']}")
    print(f"  spglib 推断一致: {details['inferred_space_group_equal']}")
    if details["issues"]:
        print("差异原因:")
        for issue in details["issues"]:
            print(f"  - {issue}")
    passed = result.passed and (not strict or result.byte_exact)
    print("结论: " + ("PASS" if passed else "FAIL"))
    if result.structure_equal and not result.byte_exact:
        print("提示: 结构和关键元数据一致，但 CIF 排版或非关键文本不同。")


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if bool(args.local) != bool(args.reference) or (args.local and args.relative_path):
        parser.error("use --local and --reference together, without relative_path")
    try:
        if args.local:
            local_cif, reference_cif = args.local, args.reference
        else:
            cpaths.warn_unpaired_filenames()
            relative = args.relative_path or default_relative_path()
            local_cif, reference_cif = cpaths.resolve_pair(relative)
        result = compare_cif(
            local_cif,
            reference_cif,
            args.lattice_tol,
            args.coord_tol,
            args.scalar_tol,
            args.ignore_atom_order,
            args.reference_sha256,
        )
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2

    passed = result.passed and (not args.strict or result.byte_exact)
    if args.json:
        payload = asdict(result)
        payload["passed"] = passed
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        _print_result(result, args.strict)
    return 0 if passed else 1
