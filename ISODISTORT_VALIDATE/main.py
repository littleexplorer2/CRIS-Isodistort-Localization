"""Interactive and CLI entry point for ISODISTORT_VALIDATE."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from isodistort_validate import compare_paths as cpaths
from isodistort_validate.batch_compare import main as batch_cli_main
from isodistort_validate.batch_compare import run_batch
from isodistort_validate.compare_cif import _print_result, compare_cif
from isodistort_validate.compare_cif import main as compare_cli_main

CLI_USAGE = """无参数时进入交互菜单。命令行：
  python main.py compare [相对路径] [选项]
  python main.py batch [选项]
"""


def _line() -> None:
    print("-" * 72)


def _prompt(text: str, default: str | None = None, allow_empty: bool = False) -> str:
    hint = f" [{default}]" if default is not None else ""
    while True:
        value = input(f"{text}{hint}: ").strip()
        if value:
            return value
        if default is not None:
            return default
        if allow_empty:
            return ""
        print("输入不能为空。")


def _prompt_float(text: str, default: float) -> float:
    while True:
        try:
            value = float(_prompt(text, str(default)))
            if value < 0:
                raise ValueError
            return value
        except ValueError:
            print("请输入大于等于 0 的数字。")


def _prompt_yes_no(text: str, default: bool = False) -> bool:
    default_text = "y" if default else "n"
    while True:
        value = _prompt(f"{text} (y/n)", default_text).lower()
        if value in {"y", "yes", "是"}:
            return True
        if value in {"n", "no", "否"}:
            return False
        print("请输入 y 或 n。")


def _prompt_common_options() -> dict[str, object]:
    options: dict[str, object] = {
        "lattice_tolerance": _prompt_float("晶格容差", 1e-5),
        "coordinate_tolerance": _prompt_float("分数坐标容差", 1e-5),
        "scalar_tolerance": _prompt_float("占据率/磁矩容差", 1e-5),
    }
    print(cpaths.IGNORE_ATOM_ORDER_HELP)
    options["ignore_atom_order"] = _prompt_yes_no("是否忽略原子排列顺序", False)
    options["strict"] = _prompt_yes_no("是否启用严格字节比较", False)
    return options


def _print_compare_dirs() -> None:
    cpaths.ensure_compare_dirs()
    print(f"待验证 CIF: {cpaths.ITEM_DIR}")
    print(f"官网参考 CIF: {cpaths.TRUE_DIR}")


def _print_batch_hint() -> None:
    print(cpaths.BATCH_PAIRING_HINT)


def _choose_relative_path() -> str | None:
    cpaths.ensure_compare_dirs()
    paired, item_only, true_only = cpaths.pairing_status()
    warning = cpaths.format_unpaired_warning(item_only, true_only)
    if warning:
        print(warning)
    if not paired and not item_only and not true_only:
        print("compare/item 与 compare/true 中没有 CIF。")
        print("请把官网导出的参考 CIF 放入 compare/true/，把本地 CIF 放入 compare/item/。")
        print("两侧文件必须使用相同的相对路径（含子目录）。")
        return None
    if not paired:
        print("目前没有可以比较的成对文件。请先按上面的提示修改文件名。")
        return None
    if len(paired) == 1:
        print(f"将比较: {paired[0]}")
        return paired[0]
    print("已配对、可比较的 CIF 相对路径:")
    for index, name in enumerate(paired, start=1):
        print(f"  {index}. {name}  [item/true]")
    while True:
        raw = _prompt("选择序号（确认要比较的这一对 CIF）", "1")
        try:
            index = int(raw)
        except ValueError:
            print("请输入数字序号。")
            continue
        if 1 <= index <= len(paired):
            return paired[index - 1]
        print("序号超出范围。")


def _run_single() -> None:
    _line()
    print("单对 CIF 比较（固定目录 compare/item ↔ compare/true）")
    _print_compare_dirs()
    relative = _choose_relative_path()
    if not relative:
        return
    options = _prompt_common_options()
    hash_value = _prompt("参考文件 SHA-256（可留空）", allow_empty=True)
    try:
        local_path, reference_path = cpaths.resolve_pair(relative)
        result = compare_cif(
            local_path,
            reference_path,
            options["lattice_tolerance"],
            options["coordinate_tolerance"],
            options["scalar_tolerance"],
            options["ignore_atom_order"],
            hash_value or None,
        )
        _print_result(result, options["strict"])
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}")


def _run_batch() -> None:
    _line()
    print("批量 CIF 回归（固定目录 compare/item ↔ compare/true）")
    _print_compare_dirs()
    paired, item_only, true_only = cpaths.pairing_status()
    warning = cpaths.format_unpaired_warning(item_only, true_only)
    if warning:
        print(warning)
    if not paired:
        if item_only or true_only:
            print("目前没有可以比较的成对文件。请先按上面的提示修改文件名。")
        else:
            print("没有可配对的 CIF。请把文件放入 compare/item 与 compare/true。")
            _print_batch_hint()
        return
    _print_batch_hint()
    pattern = _prompt("文件匹配模式", "*.cif")
    options = _prompt_common_options()
    manifest_text = _prompt("SHA-256 manifest 路径（可留空）", allow_empty=True)
    json_output = _prompt_yes_no("是否输出 JSON 报告", False)
    args = argparse.Namespace(
        local_dir=str(cpaths.ITEM_DIR),
        reference_dir=str(cpaths.TRUE_DIR),
        pattern=pattern,
        lattice_tol=options["lattice_tolerance"],
        coord_tol=options["coordinate_tolerance"],
        scalar_tol=options["scalar_tolerance"],
        ignore_atom_order=options["ignore_atom_order"],
        hash_manifest=Path(manifest_text) if manifest_text else None,
        strict=options["strict"],
        json=json_output,
    )
    try:
        results, failed = run_batch(args)
    except (OSError, ValueError, TypeError) as exc:
        print(f"错误: {exc}")
        return

    if json_output:
        import json

        summary = {"total": len(results), "passed": len(results) - failed,
                   "failed": failed, "results": results}
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return

    if not results:
        print("没有可配对的 CIF。请把文件放入 compare/item 与 compare/true。")
        _print_batch_hint()
        return
    for result in results:
        status = "PASS" if result["passed"] else "FAIL"
        print(f"{status:4s} {result['relative_path']}")
        for issue in result["issues"]:
            print(f"      - {issue}")
    print(f"Summary: total={len(results)}, passed={len(results) - failed}, failed={failed}")


def _show_help() -> None:
    _line()
    print("ISODISTORT_VALIDATE 用于比较本地 ISODISTORT 输出与官网参考 CIF。")
    print("比较目录固定为（不要改用其它路径）：")
    print(f"  待验证: {cpaths.ITEM_DIR}")
    print(f"  标准答案: {cpaths.TRUE_DIR}")
    print("两侧按相对路径配对，例如：")
    print("  compare/item/LD1 C1/subgroup.cif")
    print("  compare/true/LD1 C1/subgroup.cif")
    _print_batch_hint()
    print("单文件模式检查晶胞、原子、坐标、占据率、磁矩和空间群。")
    print("批量模式比较两个文件夹中的全部配对，并报告缺失和失败用例。")
    print("默认 n：按 CIF 原子行号一一对应。选 y：按同种元素和分数坐标配对，可忽略行顺序，")
    print("但不忽略坐标或元素种类的真实差异。")
    print()
    print(CLI_USAGE)


def _interactive_menu() -> int:
    while True:
        _line()
        print("ISODISTORT_VALIDATE 验证入口")
        cpaths.warn_unpaired_filenames(file=sys.stdout)
        print("  1. 比较一对 CIF")
        print("  2. 批量回归验证")
        print("     （true/ 中官网文件须改名，与 item/ 一一对应）")
        print("  3. 查看验证说明")
        print("  0. 退出")
        choice = _prompt("请选择", "1")
        if choice == "1":
            _run_single()
        elif choice == "2":
            _run_batch()
        elif choice == "3":
            _show_help()
        elif choice == "0":
            print("已退出。")
            return 0
        else:
            print("无效选择，请输入 0、1、2 或 3。")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    cpaths.ensure_compare_dirs()
    if not argv:
        return _interactive_menu()
    if argv[0] in {"-h", "--help", "help"}:
        print("比较 compare/item 与 compare/true 中的 CIF。")
        print(CLI_USAGE)
        _print_batch_hint()
        print("查看子命令帮助：")
        print("  python main.py compare --help")
        print("  python main.py batch --help")
        return 0
    if argv[0] == "compare":
        return compare_cli_main(argv[1:])
    if argv[0] == "batch":
        return batch_cli_main(argv[1:])
    print(f"未知命令: {argv[0]}", file=sys.stderr)
    print(CLI_USAGE, file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
