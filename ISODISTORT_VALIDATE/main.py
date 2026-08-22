"""Interactive terminal entry point for ISODISTORT_VALIDATE."""

from __future__ import annotations

import argparse
from pathlib import Path

from batch_compare import run_batch
from compare_cif import _print_result, compare_cif


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
    return {
        "lattice_tolerance": _prompt_float("晶格容差", 1e-5),
        "coordinate_tolerance": _prompt_float("分数坐标容差", 1e-5),
        "scalar_tolerance": _prompt_float("占据率/磁矩容差", 1e-5),
        "ignore_atom_order": _prompt_yes_no("是否忽略原子排列顺序", False),
        "strict": _prompt_yes_no("是否启用严格字节比较", False),
    }


def _run_single() -> None:
    _line()
    print("单对 CIF 比较")
    local_path = _prompt("本地 ISODISTORT CIF 路径")
    reference_path = _prompt("官网/参考 CIF 路径")
    options = _prompt_common_options()
    hash_value = _prompt("参考文件 SHA-256（可留空）", allow_empty=True)
    try:
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
    print("批量 CIF 回归")
    local_dir = _prompt("本地 CIF 目录")
    reference_dir = _prompt("参考 CIF 目录")
    pattern = _prompt("文件匹配模式", "*.cif")
    options = _prompt_common_options()
    manifest_text = _prompt("SHA-256 manifest 路径（可留空）", allow_empty=True)
    json_output = _prompt_yes_no("是否输出 JSON 报告", False)
    args = argparse.Namespace(
        local_dir=local_dir,
        reference_dir=reference_dir,
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

    for result in results:
        status = "PASS" if result["passed"] else "FAIL"
        print(f"{status:4s} {result['relative_path']}")
        for issue in result["issues"]:
            print(f"      - {issue}")
    print(f"Summary: total={len(results)}, passed={len(results) - failed}, failed={failed}")


def _show_help() -> None:
    _line()
    print("ISODISTORT_VALIDATE 用于比较本地 ISODISTORT 输出与参考 CIF。")
    print("单文件模式检查晶胞、原子、坐标、占据率、磁矩和空间群。")
    print("批量模式按相对路径配对多个 CIF，并报告缺失和失败用例。")
    print("默认按语义结构判断，通过排版差异不会导致失败；严格模式要求字节一致。")


def main() -> int:
    while True:
        _line()
        print("ISODISTORT_VALIDATE 验证入口")
        print("  1. 比较一对 CIF")
        print("  2. 批量回归验证")
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


if __name__ == "__main__":
    raise SystemExit(main())
