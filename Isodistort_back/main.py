"""交互式主入口：引导用户完成一轮完整畸变计算。"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import List

from isocore.api import IsoDistort


def _prompt(text: str, default: str | None = None) -> str:
    """Relative path: main.py"""

    hint = f" [{default}]" if default is not None else ""
    while True:
        val = input(f"{text}{hint}: ").strip()
        if val:
            return val
        if default is not None:
            return default


def _prompt_float(text: str, default: float) -> float:
    """Relative path: main.py"""

    while True:
        raw = _prompt(text, str(default))
        try:
            return float(raw)
        except ValueError:
            print("输入必须是数字，请重试。")


def _prompt_int(text: str, default: int) -> int:
    """Relative path: main.py"""

    while True:
        raw = _prompt(text, str(default))
        try:
            return int(raw)
        except ValueError:
            print("输入必须是整数，请重试。")


def _prompt_supercell(default_cell: List[int]) -> List[int]:
    """Relative path: main.py"""

    default_str = ",".join(str(x) for x in default_cell)
    while True:
        raw = _prompt("请输入超胞参数 a,b,c", default_str)
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        if len(parts) != 3:
            print("需要输入 3 个整数，例如 1,1,1")
            continue
        try:
            cell = [int(x) for x in parts]
        except ValueError:
            print("超胞参数必须是整数。")
            continue
        if any(x < 1 for x in cell):
            print("超胞参数必须 >= 1。")
            continue
        return cell


def _find_cif_candidates(root: Path, limit: int = 20) -> List[Path]:
    """Relative path: main.py"""

    candidates: List[Path] = []
    for p in root.rglob("*.cif"):
        if p.is_file():
            candidates.append(p)
            if len(candidates) >= limit:
                break
    return sorted(candidates)


def _choose_cif(project_root: Path) -> str:
    """Relative path: main.py"""

    print("\n=== 步骤 1：选择 CIF 文件 ===")
    candidates = _find_cif_candidates(project_root)

    if candidates:
        print("检测到以下 CIF 文件：")
        for idx, path in enumerate(candidates, start=1):
            try:
                rel = path.relative_to(project_root)
                show = str(rel)
            except ValueError:
                show = str(path)
            print(f"  {idx}. {show}")
        print("  0. 手动输入路径")

        while True:
            choice = _prompt_int("请选择 CIF 编号", 1)
            if choice == 0:
                break
            if 1 <= choice <= len(candidates):
                return str(candidates[choice - 1])
            print("编号超出范围，请重试。")

    manual = _prompt("请输入 CIF 文件路径（可相对当前目录）")
    return manual


def _choose_distortion_type() -> str:
    """Relative path: main.py"""

    print("\n=== 步骤 2：选择畸变类型 ===")
    options = {
        1: "displacement",
        2: "order",
        3: "strain",
        4: "magnetic",
    }
    print("1. displacement（原子位移）")
    print("2. order（原子有序）")
    print("3. strain（晶格应变）")
    print("4. magnetic（磁矩）")
    choice = _prompt_int("请选择类型编号", 1)
    return options.get(choice, "displacement")


def _choose_subgroup_idx(subgroups) -> int:
    """Relative path: main.py"""

    print("\n=== 步骤 3：选择子群路径 ===")
    if not subgroups:
        raise RuntimeError("未找到可选子群。")

    print("可选子群（最多显示前 20 个）：")
    for sg in subgroups[:20]:
        print(f"  idx={sg.index:3d} | SG #{sg.space_group_number:<3d} | {sg.space_group_symbol}")
    if len(subgroups) > 20:
        print(f"  ... 其余 {len(subgroups) - 20} 个请按 idx 输入")

    default_idx = subgroups[0].index
    return _prompt_int("请输入要选择的 subgroup_idx", default_idx)


def _choose_irrep_label(iso: IsoDistort) -> str | None:
    """Relative path: main.py"""

    modes = iso.distortion_modes
    if not modes:
        return None

    print("\n=== 步骤 4：选择模式（irrep） ===")
    for i, mode in enumerate(modes, start=1):
        print(f"  {i}. {mode.irrep_label} ({mode.dimension}D)")
    print("直接回车则使用第一个模式。")

    raw = _prompt("请输入 irrep 标签（如 GM1+）或序号", "")
    if not raw:
        return None
    if raw.isdigit():
        idx = int(raw)
        if 1 <= idx <= len(modes):
            return modes[idx - 1].irrep_label
    return raw


def _choose_export_formats() -> List[str]:
    """Relative path: main.py"""

    print("\n=== 步骤 6：导出格式 ===")
    print("可选格式：cif, poscar")
    raw = _prompt("请输入导出格式（逗号分隔）", "cif")
    formats = [x.strip().lower() for x in raw.split(",") if x.strip()]
    return formats or ["cif"]


def run_wizard() -> int:
    """Relative path: main.py"""

    print("isodistort - Interactive Workflow")
    print("将按向导完成：加载 CIF -> 枚举子群 -> 选择路径 -> 生成畸变 -> 导出")

    project_root = Path(__file__).resolve().parent
    iso = IsoDistort()

    # 1) CIF
    cif_path = _choose_cif(project_root)
    iso.load_structure(cif_path)

    # 2) 畸变类型 + 子群
    distortion_type = _choose_distortion_type()
    subgroups = iso.list_subgroups(distortion_type=distortion_type)
    subgroup_idx = _choose_subgroup_idx(subgroups)
    iso.select_path(subgroup_idx=subgroup_idx, distortion_type=distortion_type)

    # 3) 模式 + 幅度 + 超胞
    irrep = _choose_irrep_label(iso)
    default_amp = float(iso.cfg.defaults.get("default_amplitude", 1.0))
    default_cell = list(iso.cfg.defaults.get("default_supercell", [1, 1, 1]))
    amplitude = _prompt_float("请输入幅度 amplitude", default_amp)
    supercell = _prompt_supercell(default_cell)

    print("\n=== 步骤 5：生成畸变结构 ===")
    distorted = iso.generate_distortion(
        irrep_label=irrep,
        amplitude=amplitude,
        supercell=supercell,
    )
    print(f"畸变后结构原子数: {len(distorted)}")

    # 4) 导出
    formats = _choose_export_formats()
    name = _prompt("请输入导出文件名前缀", "distorted_output")
    out_paths = iso.export(name, formats=formats)

    print("\n=== 完成 ===")
    print("导出文件：")
    for p in out_paths:
        print(f"  {p}")

    return 0


def main() -> int:
    """Relative path: main.py"""
    
    try:
        return run_wizard()
    except KeyboardInterrupt:
        print("\n用户中断。")
        return 130
    except Exception as exc:
        print(f"\n运行失败: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())