"""统一交互入口：以终端菜单方式复现 ISODISTORT 网站核心交互流程。"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Dict, List, Sequence

from isocore.api import IsoDistort


DISTORTION_TYPE_MAP = {
    1: "displacement",
    2: "order",
    3: "strain",
    4: "magnetic",
}


def _line() -> None:
    print("-" * 78)


def _prompt(text: str, default: str | None = None) -> str:
    hint = f" [{default}]" if default is not None else ""
    while True:
        value = input(f"{text}{hint}: ").strip()
        if value:
            return value
        if default is not None:
            return default


def _prompt_int(text: str, default: int | None = None) -> int:
    while True:
        raw = _prompt(text, str(default) if default is not None else None)
        try:
            return int(raw)
        except ValueError:
            print("输入必须是整数，请重试。")


def _prompt_float(text: str, default: float | None = None) -> float:
    while True:
        raw = _prompt(text, str(default) if default is not None else None)
        try:
            return float(raw)
        except ValueError:
            print("输入必须是数字，请重试。")


def _prompt_yes_no(text: str, default_yes: bool = False) -> bool:
    default = "y" if default_yes else "n"
    while True:
        raw = _prompt(f"{text} (y/n)", default).lower()
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False
        print("请输入 y 或 n。")


def _find_cif_candidates(root: Path, limit: int = 30) -> List[Path]:
    found: List[Path] = []
    for path in root.rglob("*.cif"):
        if path.is_file():
            found.append(path)
            if len(found) >= limit:
                break
    return sorted(found)


def _choose_cif(project_root: Path, title: str) -> str:
    _line()
    print(title)
    candidates = _find_cif_candidates(project_root)

    if candidates:
        print("检测到 CIF 文件：")
        for idx, path in enumerate(candidates, start=1):
            try:
                display = str(path.relative_to(project_root))
            except ValueError:
                display = str(path)
            print(f"  {idx:2d}. {display}")
        print("   0. 手动输入路径")

        while True:
            choice = _prompt_int("请选择编号", 1)
            if choice == 0:
                break
            if 1 <= choice <= len(candidates):
                return str(candidates[choice - 1])
            print("编号超出范围，请重试。")

    return _prompt("请输入 CIF 文件路径")


def _prompt_distortion_types(default: Sequence[str]) -> List[str]:
    _line()
    print("Distortion Types（对应网页上的勾选框）")
    print("  1. displacement")
    print("  2. order")
    print("  3. strain")
    print("  4. magnetic")
    default_raw = ",".join(default)
    raw = _prompt(
        "请输入要启用的类型（编号或名称，逗号分隔，如 1,3 或 displacement,strain）",
        default_raw,
    )

    labels: List[str] = []
    for token in [x.strip().lower() for x in raw.split(",") if x.strip()]:
        if token.isdigit() and int(token) in DISTORTION_TYPE_MAP:
            labels.append(DISTORTION_TYPE_MAP[int(token)])
        elif token in DISTORTION_TYPE_MAP.values():
            labels.append(token)

    dedup: List[str] = []
    seen = set()
    for label in labels:
        if label not in seen:
            dedup.append(label)
            seen.add(label)

    if not dedup:
        return ["displacement", "strain"]
    return dedup


def _prompt_supercell(default_cell: Sequence[int]) -> List[int]:
    raw = _prompt("请输入超胞 (a,b,c)", ",".join(str(x) for x in default_cell))
    parts = [x.strip() for x in raw.split(",") if x.strip()]
    if len(parts) != 3:
        raise ValueError("超胞必须输入 3 个整数")
    values = [int(x) for x in parts]
    if any(v < 1 for v in values):
        raise ValueError("超胞参数必须 >= 1")
    return values


def _prompt_basis_matrix() -> List[List[str]]:
    print("请输入 3x3 基矢矩阵（每行 3 个值，可用分数如 1/2，用空格分隔）")
    rows: List[List[str]] = []
    for i in range(3):
        row = _prompt(f"第 {i + 1} 行", "1 0 0" if i == 0 else "0 1 0" if i == 1 else "0 0 1")
        parts = [x for x in row.split() if x]
        if len(parts) != 3:
            raise ValueError("每一行必须有 3 个值")
        rows.append(parts)
    return rows


class IsoDistortConsoleApp:
    def __init__(self) -> None:
        self.project_root = Path(__file__).resolve().parent
        self.iso = IsoDistort()
        self.distortion_types: List[str] = ["displacement", "strain"]
        self.last_method1 = []
        self.last_method2 = None
        self.last_method3 = []

    def run(self) -> int:
        self._banner()
        self._load_parent_structure()
        while True:
            choice = self._main_menu()
            if choice == 0:
                print("已退出。")
                return 0
            if choice == 1:
                self._load_parent_structure()
            elif choice == 2:
                self._set_distortion_types()
            elif choice == 3:
                self._run_method_1()
            elif choice == 4:
                self._run_method_2()
            elif choice == 5:
                self._run_method_3()
            elif choice == 6:
                self._run_method_4()
            elif choice == 7:
                self._distortion_page_menu()
            elif choice == 8:
                self._show_state()

    def _banner(self) -> None:
        """打印使用前说明"""
        _line()
        print("ISODISTORT Local Console")
        print("按官网 Search/Distortion 结构进行交互操作，以终端菜单方式完成 ISODISTORT 网站核心交互流程。")
        _line()

    def _load_parent_structure(self) -> None:
        cif_path = _choose_cif(self.project_root, "Parent Structure: 请选择母相 CIF")
        self.iso.load_structure(cif_path)
        self.last_method1 = []
        self.last_method2 = None
        self.last_method3 = []

    def _set_distortion_types(self) -> None:
        self.distortion_types = _prompt_distortion_types(self.distortion_types)
        print(f"当前启用类型: {', '.join(self.distortion_types)}")

    def _main_menu(self) -> int:
        _line()
        print("Search Page")
        print("  1. 重新加载 Parent CIF")
        print("  2. 设置 Distortion Types")
        print("  3. Method 1: Search over all special k points")
        print("  4. Method 2: General method - specific k points")
        print("  5. Method 3: Arbitrary k + specified point-group/space-group + supercell")
        print("  6. Method 4: Mode decomposition of distorted structure")
        print("  7. 进入 Distortion Page（模式生成/导出/畴）")
        print("  8. 查看当前状态")
        print("  0. 退出")
        return _prompt_int("请选择", 3)

    def _run_method_1(self) -> None:
        _line()
        print("Method 1: Search over all special k points")

        crystal_system = _prompt(
            "可选 crystal system（triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic；留空表示不过滤）",
            "",
        ).strip().lower()
        crystal_system = crystal_system or None

        sg_raw = _prompt("可选 subgroup space group number（留空表示不过滤）", "").strip()
        subgroup_space_group = int(sg_raw) if sg_raw else None

        maximal_only = _prompt_yes_no("是否只保留 maximal subgroup", False)

        direct_sublattice = None
        if _prompt_yes_no("是否输入 direct sublattice 过滤 (a,b,c)", False):
            direct_sublattice = _prompt_supercell([1, 1, 1])

        result = self.iso.search_method_1(
            distortion_types=self.distortion_types,
            crystal_system=crystal_system,
            subgroup_space_group=subgroup_space_group,
            direct_sublattice=direct_sublattice,
            maximal_subgroup_only=maximal_only,
        )
        self.last_method1 = result

        print(f"Method 1 返回 {len(result)} 条候选（显示前 30 条）")
        for item in result[:30]:
            sg = item.subgroup
            print(
                f"  idx={sg.index:3d} | SG #{sg.space_group_number:<3d} {sg.space_group_symbol:<12s} "
                f"| crystal_system={item.crystal_system:<12s} | maximal={item.is_maximal}"
            )

    def _run_method_2(self) -> None:
        _line()
        print("Method 2: General method - search over specific k points")

        subgroup_idx = _prompt_int("请输入 subgroup_idx")
        dtype = _prompt(
            "distortion_type (displacement/order/strain/magnetic)",
            self.distortion_types[0] if self.distortion_types else "displacement",
        ).strip().lower()

        k_label = _prompt("k-point label（可留空）", "").strip() or None
        k_coords_raw = _prompt("k-point coordinates（可留空；格式如 1/2,0,0）", "").strip()
        k_coords = [x.strip() for x in k_coords_raw.split(",") if x.strip()] if k_coords_raw else None

        k_params: Dict[str, str] = {}
        if _prompt_yes_no("是否输入 k 参数（如 a,b,g）", False):
            while True:
                pair = _prompt("输入参数 key=value（回车结束）", "").strip()
                if not pair:
                    break
                if "=" not in pair:
                    print("格式错误，应为 key=value")
                    continue
                key, value = pair.split("=", 1)
                k_params[key.strip()] = value.strip()

        n_mod = _prompt_int("independent modulations 数（incommensurate 时常用）", 0)
        n_superposed = _prompt_int("superposed IR 数", 1)
        opd = _prompt("指定 OPD（可留空）", "").strip() or None

        result = self.iso.search_method_2(
            subgroup_idx=subgroup_idx,
            distortion_type=dtype,
            k_point_label=k_label,
            k_point_coordinates=k_coords,
            k_parameters=k_params,
            number_of_independent_modulations=n_mod,
            number_of_superposed_irs=n_superposed,
            specified_opd=opd,
        )
        self.last_method2 = result

        print(f"Method 2: 模式数 {len(result.modes)}")
        for mode in result.modes[:40]:
            print(f"  {mode.irrep_label:<10s} dim={mode.dimension}")

        if _prompt_yes_no("是否立即进入 Distortion 生成", True):
            self._generate_single_mode_flow()

    def _run_method_3(self) -> None:
        _line()
        print("Method 3: arbitrary k for specified point group and supercell")

        point_group = _prompt("point group（留空可不填）", "").strip() or None

        sg_raw = _prompt("space-group type（留空可不填）", "").strip()
        space_group_type = int(sg_raw) if sg_raw else None

        basis = None
        if _prompt_yes_no("是否输入 supercell basis matrix", True):
            basis = _prompt_basis_matrix()

        centering = _prompt("direct sublattice centering（留空可不填）", "").strip() or None

        result = self.iso.search_method_3(
            distortion_types=self.distortion_types,
            point_group=point_group,
            space_group_type=space_group_type,
            supercell_basis=basis,
            direct_sublattice_centering=centering,
        )
        self.last_method3 = result

        print(f"Method 3 返回 {len(result)} 条候选（显示前 30 条）")
        for item in result[:30]:
            sg = item.subgroup
            print(
                f"  idx={sg.index:3d} | SG #{sg.space_group_number:<3d} {sg.space_group_symbol:<12s} "
                f"| point_group={item.point_group}"
            )

    def _run_method_4(self) -> None:
        _line()
        print("Method 4: Mode decomposition of a distorted structure")
        daughter_cif = _choose_cif(self.project_root, "请选择 Daughter CIF")

        matching = _prompt("atom matching method (nearest-site/robust)", "nearest-site").strip()
        threshold = _prompt_float("robust threshold（分数坐标距离）", 0.25)

        result = self.iso.search_method_4(
            distorted_cif_path=daughter_cif,
            atom_matching_method=matching,
            robust_distance_threshold=threshold,
            provided_origin_shift=None,
        )

        print("分解结果（按绝对幅度排序，前 20 项）")
        ranked = sorted(result.amplitudes.items(), key=lambda kv: abs(kv[1]), reverse=True)
        for label, amp in ranked[:20]:
            print(f"  {label:<12s} {amp:+.8f}")
        print(f"RMS residual: {result.rms_residual:.8e}")
        print(f"Max residual: {result.max_abs_residual:.8e}")

    def _distortion_page_menu(self) -> None:
        while True:
            _line()
            print("Distortion Page")
            print("  1. Generate single-mode distortion")
            print("  2. Generate mixed-mode distortion")
            print("  3. Export current distorted structure")
            print("  4. Generate domains")
            print("  0. 返回 Search Page")
            choice = _prompt_int("请选择", 1)

            if choice == 0:
                return
            if choice == 1:
                self._generate_single_mode_flow()
            elif choice == 2:
                self._generate_mixed_mode_flow()
            elif choice == 3:
                self._export_flow()
            elif choice == 4:
                self._domains_flow()

    def _generate_single_mode_flow(self) -> None:
        if not self.iso.mode_displacements:
            print("当前没有可用模式，请先执行 Method 2 或 select_path。")
            return

        labels = list(self.iso.mode_displacements.keys())
        print("可用模式：")
        for i, label in enumerate(labels, start=1):
            print(f"  {i:2d}. {label}")

        raw = _prompt("选择模式（序号或标签）", "1")
        if raw.isdigit():
            idx = int(raw)
            if 1 <= idx <= len(labels):
                irrep = labels[idx - 1]
            else:
                print("模式序号超出范围")
                return
        else:
            irrep = raw.strip()

        amplitude = _prompt_float(
            "amplitude",
            float(self.iso.cfg.defaults.get("default_amplitude", 1.0)),
        )
        supercell = _prompt_supercell(self.iso.cfg.defaults.get("default_supercell", [1, 1, 1]))

        distorted = self.iso.generate_distortion(irrep_label=irrep, amplitude=amplitude, supercell=supercell)
        print(f"已生成结构，原子数: {len(distorted)}")

    def _generate_mixed_mode_flow(self) -> None:
        if not self.iso.mode_displacements:
            print("当前没有可用模式，请先执行 Method 2 或 select_path。")
            return

        print("请输入模式贡献，格式 irrep=amplitude，逐行输入，空行结束")
        contributions: Dict[str, float] = {}
        while True:
            row = _prompt("贡献", "").strip()
            if not row:
                break
            if "=" not in row:
                print("格式错误，应为 irrep=amplitude")
                continue
            key, value = row.split("=", 1)
            key = key.strip()
            try:
                contributions[key] = float(value.strip())
            except ValueError:
                print("amplitude 必须是数字")

        if not contributions:
            print("未输入任何贡献，取消。")
            return

        supercell = _prompt_supercell(self.iso.cfg.defaults.get("default_supercell", [1, 1, 1]))
        distorted = self.iso.generate_mixed_distortion(contributions=contributions, supercell=supercell)
        print(f"已生成混合畸变结构，原子数: {len(distorted)}")

    def _export_flow(self) -> None:
        if self.iso.distorted_structure is None:
            print("当前没有可导出的畸变结构，请先生成。")
            return

        name = _prompt("导出文件名前缀", "distorted_output")
        formats_raw = _prompt("导出格式（逗号分隔，如 cif,poscar）", "cif")
        formats = [x.strip().lower() for x in formats_raw.split(",") if x.strip()]
        paths = self.iso.export(name, formats=formats)

        print("导出完成：")
        for path in paths:
            print(f"  {path}")

    def _domains_flow(self) -> None:
        try:
            domains = self.iso.generate_domains()
        except Exception as exc:
            print(f"生成畴失败: {exc}")
            return
        print(f"生成畴数: {len(domains)}")

    def _show_state(self) -> None:
        _line()
        print("Current Session State")
        if self.iso.structure is None:
            print("  Parent structure: 未加载")
        else:
            sg = self.iso.symmetry_info["space_group_number"]
            sym = self.iso.symmetry_info["space_group_symbol"]
            n_atoms = len(self.iso.structure)
            print(f"  Parent structure: SG #{sg} ({sym}), atoms={n_atoms}")

        print(f"  Distortion types: {', '.join(self.distortion_types)}")
        print(f"  Last Method1 count: {len(self.last_method1)}")
        print(f"  Last Method3 count: {len(self.last_method3)}")

        mode_count = len(self.iso.mode_displacements)
        print(f"  Available mapped modes: {mode_count}")

        if self.iso.distorted_structure is not None:
            print(f"  Distorted structure atoms: {len(self.iso.distorted_structure)}")
        else:
            print("  Distorted structure: 未生成")


def main() -> int:
    try:
        return IsoDistortConsoleApp().run()
    except KeyboardInterrupt:
        print("\n用户中断。")
        return 130
    except Exception as exc:
        print(f"\n运行失败: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
