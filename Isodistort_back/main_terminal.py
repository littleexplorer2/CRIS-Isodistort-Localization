"""统一交互入口：以终端菜单方式复现 ISODISTORT 网站核心交互流程。

对应官网页面：
- Search Page：加载母相 CIF、选择 Distortion Types、Method 1-4
- Distortion Page：单模式/多模式畸变生成、导出、畴

说明：
- 本程序依赖 WSL（isobyu 中为 Linux ELF 二进制），首次运行会自动
  在 WSL 侧建立短路径暂存目录与 ISODATA 符号链接。
- Method 1 枚举全部特殊 k 点的各向同性子群，候选较多时计算约需
  数秒到数十秒（与官网“数据库查询”一致，属正常等待）。
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

from isocore.api import IsoDistort
from isocore.i18n import get_language, set_language, t
from isocore.utils import IsodistortError

DISTORTION_TYPE_MAP = {
    1: "displacement",
    2: "order",
    3: "strain",
    4: "magnetic",
}


def _line() -> None:
    """打印分隔线"""
    print("-" * 78)


def _prompt(text: str, default: str | None = None) -> str:
    """通用的读取用户输入的方法：提示用户输入，支持默认值"""
    hint = f" [{default}]" if default is not None else ""
    while True:  # 循环直到用户输入有效值（非空）
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
            print(t("ui.prompt.int"))


def _prompt_float(text: str, default: float | None = None) -> float:
    while True:
        raw = _prompt(text, str(default) if default is not None else None)
        try:
            return float(raw)
        except ValueError:
            print(t("ui.prompt.float"))


def _prompt_yes_no(text: str, default_yes: bool = False) -> bool:
    default = "y" if default_yes else "n"
    while True:
        raw = _prompt(f"{text} (y/n)", default).lower()
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False
        print(t("ui.prompt.yes_no"))


def _find_cif_candidates(root: Path, limit: int = 30) -> list[Path]:
    """在指定目录下查找 CIF 文件，返回前 limit 个结果，按路径排序。"""
    found: list[Path] = []
    for path in root.rglob("*.cif"):
        if path.is_file():
            found.append(path)
            if len(found) >= limit:
                break
    return sorted(found)


def _choose_cif(project_root: Path, title: str) -> str:
    """交互式选择 CIF 文件（列出候选或手动输入路径）。"""
    _line()
    print(title)
    candidates = _find_cif_candidates(project_root)

    if candidates:
        print(t("ui.msg.cif_found"))
        for idx, path in enumerate(candidates, start=1):
            try:
                display = str(path.relative_to(project_root))
            except ValueError:
                display = str(path)
            print(f"  {idx:2d}. {display}")
        print("   0. " + t("ui.msg.cif_manual"))

        while True:
            choice = _prompt_int(t("ui.msg.cif_choice"), 1)
            if choice == 0:
                break
            if 1 <= choice <= len(candidates):
                return str(candidates[choice - 1])
            print(t("ui.msg.cif_range"))

    return _prompt(t("ui.prompt.cif_manual"))


def _prompt_distortion_types(default: Sequence[str]) -> list[str]:
    """选择 Distortion Types（对应官网勾选框，逻辑 OR）。"""
    _line()
    print("Distortion Types")
    print("  1. displacement（" + t("ui.type.displacement") + "）")
    print("  2. order（" + t("ui.type.order") + "）")
    print("  3. strain（" + t("ui.type.strain") + "）")
    print("  4. magnetic（" + t("ui.type.magnetic") + "）")
    default_raw = ",".join(default)
    raw = _prompt(t("ui.prompt.types"), default_raw)

    labels: list[str] = []
    for token in [x.strip().lower() for x in raw.split(",") if x.strip()]:
        if token.isdigit() and int(token) in DISTORTION_TYPE_MAP:
            labels.append(DISTORTION_TYPE_MAP[int(token)])
        elif token in DISTORTION_TYPE_MAP.values():
            labels.append(token)

    dedup: list[str] = []
    seen = set()
    for label in labels:
        if label not in seen:
            dedup.append(label)
            seen.add(label)

    if not dedup:
        return ["displacement", "strain"]
    return dedup


def _prompt_supercell(default_cell: Sequence[int]) -> list[int]:
    raw = _prompt(t("ui.prompt.supercell"), ",".join(str(x) for x in default_cell))
    parts = [x.strip() for x in raw.split(",") if x.strip()]
    if len(parts) != 3:
        raise ValueError("超胞必须输入 3 个整数")
    values = [int(x) for x in parts]
    if any(v < 1 for v in values):
        raise ValueError("超胞参数必须 >= 1")
    return values


def _prompt_basis_matrix() -> list[list[str]]:
    print("请输入 3x3 基矢矩阵（每行 3 个值，可用分数如 1/2，用空格分隔）")
    rows: list[list[str]] = []
    for i in range(3):
        row = _prompt(f"第 {i + 1} 行", "1 0 0" if i == 0 else "0 1 0" if i == 1 else "0 0 1")
        parts = [x for x in row.split() if x]
        if len(parts) != 3:
            raise ValueError("每一行必须有 3 个值")
        rows.append(parts)
    return rows


class IsoDistortConsoleApp:
    """终端交互应用：Search Page + Distortion Page。"""

    def __init__(self) -> None:
        self.project_root = Path(__file__).resolve().parent
        self.iso = IsoDistort()
        self.distortion_types: list[str] = ["displacement", "strain"]
        self.last_method1: list = []  # Method 1 候选（Method1ResultItem 列表）
        self.last_method2 = None      # Method 2 结果
        self.last_method3: list = []  # Method 3 候选

    def run(self) -> int:
        self._banner()
        self._load_parent_structure()
        while True:
            choice = self._main_menu()
            if choice == 0:
                print(t("ui.exit.done"))
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
            elif choice == 9:
                self._switch_language()

    # ----------------------------------------------------------------
    # 基础流程
    # ----------------------------------------------------------------

    def _banner(self) -> None:
        """打印使用前说明"""
        _line()
        print(t("ui.banner.title"))
        print(t("ui.banner.desc"))
        _line()

    def _load_parent_structure(self) -> None:
        cif_path = _choose_cif(self.project_root, t("ui.prompt.cif"))
        try:
            self.iso.load_structure(cif_path)
        except Exception as exc:  # noqa: BLE001 - CLI 边界：任何加载错误都转为提示
            print(t("load.failed", err=exc))
            return
        self.last_method1 = []
        self.last_method2 = None
        self.last_method3 = []

    def _set_distortion_types(self) -> None:
        self.distortion_types = _prompt_distortion_types(self.distortion_types)
        print(f"当前启用类型: {', '.join(self.distortion_types)}")

    def _switch_language(self) -> None:
        """切换界面语言（zh / en / mixed）。"""
        print(t("ui.lang.current", lang=get_language()))
        raw = _prompt(t("ui.lang.prompt"), "").strip().lower()
        if raw in ("zh", "en", "mixed"):
            set_language(raw)
            print(t("ui.lang.current", lang=get_language()))
        else:
            print(t("ui.lang.invalid"))

    def _main_menu(self) -> int:
        """打印 Search Page 主菜单并获取用户选择"""
        _line()
        print(t("ui.search_page"))
        print(t("ui.menu.reload"))
        print(t("ui.menu.types"))
        print(t("ui.menu.method1"))
        print(t("ui.menu.method2"))
        print(t("ui.menu.method3"))
        print(t("ui.menu.method4"))
        print(t("ui.menu.distortion"))
        print(t("ui.menu.state"))
        print(t("ui.menu.language"))
        print(t("ui.menu.exit"))
        return _prompt_int(t("ui.prompt.action"), 3)

    # ----------------------------------------------------------------
    # Method 1-4
    # ----------------------------------------------------------------

    def _run_method_1(self) -> None:
        _line()
        print("Method 1: Search over all special k points")

        crystal_system = _prompt(
            t("m1.cs"),
            "",
        ).strip().lower()
        crystal_system = crystal_system or None

        sg_raw = _prompt(t("m1.sg"), "").strip()
        subgroup_space_group = int(sg_raw) if sg_raw else None

        maximal_only = _prompt_yes_no(t("m1.maximal"), False)

        direct_sublattice = None
        if _prompt_yes_no(t("m1.sublattice"), False):
            direct_sublattice = _prompt_supercell([1, 1, 1])

        print(t("m1.wait"))
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
                f"  idx={sg.index:3d} | SG #{sg.space_group_number:<3d} "
                f"{sg.space_group_symbol:<12s} | k={sg.k_point_label:<4s} "
                f"IR={sg.irrep_label:<6s} OPD={sg.opd_symbol:<4s} "
                f"| crystal_system={item.crystal_system:<12s} "
                f"| maximal={item.is_maximal}"
            )
        if len(result) > 30:
            print(f"  ... 还有 {len(result) - 30} 条")

    def _run_method_2(self) -> None:
        _line()
        print("Method 2: General method - search over specific k points")
        print(t("m2.mode1"))
        print(t("m2.mode2"))
        mode = _prompt_int(t("m2.choose_mode"), 1)

        if mode == 2:
            subgroup_idx = self._direct_kpoint_search()
        elif not self.last_method1:
            print(t("m2.no_candidates"))
            return
        else:
            subgroup_idx = _prompt_int(
                f"subgroup_idx（Method 1 候选序号，"
                f"范围 0-{len(self.last_method1) - 1}）"
            )
        if subgroup_idx is None:
            return

        dtype = _prompt(
            "distortion_type (displacement/order/strain/magnetic)",
            self.distortion_types[0] if self.distortion_types else "displacement",
        ).strip().lower()

        # k 点 / OPD 信息（官网在 Method 2 页选择 k 点、IR、OPD；
        # 本地已由枚举确定，此处仅作展示与确认）
        target = next(
            (sg for sg in self.iso.subgroups
             if sg.index == subgroup_idx), None
        )
        if target is None:
            print(t("m2.idx_range"))
            return
        print(
            f"已选子群: {target.describe()} "
            f"(index={target.subgroup_index}, size={target.size})"
        )

        n_mod = _prompt_int("independent modulations 数（incommensurate 时常用）", 0)
        n_superposed = _prompt_int("superposed IR 数", 1)
        opd = _prompt("指定 OPD（留空使用子群默认 OPD）", "").strip() or None

        result = self.iso.search_method_2(
            subgroup_idx=subgroup_idx,
            distortion_type=dtype,
            k_point_label=target.k_point_label,
            k_parameters={},
            number_of_independent_modulations=n_mod,
            number_of_superposed_irs=n_superposed,
            specified_opd=opd,
        )
        self.last_method2 = result

        print(f"Method 2: 模式数 {len(result.modes)}")
        for mode in result.modes:
            sites = sorted({b.wyckoff_letter for b in mode.bush_modes})
            print(
                f"  {mode.irrep_label:<8s} OPD={mode.opd_symbol:<6s} "
                f"dim={mode.dimension:<2d} 位点={sites}"
            )

        if _prompt_yes_no(t("m2.enter_dist"), True):
            self._generate_single_mode_flow()

    def _direct_kpoint_search(self) -> int | None:
        """直接 k 点搜索：k 点 -> IR -> OPD -> 子群（对齐官网 Method 2 流程）。"""
        print("\n--- 选择 k 点 / Choose a k point ---")
        kpoints = self.iso.list_k_points()
        for i, kp in enumerate(kpoints, start=1):
            params = f"（参数: {','.join(kp.parameters)}）" if kp.parameters else ""
            print(f"  {i:2d}. {kp.label:<4s} {kp.coordinates}{params}")
        choice = _prompt_int(t("m2.choose_k"), 1)
        if not (1 <= choice <= len(kpoints)):
            print(t("m2.range"))
            return None
        kp = kpoints[choice - 1]

        k_parameters: list[str] | None = None
        if kp.parameters:
            print(
                f"k 点 {kp.label} 坐标为 {kp.coordinates}，"
                f"需设置参数（按顺序输入，官网要求有理数如 1/2）"
            )
            print("提示：iso 的参数约定与官网可能差整数倍（如 iso 用 2a、官网用 g），"
                  "请根据 k 点坐标形式自行换算。")
            values: list[str] = []
            for p in kp.parameters:
                values.append(_prompt(t("m2.param_value", p=p), "").strip())
            k_parameters = values if all(values) else None
            if k_parameters is None:
                print(t("m2.param_empty"))
                return None

        print("\n--- 选择不可约表示 / Choose an IR ---")
        irreps = self.iso.list_irreps(kp.label, k_parameters)
        for i, ir in enumerate(irreps, start=1):
            print(f"  {i:2d}. {ir.label:<8s} dim={ir.dimension}")
        choice = _prompt_int(t("m2.choose_ir"), 1)
        if not (1 <= choice <= len(irreps)):
            print(t("m2.range"))
            return None
        irrep = irreps[choice - 1]

        print(t("m2.enum_wait"))
        try:
            subs = self.iso.list_subgroups_at(
                kp.label, irrep.label, k_parameters=k_parameters,
            )
        except IsodistortError as exc:
            if "需要在线生成" in str(exc) and _prompt_yes_no(
                "子群数据库不存在，是否在线生成（官网 “Generate isotropy "
                "subgroups”，可能耗时数分钟到数小时）",
                False,
            ):
                subs = self.iso.list_subgroups_at(
                    kp.label, irrep.label, k_parameters=k_parameters,
                    generate_if_missing=True,
                )
            else:
                print(f"枚举失败: {exc}")
                return None
        if not subs:
            print(t("m2.no_subs"))
            return None
        for i, sg in enumerate(subs, start=1):
            print(
                f"  {i:2d}. idx={sg.index} SG {sg.space_group_number} "
                f"{sg.space_group_symbol:<10s} OPD={sg.opd_symbol:<4s} "
                f"index={sg.subgroup_index} size={sg.size}"
            )
        choice = _prompt_int(t("m2.choose_sg"), 1)
        if not (1 <= choice <= len(subs)):
            print(t("m2.range"))
            return None
        return subs[choice - 1].index

    def _run_method_3(self) -> None:
        _line()
        print("Method 3: Search over arbitrary k for point group and supercell")

        point_group = _prompt(t("m3.pg"), "").strip() or None

        sg_raw = _prompt(t("m3.sg"), "").strip()
        space_group_type = int(sg_raw) if sg_raw else None

        basis = None
        if _prompt_yes_no(t("m3.basis_q"), True):
            basis = _prompt_basis_matrix()

        centering = _prompt(t("m3.centering"), "").strip() or None

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
                f"  idx={sg.index:3d} | SG #{sg.space_group_number:<3d} "
                f"{sg.space_group_symbol:<12s} | k={sg.k_point_label:<4s} "
                f"IR={sg.irrep_label:<6s} | point_group={item.point_group}"
            )
        if len(result) > 30:
            print(f"  ... 还有 {len(result) - 30} 条")

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

        print(t("m4.result_title"))
        ranked = sorted(result.amplitudes.items(), key=lambda kv: abs(kv[1]), reverse=True)
        for label, amp in ranked[:20]:
            print(f"  {label:<12s} {amp:+.8f}")
        print(f"RMS residual: {result.rms_residual:.8e}")
        print(f"Max residual: {result.max_abs_residual:.8e}")

    # ----------------------------------------------------------------
    # Distortion Page
    # ----------------------------------------------------------------

    def _distortion_page_menu(self) -> None:
        while True:
            _line()
            print(t("ui.distortion_page"))
            print(t("ui.dist.single"))
            print(t("ui.dist.mixed"))
            print(t("ui.dist.export"))
            print(t("ui.dist.domains"))
            print(t("ui.dist.back"))
            choice = _prompt_int(t("ui.dist.choice"), 1)

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
            print(t("ui.msg.no_modes"))
            return

        labels = list(self.iso.mode_displacements.keys())
        print("可用模式 / Available modes:")
        for i, label in enumerate(labels, start=1):
            print(f"  {i:2d}. {label}")

        raw = _prompt(t("dist.choose_mode"), "1")
        if raw.isdigit():
            idx = int(raw)
            if 1 <= idx <= len(labels):
                irrep = labels[idx - 1]
            else:
                print(t("dist.mode_range"))
                return
        else:
            irrep = raw.strip()

        amplitude = _prompt_float(
            "amplitude",
            float(self.iso.cfg.defaults.get("default_amplitude", 1.0)),
        )
        supercell = None
        if self.iso.phase_path is not None and self.iso.phase_path.basis_vectors:
            print("将使用所选子群的超胞基矢 / Using the subgroup supercell basis.")
        elif _prompt_yes_no(t("dist.supercell_q"), False):
            supercell = _prompt_supercell(
                self.iso.cfg.defaults.get("default_supercell", [1, 1, 1])
            )

        distorted = self.iso.generate_distortion(
            irrep_label=irrep, amplitude=amplitude, supercell=supercell
        )
        print(t("ui.msg.generated", n=len(distorted)))

    def _generate_mixed_mode_flow(self) -> None:
        if not self.iso.mode_displacements:
            print(t("ui.msg.no_modes"))
            return

        print(t("dist.contrib"))
        contributions: dict[str, float] = {}
        while True:
            row = _prompt(t("dist.contrib_prompt"), "").strip()
            if not row:
                break
            if "=" not in row:
                print(t("dist.format_err"))
                continue
            key, value = row.split("=", 1)
            key = key.strip()
            try:
                contributions[key] = float(value.strip())
            except ValueError:
                print(t("dist.amp_num"))

        if not contributions:
            print(t("dist.no_contrib"))
            return

        supercell = None
        if self.iso.phase_path is not None and self.iso.phase_path.basis_vectors:
            print("将使用所选子群的超胞基矢 / Using the subgroup supercell basis.")
        elif _prompt_yes_no(t("dist.supercell_q"), False):
            supercell = _prompt_supercell(
                self.iso.cfg.defaults.get("default_supercell", [1, 1, 1])
            )

        distorted = self.iso.generate_mixed_distortion(
            contributions=contributions, supercell=supercell
        )
        print(t("ui.msg.generated_mixed", n=len(distorted)))

    def _export_flow(self) -> None:
        if self.iso.distorted_structure is None:
            print(t("ui.msg.no_export"))
            return

        name = _prompt("导出文件名前缀 / File name prefix", "distorted_output")
        formats_raw = _prompt("导出格式（逗号分隔，如 cif,poscar）", "cif")
        formats = [x.strip().lower() for x in formats_raw.split(",") if x.strip()]
        paths = self.iso.export(name, formats=formats)

        print(t("ui.msg.export_done"))
        for path in paths:
            print(f"  {path}")

    def _domains_flow(self) -> None:
        try:
            domains = self.iso.generate_domains()
        except IsodistortError as exc:
            print(f"生成畴失败 / Domains failed: {exc}")
            return
        print(t("ui.msg.domains_done", n=len(domains)))
        for d in domains[:20]:
            print(
                f"  domain {d.domain_number:3d} | gen={d.generator:<14s} "
                f"SG {d.space_group_number} {d.space_group_symbol}"
            )
        if len(domains) > 20:
            print(f"  ... 还有 {len(domains) - 20} 个")

    # ----------------------------------------------------------------
    # 状态
    # ----------------------------------------------------------------

    def _show_state(self) -> None:
        _line()
        print(t("ui.msg.state"))
        if self.iso.structure is None:
            print(f"  Parent structure: {t('state.not_loaded')}")
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
            print(f"  Distorted structure: {t('state.not_generated')}")


def main() -> int:
    # Windows 控制台默认代码页可能不是 UTF-8，统一重配置避免中文乱码
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        return IsoDistortConsoleApp().run()
    except KeyboardInterrupt:
        print(f"\n{t('ui.exit.interrupt')}")
        return 130
    except IsodistortError as exc:
        print(f"\n运行失败: {exc}")
        return 1
    except Exception as exc:  # noqa: BLE001 - CLI 入口兜底，防止未捕获异常直接崩溃
        print(f"\n运行失败: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
