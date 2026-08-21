"""
核心 Python API - 一站式 ISODISTORT 使用入口

对应全阶段：封装官网 Search Page → Distortion Page 的完整交互流程。

使用示例：
    from isocore.api import IsoDistort

    iso = IsoDistort()
    # 1. 从 CIF 加载结构
    iso.load_structure("input.cif")
    # 2. Method 1：枚举全部特殊 k 点子群并过滤
    candidates = iso.search_method_1(crystal_system="tetragonal")
    # 3. Method 2：选择子群并计算畸变模式
    result = iso.search_method_2(subgroup_idx=candidates[0].subgroup.index)
    # 4. Distortion Page：生成畸变结构
    distorted = iso.generate_distortion(amplitude=0.1)
    # 5. 导出
    iso.export("output", formats=["cif", "poscar"])
"""
import threading
from fractions import Fraction
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup

from ..backend import (
    DistortionMode,
    FindsymWrapper,
    IsoWrapper,
    SubgroupInfo,
)
from ..distortion import (
    DISTORTION_TYPES,
    DistortionEngine,
    DistortionMapper,
    DomainGenerator,
    IsoSearchEngine,
    Method1Query,
    Method2Query,
    Method3Query,
    Method4Query,
    OccupationalModeGenerator,
    PhasePath,
    normalize_distortion_types,
)
from ..i18n import get_language, set_language, t
from ..io import StructureExporter
from ..structure import (
    SymmetryValidator,
    read_cif,
    read_structure,
)
from ..utils import IsodistortError, get_config
from ..utils.text_parser import parse_fraction


class IsoDistort:
    """
    ISODISTORT 主入口类

    封装完整工作流：
        加载结构 → 识别对称 → 枚举子群（Method 1）→ 选择路径（Method 2）
        → 计算畸变模式 → 生成畸变结构 → 导出/畴

    语言支持：``IsoDistort(language="en")`` 或 ``iso.set_language("zh")``
    可随时切换控制台输出语言（终端/网页端通过 isocore.i18n 全局切换）。
    """

    def __init__(self, language: str | None = None) -> None:
        self.cfg = get_config()

        # 界面语言（None 时取配置 runtime.language，默认 zh）
        if language is not None:
            set_language(language)
        else:
            set_language(self.cfg.language)

        # 底层封装
        self._findsym = FindsymWrapper()
        self._iso = IsoWrapper()

        # 业务层
        self._sym_val = SymmetryValidator()
        self._dist_mapper = DistortionMapper()
        self._dist_engine = DistortionEngine(self._dist_mapper)
        self._domain_gen = DomainGenerator(self._iso)
        self._search = IsoSearchEngine(self._iso)

        # 输出层
        self._exporter = StructureExporter()

        # 状态
        self.structure: Structure | None = None
        self.symmetry_info: dict | None = None
        self.subgroups: list[SubgroupInfo] = []
        self.phase_path: PhasePath | None = None
        self.distortion_modes: list[DistortionMode] = []
        self.mode_displacements: dict = {}
        self.mode_occupancies: dict = {}          # occupational 模式（占据率调制）
        self.distorted_structure: Structure | None = None

        # 畸变类型作用域（对齐官网 per-species 复选框）：type -> 物种列表（"*"=全部）
        self.distortion_scope: dict[str, list[str]] = {}
        self._special_subgroups_cache: list | None = None
        self._special_subgroups_lock = threading.Lock()
        self._conv_to_prim_cache: np.ndarray | None = None

    def set_language(self, language: str) -> None:
        """切换控制台/界面语言（"zh" 中文 / "en" English）。"""
        set_language(language)

    def get_language(self) -> str:
        """当前语言（"zh" 或 "en"）。"""
        return get_language()

    # ================================================================
    # 阶段一：结构输入与对称识别
    # ================================================================

    def load_structure(self, cif_path: str | Path) -> Structure:
        """
        步骤1-3：加载结构文件，识别对称性

        支持 CIF / VASP POSCAR / xyz（按扩展名自动识别，见 read_structure）。

        Args:
            cif_path: 结构文件路径（CIF / POSCAR / xyz）

        Returns:
            Structure: 加载后的结构对象
        """
        path = Path(cif_path)
        self.structure = (read_cif(path) if path.suffix.lower() == ".cif"
                          else read_structure(path))
        self.symmetry_info = self._sym_val.validate(self.structure)
        self._reset_derived_state()

        sg_num = self.symmetry_info["space_group_number"]
        sg_sym = self.symmetry_info["space_group_symbol"]
        n_atoms = len(self.structure)

        print(t("load.done", sg=sg_num, sym=sg_sym, n=n_atoms))
        return self.structure

    def set_structure(self, structure: Structure) -> Structure:
        """直接设置 Structure 对象"""
        self.structure = structure
        self.symmetry_info = self._sym_val.validate(self.structure)
        self._reset_derived_state()
        return self.structure

    def _reset_derived_state(self) -> None:
        """加载新结构后清空所有派生状态。"""
        self.subgroups = []
        self.phase_path = None
        self.distortion_modes = []
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.distorted_structure = None
        self._special_subgroups_cache = None
        self._conv_to_prim_cache = None
        self.distorted_structure = None
        self._preferences = None

    def _wyckoff_letters(self) -> list[str]:
        """当前结构各 Wyckoff 位置字母（用于 iso 模式计算）。"""
        if not self.symmetry_info:
            return []
        return [s["wyckoff_letter"] for s in self.symmetry_info["wyckoff_sites"]]

    # ================================================================
    # 畸变类型作用域（对齐官网 Types 面板的 per-species 复选框）
    # ================================================================

    def species(self) -> list[str]:
        """当前结构包含的元素符号（去重、排序）。"""
        if self.structure is None:
            return []
        return sorted({s.species_string for s in self.structure})

    def set_distortion_scope(self, scope: dict | None) -> None:
        """
        设置各畸变类型的作用域物种（官网 Displacive/Occupational/Magnetic/
        Rotational 行内的 all/none/Eu/Al 复选框）。

        Args:
            scope: {类型名: "*"（全部）或 ["Eu", "Al", ...]}；
                不在 scope 中的类型默认作用于全部物种
        """
        self.distortion_scope = {}
        if not scope:
            return
        valid = set(DISTORTION_TYPES)
        for tp, val in scope.items():
            if tp not in valid:
                continue
            if isinstance(val, str) and val.strip().lower() in ("*", "all", "全部"):
                self.distortion_scope[tp] = ["*"]
            elif isinstance(val, (list, tuple)):
                species = [str(s) for s in val if str(s) != "*"]
                self.distortion_scope[tp] = species if species else ["*"]

    def _scope_species(self, type_name: str) -> set[str]:
        """某畸变类型作用域内的物种集合（未设置时默认全部物种）。"""
        all_species = set(self.species())
        val = self.distortion_scope.get(type_name)
        if not val or "*" in val:
            return all_species
        return {s for s in val if s in all_species}

    def _letters_for_species(self, species_set) -> list[str]:
        """物种集合对应的 Wyckoff 位置字母（去重）。"""
        if not self.symmetry_info:
            return []
        letters: list[str] = []
        for site in self.symmetry_info["wyckoff_sites"]:
            if site["species"] in species_set and site["wyckoff_letter"] not in letters:
                letters.append(site["wyckoff_letter"])
        return letters

    # ================================================================
    # Method 1 下拉数据（对齐官网：可达子群空间群 + conventional/primitive lattice）
    # ================================================================

    def space_group_preferences(self) -> str:
        """官网页头 “Default space-group preferences: ...” 行。

        本地引擎采用官网默认值（国际取位/标准设定），与实际计算行为一致。
        """
        return ("monoclinic axes a(b)c, monoclinic cell choice 1, "
                "orthorhombic axes abc, origin choice 2, hexagonal axes, "
                "SSG standard setting")

    def set_space_group_preferences(self, preferences: dict | None = None) -> dict:
        """记录 Space-Group Preferences 选择（对齐官网 settings 面板交互）。

        官网允许用户修改这些偏好并影响后续子群枚举的取位。本地 iso 二进制
        固定使用国际标准取位（International, new ed. with conventional basis
        vectors），不支持自定义 setting（Syntax error），因此**计算仍采用
        官网默认值**；此处仅记录用户选择供界面展示与状态查询，并返回当前
        生效的偏好说明。

        Args:
            preferences: 可选 dict，键为 settingaxesm / settingcell /
                settingaxeso / settingaxesh / settingorigin / settingssg /
                parentlike（与官网表单字段名一致）

        Returns:
            dict: 当前生效的偏好（默认值 + 用户选择记录）
        """
        defaults = {
            "settingaxesm": "a(b)c", "settingcell": "1",
            "settingaxeso": "abc", "settingaxesh": "h",
            "settingorigin": "2", "settingssg": "standard",
            "parentlike": False,
        }
        self._preferences = dict(defaults)
        if preferences:
            for k in defaults:
                if k in preferences and preferences[k] is not None:
                    self._preferences[k] = preferences[k]
        # 本地 iso 固定国际标准取位：自定义偏好仅记录，不影响计算
        self._preferences["effective"] = "international (default)"
        return dict(self._preferences)

    def _ensure_special_subgroups(self) -> list:
        """枚举全部特殊 k 点子群（线程安全缓存，Method 1 下拉与搜索共用）。"""
        if self.structure is None:
            raise RuntimeError(t("err.load_first"))
        with self._special_subgroups_lock:
            if self._special_subgroups_cache is None:
                sg_num = self.symmetry_info["space_group_number"]
                self._special_subgroups_cache = self._iso.enumerate_all_special_subgroups(
                    sg_num
                )
        return self._special_subgroups_cache

    def _conv_to_prim(self) -> np.ndarray:
        """母相惯用格子 -> 原胞格子的变换矩阵 T（L_prim = L_conv @ T）。"""
        if self._conv_to_prim_cache is None:
            sga = SpacegroupAnalyzer(self.structure)
            self._conv_to_prim_cache = np.asarray(
                sga.get_conventional_to_primitive_transformation_matrix(),
                dtype=float,
            )
        return self._conv_to_prim_cache

    def lattice_in_conventional_frame(self, matrix, frame: str = "conventional"
                                      ) -> list[list[float]]:
        """把用户选择的 lattice 基矢换算到 iso 惯用（conventional）坐标系。

        官网下拉的 Conventional lattice 选项即 iso 输出基矢（惯用坐标）；
        Primitive lattice 选项为同一超胞在原胞坐标下的表达
        （B_conv = B_prim @ T，见 _conv_to_prim）。
        """
        m = np.asarray(matrix, dtype=float)
        if frame == "primitive":
            m = m @ self._conv_to_prim()
        return m.tolist()

    @staticmethod
    def _same_lattice(a, b) -> bool:
        """两个 3x3 超胞基矢是否生成同一格点（GL(3,Z) 等价）。

        判定：行列式绝对值相等，且 B_a @ inv(B_b) 为整数矩阵
        （此时 B_a = B_b @ U，U 为幺模整数矩阵 -> 同一格点）。
        对含分数坐标的基矢（如 (-1/2,1/2,1/2)）同样成立。
        """
        a_arr = np.asarray(a, dtype=float)
        b_arr = np.asarray(b, dtype=float)
        if a_arr.shape != (3, 3) or b_arr.shape != (3, 3):
            return False
        if abs(abs(np.linalg.det(a_arr)) - abs(np.linalg.det(b_arr))) > 1e-6:
            return False
        try:
            n = a_arr @ np.linalg.inv(b_arr)
        except np.linalg.LinAlgError:
            return False
        return bool(np.allclose(n, np.round(n), atol=1e-5))

    def _distinct_lattices(self, bases, to_conventional: np.ndarray | None = None
                           ) -> list[dict]:
        """从一组超胞基矢中提取去重后的 lattice 选项（对齐官网下拉）。

        去重按“格点等价”（GL(3,Z) 幺模变换）判定——旧实现只做行排序/符号
        归一，无法合并同一格点的不同基矢表达（如行置换、幺模变换），
        导致选项数量远超官网。
        每类的显示代表取该类中“最简”的基矢（元素平方和最小，破平按字典序）。

        参数:
            to_conventional: 分类后在“原胞坐标”下进行（basis @ T⁻¹），
                显示时需转回惯用坐标（best @ to_conventional），与官网
                Primitive lattice 下拉一致（官网选项为惯用坐标表达，
                如原胞本身显示为 (-1/2,1/2,1/2),... 而非原胞坐标矩阵）。

        顺序：保持首次出现顺序（官网下拉顺序 = 子群数据库枚举顺序，
        不按行列式排序——旧实现按 det 排序会打乱官网顺序）。
        """
        classes: list[list[np.ndarray]] = []
        for b in bases:
            arr = np.asarray(b, dtype=float)
            if arr.shape != (3, 3):
                continue
            for cls in classes:
                if self._same_lattice(arr, cls[0]):
                    cls.append(arr)
                    break
            else:
                classes.append([arr])

        def _norm(row: np.ndarray) -> np.ndarray:
            r = np.round(row, 6)
            for x in r:
                if abs(x) > 1e-6:
                    if x < 0:
                        r = -r
                    break
            return r

        options = []
        for cls in classes:
            # 选“最简”代表（元素平方和最小；并列按排序后字典序）
            best = min(cls, key=lambda m: (
                float(np.sum(m * m)),
                tuple(float(x) for x in np.sort(np.round(m, 6), axis=0).flatten()),
            ))
            if to_conventional is not None:
                # 分类在原胞坐标下进行，显示前转回惯用坐标（官网同款）
                best = best @ to_conventional
            key = tuple(tuple(float(x) for x in _norm(r)) for r in sorted(best, key=tuple))
            options.append((round(abs(np.linalg.det(best)), 6), key))
        # 保持首次出现的枚举顺序（官网下拉 = 子群数据库枚举顺序，不按 det 排序）
        return [
            {
                "label": self._format_lattice(item[1]),
                "basis": [list(r) for r in item[1]],
            }
            for item in options
        ]

    @staticmethod
    def _format_lattice(key: tuple) -> str:
        """把 3x3 基矢渲染为官网风格的 "(1,0,0),(0,1,0),(0,0,1)"。"""

        def fmt_row(row):
            parts = []
            for x in row:
                f = Fraction(float(x)).limit_denominator(12)
                parts.append(str(f.numerator) if f.denominator == 1
                             else f"{f.numerator}/{f.denominator}")
            return "(" + ",".join(parts) + ")"

        return ",".join(fmt_row(r) for r in key)

    def method1_options(self) -> dict:
        """
        Method 1 下拉数据（对齐官网搜索页）：
        - space_groups：可达子群空间群（官网只列出与母相结构相容的对称性，
          不显示全部 230 个）
        - conventional_lattices / primitive_lattices：官网 Conventional lattice
          与 Primitive lattice 下拉选项（由真实枚举得到的超胞基矢生成）
        """
        subs = self._ensure_special_subgroups()
        numbers: list[int] = []
        for sg in subs:
            if sg.space_group_number not in numbers:
                numbers.append(sg.space_group_number)
        numbers.sort()  # 官网下拉按序号升序
        space_groups = [
            {"number": n, "symbol": SpaceGroup.from_int_number(n).symbol}
            for n in numbers
        ]
        conventional = self._distinct_lattices([sg.basis_vectors for sg in subs])
        t_prim = np.linalg.inv(self._conv_to_prim())
        # Primitive lattice：在原胞坐标下分类（basis @ T⁻¹，整数矩阵），
        # 显示前转回惯用坐标（best @ T）——对齐官网 isoplattice 下拉
        # （如“原胞本身”显示为 (-1/2,1/2,1/2),(1/2,-1/2,1/2),(1/2,1/2,-1/2)）。
        primitive = self._distinct_lattices(
            [np.asarray(sg.basis_vectors, dtype=float) @ t_prim for sg in subs],
            to_conventional=self._conv_to_prim(),
        )
        return {
            "space_groups": space_groups,
            "conventional_lattices": conventional,
            "primitive_lattices": primitive,
        }

    # ================================================================
    # 阶段二：子群枚举（Method 1 / Method 2 数据源）
    # ================================================================

    def list_k_points(self) -> list:
        """
        Method 2 数据源：枚举母相的全部 k 点（官网 Method 2 的 k 点下拉列表）。

        Returns:
            List[KPointInfo]
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        return self._iso.list_k_points(self.symmetry_info["space_group_number"])

    def _resolve_k_vector(self, k_point_label: str) -> list[float]:
        """把 k 点标签解析为数值坐标（母相倒格分数单位）。

        仅支持无自由参数的特殊 k 点（如 GM=(0,0,0)、M=(1/2,1/2,0)）；
        带参数 k 点（如 LD 的 a/b/g）无法在超胞副本间确定相位，
        返回空列表（引擎将按 k=Γ 处理，即不调制）。
        """
        if not k_point_label:
            return []
        try:
            kpoints = self._iso.list_k_points(
                self.symmetry_info["space_group_number"])
        except Exception:  # noqa: BLE001 - 解析失败按无 k 向量处理
            return []
        for kp in kpoints:
            if kp.label != k_point_label:
                continue
            if not kp.is_special:
                return []
            try:
                return [parse_fraction(c) for c in kp.coordinates]
            except (ValueError, IndexError):
                return []
        return []

    def list_irreps(self, k_point_label: str,
                    k_parameters: list | None = None) -> list:
        """
        Method 2 数据源：枚举指定 k 点下的不可约表示（官网的 IR 下拉列表）。

        Args:
            k_point_label: k 点标签（Miller-Love 记号）
            k_parameters: k 点参数值（带参数 k 点必须提供，如 ["1/6"]）

        Returns:
            List[IrrepInfo]
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        return self._iso.list_irreps(
            self.symmetry_info["space_group_number"], k_point_label, k_parameters
        )

    def list_subgroups_at(self, k_point_label: str, irrep_label: str,
                          k_parameters: list | None = None,
                          opd_symbol: str | None = None,
                          generate_if_missing: bool = False) -> list:
        """
        Method 2 数据源：枚举指定 (k 点, IR) 下的各向同性子群，
        并记录为当前候选列表（供 search_method_2 使用）。

        Args:
            k_point_label: k 点标签
            irrep_label: 不可约表示标签
            k_parameters: k 点参数（带参数 k 点必须提供）
            opd_symbol: 只返回该序参量方向对应的子群（可选）
            generate_if_missing: 子群数据库缺失时是否自动在线生成
                （默认 False；生成可能耗时数分钟到数小时，请谨慎开启）

        Returns:
            List[SubgroupInfo]
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        self.subgroups = self._iso.list_subgroups(
            self.symmetry_info["space_group_number"],
            k_point_label,
            irrep_label,
            k_parameters=k_parameters,
            opd_symbol=opd_symbol,
            generate_if_missing=generate_if_missing,
        )
        return self.subgroups

    def list_subgroups_at_kpoint(self, k_point_label: str,
                                 k_parameters: list | None = None,
                                 generate_if_missing: bool = False) -> list:
        """
        Method 2 数据源：枚举指定 k 点下**全部不可约表示**的各向同性子群。

        对齐官网 Method 2（General method - search over specific k points）：
        官网只选择 k 点（及其参数 a/b/g）后提交，即返回该 k 点全部 IR 的
        子群列表（不预先选择 IR / OPD）。

        Args:
            k_point_label: k 点标签
            k_parameters: k 点参数（带参数 k 点必须提供）
            generate_if_missing: 子群数据库缺失时是否自动在线生成

        Returns:
            List[SubgroupInfo]（按 k 点下 IR 的枚举顺序）
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        irreps = self._iso.list_irreps(
            self.symmetry_info["space_group_number"], k_point_label, k_parameters
        )
        merged: list = []
        for ir in irreps:
            try:
                subs = self._iso.list_subgroups(
                    self.symmetry_info["space_group_number"],
                    k_point_label, ir.label,
                    k_parameters=k_parameters,
                    opd_symbol=None,
                    generate_if_missing=generate_if_missing,
                )
            except IsodistortError:
                continue  # 无子群的 IR（数据库缺失等）跳过
            merged.extend(subs)
        self.subgroups = merged
        return merged

    def list_subgroups(self, distortion_type: str | list[str] | None = None
                       ) -> list[SubgroupInfo]:
        """
        步骤4：枚举母相全部特殊 k 点的各向同性子群（Method 1 数据源）。

        结果按会话缓存（Method 1 下拉与搜索共用，避免重复枚举）。

        Args:
            distortion_type: 畸变类型（保留参数；类型过滤在模式计算阶段执行）

        Returns:
            List[SubgroupInfo]
        """
        if self.structure is None:
            raise RuntimeError(t("err.load_first"))

        self.subgroups = self._ensure_special_subgroups()

        print(t("subgroups.found", n=len(self.subgroups)))
        for sg in self.subgroups[:10]:
            print(f"  {sg.describe()}")
        if len(self.subgroups) > 10:
            print(t("subgroups.more", n=len(self.subgroups) - 10))
        return self.subgroups

    # ================================================================
    # 阶段三：路径选择与畸变模式计算（Method 2）
    # ================================================================

    def select_path(self, subgroup_idx: int,
                    distortion_type: str | list[str] | None = None) -> PhasePath:
        """
        步骤5-6：选择相变路径，计算畸变模式（Mode Basis）。

        Args:
            subgroup_idx: 子群序号（来自 Method 1 候选）
            distortion_type: 畸变类型（单个或列表，默认 displacive/strain；
                按 self.distortion_scope 限定物种作用域）

        Returns:
            PhasePath
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        if not self.subgroups:
            self.list_subgroups(distortion_type)

        if distortion_type is None:
            distortion_type = "displacive"

        target = next((sg for sg in self.subgroups if sg.index == subgroup_idx), None)
        if target is None:
            raise ValueError(t("subgroup.not_found", idx=subgroup_idx))

        types = normalize_distortion_types(distortion_type)
        self.phase_path = PhasePath.from_subgroup(
            self.symmetry_info["space_group_number"],
            target,
            types,
        )
        self.phase_path.validate()

        print(t("path.selected", desc=self.phase_path.describe()))

        # 计算畸变模式（按类型 + 物种作用域）
        self.distortion_modes = self._compute_scoped_modes(
            self.symmetry_info["space_group_number"], target, types
        )

        print(t("modes.found", n=len(self.distortion_modes)
                + len(self.mode_occupancies)))
        for m in self.distortion_modes:
            n_sites = len({b.wyckoff_letter for b in m.bush_modes})
            print(t("mode.sites", irrep=m.irrep_label, opd=m.opd_symbol, n=n_sites))

        # 映射到原子位移（仅 displacive/rotational 位移模式）
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )

        return self.phase_path

    def _union_scope_species(self, types: list[str]) -> set[str]:
        """全部启用类型（除 strain 外）作用域物种的并集。"""
        scoped: set[str] = set()
        for tp in types:
            if tp == "strain":
                continue
            scoped |= self._scope_species(tp)
        return scoped

    def _compute_scoped_modes(self, parent_sg: int, target: SubgroupInfo,
                              types: list[str],
                              raw_modes: list[DistortionMode] | None = None
                              ) -> list[DistortionMode]:
        """
        按畸变类型 + 物种作用域整理子群路径的模式：
        - displacive / rotational：BUSH 位移模式（raw_modes 已按作用域
          Wyckoff 位置计算），再按各类型作用域过滤并标注 mode_type；
        - occupational：本地占据率模式生成器（+1/-1 交替占据，见
          OccupationalModeGenerator），结果存入 self.mode_occupancies；
        - strain / magnetic：本地引擎暂不产生对应模式（见 README 已知差异）。
        """
        self.mode_occupancies = {}
        bush_types = [tp for tp in types if tp in ("displacive", "rotational")]

        modes: list[DistortionMode] = []
        if raw_modes is None and bush_types:
            letters = self._letters_for_species(self._union_scope_species(types))
            if letters:
                raw_modes = self._iso.calc_distortion_modes(parent_sg, target, letters)
        if raw_modes:
            allowed_letters: set[str] = set()
            for tp in bush_types:
                allowed_letters |= set(self._letters_for_species(self._scope_species(tp)))
            for m in raw_modes:
                mode_letters = {b.wyckoff_letter for b in m.bush_modes}
                if mode_letters & allowed_letters:
                    m.mode_type = "displacive" if "displacive" in bush_types else "rotational"
                    modes.append(m)

        if "occupational" in types:
            generator = OccupationalModeGenerator()
            occ_modes = generator.generate(
                self.structure,
                self.symmetry_info["wyckoff_sites"],
                target,
                self._scope_species("occupational"),
            )
            for om in occ_modes:
                self.mode_occupancies[om.label] = {
                    "mode": om,
                    "pattern": om.pattern,
                    "basis": om.basis_vectors,
                    "validated": om.validated,
                    "note": om.note,
                }
        return modes

    # ================================================================
    # 阶段四-五：生成畸变结构
    # ================================================================

    def generate_distortion(self, irrep_label: str | None = None,
                            amplitude: float | None = None,
                            supercell: list | None = None) -> Structure:
        """
        步骤8-9：生成畸变结构

        Args:
            irrep_label: 不可约表示标号，为 None 则用第一个模式
            amplitude: 畸变幅度（位移向量最大分量按该幅度缩放）
            supercell: 超胞规格；None 时使用所选子群的基矢（3x3 矩阵）

        Returns:
            Structure: 畸变后的结构
        """
        if not self.mode_displacements and not self.mode_occupancies:
            raise RuntimeError(t("err.select_path_first"))

        if irrep_label is None:
            irrep_label = next(iter({**self.mode_displacements, **self.mode_occupancies}))

        if supercell is None and self.phase_path is not None:
            # 默认使用子群超胞基矢
            supercell = self.phase_path.supercell_basis()

        # occupational 模式：占据率调制（+1 类全占据，-1 类 1-amplitude）
        if irrep_label in self.mode_occupancies:
            entry = self.mode_occupancies[irrep_label]
            basis = supercell or entry["basis"]
            if not np.allclose(np.asarray(basis, dtype=float),
                               np.asarray(entry["basis"], dtype=float), atol=1e-8):
                raise ValueError(
                    "occupational 模式必须使用与生成时一致的子群超胞基矢"
                )
            amp = (self._dist_engine.default_amplitude
                   if amplitude is None else amplitude)
            self.distorted_structure = self._dist_engine.generate_modes(
                self.structure, basis,
                parent_displacements=None,
                occupancy_patterns=[(entry["pattern"], amp)],
            )
            print(t("distortion.generated", irrep=irrep_label, amp=amp,
                    n1=len(self.structure), n2=len(self.distorted_structure), r=1))
            fname = f"distorted_{irrep_label}"
            paths = self._exporter.auto_export(self.distorted_structure, fname,
                                               formats=["cif"])
            if paths:
                print(t("export.default", path=paths[0]))
            return self.distorted_structure

        if irrep_label not in self.mode_displacements:
            raise ValueError(t("mode.invalid", label=irrep_label))

        disp = self.mode_displacements[irrep_label]["displacements"]
        k_vector = (self.phase_path.k_vector if self.phase_path is not None
                    else None)
        self.distorted_structure = self._dist_engine.generate_single_mode(
            self.structure, disp, amplitude, supercell, k_vector=k_vector
        )

        n_ratio = len(self.distorted_structure) / len(self.structure)
        print(t("distortion.generated", irrep=irrep_label, amp=amplitude,
                n1=len(self.structure), n2=len(self.distorted_structure), r=n_ratio))

        # 默认导出畸变后的 CIF 文件
        fname = f"distorted_{irrep_label}" if irrep_label else "distorted"
        if amplitude is not None:
            amp_str = str(amplitude).replace(".", "p")
            fname = f"{fname}_a{amp_str}"
        paths = self._exporter.auto_export(self.distorted_structure, fname, formats=["cif"])
        if paths:
            print(t("export.default", path=paths[0]))
        return self.distorted_structure

    def generate_mixed_distortion(self, contributions: dict[str, float],
                                  supercell: list | None = None) -> Structure:
        """生成多模式混合畸变（可同时包含位移模式与 occupational 占据率模式）

        Args:
            contributions: {irrep_label: amplitude} 各模式贡献
            supercell: 超胞规格；None 时使用所选子群的基矢
        """
        if supercell is None and self.phase_path is not None:
            supercell = self.phase_path.supercell_basis()

        total_disp: np.ndarray | None = None
        occ_patterns: list[tuple[np.ndarray, float]] = []
        for label, amp in contributions.items():
            if label in self.mode_occupancies:
                entry = self.mode_occupancies[label]
                if supercell is None:
                    supercell = entry["basis"]
                if not np.allclose(np.asarray(supercell, dtype=float),
                                   np.asarray(entry["basis"], dtype=float), atol=1e-8):
                    raise ValueError(
                        "occupational 模式必须使用与生成时一致的子群超胞基矢"
                    )
                occ_patterns.append((entry["pattern"], float(amp)))
            elif label in self.mode_displacements:
                contribution = float(amp) * np.asarray(
                    self.mode_displacements[label]["displacements"], dtype=float
                )
                total_disp = contribution if total_disp is None else total_disp + contribution

        if total_disp is None and not occ_patterns:
            raise ValueError("未提供任何有效的模式贡献（位移或占据率）")

        self.distorted_structure = self._dist_engine.generate_modes(
            self.structure, supercell,
            parent_displacements=total_disp,
            occupancy_patterns=occ_patterns or None,
            k_vector=(self.phase_path.k_vector
                      if self.phase_path is not None else None),
        )
        # 默认导出混合畸变为 CIF
        label = "mixed"
        keys = "+".join(sorted(contributions.keys()))
        if keys:
            label = f"mixed_{keys}"
        paths = self._exporter.auto_export(self.distorted_structure, label, formats=["cif"])
        if paths:
            print(t("export.default", path=paths[0]))
        return self.distorted_structure

    # ================================================================
    # 阶段六：导出
    # ================================================================

    def export(self, filename: str, formats: list | None = None) -> list:
        """
        步骤11：导出畸变结构

        Args:
            filename: 文件名（不含后缀）
            formats: 导出格式列表，默认 ["cif"]

        Returns:
            list of Path: 导出文件路径
        """
        if self.distorted_structure is None:
            raise RuntimeError(t("err.generate_first"))

        paths = self._exporter.auto_export(
            self.distorted_structure, filename, formats
        )
        print(t("export.done", n=len(paths)))
        for p in paths:
            print(f"  {p}")
        return paths

    # ================================================================
    # 畴变体
    # ================================================================

    def generate_domains(self) -> list:
        """生成所有畴变体描述（官网 Domains 输出）。

        畴总数 = 子群在母相中的指数；需要先选择路径（select_path / Method 2）。
        """
        if self.phase_path is None or not self.subgroups:
            raise RuntimeError(t("err.domains_need_path"))

        target = next(
            (sg for sg in self.subgroups
             if sg.index == self.phase_path.subgroup_index), None
        )
        if target is None:
            raise RuntimeError(t("err.domains_not_in_list"))

        domains = self._domain_gen.generate_domains(self.phase_path, target)
        print(t("domains.found", n=len(domains)))
        return domains

    # ================================================================
    # ISODISTORT Search Method 1-4
    # ================================================================

    def search_method_1(self,
                        distortion_types: str | list[str] | None = None,
                        crystal_system: str | None = None,
                        subgroup_space_group: int | None = None,
                        direct_sublattice: list[int] | None = None,
                        lattice: list[list[float]] | None = None,
                        maximal_subgroup_only: bool = False):
        """
        Method 1: Search over all special k points.

        支持多条件同时过滤（逻辑 AND，与官网一致）：
        - lattice：官网 Conventional lattice / Primitive lattice 下拉所选
          子格（3x3 矩阵，惯用坐标；Primitive 选项请先用
          lattice_in_conventional_frame 换算）
        - direct_sublattice：旧版直接子格 (a,b,c)（对角阵，兼容保留）

        枚举结果按会话缓存，重复调用秒回。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        query = Method1Query(
            distortion_types=distortion_types,
            crystal_system=crystal_system,
            subgroup_space_group=subgroup_space_group,
            direct_sublattice=direct_sublattice,
            lattice=lattice,
            maximal_subgroup_only=maximal_subgroup_only,
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_1_search(
            parent_sg, query, subgroups=self._ensure_special_subgroups()
        )

        # 记录过滤后的候选，供 Method 2 使用
        self.subgroups = [item.subgroup for item in result]
        print(t("method1.result", n=len(result)))
        return result

    def search_method_2(self,
                        subgroup_idx: int,
                        distortion_type: str | list[str] = "displacive",
                        k_point_label: str | None = None,
                        k_point_coordinates: list | None = None,
                        k_parameters: dict[str, str | int | float] | None = None,
                        number_of_independent_modulations: int = 0,
                        number_of_superposed_irs: int = 1,
                        specified_opd: str | None = None):
        """
        Method 2: General method - search over specific k points.

        在 Method 1 候选（或 list_subgroups 枚举）中按序号选择子群，
        通过真实 iso（DISPLAY BUSH）计算其畸变模式基矢；
        按 self.distortion_scope 限定物种作用域（displacive/occupational 等），
        occupational 模式由本地生成器产生（存入 self.mode_occupancies）。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        types = normalize_distortion_types(distortion_type)
        if not self.subgroups:
            self.list_subgroups(distortion_type=distortion_type)

        query = Method2Query(
            subgroup_idx=subgroup_idx,
            distortion_type=types,
            k_point_label=k_point_label,
            k_point_coordinates=k_point_coordinates,
            k_parameters=k_parameters or {},
            number_of_independent_modulations=number_of_independent_modulations,
            number_of_superposed_irs=number_of_superposed_irs,
            specified_opd=specified_opd,
        )

        parent_sg = self.symmetry_info["space_group_number"]
        # 按作用域限制 BUSH 的 Wyckoff 位置（避免重复计算）
        scoped_letters = self._letters_for_species(self._union_scope_species(types))
        result = self._search.method_2_search(
            parent_sg, self.subgroups, query,
            wyckoff_letters=scoped_letters,
        )

        # 记录路径与模式，供 Distortion Page 使用
        self.phase_path = PhasePath.from_subgroup(
            parent_sg, result.subgroup, types
        )
        # 解析 k 点坐标（Bloch 相位调制用；仅特殊 k 点可直接求值）
        self.phase_path.k_vector = self._resolve_k_vector(
            result.subgroup.k_point_label)
        self.phase_path.validate()
        self.distortion_modes = self._compute_scoped_modes(
            parent_sg, result.subgroup, types, raw_modes=result.modes
        )
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )
        print(t("method2.result", idx=subgroup_idx,
                n=len(result.modes) + len(self.mode_occupancies)))
        return result

    def search_method_3(self,
                        distortion_types: str | list[str] | None = None,
                        point_group: str | None = None,
                        space_group_type: int | None = None,
                        supercell_basis: list[list[str | int | float]] | None = None,
                        direct_sublattice_centering: str | None = None,
                        lattice_type: str = "direct"):
        """
        Method 3: Search over arbitrary k points for a specified point group and supercell.

        若 point_group 与 space_group_type 同时提供，按官网规则优先采用
        space_group_type。lattice_type 为官网 radio（direct/reciprocal）；
        本地引擎暂不支持 reciprocal（倒易超格）模式，会给出明确错误。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        if lattice_type != "direct":
            raise ValueError(
                "本地引擎暂不支持 reciprocal（倒易空间超格）模式，"
                "请使用 direct（实空间子格）。官网该选项的完整实现在后续版本中支持。"
            )

        query = Method3Query(
            distortion_types=distortion_types,
            point_group=point_group,
            space_group_type=space_group_type,
            supercell_basis=supercell_basis,
            direct_sublattice_centering=direct_sublattice_centering,
            lattice_type=lattice_type,
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_3_search(parent_sg, query)

        # 记录过滤后的候选，供 Method 2（search_method_2）使用（与 Method 1 一致）
        self.subgroups = [item.subgroup for item in result]
        print(t("method3.result", n=len(result)))
        return result

    def search_method_4(self,
                        distorted_cif_path: str | Path,
                        atom_matching_method: str = "nearest-site",
                        robust_distance_threshold: float = 0.25,
                        provided_origin_shift: list[float] | None = None):
        """
        Method 4: Mode decomposition of a distorted structure.

        要求已经通过 Method 2 或 select_path 计算出可用模式。
        """
        if self.structure is None:
            raise RuntimeError("请先加载母相结构 (load_structure)")
        if not self.mode_displacements:
            raise RuntimeError("请先通过 select_path 或 search_method_2 计算模式")

        distorted_structure = read_cif(distorted_cif_path)
        mode_disp = {k: v["displacements"] for k, v in self.mode_displacements.items()}
        query = Method4Query(
            atom_matching_method=atom_matching_method,
            robust_distance_threshold=robust_distance_threshold,
            provided_origin_shift=provided_origin_shift,
        )

        result = self._search.method_4_decompose(
            self.structure,
            distorted_structure,
            mode_disp,
            query,
        )

        print(t("method4.result", n=len(result.amplitudes), rms=result.rms_residual))
        return result
