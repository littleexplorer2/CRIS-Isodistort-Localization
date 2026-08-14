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
from pathlib import Path

from pymatgen.core import Structure

from ..backend import (
    DistortionMode,
    FindsymWrapper,
    IsoWrapper,
    SubgroupInfo,
)
from ..distortion import (
    DistortionEngine,
    DistortionMapper,
    DomainGenerator,
    IsoSearchEngine,
    Method1Query,
    Method2Query,
    Method3Query,
    Method4Query,
    PhasePath,
)
from ..io import StructureExporter
from ..structure import (
    SymmetryValidator,
    read_cif,
)
from ..utils import get_config


class IsoDistort:
    """
    ISODISTORT 主入口类

    封装完整工作流：
        加载结构 → 识别对称 → 枚举子群（Method 1）→ 选择路径（Method 2）
        → 计算畸变模式 → 生成畸变结构 → 导出/畴
    """

    def __init__(self) -> None:
        self.cfg = get_config()

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
        self.distorted_structure: Structure | None = None

    # ================================================================
    # 阶段一：结构输入与对称识别
    # ================================================================

    def load_structure(self, cif_path: str | Path) -> Structure:
        """
        步骤1-3：加载 CIF，识别对称性

        Args:
            cif_path: CIF 文件路径

        Returns:
            Structure: 加载后的结构对象
        """
        self.structure = read_cif(cif_path)
        self.symmetry_info = self._sym_val.validate(self.structure)
        self._reset_derived_state()

        sg_num = self.symmetry_info["space_group_number"]
        sg_sym = self.symmetry_info["space_group_symbol"]
        n_atoms = len(self.structure)

        print(f"[加载完成] 空间群 #{sg_num} ({sg_sym}), {n_atoms} 个原子")
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
        self.distorted_structure = None

    def _wyckoff_letters(self) -> list[str]:
        """当前结构各 Wyckoff 位点字母（用于 iso 模式计算）。"""
        if not self.symmetry_info:
            return []
        return [s["wyckoff_letter"] for s in self.symmetry_info["wyckoff_sites"]]

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

    def list_subgroups(self, distortion_type: str | list[str] | None = None
                       ) -> list[SubgroupInfo]:
        """
        步骤4：枚举母相全部特殊 k 点的各向同性子群（Method 1 数据源）。

        Args:
            distortion_type: 畸变类型（保留参数；类型过滤在模式计算阶段执行）

        Returns:
            List[SubgroupInfo]
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        sg_num = self.symmetry_info["space_group_number"]
        self.subgroups = self._iso.enumerate_all_special_subgroups(
            sg_num, distortion_type
        )

        print(f"[子群枚举] 共找到 {len(self.subgroups)} 个各向同性子群")
        for sg in self.subgroups[:10]:
            print(f"  {sg.describe()}")
        if len(self.subgroups) > 10:
            print(f"  ... 还有 {len(self.subgroups) - 10} 个")
        return self.subgroups

    # ================================================================
    # 阶段三：路径选择与畸变模式计算（Method 2）
    # ================================================================

    def select_path(self, subgroup_idx: int,
                    distortion_type: str | None = None) -> PhasePath:
        """
        步骤5-6：选择相变路径，计算畸变模式（Mode Basis）。

        Args:
            subgroup_idx: 子群序号（来自 Method 1 候选）
            distortion_type: 畸变类型（默认 displacement）

        Returns:
            PhasePath
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        if not self.subgroups:
            self.list_subgroups(distortion_type)

        if distortion_type is None:
            distortion_type = "displacement"

        target = next((sg for sg in self.subgroups if sg.index == subgroup_idx), None)
        if target is None:
            raise ValueError(f"子群序号 {subgroup_idx} 不存在")

        self.phase_path = PhasePath.from_subgroup(
            self.symmetry_info["space_group_number"],
            target,
            distortion_type,
        )
        self.phase_path.validate()

        print(f"[路径选择] {self.phase_path.describe()}")

        # 计算畸变模式基矢（DISPLAY BUSH）
        self.distortion_modes = self._iso.calc_distortion_modes(
            self.symmetry_info["space_group_number"],
            target,
            wyckoff_letters=self._wyckoff_letters(),
        )

        print(f"[模式计算] 共 {len(self.distortion_modes)} 个畸变模式")
        for m in self.distortion_modes:
            n_sites = len({b.wyckoff_letter for b in m.bush_modes})
            print(f"  {m.irrep_label:8s}  {m.opd_symbol:6s} 涉及位点 {n_sites} 个")

        # 映射到原子位移
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )

        return self.phase_path

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
        if not self.mode_displacements:
            raise RuntimeError("请先选择相变路径 (select_path) 或执行 Method 2")

        if irrep_label is None:
            irrep_label = next(iter(self.mode_displacements))

        if irrep_label not in self.mode_displacements:
            raise ValueError(f"模式 {irrep_label} 不存在")

        disp = self.mode_displacements[irrep_label]["displacements"]
        if supercell is None and self.phase_path is not None:
            # 默认使用子群超胞基矢
            supercell = self.phase_path.supercell_basis()

        self.distorted_structure = self._dist_engine.generate_single_mode(
            self.structure, disp, amplitude, supercell
        )

        n_ratio = len(self.distorted_structure) / len(self.structure)
        print(f"[畸变生成] 模式 {irrep_label}, 幅度 {amplitude}, "
              f"原子数 {len(self.structure)} -> {len(self.distorted_structure)} "
              f"(超胞倍数 {n_ratio:g})")

        # 默认导出畸变后的 CIF 文件
        fname = f"distorted_{irrep_label}" if irrep_label else "distorted"
        if amplitude is not None:
            amp_str = str(amplitude).replace(".", "p")
            fname = f"{fname}_a{amp_str}"
        paths = self._exporter.auto_export(self.distorted_structure, fname, formats=["cif"])
        if paths:
            print(f"[默认导出] 已生成 CIF: {paths[0]}")
        return self.distorted_structure

    def generate_mixed_distortion(self, contributions: dict[str, float],
                                  supercell: list | None = None) -> Structure:
        """生成多模式混合畸变

        Args:
            contributions: {irrep_label: amplitude} 各模式贡献
            supercell: 超胞规格；None 时使用所选子群的基矢
        """
        all_disp = {k: v["displacements"] for k, v in self.mode_displacements.items()}
        if supercell is None and self.phase_path is not None:
            supercell = self.phase_path.supercell_basis()
        self.distorted_structure = self._dist_engine.generate_mixed_mode(
            self.structure, contributions, all_disp, supercell
        )
        # 默认导出混合畸变为 CIF
        label = "mixed"
        keys = "+".join(sorted(contributions.keys()))
        if keys:
            label = f"mixed_{keys}"
        paths = self._exporter.auto_export(self.distorted_structure, label, formats=["cif"])
        if paths:
            print(f"[默认导出] 已生成 CIF: {paths[0]}")
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
            raise RuntimeError("请先生成畸变结构 (generate_distortion)")

        paths = self._exporter.auto_export(
            self.distorted_structure, filename, formats
        )
        print(f"[导出完成] 共 {len(paths)} 个文件:")
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
            raise RuntimeError("请先选择相变路径 (select_path)")

        target = next(
            (sg for sg in self.subgroups
             if sg.index == self.phase_path.subgroup_index), None
        )
        if target is None:
            raise RuntimeError("当前路径对应的子群不在候选列表中")

        domains = self._domain_gen.generate_domains(self.phase_path, target)
        print(f"[畴生成] 共 {len(domains)} 个畴变体")
        return domains

    # ================================================================
    # ISODISTORT Search Method 1-4
    # ================================================================

    def search_method_1(self,
                        distortion_types: str | list[str] | None = None,
                        crystal_system: str | None = None,
                        subgroup_space_group: int | None = None,
                        direct_sublattice: list[int] | None = None,
                        maximal_subgroup_only: bool = False):
        """
        Method 1: Search over all special k points.

        支持多条件同时过滤（逻辑 AND，与官网一致）。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        query = Method1Query(
            distortion_types=distortion_types,
            crystal_system=crystal_system,
            subgroup_space_group=subgroup_space_group,
            direct_sublattice=direct_sublattice,
            maximal_subgroup_only=maximal_subgroup_only,
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_1_search(parent_sg, query)

        # 记录过滤后的候选，供 Method 2 使用
        self.subgroups = [item.subgroup for item in result]
        print(f"[Method 1] 过滤后得到 {len(result)} 条候选")
        return result

    def search_method_2(self,
                        subgroup_idx: int,
                        distortion_type: str = "displacement",
                        k_point_label: str | None = None,
                        k_point_coordinates: list | None = None,
                        k_parameters: dict[str, str | int | float] | None = None,
                        number_of_independent_modulations: int = 0,
                        number_of_superposed_irs: int = 1,
                        specified_opd: str | None = None):
        """
        Method 2: General method - search over specific k points.

        在 Method 1 候选（或 list_subgroups 枚举）中按序号选择子群，
        通过真实 iso（DISPLAY BUSH）计算其畸变模式基矢。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        if not self.subgroups:
            self.list_subgroups(distortion_type=distortion_type)

        query = Method2Query(
            subgroup_idx=subgroup_idx,
            distortion_type=distortion_type,
            k_point_label=k_point_label,
            k_point_coordinates=k_point_coordinates,
            k_parameters=k_parameters or {},
            number_of_independent_modulations=number_of_independent_modulations,
            number_of_superposed_irs=number_of_superposed_irs,
            specified_opd=specified_opd,
        )

        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_2_search(
            parent_sg, self.subgroups, query,
            wyckoff_letters=self._wyckoff_letters(),
        )

        # 记录路径与模式，供 Distortion Page 使用
        self.phase_path = PhasePath.from_subgroup(
            parent_sg, result.subgroup, distortion_type
        )
        self.phase_path.validate()
        self.distortion_modes = result.modes
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )
        print(f"[Method 2] 子群 #{subgroup_idx} 共得到 {len(result.modes)} 个模式")
        return result

    def search_method_3(self,
                        distortion_types: str | list[str] | None = None,
                        point_group: str | None = None,
                        space_group_type: int | None = None,
                        supercell_basis: list[list[str | int | float]] | None = None,
                        direct_sublattice_centering: str | None = None):
        """
        Method 3: Search over arbitrary k points for a specified point group and supercell.

        若 point_group 与 space_group_type 同时提供，按官网规则优先采用
        space_group_type。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        query = Method3Query(
            distortion_types=distortion_types,
            point_group=point_group,
            space_group_type=space_group_type,
            supercell_basis=supercell_basis,
            direct_sublattice_centering=direct_sublattice_centering,
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_3_search(parent_sg, query)
        print(f"[Method 3] 约束搜索得到 {len(result)} 条候选")
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

        print(
            f"[Method 4] 分解完成，{len(result.amplitudes)} 个模式，"
            f"RMS 残差 = {result.rms_residual:.6f}"
        )
        return result
