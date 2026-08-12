"""
核心 Python API - 一站式 ISODISTORT 使用入口

对应全阶段：封装 7 大阶段 12 个步骤，提供面向用户的高级接口
实现方式：❌ 自研（上层业务编排）

使用示例：
    from isodistort.api import IsoDistort

    iso = IsoDistort()
    # 1. 从 CIF 加载结构
    iso.load_structure("input.cif")
    # 2. 列出所有子群
    subgroups = iso.list_subgroups()
    # 3. 选择相变路径
    iso.select_path(subgroup_idx=3, distortion_type="displacement")
    # 4. 生成畸变结构
    distorted = iso.generate_distortion(amplitude=0.1)
    # 5. 导出
    iso.export("output", formats=["cif", "poscar"])
"""
from pathlib import Path
from typing import List, Dict, Optional

from pymatgen.core import Structure

from ..backend import FindsymWrapper, IsoWrapper, SubgroupInfo, DistortionMode
from ..structure import (
    read_cif, SymmetryValidator, SiteMapper,
)
from ..distortion import (
    PhasePath, DistortionMapper, DistortionEngine, DomainGenerator,
    IsoSearchEngine,
    Method1Query,
    Method2Query,
    Method3Query,
    Method4Query,
)
from ..io import StructureExporter, ResultSerializer
from ..utils import get_config


class IsoDistort:
    """
    ISODISTORT 主入口类

    封装完整工作流：
        加载结构 → 识别对称 → 枚举子群 → 选择路径 →
        计算畸变模式 → 生成畸变结构 → 导出/可视化
    """

    def __init__(self):
        """Relative path: isocore/api/core_api.py"""
        
        self.cfg = get_config()

        # 底层封装
        self._findsym = FindsymWrapper()
        self._iso = IsoWrapper()

        # 业务层
        self._sym_val = SymmetryValidator()
        self._site_mapper = SiteMapper()
        self._dist_mapper = DistortionMapper(self._site_mapper)
        self._dist_engine = DistortionEngine(self._dist_mapper)
        self._domain_gen = DomainGenerator(self._iso)
        self._search = IsoSearchEngine(self._iso)

        # 输出层
        self._exporter = StructureExporter()
        self._serializer = ResultSerializer()

        # 状态
        self.structure: Optional[Structure] = None
        self.symmetry_info: Optional[Dict] = None
        self.subgroups: List[SubgroupInfo] = []
        self.phase_path: Optional[PhasePath] = None
        self.distortion_modes: List[DistortionMode] = []
        self.mode_displacements: Dict = {}
        self.distorted_structure: Optional[Structure] = None

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
        
        Relative path: isocore/api/core_api.py"""

        self.structure = read_cif(cif_path)
        self.symmetry_info = self._sym_val.validate(self.structure)

        sg_num = self.symmetry_info["space_group_number"]
        sg_sym = self.symmetry_info["space_group_symbol"]
        n_atoms = len(self.structure)

        print(f"[加载完成] 空间群 #{sg_num} ({sg_sym}), {n_atoms} 个原子")
        return self.structure

    def set_structure(self, structure: Structure):
        """直接设置 Structure 对象

        Relative path: isocore/api/core_api.py"""

        self.structure = structure
        self.symmetry_info = self._sym_val.validate(self.structure)

    # ================================================================
    # 阶段二：子群枚举
    # ================================================================

    def list_subgroups(self, distortion_type: str | list[str] | None = None
                    ) -> List[SubgroupInfo]:
        """
        步骤4：枚举母相的所有各向同性子群。

        默认同时包含 displacement 与 strain 两种模式。
        """

        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        sg_num = self.symmetry_info["space_group_number"]
        self.subgroups = self._iso.list_subgroups(sg_num, distortion_type)

        print(f"[子群枚举] 共找到 {len(self.subgroups)} 个各向同性子群")
        for sg in self.subgroups[:10]:  # 只打印前10个
            print(f"  #{sg.index:2d}  空间群 #{sg.space_group_number}  "
                f"{sg.space_group_symbol}")
        if len(self.subgroups) > 10:
            print(f"  ... 还有 {len(self.subgroups) - 10} 个")

        return self.subgroups

    # ================================================================
    # 阶段三：畸变模式计算
    # ================================================================

    def select_path(self, subgroup_idx: int,
                    distortion_type: str | None = None) -> PhasePath:
        """
        步骤5-6：选择相变路径，计算畸变模式。

        若不指定 distortion_type，则按默认 displacement + strain 取第一个
        能命中该子群的模式。
        """

        if not self.subgroups:
            self.list_subgroups(distortion_type)

        if distortion_type is None:
            distortion_type = "displacement"

        # 找到对应子群
        target = None
        for sg in self.subgroups:
            if sg.index == subgroup_idx:
                target = sg
                break
        if target is None:
            raise ValueError(f"子群序号 {subgroup_idx} 不存在")

        self.phase_path = PhasePath.from_subgroup(
            self.symmetry_info["space_group_number"],
            target,
            distortion_type,
        )
        self.phase_path.validate()

        print(f"[路径选择] {self.phase_path.describe()}")

        # 计算畸变模式
        sg_num = self.symmetry_info["space_group_number"]
        self.distortion_modes = self._iso.calc_distortion_modes(
            sg_num, subgroup_idx, distortion_type
        )

        print(f"[模式计算] 共 {len(self.distortion_modes)} 个畸变模式")
        for m in self.distortion_modes:
            print(f"  {m.irrep_label:8s}  {m.dimension} 维")

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

    def generate_distortion(self, irrep_label: str = None,
                            amplitude: float = None,
                            supercell: list = None) -> Structure:
        """
        步骤8-9：生成畸变结构

        Args:
            irrep_label: 不可约表示标号，为 None 则用第一个模式
            amplitude: 畸变幅度
            supercell: 超胞大小

        Returns:
            Structure: 畸变后的结构
        
        Relative path: isocore/api/core_api.py"""

        if not self.mode_displacements:
            raise RuntimeError("请先选择相变路径 (select_path)")

        if irrep_label is None:
            irrep_label = list(self.mode_displacements.keys())[0]

        if irrep_label not in self.mode_displacements:
            raise ValueError(f"模式 {irrep_label} 不存在")

        disp = self.mode_displacements[irrep_label]["displacements"]
        self.distorted_structure = self._dist_engine.generate_single_mode(
            self.structure, disp, amplitude, supercell
        )

        vol_change = (self.distorted_structure.volume - self.structure.volume) \
                     / self.structure.volume * 100
        print(f"[畸变生成] 模式 {irrep_label}, 幅度 {amplitude}, "
            f"体积变化 {vol_change:+.4f}%")

        # 默认导出畸变后的 CIF 文件
        fname = f"distorted_{irrep_label}" if irrep_label else "distorted"
        if amplitude is not None:
            try:
                amp_str = str(amplitude).replace(".", "p")
                fname = f"{fname}_a{amp_str}"
            except Exception:
                pass
        paths = self._exporter.auto_export(self.distorted_structure, fname, formats=["cif"])
        if paths:
            print(f"[默认导出] 已生成 CIF: {paths[0]}")
        return self.distorted_structure

    def generate_mixed_distortion(self, contributions: Dict[str, float],
                                supercell: list = None) -> Structure:
        """生成多模式混合畸变

        Relative path: isocore/api/core_api.py"""

        all_disp = {k: v["displacements"] for k, v in self.mode_displacements.items()}
        self.distorted_structure = self._dist_engine.generate_mixed_mode(
            self.structure, contributions, all_disp, supercell
        )
        # 默认导出混合畸变为 CIF
        label = "mixed"
        try:
            keys = "+".join(sorted(contributions.keys()))
            label = f"mixed_{keys}"
        except Exception:
            pass
        paths = self._exporter.auto_export(self.distorted_structure, label, formats=["cif"])
        if paths:
            print(f"[默认导出] 已生成 CIF: {paths[0]}")
        return self.distorted_structure

    # ================================================================
    # 阶段六：导出
    # ================================================================

    def export(self, filename: str, formats: list = None) -> list:
        """
        步骤11：导出畸变结构

        Args:
            filename: 文件名（不含后缀）
            formats: 导出格式列表，默认 ["cif"]

        Returns:
            list of Path: 导出文件路径
        
        Relative path: isocore/api/core_api.py"""

        if self.distorted_structure is None:
            raise RuntimeError("请先生成畸变结构 (generate_distortion)")

        paths = self._exporter.auto_export(
            self.distorted_structure, filename, formats
        )
        print(f"[导出完成] 共 {len(paths)} 个文件:")
        for p in paths:
            print(f"  {p}")
        return paths

    # 阶段七：可视化（已删除 vis 依赖；可使用导出的 CIF 用外部程序查看）

    # ================================================================
    # 畴变体
    # ================================================================

    def generate_domains(self) -> list:
        """生成所有畴变体结构

        Relative path: isocore/api/core_api.py"""

        if self.distorted_structure is None or self.phase_path is None:
            raise RuntimeError("请先生成畸变结构并选择路径")

        domains = self._domain_gen.generate_domains(
            self.distorted_structure, self.phase_path
        )
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

        支持多条件同时过滤（逻辑 AND）。
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
        print(f"[Method 1] 过滤后得到 {len(result)} 条候选")
        return result

    def search_method_2(self,
                        subgroup_idx: int,
                        distortion_type: str = "displacement",
                        k_point_label: str | None = None,
                        k_point_coordinates: list | None = None,
                        k_parameters: Dict[str, str | int | float] | None = None,
                        number_of_independent_modulations: int = 0,
                        number_of_superposed_irs: int = 1,
                        specified_opd: str | None = None):
        """
        Method 2: General method - search over specific k points.

        当前本地实现会保留 k 点与 OPD 参数，并基于所选子群计算模式列表。
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
        result = self._search.method_2_search(parent_sg, self.subgroups, query)

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

        若 point_group 与 space_group_type 同时提供，将按官网规则优先采用
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
