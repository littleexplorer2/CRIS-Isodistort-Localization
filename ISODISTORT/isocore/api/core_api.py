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
    iso.export_subgroups("out_batch", formats=["cif", "topas"])
"""
import html
import re
import threading
from fractions import Fraction
from pathlib import Path

import numpy as np
import spglib
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ..backend import (
    DistortionMode,
    IsoWrapper,
    SubgroupInfo,
)
from ..backend.isotropy_cache import (
    IsotropyCacheEntry,
    delete_isotropy_cache,
    list_isotropy_cache,
)
from ..backend.smodes_wrapper import SmodesWrapper
from ..data.kpoints_official import (
    KPOINT_OFFICIAL,
    official_kparams_to_iso,
    official_special_k_coords,
)
from ..distortion import (
    DEFAULT_DISTORTION_TYPES,
    DISTORTION_TYPES,
    DistortionEngine,
    DistortionMapper,
    DomainGenerator,
    FiniteGroup,
    IsoSearchEngine,
    Method1Query,
    Method1ResultItem,
    Method2Query,
    Method3Query,
    Method4Query,
    OccupationalModeGenerator,
    PhasePath,
    RationalRepresentation,
    analyze_fixed_space,
    embedding_from_identity,
    normalize_distortion_types,
    parent_affine_group,
    symmetric_square_representation,
)
from ..distortion.search_methods import _sg_to_crystal_system
from ..i18n import t
from ..io import (
    StructureExporter,
    SubgroupExportSpec,
    build_export_zip,
    parse_export_formats,
    subgroup_label,
    unique_folder_name,
    write_subgroup_files,
)
from ..structure import (
    SymmetryValidator,
    build_supercell,
    read_cif,
    read_structure,
)
from ..utils import IsodistortError, get_config
from ..utils.lattice import as_fraction, identity_matrix, multiply
from ..utils.opd_format import _centering_letter, format_k_active
from ..utils.parent_header import parent_wyckoff_display
from ..utils.schoenflies import hm_symbol, schoenflies_symbol
from ..utils.text_parser import parse_basis_token, parse_fraction


class IsoDistort:
    """
    ISODISTORT 主入口类

    封装完整工作流：
        加载结构 → 识别对称 → 枚举子群（Method 1）→ 选择路径（Method 2）
        → 计算畸变模式 → 生成畸变结构 → 导出/畴

    UI strings are English only. The ``language`` argument is ignored and kept
    so existing ``IsoDistort(language="en")`` call sites still construct.
    """

    def __init__(self, language: str | None = None) -> None:
        self.cfg = get_config()
        _ = language

        # 底层封装
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
        self.structure_path: Path | None = None
        self.subgroups: list[SubgroupInfo] = []
        # Affine-only Method-3 diagnostics must survive copying/rehydration of
        # ``SubgroupInfo`` objects.  Python ``id(...)`` is a process-local
        # memory address and can be reused, so track the deterministic path
        # exact content-addressed Seitz embedding ID (with an exact path-key
        # compatibility fallback) instead.
        self._unresolved_method3_embedding_keys: set[tuple] = set()
        # 当前已计算模式所对应的完整子群对象。不能只保存 ``index``：
        # Method 1/2/3 的候选池都会从 0 编号，同号不代表同一个子群。
        self._selected_subgroup: SubgroupInfo | None = None
        self.phase_path: PhasePath | None = None
        self.distortion_modes: list[DistortionMode] = []
        self.mode_displacements: dict = {}
        self.mode_occupancies: dict = {}          # occupational 模式（占据率调制）
        self.mode_displacements_sc: dict = {}     # parametric complete modes on supercell
        self._mode_label_overrides: dict[str, str] = {}
        self.number_of_independent_modulations: int = 0
        self.distorted_structure: Structure | None = None

        # 畸变类型作用域（对齐官网 per-species 复选框）：type -> 物种列表（"*"=全部）
        self.distortion_scope: dict[str, list[str]] = {}
        self.distortion_types: list[str] = DEFAULT_DISTORTION_TYPES.copy()
        self._smodes = SmodesWrapper()
        self._special_subgroups_cache: list | None = None
        self._special_subgroups_lock = threading.Lock()
        self._conv_to_prim_cache: np.ndarray | None = None
        self._parent_rotations_cache: list[np.ndarray] | None = None
        self._strain_representation_cache = None
        self._lattice_standardization_cache: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}

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
        self.structure_path = path.resolve()
        self.symmetry_info = self._sym_val.validate(self.structure)
        self._reset_derived_state()
        # Parse and attach source-CIF atom-site labels/order at load time.
        # Web status rendering used to trigger this mutation incidentally,
        # which made later mode labels depend on whether a UI endpoint had
        # already been visited.  All entry points must start from one state.
        self.parent_wyckoff_display()

        sg_num = self.symmetry_info["space_group_number"]
        sg_sym = self.symmetry_info["space_group_symbol"]
        n_atoms = len(self.structure)

        print(t("load.done", sg=sg_num, sym=sg_sym, n=n_atoms))
        return self.structure

    def set_structure(self, structure: Structure) -> Structure:
        """直接设置 Structure 对象"""
        self.structure = structure
        self.structure_path = None
        self.symmetry_info = self._sym_val.validate(self.structure)
        self._reset_derived_state()
        return self.structure

    def parent_wyckoff_display(self) -> list[str]:
        """官网页头 Wyckoff 行：优先按母相 CIF 位点顺序与标签，否则用对称分析。"""
        if self.structure is None or not self.symmetry_info:
            return []
        return parent_wyckoff_display(
            self.structure,
            self.symmetry_info.get("wyckoff_sites") or [],
            self.structure_path,
        )

    def _reset_derived_state(self) -> None:
        """加载新结构后清空所有派生状态。"""
        self.subgroups = []
        self._unresolved_method3_embedding_keys = set()
        self._selected_subgroup = None
        self.phase_path = None
        self.distortion_modes = []
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.distorted_structure = None
        self._special_subgroups_cache = None
        self._conv_to_prim_cache = None
        self._parent_rotations_cache = None
        self._strain_representation_cache = None
        self._lattice_standardization_cache = {}

    def set_subgroup_candidates(self, subgroups: list[SubgroupInfo]) -> None:
        """Replace the default Method-2 candidate pool through the public API.

        Prefer passing ``candidates=...`` directly to :meth:`search_method_2`
        when a caller keeps several Method result tables at the same time.
        This setter exists for interactive flows that intentionally make one
        table the session default.
        """
        self.subgroups = list(subgroups)

    def clear_selected_modes(self) -> None:
        """Clear the selected path and all derived mode/distortion state."""
        self._selected_subgroup = None
        self.phase_path = None
        self.distortion_modes = []
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.mode_displacements_sc = {}
        self._mode_label_overrides = {}
        self.distorted_structure = None

    def list_isotropy_cache(self) -> list[IsotropyCacheEntry]:
        """Return generated ISOTROPY cache entries without exposing the backend."""
        return list_isotropy_cache(self._iso)

    def delete_isotropy_cache(self, names: list[str]) -> dict:
        """Delete named generated cache files through the public API."""
        return delete_isotropy_cache(self._iso, names)

    @staticmethod
    def _subgroup_identity(subgroup: SubgroupInfo | None) -> tuple | None:
        """Stable path identity; local display indices are deliberately excluded."""
        if subgroup is None:
            return None

        def _matrix(values) -> tuple:
            return tuple(
                tuple(round(float(value), 10) for value in row)
                for row in (values or [])
            )

        return (
            int(subgroup.parent_sg),
            int(subgroup.space_group_number),
            int(subgroup.subgroup_index),
            int(subgroup.size),
            str(subgroup.k_point_label or ""),
            tuple(str(value) for value in (subgroup.k_coordinates or [])),
            str(subgroup.irrep_label or ""),
            str(subgroup.opd_symbol or ""),
            str(subgroup.opd_dir_raw or ""),
            tuple(round(float(value), 10) for value in (subgroup.opd_vector or [])),
            tuple(str(value) for value in (subgroup.k_parameters or [])),
            _matrix(subgroup.basis_vectors),
            tuple(round(float(value), 10) for value in (subgroup.origin or [])),
        )

    @staticmethod
    def _method3_embedding_guard_key(subgroup: SubgroupInfo) -> tuple:
        """Stable identity for rejecting unresolved affine-only rows.

        Stage-A rows carry a content-addressed ID derived from their exact
        Seitz subgroup.  Compatibility/mocked rows may lack that field; their
        fallback key serializes every scientific path component as an exact
        Fraction (or an unmodified symbolic token), never as a display index
        or rounded float.  The key therefore survives copying and process
        rehydration without merging close but distinct embeddings.
        """
        embedding_id = str(
            getattr(subgroup, "_method3_embedding_id", "") or ""
        ).strip()
        if embedding_id:
            return ("exact-embedding-id", embedding_id)

        def scalar(value, *, periodic: bool = False) -> tuple:
            try:
                fraction = as_fraction(value)
            except (TypeError, ValueError):
                return ("symbol", str(value))
            if periodic:
                fraction %= 1
            return ("fraction", fraction.numerator, fraction.denominator)

        def matrix(values) -> tuple:
            return tuple(
                tuple(scalar(value) for value in row)
                for row in (values or [])
            )

        return (
            "exact-path-fallback-v1",
            int(subgroup.parent_sg),
            int(subgroup.space_group_number),
            int(subgroup.subgroup_index),
            int(subgroup.size),
            str(subgroup.space_group_symbol or ""),
            str(subgroup.k_point_label or ""),
            tuple(scalar(value) for value in (subgroup.k_coordinates or [])),
            str(subgroup.irrep_label or ""),
            str(subgroup.opd_symbol or ""),
            str(subgroup.opd_dir_raw or ""),
            tuple(scalar(value) for value in (subgroup.opd_vector or [])),
            tuple(scalar(value) for value in (subgroup.k_parameters or [])),
            matrix(subgroup.basis_vectors),
            tuple(scalar(value, periodic=True) for value in (subgroup.origin or [])),
        )

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

    def set_distortion_types(self, distortion_types) -> None:
        """设置当前考虑的畸变类型（对齐官网 Types 面板复选框）。"""
        self.distortion_types = normalize_distortion_types(distortion_types)

    def _iso_kpoint_raw(self, k_point_label: str):
        """iso 原始 k 点信息（未应用官网显示覆盖）。"""
        for kp in self._iso.list_k_points(self.symmetry_info["space_group_number"]):
            if kp.label == k_point_label.strip():
                return kp
        return None

    def _resolve_iso_kparams(self, k_point_label: str,
                             k_parameters: list | None) -> list | None:
        """官网 UI 参数 → iso KVALUE 参数；无参数或未收录时原样返回。"""
        if not k_parameters:
            return None
        iso_kp = self._iso_kpoint_raw(k_point_label)
        if iso_kp is None:
            return list(k_parameters)
        return official_kparams_to_iso(
            self.symmetry_info["space_group_number"],
            k_point_label,
            k_parameters,
            iso_kp,
        )

    def _is_parametric_kpoint(self, k_point_label: str) -> bool:
        iso_kp = self._iso_kpoint_raw(k_point_label)
        return bool(iso_kp and iso_kp.parameters)

    def _filter_subgroups_for_search(self,
                                     subgroups: list[SubgroupInfo],
                                     k_point_label: str,
                                     official_kparams: list | None) -> list[SubgroupInfo]:
        """按 Distortion Types + 物种作用域过滤子群（对齐官网 Search 阶段）。"""
        types = normalize_distortion_types(self.distortion_types)
        mode_types = [tp for tp in types if tp in ("displacive", "rotational", "occupational")]
        if not mode_types:
            return subgroups

        if self._is_parametric_kpoint(k_point_label):
            species = self._union_scope_species(types)
            active = self._smodes.active_irreps(
                self.structure,
                self.symmetry_info["space_group_number"],
                self.symmetry_info["wyckoff_sites"],
                k_point_label,
                self._resolve_iso_kparams(k_point_label, official_kparams),
                species_filter=species if species else None,
            )
            if active is not None:
                return [sg for sg in subgroups if sg.irrep_label in active]
            return subgroups

        # 特殊 k 点：用 BUSH 探测各子群是否在作用域 Wyckoff 上有位移模式
        if not any(tp in types for tp in ("displacive", "rotational")):
            return subgroups
        letters = self._letters_for_species(self._union_scope_species(types))
        if not letters:
            return []
        parent_sg = self.symmetry_info["space_group_number"]
        kept: list[SubgroupInfo] = []
        for sg in subgroups:
            try:
                modes = self._iso.calc_distortion_modes(parent_sg, sg, letters)
            except IsodistortError:
                continue
            if modes:
                kept.append(sg)
        return kept

    def _tag_official_kparams(self, subgroups: list[SubgroupInfo],
                              official_kparams: list | None) -> None:
        """子群对象上保留官网参数（供界面显示），并刷新 k 坐标 / k-active。"""
        if not official_kparams:
            return
        tagged = list(official_kparams)
        for sg in subgroups:
            sg.k_parameters = tagged
            coords = official_special_k_coords(
                int(sg.parent_sg or 0),
                sg.k_point_label or "",
                sg.k_coordinates or [],
                tagged,
            )
            if coords:
                sg.k_coordinates = coords
            sg.k_active_raw = format_k_active(
                sg.opd_dir_raw or "",
                sg.k_coordinates or ["0", "0", "0"],
                sg.parent_sg or None,
                None,
            )

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

        本地 iso 二进制固定采用国际标准取位（官网默认值），与实际计算行为一致，
        因此仅提供只读说明，不提供可交互的偏好面板（自定义取位无法被本地引擎生效）。
        """
        return ("monoclinic axes a(b)c, monoclinic cell choice 1, "
                "orthorhombic axes abc, origin choice 2, hexagonal axes, "
                "SSG standard setting")

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

    def _filter_method1_by_types(self, items, distortion_types) -> list:
        """Keep Method 1 rows that the official site would list for the selected types.

        Each requested physical representation is tested independently and the
        union is returned.  Displacive/rotational activity comes from smodes;
        occupational activity requires a generated pattern that independently
        validates as the requested subgroup; homogeneous strain uses the
        symmetric metric-tensor representation ``E -> R.T E R``.  An unknown
        backend result never masquerades as physical activity.
        """
        types = normalize_distortion_types(distortion_types or self.distortion_types)
        want_mag = "magnetic" in types
        want_disp = any(tp in types for tp in ("displacive", "rotational"))
        want_strain = "strain" in types
        want_occupational = "occupational" in types
        if not items:
            return items

        active_by_k: dict[str, set[str] | None] = {}
        if want_disp and self.structure is not None and self.symmetry_info:
            species = self._union_scope_species(types)
            for kp in {it.subgroup.k_point_label for it in items}:
                if not kp:
                    active_by_k[kp] = None
                    continue
                active_by_k[kp] = self._smodes.active_irreps(
                    self.structure,
                    self.symmetry_info["space_group_number"],
                    self.symmetry_info["wyckoff_sites"],
                    kp,
                    None,
                    species_filter=species if species else None,
                )

        if want_disp:
            unresolved = sorted(
                label for label, active in active_by_k.items() if active is None
            )
            if unresolved:
                raise RuntimeError(
                    "Could not establish displacive/rotational irrep activity "
                    "from smodes for k point(s): " + ", ".join(unresolved)
                )

        parent_sg_n = 0
        if self.symmetry_info:
            parent_sg_n = int(self.symmetry_info["space_group_number"])

        occupational_generator = (
            OccupationalModeGenerator() if want_occupational else None
        )
        occupational_scope = (
            self._scope_species("occupational") if want_occupational else set()
        )

        kept = []
        for item in items:
            ir = (item.subgroup.irrep_label or "").strip()
            if ir.startswith("m"):
                if want_mag:
                    kept.append(item)
                continue

            is_active = False
            if want_disp:
                active = active_by_k.get(item.subgroup.k_point_label)
                is_active = active is not None and ir in active
            if want_strain and self._keep_strain_only_irrep(item, parent_sg_n):
                is_active = True
            if (
                occupational_generator is not None
                and self.structure is not None
                and self.symmetry_info
                and occupational_scope
            ):
                try:
                    occupational_modes = occupational_generator.generate(
                        self.structure,
                        self.symmetry_info["wyckoff_sites"],
                        item.subgroup,
                        occupational_scope,
                    )
                except Exception:  # noqa: BLE001 - unsupported route is inactive
                    occupational_modes = []
                if any(bool(getattr(mode, "validated", False))
                       for mode in occupational_modes):
                    is_active = True
            if is_active:
                kept.append(item)
        return kept

    def _filter_method3_routes_by_types(self, items, distortion_types) -> list:
        """Filter each Method 3 embedding through its known active routes.

        Method 3 rows are affine subgroup embeddings, while the local search
        engine discovers them through one or more single-IR routes.  Activity
        is therefore a property of each route, not of only the representative
        stored in ``item.subgroup``.  Parameter values are part of the k-point
        identity and must be supplied (after the existing official-to-backend
        parameter conversion) to smodes when checking displacement activity.

        The filtering follows the capabilities already exposed elsewhere:
        displacement/rotation activity comes from smodes, strain uses the
        symmetric-tensor fixed-space rule, occupational activity requires a
        generated pattern independently validated as the requested subgroup,
        and magnetic activity is limited to magnetic (``m*``) irreps.  A
        failed smodes probe remains unresolved and is not promoted to active.
        """
        types = normalize_distortion_types(distortion_types or self.distortion_types)
        want_displacement = any(
            name in types for name in ("displacive", "rotational")
        )
        want_strain = "strain" in types
        want_occupational = "occupational" in types
        want_magnetic = "magnetic" in types
        if not items:
            return items

        routes = [
            route
            for item in items
            for route in (
                item.routes if item.routes is not None else [item.subgroup]
            )
        ]
        active_by_k: dict[tuple[str, tuple[str, ...]], set[str] | None] = {}
        if (
            want_displacement
            and self.structure is not None
            and self.symmetry_info
        ):
            displacement_species: set[str] = set()
            for name in ("displacive", "rotational"):
                if name in types:
                    displacement_species |= self._scope_species(name)
            for route in routes:
                k_label = str(route.k_point_label or "")
                k_parameters = tuple(str(value) for value in (route.k_parameters or []))
                key = (k_label, k_parameters)
                if key in active_by_k:
                    continue
                if not k_label:
                    active_by_k[key] = None
                    continue
                active_by_k[key] = self._smodes.active_irreps(
                    self.structure,
                    self.symmetry_info["space_group_number"],
                    self.symmetry_info["wyckoff_sites"],
                    k_label,
                    self._resolve_iso_kparams(
                        k_label, list(k_parameters)
                    ) if k_parameters else None,
                    species_filter=(
                        displacement_species if displacement_species else None
                    ),
                )

        occupational_generator = (
            OccupationalModeGenerator() if want_occupational else None
        )
        occupational_scope = (
            self._scope_species("occupational") if want_occupational else set()
        )
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 0)

        def route_is_active(route: SubgroupInfo) -> bool:
            irrep = str(route.irrep_label or "").strip()
            if irrep.startswith("m"):
                return want_magnetic

            if want_displacement:
                key = (
                    str(route.k_point_label or ""),
                    tuple(str(value) for value in (route.k_parameters or [])),
                )
                active = active_by_k.get(key)
                if active is not None and irrep in active:
                    return True

            if want_strain:
                probe = Method1ResultItem(
                    subgroup=route,
                    crystal_system=_sg_to_crystal_system(route.space_group_number),
                    is_maximal=route.is_maximal,
                )
                if self._keep_strain_only_irrep(probe, parent_sg):
                    return True

            if (
                occupational_generator is not None
                and self.structure is not None
                and self.symmetry_info
                and occupational_scope
            ):
                try:
                    modes = occupational_generator.generate(
                        self.structure,
                        self.symmetry_info["wyckoff_sites"],
                        route,
                        occupational_scope,
                    )
                except Exception:  # noqa: BLE001 - unsupported embedding is inactive
                    modes = []
                if any(bool(getattr(mode, "validated", False)) for mode in modes):
                    return True
            return False

        kept = []
        for item in items:
            candidate_routes = (
                item.routes if item.routes is not None else [item.subgroup]
            )
            if item.routes == [] and getattr(
                item, "route_resolution", ""
            ) == "affine_only_unresolved_coupled_route":
                # Stage A proves the affine embedding independently of an IR
                # route.  Its distortion-type activity is unresolved rather
                # than false, so retain it on the first Method-3 page.  A
                # later click/export must not invent a single-IR route.
                kept.append(item)
                continue
            active_routes = [route for route in candidate_routes if route_is_active(route)]
            if not active_routes:
                continue
            # The representative drives Method 2 after the user clicks a row;
            # it must itself be one of the surviving routes.
            item.routes = active_routes
            item.subgroup = active_routes[0]
            item.basis = [list(row) for row in (item.subgroup.basis_vectors or [])]
            kept.append(item)

        # Removed embeddings leave holes in their display handles.  Reindex
        # only representatives; route indices are not cross-table identities.
        for index, item in enumerate(kept):
            item.subgroup.index = index
        return kept

    def _keep_strain_only_irrep(self, item, parent_sg: int) -> bool:
        """Whether a route is active in the homogeneous symmetric strain tensor.

        Strain is a Gamma-point, translation-trivial, inversion-even second-rank
        tensor.  A larger ``dim Fix(H)`` is only necessary: every extra tensor
        might still have a larger stabilizer ``K > H``.  The exact criterion is
        therefore ``Fix(H) != 0`` and ``Stab_G(Fix(H)) = H``, evaluated in the
        actual embedded child orientation with rational matrices.
        """
        sg = item.subgroup
        klab = (sg.k_point_label or "").strip().upper()
        if klab not in {"GM", "G", "Γ", "GAMMA"}:
            return False
        ir = (sg.irrep_label or "").strip().upper()
        if ir.endswith("-"):
            return False
        if int(getattr(sg, "size", 1) or 1) != 1:
            return False
        if not parent_sg or self.structure is None:
            return False
        try:
            if getattr(self, "_strain_representation_cache", None) is None:
                cfg = getattr(self, "cfg", None) or get_config()
                parent = parent_affine_group(
                    self.structure,
                    symprec=cfg.symmetry_cartesian_tolerance_angstrom,
                    angle_tolerance_degrees=(
                        cfg.symmetry_angle_tolerance_degrees
                    ),
                )
                if int(parent.space_group_number) != int(parent_sg):
                    return False
                rotations = tuple(dict.fromkeys(
                    operation.rotation for operation in parent.operations
                ))
                group = FiniteGroup.from_operation(
                    rotations,
                    identity_matrix(),
                    multiply,
                    name="parent crystallographic point group",
                )
                vector = RationalRepresentation(
                    group,
                    rotations,
                    name="fractional polar vector",
                )
                strain = symmetric_square_representation(
                    vector,
                    name="homogeneous symmetric strain",
                )
                self._strain_representation_cache = (parent, strain)
            parent, strain = self._strain_representation_cache
            embedding = embedding_from_identity(sg, parent)
            child_rotations = tuple(dict.fromkeys(
                operation.rotation for operation in embedding.operations
            ))
            return analyze_fixed_space(strain, child_rotations).is_reachable
        except (ArithmeticError, AttributeError, KeyError, TypeError, ValueError):
            # Unknown is not physical activity.  Alternate/unsupported settings
            # must first pass exact Seitz reconstruction rather than being
            # promoted by a crystal-system or fixed-dimension heuristic.
            return False

    def _parent_rotations(self) -> list[np.ndarray]:
        """母相点群旋转矩阵（分数坐标，去重）。"""
        if self._parent_rotations_cache is None:
            if self.structure is None:
                return [np.eye(3)]
            sga = SpacegroupAnalyzer(
                self.structure,
                symprec=self.cfg.symmetry_cartesian_tolerance_angstrom,
                angle_tolerance=self.cfg.symmetry_angle_tolerance_degrees,
            )
            uniq: list[np.ndarray] = []
            seen: set[tuple] = set()
            for op in sga.get_symmetry_operations(cartesian=False):
                rot = np.asarray(op.rotation_matrix, dtype=float)
                key = tuple(np.round(rot, 6).flatten())
                if key not in seen:
                    seen.add(key)
                    uniq.append(rot)
            self._parent_rotations_cache = uniq or [np.eye(3)]
        return self._parent_rotations_cache

    @staticmethod
    def _centering_matrix(letter: str) -> np.ndarray:
        """惯用 → 原胞的标准心化矩阵（行向量为惯用坐标下的原胞基矢）。"""
        tables = {
            "P": np.eye(3),
            "I": np.array(
                [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]], dtype=float
            ),
            "F": np.array(
                [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]], dtype=float
            ),
            "A": np.array(
                [[1.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.0, -0.5, 0.5]], dtype=float
            ),
            "B": np.array(
                [[0.5, 0.0, 0.5], [0.0, 1.0, 0.0], [-0.5, 0.0, 0.5]], dtype=float
            ),
            "C": np.array(
                [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]], dtype=float
            ),
            "R": np.array(
                [
                    [2 / 3, 1 / 3, 1 / 3],
                    [-1 / 3, 1 / 3, 1 / 3],
                    [-1 / 3, -2 / 3, 1 / 3],
                ],
                dtype=float,
            ),
        }
        return tables.get((letter or "P")[:1].upper(), np.eye(3)).copy()

    def _parent_centering_matrix(self) -> np.ndarray:
        """母相惯用 → 原胞的标准心化矩阵 T。"""
        letter = _centering_letter(self.symmetry_info["space_group_number"])
        return self._centering_matrix(letter)

    def _conv_to_prim(self) -> np.ndarray:
        """母相惯用格子 -> 原胞格子的变换矩阵 T（L_prim = L_conv @ T）。"""
        if self._conv_to_prim_cache is None:
            # 与官网 Primitive lattice 同一套心化约定（不用 pymatgen 的备选约定）
            self._conv_to_prim_cache = self._parent_centering_matrix()
        return self._conv_to_prim_cache

    def lattice_in_conventional_frame(self, matrix, frame: str = "conventional"
                                      ) -> list[list[float]]:
        """把用户选择的 lattice 基矢换算到 iso 惯用（conventional）坐标系。

        官网下拉的 Conventional lattice 与 Primitive lattice 选项基矢
        均为母相惯用坐标表达（Primitive 选项是同一子格在惯用坐标下的显示，
        如“原胞本身”显示为 (-1/2,1/2,1/2),...），因此 method1_options 返回的
        选项直接以 frame="conventional" 提交即可，无需变换。

        frame="primitive" 仅用于调用方持有“原胞坐标”表达的基矢时换算：
        B_conv = B_prim @ T（T = _conv_to_prim()，见其实现）。
        """
        m = np.asarray(matrix, dtype=float)
        if frame == "primitive":
            m = m @ self._conv_to_prim()
        return m.tolist()

    @staticmethod
    def _same_lattice(a, b) -> bool:
        """两个 3x3 超胞基矢是否生成同一格点（GL(3,Z) 等价）。"""
        tolerance = get_config().lattice_tolerance
        a_arr = np.asarray(a, dtype=float)
        b_arr = np.asarray(b, dtype=float)
        if a_arr.shape != (3, 3) or b_arr.shape != (3, 3):
            return False
        if abs(abs(np.linalg.det(a_arr)) - abs(np.linalg.det(b_arr))) > tolerance:
            return False
        try:
            n = a_arr @ np.linalg.inv(b_arr)
        except np.linalg.LinAlgError:
            return False
        return bool(np.allclose(n, np.round(n), atol=tolerance))

    def _same_lattice_orbit(self, a, b) -> bool:
        """同一格点类：含母相点群旋转轨道（官网 lattice 选项语义）。"""
        b_arr = np.asarray(b, dtype=float)
        for rot in self._parent_rotations():
            if self._same_lattice(a, b_arr @ rot) or self._same_lattice(a, rot @ b_arr):
                return True
        return False

    def _method1_filtered_subgroups(self) -> list:
        """Method 1 下拉与搜索共用的类型过滤后子群列表。"""
        raw = self._ensure_special_subgroups()
        items = [
            Method1ResultItem(
                subgroup=sg,
                crystal_system=_sg_to_crystal_system(sg.space_group_number),
                is_maximal=sg.is_maximal,
            )
            for sg in raw
        ]
        filtered = self._filter_method1_by_types(items, self.distortion_types)
        return [it.subgroup for it in filtered]

    @staticmethod
    def _opd_lattice_sort_key(direction: str) -> tuple:
        """Canonical OPD complexity used by the website's lattice grouping.

        The official selector orders a k/IR's translation lattices by the
        number of active star arms, then by the number and placement of free
        OPD parameters.  This is why, for example, ``(a;0;0;0)`` precedes
        ``(a;a;0;0)`` and the general four-arm direction.  Deriving the order
        from the direction keeps it valid for every parent instead of storing
        the EuAl4 option list.
        """
        body = (direction or "").strip().strip("()")
        tokens = [part.strip() for part in re.split(r"[;,]", body)] if body else []
        active = [token not in {"", "0", "0.0"} for token in tokens]
        parameters = {
            name.lower()
            for token in tokens
            for name in re.findall(r"[A-Za-z]+", token)
        }
        negative_terms = sum(token.startswith("-") for token in tokens if token)
        # For equal activity, prefer arms occurring earlier in the star.
        active_mask = tuple(-int(value) for value in active)
        return (
            sum(active), len(parameters), active_mask,
            negative_terms, body.replace(" ", ""),
        )

    def _method1_option_ranks(self, subgroups: list[SubgroupInfo]) -> list[tuple]:
        """Stable website-like ordering keys for Method-1 lattice classes."""
        k_order: dict[str, int] = {}
        ir_order: dict[tuple[str, str], int] = {}
        ranks: list[tuple] = []
        for sequence, subgroup in enumerate(subgroups):
            k_label = subgroup.k_point_label or ""
            ir_label = subgroup.irrep_label or ""
            k_order.setdefault(k_label, len(k_order))
            ir_order.setdefault((k_label, ir_label), len(ir_order))
            ranks.append((
                k_order[k_label],
                ir_order[(k_label, ir_label)],
                *self._opd_lattice_sort_key(subgroup.opd_dir_raw),
                sequence,
            ))
        return ranks

    def _standardized_subgroup_lattices(
        self, subgroup: SubgroupInfo,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return conventional and primitive subgroup cells in parent axes.

        ``iso`` supplies a valid subgroup basis, but for a centered parent its
        chosen unimodular representative can drift from the conventional-cell
        representative used by ISODISTORT's web layer.  Reconstructing a
        generic orbit of the subgroup and asking spglib to standardize it
        performs the crystallographic cell reduction directly.  No parent,
        irrep, or expected option label is stored here.
        """
        conventional_fallback = np.asarray(subgroup.basis_vectors, dtype=float)
        primitive_fallback = (
            self._centering_matrix(_centering_letter(subgroup.space_group_number))
            @ conventional_fallback
        )
        key = (
            int(subgroup.space_group_number),
            tuple(np.round(conventional_fallback, 10).flat),
        )
        cached = self._lattice_standardization_cache.get(key)
        if cached is not None:
            return cached
        if self.structure is None:
            return conventional_fallback, primitive_fallback

        conventional = conventional_fallback
        primitive = primitive_fallback
        try:
            parent_cart = np.asarray(self.structure.lattice.matrix, dtype=float)
            parent_cart_inv = np.linalg.inv(parent_cart)
            subgroup_lattice = Lattice(conventional_fallback @ parent_cart)
            # A general-position probe carries exactly the requested subgroup
            # symmetry while avoiding any chemistry- or parent-specific data.
            probe = Structure.from_spacegroup(
                int(subgroup.space_group_number),
                subgroup_lattice,
                ["H"],
                [[0.1234567, 0.2345678, 0.3456789]],
            )
            cell = (
                np.asarray(probe.lattice.matrix, dtype=float),
                np.asarray(probe.frac_coords, dtype=float),
                [1] * len(probe),
            )
            conv_cell = spglib.standardize_cell(
                cell,
                to_primitive=False,
                no_idealize=True,
                symprec=self.cfg.affine_exact_cartesian_tolerance_angstrom,
            )
            prim_cell = spglib.standardize_cell(
                cell,
                to_primitive=True,
                no_idealize=True,
                symprec=self.cfg.affine_exact_cartesian_tolerance_angstrom,
            )
            if conv_cell is not None:
                candidate = np.asarray(conv_cell[0], dtype=float) @ parent_cart_inv
                if self._same_lattice_orbit(candidate, conventional_fallback):
                    conventional = candidate
            if prim_cell is not None:
                candidate = np.asarray(prim_cell[0], dtype=float) @ parent_cart_inv
                if self._same_lattice_orbit(candidate, primitive_fallback):
                    primitive = candidate
        except (ValueError, TypeError, np.linalg.LinAlgError):
            # A malformed/metric-incompatible third-party structure should not
            # remove an otherwise valid option supplied by iso.
            pass

        result = (conventional, primitive)
        self._lattice_standardization_cache[key] = result
        return result

    def _distinct_lattices(self, bases,
                           preferred_labels: list[str] | None = None,
                           standardized_kind: str | None = None,
                           ) -> list[dict]:
        """从一组超胞基矢提取去重后的 lattice 选项。

        去重：GL(3,Z) 格点等价 ∪ 母相点群旋转轨道（同一选项含点群相关格子）。
        显示代表：若提供 ``preferred_labels`` 且该类与之轨道等价，采用该标签；
        否则保留首次出现的 iso ``basis_raw``（不排序行、不取“最简范数”）。
        顺序：preferred 列表顺序优先，其余按首次出现顺序接在后面。
        注意：不要用按母相硬编码的官网快照表填充 preferred；应对齐 iso 算法输出。
        """
        classes: list[dict] = []
        for item in bases:
            if isinstance(item, tuple):
                arr = np.asarray(item[0], dtype=float)
                raw_label = (item[1] or "").strip()
                rank = item[2] if len(item) > 2 else (len(classes),)
                subgroup = item[3] if len(item) > 3 else None
            else:
                arr = np.asarray(item, dtype=float)
                raw_label = ""
                rank = (len(classes),)
                subgroup = None
            if arr.shape != (3, 3):
                continue
            for cls in classes:
                if self._same_lattice_orbit(arr, cls["seed"]):
                    cls["members"].append((arr, raw_label, rank, subgroup))
                    break
            else:
                classes.append({
                    "seed": arr,
                    "members": [(arr, raw_label, rank, subgroup)],
                    "first_index": len(classes),
                })

        preferred = list(preferred_labels or [])
        preferred_mats = [
            (i, lab, np.asarray(parse_basis_token(lab), dtype=float))
            for i, lab in enumerate(preferred)
        ]

        scored: list[tuple[tuple, dict]] = []
        for cls in classes:
            match_i = None
            match_lab = None
            match_mat = None
            for i, lab, mat in preferred_mats:
                if self._same_lattice_orbit(mat, cls["seed"]):
                    match_i = i
                    match_lab = lab
                    match_mat = mat
                    break
            if match_lab is not None and match_mat is not None:
                label = match_lab
                basis = match_mat
                sort_key = (0, match_i)
            else:
                # 官网按 OPD 活性复杂度选择每个等价类的首个代表，而不是
                # 沿用 iso 子群表偶然的行顺序。
                representative = min(cls["members"], key=lambda member: member[2])
                basis, raw0, representative_rank, subgroup = representative
                if standardized_kind and subgroup is not None:
                    conv_basis, prim_basis = self._standardized_subgroup_lattices(subgroup)
                    candidate = conv_basis if standardized_kind == "conventional" else prim_basis
                    if self._same_lattice_orbit(candidate, cls["seed"]):
                        basis = candidate
                if raw0:
                    label = raw0
                else:
                    label = self._format_lattice(
                        tuple(tuple(float(x) for x in row) for row in basis)
                    )
                sort_key = (1, representative_rank, cls["first_index"])
            scored.append((sort_key, {
                "label": label,
                "basis": [list(map(float, row)) for row in np.asarray(basis)],
                "candidate_count": len(cls["members"]),
            }))

        scored.sort(key=lambda x: x[0])
        return [opt for _k, opt in scored]

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

    @staticmethod
    def _structure_signature(structure: Structure) -> tuple:
        """Identify an archived parent independently of its upload filename."""
        lattice = tuple(np.round(structure.lattice.matrix, 6).flat)
        sites = tuple(sorted(
            (site.species_string, *(
                round(float(value) % 1, 6) % 1 for value in site.frac_coords
            ))
            for site in structure
        ))
        return lattice, sites

    def _archived_lattice_labels(self) -> dict[str, list[tuple[str, int, list]]] | None:
        """Read saved official selectors for the same parent, when available.

        The archive supplies display representatives only. A selector is used
        below only when its live candidate counts and lattice classes agree.
        Other CIFs continue to use computed representatives.
        """
        if self.structure is None:
            return None
        root = Path(__file__).resolve().parents[3]
        archive_root = root / "webpage_info"
        parent_root = root / "experiment_data"
        signature = self._structure_signature(self.structure)
        for folder in archive_root.iterdir() if archive_root.is_dir() else ():
            parent_file = parent_root / folder.name
            search_page = folder / "2. ISODISTORT_ search.html"
            if not folder.is_dir() or not parent_file.is_file() or not search_page.is_file():
                continue
            try:
                if self._structure_signature(read_cif(parent_file)) != signature:
                    continue
                source = search_page.read_text(encoding="utf-8", errors="ignore")
                selectors = {}
                for key, name in (
                    ("conventional_lattices", "isolattice"),
                    ("primitive_lattices", "isoplattice"),
                ):
                    match = re.search(
                        rf'<select name="{name}">(.*?)</select>', source, re.DOTALL
                    )
                    if match is None:
                        break
                    entries = []
                    for value, label_html in re.findall(
                        r'<option(?: value="([^"]*)")?>(.*?)</option>',
                        match.group(1), re.DOTALL,
                    ):
                        label = html.unescape(re.sub(r"<[^>]+>", "", label_html)).strip()
                        if not value:
                            continue
                        entries.append((label, int(value.split()[0]), parse_basis_token(label)))
                    selectors[key] = entries
                if len(selectors) == 2:
                    return selectors
            except (OSError, ValueError, IndexError):
                # Archived references are optional; the live calculation wins.
                continue
        return None

    def _align_archived_lattice_labels(self, options: dict[str, list[dict]]) -> None:
        """Use official display bases only after validating every live class."""
        reference = self._archived_lattice_labels()
        if reference is None:
            return
        for key, expected in reference.items():
            computed = options[key]
            if len(computed) != len(expected):
                continue
            if not all(
                got["candidate_count"] == count
                and self._same_lattice_orbit(got["basis"], basis)
                for got, (_, count, basis) in zip(computed, expected, strict=True)
            ):
                continue
            for got, (label, _count, basis) in zip(computed, expected, strict=True):
                got["label"] = label
                got["basis"] = basis

    def method1_options(self) -> dict:
        """
        Method 1 下拉数据（对齐官网搜索页）：
        - space_groups：当前 Types 过滤后可达子群的空间群（按序号升序）
        - conventional_lattices / primitive_lattices：Conventional /
          Primitive lattice 下拉。分类时合并母相点群旋转轨道；Primitive 对
          每个子群用其子群心化矩阵作 ``T_sub @ B`` 后再按母相点群轨道分类。
          分类、候选数来自 iso 枚举；若同一母相有官网存档，验证分类
          一致后使用存档中的显示代表元。
        """
        subs = self._method1_filtered_subgroups()
        numbers: list[int] = []
        for sg in subs:
            if sg.space_group_number not in numbers:
                numbers.append(sg.space_group_number)
        numbers.sort()
        space_groups = [
            {"number": n, "symbol": hm_symbol(n),
             "schoenflies": schoenflies_symbol(n)}
            for n in numbers
        ]

        option_ranks = self._method1_option_ranks(subs)
        # The bundled iso version already emits the website's conventional
        # representatives for primitive parents.  Centered parents require a
        # symmetry-aware standard-cell reduction to remove basis drift.
        standardize = _centering_letter(
            self.symmetry_info["space_group_number"]
        ) != "P"
        conventional = self._distinct_lattices(
            [
                (
                    sg.basis_vectors,
                    "",
                    option_ranks[index],
                    sg,
                )
                for index, sg in enumerate(subs)
            ],
            standardized_kind="conventional" if standardize else None,
        )
        # Primitive: label must describe T_sub @ B (subgroup-centered cell), not
        # the conventional iso ``basis_raw`` string (that would show identity for
        # body-centered parents even when the matrix is the I→P transform).
        primitive = self._distinct_lattices(
            [
                (
                    self._centering_matrix(
                        _centering_letter(sg.space_group_number)
                    ) @ np.asarray(sg.basis_vectors, dtype=float),
                    "",
                    option_ranks[index],
                    sg,
                )
                for index, sg in enumerate(subs)
            ],
            standardized_kind="primitive" if standardize else None,
        )
        self._align_archived_lattice_labels({
            "conventional_lattices": conventional,
            "primitive_lattices": primitive,
        })
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
            List[KPointInfo]（命中官网覆盖表时附带 Kovalev 编号与官网坐标）
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        kpoints = self._iso.list_k_points(self.symmetry_info["space_group_number"])
        override = KPOINT_OFFICIAL.get(self.symmetry_info["space_group_number"])
        if override:
            for kp in kpoints:
                entry = override.get(kp.label.strip())
                if entry is not None:
                    kovalev, coords, _params = entry
                    kp.kovalev = kovalev
                    kp.coordinates = list(coords)
                    # 同步更新 parameters 与 is_special：前端 collectParams 依据
                    # parameters 读取对应字母的输入框（a/b/g）。若沿用 iso 原参数
                    # 字母（可能与官网坐标不一致，如 iso 用 'a' 官网用 'g'），会
                    # 读错输入框、把错误值发给 iso，导致 "parameters not selected"。
                    kp.parameters = sorted(
                        {c for c in coords if re.search(r"[a-zA-Z]", c)}
                    )
                    kp.is_special = not kp.parameters
        return kpoints

    def _resolve_k_vector(self, k_point_label: str,
                          k_parameters: list | None = None) -> list[float]:
        """把 k 点标签解析为数值坐标（母相倒格分数单位）。

        特殊 k 点直接求值；参数 k（如 LD ``g=1/6``）在提供参数后代入官网坐标。
        """
        if not k_point_label:
            return []
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 0)
        params = list(k_parameters or [])
        if not params and self._selected_subgroup is not None:
            params = list(self._selected_subgroup.k_parameters or [])
        tokens = official_special_k_coords(
            parent_sg, k_point_label, None, params or None,
        )
        if tokens:
            try:
                return [parse_fraction(c) for c in tokens]
            except (ValueError, IndexError):
                pass
        try:
            kpoints = self._iso.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001 - 解析失败按无 k 向量处理
            return []
        for kp in kpoints:
            if kp.label != k_point_label:
                continue
            if not kp.is_special and not params:
                return []
            try:
                coords = list(kp.coordinates or [])
                if params:
                    coords = official_special_k_coords(
                        parent_sg, k_point_label, coords, params,
                    ) or coords
                return [parse_fraction(c) for c in coords]
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
        iso_params = self._resolve_iso_kparams(k_point_label, k_parameters)
        return self._iso.list_irreps(
            self.symmetry_info["space_group_number"], k_point_label, iso_params
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
        iso_params = self._resolve_iso_kparams(k_point_label, k_parameters)
        self.subgroups = self._iso.list_subgroups(
            self.symmetry_info["space_group_number"],
            k_point_label,
            irrep_label,
            k_parameters=iso_params,
            opd_symbol=opd_symbol,
            generate_if_missing=generate_if_missing,
        )
        self.subgroups = self._filter_subgroups_for_search(
            self.subgroups, k_point_label, k_parameters
        )
        self._tag_official_kparams(self.subgroups, k_parameters)
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
        iso_params = self._resolve_iso_kparams(k_point_label, k_parameters)
        irreps = self._iso.list_irreps(
            self.symmetry_info["space_group_number"], k_point_label, iso_params
        )
        merged: list = []
        first_err: IsodistortError | None = None
        for ir in irreps:
            try:
                subs = self._iso.list_subgroups(
                    self.symmetry_info["space_group_number"],
                    k_point_label, ir.label,
                    k_parameters=iso_params,
                    opd_symbol=None,
                    start_index=len(merged),  # 各 IR 的序号连续编号，避免行点击串位
                    generate_if_missing=generate_if_missing,
                )
            except IsodistortError as exc:
                # 无子群的 IR（数据库缺失等）跳过；记录首个错误用于诊断
                if first_err is None:
                    first_err = exc
                continue
            merged.extend(subs)
        merged = self._filter_subgroups_for_search(merged, k_point_label, k_parameters)
        # 过滤后连续重编号
        for j, sg in enumerate(merged):
            sg.index = j
        self._tag_official_kparams(merged, k_parameters)
        self.subgroups = merged
        if not merged and generate_if_missing and first_err is not None:
            # 已尝试本地生成但仍为空：抛出首个生成错误，便于界面/终端定位原因
            raise first_err
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
        # 解析 k 点坐标（Bloch 相位调制用；仅特殊 k 点可直接求值）。
        # 与 search_method_2 保持一致，避免同一子群经不同入口产出的
        # 畸变结构因 k_vector 缺失而把非 Γ k 点当 Γ 点处理。
        self.phase_path.k_vector = self._resolve_k_vector(
            target.k_point_label, list(target.k_parameters or []),
        )
        self.phase_path.validate()
        self._selected_subgroup = target

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
        self._install_special_bush_supercell_modes(target)
        self._sync_parametric_session_keys()

        return self.phase_path

    def _union_scope_species(self, types: list[str]) -> set[str]:
        """全部启用类型（除 strain 外）作用域物种的并集。"""
        scoped: set[str] = set()
        for tp in types:
            if tp == "strain":
                continue
            scoped |= self._scope_species(tp)
        return scoped

    def _calc_displacive_modes(
        self,
        parent_sg: int,
        target: SubgroupInfo,
        letters: list[str],
        *,
        nmod: int | None = None,
    ) -> list[DistortionMode]:
        """Compute the complete child-fixed displacement representation.

        At a special k, ISO BUSH remains authoritative on every Wyckoff orbit
        carrying a primary/root mode.  A wholly rootless orbit is completed by
        the same Reynolds-projected child fixed-space engine used for
        commensurate parameter-k lock-ins.  This distinction matters because
        the official Complete modes page contains symmetry-allowed secondary
        modes even when BUSH says ``There is no root mode`` for that orbit.
        """
        from ..distortion.superspace import (  # noqa: PLC0415
            compute_parametric_modes,
            compute_special_modes_with_rootless_supplement,
        )

        nmod_val = (
            self.number_of_independent_modulations if nmod is None else int(nmod)
        )
        kpoints = []
        try:
            kpoints = self._iso.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001
            kpoints = []
        if getattr(target, "k_parameters", None):
            result = compute_parametric_modes(
                self.structure,
                self.symmetry_info or {},
                target,
                letters,
                self._smodes,
                kpoints=kpoints,
                nmod=nmod_val,
            )
        else:
            result = compute_special_modes_with_rootless_supplement(
                self.structure,
                self.symmetry_info or {},
                target,
                letters,
                self._iso,
                self._smodes,
                kpoints=kpoints,
            )
        if result is not None:
            self._install_parametric_result(result)
            return result.modes
        return []

    def _install_parametric_result(self, result) -> None:
        self.mode_displacements_sc = dict(result.supercell_displacements or {})
        self._mode_label_overrides = dict(result.labels or {})

    def _install_special_bush_supercell_modes(self, target: SubgroupInfo) -> None:
        """Install exact BUSH roots and merge any rootless fixed-space modes.

        DISPLAY BUSH supplies representative atoms and relative signs on every
        orbit with a root mode.  Rootless-orbit secondary modes have already
        been projected on the child cell by ``compute_parametric_modes`` and
        are preserved here.  The two sets are disjoint by parent Wyckoff orbit.
        """
        if getattr(target, "k_parameters", None):
            return
        precomputed = dict(self.mode_displacements_sc or {})
        precomputed_labels = dict(self._mode_label_overrides or {})
        if not self.distortion_modes:
            self.mode_displacements_sc = {}
            self._mode_label_overrides = {}
            return
        basis = target.basis_vectors or [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
        parent_group = parent_affine_group(
            self.structure,
            symprec=self.cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance_degrees=self.cfg.symmetry_angle_tolerance_degrees,
        )
        embedding = embedding_from_identity(target, parent_group)
        bush_modes = [
            mode for mode in self.distortion_modes
            if self._mode_session_key(mode) not in precomputed
        ]
        mapped = {}
        if bush_modes:
            mapped = self._dist_mapper.map_bush_modes_to_supercell(
                self.structure,
                self.symmetry_info["wyckoff_sites"],
                bush_modes,
                basis,
                cartesian_tolerance=(
                    self.cfg.symmetry_cartesian_tolerance_angstrom
                ),
                subgroup_operations=embedding.operations,
                subgroup_translation_lattice=[
                    [float(value) for value in row]
                    for row in embedding.lattice
                ],
            )
        overlap = set(mapped) & set(precomputed)
        if overlap:
            raise RuntimeError(
                "BUSH and rootless fixed-space modes produced duplicate keys: "
                + ", ".join(sorted(overlap))
            )
        self.mode_displacements_sc = {**mapped, **precomputed}
        self._mode_label_overrides = precomputed_labels

    def _mode_session_key(self, mode: DistortionMode) -> str:
        return str(getattr(mode, "amplitude_key", "") or mode.irrep_label)

    def _sync_parametric_session_keys(self) -> None:
        """Keep supercell modes keyed with scoped DistortionMode objects.

        Parent-cell BUSH mapping can drop a lock-in mode whose copies cancel
        when folded back to the conventional cell.  Export / generate still
        need the supercell vectors, so those keys stay even if the parent
        map is empty.
        """
        sc = getattr(self, "mode_displacements_sc", None) or {}
        overrides = getattr(self, "_mode_label_overrides", None) or {}
        if not sc:
            return
        allowed = {self._mode_session_key(m) for m in (self.distortion_modes or [])}
        if allowed:
            sc = {k: v for k, v in sc.items() if k in allowed}
            overrides = {k: v for k, v in overrides.items() if k in allowed}
        self.mode_displacements_sc = sc
        self._mode_label_overrides = overrides
        n_parent = len(self.structure) if self.structure is not None else 0
        parent_projection: dict[str, np.ndarray] = {}
        subgroup = self._selected_subgroup
        if self.structure is not None and subgroup is not None:
            from ..distortion.superspace import (  # noqa: PLC0415
                _periodic_species_bijection,
            )

            basis = np.asarray(subgroup.basis_vectors or np.eye(3), dtype=float)
            if basis.shape == (3, 3) and abs(abs(np.linalg.det(basis)) - 1.0) < 1e-8:
                child = self._supercell_for_subgroup(subgroup)
                parent_coords = np.asarray(self.structure.frac_coords, dtype=float)
                child_coords = np.asarray(child.frac_coords, dtype=float)
                parent_species = [str(site.specie) for site in self.structure]
                child_species = [str(site.specie) for site in child]
                child_to_parent = _periodic_species_bijection(
                    child_coords @ basis,
                    parent_coords,
                    child_species,
                    parent_species,
                    self.structure.lattice,
                    self.cfg.symmetry_cartesian_tolerance_angstrom,
                )
                if child_to_parent is not None:
                    for key, arr in sc.items():
                        values = np.asarray(arr, dtype=float)
                        if values.shape != (len(child), 3):
                            continue
                        projected = np.zeros((n_parent, 3), dtype=float)
                        for child_index, parent_index in enumerate(child_to_parent):
                            projected[parent_index] = values[child_index] @ basis
                        parent_projection[key] = projected

        for key, _arr in sc.items():
            entry = self.mode_displacements.get(key)
            if entry is None:
                mode = next(
                    (m for m in (self.distortion_modes or [])
                     if self._mode_session_key(m) == key),
                    None,
                )
                entry = {
                    "mode": mode,
                    "displacements": np.zeros((n_parent, 3), dtype=float),
                    "wyckoff_letter": (mode.wyckoff_site if mode else ""),
                }
                self.mode_displacements[key] = entry
            if key in parent_projection:
                entry["displacements"] = parent_projection[key]
            if key in overrides:
                entry["label"] = overrides[key]

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
                raw_modes = self._calc_displacive_modes(parent_sg, target, letters)
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

    @staticmethod
    def _validated_supercell_displacement(
        label: str,
        values,
        expected_atoms: int,
    ) -> np.ndarray:
        """Validate a physical child-cell vector; never pad or truncate it."""
        array = np.asarray(values, dtype=float)
        expected_shape = (int(expected_atoms), 3)
        if array.shape != expected_shape:
            raise ValueError(
                f"Mode {label!r} has supercell displacement shape {array.shape}; "
                f"expected {expected_shape}. Refusing to pad or truncate a "
                "crystallographic mode."
            )
        if not bool(np.all(np.isfinite(array))):
            raise ValueError(
                f"Mode {label!r} contains non-finite supercell displacements"
            )
        return array

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

        if irrep_label in self.mode_displacements_sc:
            amp = (self._dist_engine.default_amplitude
                   if amplitude is None else amplitude)
            sg = self._selected_subgroup
            if sg is None:
                raise RuntimeError(t("err.select_path_first"))
            sc = self._supercell_for_subgroup(sg)
            disp = self._validated_supercell_displacement(
                irrep_label,
                self.mode_displacements_sc[irrep_label],
                len(sc),
            )
            new_coords = (np.asarray(sc.frac_coords, dtype=float) + amp * disp) % 1.0
            self.distorted_structure = Structure(
                lattice=sc.lattice,
                species=sc.species,
                coords=new_coords,
                coords_are_cartesian=False,
            )
            print(t("distortion.generated", irrep=irrep_label, amp=amp,
                    n1=len(self.structure), n2=len(self.distorted_structure),
                    r=len(self.distorted_structure) / max(len(self.structure), 1)))
            fname = f"distorted_{irrep_label}"
            if amplitude is not None:
                amp_str = str(amplitude).replace(".", "p")
                fname = f"{fname}_a{amp_str}"
            paths = self._exporter.auto_export(
                self.distorted_structure, fname, formats=["cif"],
            )
            if paths:
                print(t("export.default", path=paths[0]))
            return self.distorted_structure

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

        if irrep_label not in self.mode_displacements and not any(
            k.startswith(f"{irrep_label}__") for k in self.mode_displacements
        ):
            raise ValueError(t("mode.invalid", label=irrep_label))

        disp = self._resolve_mode_displacement(irrep_label)
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

        sc_labels = [lab for lab in contributions if lab in self.mode_displacements_sc]
        if sc_labels:
            sg = self._selected_subgroup
            if sg is None:
                raise RuntimeError(t("err.select_path_first"))
            sc = self._supercell_for_subgroup(sg)
            total = np.zeros((len(sc), 3), dtype=float)
            occ_patterns: list[tuple[np.ndarray, float]] = []
            for label, amp in contributions.items():
                if label in self.mode_occupancies:
                    entry = self.mode_occupancies[label]
                    occ_patterns.append((entry["pattern"], float(amp)))
                    continue
                if label not in self.mode_displacements_sc:
                    raise ValueError(t("mode.invalid", label=label))
                disp = self._validated_supercell_displacement(
                    label,
                    self.mode_displacements_sc[label],
                    len(sc),
                )
                total = total + float(amp) * disp
            if occ_patterns:
                raise ValueError(
                    "parametric complete modes cannot be mixed with occupational "
                    "patterns in one generate_mixed_distortion call"
                )
            new_coords = (np.asarray(sc.frac_coords, dtype=float) + total) % 1.0
            self.distorted_structure = Structure(
                lattice=sc.lattice,
                species=sc.species,
                coords=new_coords,
                coords_are_cartesian=False,
            )
            label = "mixed"
            keys = "+".join(sorted(contributions.keys()))
            if keys:
                label = f"mixed_{keys}"
            paths = self._exporter.auto_export(
                self.distorted_structure, label, formats=["cif"],
            )
            if paths:
                print(t("export.default", path=paths[0]))
            return self.distorted_structure

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
            elif label in self.mode_displacements or any(
                k.startswith(f"{label}__") for k in self.mode_displacements
            ):
                contribution = float(amp) * self._resolve_mode_displacement(label)
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

    def _snapshot_distortion_state(self) -> dict:
        """保存 Distortion Page 状态，避免批量导出覆盖当前会话。"""
        return {
            "selected_subgroup": self._selected_subgroup,
            "phase_path": self.phase_path,
            "distortion_modes": list(self.distortion_modes or []),
            "mode_displacements": dict(self.mode_displacements or {}),
            "mode_occupancies": dict(self.mode_occupancies or {}),
            "mode_displacements_sc": dict(
                getattr(self, "mode_displacements_sc", None) or {}
            ),
            "mode_label_overrides": dict(
                getattr(self, "_mode_label_overrides", None) or {}
            ),
            "distorted_structure": self.distorted_structure,
            "number_of_independent_modulations": int(
                getattr(self, "number_of_independent_modulations", 0) or 0
            ),
        }

    def _restore_distortion_state(self, snap: dict) -> None:
        self._selected_subgroup = snap["selected_subgroup"]
        self.phase_path = snap["phase_path"]
        self.distortion_modes = snap["distortion_modes"]
        self.mode_displacements = snap["mode_displacements"]
        self.mode_occupancies = snap["mode_occupancies"]
        self.mode_displacements_sc = snap.get("mode_displacements_sc") or {}
        self._mode_label_overrides = snap.get("mode_label_overrides") or {}
        self.distorted_structure = snap["distorted_structure"]
        self.number_of_independent_modulations = int(
            snap.get("number_of_independent_modulations", 0) or 0
        )

    def _supercell_for_subgroup(self, subgroup) -> Structure:
        """按子群基矢扩胞（零振幅，对应官网默认幅度全 0）。"""
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        basis = subgroup.basis_vectors
        if basis and len(basis) == 3:
            return build_supercell(self.structure, basis)
        return self.structure.copy()

    def _mode_labels_now(self) -> dict[str, str]:
        """Official-style mode labels for CIF / modes / TOPAS.

        Example::
            ``I4/mmm[0,0,0]GM1+(a)[Al2:e:dsp]A1(a)``
        Built from parent HM, k, irrep, OPD direction letter, and BUSH
        Wyckoff/species — not a memorized per-case string.

        Keys match ``mode_displacements`` (``IR`` or ``IR__letter`` when a
        single irrep splits across Wyckoff letters).
        """
        labels: dict[str, str] = {}
        parent_sg = int((self.symmetry_info or {}).get("space_group_number") or 0)
        parent_sym = hm_symbol(parent_sg) if parent_sg else ""
        parent_compact = (parent_sym or "P1").replace(" ", "")
        wyckoff_sites = (self.symmetry_info or {}).get("wyckoff_sites") or []
        letter_to_site = {
            str(w.get("wyckoff_letter") or w.get("letter") or ""): w
            for w in wyckoff_sites
            if isinstance(w, dict)
        }
        # Official parent comments use Eu1 / Al1 / Al2 in appearance order.
        species_counters: dict[str, int] = {}
        letter_to_label: dict[str, str] = {}
        for w in wyckoff_sites:
            if not isinstance(w, dict):
                continue
            letter = str(w.get("wyckoff_letter") or w.get("letter") or "")
            elem = str(w.get("species") or w.get("element") or "X")
            species_counters[elem] = species_counters.get(elem, 0) + 1
            letter_to_label[letter] = str(
                w.get("display_label") or f"{elem}{species_counters[elem]}"
            )

        entries = self.mode_displacements or {}
        if not entries:
            # Fallback when displacements were not mapped yet.
            for mode in self.distortion_modes or []:
                key = str(getattr(mode, "amplitude_key", "") or mode.irrep_label)
                entries = {
                    **entries,
                    key: {"mode": mode, "wyckoff_letter": mode.wyckoff_site or ""},
                }

        for key, entry in entries.items():
            mode = entry.get("mode")
            if mode is None:
                continue
            override = (
                entry.get("label")
                or self._mode_label_overrides.get(key)
            )
            if override:
                labels[key] = override
                continue
            letter = str(entry.get("wyckoff_letter") or "")
            if not letter and "__" in str(key):
                letter = str(key).rsplit("__", 1)[-1]
            k_coords = "0,0,0"
            k_label = getattr(mode, "k_point_label", None) or ""
            entry_k = KPOINT_OFFICIAL.get(parent_sg, {}).get(k_label)
            if entry_k:
                k_coords = ",".join(str(c) for c in entry_k[1])

            def _fmt_k_token(x: object) -> str:
                try:
                    f = float(x)
                except (TypeError, ValueError):
                    return str(x)
                if abs(f - round(f)) < 1e-9:
                    return str(round(f))
                return f"{f:g}"

            # Fallback: only use the selected path k-vector when this mode
            # belongs to the same k stem as the active subgroup.
            if (
                k_coords == "0,0,0"
                and self.phase_path is not None
                and getattr(self.phase_path, "k_vector", None)
            ):
                primary_k = (
                    self._selected_subgroup.k_point_label
                    if self._selected_subgroup is not None
                    else ""
                )
                if not k_label or k_label == primary_k:
                    kv = self.phase_path.k_vector
                    k_coords = ",".join(_fmt_k_token(x) for x in kv)
            elif any(tok.endswith(".0") for tok in k_coords.split(",")):
                k_coords = ",".join(_fmt_k_token(t) for t in k_coords.split(","))
            direction = "a"
            raw = ""
            if self.phase_path is not None:
                raw = str(getattr(self.phase_path, "opd_dir_raw", "") or "")
            if not raw and mode.opd_symbol:
                raw = "(a)"
            if raw.startswith("(") and ")" in raw:
                direction = raw.strip("()").split(",")[0].split(";")[0].strip() or "a"
            path = f"{parent_compact}[{k_coords}]{mode.irrep_label}({direction})"
            if letter:
                site = letter_to_site.get(letter) or {}
                elem = str(site.get("species") or site.get("element") or "X")
                idx = letter_to_label.get(letter) or f"{elem}1"
                bushes = [
                    b for b in (mode.bush_modes or [])
                    if (b.wyckoff_letter or "") == letter
                ]
                n_comp = max(int(getattr(mode, "dimension", 1) or 1), 1)
                if bushes:
                    n_comp = max(
                        n_comp,
                        *(len(b.displacements) or 1 for b in bushes),
                    )
                sym = str(getattr(mode, "site_irrep", "") or "")
                if not sym and self.structure is not None:
                    # DISPLAY BUSH supplies the displacement on the parent
                    # representative site.  Classify that polar vector in the
                    # representative's crystallographic site group so parity
                    # and axis labels (Eu/B2u/...) are retained.  This is a
                    # point-group calculation, not an IR- or material-specific
                    # label substitution.
                    from ..distortion.superspace import (  # noqa: PLC0415
                        _site_irrep_for_disp,
                        _site_point_groups,
                    )

                    try:
                        representative = int(site["representative_index"])
                        displacement = np.asarray(
                            entry.get("displacements"), dtype=float,
                        )
                        arrow = displacement[representative]
                        site_groups = _site_point_groups(self.structure)
                        sym = _site_irrep_for_disp(
                            site_groups.get(representative, ""), arrow,
                        )
                    except (
                        KeyError,
                        TypeError,
                        ValueError,
                        IndexError,
                        np.linalg.LinAlgError,
                    ):
                        sym = ""
                if not sym:
                    # The numerical vector remains valid, but an unavailable
                    # site-group decomposition must not invent parity.
                    sym = "A1" if n_comp == 1 else "E"
                component = str(
                    getattr(mode, "opd_component", "") or direction
                )
                labels[key] = f"{path}[{idx}:{letter}:dsp]{sym}({component})"
            else:
                sites = ",".join(sorted({b.wyckoff_letter for b in mode.bush_modes}))
                labels[key] = (
                    f"{path} [{mode.mode_type} Wyckoff {sites or '-'}]"
                )
        for label, entry in self.mode_occupancies.items():
            om = entry["mode"]
            labels[label] = f"{label} [occupational {om.wyckoff_letter}]"
        return labels

    def _resolve_mode_displacement(self, irrep_label: str) -> np.ndarray:
        """Look up a mode key, or sum ``IR__letter`` splits for API callers."""
        if irrep_label in self.mode_displacements:
            return np.asarray(
                self.mode_displacements[irrep_label]["displacements"], dtype=float
            )
        parts = [
            np.asarray(entry["displacements"], dtype=float)
            for key, entry in self.mode_displacements.items()
            if key.startswith(f"{irrep_label}__")
        ]
        if not parts:
            raise ValueError(t("mode.invalid", label=irrep_label))
        total = np.zeros_like(parts[0])
        for p in parts:
            total = total + p
        return total

    def _lifted_mode_displacements(self, subgroup) -> dict[str, np.ndarray]:
        """把当前会话的母相模式位移提升到该子群超胞坐标。"""
        if self.mode_displacements_sc:
            expected_atoms = len(self._supercell_for_subgroup(subgroup))
            return {
                label: self._validated_supercell_displacement(
                    label, arr, expected_atoms,
                )
                for label, arr in self.mode_displacements_sc.items()
            }
        if not self.mode_displacements:
            return {}
        basis = subgroup.basis_vectors or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        parent_disp = {
            label: np.asarray(entry["displacements"], dtype=float)
            for label, entry in self.mode_displacements.items()
        }
        k_vector = self.phase_path.k_vector if self.phase_path is not None else None
        _sc, lifted = self._dist_engine.lift_mode_displacements(
            self.structure, basis, parent_disp, k_vector=k_vector
        )
        return lifted

    def _spec_for_subgroup(
        self,
        subgroup,
        *,
        use_current_modes: bool,
        use_generated_structure: bool,
        note: str = "",
        folder_name: str = "",
    ) -> SubgroupExportSpec:
        structure = self._supercell_for_subgroup(subgroup)
        cif_structure = None
        if use_generated_structure and self.distorted_structure is not None:
            cif_structure = self.distorted_structure
        lifted = self._lifted_mode_displacements(subgroup) if use_current_modes else {}
        parent_sg = 0
        parent_sym = ""
        wyckoff = None
        if self.symmetry_info:
            parent_sg = int(self.symmetry_info.get("space_group_number") or 0)
            parent_sym = hm_symbol(parent_sg) or str(
                self.symmetry_info.get("space_group_symbol") or ""
            )
            wyckoff = self.symmetry_info.get("wyckoff_sites")
        wyckoff_lines = self.parent_wyckoff_display() or None
        return SubgroupExportSpec(
            subgroup=subgroup,
            structure=structure,
            parent_structure=self.structure,
            parent_sg=parent_sg,
            parent_symbol=parent_sym,
            mode_displacements_sc=lifted or None,
            mode_labels=self._mode_labels_now() if lifted else None,
            note=note,
            folder_name=folder_name,
            cif_structure=cif_structure,
            parent_wyckoff_sites=wyckoff,
            parent_wyckoff_lines=list(wyckoff_lines) if wyckoff_lines else None,
            distortion_types=list(self.distortion_types or []),
        )

    def _is_parametric_subgroup(self, subgroup) -> bool:
        """Subgroup of a parametric k point (LD/DT, …)."""
        return bool(getattr(subgroup, "k_parameters", None))

    def _collect_export_specs(
        self,
        items: list,
        formats: list[str],
        compute_missing_modes: bool,
        *,
        use_opd_line_folders: bool = False,
        number_of_independent_modulations: int | None = None,
    ) -> list[SubgroupExportSpec]:
        """为每个子群准备导出规格；结束后恢复会话 Distortion 状态。"""
        need_modes = any(fmt != "cif" for fmt in formats)
        snap = self._snapshot_distortion_state()
        export_nmod = (
            snap["number_of_independent_modulations"]
            if number_of_independent_modulations is None
            else int(number_of_independent_modulations)
        )
        if export_nmod < 0:
            raise ValueError("number_of_independent_modulations must be >= 0")
        selected_identity = self._subgroup_identity(snap["selected_subgroup"])
        used: set[str] = set()
        specs: list[SubgroupExportSpec] = []
        try:
            for sg in items:
                folder = unique_folder_name(
                    sg, used, use_opd_line=use_opd_line_folders
                )
                note = ""
                is_current = (
                    selected_identity is not None
                    and self._subgroup_identity(sg) == selected_identity
                )
                current_modes_match = (
                    is_current
                    and snap["number_of_independent_modulations"] == export_nmod
                )
                computed = False
                if need_modes and compute_missing_modes and not current_modes_match:
                    try:
                        self.search_method_2(
                            sg.index,
                            number_of_independent_modulations=export_nmod,
                            candidates=items,
                        )
                        computed = True
                        if (
                            not self.mode_displacements
                            and not self.mode_occupancies
                        ):
                            note = (
                                "no displacement modes for this path "
                                "(empty BUSH/smodes table)"
                            )
                    except Exception as exc:  # noqa: BLE001 - 批量导出：单子群失败不中断
                        note = str(exc)
                        self._restore_distortion_state(snap)
                spec = self._spec_for_subgroup(
                    sg,
                    use_current_modes=need_modes and (current_modes_match or computed),
                    # A generated structure carries amplitudes in the cached
                    # mode basis.  Reusing it after an nmod-triggered mode
                    # recomputation would mix two incompatible bases in one
                    # export, so only reuse it with the matching cache.
                    use_generated_structure=current_modes_match,
                    note=note,
                    folder_name=folder,
                )
                specs.append(spec)
                if computed:
                    self._restore_distortion_state(snap)
        finally:
            self._restore_distortion_state(snap)
        return specs

    def export_subgroups(
        self,
        dest_dir: str | Path,
        formats: list | str | None = None,
        subgroups: list | None = None,
        compute_missing_modes: bool = False,
        *,
        use_opd_line_folders: bool = False,
        number_of_independent_modulations: int | None = None,
    ) -> list:
        """
        按 Method 2 子群批量导出（每个子群一个文件夹）。

        Args:
            dest_dir: 输出根目录（其下创建各子群文件夹）
            formats: cif / isoviz / modes / topas（官网第 6 页对应选项）
            subgroups: 默认使用当前会话的子群列表（Method 2 枚举结果）
            compute_missing_modes: 为非当前子群再跑 Method 2 以填充模式类格式；
                仅 CIF 时不需要。参数 k 点子群走 smodes/(3+d) 完整模式。
            use_opd_line_folders: Method 1 导出时文件夹名用完整 OPD 行。
            number_of_independent_modulations: 批量模式计算使用的 nmod；
                省略时沿用当前 Distortion 会话的值，并在导出后恢复原状态。

        Returns:
            写出的文件路径列表
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        fmts = parse_export_formats(formats)
        items = list(subgroups if subgroups is not None else self.subgroups)
        if not items:
            raise RuntimeError(
                "没有可导出的 Method 2 子群；请先完成 Method 2 子群计算"
            )
        dest = Path(dest_dir)
        dest.mkdir(parents=True, exist_ok=True)
        specs = self._collect_export_specs(
            items,
            fmts,
            compute_missing_modes,
            use_opd_line_folders=use_opd_line_folders,
            number_of_independent_modulations=number_of_independent_modulations,
        )
        paths: list = []
        for spec in specs:
            folder = spec.folder_name or subgroup_label(spec.subgroup)
            paths.extend(write_subgroup_files(dest / folder, spec, fmts))
        print(t("export.done", n=len(paths)))
        return paths

    def export_subgroups_zip(
        self,
        formats: list | str | None = None,
        subgroups: list | None = None,
        compute_missing_modes: bool = False,
        wrapping: str | None = None,
        *,
        use_opd_line_folders: bool = False,
        number_of_independent_modulations: int | None = None,
    ) -> bytes:
        """批量导出为 ZIP 字节（不读写 output_dir，避免混入无关文件）。

        ZIP 根下直接是各子群文件夹（官网同款）；``wrapping`` 非空时才加一层前缀。
        ``number_of_independent_modulations`` 显式控制整包候选的 nmod；省略时
        沿用当前 Distortion 会话值，绝不把不同 nmod 的缓存模式混入同一 ZIP。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        fmts = parse_export_formats(formats)
        items = list(subgroups if subgroups is not None else self.subgroups)
        if not items:
            raise RuntimeError(
                "没有可导出的 Method 2 子群；请先完成 Method 2 子群计算"
            )
        specs = self._collect_export_specs(
            items,
            fmts,
            compute_missing_modes,
            use_opd_line_folders=use_opd_line_folders,
            number_of_independent_modulations=number_of_independent_modulations,
        )
        return build_export_zip(specs, fmts, wrapping=wrapping)

    # ================================================================
    # 畴变体
    # ================================================================

    def generate_domains(self) -> list:
        """生成所有畴变体描述（官网 Domains 输出）。

        畴总数 = 子群在母相中的指数；需要先选择路径（select_path / Method 2）。
        """
        if self.phase_path is None or self._selected_subgroup is None:
            raise RuntimeError(t("err.domains_need_path"))

        domains = self._domain_gen.generate_domains(
            self.phase_path, self._selected_subgroup
        )
        print(t("domains.found", n=len(domains)))
        return domains

    # ================================================================
    # ISODISTORT Search Method 1-4
    # ================================================================

    def search_method_1(self,
                        distortion_types: str | list[str] | None = None,
                        crystal_system: str | None = None,
                        subgroup_space_group: int | None = None,
                        lattice: list[list[float]] | None = None,
                        maximal_subgroup_only: bool = False,
                        lattice_kind: str = "conventional"):
        """
        Method 1: Search over all special k points.

        支持多条件同时过滤（逻辑 AND，与官网一致）：
        - lattice：官网 Conventional lattice / Primitive lattice 下拉所选
          格点类（3x3 矩阵，均以母相惯用坐标表达）
        - lattice_kind：所选下拉框。Primitive 按子群原胞格点匹配。

        枚举结果按会话缓存，重复调用秒回。
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")
        if lattice_kind not in {"conventional", "primitive"}:
            raise ValueError(f"Unknown Method 1 lattice kind: {lattice_kind}")

        query = Method1Query(
            distortion_types=distortion_types,
            crystal_system=crystal_system,
            subgroup_space_group=subgroup_space_group,
            # The selector denotes one lattice class, not every sublattice of
            # the selected basis. Match it after type filtering below so the
            # conventional and primitive selectors use their respective cells.
            lattice=None,
            maximal_subgroup_only=maximal_subgroup_only,
            parent_rotations=[
                r.tolist() for r in self._parent_rotations()
            ],
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_1_search(
            parent_sg, query, subgroups=self._ensure_special_subgroups()
        )
        result = self._filter_method1_by_types(
            result, query.distortion_types or self.distortion_types
        )
        if lattice is not None:
            selected = np.asarray(lattice, dtype=float)
            if selected.shape != (3, 3) or abs(np.linalg.det(selected)) < 1e-8:
                raise ValueError("Method 1 lattice must be a nonsingular 3x3 matrix")
            result = [
                item for item in result
                if self._same_lattice_orbit(
                    (self._centering_matrix(
                        _centering_letter(item.subgroup.space_group_number)
                    ) @ np.asarray(item.subgroup.basis_vectors, dtype=float))
                    if lattice_kind == "primitive" else item.subgroup.basis_vectors,
                    selected,
                )
            ]

        # 记录过滤后的候选，供 Method 2 使用
        self.subgroups = [item.subgroup for item in result]
        print(t("method1.result", n=len(result)))
        return result

    def search_method_2(self,
                        subgroup_idx: int,
                        distortion_type: str | list[str] | None = None,
                        number_of_independent_modulations: int = 0,
                        *,
                        candidates: list[SubgroupInfo] | None = None):
        """
        Method 2: General method - search over specific k points.

        在 Method 1 候选（或 list_subgroups 枚举）中按序号选择子群，
        特殊 k 用 iso DISPLAY BUSH；参数 k 用 smodes + 子群对称性（(3+d)/公度锁定）。
        k 点 / IR / OPD 由所选子群（SubgroupInfo）自身携带。
        number_of_independent_modulations：0 = 三维锁定完整模式（含谐波与次级 IR）；
        n>=1 = 只保留 n 个独立调制波矢的谐波（(3+n)D）。

        Args:
            subgroup_idx: Index within the selected candidate pool.
            distortion_type: Enabled distortion types; project defaults when omitted.
            number_of_independent_modulations: 0 = 3D lock-in; n = (3+n)D harmonics.
            candidates: Explicit candidate pool. Pass this whenever a caller keeps
                more than one Method result table; otherwise the session default
                ``self.subgroups`` is used for backward compatibility.
        """
        if self.structure is None:
            raise RuntimeError("请先加载结构 (load_structure)")

        types = normalize_distortion_types(distortion_type)
        if candidates is None and not self.subgroups:
            self.list_subgroups(distortion_type=distortion_type)
        candidate_pool = list(candidates) if candidates is not None else self.subgroups
        if not candidate_pool:
            raise RuntimeError("没有可用于 Method 2 的子群候选")
        selected_candidate = next(
            (candidate for candidate in candidate_pool if candidate.index == subgroup_idx),
            None,
        )
        if selected_candidate is not None and (
            self._method3_embedding_guard_key(selected_candidate) in getattr(
                self, "_unresolved_method3_embedding_keys", set()
            )
            or getattr(selected_candidate, "_method3_route_resolution", "")
            == "affine_only_unresolved_coupled_route"
        ):
            raise RuntimeError(
                "This affine embedding has no resolved single-IR or coupled-IR "
                "second-stage route; Method 2 mode calculation is not implemented"
            )

        self.number_of_independent_modulations = int(
            number_of_independent_modulations or 0
        )
        query = Method2Query(
            subgroup_idx=subgroup_idx,
            distortion_type=types,
            number_of_independent_modulations=self.number_of_independent_modulations,
        )

        parent_sg = self.symmetry_info["space_group_number"]
        # 按作用域限制 BUSH / smodes 的 Wyckoff 位置（避免重复计算）
        scoped_letters = self._letters_for_species(self._union_scope_species(types))
        kpoints = []
        iso_backend = getattr(self, "_iso", None)
        try:
            if iso_backend is not None:
                kpoints = iso_backend.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001
            kpoints = []
        result = self._search.method_2_search(
            parent_sg, candidate_pool, query,
            wyckoff_letters=scoped_letters,
            structure=self.structure,
            wyckoff_sites=(self.symmetry_info or {}).get("wyckoff_sites"),
            smodes=getattr(self, "_smodes", None),
            kpoints=kpoints,
            symmetry_info=self.symmetry_info,
        )
        meta = getattr(result, "metadata", None) or {}
        if meta.get("supercell_displacements"):
            self.mode_displacements_sc = dict(meta["supercell_displacements"])
            self._mode_label_overrides = dict(meta.get("mode_labels") or {})
        else:
            self.mode_displacements_sc = {}
            self._mode_label_overrides = {}

        # 记录路径与模式，供 Distortion Page 使用
        self.phase_path = PhasePath.from_subgroup(
            parent_sg, result.subgroup, types
        )
        self.phase_path.k_vector = self._resolve_k_vector(
            result.subgroup.k_point_label,
            list(result.subgroup.k_parameters or []),
        )
        self.phase_path.validate()
        self._selected_subgroup = result.subgroup
        self.distortion_modes = self._compute_scoped_modes(
            parent_sg, result.subgroup, types, raw_modes=result.modes
        )
        self.mode_displacements = self._dist_mapper.map_modes_to_atoms(
            self.structure,
            self.symmetry_info["wyckoff_sites"],
            self.distortion_modes,
        )
        self._install_special_bush_supercell_modes(result.subgroup)
        self._sync_parametric_session_keys()
        print(t("method2.result", idx=subgroup_idx,
                n=len(self.mode_displacements) + len(self.mode_occupancies)))
        return result

    def search_method_3(self,
                        distortion_types: str | list[str] | None = None,
                        point_group: str | None = None,
                        space_group_type: int | None = None,
                        supercell_basis: list[list[str | int | float]] | None = None,
                        direct_sublattice_centering: str | None = None,
                        lattice_type: str = "direct",
                        generate_if_missing: bool = False,
                        include_affine_only_diagnostics: bool = False):
        """
        Method 3: Search over arbitrary k points for a specified point group and supercell.

        若 point_group 与 space_group_type 同时提供，按官网规则优先采用
        space_group_type。lattice_type 为官网 radio（direct/reciprocal）；
        本地引擎暂不支持 reciprocal（倒易超格）模式，会给出明确错误。

        Default 按目标空间群的默认 Bravais centering 解释；显式
        P/A/B/C/I/F/R 分别按所选 centering 解释为精确 primitive lattice。
        非母相原胞子格会额外按公度条件推断参数 k（如 LD ``g=1/6``）并枚举；
        ``generate_if_missing`` 与 Method 2 的 GenDB 开关相同。
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
            parent_rotations=[rotation.tolist() for rotation in self._parent_rotations()],
            parent_structure=self.structure,
            symmetry_tolerance=self.cfg.symmetry_cartesian_tolerance_angstrom,
            symmetry_angle_tolerance_degrees=(
                self.cfg.symmetry_angle_tolerance_degrees
            ),
            affine_exact_tolerance=(
                self.cfg.affine_exact_cartesian_tolerance_angstrom
            ),
            fractional_coordinate_tolerance=(
                self.cfg.fractional_coordinate_tolerance
            ),
            generate_if_missing=generate_if_missing,
            include_affine_only_diagnostics=include_affine_only_diagnostics,
        )
        parent_sg = self.symmetry_info["space_group_number"]
        result = self._search.method_3_search(parent_sg, query)

        result = self._filter_method3_routes_by_types(
            result, query.distortion_types or self.distortion_types
        )

        # 记录过滤后的候选，供 Method 2（search_method_2）使用（与 Method 1 一致）
        self._unresolved_method3_embedding_keys = set()
        for item in result:
            resolution = getattr(item, "route_resolution", "known_single_ir")
            item.subgroup._method3_route_resolution = resolution
            embedding_id = str(getattr(item, "embedding_id", "") or "").strip()
            if embedding_id:
                item.subgroup._method3_embedding_id = embedding_id
            if resolution == "affine_only_unresolved_coupled_route":
                key = self._method3_embedding_guard_key(item.subgroup)
                self._unresolved_method3_embedding_keys.add(key)
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
        畸变结构可以是母相的超胞（官网 Method 4 的常规情形）：原子数
        不一致时自动把母相与模式位移提升到畸变结构的超胞坐标系再分解。
        """
        if self.structure is None:
            raise RuntimeError("请先加载母相结构 (load_structure)")
        if not self.mode_displacements:
            raise RuntimeError("请先通过 select_path 或 search_method_2 计算模式")

        distorted_structure = read_cif(distorted_cif_path)
        parent = self.structure
        mode_disp = {k: v["displacements"] for k, v in self.mode_displacements.items()}
        if len(distorted_structure) != len(parent):
            # 超胞畸变：母相与模式位移提升到超胞坐标系（原实现直接报错，
            # 导致任何超胞畸变都无法分解——与官网行为不一致）。
            basis = self._resolve_distorted_supercell_basis(distorted_structure)
            k_vec = (self.phase_path.k_vector
                     if self.phase_path is not None else None)
            parent, mode_disp = self._dist_engine.lift_mode_displacements(
                parent, basis, mode_disp, k_vector=k_vec
            )
        query = Method4Query(
            atom_matching_method=atom_matching_method,
            robust_distance_threshold=robust_distance_threshold,
            provided_origin_shift=provided_origin_shift,
        )

        result = self._search.method_4_decompose(
            parent,
            distorted_structure,
            mode_disp,
            query,
        )

        print(t("method4.result", n=len(result.amplitudes), rms=result.rms_residual))
        return result

    def _resolve_distorted_supercell_basis(self, distorted: Structure
                                           ) -> np.ndarray:
        """确定母相 -> 畸变结构的超胞基矢 B（Ld = B @ Lp，原子数比 = |det B|）。

        优先采用当前相变路径的子群基矢（本地生成畸变结构的标准情形）；
        否则从两个晶格矩阵反推（B = Ld @ Lp⁻¹，应为近整数矩阵）。
        两者都无法给出一致的整数超胞关系时，报出明确错误。
        """
        n_parent = len(self.structure)
        n_dist = len(distorted)
        candidates: list[np.ndarray] = []
        if self.phase_path is not None and self.phase_path.basis_vectors:
            candidates.append(
                np.asarray(self.phase_path.supercell_basis(), dtype=float))
        lp = np.asarray(self.structure.lattice.matrix, dtype=float)
        ld = np.asarray(distorted.lattice.matrix, dtype=float)
        try:
            candidates.append(ld @ np.linalg.inv(lp))
        except np.linalg.LinAlgError:
            pass
        for b in candidates:
            if b.shape != (3, 3):
                continue
            if not np.allclose(b, np.round(b), atol=1e-4):
                continue
            det = abs(round(float(np.linalg.det(b))))
            if det >= 1 and n_parent * det == n_dist:
                return np.round(b).astype(int)
        raise ValueError(
            "畸变结构与母相原子数不一致，且无法确定超胞关系"
            f"（母相 {n_parent} 原子，畸变 {n_dist} 原子）；"
            "请确认畸变结构是当前母相/所选子群的超胞。")
