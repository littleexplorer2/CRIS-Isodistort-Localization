"""Search workflows aligned with ISODISTORT Method 1-4（基于真实 iso 枚举）。"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.groups import SpaceGroup

from ..backend import DistortionMode, IsoWrapper, KPointInfo, SubgroupInfo
from ..data.kpoints_official import (
    KPOINT_OFFICIAL,
    _param_coeff_in_component,
    official_kparams_to_iso,
    official_special_k_coords,
)
from .phase_path import normalize_distortion_types

CRYSTAL_SYSTEMS = {
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic",
}


def _sg_to_crystal_system(space_group_number: int) -> str:
    if 1 <= space_group_number <= 2:
        return "triclinic"
    if 3 <= space_group_number <= 15:
        return "monoclinic"
    if 16 <= space_group_number <= 74:
        return "orthorhombic"
    if 75 <= space_group_number <= 142:
        return "tetragonal"
    if 143 <= space_group_number <= 167:
        return "trigonal"
    if 168 <= space_group_number <= 194:
        return "hexagonal"
    if 195 <= space_group_number <= 230:
        return "cubic"
    return "unknown"


def _to_float(value: str | float) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if "/" in text:
        return float(Fraction(text))
    return float(text)


def _parse_basis_rows(rows: Sequence[Sequence[str | int | float]]) -> list[list[float]]:
    basis = []
    for row in rows:
        if len(row) != 3:
            raise ValueError("Each basis row must contain 3 values")
        basis.append([_to_float(v) for v in row])
    if len(basis) != 3:
        raise ValueError("Basis matrix must be 3x3")
    return basis


def _lattice_equivalent(a: Sequence[Sequence[float]],
                        b: Sequence[Sequence[float]]) -> bool:
    """两个 3x3 超胞基矢是否生成同一格点（GL(3,Z) 幺模等价）。

    Method 3 在用户指定**非恒等**目标子格时用此判定（对齐官网「指定格子」）。
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


def _is_identity_basis(basis: Sequence[Sequence[float]], atol: float = 1e-6) -> bool:
    arr = np.asarray(basis, dtype=float)
    return arr.shape == (3, 3) and bool(np.allclose(arr, np.eye(3), atol=atol))


def _integer_basis_matrix(basis: Sequence[Sequence[float]]) -> np.ndarray | None:
    """近整数的 3x3 超胞矩阵；否则 None。"""
    arr = np.asarray(basis, dtype=float)
    if arr.shape != (3, 3):
        return None
    rounded = np.rint(arr)
    if not np.allclose(arr, rounded, atol=1e-5):
        return None
    if abs(float(np.linalg.det(rounded))) < 1e-6:
        return None
    return rounded.astype(int)


def _validate_centering(value: str | None) -> str:
    """校验 Method 3 的 direct_sublattice_centering（官网 d/P/A/B/C/I/F/R radio）。

    官网帮助：未选带心时按所选空间群默认带心，点群则按 primitive。
    表单选 **P** 表示 primitive / no centering（结果摘要常写 ``no centering``），
    与 Default(``d``) 在本地均可接受（不另建带心数据库）。

    A/B/C/I/F/R 表示以该心化解读所填基矢；本地 iso 无法按任意带心再生成
    Method 3 子群库，故明确报错。
    """
    if value is None or str(value).strip() == "":
        return "d"
    key = str(value).strip().lower()
    if key in ("d", "default", "p"):
        return "d"
    raise ValueError(
        f"本地引擎 Method 3 仅支持 Default (d) 或 P（primitive / no centering），收到 {value!r}。"
        "A/B/C/I/F/R 带心需要官网按指定子格在线生成子群数据库。 "
        f"/ Local Method 3 supports Default (d) or P (primitive/no centering); got {value!r}. "
        "A/B/C/I/F/R centering needs the website's on-demand subgroup database."
    )


def _fraction_str(value: Fraction) -> str:
    value = Fraction(value).limit_denominator(10_000)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _param_value_candidates(M: np.ndarray) -> list[Fraction]:
    """由超胞倍数推断一维参数 k 的候选（如 c'=6 → 1/6、1/3、1/2）。"""
    cands: set[Fraction] = set()
    for v in np.abs(np.diag(M)).astype(int).tolist():
        if v > 1:
            cands.add(Fraction(1, int(v)))
    det = abs(int(round(float(np.linalg.det(M)))))
    for d in range(2, det + 1):
        if det % d == 0:
            cands.add(Fraction(1, d))
    return sorted(cands, key=lambda f: (f.denominator, f.numerator))


def _k_compatible_with_supercell(kvec: np.ndarray, M: np.ndarray) -> bool:
    """k 是否为超胞倒格点（M^T k ∈ ℤ³），即公度锁定在该子格上。"""
    v = np.asarray(M, dtype=float).T @ np.asarray(kvec, dtype=float)
    return bool(np.allclose(v, np.round(v), atol=1e-5))


def _eval_simple_k_component(comp: str, param_name: str, param_val: Fraction) -> Fraction | None:
    """解析 ``0`` / ``g`` / ``2a`` 形式；含 ``+``/``-`` 混杂的复杂式返回 None。"""
    text = str(comp).strip()
    if re.fullmatch(r"-?\d+(/\d+)?", text):
        return Fraction(text)
    coeff = _param_coeff_in_component(text, param_name)
    if coeff is None:
        return None
    if text not in (param_name, f"-{param_name}") and not re.fullmatch(
        rf"-?\d+{re.escape(param_name)}", text
    ):
        return None
    return coeff * param_val


def _kvec_from_template(
    coords: Sequence[str], param_name: str, param_val: Fraction
) -> np.ndarray | None:
    comps: list[float] = []
    for comp in coords:
        val = _eval_simple_k_component(comp, param_name, param_val)
        if val is None:
            return None
        comps.append(float(val))
    return np.asarray(comps, dtype=float)


def _line_kpoint_template(parent_sg: int, kp: KPointInfo) -> tuple[list[str], str] | None:
    """单参数线 k 点的坐标模板与参数名（优先官网表）。"""
    override = KPOINT_OFFICIAL.get(int(parent_sg), {}).get(kp.label.strip())
    if override is not None:
        _kov, coords, params = override
        if len(params) == 1:
            return [str(c) for c in coords], str(params[0])
        return None
    params = list(kp.parameters or [])
    if len(params) != 1:
        return None
    return [str(c) for c in (kp.coordinates or [])], str(params[0])


def _space_group_to_point_group(space_group_number: int) -> str:
    try:
        return SpaceGroup.from_int_number(space_group_number).point_group
    except (ValueError, KeyError, TypeError):
        return ""


def _basis_is_sublattice_of(basis: Sequence[Sequence[float]],
                            sublattice: Sequence[Sequence[float]],
                            parent_rotations: Sequence[np.ndarray] | None = None
                            ) -> bool:
    """
    判断子群超胞基矢 B 的格点是否是被选子格 S 的子格。

    对应官网 Method 1 的 direct sublattice / Conventional lattice /
    Primitive lattice 过滤：B 的每一行必须是 S 的整系数线性组合，
    即 N = B @ inv(S) 的元素全部为整数（S 可为对角阵或任意 3x3 矩阵）。

    官网同一 lattice 选项还包含母相点群任意旋转得到的等价子格；若提供
    ``parent_rotations``，则对 S' = S @ R 与 R @ S 一并判定。
    """
    b = np.asarray(basis, dtype=float)
    s0 = np.asarray(sublattice, dtype=float)
    if b.shape != (3, 3) or s0.shape != (3, 3):
        return False

    candidates: list[np.ndarray] = [s0]
    if parent_rotations:
        for rot in parent_rotations:
            r = np.asarray(rot, dtype=float)
            if r.shape != (3, 3):
                continue
            candidates.append(s0 @ r)
            candidates.append(r @ s0)

    for s in candidates:
        try:
            n = b @ np.linalg.inv(s)
        except np.linalg.LinAlgError:
            continue
        if np.allclose(n, np.round(n), atol=1e-6):
            return True
    return False


@dataclass
class Method1Query:
    """Method 1: search over all special k points."""

    distortion_types: str | Sequence[str] | None = None
    crystal_system: str | Sequence[str] | None = None  # 单个或列表（多选=OR）
    subgroup_space_group: int | None = None
    lattice: Sequence[Sequence[float]] | None = None  # 官网 conventional/primitive lattice（3x3 子格矩阵）
    maximal_subgroup_only: bool = False
    # 母相点群旋转（分数坐标）；lattice 过滤时与官网一样合并点群轨道
    parent_rotations: Sequence[Sequence[Sequence[float]]] | None = None


@dataclass
class Method1ResultItem:
    subgroup: SubgroupInfo
    crystal_system: str
    is_maximal: bool


@dataclass
class Method2Query:
    """Method 2: general method over a selected subgroup (k/IR/OPD 由子群自身携带)。"""

    subgroup_idx: int
    distortion_type: str | Sequence[str] = "displacive"
    number_of_independent_modulations: int = 0  # 仅支持 0（公度）；非 0 报错


@dataclass
class Method2Result:
    subgroup: SubgroupInfo
    modes: list[DistortionMode]
    metadata: dict[str, object]


@dataclass
class Method3Query:
    """Method 3: search over arbitrary k for point/space group + supercell."""

    distortion_types: str | Sequence[str] | None = None
    point_group: str | None = None
    space_group_type: int | None = None
    supercell_basis: Sequence[Sequence[str | int | float]] | None = None
    direct_sublattice_centering: str | None = None
    lattice_type: str = "direct"  # 官网 radio：direct（实空间子格）/ reciprocal（倒易超格）
    # 参数 k 点子群库缺失时是否在线生成（与 Method 2 GenDB 同一开关）
    generate_if_missing: bool = False


@dataclass
class Method3ResultItem:
    subgroup: SubgroupInfo
    point_group: str
    basis: list[list[float]]


@dataclass
class Method4Query:
    """Method 4: mode decomposition of distorted structure."""

    atom_matching_method: str = "nearest-site"
    robust_distance_threshold: float = 0.25
    provided_origin_shift: Sequence[float] | None = None


@dataclass
class Method4Result:
    amplitudes: dict[str, float]
    rms_residual: float
    max_abs_residual: float
    assignments: list[int]
    metadata: dict[str, object]


class IsoSearchEngine:
    """Implements local Method 1-4 search workflows on top of the real iso binary."""

    def __init__(self, iso_wrapper: IsoWrapper) -> None:
        self._iso = iso_wrapper

    # ----------------------------------------------------------------
    # Method 1：全特殊 k 点搜索 + 客户端过滤（与官网逻辑 AND 语义一致）
    # ----------------------------------------------------------------

    def method_1_search(self, parent_sg: int, query: Method1Query,
                        distortion_types=None,
                        subgroups: Sequence[SubgroupInfo] | None = None
                        ) -> list[Method1ResultItem]:
        """
        官网 Method 1：遍历全部特殊 k 点，得到子群候选后按用户条件过滤。

        过滤条件（多条件同时生效，逻辑 AND；同一条件内多选为 OR）：
        - crystal system：子群所属晶系（单个或列表，列表任中其一即通过）
        - subgroup space group：子群空间群号
        - maximal subgroup only：仅保留 maximal 子群
        - lattice：超胞格是否为所选子格的子格
          （官网 Conventional lattice / Primitive lattice）

        Args:
            parent_sg: 母相空间群号
            query: Method 1 查询
            distortion_types: 畸变类型（保留参数；类型过滤在模式计算阶段执行）
            subgroups: 预枚举的子群候选（缓存复用）；None 时现场枚举
        """
        _ = distortion_types
        if subgroups is None:
            subgroups = self._iso.enumerate_all_special_subgroups(
                parent_sg, query.distortion_types
            )

        crystal_systems: set[str] | None = None
        if query.crystal_system:
            raw = query.crystal_system
            if isinstance(raw, str):
                raw = [raw]
            crystal_systems = {x.strip().lower() for x in raw if x.strip()}

        result: list[Method1ResultItem] = []
        for sg in subgroups:
            crystal_system = _sg_to_crystal_system(sg.space_group_number)
            item = Method1ResultItem(
                subgroup=sg,
                crystal_system=crystal_system,
                is_maximal=sg.is_maximal,
            )
            if crystal_systems and crystal_system not in crystal_systems:
                continue
            if query.subgroup_space_group and sg.space_group_number != query.subgroup_space_group:
                continue
            if query.maximal_subgroup_only and not sg.is_maximal:
                continue
            if query.lattice is not None:
                rots = None
                if query.parent_rotations:
                    rots = [
                        np.asarray(r, dtype=float) for r in query.parent_rotations
                    ]
                if not _basis_is_sublattice_of(
                    sg.basis_vectors, query.lattice, parent_rotations=rots
                ):
                    continue
            result.append(item)
        return result

    # ----------------------------------------------------------------
    # Method 2：指定 k 点/IR/OPD 的模式计算
    # ----------------------------------------------------------------

    def method_2_search(self, parent_sg: int, subgroups: Sequence[SubgroupInfo],
                        query: Method2Query,
                        wyckoff_letters: Sequence[str] | None = None) -> Method2Result:
        """
        官网 Method 2：在已枚举子群中按序号选择目标子群，计算其畸变模式。

        模式基矢由真实 iso 的 DISPLAY BUSH 计算（需要母相 Wyckoff 位置）。

        Args:
            parent_sg: 母相空间群号
            subgroups: 子群候选列表（来自 Method 1 或 list_subgroups）
            query: Method 2 查询参数（subgroup_idx 必填）
            wyckoff_letters: 母相结构各原子的 Wyckoff 位置字母

        Returns:
            Method2Result
        """
        target = next((s for s in subgroups if s.index == query.subgroup_idx), None)
        if target is None:
            raise ValueError(
                f"Subgroup index {query.subgroup_idx} not found; "
                "请先执行 Method 1 或 list_subgroups 获得候选列表"
            )

        # 官网 nmod（# of independent incommensurate modulations）仅对非公度
        # （参数 k 点）有意义；本地引擎只支持公度特殊 k 点（nmod=0），
        # 非零值明确报错而非静默忽略。
        if query.number_of_independent_modulations:
            raise ValueError(
                "本地引擎不支持非公度调制叠加"
                "（number_of_independent_modulations 必须为 0）；"
                "该参数仅对官网参数 k 点的 (3+d) 维超空间机制有意义，"
                "本地 iso 二进制无法完成。"
                "/ The local engine does not support incommensurate modulation superposition "
                "(number_of_independent_modulations must be 0)."
            )

        # wyckoff_letters 为 None 表示调用方未提供（误用）；空列表表示
        # 作用域内无 Wyckoff 位置（如全部类型选 none）-> 直接返回空模式
        if wyckoff_letters is None:
            raise ValueError(
                "Method 2 计算模式需要母相结构的 Wyckoff 位置信息，"
                "请先加载结构（load_structure）"
            )
        modes = []
        if wyckoff_letters:
            modes = self._iso.calc_distortion_modes(
                parent_sg, target, wyckoff_letters=wyckoff_letters
            )

        metadata = {
            # k 点 / IR / OPD 来自所选子群自身（查询时已确定，无需重复传参）
            "k_point_label": target.k_point_label,
            "irrep_label": target.irrep_label,
            "opd_symbol": target.opd_symbol,
            "k_parameters": list(target.k_parameters),
            "number_of_independent_modulations": query.number_of_independent_modulations,
        }

        return Method2Result(subgroup=target, modes=modes, metadata=metadata)

    # ----------------------------------------------------------------
    # Method 3：指定点群/空间群 + 超胞
    # ----------------------------------------------------------------

    def method_3_search(self, parent_sg: int, query: Method3Query) -> list[Method3ResultItem]:
        """
        官网 Method 3 的本地实现（特殊 k + 由超胞推断的公度参数 k）：

        - 若同时提供 point_group 与 space_group_type，空间群选择优先；
        - 带心：Default(``d``) 与 **P**（primitive / no centering）可接受；
          A/B/C/I/F/R 明确报错；
        - 基矢：恒等（网页默认）用子格包容过滤（便于浏览）；非恒等目标子格
          用 ``_lattice_equivalent``（对齐官网指定格子）；
        - 非恒等整数超胞时，由 M 推断公度参数 k（如 ``(0,0,6)`` → LD ``g=1/6``），
          再枚举该 k 上各 IR 子群（可 ``generate_if_missing``），以覆盖官网
          Method 3 对参数 k 的搜索（如 EuAl4 → 99 P4mm, s=12, i=24）。

        已知限制：reciprocal 模式不支持；多参数平面/一般 k（GP）未自动推断；
        参数 k 点位移模式仍不能本地计算（与 Method 2 相同）。
        """
        distortion_types = normalize_distortion_types(query.distortion_types)
        _validate_centering(query.direct_sublattice_centering)

        if query.point_group and query.space_group_type:
            point_group_filter = None
        else:
            point_group_filter = query.point_group

        basis: list[list[float]] | None = None
        if query.supercell_basis:
            basis = _parse_basis_rows(query.supercell_basis)

        # 1) 特殊 k 点候选（Method 1 同源）
        subgroups = list(
            self._iso.enumerate_all_special_subgroups(parent_sg, distortion_types)
        )

        # 2) 非恒等超胞 → 公度参数 k 回退（官网 Method 3 的任意 k 近似）
        if basis is not None and not _is_identity_basis(basis):
            parametric = self._method3_parametric_subgroups(
                parent_sg,
                basis,
                distortion_types,
                generate_if_missing=bool(query.generate_if_missing),
            )
            subgroups.extend(parametric)

        matched: list[SubgroupInfo] = []
        for sg in subgroups:
            point_group = _space_group_to_point_group(sg.space_group_number)
            if query.space_group_type and sg.space_group_number != query.space_group_type:
                continue
            if point_group_filter and point_group != point_group_filter:
                continue
            if basis is not None:
                sg_basis = sg.basis_vectors or [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                if _is_identity_basis(basis):
                    if not _basis_is_sublattice_of(sg_basis, basis):
                        continue
                elif not _lattice_equivalent(sg_basis, basis):
                    continue
            matched.append(sg)

        result: list[Method3ResultItem] = []
        for i, sg in enumerate(matched):
            sg.index = i
            result.append(Method3ResultItem(
                subgroup=sg,
                point_group=_space_group_to_point_group(sg.space_group_number),
                basis=[list(row) for row in (sg.basis_vectors or [])],
            ))
        return result

    def _method3_parametric_subgroups(
        self,
        parent_sg: int,
        basis: Sequence[Sequence[float]],
        distortion_types,
        generate_if_missing: bool = False,
    ) -> list[SubgroupInfo]:
        """由超胞基矢推断公度线 k 点并枚举子群（Method 3 参数 k 回退）。"""
        M = _integer_basis_matrix(basis)
        if M is None or abs(int(round(float(np.linalg.det(M))))) <= 1:
            return []

        try:
            kpoints = self._iso.list_k_points(parent_sg)
        except Exception:  # noqa: BLE001 - Method 3 回退失败时仍返回特殊 k 结果
            return []

        param_cands = _param_value_candidates(M)
        searches: list[tuple[KPointInfo, list[str], list[str]]] = []
        seen: set[tuple[str, tuple[str, ...]]] = set()
        for kp in kpoints:
            if kp.is_special:
                continue
            tmpl = _line_kpoint_template(parent_sg, kp)
            if tmpl is None:
                continue
            coords, pname = tmpl
            for pval in param_cands:
                kvec = _kvec_from_template(coords, pname, pval)
                if kvec is None or not _k_compatible_with_supercell(kvec, M):
                    continue
                if np.allclose(kvec, 0.0):
                    continue
                official = [_fraction_str(pval)]
                key = (kp.label, tuple(official))
                if key in seen:
                    continue
                seen.add(key)
                iso_vals = official_kparams_to_iso(parent_sg, kp.label, official, kp)
                searches.append((kp, official, iso_vals))

        out: list[SubgroupInfo] = []
        for kp, official_params, iso_params in searches:
            try:
                irreps = self._iso.list_irreps(
                    parent_sg, kp.label, k_parameters=iso_params
                )
            except Exception:  # noqa: BLE001
                continue
            irreps = [
                ir for ir in irreps
                if self._iso._include_irrep(ir, distortion_types)
            ]
            for ir in irreps:
                try:
                    batch = self._iso.list_subgroups(
                        parent_sg,
                        kp.label,
                        ir.label,
                        k_parameters=iso_params,
                        generate_if_missing=generate_if_missing,
                        start_index=len(out),
                    )
                except Exception:  # noqa: BLE001
                    continue
                for sg in batch:
                    sg.k_parameters = list(official_params)
                    sg.k_coordinates = official_special_k_coords(
                        parent_sg, kp.label, None, official_params
                    )
                    out.append(sg)
        return out

    # ----------------------------------------------------------------
    # Method 4：模式分解（自研最小二乘拟合）
    # ----------------------------------------------------------------

    def method_4_decompose(self,
                           parent_structure: Structure,
                           distorted_structure: Structure,
                           mode_displacements: dict[str, np.ndarray],
                           query: Method4Query) -> Method4Result:
        if len(parent_structure) != len(distorted_structure):
            raise ValueError(
                "Mode decomposition currently requires parent and distorted structures "
                "to have the same atom count"
            )

        assignments = self._match_atoms(parent_structure, distorted_structure, query)
        delta = self._build_delta_vector(parent_structure, distorted_structure, assignments)

        mode_labels = list(mode_displacements.keys())
        if not mode_labels:
            raise ValueError("No mode displacements are available for decomposition")

        columns = []
        for label in mode_labels:
            vec = np.asarray(mode_displacements[label], dtype=float).reshape(-1)
            if vec.size != delta.size:
                raise ValueError(
                    f"Mode {label} has incompatible size {vec.size}, expected {delta.size}"
                )
            columns.append(vec)
        a_matrix = np.column_stack(columns)

        coeffs, _, _, _ = np.linalg.lstsq(a_matrix, delta, rcond=None)
        reconstructed = a_matrix @ coeffs
        residual = delta - reconstructed

        amplitudes = {mode_labels[i]: float(coeffs[i]) for i in range(len(mode_labels))}
        rms = float(np.sqrt(np.mean(residual ** 2)))
        max_abs = float(np.max(np.abs(residual)))

        return Method4Result(
            amplitudes=amplitudes,
            rms_residual=rms,
            max_abs_residual=max_abs,
            assignments=assignments,
            metadata={
                "atom_matching_method": query.atom_matching_method,
                "provided_origin_shift": list(query.provided_origin_shift)
                if query.provided_origin_shift is not None else None,
            },
        )

    @staticmethod
    def _fractional_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        delta = b - a
        delta -= np.round(delta)
        return delta

    def _match_atoms(self,
                     parent_structure: Structure,
                     distorted_structure: Structure,
                     query: Method4Query) -> list[int]:
        if query.atom_matching_method not in {"nearest-site", "robust"}:
            raise ValueError("atom_matching_method must be 'nearest-site' or 'robust'")

        assignments: list[int] = []
        used: set[int] = set()

        for i, site in enumerate(parent_structure):
            species = site.species_string
            parent_coord = np.asarray(site.frac_coords)
            candidates = []
            for j, dst in enumerate(distorted_structure):
                if j in used:
                    continue
                if dst.species_string != species:
                    continue
                delta = self._fractional_delta(parent_coord, np.asarray(dst.frac_coords))
                dist = float(np.linalg.norm(delta))
                if query.atom_matching_method == "nearest-site":
                    candidates.append((dist, j))
                else:
                    if dist <= query.robust_distance_threshold:
                        candidates.append((dist, j))

            if not candidates:
                raise ValueError(f"Cannot match parent atom index {i} ({species})")

            candidates.sort(key=lambda t: t[0])
            chosen = candidates[0][1]
            assignments.append(chosen)
            used.add(chosen)

        return assignments

    def _build_delta_vector(self,
                            parent_structure: Structure,
                            distorted_structure: Structure,
                            assignments: Sequence[int]) -> np.ndarray:
        delta_rows = []
        for i, j in enumerate(assignments):
            p = np.asarray(parent_structure[i].frac_coords)
            d = np.asarray(distorted_structure[j].frac_coords)
            delta_rows.append(self._fractional_delta(p, d))
        return np.asarray(delta_rows, dtype=float).reshape(-1)
