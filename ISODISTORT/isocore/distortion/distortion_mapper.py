"""
畸变模式映射器 - 将 iso（DISPLAY BUSH）输出的模式基矢映射到每个原子的位移向量

对应阶段四，步骤8：畸变基矢到原子坐标的映射

映射规则（基于实测 DISPLAY BUSH 输出）：
- BUSH 每个位移列是子群固定子空间的一条独立基矢。
  ``IsoWrapper`` 已按列拆分，因此这里的每个 ``BushMode`` 恰含一个向量。
- 若某位点只有一个代表原子行（常见于 Gamma 点均匀模式），则将该位移
  均匀作用于该位点的全部等效原子。
- 若某位点有多个代表原子行，每个代表对应模式分裂出的一个子轨道：
  - 先按**结构真实坐标**解析代表点中的自由参数（如 ``(0,0,z)`` 中的 z
    取该位点原子的实际坐标，否则两个相反代表点会被解析成同一点）；
  - 再按“模母相格点等价”（格点平移副本取同一位移，如带心平移）把位点
    原子分配到各代表；无法匹配的原子按周期最近邻回退。

已知差异（见 README）：
- 位移向量为超胞分数坐标单位、最大分量为 1 的约定（与官网一致），
  但尚未实现官网 As/Ap 振幅的 normfactor 归一化换算。
"""

import re
from collections.abc import Sequence
from fractions import Fraction

import numpy as np
from pymatgen.core import Lattice, Structure

from ..backend import BushMode, DistortionMode
from ..structure import build_supercell


class DistortionMapper:
    """
    畸变模式 -> 原子位移 映射器

    负责将 iso 计算出的抽象模式基矢，转化为每个原子的实际位移向量。
    """

    def __init__(self) -> None:
        self._prim_matrix: np.ndarray | None = None
        self._conv_matrix: np.ndarray | None = None

    def map_modes_to_atoms(self, structure: Structure,
                           wyckoff_sites: list[dict],
                           modes: list[DistortionMode]) -> dict:
        """
        将畸变模式映射到结构中每个原子

        Args:
            structure: 母相晶体结构
            wyckoff_sites: 结构的 Wyckoff 位置分组信息（来自 SymmetryValidator）
            modes: iso 计算出的畸变模式列表（含 bush_modes）

        Returns:
            dict: 每个模式对应的原子位移向量
                {
                    "irrep_label": {
                        "mode": DistortionMode,
                        "displacements": np.ndarray (N_atoms x 3),
                    }
                }
        """
        n_atoms = len(structure)
        self._lattice_setup(structure)

        # 按 (字母, 物种) 分组：同一字母可能被多个物种共用（如 Fe/O 同在 k）
        groups: list[tuple[str, list[int]]] = []
        for site_info in wyckoff_sites:
            groups.append((site_info["wyckoff_letter"],
                           list(site_info["equivalent_indices"])))

        result: dict[str, dict] = {}
        for mode in modes:
            # 按位点字母分组 BUSH 行；多位点时拆成独立幅度键（对齐官网
            # 每个 [Site:letter:dsp]… 各算一个 displacive mode）。
            by_letter: dict[str, list[BushMode]] = {}
            for bush in mode.bush_modes:
                by_letter.setdefault(bush.wyckoff_letter, []).append(bush)

            letter_disps: dict[str, np.ndarray] = {}
            for letter, indices in groups:
                bushes = by_letter.get(letter, [])
                if not bushes or not indices:
                    continue
                displacements = np.zeros((n_atoms, 3))

                # 解析每个代表点（含自由参数 -> 结构真实坐标）。
                # 列在 wrapper 层已拆分；如果再收到多列数据，说明上游
                # 违反模式身份约定，必须显式报错而不是静默相加。
                reps: list[tuple[np.ndarray, np.ndarray]] = []
                for bush in bushes:
                    if not bush.displacements:
                        continue
                    if len(bush.displacements) != 1:
                        raise ValueError(
                            "DISPLAY BUSH displacement columns must be split "
                            "into independent DistortionMode objects before mapping"
                        )
                    point = self._resolve_rep_point(bush, structure, indices)
                    vec = np.asarray(bush.displacements[0], dtype=float)
                    reps.append((point, vec))
                if not reps:
                    continue

                # 只保留能“模格点”匹配到本组原子的代表（同一字母不同物种时，
                # 彼此的符号点取位约定可能不同，无法匹配的物种保持不动，
                # 避免错误回退破坏对称性）
                matched_reps = [
                    (p, v) for p, v in reps
                    if any(self._lattice_equiv(
                        np.asarray(structure[i].frac_coords), p) for i in indices)
                ]
                if not matched_reps:
                    continue  # 本组原子保持不动
                if len(matched_reps) == 1:
                    # 均匀模式：同一位移作用于该位点的全部等效原子
                    disp = matched_reps[0][1]
                    for idx in indices:
                        displacements[idx] = disp
                else:
                    # 多代表：优先按“周期等价”（模整数平移）把代表点分配给
                    # 原子。非 Γ k 点模式中，BUSH 会为带心副本给出相位反号的
                    # 独立代表（如 Im-3m H4-：`(0,0,0)→(1,0,1)` 与
                    # `(-1/2,1/2,1/2)→(-1,0,-1)`），这两个代表相差原始格点
                    # 平移、互为“模格点等价”，若按旧逻辑“首个代表胜出”会把
                    # 全部原子分给同一代表、丢失 k 点相位（畸变退化为刚性平移）。
                    # 代表点与原子分数坐标 mod 1 重合时即视为同一原子，直接
                    # 采用该代表的位移；无法周期匹配时回退“模母相格点等价”
                    # （子轨道/带心副本同位移，Γ 点均匀语义），最后退周期最近邻。
                    for idx in indices:
                        coord = np.asarray(structure[idx].frac_coords)
                        assigned = False
                        for point, disp in matched_reps:
                            if self._periodic_equiv(coord, point):
                                displacements[idx] = disp
                                assigned = True
                                break
                        if not assigned:
                            for point, disp in matched_reps:
                                if self._lattice_equiv(coord, point):
                                    displacements[idx] = disp
                                    assigned = True
                                    break
                        if not assigned:
                            # 回退：周期最近邻（仅本组代表内）
                            best_disp = matched_reps[0][1]
                            best_dist = float("inf")
                            for point, disp in matched_reps:
                                d = self._periodic_distance(coord, point)
                                if d < best_dist:
                                    best_dist = d
                                    best_disp = disp
                            displacements[idx] = best_disp

                max_norm = float(np.max(np.linalg.norm(displacements, axis=1)))
                if max_norm > 1e-12:
                    letter_disps[letter] = displacements / max_norm

            if not letter_disps:
                continue
            unique = str(getattr(mode, "amplitude_key", "") or "")
            if unique:
                letter, displacements = next(iter(letter_disps.items()))
                result[unique] = {
                    "mode": mode,
                    "displacements": displacements,
                    "wyckoff_letter": letter,
                }
                continue
            if len(letter_disps) == 1:
                letter, displacements = next(iter(letter_disps.items()))
                result[mode.irrep_label] = {
                    "mode": mode,
                    "displacements": displacements,
                    "wyckoff_letter": letter,
                }
            else:
                for letter, displacements in letter_disps.items():
                    key = f"{mode.irrep_label}__{letter}"
                    result[key] = {
                        "mode": mode,
                        "displacements": displacements,
                        "wyckoff_letter": letter,
                    }

        return result

    def map_bush_modes_to_supercell(
        self,
        structure: Structure,
        wyckoff_sites: list[dict],
        modes: list[DistortionMode],
        basis_vectors: list[list[float]],
        *,
        cartesian_tolerance: float = 1e-4,
        subgroup_operations: Sequence[object] | None = None,
        subgroup_translation_lattice: Sequence[Sequence[float]] | None = None,
    ) -> dict[str, np.ndarray]:
        """Map complete DISPLAY BUSH basis vectors directly to the child cell.

        BUSH continuation rows already contain the relative signs/phases on
        every representative atom of the selected isotropy subgroup.  For a
        child basis ``B`` (rows in parent fractional coordinates), a child
        atom at parent coordinate ``x`` belongs to a BUSH representative
        ``p`` precisely when ``x-p`` belongs to the subgroup's *primitive*
        translation lattice ``T_H``.  Its arrow is transformed to the printed
        conventional child basis as ``du_child = du_parent @ inv(B)``.  The
        distinction matters for A/B/C/I/F/R centered child groups.

        This is deliberately different from reconstructing every mode with a
        cosine of the *primary* k vector: a single BUSH table can also contain
        secondary irreps at GM/X/N/etc., each with its own phase pattern.
        Missing or conflicting representatives are treated as invalid output
        rather than silently replaced by a guessed Bloch phase.
        """
        if not modes:
            return {}
        basis = np.asarray(basis_vectors, dtype=float)
        if basis.shape != (3, 3) or not np.all(np.isfinite(basis)):
            raise ValueError("subgroup basis must be a finite 3x3 matrix")
        det = float(np.linalg.det(basis))
        if abs(det) < 1e-12:
            raise ValueError("subgroup basis vectors must be linearly independent")
        if cartesian_tolerance <= 0 or not np.isfinite(cartesian_tolerance):
            raise ValueError("cartesian_tolerance must be finite and positive")

        child = build_supercell(structure, basis)
        inverse = np.linalg.inv(basis)
        translation_basis = np.asarray(
            subgroup_translation_lattice
            if subgroup_translation_lattice is not None
            else basis,
            dtype=float,
        )
        if (
            translation_basis.shape != (3, 3)
            or not np.all(np.isfinite(translation_basis))
            or abs(float(np.linalg.det(translation_basis))) < 1e-12
        ):
            raise ValueError(
                "subgroup primitive translation lattice must be a finite "
                "nonsingular 3x3 matrix"
            )
        translation_inverse = np.linalg.inv(translation_basis)
        parent_lattice = np.asarray(structure.lattice.matrix, dtype=float)
        translation_metric = Lattice(translation_basis @ parent_lattice)
        parent_coords = np.asarray(structure.frac_coords, dtype=float)
        child_as_parent = np.asarray(child.frac_coords, dtype=float) @ basis

        # Map each generated child atom back to one actual parent-cell site.
        # The species check prevents coincident special positions belonging to
        # different chemical orbits from being conflated.
        child_parent_index: list[int] = []
        for child_index, coord in enumerate(child_as_parent):
            best_index = -1
            best_distance = float("inf")
            for parent_index, parent_coord in enumerate(parent_coords):
                if child[child_index].species != structure[parent_index].species:
                    continue
                # Component-wise wrapping of fractional coordinates is not a
                # nearest-image construction in a skew cell.  The direct-
                # lattice metric must choose the integer image that minimizes
                # ``||(coord-parent_coord-n) B||``.
                distance = float(
                    structure.lattice.get_distance_and_image(
                        parent_coord, coord,
                    )[0]
                )
                if distance < best_distance:
                    best_index = parent_index
                    best_distance = distance
            if best_index < 0 or best_distance > cartesian_tolerance:
                raise ValueError(
                    "cannot map child atom to a chemically matching parent site: "
                    f"child index {child_index}, residual {best_distance:.6g} Å"
                )
            child_parent_index.append(best_index)

        parent_groups: list[tuple[str, set[int], list[int]]] = []
        for site_info in wyckoff_sites:
            indices = [int(value) for value in site_info["equivalent_indices"]]
            parent_groups.append((
                str(site_info["wyckoff_letter"]),
                set(indices),
                indices,
            ))

        operation_data: list[tuple[np.ndarray, np.ndarray]] = []
        for operation in subgroup_operations or ():
            rotation = np.asarray(operation.rotation, dtype=float)
            translation = np.asarray(
                operation.translation, dtype=float,
            )
            if rotation.shape != (3, 3) or translation.shape != (3,):
                raise ValueError("subgroup Seitz operations must contain R(3x3), t(3)")
            operation_data.append((rotation, translation))
        if not operation_data:
            operation_data.append((np.eye(3), np.zeros(3)))

        result: dict[str, np.ndarray] = {}
        for mode in modes:
            key = str(mode.amplitude_key or mode.irrep_label)
            if not key:
                raise ValueError("DISPLAY BUSH mode has no stable amplitude key")
            parent_arrows = np.zeros((len(child), 3), dtype=float)
            assigned = np.zeros(len(child), dtype=bool)
            target_letter = str(mode.wyckoff_site or "")

            for letter, parent_index_set, parent_indices in parent_groups:
                if target_letter and letter != target_letter:
                    continue
                bushes = [
                    bush for bush in mode.bush_modes
                    if str(bush.wyckoff_letter) == letter
                ]
                if not bushes:
                    continue
                reps: list[tuple[np.ndarray, np.ndarray]] = []
                for bush in bushes:
                    if len(bush.displacements) != 1:
                        raise ValueError(
                            "DISPLAY BUSH displacement columns must be split "
                            "before child-cell mapping"
                        )
                    point = self._resolve_rep_point(
                        bush, structure, parent_indices,
                    )
                    arrow = np.asarray(bush.displacements[0], dtype=float)
                    if point.shape != (3,) or arrow.shape != (3,):
                        raise ValueError("BUSH point and displacement must be 3-vectors")
                    for rotation, translation in operation_data:
                        # Seitz convention is x' = R x + t for column
                        # coordinates.  Points/arrows here are row vectors.
                        reps.append((
                            point @ rotation.T + translation,
                            arrow @ rotation.T,
                        ))

                group_children = [
                    index for index, parent_index in enumerate(child_parent_index)
                    if parent_index in parent_index_set
                ]
                for child_index in group_children:
                    coord = child_as_parent[child_index]
                    matches: list[np.ndarray] = []
                    for point, arrow in reps:
                        translation_delta = (coord - point) @ translation_inverse
                        # ``translation_delta`` is fractional in an arbitrary
                        # primitive basis of the subgroup translation lattice.
                        # Rounding its components is not a nearest-image
                        # construction when that basis is skew or unreduced;
                        # it would make BUSH coverage depend on the chosen
                        # (GL(3,Z)-equivalent) lattice basis.  Minimize in the
                        # actual direct-lattice metric instead.
                        distance = float(
                            translation_metric.get_distance_and_image(
                                translation_delta, np.zeros(3),
                            )[0]
                        )
                        if distance <= cartesian_tolerance:
                            matches.append(arrow)
                    if not matches:
                        raise ValueError(
                            "DISPLAY BUSH table does not cover child atom "
                            f"{child_index} in Wyckoff orbit {letter!r}"
                        )
                    reference = matches[0]
                    if any(
                        float(np.linalg.norm((other - reference) @ parent_lattice))
                        > cartesian_tolerance
                        for other in matches[1:]
                    ):
                        raise ValueError(
                            "DISPLAY BUSH representatives assign conflicting arrows "
                            f"to child atom {child_index}"
                        )
                    parent_arrows[child_index] = reference
                    assigned[child_index] = True

            if not np.any(assigned):
                raise ValueError(
                    f"DISPLAY BUSH mode {key!r} did not map to any child atom"
                )
            # Preserve the established UI amplitude convention: the largest
            # parent-fractional arrow norm is one.  The physically essential
            # direction/phase is exact before the child-coordinate transform.
            scale = float(np.max(np.linalg.norm(parent_arrows, axis=1)))
            if scale <= 1e-12:
                raise ValueError(f"DISPLAY BUSH mode {key!r} is identically zero")
            result[key] = (parent_arrows / scale) @ inverse
        return result

    # ----------------------------------------------------------------
    # 代表点解析 / 格点等价
    # ----------------------------------------------------------------

    def _lattice_setup(self, structure: Structure) -> None:
        """缓存母相惯用格矩阵与原胞格矩阵（模格点等价判定用）。

        用 structure.get_primitive_structure()（纯胞约化，与母相格点同一格）
        而非 get_primitive_standard_structure()（可能给出旋转后的不同格点，
        实测对 R3m 等带心结构会把带心平移判成非格点）。
        """
        prim = structure.get_primitive_structure()
        self._prim_matrix = np.asarray(prim.lattice.matrix, dtype=float)
        self._conv_matrix = np.asarray(structure.lattice.matrix, dtype=float)

    def _resolve_rep_point(self, bush: BushMode, structure: Structure,
                           indices: list[int]) -> np.ndarray:
        """把 BUSH 代表点解析为具体分数坐标。

        带自由参数的代表点（如 ``(0,0,z)``、``(x,2x,z)``、``(x+1/2,-y+1,-z)``）
        中，字母参数按其表示的坐标分量取值：x/y/z（或 a/b/g）分别取位点
        第一个原子的第 0/1/2 个分数坐标分量。否则相反符号的代表点
        （``(0,0,z)`` 与 ``(0,0,-z)``）会解析成同一点，导致相反位移的
        代表无法区分。
        """
        _LETTER_AXIS = {"x": 0, "y": 1, "z": 2, "a": 0, "b": 1, "g": 2}
        tokens = [str(t).strip() for t in (bush.point_raw or [])]
        if not tokens or not any(_has_letter(t) for t in tokens):
            return np.asarray(bush.point, dtype=float)
        atom0 = np.asarray(structure[indices[0]].frac_coords, dtype=float)
        letters: dict[str, float] = {}
        for tok in tokens:
            m = re.search(r"[a-zA-Z]", tok)
            if m and m.group() not in letters:
                letters[m.group()] = float(
                    atom0[_LETTER_AXIS.get(m.group(), 0)])
        point = self._eval_tokens(tokens, letters)
        return point if len(point) == 3 else np.asarray(bush.point, dtype=float)

    @staticmethod
    def _eval_token(tok: str, letters: dict[str, float]) -> float:
        """求值单个坐标 token（可含常数项与线性组合）。

        支持 ``1/2``、``x``、``-x``、``2x``、``-y+1``、``x-1/2``、
        ``-x+3/2``、``y+1/2``、``-2a+1`` 等形式。
        """
        s = tok.replace(" ", "")
        total = 0.0
        for m in re.finditer(r"([+-]?)(\d+(?:\.\d+)?(?:/\d+)?)?([a-zA-Z]*)", s):
            if not m.group(0):
                continue
            sign, coeff, alpha = m.group(1), m.group(2), m.group(3)
            if alpha:
                base = letters.get(alpha, 0.0) * (
                    float(Fraction(coeff)) if coeff else 1.0)
            else:
                base = float(Fraction(coeff)) if coeff else 0.0
            total += base if sign != "-" else -base
        return total

    @staticmethod
    def _eval_tokens(tokens: list[str], letters: dict[str, float]) -> np.ndarray:
        """按字母参数表求值坐标 token 列表。"""
        return np.asarray(
            [DistortionMapper._eval_token(t, letters) for t in tokens],
            dtype=float,
        )

    def _lattice_equiv(self, a: np.ndarray, b: np.ndarray,
                       atol: float = 1e-3) -> bool:
        """两个分数坐标是否相差一个母相格点（a - b ∈ L）。"""
        delta_cart = (a - b) @ self._conv_matrix
        n = delta_cart @ np.linalg.inv(self._prim_matrix)
        return bool(np.allclose(n, np.round(n), atol=atol))

    @staticmethod
    def _periodic_equiv(a: np.ndarray, b: np.ndarray,
                        atol: float = 1e-3) -> bool:
        """两个分数坐标是否周期等价（相差整数平移，mod 1 相等）。

        用于把 BUSH 代表点匹配到具体原子：代表点与原子分数坐标在 [0,1) 内
        重合即视为同一原子（如带心副本相位反号代表 `(-1/2,1/2,1/2)` 与
        体心原子 `(1/2,1/2,1/2)` mod 1 相等）。
        """
        delta = a - b
        return bool(np.allclose(delta - np.round(delta), 0.0, atol=atol))

    def _periodic_distance(self, a: np.ndarray, b: np.ndarray) -> float:
        """Return the shortest periodic distance in the direct-lattice metric.

        Fractional Euclidean distance has no crystallographic meaning unless
        the basis is orthonormal.  ``get_distance_and_image`` minimizes over
        lattice translations using the actual parent lattice, including skew
        and anisotropic cells.
        """
        if self._conv_matrix is None:
            raise RuntimeError("parent lattice is not initialized")
        return float(Lattice(self._conv_matrix).get_distance_and_image(a, b)[0])


def _has_letter(token: str) -> bool:
    """token 是否含字母（自由参数标记，如 x/y/z）。"""
    return bool(re.search(r"[a-zA-Z]", str(token)))
