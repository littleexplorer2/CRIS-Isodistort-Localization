"""超空间群与对称操作（IT-C 标准取位；nmod = d）。

三维空间群操作 ``x' = R x + v`` 嵌入 (3+nmod) 维：
    M = [[R, 0], [0, ε]],  t = (v, 0)
内部块 ε ∈ GL(nmod, Z) 由调制波矢在点群下的封闭性决定（IT-C 基本嵌入，
内部平移取 0，保证与惯用胞带心代表元的群运算闭合）。
``nmod=0`` 退化为普通三维空间群；``nmod=1`` 即 (3+1)。
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import permutations

import numpy as np
from pymatgen.symmetry.groups import SpaceGroup

from ..utils import DimensionMismatchError, NumericalSingularError, SymmetryIncompatibleError
from ..utils.opd_format import _centering_letter, _g_allowed
from ..utils.schoenflies import hm_symbol, schoenflies_symbol
from .constants import EPS, PHYSICAL_DIM, eps
from .validate import (
    dim_3pd,
    require_invertible,
    validate_matrix_shape,
    validate_nmod,
    validate_space_group_number,
    validate_vector_dim,
)


def wrap_translation(vec: np.ndarray, tol: float | None = None) -> np.ndarray:
    """分数平移归约到 [0, 1)。"""
    tol = EPS if tol is None else tol
    wrapped = np.mod(np.asarray(vec, dtype=float), 1.0)
    wrapped[np.abs(wrapped - 1.0) < tol] = 0.0
    wrapped[np.abs(wrapped) < tol] = 0.0
    return wrapped


@dataclass
class SuperspaceOperation:
    """一条 (3+nmod) 维超空间对称操作（旋转块 + 平移）。"""

    rotation: np.ndarray
    translation: np.ndarray
    nmod: int = 0

    def __post_init__(self) -> None:
        self.nmod = validate_nmod(self.nmod)
        dim = dim_3pd(self.nmod)
        self.rotation = validate_matrix_shape(
            self.rotation, (dim, dim), name="superspace rotation"
        )
        self.translation = wrap_translation(
            validate_vector_dim(self.translation, dim, name="superspace translation")
        )

    @property
    def dim(self) -> int:
        return dim_3pd(self.nmod)

    def compose(self, other: SuperspaceOperation) -> SuperspaceOperation:
        """self ∘ other，即先 other 再 self。"""
        if self.nmod != other.nmod:
            raise SymmetryIncompatibleError(
                f"Cannot compose superspace operations with nmod={self.nmod} and nmod={other.nmod}"
            )
        rot = self.rotation @ other.rotation
        trans = wrap_translation(self.translation + self.rotation @ other.translation)
        return SuperspaceOperation(rotation=rot, translation=trans, nmod=self.nmod)

    def inverse(self) -> SuperspaceOperation:
        inv = np.linalg.inv(require_invertible(self.rotation, name="superspace rotation"))
        trans = wrap_translation(-inv @ self.translation)
        return SuperspaceOperation(rotation=inv, translation=trans, nmod=self.nmod)

    def apply(self, vec: Sequence[float]) -> np.ndarray:
        v = validate_vector_dim(vec, self.dim, name="superspace vector")
        return wrap_translation(self.rotation @ v + self.translation)

    def is_identity(self, tol: float | None = None) -> bool:
        tol = eps() if tol is None else tol
        dim = self.dim
        if not np.allclose(self.rotation, np.eye(dim), atol=tol):
            return False
        delta = wrap_translation(self.translation, tol)
        return bool(np.allclose(delta, 0.0, atol=tol))

    def equivalent(self, other: SuperspaceOperation, tol: float | None = None) -> bool:
        tol = eps() if tol is None else tol
        if self.nmod != other.nmod:
            return False
        if not np.allclose(self.rotation, other.rotation, atol=tol):
            return False
        delta = wrap_translation(self.translation - other.translation, tol)
        return bool(np.allclose(delta, 0.0, atol=tol))

    def physical_rotation(self) -> np.ndarray:
        return np.asarray(self.rotation[:PHYSICAL_DIM, :PHYSICAL_DIM], dtype=float)

    def internal_rotation(self) -> np.ndarray:
        if self.nmod == 0:
            return np.zeros((0, 0))
        return np.asarray(
            self.rotation[PHYSICAL_DIM:, PHYSICAL_DIM:], dtype=float
        )

    def to_dict(self) -> dict:
        return {
            "nmod": self.nmod,
            "rotation": self.rotation.tolist(),
            "translation": self.translation.tolist(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> SuperspaceOperation:
        return cls(
            rotation=data["rotation"],
            translation=data["translation"],
            nmod=data.get("nmod", 0),
        )


def _signed_permutation_matrices(d: int) -> list[np.ndarray]:
    if d <= 0:
        return [np.zeros((0, 0))]
    out: list[np.ndarray] = []
    eye = np.eye(d, dtype=float)
    for perm in permutations(range(d)):
        p = eye[:, list(perm)]
        for signs in range(2 ** d):
            s = np.diag([1.0 if (signs >> i) & 1 == 0 else -1.0 for i in range(d)])
            out.append(s @ p)
    return out


def _internal_epsilon(
    rot3: np.ndarray,
    q_matrix: np.ndarray,
    centering: str,
    tol: float,
) -> np.ndarray:
    """求 ε ∈ GL(nmod,Z)，使 R^{-T} Q = Q ε + G（G 为允许的倒格矢）。"""
    d = q_matrix.shape[1]
    if d == 0:
        return np.zeros((0, 0))
    if np.allclose(q_matrix, 0.0, atol=tol):
        return np.eye(d)
    try:
        r_inv_t = np.linalg.inv(rot3).T
    except np.linalg.LinAlgError as exc:
        raise NumericalSingularError("3D rotation matrix is singular") from exc
    q_prime = r_inv_t @ q_matrix
    for eps_mat in _signed_permutation_matrices(d):
        residual = q_prime - q_matrix @ eps_mat
        ok = True
        for j in range(d):
            g = residual[:, j]
            rounded = np.rint(g)
            if not np.allclose(g, rounded, atol=tol):
                ok = False
                break
            h, k, ell = int(rounded[0]), int(rounded[1]), int(rounded[2])
            if not _g_allowed(h, k, ell, centering):
                ok = False
                break
        if ok and abs(round(float(np.linalg.det(eps_mat)))) == 1:
            return eps_mat
    raise SymmetryIncompatibleError(
        "Modulation wavevectors are not closed under the parent point group; "
        "cannot lift the 3D operation into (3+nmod) superspace (IT-C)."
    )


def _block_matrix(rot3: np.ndarray, epsilon: np.ndarray, nmod: int) -> np.ndarray:
    dim = dim_3pd(nmod)
    mat = np.eye(dim)
    mat[:PHYSICAL_DIM, :PHYSICAL_DIM] = rot3
    if nmod:
        mat[PHYSICAL_DIM:, PHYSICAL_DIM:] = epsilon
    return mat


def lift_operation(
    rot3: np.ndarray,
    trans3: np.ndarray,
    q_matrix: np.ndarray,
    nmod: int,
    centering: str,
    tol: float | None = None,
) -> SuperspaceOperation:
    """把三维操作 {R|v} 嵌入 (3+nmod)。q_matrix 为 3×nmod（列=调制波矢）。"""
    tol = eps() if tol is None else tol
    nmod = validate_nmod(nmod)
    rot3 = validate_matrix_shape(rot3, (3, 3), name="3D rotation")
    trans3 = validate_vector_dim(trans3, 3, name="3D translation")
    q_matrix = np.asarray(q_matrix, dtype=float).reshape(3, nmod) if nmod else np.zeros((3, 0))
    epsilon = _internal_epsilon(rot3, q_matrix, centering, tol)
    rot_s = _block_matrix(rot3, epsilon, nmod)
    trans = np.zeros(dim_3pd(nmod))
    trans[:PHYSICAL_DIM] = trans3
    return SuperspaceOperation(rotation=rot_s, translation=trans, nmod=nmod)


def lift_operation_3p1(
    rot3: np.ndarray,
    trans3: np.ndarray,
    q_vector: Sequence[float],
    centering: str,
    tol: float | None = None,
) -> SuperspaceOperation:
    """专用 (3+1) 提升：4×4 矩阵，ε=±1。与 ``lift_operation(..., nmod=1)`` 数值对齐。"""
    tol = eps() if tol is None else tol
    rot3 = validate_matrix_shape(rot3, (3, 3), name="3D rotation")
    trans3 = validate_vector_dim(trans3, 3, name="3D translation")
    q = validate_vector_dim(q_vector, 3, name="q_vector")
    try:
        q_prime = np.linalg.inv(rot3).T @ q
    except np.linalg.LinAlgError as exc:
        raise NumericalSingularError("3D rotation matrix is singular") from exc
    epsilon = None
    for eps_val in (1.0, -1.0):
        g = q_prime - eps_val * q
        rounded = np.rint(g)
        if not np.allclose(g, rounded, atol=tol):
            continue
        h, k, ell = int(rounded[0]), int(rounded[1]), int(rounded[2])
        if _g_allowed(h, k, ell, centering):
            epsilon = eps_val
            break
    if epsilon is None:
        raise SymmetryIncompatibleError(
            "q-vector is not closed under this 3D rotation; cannot build (3+1) superspace."
        )
    rot4 = np.eye(4)
    rot4[:3, :3] = rot3
    rot4[3, 3] = epsilon
    trans = np.zeros(4)
    trans[:3] = trans3
    return SuperspaceOperation(rotation=rot4, translation=trans, nmod=1)


def _parent_ops(space_group_number: int) -> list[tuple[np.ndarray, np.ndarray]]:
    sg = SpaceGroup.from_int_number(int(space_group_number))
    out: list[tuple[np.ndarray, np.ndarray]] = []
    for op in sg.symmetry_ops:
        rot = np.asarray(op.rotation_matrix, dtype=float)
        trans = np.asarray(op.translation_vector, dtype=float).reshape(3)
        out.append((rot, trans))
    return out


@dataclass
class SuperspaceGroup:
    """(3+nmod) 超空间群（由三维空间群 + nmod 个调制波矢构造）。"""

    space_group_number: int
    nmod: int
    operations: list[SuperspaceOperation] = field(default_factory=list)
    q_vectors: list[list[float]] = field(default_factory=list)
    space_group_symbol: str = ""
    schoenflies: str = ""
    centering: str = "P"
    setting: str = "standard (IT-C)"

    def __post_init__(self) -> None:
        self.nmod = validate_nmod(self.nmod)
        self.space_group_number = validate_space_group_number(self.space_group_number)
        if not self.space_group_symbol:
            self.space_group_symbol = hm_symbol(self.space_group_number)
        if not self.schoenflies:
            self.schoenflies = schoenflies_symbol(self.space_group_number)
        if not self.centering:
            self.centering = _centering_letter(self.space_group_number)

    @property
    def dim(self) -> int:
        return dim_3pd(self.nmod)

    @property
    def order(self) -> int:
        return len(self.operations)

    @classmethod
    def from_space_group(
        cls,
        space_group_number: int,
        nmod: int = 0,
        q_vectors: Sequence[Sequence[float]] | None = None,
        *,
        check_closed: bool = True,
    ) -> SuperspaceGroup:
        """从三维空间群构造 (3+nmod) 超空间群。

        ``q_vectors`` 为 None 或全零时用平凡嵌入（ε=I，τ=0），始终合法。
        给定与点群封闭的调制波矢时按 IT-C 基本嵌入提升。
        """
        nmod = validate_nmod(nmod)
        sg_num = validate_space_group_number(space_group_number)
        centering = _centering_letter(sg_num)
        q_matrix = np.zeros((3, nmod))
        q_list: list[list[float]] = []
        if nmod and q_vectors:
            if len(q_vectors) != nmod:
                raise DimensionMismatchError(
                    f"Need {nmod} modulation q-vectors for nmod={nmod}; got {len(q_vectors)}"
                )
            for i, q in enumerate(q_vectors):
                col = validate_vector_dim(q, 3, name=f"q_vectors[{i}]")
                q_matrix[:, i] = col
                q_list.append(col.tolist())
        else:
            q_list = [[0.0, 0.0, 0.0] for _ in range(nmod)]

        ops: list[SuperspaceOperation] = []
        for rot, trans in _parent_ops(sg_num):
            ops.append(lift_operation(rot, trans, q_matrix, nmod, centering))
        group = cls(
            space_group_number=sg_num,
            nmod=nmod,
            operations=ops,
            q_vectors=q_list,
            centering=centering,
        )
        if check_closed:
            group.assert_closed()
        return group

    @classmethod
    def from_space_group_3p1(
        cls,
        space_group_number: int,
        q_vector: Sequence[float] | None = None,
        *,
        check_closed: bool = True,
    ) -> SuperspaceGroup:
        """(3+1) 专用构造（与 ``from_space_group(..., nmod=1)`` 对照用）。"""
        q = [0.0, 0.0, 0.0] if q_vector is None else list(q_vector)
        sg_num = validate_space_group_number(space_group_number)
        centering = _centering_letter(sg_num)
        ops = [
            lift_operation_3p1(rot, trans, q, centering)
            for rot, trans in _parent_ops(sg_num)
        ]
        group = cls(
            space_group_number=sg_num,
            nmod=1,
            operations=ops,
            q_vectors=[list(validate_vector_dim(q, 3, name="q_vector"))],
            centering=centering,
        )
        if check_closed:
            group.assert_closed()
        return group

    def identity(self) -> SuperspaceOperation:
        dim = self.dim
        return SuperspaceOperation(
            rotation=np.eye(dim), translation=np.zeros(dim), nmod=self.nmod
        )

    def find_inverse(self, op: SuperspaceOperation) -> SuperspaceOperation:
        inv = op.inverse()
        for cand in self.operations:
            if cand.equivalent(inv):
                return cand
        raise SymmetryIncompatibleError("Inverse operation is not in the superspace group")

    def contains(self, op: SuperspaceOperation, tol: float | None = None) -> bool:
        return any(op.equivalent(other, tol) for other in self.operations)

    def assert_closed(self, tol: float | None = None) -> None:
        """乘法闭合 + 每元有逆。"""
        tol = eps() if tol is None else tol
        ident = self.identity()
        if not any(op.equivalent(ident, tol) for op in self.operations):
            raise SymmetryIncompatibleError("Superspace group is missing the identity")
        for op in self.operations:
            inv = op.inverse()
            if not any(inv.equivalent(other, tol) for other in self.operations):
                raise SymmetryIncompatibleError(
                    "Superspace group is not closed under inversion of operations"
                )
        for a in self.operations:
            for b in self.operations:
                prod = a.compose(b)
                if not any(prod.equivalent(other, tol) for other in self.operations):
                    raise SymmetryIncompatibleError(
                        "Superspace operations are not closed under composition"
                    )

    def check_custom_operations_closed(
        self,
        operations: Sequence[SuperspaceOperation],
        tol: float | None = None,
    ) -> None:
        """校验一组自定义超空间操作是否构成群（闭合 / 有逆）。"""
        tmp = SuperspaceGroup(
            space_group_number=self.space_group_number,
            nmod=self.nmod,
            operations=list(operations),
            q_vectors=self.q_vectors,
            centering=self.centering,
        )
        tmp.assert_closed(tol)

    def to_dict(self) -> dict:
        return {
            "kind": "superspace",
            "nmod": self.nmod,
            "space_group_number": self.space_group_number,
            "space_group_symbol": self.space_group_symbol,
            "schoenflies": self.schoenflies,
            "centering": self.centering,
            "setting": self.setting,
            "q_vectors": self.q_vectors,
            "operations": [op.to_dict() for op in self.operations],
        }

    @classmethod
    def from_dict(cls, data: dict) -> SuperspaceGroup:
        nmod = validate_nmod(data.get("nmod", data.get("d", 0)))
        ops = [SuperspaceOperation.from_dict(item) for item in data.get("operations", [])]
        return cls(
            space_group_number=int(data.get("space_group_number", 1)),
            nmod=nmod,
            operations=ops,
            q_vectors=[list(map(float, q)) for q in data.get("q_vectors", [])],
            space_group_symbol=str(data.get("space_group_symbol", "")),
            schoenflies=str(data.get("schoenflies", "")),
            centering=str(data.get("centering", "")),
            setting=str(data.get("setting", "standard (IT-C)")),
        )
