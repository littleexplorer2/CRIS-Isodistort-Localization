"""(3+nmod) 超空间波矢、点阵度量、倒格子、3D ↔ 超空间投影。"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from pymatgen.core import Lattice

from ..utils import DimensionMismatchError, NumericalSingularError, SymmetryIncompatibleError
from ..utils.opd_format import _canonical_k, _centering_letter, _g_allowed, _k_equivalent
from .constants import PHYSICAL_DIM, eps
from .validate import dim_3pd, require_invertible, validate_nmod, validate_vector_dim


def default_lattice_matrix(space_group_number: int) -> np.ndarray:
    """无母相 CIF 时用的常规晶胞（边长 1；正方/六方取对应晶系默认角）。"""
    n = int(space_group_number)
    if 168 <= n <= 194 or 143 <= n <= 167:
        return np.asarray(Lattice.hexagonal(1.0, 1.0).matrix, dtype=float)
    if 75 <= n <= 142:
        return np.asarray(Lattice.tetragonal(1.0, 1.0).matrix, dtype=float)
    if 195 <= n <= 230:
        return np.asarray(Lattice.cubic(1.0).matrix, dtype=float)
    if 16 <= n <= 74:
        return np.asarray(Lattice.orthorhombic(1.0, 1.0, 1.0).matrix, dtype=float)
    if 3 <= n <= 15:
        return np.asarray(Lattice.monoclinic(1.0, 1.0, 1.0, 90.0).matrix, dtype=float)
    return np.asarray(Lattice.triclinic(1.0, 1.0, 1.0, 90.0, 90.0, 90.0).matrix, dtype=float)


@dataclass
class SuperspaceLattice:
    """(3+nmod) 正格子度量、倒格子、子空间投影。内部维取单位度量（IT-C 分数坐标）。"""

    lattice_3d: np.ndarray
    nmod: int
    space_group_number: int = 1

    def __post_init__(self) -> None:
        self.nmod = validate_nmod(self.nmod)
        self.lattice_3d = np.asarray(self.lattice_3d, dtype=float)
        if self.lattice_3d.shape != (3, 3):
            raise DimensionMismatchError(
                f"lattice_3d must be 3x3; got {self.lattice_3d.shape}"
            )
        require_invertible(self.lattice_3d, name="lattice_3d")

    @property
    def dim(self) -> int:
        return dim_3pd(self.nmod)

    @property
    def metric(self) -> np.ndarray:
        g = np.eye(self.dim)
        g[:PHYSICAL_DIM, :PHYSICAL_DIM] = self.lattice_3d @ self.lattice_3d.T
        return g

    @property
    def reciprocal_metric(self) -> np.ndarray:
        try:
            return np.linalg.inv(self.metric)
        except np.linalg.LinAlgError as exc:
            raise NumericalSingularError("superspace metric is singular") from exc

    @property
    def reciprocal_lattice_3d(self) -> np.ndarray:
        return np.linalg.inv(self.lattice_3d).T

    def embed_3d(self, vec3: Sequence[float], internal: Sequence[float] | None = None) -> np.ndarray:
        v = validate_vector_dim(vec3, 3, name="3D vector")
        if self.nmod == 0:
            return v
        if internal is None:
            extra = np.zeros(self.nmod)
        else:
            extra = validate_vector_dim(internal, self.nmod, name="internal vector")
        return np.concatenate([v, extra])

    def project_3d(self, vec: Sequence[float]) -> np.ndarray:
        v = validate_vector_dim(vec, self.dim, name="superspace vector")
        return v[:PHYSICAL_DIM].copy()

    def to_dict(self) -> dict:
        return {
            "nmod": self.nmod,
            "space_group_number": self.space_group_number,
            "lattice_3d": self.lattice_3d.tolist(),
            "metric": self.metric.tolist(),
            "reciprocal_metric": self.reciprocal_metric.tolist(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> SuperspaceLattice:
        return cls(
            lattice_3d=data["lattice_3d"],
            nmod=data.get("nmod", 0),
            space_group_number=int(data.get("space_group_number", 1)),
        )


def _wrap_ks_internal(coords: np.ndarray, nmod: int, tol: float) -> np.ndarray:
    out = coords.copy()
    if nmod:
        extra = np.mod(out[PHYSICAL_DIM:], 1.0)
        extra[np.abs(extra - 1.0) < tol] = 0.0
        extra[np.abs(extra) < tol] = 0.0
        out[PHYSICAL_DIM:] = extra
    return out


@dataclass
class KsVector:
    """超空间波矢 k_s ∈ R^{3+nmod}（分数坐标，与 Method 2 k 点同一套分数/分数记号）。"""

    coords: np.ndarray
    nmod: int
    k_point_label: str = ""
    space_group_number: int = 1

    def __post_init__(self) -> None:
        self.nmod = validate_nmod(self.nmod)
        self.coords = validate_vector_dim(
            self.coords, dim_3pd(self.nmod), name="k_s"
        )

    @property
    def dim(self) -> int:
        return dim_3pd(self.nmod)

    @property
    def external(self) -> np.ndarray:
        return self.coords[:PHYSICAL_DIM].copy()

    @property
    def internal(self) -> np.ndarray:
        if self.nmod == 0:
            return np.zeros(0)
        return self.coords[PHYSICAL_DIM:].copy()

    def reduce(self, centering: str | None = None, tol: float | None = None) -> KsVector:
        """约化到惯用倒胞代表（外部分量沿用 Method 1 k-star 的 Bravais 规则）。"""
        tol = eps() if tol is None else tol
        letter = centering or _centering_letter(self.space_group_number)
        ext = _canonical_k(self.external, letter)
        coords = np.concatenate([ext, self.internal]) if self.nmod else ext
        coords = _wrap_ks_internal(coords, self.nmod, tol)
        return KsVector(
            coords=coords,
            nmod=self.nmod,
            k_point_label=self.k_point_label,
            space_group_number=self.space_group_number,
        )

    def equivalent(self, other: KsVector, centering: str | None = None, tol: float | None = None) -> bool:
        if self.nmod != other.nmod:
            return False
        tol = eps() if tol is None else tol
        letter = centering or _centering_letter(self.space_group_number)
        if not _k_equivalent(self.external, other.external, letter):
            return False
        if self.nmod == 0:
            return True
        delta = self.internal - other.internal
        rounded = np.rint(delta)
        return bool(np.allclose(delta, rounded, atol=tol))

    def project_to_3d(self) -> np.ndarray:
        return self.external

    @classmethod
    def from_3d(
        cls,
        k3: Sequence[float],
        nmod: int,
        *,
        satellite: Sequence[float] | None = None,
        k_point_label: str = "",
        space_group_number: int = 1,
    ) -> KsVector:
        """三维 k 嵌入超空间。nmod≥1 时默认第一内部坐标为 1（一阶卫星）。"""
        nmod = validate_nmod(nmod)
        ext = validate_vector_dim(k3, 3, name="k")
        if nmod == 0:
            extra = np.zeros(0)
        elif satellite is None:
            extra = np.zeros(nmod)
            extra[0] = 1.0
        else:
            extra = validate_vector_dim(satellite, nmod, name="k_s internal")
        return cls(
            coords=np.concatenate([ext, extra]) if nmod else ext,
            nmod=nmod,
            k_point_label=k_point_label,
            space_group_number=space_group_number,
        )

    def star(self, rotations_3d: Sequence[np.ndarray], centering: str | None = None) -> list[KsVector]:
        """三维点群作用下的 k 星（外部分量），内部坐标随 ε=±1 的平凡复制。"""
        letter = centering or _centering_letter(self.space_group_number)
        arms: list[KsVector] = [self.reduce(letter)]
        k0 = self.external
        for rot in rotations_3d:
            r = np.asarray(rot, dtype=float)
            try:
                kp = np.linalg.inv(r).T @ k0
            except np.linalg.LinAlgError:
                continue
            cand = KsVector(
                coords=np.concatenate([kp, self.internal]) if self.nmod else kp,
                nmod=self.nmod,
                k_point_label=self.k_point_label,
                space_group_number=self.space_group_number,
            ).reduce(letter)
            if not any(cand.equivalent(seen, letter) for seen in arms):
                arms.append(cand)
        return arms

    def to_dict(self) -> dict:
        return {
            "nmod": self.nmod,
            "coords": self.coords.tolist(),
            "k_point_label": self.k_point_label,
            "space_group_number": self.space_group_number,
        }

    @classmethod
    def from_dict(cls, data: dict) -> KsVector:
        return cls(
            coords=data["coords"],
            nmod=data.get("nmod", 0),
            k_point_label=str(data.get("k_point_label", "")),
            space_group_number=int(data.get("space_group_number", 1)),
        )


def assert_ks_compatible(
    ks: KsVector,
    q_vectors: Sequence[Sequence[float]],
    tol: float | None = None,
) -> None:
    """波矢外部分量须落在调制波矢 Q 张成的模格子上（加倒格矢）。"""
    tol = eps() if tol is None else tol
    if ks.nmod == 0:
        return
    q = np.asarray(q_vectors, dtype=float)
    if q.size == 0 or np.allclose(q, 0.0, atol=tol):
        return
    if q.shape != (ks.nmod, 3) and q.shape != (3, ks.nmod):
        if q.ndim == 2 and q.shape[1] == 3:
            qmat = q.T
        else:
            raise SymmetryIncompatibleError("q_vectors shape is incompatible with k_s")
    else:
        qmat = q if q.shape[0] == 3 else q.T
    # k_ext ≈ Q @ alpha + G
    try:
        alpha, *_ = np.linalg.lstsq(qmat, ks.external, rcond=None)
        recon = qmat @ alpha
    except np.linalg.LinAlgError as exc:
        raise NumericalSingularError("cannot project k_s onto modulation q-vectors") from exc
    residual = ks.external - recon
    rounded = np.rint(residual)
    if not np.allclose(residual, rounded, atol=tol):
        raise SymmetryIncompatibleError(
            "k_s external components are not compatible with the superspace modulation q-vectors"
        )
    letter = _centering_letter(ks.space_group_number)
    h, k, ell = int(rounded[0]), int(rounded[1]), int(rounded[2])
    if not _g_allowed(h, k, ell, letter):
        raise SymmetryIncompatibleError(
            "k_s differs from the modulation span by a reciprocal vector forbidden by Bravais centering"
        )
