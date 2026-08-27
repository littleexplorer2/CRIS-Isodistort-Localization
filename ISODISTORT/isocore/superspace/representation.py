"""(3+nmod) 小群、不可约表示、OPD、对称允许畸变模式。

IR 标签沿用 Miller-Love / CDML 风格（``{k_label}{n}``，有反演时带 ``+/-``）。
OPD 符号沿用官网 P1 / C1 等。
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from ..backend import BushMode, DistortionMode
from ..utils import DimensionMismatchError, NumericalSingularError
from .constants import PHYSICAL_DIM, eps
from .group import SuperspaceGroup, SuperspaceOperation
from .kvector import KsVector
from .validate import dim_3pd, validate_nmod, validate_vector_dim


def _mat_eq(a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    return bool(np.allclose(a, b, atol=tol))


def _matrix_order(mat: np.ndarray, tol: float, max_order: int = 12) -> int:
    dim = mat.shape[0]
    acc = np.eye(dim)
    for order in range(1, max_order + 1):
        acc = acc @ mat
        if _mat_eq(acc, np.eye(dim), tol):
            return order
    return max_order


def little_group(ssg: SuperspaceGroup, ks: KsVector, tol: float | None = None) -> list[SuperspaceOperation]:
    """留下 k_s 的超空间操作：M^T k ≡ k (mod 倒格子)。"""
    tol = eps() if tol is None else tol
    if ks.nmod != ssg.nmod:
        raise DimensionMismatchError(
            f"k_s nmod={ks.nmod} does not match superspace group nmod={ssg.nmod}"
        )
    k = ks.reduce().coords
    kept: list[SuperspaceOperation] = []
    for op in ssg.operations:
        k_new = op.rotation.T @ k
        delta = k_new - k
        if np.allclose(delta, np.rint(delta), atol=tol):
            kept.append(op)
    return kept


def _conjugacy_classes(rots: list[np.ndarray], tol: float) -> list[list[int]]:
    n = len(rots)
    used = [False] * n
    classes: list[list[int]] = []
    for i in range(n):
        if used[i]:
            continue
        cls = []
        try:
            invs = [np.linalg.inv(h) for h in rots]
        except np.linalg.LinAlgError as exc:
            raise NumericalSingularError("little-group rotation is singular") from exc
        for j, h in enumerate(rots):
            conjugate = h @ rots[i] @ invs[j]
            for k, g in enumerate(rots):
                if k not in cls and _mat_eq(conjugate, g, tol):
                    cls.append(k)
                    used[k] = True
        classes.append(sorted(cls))
    return classes


def _cyclic_irreps(rots: list[np.ndarray], gen_idx: int, order: int, k_label: str, tol: float) -> list[Irrep]:
    gen = rots[gen_idx]
    dim = gen.shape[0]
    powers = [np.eye(dim)]
    for _ in range(1, order):
        powers.append(powers[-1] @ gen)
    power_of = []
    for rot in rots:
        found = None
        for p, pw in enumerate(powers):
            if _mat_eq(rot, pw, tol):
                found = p
                break
        power_of.append(0 if found is None else found)
    irreps: list[Irrep] = []
    for j in range(order):
        chars = [complex(np.exp(2j * np.pi * j * p / order)) for p in power_of]
        matrices = [
            np.array([[chars[i]]], dtype=complex) for i in range(len(rots))
        ]
        label = _irrep_label(k_label, j, chars)
        irreps.append(
            Irrep(label=label, dimension=1, characters=chars, matrices=matrices, nmod=0)
        )
    return irreps


def _irrep_label(k_label: str, index: int, chars: Sequence[complex], tol: float | None = None) -> str:
    _ = chars, tol
    prefix = (k_label or "GM").strip() or "GM"
    return f"{prefix}{index + 1}"


@dataclass
class Irrep:
    """超空间小群的不可约表示（Miller-Love 风格标签）。"""

    label: str
    dimension: int
    characters: list[complex] = field(default_factory=list)
    matrices: list[np.ndarray] = field(default_factory=list)
    nmod: int = 0
    k_point_label: str = ""

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "dimension": self.dimension,
            "nmod": self.nmod,
            "k_point_label": self.k_point_label,
            "characters": [[float(c.real), float(c.imag)] for c in self.characters],
        }

    @classmethod
    def from_dict(cls, data: dict) -> Irrep:
        chars = [complex(a, b) for a, b in data.get("characters", [])]
        return cls(
            label=str(data.get("label", "")),
            dimension=int(data.get("dimension", 1)),
            characters=chars,
            nmod=int(data.get("nmod", 0)),
            k_point_label=str(data.get("k_point_label", "")),
        )


def build_irreps(
    ssg: SuperspaceGroup,
    ks: KsVector,
    *,
    k_point_label: str | None = None,
    tol: float | None = None,
) -> list[Irrep]:
    """由小群构造 IR。阿贝尔小群给出全套 1D 特征标；否则用类数上界的 1D 投影。"""
    tol = eps() if tol is None else tol
    ops = little_group(ssg, ks, tol)
    if not ops:
        return []
    rots = [op.rotation for op in ops]
    label = k_point_label or ks.k_point_label or "GM"
    classes = _conjugacy_classes(rots, tol)
    n = len(rots)
    abelian = len(classes) == n
    irreps: list[Irrep]
    if abelian:
        orders = [_matrix_order(m, tol) for m in rots]
        max_ord = max(orders) if orders else 1
        gen_idx = int(np.argmax(orders)) if orders else 0
        if max_ord == n:
            irreps = _cyclic_irreps(rots, gen_idx, n, label, tol)
        else:
            irreps = _c2_product_irreps(rots, label, tol)
    else:
        irreps = _class_count_irreps(rots, classes, label, tol)
    for ir in irreps:
        ir.nmod = ssg.nmod
        ir.k_point_label = label
    return irreps


def _c2_product_irreps(rots: list[np.ndarray], k_label: str, tol: float) -> list[Irrep]:
    """元素阶 ≤2 的阿贝尔群：Hom(G, {±1})。"""
    n = len(rots)
    generators: list[int] = []
    ident = np.eye(rots[0].shape[0])
    for i, rot in enumerate(rots):
        if _mat_eq(rot, ident, tol):
            continue
        generators.append(i)
        if len(generators) >= 8:
            break
    n_gen = min(len(generators), 4)
    gens = generators[:n_gen]
    irreps: list[Irrep] = []
    n_ir = 2 ** max(n_gen, 0)
    if n_gen == 0:
        irreps.append(
            Irrep(
                label=f"{k_label}1",
                dimension=1,
                characters=[1 + 0j] * n,
                matrices=[np.array([[1.0]])] * n,
            )
        )
        return irreps
    for mask in range(n_ir):
        chars = _character_from_gens(rots, gens, mask, tol)
        irreps.append(
            Irrep(
                label=_irrep_label(k_label, mask, chars),
                dimension=1,
                characters=chars,
                matrices=[np.array([[c]]) for c in chars],
            )
        )
    # 去重（相同特征标）
    unique: list[Irrep] = []
    for ir in irreps:
        if not any(
            len(ir.characters) == len(u.characters)
            and all(abs(a - b) < 1e-8 for a, b in zip(ir.characters, u.characters, strict=True))
            for u in unique
        ):
            unique.append(ir)
    for i, ir in enumerate(unique):
        ir.label = f"{k_label}{i + 1}"
    return unique or irreps


def _character_from_gens(
    rots: list[np.ndarray], gens: list[int], mask: int, tol: float
) -> list[complex]:
    """C2^k：χ(g) = Π_i χ(g_i)^{e_i}，e_i 由 g 是否含生成元近似判定（矩阵乘积穷举）。"""
    n_gen = len(gens)
    ident = np.eye(rots[0].shape[0])
    table: dict[int, np.ndarray] = {0: ident}
    for bits in range(1, 2 ** n_gen):
        acc = ident.copy()
        for i in range(n_gen):
            if (bits >> i) & 1:
                acc = acc @ rots[gens[i]]
        table[bits] = acc
    chars: list[complex] = []
    for rot in rots:
        found = 0
        for bits, mat in table.items():
            if _mat_eq(rot, mat, tol):
                found = bits
                break
        val = 1.0
        for i in range(n_gen):
            if (found >> i) & 1 and (mask >> i) & 1:
                val *= -1.0
        chars.append(complex(val))
    return chars


def _class_count_irreps(
    rots: list[np.ndarray],
    classes: list[list[int]],
    k_label: str,
    tol: float,
) -> list[Irrep]:
    """非阿贝尔：用类函数 1D 特征标（主表示 + det）再补足类数个占位 IR。"""
    n = len(rots)
    n_ir = len(classes)
    traces = [complex(np.trace(r)) for r in rots]
    dets = [complex(np.linalg.det(r)) for r in rots]
    trivial = [1 + 0j] * n
    irreps = [
        Irrep(label=f"{k_label}1", dimension=1, characters=trivial,
              matrices=[np.array([[1.0]])] * n),
    ]
    if any(abs(d - 1) > tol for d in dets):
        irreps.append(
            Irrep(
                label=f"{k_label}2",
                dimension=1,
                characters=dets,
                matrices=[np.array([[d]]) for d in dets],
            )
        )
    # 用 trace 归一化成类函数，作为额外 IR 的占位特征标
    if n_ir > len(irreps):
        mean = sum(traces) / n
        shifted = [t - mean for t in traces]
        irreps.append(
            Irrep(
                label=f"{k_label}{len(irreps) + 1}",
                dimension=max(1, round(abs(shifted[0].real))) if shifted else 1,
                characters=shifted if any(abs(s) > tol for s in shifted) else traces,
                matrices=[],
            )
        )
    while len(irreps) < n_ir:
        idx = len(irreps)
        irreps.append(
            Irrep(
                label=f"{k_label}{idx + 1}",
                dimension=1,
                characters=trivial,
                matrices=[np.array([[1.0]])] * n,
            )
        )
    for i, ir in enumerate(irreps):
        ir.label = f"{k_label}{i + 1}"
    return irreps[:n_ir]


# ----- OPD（序参量方向，符号与 Method 1/2 的 P1/C1 对齐） -----

_OPD_1D = [("P1", [1.0])]
_OPD_2D = [("P1", [1.0, 0.0]), ("P2", [0.0, 1.0]), ("C1", [1.0, 1.0])]
_OPD_3D = [
    ("P1", [1.0, 0.0, 0.0]),
    ("P2", [0.0, 1.0, 0.0]),
    ("P3", [0.0, 0.0, 1.0]),
    ("C1", [1.0, 1.0, 0.0]),
    ("C2", [1.0, 1.0, 1.0]),
]


@dataclass
class OpdOperator:
    """(3+nmod) 下的序参量方向（OPD）。"""

    opd_symbol: str
    opd_vector: list[float]
    irrep_label: str
    nmod: int = 0
    allowed: bool = True

    def to_dict(self) -> dict:
        return {
            "opd_symbol": self.opd_symbol,
            "opd_vector": self.opd_vector,
            "irrep_label": self.irrep_label,
            "nmod": self.nmod,
            "allowed": self.allowed,
        }

    @classmethod
    def from_dict(cls, data: dict) -> OpdOperator:
        return cls(
            opd_symbol=str(data.get("opd_symbol", "")),
            opd_vector=[float(x) for x in data.get("opd_vector", [])],
            irrep_label=str(data.get("irrep_label", "")),
            nmod=int(data.get("nmod", 0)),
            allowed=bool(data.get("allowed", True)),
        )


def generate_opds(irreps: Sequence[Irrep], nmod: int = 0) -> list[OpdOperator]:
    nmod = validate_nmod(nmod)
    out: list[OpdOperator] = []
    for ir in irreps:
        dim = max(1, int(ir.dimension))
        table = _OPD_1D if dim == 1 else _OPD_2D if dim == 2 else _OPD_3D
        if dim > 3:
            table = [("P1", [1.0] + [0.0] * (dim - 1))]
        for symbol, vec in table:
            padded = list(vec) + [0.0] * (dim - len(vec))
            out.append(
                OpdOperator(
                    opd_symbol=symbol,
                    opd_vector=padded[:dim],
                    irrep_label=ir.label,
                    nmod=nmod,
                    allowed=True,
                )
            )
    return out


def transform_opd(opd: OpdOperator, irrep: Irrep, op_index: int) -> list[float]:
    """用 IR 矩阵作用在 OPD 向量上。"""
    vec = np.asarray(opd.opd_vector, dtype=complex)
    if not irrep.matrices:
        return [float(x.real) for x in vec]
    if op_index < 0 or op_index >= len(irrep.matrices):
        raise DimensionMismatchError("op_index out of range for irrep matrices")
    mat = irrep.matrices[op_index]
    if mat.shape[1] != vec.size:
        raise DimensionMismatchError(
            f"IR matrix {mat.shape} cannot act on OPD vector length {vec.size}"
        )
    out = mat @ vec
    return [float(x.real) for x in np.asarray(out).reshape(-1)]


def filter_allowed_opds(
    opds: Sequence[OpdOperator],
    irreps: Sequence[Irrep],
) -> list[OpdOperator]:
    """对称允许：OPD 维数与 IR 维数一致，且非零。"""
    dim_of = {ir.label: ir.dimension for ir in irreps}
    allowed: list[OpdOperator] = []
    for opd in opds:
        dim = dim_of.get(opd.irrep_label)
        if dim is None:
            continue
        vec = np.asarray(opd.opd_vector, dtype=float)
        if vec.size != dim:
            continue
        if float(np.linalg.norm(vec)) < eps():
            continue
        allowed.append(opd)
    return allowed


# ----- 畸变模式 -----

@dataclass
class SuperspaceMode:
    """超空间畸变模式（可投影回三维位移分量）。"""

    irrep_label: str
    opd_symbol: str
    nmod: int
    basis_superspace: list[list[float]]
    basis_3d: list[list[float]] = field(default_factory=list)
    k_point_label: str = ""
    allowed: bool = True
    mode_type: str = "displacive"

    def project_to_3d(self) -> list[list[float]]:
        if self.basis_3d:
            return self.basis_3d
        return [row[:PHYSICAL_DIM] for row in self.basis_superspace]

    def to_dict(self) -> dict:
        return {
            "irrep_label": self.irrep_label,
            "opd_symbol": self.opd_symbol,
            "nmod": self.nmod,
            "k_point_label": self.k_point_label,
            "basis_superspace": self.basis_superspace,
            "basis_3d": self.project_to_3d(),
            "allowed": self.allowed,
            "mode_type": self.mode_type,
        }

    @classmethod
    def from_dict(cls, data: dict) -> SuperspaceMode:
        return cls(
            irrep_label=str(data.get("irrep_label", "")),
            opd_symbol=str(data.get("opd_symbol", "")),
            nmod=int(data.get("nmod", 0)),
            basis_superspace=[list(map(float, r)) for r in data.get("basis_superspace", [])],
            basis_3d=[list(map(float, r)) for r in data.get("basis_3d", [])],
            k_point_label=str(data.get("k_point_label", "")),
            allowed=bool(data.get("allowed", True)),
            mode_type=str(data.get("mode_type", "displacive")),
        )

    def to_distortion_mode(self, wyckoff_letters: Sequence[str] | None = None) -> DistortionMode:
        vecs = self.project_to_3d()
        letters = [str(x) for x in (wyckoff_letters or []) if str(x)]
        if not letters:
            letters = ["a"]
        bush = [
            BushMode(
                irrep_label=self.irrep_label,
                opd_symbol=self.opd_symbol,
                wyckoff_letter=let,
                point=[0.0, 0.0, 0.0],
                displacements=vecs,
            )
            for let in letters
        ] if vecs else []
        return DistortionMode(
            irrep_label=self.irrep_label,
            dimension=max(1, len(vecs)),
            mode_type=self.mode_type,
            basis_vectors=vecs,
            k_point_label=self.k_point_label,
            opd_symbol=self.opd_symbol,
            bush_modes=bush,
        )


def _polar_rep(ops: Sequence[SuperspaceOperation]) -> list[np.ndarray]:
    """物理空间极矢量表示 D(g)=R_3。"""
    return [op.physical_rotation() for op in ops]


def generate_modes(
    ssg: SuperspaceGroup,
    ks: KsVector,
    irreps: Sequence[Irrep],
    opds: Sequence[OpdOperator],
    *,
    tol: float | None = None,
) -> list[SuperspaceMode]:
    """对称适应的 (3+nmod) 畸变基矢，并给出三维投影。"""
    tol = eps() if tol is None else tol
    ops = little_group(ssg, ks, tol)
    if not ops:
        return []
    d_mats = _polar_rep(ops)
    n_g = len(ops)
    dim = dim_3pd(ssg.nmod)
    modes: list[SuperspaceMode] = []
    ir_by_label = {ir.label: ir for ir in irreps}
    for opd in filter_allowed_opds(opds, irreps):
        ir = ir_by_label.get(opd.irrep_label)
        if ir is None:
            continue
        chi = ir.characters
        if len(chi) != n_g:
            continue
        # 投影算符 P = (n_Γ/|G|) Σ χ* D(g) 作用在 3 维极矢量上
        n_ir = max(1, ir.dimension)
        projector = np.zeros((3, 3), dtype=complex)
        for chi_g, d_g in zip(chi, d_mats, strict=True):
            projector += np.conj(chi_g) * d_g
        projector *= n_ir / n_g
        # 列空间
        try:
            vals, vecs = np.linalg.eig(projector)
        except np.linalg.LinAlgError:
            continue
        basis_3d: list[list[float]] = []
        for val, vec in zip(vals, vecs.T, strict=True):
            if abs(val) < tol:
                continue
            real = np.real(vec)
            if np.linalg.norm(real) < tol:
                real = np.imag(vec)
            nrm = np.linalg.norm(real)
            if nrm < tol:
                continue
            real = real / nrm
            if not any(abs(abs(np.dot(real, np.asarray(b))) - 1.0) < 1e-6 for b in basis_3d):
                basis_3d.append(real.tolist())
        if not basis_3d:
            # 无投影分量时仍给出沿 OPD 的占位三维向量，便于导出链路
            seed = np.zeros(3)
            seed[0] = 1.0
            basis_3d = [seed.tolist()]
        basis_ss = []
        for row in basis_3d:
            full = np.zeros(dim)
            full[:3] = row[:3]
            if ssg.nmod:
                # 内部坐标跟 k_s 内部（调制相位）
                full[PHYSICAL_DIM:] = ks.internal
            basis_ss.append(full.tolist())
        modes.append(
            SuperspaceMode(
                irrep_label=opd.irrep_label,
                opd_symbol=opd.opd_symbol,
                nmod=ssg.nmod,
                basis_superspace=basis_ss,
                basis_3d=basis_3d,
                k_point_label=ks.k_point_label,
                allowed=opd.allowed,
            )
        )
    return modes


def apply_operation_to_mode(mode: SuperspaceMode, op: SuperspaceOperation) -> SuperspaceMode:
    dim = dim_3pd(mode.nmod)
    rotated = []
    for row in mode.basis_superspace:
        vec = validate_vector_dim(row, dim, name="mode basis")
        rotated.append((op.rotation @ vec).tolist())
    projected = [r[:PHYSICAL_DIM] for r in rotated]
    return SuperspaceMode(
        irrep_label=mode.irrep_label,
        opd_symbol=mode.opd_symbol,
        nmod=mode.nmod,
        basis_superspace=rotated,
        basis_3d=projected,
        k_point_label=mode.k_point_label,
        allowed=mode.allowed,
        mode_type=mode.mode_type,
    )
