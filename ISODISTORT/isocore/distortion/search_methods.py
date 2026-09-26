"""Search workflows aligned with ISODISTORT Method 1-4（基于真实 iso 枚举）。"""

from __future__ import annotations

import copy
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from heapq import merge

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.groups import SpaceGroup

from ..backend import DistortionMode, IsoWrapper, KPointInfo, SubgroupInfo
from ..data.kpoints_official import (
    KPOINT_OFFICIAL,
    official_kparams_to_iso,
)
from ..utils.config_loader import get_config
from ..utils.lattice import (
    as_fraction,
    centering_primitive_matrix,
    determinant,
    identity_matrix,
    inverse,
    is_sublattice,
    matrix_key,
    multiply,
    rational_matrix,
    same_lattice_in_point_group_orbit,
    transpose,
)
from ..utils.opd_format import _centering_letter, format_k_active
from .affine_embeddings import (
    AffineEmbedding,
    ParentAffineGroup,
    affine_equivalence,
    embedding_from_identity,
    enumerate_target_affine_embeddings,
    parent_affine_group,
    standardized_embedding_identity,
)
from .inverse_landau_adapter import stable_embedding_id
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


def _parse_basis_rows(
    rows: Sequence[Sequence[str | int | float | Fraction]],
) -> list[list[Fraction]]:
    """Parse a basis once and retain its exact rational entries."""
    return [list(row) for row in rational_matrix(rows)]


def _lattice_equivalent(a: Sequence[Sequence[float]],
                        b: Sequence[Sequence[float]]) -> bool:
    """两个 3x3 超胞基矢是否生成同一格点（GL(3,Z) 幺模等价）。

    Method 3 在用户指定**非恒等**目标子格时用此判定（对齐官网「指定格子」）。
    """
    try:
        return same_lattice_in_point_group_orbit(
            rational_matrix(a), rational_matrix(b), None
        )
    except ValueError:
        return False


def _is_identity_basis(basis: Sequence[Sequence[float]], atol: float = 1e-6) -> bool:
    arr = np.asarray(basis, dtype=float)
    return arr.shape == (3, 3) and bool(np.allclose(arr, np.eye(3), atol=atol))


def _integer_basis_matrix(basis: Sequence[Sequence[float]]) -> np.ndarray | None:
    """精确整数的 3x3 超胞矩阵；否则 None。"""
    try:
        exact = rational_matrix(basis)
    except ValueError:
        return None
    if any(value.denominator != 1 for row in exact for value in row):
        return None
    return np.asarray([[int(value) for value in row] for row in exact], dtype=int)


def _validate_centering(value: str | None) -> str:
    """校验 Method 3 的 direct_sublattice_centering（官网 d/P/A/B/C/I/F/R radio）。

    官网帮助：未选带心时按所选空间群默认带心，点群则按 primitive。
    表单选 **P** 表示 primitive / no centering（结果摘要常写 ``no centering``），
    与 Default(``d``) 的格子语义不同。A/B/C/I/F/R 同样只规定如何把所填
    conventional basis 还原成唯一 primitive translation lattice；它们不是
    另一套子群数据库开关。
    """
    if value is None or str(value).strip() == "":
        return "d"
    text = str(value).strip()
    key = text.lower()
    if key in ("d", "default"):
        return "d"
    letter = text.upper()
    if letter in {"P", "A", "B", "C", "I", "F", "R"}:
        return letter
    raise ValueError(f"Unknown Method 3 centering: {value!r}")


def _selected_primitive_lattice(
    parent_sg: int,
    query: Method3Query,
    basis: Sequence[Sequence[float]],
    *,
    parent_primitive=None,
    standard_basis=None,
):
    """Interpret the website's basis+centering controls as one exact lattice.

    Website/ISO basis rows are expressed in the standard conventional parent
    frame.  Exact affine operations, however, live in the fractional axes of
    the uploaded CIF.  ``standard_basis`` is the exact row-basis bridge from
    the former to the latter.  ``parent_primitive`` is reconstructed from the
    uploaded structure's identity-rotation translations.  The legacy HM-only
    fallback is retained solely for low-level callers that provide no parent
    structure.
    """
    centering = _validate_centering(query.direct_sublattice_centering)
    if centering == "d":
        centering = (
            _centering_letter(query.space_group_type)
            if query.space_group_type is not None
            else "P"
        )
    selected_standard = multiply(
        centering_primitive_matrix(centering), rational_matrix(basis)
    )
    if standard_basis is None:
        standard_basis = identity_matrix()
    selected = multiply(selected_standard, rational_matrix(standard_basis))
    parent = (
        rational_matrix(parent_primitive)
        if parent_primitive is not None
        else centering_primitive_matrix(_centering_letter(parent_sg))
    )
    if not is_sublattice(selected, parent):
        raise ValueError(
            "Basis vectors must be lattice vectors of the parent. "
            "The selected centering and basis do not define a parent-lattice sublattice."
        )
    return selected, parent


def _subgroup_primitive_lattice(subgroup: SubgroupInfo, standard_basis=None):
    # A subgroup row can use a non-standard setting whose centering differs
    # from pymatgen's standard symbol for the international number.  The HM
    # symbol printed by iso describes that actual row and is therefore the
    # primary source; the number is only a fallback for malformed/empty data.
    symbol = str(subgroup.space_group_symbol or "").strip()
    match = re.match(r"([PABCIFR])", symbol, flags=re.IGNORECASE)
    centering = (
        match.group(1).upper()
        if match is not None
        else _centering_letter(subgroup.space_group_number)
    )
    basis = subgroup.basis_vectors or identity_matrix()
    primitive_standard = multiply(
        centering_primitive_matrix(centering),
        rational_matrix(basis),
    )
    frame = (
        identity_matrix()
        if standard_basis is None
        else rational_matrix(standard_basis)
    )
    return multiply(
        primitive_standard,
        frame,
    )


def _fraction_str(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _param_value_candidates(M: np.ndarray) -> list[Fraction]:
    """Return the finite translation quotient's complete scalar grid.

    This compatibility helper is used by the superspace code.  It derives the
    exponent (largest Smith invariant) from exact determinants/minors instead
    of scanning every divisor up to ``det(M)`` or computing a float
    determinant.  Method 3 line searches use
    :func:`_commensurate_line_parameters` below, which is more selective and
    solves the actual affine line rather than guessing ``1/d`` values.
    """
    exact = rational_matrix(M)
    if any(value.denominator != 1 for row in exact for value in row):
        raise ValueError("Supercell matrix must contain exact integers")
    det = abs(int(determinant(exact)))
    if det <= 1:
        return []

    entries = [[int(value) for value in row] for row in exact]
    minors: list[int] = []
    for row_a in range(3):
        for row_b in range(row_a + 1, 3):
            for col_a in range(3):
                for col_b in range(col_a + 1, 3):
                    minors.append(abs(
                        entries[row_a][col_a] * entries[row_b][col_b]
                        - entries[row_a][col_b] * entries[row_b][col_a]
                    ))
    gcd_minors = 0
    for value in minors:
        gcd_minors = math.gcd(gcd_minors, value)
    exponent = det // gcd_minors if gcd_minors else det
    return sorted({Fraction(n, exponent) for n in range(1, exponent)})


def _k_compatible_with_supercell(
    kvec: Sequence[str | int | float | Fraction],
    translations: Sequence[Sequence[str | int | float | Fraction]],
) -> bool:
    """k 是否为超胞倒格点（M k ∈ ℤ³）。

    ``basis_vectors`` 与 :func:`build_supercell` 一致，每一行是一条用
    母相基矢表示的子格基矢。因此子格平移 ``m_i`` 与母相
    倒格坐标 ``k`` 的相位是第 ``i`` 行与 ``k`` 的点积，即
    ``M @ k``。转置只会在对角超胞上偶然给出相同结果，轴置换或
    剪切子格会错误丢弃本应折叠到 Gamma 的波矢。
    """
    if len(kvec) != 3:
        return False
    try:
        matrix = rational_matrix(translations)
        vector = tuple(as_fraction(value) for value in kvec)
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    return all(
        sum(row[col] * vector[col] for col in range(3)).denominator == 1
        for row in matrix
    )


_RATIONAL_LITERAL = r"(?:\d+(?:/\d+)?|(?:\d+\.\d*|\.\d+))"


def _parse_affine_k_component(
    component: str,
    param_name: str,
) -> tuple[Fraction, Fraction] | None:
    """Parse a strict one-variable affine coordinate ``slope*p + offset``.

    Accepted terms are exact rational constants and rational multiples of the
    selected parameter, with optional ``*`` and arbitrary term order.  Thus
    website forms such as ``-a+1``, ``1-a`` and ``a+1/2`` are handled without
    ``eval`` or floating-point arithmetic.  Other variables, products and
    nonlinear expressions are rejected.
    """
    name = str(param_name).strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
        return None
    text = re.sub(r"\s+", "", str(component))
    if not text:
        return None
    signed = text if text[0] in "+-" else f"+{text}"
    terms = re.findall(r"[+-][^+-]+", signed)
    if not terms or "".join(terms) != signed:
        return None

    parameter_term = re.compile(
        rf"(?:(?P<coefficient>{_RATIONAL_LITERAL})\*?)?{re.escape(name)}"
    )
    constant_term = re.compile(_RATIONAL_LITERAL)
    slope = Fraction(0)
    offset = Fraction(0)
    for term in terms:
        sign = -1 if term[0] == "-" else 1
        body = term[1:]
        match = parameter_term.fullmatch(body)
        if match is not None:
            raw_coefficient = match.group("coefficient")
            coefficient = Fraction(raw_coefficient) if raw_coefficient else Fraction(1)
            slope += sign * coefficient
            continue
        if constant_term.fullmatch(body) is not None:
            offset += sign * Fraction(body)
            continue
        return None
    return slope, offset


def _eval_simple_k_component(
    comp: str,
    param_name: str,
    param_val: Fraction,
) -> Fraction | None:
    affine = _parse_affine_k_component(comp, param_name)
    if affine is None:
        return None
    slope, offset = affine
    return slope * Fraction(param_val) + offset


def _kvec_from_template(
    coords: Sequence[str], param_name: str, param_val: Fraction
) -> tuple[Fraction, Fraction, Fraction] | None:
    comps: list[Fraction] = []
    for comp in coords:
        val = _eval_simple_k_component(comp, param_name, param_val)
        if val is None:
            return None
        comps.append(val)
    if len(comps) != 3:
        return None
    return comps[0], comps[1], comps[2]


def _matrix_vector(matrix, vector: Sequence[Fraction]) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(
        sum(matrix[row][col] * vector[col] for col in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _parent_standard_basis(parent: ParentAffineGroup | None):
    """Exact standard-conventional direct basis in uploaded parent axes."""

    if parent is None or parent.standard_basis is None:
        return identity_matrix()
    return parent.standard_basis


def _standard_k_to_parent(
    vector: Sequence[Fraction], standard_basis
) -> tuple[Fraction, Fraction, Fraction]:
    """Transform a reciprocal vector from standard to uploaded coordinates.

    Direct row coordinates obey ``v_input = v_std S``.  Phase invariance then
    gives ``k_input = S^-1 k_std`` for reciprocal-coordinate columns.
    """

    return _matrix_vector(inverse(rational_matrix(standard_basis)), vector)


def _subgroup_identity_in_parent_axes(
    subgroup: SubgroupInfo,
    parent: ParentAffineGroup,
) -> dict[str, object]:
    """Convert one ISO standard-frame subgroup identity to uploaded axes."""

    standard_basis = _parent_standard_basis(parent)
    basis = multiply(
        rational_matrix(subgroup.basis_vectors or identity_matrix()),
        standard_basis,
    )
    standard_origin_values = tuple(
        as_fraction(value) for value in (subgroup.origin or [0, 0, 0])
    )
    if len(standard_origin_values) != 3:
        raise ValueError("embedding origin must contain three coordinates")
    standard_origin = standard_origin_values  # type: ignore[assignment]
    parent_origin = tuple(parent.standard_origin)
    uploaded_origin = tuple(
        left + right
        for left, right in zip(
            _matrix_vector(transpose(standard_basis), standard_origin),
            parent_origin,
            strict=True,
        )
    )
    return {
        "space_group_number": int(subgroup.space_group_number),
        "symbol": str(subgroup.space_group_symbol),
        "basis": basis,
        "origin": uploaded_origin,
    }


def _embedding_identity_in_standard_axes(
    identity: dict[str, object],
    parent: ParentAffineGroup,
) -> dict[str, object]:
    """Convert an exact uploaded-axis identity back to the ISO/UI frame."""

    standard_basis = _parent_standard_basis(parent)
    uploaded_basis = rational_matrix(identity["basis"])  # type: ignore[arg-type]
    standard_conventional_basis = multiply(
        uploaded_basis, inverse(standard_basis)
    )
    uploaded_origin_values = tuple(
        as_fraction(value) for value in identity["origin"]  # type: ignore[union-attr]
    )
    if len(uploaded_origin_values) != 3:
        raise ValueError("embedding origin must contain three coordinates")
    relative_origin = tuple(
        uploaded_origin_values[index] - parent.standard_origin[index]
        for index in range(3)
    )
    standard_origin = _matrix_vector(
        inverse(transpose(standard_basis)),
        relative_origin,
    )
    return {
        **identity,
        "basis": standard_conventional_basis,
        "origin": standard_origin,
    }


def _affine_k_star_arms(
    slopes: Sequence[Fraction],
    offsets: Sequence[Fraction],
    parent_rotations: Sequence[Sequence[Sequence[str | int | float | Fraction]]] | None,
) -> list[tuple[tuple[Fraction, Fraction, Fraction], tuple[Fraction, Fraction, Fraction]]]:
    """Return exact reciprocal-star affine arms ``k(t)=s*t+q``."""
    base_slope = (slopes[0], slopes[1], slopes[2])
    base_offset = (offsets[0], offsets[1], offsets[2])
    arms = [(base_slope, base_offset)]
    seen = {(base_slope, base_offset)}
    for raw_rotation in parent_rotations or []:
        rotation = rational_matrix(raw_rotation)
        reciprocal_rotation = transpose(inverse(rotation))
        arm = (
            _matrix_vector(reciprocal_rotation, base_slope),
            _matrix_vector(reciprocal_rotation, base_offset),
        )
        if arm not in seen:
            seen.add(arm)
            arms.append(arm)
    return arms


def _reciprocal_star_vectors(
    vector: Sequence[Fraction],
    parent_rotations: Sequence[
        Sequence[Sequence[str | int | float | Fraction]]
    ] | None,
) -> set[tuple[Fraction, Fraction, Fraction]]:
    """Return exact reciprocal-space star vectors before lattice reduction."""
    zero = (Fraction(0), Fraction(0), Fraction(0))
    return {
        arm_vector
        for _slope, arm_vector in _affine_k_star_arms(
            zero,
            vector,
            parent_rotations,
        )
    }


def _special_k_vector(parent_sg: int, kp: KPointInfo):
    """Return one exact special-k representative, preferring website metadata."""
    override = KPOINT_OFFICIAL.get(int(parent_sg), {}).get(kp.label.strip())
    if override is not None:
        _kovalev, coords, params = override
        if params:
            return None
    else:
        if not kp.is_special:
            return None
        coords = list(kp.coordinates or [])
    if len(coords) != 3:
        return None
    try:
        values = tuple(as_fraction(value) for value in coords)
    except (ValueError, ZeroDivisionError):
        return None
    return values


def _special_k_star_vectors(
    parent_sg: int,
    kpoints: Sequence[KPointInfo],
    parent_rotations: Sequence[
        Sequence[Sequence[str | int | float | Fraction]]
    ] | None,
    standard_basis=None,
) -> set[tuple[Fraction, Fraction, Fraction]]:
    """All known special points and their exact reciprocal-space stars."""
    # Gamma is universal and is part of the special-k pass even if a minimal
    # test/backend k-point table omits its explicit row.
    keys: set[tuple[Fraction, Fraction, Fraction]] = {
        (Fraction(0), Fraction(0), Fraction(0))
    }
    for kp in kpoints:
        vector = _special_k_vector(parent_sg, kp)
        if vector is not None:
            if standard_basis is not None:
                vector = _standard_k_to_parent(vector, standard_basis)
            keys.update(_reciprocal_star_vectors(vector, parent_rotations))
    return keys


def _reciprocal_vectors_equivalent(
    left: Sequence[Fraction],
    right: Sequence[Fraction],
    parent_primitive,
) -> bool:
    """Whether ``left-right`` belongs to the parent's reciprocal lattice.

    The k coordinates use the parent conventional reciprocal basis.  For a
    centered parent, reducing each coordinate modulo one is wrong: e.g.
    I-centered ``(1,1,1)`` is not Gamma.  A reciprocal shift ``delta`` is
    allowed exactly when every parent primitive direct translation has an
    integer phase, ``P @ delta in Z^3``.
    """
    delta = tuple(left[index] - right[index] for index in range(3))
    phases = _matrix_vector(parent_primitive, delta)
    return all(value.denominator == 1 for value in phases)


def _is_special_k_equivalent(
    vector: Sequence[Fraction],
    special_star_vectors: set[tuple[Fraction, Fraction, Fraction]],
    parent_primitive,
    parent_rotations: Sequence[
        Sequence[Sequence[str | int | float | Fraction]]
    ] | None,
) -> bool:
    """Check one parameter value against all special-k stars exactly."""
    return any(
        _reciprocal_vectors_equivalent(arm, special, parent_primitive)
        for arm in _reciprocal_star_vectors(vector, parent_rotations)
        for special in special_star_vectors
    )


def _same_reciprocal_k_star(
    left: Sequence[Fraction],
    right: Sequence[Fraction],
    parent_primitive,
    parent_rotations: Sequence[
        Sequence[Sequence[str | int | float | Fraction]]
    ]
    | None,
) -> bool:
    """Whether two representatives describe one physical parent k-star.

    ``left`` is expanded by the parent point group and compared with ``right``
    modulo the reciprocal lattice dual to the *actual uploaded* parent
    translation lattice.  This must be applied before querying the parameter-k
    backend: otherwise symmetry-related values such as ``1/3`` and ``2/3``
    can acquire different printed origins and survive a later embedding
    deduplication even though they belong to one physical star.
    """

    return any(
        _reciprocal_vectors_equivalent(arm, right, parent_primitive)
        for arm in _reciprocal_star_vectors(left, parent_rotations)
    )


def _ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _iter_line_congruence_solutions(
    translations,
    slopes: Sequence[Fraction],
    offsets: Sequence[Fraction],
):
    """Yield ordered solutions without allocating a superlattice-sized list."""
    matrix = rational_matrix(translations)
    slope_terms = _matrix_vector(matrix, slopes)
    offset_terms = _matrix_vector(matrix, offsets)
    equations = list(zip(slope_terms, offset_terms, strict=True))

    for slope, offset in equations:
        if slope == 0 and offset.denominator != 1:
            return
    finite = [(slope, offset) for slope, offset in equations if slope != 0]
    if not finite:
        # A nonsingular translation matrix and a genuine parameter direction
        # cannot reach this branch.  Returning no discrete values is safer
        # than pretending that an unconstrained continuum is enumerable.
        return

    def root_count(equation: tuple[Fraction, Fraction]) -> int:
        slope, offset = equation
        if slope > 0:
            return max(0, _ceil_fraction(offset + slope) - _ceil_fraction(offset))
        return max(0, offset.numerator // offset.denominator - (
            (offset + slope).numerator // (offset + slope).denominator
        ))

    seed_slope, seed_offset = min(finite, key=root_count)
    if seed_slope > 0:
        integer_values = range(
            _ceil_fraction(seed_offset),
            _ceil_fraction(seed_offset + seed_slope),
        )
    else:
        lower = (seed_offset + seed_slope).numerator // (seed_offset + seed_slope).denominator
        upper = seed_offset.numerator // seed_offset.denominator
        # A negative slope reverses the integer-to-parameter ordering.
        integer_values = range(upper, lower, -1)

    for integer in integer_values:
        value = (Fraction(integer) - seed_offset) / seed_slope
        if 0 <= value < 1 and all(
            (slope * value + offset).denominator == 1
            for slope, offset in equations
        ):
            yield value


def _line_congruence_solutions(
    translations,
    slopes: Sequence[Fraction],
    offsets: Sequence[Fraction],
) -> list[Fraction]:
    """Solve ``T @ (slopes*t + offsets) in Z^3`` for every ``0 <= t < 1``."""
    return list(_iter_line_congruence_solutions(translations, slopes, offsets))


def _iter_commensurate_line_parameters(
    translations,
    coords: Sequence[str],
    param_name: str,
    parent_rotations: Sequence[
        Sequence[Sequence[str | int | float | Fraction]]
    ] | None = None,
    standard_basis=None,
):
    """Yield the sorted union of exact parameter values over all star arms."""
    parsed = [_parse_affine_k_component(component, param_name) for component in coords]
    if len(parsed) != 3 or any(item is None for item in parsed):
        return
    affine = [item for item in parsed if item is not None]
    slopes = tuple(item[0] for item in affine)
    offsets = tuple(item[1] for item in affine)
    if standard_basis is not None:
        slopes = _standard_k_to_parent(slopes, standard_basis)
        offsets = _standard_k_to_parent(offsets, standard_basis)
    arm_solutions = (
        _iter_line_congruence_solutions(translations, arm_slopes, arm_offsets)
        for arm_slopes, arm_offsets in _affine_k_star_arms(
            slopes,
            offsets,
            parent_rotations,
        )
    )
    previous: Fraction | None = None
    for value in merge(*arm_solutions):
        if value != previous:
            yield value
            previous = value


def _commensurate_line_parameters(
    translations,
    coords: Sequence[str],
    param_name: str,
    parent_rotations: Sequence[Sequence[Sequence[str | int | float | Fraction]]] | None = None,
    standard_basis=None,
) -> list[Fraction]:
    """Solve all commensurate values on every parent-point-group star arm.

    ``coords`` comes from the ISO standard conventional reciprocal frame.
    When the uploaded CIF uses primitive or nonstandard axes, transform both
    the affine slope and offset into those axes before applying the actual
    parent rotations and child translations.
    """
    return list(
        _iter_commensurate_line_parameters(
            translations,
            coords,
            param_name,
            parent_rotations,
            standard_basis,
        )
    )


def _line_kpoint_template(parent_sg: int, kp: KPointInfo) -> tuple[list[str], str] | None:
    """单参数线 k 点的坐标模板与参数名（优先官网表）。"""
    override = KPOINT_OFFICIAL.get(int(parent_sg), {}).get(kp.label.strip())
    if override is not None:
        _kov, coords, params = override
        if len(params) == 1:
            return [str(c) for c in coords], str(params[0])
        return None
    coords = [str(c) for c in (kp.coordinates or [])]
    # Older parser rows stored whole coordinate expressions (``"2a"`` or
    # ``"-a+1"``) in ``parameters``.  Derive identifiers from the coordinate
    # grammar itself so arbitrary parent groups get the same affine handling
    # as the curated website display table.
    params = sorted({
        name
        for component in coords
        for name in re.findall(r"[A-Za-z][A-Za-z0-9_]*", component)
    })
    if len(params) != 1:
        return None
    return coords, params[0]


def _space_group_to_point_group(space_group_number: int) -> str:
    try:
        return SpaceGroup.from_int_number(space_group_number).point_group
    except (ValueError, KeyError, TypeError):
        return ""


def _same_method3_route(
    left: SubgroupInfo,
    right: SubgroupInfo,
    parent_primitive,
    parent_rotations,
    standard_basis=None,
) -> bool:
    """Deduplicate one physical route, including equivalent k-star arms.

    Parameter values are coordinates on a named line, not physical identities:
    two values can describe the same star through different representatives.
    We therefore compare the exact physical k vectors modulo the parent
    reciprocal lattice and all supplied parent point-group rotations.  IR and
    OPD labels must still agree, so genuinely different second-stage routes
    into one affine embedding remain available.
    """

    labels_left = (
        left.irrep_label,
        left.opd_symbol,
    )
    labels_right = (
        right.irrep_label,
        right.opd_symbol,
    )
    if labels_left != labels_right:
        return False
    if (
        left.k_point_label,
        tuple(left.k_parameters or []),
    ) == (
        right.k_point_label,
        tuple(right.k_parameters or []),
    ):
        return True
    try:
        left_vector = tuple(as_fraction(value) for value in (left.k_coordinates or []))
        right_vector = tuple(as_fraction(value) for value in (right.k_coordinates or []))
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    if len(left_vector) != 3 or len(right_vector) != 3:
        return False
    if standard_basis is not None:
        left_vector = _standard_k_to_parent(left_vector, standard_basis)
        right_vector = _standard_k_to_parent(right_vector, standard_basis)
    return _same_reciprocal_k_star(
        left_vector,
        right_vector,
        parent_primitive,
        parent_rotations,
    )


def _basis_is_sublattice_of(basis: Sequence[Sequence[float]],
                            sublattice: Sequence[Sequence[float]],
                            parent_rotations: Sequence[np.ndarray] | None = None
                            ) -> bool:
    """
    判断子群超胞基矢 B 的格点是否是被选子格 S 的子格。

    用于判断 direct sublattice 包含关系：B 的每一行必须是 S 的整系数线性组合，
    即 N = B @ inv(S) 的元素全部为整数（S 可为对角阵或任意 3x3 矩阵）。

    母相点群任意旋转也可生成等价子格；若提供
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


def _basis_is_same_lattice(basis: Sequence[Sequence[float]],
                           selected: Sequence[Sequence[float]],
                           parent_rotations: Sequence[np.ndarray] | None = None
                           ) -> bool:
    """Match one official lattice selector class, including parent rotations."""
    b = np.asarray(basis, dtype=float)
    s = np.asarray(selected, dtype=float)
    if b.shape != (3, 3) or s.shape != (3, 3):
        return False
    if not np.isclose(abs(np.linalg.det(b)), abs(np.linalg.det(s)), atol=1e-6):
        return False
    return _basis_is_sublattice_of(b, s, parent_rotations)


@dataclass
class Method1Query:
    """Method 1: search over all special k points."""

    distortion_types: str | Sequence[str] | None = None
    crystal_system: str | Sequence[str] | None = None  # 单个或列表（多选=OR）
    subgroup_space_group: int | None = None
    lattice: Sequence[Sequence[float]] | None = None  # conventional 格点类（3x3 矩阵）
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
    number_of_independent_modulations: int = 0  # 0 = 3D lock-in; n = (3+n)D harmonics


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
    # 母相点群旋转；指定格子按完整点群轨道匹配，而不是只比输入代表基矢
    parent_rotations: Sequence[Sequence[Sequence[float]]] | None = None
    # The actual uploaded structure is required to reconstruct nonsymmorphic
    # parent operations and compare candidate embeddings exactly.  Direct
    # IsoSearchEngine callers may omit it, in which case the legacy textual
    # identity remains a conservative fallback (no candidates are merged).
    parent_structure: Structure | None = None
    # spglib/pymatgen parent detection: absolute Cartesian distance in Å.
    # ``None`` resolves from settings.yaml rather than inventing a fractional
    # or lattice-matrix tolerance at this call site.
    symmetry_tolerance: float | None = None
    symmetry_angle_tolerance_degrees: float | None = None
    # Stage-A data are exact rational operations converted to float only at
    # the spglib boundary, so they use a separate, stricter Å error budget.
    affine_exact_tolerance: float | None = None
    # Dimensionless residual for recovering spglib transforms/origins as
    # rationals; exact Seitz reconstruction is still required afterwards.
    fractional_coordinate_tolerance: float | None = None
    # Diagnostic only until a fixed-subspace/inverse-Landau calculation proves
    # activity for the selected distortion types and attaches a second-stage
    # route.  User-facing web/terminal calls leave this disabled.
    include_affine_only_diagnostics: bool = False
    # 参数 k 点子群库缺失时是否在线生成（与 Method 2 GenDB 同一开关）
    generate_if_missing: bool = False


@dataclass
class Method3ResultItem:
    """One affine subgroup embedding plus the single-IR routes known locally.

    The website's first Method 3 table identifies the embedding, not one
    particular IR/OPD.  ``subgroup`` remains as a compatibility representative
    while ``routes`` prevents distinct known directions from being discarded.
    """

    subgroup: SubgroupInfo
    point_group: str
    basis: list[list[float]]
    routes: list[SubgroupInfo] | None = None
    affine_embedding: AffineEmbedding | None = None
    route_resolution: str = "known_single_ir"
    embedding_id: str | None = None


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
        # Method 3 repeatedly applies different SG/lattice filters to the same
        # complete special-k table.  Cache that expensive WSL enumeration per
        # parent/type set; returned SubgroupInfo objects are only re-indexed for
        # the current result and their scientific identity is unchanged.
        self._method3_special_cache: dict[tuple[int, tuple[str, ...]], list[SubgroupInfo]] = {}

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
        - lattice：超胞格是否与所选格子相同（含母相点群旋转轨道）
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
                if not _basis_is_same_lattice(
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
                        wyckoff_letters: Sequence[str] | None = None,
                        *,
                        structure: Structure | None = None,
                        wyckoff_sites: Sequence[dict] | None = None,
                        smodes: object | None = None,
                        kpoints: Sequence | None = None,
                        symmetry_info: dict | None = None,
                        ) -> Method2Result:
        """
        官网 Method 2：在已枚举子群中按序号选择目标子群，计算其畸变模式。

        特殊 k：iso DISPLAY BUSH 根模 + rootless Wyckoff 轨道的子群固定
        子空间。参数 k：smodes + 子群对称性（(3+d)/公度锁定）。

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

        requested_nmod = int(query.number_of_independent_modulations or 0)
        if requested_nmod < 0:
            raise ValueError(
                "number_of_independent_modulations must be >= 0 "
                "(0 = 3D lock-in complete modes; n = (3+n)D harmonics of n q-vectors)."
            )
        if requested_nmod > 3:
            requested_nmod = 3
        # nmod is meaningful only for parameterized k vectors.  Preserve the
        # caller-visible request in metadata for API compatibility, while the
        # special-k computation itself always uses the ordinary 3D lock-in.
        effective_nmod = requested_nmod if target.k_parameters else 0

        # wyckoff_letters 为 None 表示调用方未提供（误用）；空列表表示
        # 作用域内无 Wyckoff 位置（如全部类型选 none）-> 直接返回空模式
        if wyckoff_letters is None:
            raise ValueError(
                "Method 2 计算模式需要母相结构的 Wyckoff 位置信息，"
                "请先加载结构（load_structure）"
            )
        modes = []
        extra_meta: dict = {}
        if wyckoff_letters:
            if target.k_parameters and structure is not None and smodes is not None:
                from .superspace import compute_parametric_modes  # noqa: PLC0415

                info = dict(symmetry_info or {})
                info.setdefault("space_group_number", parent_sg)
                if wyckoff_sites is not None:
                    info["wyckoff_sites"] = list(wyckoff_sites)
                complete_result = compute_parametric_modes(
                    structure,
                    info,
                    target,
                    wyckoff_letters,
                    smodes,  # type: ignore[arg-type]
                    kpoints=kpoints,
                    nmod=effective_nmod,
                )
                modes = complete_result.modes
                extra_meta["supercell_displacements"] = complete_result.supercell_displacements
                extra_meta["mode_labels"] = complete_result.labels
                extra_meta["parametric_note"] = complete_result.note
            elif structure is not None and smodes is not None:
                from .superspace import (  # noqa: PLC0415
                    compute_special_modes_with_rootless_supplement,
                )

                info = dict(symmetry_info or {})
                info.setdefault("space_group_number", parent_sg)
                if wyckoff_sites is not None:
                    info["wyckoff_sites"] = list(wyckoff_sites)
                complete_result = compute_special_modes_with_rootless_supplement(
                    structure,
                    info,
                    target,
                    wyckoff_letters,
                    self._iso,
                    smodes,  # type: ignore[arg-type]
                    kpoints=kpoints,
                )
                modes = complete_result.modes
                extra_meta["supercell_displacements"] = (
                    complete_result.supercell_displacements
                )
                extra_meta["mode_labels"] = complete_result.labels
                extra_meta["special_mode_note"] = complete_result.note
            else:
                # Backend-only callers have no structure on which to build the
                # child displacement representation.  They still receive the
                # exact BUSH root modes; the high-level API always supplies the
                # structure and therefore also completes rootless orbits.
                modes = self._iso.calc_distortion_modes(
                    parent_sg, target, wyckoff_letters=wyckoff_letters
                )

        metadata = {
            # k 点 / IR / OPD 来自所选子群自身（查询时已确定，无需重复传参）
            "k_point_label": target.k_point_label,
            "irrep_label": target.irrep_label,
            "opd_symbol": target.opd_symbol,
            "k_parameters": list(target.k_parameters),
            "number_of_independent_modulations": requested_nmod,
            **extra_meta,
        }

        return Method2Result(subgroup=target, modes=modes, metadata=metadata)

    # ----------------------------------------------------------------
    # Method 3：指定点群/空间群 + 超胞
    # ----------------------------------------------------------------

    def method_3_search(self, parent_sg: int, query: Method3Query) -> list[Method3ResultItem]:
        """
        官网 Method 3 的本地单-IR 子集实现（特殊 k + 公度参数 k）：

        - 若同时提供 point_group 与 space_group_type，空间群选择优先；
        - Default 按目标空间群的 Bravais centering 解释；P/A/B/C/I/F/R
          分别按显式 centering 解释，全部转换为精确 primitive lattice；
        - 所有基矢（包括恒等矩阵）都只匹配同一 lattice class，并遍历母相
          点群轨道；不再错误接纳任意更大超胞；
        - 已知单-IR 路径按 ``(SG,basis,origin)`` 去重；不能把物理畴直接
          当成首屏 embedding 候选，因为官网还会按母群仿射共轭分类。
        - 非母相原胞子格时，在每条单参数仿射 k 线及其完整 k-star 上精确解
          ``T k(p) in Z^3``；枚举 ``0 <= p < 1`` 的全部公度解，而非只猜
          ``1/d``，再查询各 IR 子群（可 ``generate_if_missing``）。

        已知限制：这仍不是完整的 inverse-Landau/COPL 搜索；仅由单-IR 路径
        出发，可能遗漏只有 coupled IR 才产生的嵌入。reciprocal 模式不支持；
        多参数平面/一般 k（GP）未自动推断。
        参数 k 的位移模式与 Method 2 相同，在选中子群后走 smodes/(3+d)。
        Method 2 只携带一个主波矢，因此 nmod=1、2、3 都只保留该 q 的谐波
        （含 Γ）；不会再引入第二个独立调制波矢。nmod=0 才保留全部折叠 k。
        """
        distortion_types = normalize_distortion_types(query.distortion_types)
        if not query.space_group_type and not query.point_group:
            raise ValueError("Method 3 requires a space-group type or point group")
        if str(query.lattice_type or "direct").lower() != "direct":
            raise ValueError("Method 3 reciprocal-space superlattices are not implemented")

        if query.point_group and query.space_group_type:
            point_group_filter = None
        else:
            point_group_filter = query.point_group

        basis = _parse_basis_rows(
            query.supercell_basis
            or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        )
        parent_model = None
        if query.parent_structure is not None:
            parent_model = parent_affine_group(
                query.parent_structure,
                query.symmetry_tolerance,
                query.symmetry_angle_tolerance_degrees,
            )
            if parent_model.space_group_number != int(parent_sg):
                raise ValueError(
                    "Method 3 parent structure symmetry does not match parent_sg: "
                    f"detected {parent_model.space_group_number}, requested {parent_sg}"
                )
        standard_basis = _parent_standard_basis(parent_model)
        selected_primitive, parent_primitive = _selected_primitive_lattice(
            parent_sg,
            query,
            basis,
            parent_primitive=(parent_model.lattice if parent_model is not None else None),
            standard_basis=standard_basis,
        )

        # 1) 特殊 k 点候选（Method 1 同源）
        special_key = (parent_sg, tuple(sorted(distortion_types)))
        if special_key not in self._method3_special_cache:
            self._method3_special_cache[special_key] = list(
                self._iso.enumerate_all_special_subgroups(parent_sg, distortion_types)
            )
        # Result tables own their row objects because each query reassigns its
        # display indices.  Never let a later filter mutate an earlier table.
        subgroups = copy.deepcopy(self._method3_special_cache[special_key])

        # 2) 非母相原胞子格 → 公度参数 k 回退（任意 k 的有限单参数子集）
        relative = multiply(selected_primitive, inverse(parent_primitive))
        if abs(determinant(relative)) > 1:
            parametric = self._method3_parametric_subgroups(
                parent_sg,
                selected_primitive,
                parent_primitive,
                distortion_types,
                parent_rotations=query.parent_rotations,
                standard_basis=standard_basis,
                generate_if_missing=bool(query.generate_if_missing),
            )
            subgroups.extend(parametric)

        # 3) 先按目标类型收窄。不要把 SHOW DOMAIN 的物理畴直接展开成
        # Method 3 首屏候选；实测 SG12 会由官网 3 个 embedding 错增为 8 个畴。
        typed: list[SubgroupInfo] = []
        for sg in subgroups:
            point_group = _space_group_to_point_group(sg.space_group_number)
            if query.space_group_type and sg.space_group_number != query.space_group_type:
                continue
            if point_group_filter and point_group != point_group_filter:
                continue
            typed.append(sg)

        expanded: list[tuple[SubgroupInfo, AffineEmbedding | None, bool]] = [
            (subgroup, None, False) for subgroup in typed
        ]
        matched: list[SubgroupInfo] = []
        embedding_routes: list[list[SubgroupInfo]] = []
        embedding_models: list[AffineEmbedding | None] = []
        if parent_model is not None:
            if query.space_group_type and query.include_affine_only_diagnostics:
                parent_cartesian = query.parent_structure.lattice.matrix
                affine_candidates = enumerate_target_affine_embeddings(
                    parent_model,
                    selected_primitive,
                    int(query.space_group_type),
                    parent_cartesian,
                    symprec=query.affine_exact_tolerance,
                )
                translation_index = int(abs(determinant(relative)))
                for embedding in affine_candidates:
                    identity = standardized_embedding_identity(
                        embedding,
                        parent_model,
                        parent_cartesian,
                        int(query.space_group_type),
                        symprec=query.affine_exact_tolerance,
                        fractional_tolerance=(
                            query.fractional_coordinate_tolerance
                        ),
                    )
                    if identity is None:
                        raise RuntimeError(
                            "Method 3 found an exact affine embedding but could not "
                            "recover an exactly verified conventional basis/origin; "
                            "the diagnostic result was not truncated"
                        )
                    standard_identity = _embedding_identity_in_standard_axes(
                        identity, parent_model
                    )
                    basis_rows = [list(row) for row in standard_identity["basis"]]
                    origin = list(standard_identity["origin"])
                    point_order = len(embedding.operations)
                    subgroup_index = (
                        len(parent_model.operations) * translation_index
                    ) // point_order
                    generated = SubgroupInfo(
                        index=-1,
                        space_group_number=int(query.space_group_type),
                        space_group_symbol=str(identity["symbol"]),
                        subgroup_index=subgroup_index,
                        size=translation_index,
                        basis_vectors=basis_rows,
                        origin=origin,
                        parent_sg=parent_sg,
                    )
                    expanded.append((generated, embedding, True))

        for sg, supplied_model, pure_affine in expanded:
            try:
                candidate_primitive = _subgroup_primitive_lattice(
                    sg, standard_basis
                )
            except ValueError:
                continue
            if not same_lattice_in_point_group_orbit(
                candidate_primitive,
                selected_primitive,
                query.parent_rotations,
            ):
                continue
            model = supplied_model
            if model is None and parent_model is not None:
                try:
                    model = embedding_from_identity(
                        _subgroup_identity_in_parent_axes(sg, parent_model),
                        parent_model,
                    )
                except (ArithmeticError, KeyError, TypeError, ValueError):
                    # Third-party rows can use an alternate Hall setting which
                    # the preferred-Hall reconstruction cannot yet map.  Such
                    # a row must remain visible; never discard it as a failed
                    # duplicate proof.
                    model = None

            match_index: int | None = None
            if model is not None and parent_model is not None:
                for index, existing in enumerate(embedding_models):
                    relationship = (
                        affine_equivalence(
                            model,
                            existing,
                            parent_model,
                            # Parent-conjugate orientation/domain embeddings
                            # are not equality and may be separate official
                            # rows.  Merge only a proven identical operator
                            # subgroup.
                            allow_parent_conjugacy=False,
                        )
                        if existing is not None
                        else None
                    )
                    if relationship is not None:
                        match_index = index
                        break
            if match_index is None and model is None:
                origin_key = tuple(
                    (value.numerator, value.denominator)
                    for value in (
                        as_fraction(item) % 1 for item in (sg.origin or [0, 0, 0])
                    )
                )
                legacy_key = (
                    int(sg.space_group_number),
                    matrix_key(rational_matrix(sg.basis_vectors or identity_matrix())),
                    origin_key,
                )
                for index, existing in enumerate(matched):
                    existing_origin = tuple(
                        (value.numerator, value.denominator)
                        for value in (
                            as_fraction(item) % 1
                            for item in (existing.origin or [0, 0, 0])
                        )
                    )
                    existing_key = (
                        int(existing.space_group_number),
                        matrix_key(
                            rational_matrix(existing.basis_vectors or identity_matrix())
                        ),
                        existing_origin,
                    )
                    if existing_key == legacy_key:
                        match_index = index
                        break
            if match_index is None:
                match_index = len(matched)
                matched.append(sg)
                embedding_routes.append([])
                embedding_models.append(model)
            elif pure_affine is False and not embedding_routes[match_index]:
                # Prefer a real second-stage route over the route-less
                # standard-setting representative produced by stage A.
                matched[match_index] = sg
            routes = embedding_routes[match_index]
            if not pure_affine and not any(
                _same_method3_route(
                    item,
                    sg,
                    parent_primitive,
                    query.parent_rotations,
                    standard_basis,
                )
                for item in routes
            ):
                routes.append(sg)

        result: list[Method3ResultItem] = []
        for i, sg in enumerate(matched):
            sg.index = i
            route_resolution = (
                "known_single_ir" if embedding_routes[i]
                else "affine_only_unresolved_coupled_route"
            )
            embedding_id = None
            if (
                route_resolution == "affine_only_unresolved_coupled_route"
                and parent_model is not None
                and embedding_models[i] is not None
            ):
                embedding_id = stable_embedding_id(
                    parent_model, embedding_models[i]
                )
            result.append(Method3ResultItem(
                subgroup=sg,
                point_group=_space_group_to_point_group(sg.space_group_number),
                basis=[list(row) for row in (sg.basis_vectors or [])],
                routes=embedding_routes[i],
                affine_embedding=embedding_models[i],
                route_resolution=route_resolution,
                embedding_id=embedding_id,
            ))
        return result

    def _method3_parametric_subgroups(
        self,
        parent_sg: int,
        selected_primitive,
        parent_primitive,
        distortion_types,
        parent_rotations=None,
        standard_basis=None,
        generate_if_missing: bool = False,
    ) -> list[SubgroupInfo]:
        """由超胞基矢推断公度线 k 点并枚举子群（Method 3 参数 k 回退）。"""
        relative = multiply(selected_primitive, inverse(parent_primitive))
        if any(value.denominator != 1 for row in relative for value in row):
            raise ValueError("Selected primitive lattice is not a parent-lattice sublattice")
        if abs(determinant(relative)) <= 1:
            return []

        config = get_config()
        parametric_limit = config.method3_max_parametric_values
        query_limit = config.method3_max_backend_queries

        try:
            kpoints = self._iso.list_k_points(parent_sg)
        except Exception as exc:
            raise RuntimeError(
                f"Method 3 failed to list k points for parent SG {parent_sg}"
            ) from exc
        special_k_vectors = _special_k_star_vectors(
            parent_sg,
            kpoints,
            parent_rotations,
            standard_basis,
        )

        searches: list[
            tuple[
                KPointInfo,
                list[str],
                list[str],
                tuple[Fraction, Fraction, Fraction],
            ]
        ] = []
        seen: set[tuple[str, tuple[str, ...]]] = set()
        seen_parametric_stars: list[tuple[Fraction, Fraction, Fraction]] = []
        for kp in kpoints:
            if kp.is_special:
                continue
            tmpl = _line_kpoint_template(parent_sg, kp)
            if tmpl is None:
                continue
            coords, pname = tmpl
            param_cands = _iter_commensurate_line_parameters(
                selected_primitive,
                coords,
                pname,
                parent_rotations,
                standard_basis,
            )
            for pval in param_cands:
                standard_kvec = _kvec_from_template(coords, pname, pval)
                if standard_kvec is None:
                    continue
                parent_kvec = (
                    _standard_k_to_parent(standard_kvec, standard_basis)
                    if standard_basis is not None
                    else standard_kvec
                )
                # Parameter-line endpoints and boundary values often carry a
                # special-point alias (possibly only on another k-star arm).
                # They were already enumerated in the special-k pass, so do
                # not ask the parameter backend for the same physical point.
                if _is_special_k_equivalent(
                    parent_kvec,
                    special_k_vectors,
                    parent_primitive,
                    parent_rotations,
                ):
                    continue
                # A named line can contain several parameter values that are
                # related by the parent point group modulo a reciprocal-lattice
                # vector (for example Y a=1/3 and 2/3 in P4/mmm).  They are one
                # k-star and must be queried once, using the first exact value
                # in the deterministic ascending enumeration.  Deduplicating
                # only after subgroup generation is too late because the
                # backend may print different origin representatives.
                if any(
                    _same_reciprocal_k_star(
                        parent_kvec,
                        existing,
                        parent_primitive,
                        parent_rotations,
                    )
                    for existing in seen_parametric_stars
                ):
                    continue
                official = [_fraction_str(pval)]
                key = (kp.label, tuple(official))
                if key in seen:
                    continue
                seen.add(key)
                seen_parametric_stars.append(parent_kvec)
                iso_vals = official_kparams_to_iso(parent_sg, kp.label, official, kp)
                searches.append((kp, official, iso_vals, standard_kvec))
                parametric_count = len(searches)
                if parametric_limit and parametric_count > parametric_limit:
                    raise RuntimeError(
                        "Method 3 exact commensurate search requires at least "
                        f"{parametric_count} parametric values, exceeding "
                        "runtime.method3_max_parametric_values="
                        f"{parametric_limit}. Candidate generation stopped at "
                        "the first proven excess; no parametric values were "
                        "truncated and no list_irreps query was started."
                    )
                if query_limit and parametric_count > query_limit:
                    raise RuntimeError(
                        "Method 3 requires at least "
                        f"{parametric_count} backend queries "
                        f"({parametric_count} list_irreps calls), exceeding "
                        "runtime.method3_max_backend_queries="
                        f"{query_limit}. Candidate generation stopped at the "
                        "first proven excess; no list_irreps query was started."
                    )

        out: list[SubgroupInfo] = []
        parametric_count = len(searches)

        # Two phases are intentional: collect and activity-filter every IR
        # first, then know the exact list_irreps + list_subgroups cost before
        # allowing a single subgroup query.  Silent truncation would make a
        # budget failure indistinguishable from a scientifically empty search.
        resolved_searches: list[
            tuple[
                KPointInfo,
                list[str],
                list[str],
                tuple[Fraction, Fraction, Fraction],
                list,
            ]
        ] = []
        for kp, official_params, iso_params, official_kvec in searches:
            try:
                irreps = self._iso.list_irreps(
                    parent_sg, kp.label, k_parameters=iso_params
                )
            except Exception as exc:
                raise RuntimeError(
                    "Method 3 failed to list irreps for "
                    f"{kp.label} with k parameters {official_params}"
                ) from exc
            irreps = [
                ir for ir in irreps
                if self._iso._include_irrep(ir, distortion_types)
            ]
            resolved_searches.append(
                (kp, official_params, iso_params, official_kvec, irreps)
            )

        estimated_queries = parametric_count + sum(
            len(irreps) for *_, irreps in resolved_searches
        )
        if query_limit and estimated_queries > query_limit:
            raise RuntimeError(
                "Method 3 requires "
                f"{estimated_queries} backend queries "
                f"({parametric_count} list_irreps + "
                f"{estimated_queries - parametric_count} list_subgroups), "
                "exceeding runtime.method3_max_backend_queries="
                f"{query_limit}. Increase runtime.method3_max_backend_queries "
                "or set it to 0 for unlimited; no list_subgroups query was "
                "started and no IRs were truncated."
            )

        for kp, official_params, iso_params, official_kvec, irreps in resolved_searches:
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
                except Exception as exc:
                    raise RuntimeError(
                        "Method 3 failed to list subgroups for "
                        f"{kp.label}/{ir.label} with k parameters {official_params}"
                    ) from exc
                for sg in batch:
                    sg.k_parameters = list(official_params)
                    sg.k_coordinates = [_fraction_str(value) for value in official_kvec]
                    # ``iso_params`` can use a backend-specific parameter
                    # scale (LD g=1/6 is represented as 1/12 by this bundled
                    # iso build).  Never retain the backend-form k-active text
                    # after restoring the website coordinate convention.
                    sg.k_active_raw = format_k_active(
                        sg.opd_dir_raw or "",
                        sg.k_coordinates or ["0", "0", "0"],
                        parent_sg,
                        None,
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
