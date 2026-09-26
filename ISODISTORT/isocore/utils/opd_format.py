"""Format Method 1 order-parameter lines to match the official ISODISTORT OPD page.

Visible radio text on the website looks like::

    GM1+     P1   (a) 139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=(0,0,0), s=1, i=1, k-active= (0,0,0)
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from fractions import Fraction
from numbers import Integral, Real

import numpy as np
from pymatgen.symmetry.groups import SpaceGroup

from .lattice import inverse, rational_matrix, transpose

# Official radio line: irrep left-padded to 9, OPD symbol to 5.
_IRREP_WIDTH = 9
_OPD_WIDTH = 5


def format_number(value: float | int | str | Fraction) -> str:
    """Render a coordinate without inventing a nearby rational value.

    Text and :class:`Fraction` inputs carry exact intent and are preserved (or
    reduced) exactly.  Binary floats are rendered as a fraction only when
    their *exact binary value* is a reasonably small dyadic rational; other
    floats use Python's round-trippable decimal spelling.  In particular,
    ``float(1/97)`` can no longer be rounded to the false value ``0``.
    """
    if isinstance(value, str):
        text = value.strip()
        if text:
            return text
        return "0"
    if isinstance(value, Fraction):
        frac = value
    elif isinstance(value, Integral):
        return str(int(value))
    elif isinstance(value, Real):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Coordinate must be finite, got {value!r}")
        if number == 0:
            return "0"
        exact = Fraction.from_float(number)
        if exact.denominator > 65_536:
            return repr(number)
        frac = exact
    else:
        raise TypeError(f"Unsupported coordinate value: {value!r}")
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def format_tuple(
    values: Sequence[float | int | str | Fraction],
    raw: str | None = None,
) -> str:
    """``(x,y,z)`` using the iso raw token when present."""
    if raw:
        text = raw.strip()
        if text.startswith("(") and text.endswith(")"):
            return text
        if text:
            return f"({text})"
    inner = ",".join(format_number(v) for v in values)
    return f"({inner})"


def format_basis(
    vectors: Sequence[Sequence[float | int | str | Fraction]],
    raw: str | None = None,
) -> str:
    """``(1,0,0),(0,1,0),(0,0,1)`` (no extra spaces), matching the website."""
    if raw:
        return raw.strip()
    return ",".join(format_tuple(row) for row in vectors)


def _opd_arm_nonzero(arm: str) -> bool:
    """True if an OPD semicolon-group is not identically zero."""
    for part in arm.split(","):
        token = part.strip()
        if not token:
            continue
        if token in {"0", "+0", "-0", "0.0", "+0.0", "-0.0"}:
            continue
        try:
            if Fraction(token) == 0:
                continue
        except (ValueError, ZeroDivisionError):
            pass
        return True
    return False


def _centering_letter(space_group_number: int) -> str:
    symbol = SpaceGroup.from_int_number(int(space_group_number)).symbol
    return (symbol[:1] or "P").upper()


def _integer_delta(k1: np.ndarray, k2: np.ndarray) -> tuple[int, int, int] | None:
    delta = np.asarray(k1, dtype=float) - np.asarray(k2, dtype=float)
    rounded = np.rint(delta)
    if not np.allclose(delta, rounded, atol=1e-6):
        return None
    return int(rounded[0]), int(rounded[1]), int(rounded[2])


def _g_allowed(h: int, k: int, ell: int, centering: str) -> bool:
    letter = (centering or "P")[:1].upper()
    if letter == "P":
        return True
    if letter == "I":
        return (h + k + ell) % 2 == 0
    if letter == "F":
        return (h % 2) == (k % 2) == (ell % 2)
    if letter == "A":
        return (k + ell) % 2 == 0
    if letter == "B":
        return (h + ell) % 2 == 0
    if letter == "C":
        return (h + k) % 2 == 0
    return True


def _k_score(k: np.ndarray) -> tuple:
    """Lower is better: stay in [0,1] inclusive, avoid negatives, then smaller |k|."""
    out_of_unit = int(sum(1 for x in k if x < -1e-8 or x > 1.0 + 1e-8))
    negative = int(sum(1 for x in k if x < -1e-8))
    abs_sum = float(np.sum(np.abs(k)))
    return (out_of_unit, negative, abs_sum)


def _canonical_k(k: np.ndarray, centering: str, max_g: int = 2) -> np.ndarray:
    """Pick a unit-cell representative of k modulo the Bravais reciprocal lattice.

    Used for rotated star arms only. The listed special-k seed (CDML / iso
    DISPLAY KPOINT, e.g. M = (1,1,1) for I4/mmm) is kept as given: reducing it
    by an allowed G would map M onto (0,0,1), which the website never prints.
    """
    best = np.asarray(k, dtype=float)
    best_score = _k_score(best)
    for h in range(-max_g, max_g + 1):
        for ki in range(-max_g, max_g + 1):
            for ell in range(-max_g, max_g + 1):
                if not _g_allowed(h, ki, ell, centering):
                    continue
                cand = np.asarray(k, dtype=float) + np.array([h, ki, ell], dtype=float)
                score = _k_score(cand)
                if score < best_score:
                    best = cand
                    best_score = score
    return best


def _zero_one_flips(k: np.ndarray) -> list[np.ndarray]:
    """0↔1 partners on axes that sit on the unit-cell boundary."""
    out: list[np.ndarray] = []
    for i in range(3):
        if abs(float(k[i])) < 1e-8:
            kp = np.asarray(k, dtype=float).copy()
            kp[i] = 1.0
            out.append(kp)
        elif abs(float(k[i]) - 1.0) < 1e-8:
            kp = np.asarray(k, dtype=float).copy()
            kp[i] = 0.0
            out.append(kp)
    return out


def _k_equivalent(k1: np.ndarray, k2: np.ndarray, centering: str) -> bool:
    """Two conventional k vectors are the same modulo the Bravais reciprocal lattice."""
    delta = _integer_delta(k1, k2)
    if delta is None:
        return False
    h, k, ell = delta
    letter = (centering or "P")[:1].upper()
    if letter == "P":
        return True
    if letter == "I":
        return (h + k + ell) % 2 == 0
    if letter == "F":
        return (h % 2) == (k % 2) == (ell % 2)
    if letter == "A":
        return (k + ell) % 2 == 0
    if letter == "B":
        return (h + ell) % 2 == 0
    if letter == "C":
        return (h + k) % 2 == 0
    if letter == "R":
        return True
    return True


def _k_to_array(coords: Sequence[str | float]) -> np.ndarray:
    out = []
    for c in coords:
        if isinstance(c, str):
            out.append(float(Fraction(c.strip())))
        else:
            out.append(float(c))
    return np.asarray(out, dtype=float)


def _k_to_fractions(
    coords: Sequence[str | float | int | Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    """Parse one reciprocal vector without losing stated rational values."""
    if len(coords) != 3:
        raise ValueError("A reciprocal vector must contain exactly three components")
    values: list[Fraction] = []
    for value in coords:
        if isinstance(value, Fraction):
            values.append(value)
        elif isinstance(value, Integral):
            values.append(Fraction(int(value), 1))
        elif isinstance(value, str):
            values.append(Fraction(value.strip()))
        elif isinstance(value, Real):
            number = float(value)
            if not math.isfinite(number):
                raise ValueError(f"Coordinate must be finite, got {value!r}")
            # A binary float no longer carries the user's exact textual intent.
            # Preserve its round-trippable decimal value; do not guess a nearby
            # small-denominator fraction here.
            values.append(Fraction(repr(number)))
        else:
            raise TypeError(f"Unsupported reciprocal coordinate: {value!r}")
    return values[0], values[1], values[2]


def _fraction_matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )


def _k_equivalent_exact(
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
    centering: str,
) -> bool:
    """Exact reciprocal-lattice equivalence in conventional coordinates."""
    delta = tuple(a - b for a, b in zip(left, right, strict=True))
    if any(value.denominator != 1 for value in delta):
        return False
    h, k, ell = (int(value) for value in delta)
    return _g_allowed(h, k, ell, centering)


def _k_score_exact(vector: tuple[Fraction, Fraction, Fraction]) -> tuple:
    out_of_unit = sum(value < 0 or value > 1 for value in vector)
    negative = sum(value < 0 for value in vector)
    return out_of_unit, negative, sum(abs(value) for value in vector)


def _canonical_k_exact(
    vector: tuple[Fraction, Fraction, Fraction],
    centering: str,
    max_g: int = 2,
) -> tuple[Fraction, Fraction, Fraction]:
    best = vector
    best_score = _k_score_exact(best)
    for h in range(-max_g, max_g + 1):
        for k in range(-max_g, max_g + 1):
            for ell in range(-max_g, max_g + 1):
                if not _g_allowed(h, k, ell, centering):
                    continue
                candidate = (
                    vector[0] + h,
                    vector[1] + k,
                    vector[2] + ell,
                )
                score = _k_score_exact(candidate)
                if score < best_score:
                    best = candidate
                    best_score = score
    return best


def _zero_one_flips_exact(
    vector: tuple[Fraction, Fraction, Fraction],
) -> list[tuple[Fraction, Fraction, Fraction]]:
    result: list[tuple[Fraction, Fraction, Fraction]] = []
    for index, value in enumerate(vector):
        if value not in {Fraction(0), Fraction(1)}:
            continue
        flipped = list(vector)
        flipped[index] = Fraction(1) if value == 0 else Fraction(0)
        result.append((flipped[0], flipped[1], flipped[2]))
    return result


def k_star_tuples(k_coordinates: Sequence[str | float],
                  parent_sg: int) -> list[str]:
    """Star of ``k`` as ``(x,y,z)`` strings in a stable, website-like order.

    The listed seed is kept (not reduced by centering). Remaining arms are
    ordered as on the ISODISTORT OPD page: 0↔1 partners of the seed first,
    then leftover arms sorted by coordinates. For I4/mmm this yields N =
    ``(1/2,0,1/2),(1/2,1,1/2),(0,1/2,1/2),(1,1/2,1/2)``.
    """
    if not k_coordinates:
        return ["(0,0,0)"]
    centering = _centering_letter(parent_sg)
    k0 = _k_to_fractions(k_coordinates)
    sg = SpaceGroup.from_int_number(int(parent_sg))
    arms: list[tuple[Fraction, Fraction, Fraction]] = [k0]
    for op in sg.symmetry_ops:
        # Fractional direct coordinates transform as x' = R x; reciprocal
        # coordinates therefore transform by (R^-1)^T.  Space-group rotation
        # matrices are integral/unimodular, so the complete star stays exact.
        rotation = rational_matrix(np.asarray(op.rotation_matrix).tolist())
        reciprocal_rotation = transpose(inverse(rotation))
        kp = _canonical_k_exact(
            _fraction_matrix_vector(reciprocal_rotation, k0), centering,
        )
        if any(_k_equivalent_exact(kp, seen, centering) for seen in arms):
            continue
        arms.append(kp)

    ordered: list[tuple[Fraction, Fraction, Fraction]] = []

    def _add(vec: tuple[Fraction, Fraction, Fraction]) -> None:
        if any(_k_equivalent_exact(vec, seen, centering) for seen in ordered):
            return
        for arm in arms:
            if _k_equivalent_exact(vec, arm, centering):
                ordered.append(arm)
                return

    _add(k0)
    i = 0
    while i < len(ordered):
        for flip in _zero_one_flips_exact(ordered[i]):
            _add(flip)
        i += 1
    rest = [
        arm for arm in arms
        if not any(_k_equivalent_exact(arm, seen, centering) for seen in ordered)
    ]
    rest.sort()
    ordered.extend(rest)

    result = [format_tuple(k_coordinates)]
    for arm in ordered[1:]:
        result.append(format_tuple(arm))
    return result or [format_tuple(k_coordinates)]


def format_k_active(opd_dir_raw: str,
                    k_coordinates: Sequence[str | float],
                    parent_sg: int | None = None,
                    k_active_raw: str | None = None) -> str:
    """`` (0,0,0)`` or `` (1/2,1/2,0),(1/2,1/2,1)`` (leading space, no extra gaps)."""
    if k_active_raw:
        text = k_active_raw.strip()
        return text if text.startswith(" ") else f" {text}"
    star = k_star_tuples(k_coordinates, parent_sg or 1)
    inner = (opd_dir_raw or "").strip()
    if inner.startswith("(") and inner.endswith(")"):
        inner = inner[1:-1]
    arms = inner.split(";") if inner else [""]
    if len(arms) <= 1 or len(star) == 1:
        return " " + star[0]
    if len(arms) == len(star):
        active = [star[i] for i, arm in enumerate(arms) if _opd_arm_nonzero(arm)]
        if not active:
            active = [star[0]]
        return " " + ",".join(active)
    # Unequal lengths: keep arms that look nonzero, zip with the star.
    active = []
    for i, arm in enumerate(arms):
        if i < len(star) and _opd_arm_nonzero(arm):
            active.append(star[i])
    if not active:
        active = list(star)
    return " " + ",".join(active)


def official_method1_fields(
    *,
    irrep_label: str,
    opd_symbol: str,
    opd_dir_raw: str,
    space_group_number: int,
    space_group_symbol: str,
    basis_raw: str,
    origin_raw: str,
    size: int,
    subgroup_index: int,
    k_coordinates: Sequence[str | float] | None = None,
    parent_sg: int | None = None,
    k_active_raw: str | None = None,
    basis_vectors: Sequence[Sequence[float]] | None = None,
    origin: Sequence[float] | None = None,
) -> dict[str, str | int]:
    """Split the official Method 1 radio line into table columns.

    Visible radio text on the OPD page is one long line; the local UI keeps
    the same tokens but shows them as a filterable/sortable table.
    """
    direction = opd_dir_raw.strip() if opd_dir_raw else "(a)"
    basis_inner = format_basis(basis_vectors or [], basis_raw)
    orig = format_tuple(origin or [0, 0, 0], origin_raw)
    kact = format_k_active(
        direction,
        k_coordinates or ["0", "0", "0"],
        parent_sg,
        k_active_raw,
    ).strip()
    number = int(space_group_number)
    symbol = (space_group_symbol or "").strip()
    return {
        "irrep": (irrep_label or "").strip(),
        "opd": (opd_symbol or "").strip(),
        "dir": direction,
        "space_group_number": number,
        "space_group_symbol": symbol,
        "sg": f"{number} {symbol}".strip(),
        "basis": f"{{{basis_inner}}}",
        "origin": orig,
        "s": int(size),
        "i": int(subgroup_index),
        "k_active": kact,
    }


def format_opd_line_body(
    *,
    irrep_label: str = "",
    opd_symbol: str,
    opd_dir_raw: str,
    space_group_number: int,
    space_group_symbol: str,
    basis_raw: str,
    origin_raw: str,
    size: int,
    subgroup_index: int,
    k_coordinates: Sequence[str | float] | None = None,
    parent_sg: int | None = None,
    k_active_raw: str | None = None,
    basis_vectors: Sequence[Sequence[float]] | None = None,
    origin: Sequence[float] | None = None,
    pad_opd: bool = True,
) -> str:
    """OPD body without the irrep prefix (CIF Distortion comment style).

    Method 1 radio / CIF Method-1 comments pad the OPD symbol to 5 columns
    (``P1   (a)``). Method 2 CIF comments use a single space (``C1 (a,b)``).
    """
    fields = official_method1_fields(
        irrep_label=irrep_label,
        opd_symbol=opd_symbol,
        opd_dir_raw=opd_dir_raw,
        space_group_number=space_group_number,
        space_group_symbol=space_group_symbol,
        basis_raw=basis_raw,
        origin_raw=origin_raw,
        size=size,
        subgroup_index=subgroup_index,
        k_coordinates=k_coordinates,
        parent_sg=parent_sg,
        k_active_raw=k_active_raw,
        basis_vectors=basis_vectors,
        origin=origin,
    )
    opd = (
        f"{fields['opd']!s:<{_OPD_WIDTH}s}"
        if pad_opd
        else f"{fields['opd']!s} "
    )
    return (
        f"{opd}{fields['dir']} {int(fields['space_group_number']):>3d} "
        f"{fields['space_group_symbol']}, basis={fields['basis']}, origin={fields['origin']}, "
        f"s={fields['s']}, i={fields['i']}, k-active= {fields['k_active']}"
    )


def format_opd_line(
    *,
    irrep_label: str,
    opd_symbol: str,
    opd_dir_raw: str,
    space_group_number: int,
    space_group_symbol: str,
    basis_raw: str,
    origin_raw: str,
    size: int,
    subgroup_index: int,
    k_coordinates: Sequence[str | float] | None = None,
    parent_sg: int | None = None,
    k_active_raw: str | None = None,
    basis_vectors: Sequence[Sequence[float]] | None = None,
    origin: Sequence[float] | None = None,
) -> str:
    """One visible Method 1 radio line (no maximal asterisk)."""
    fields = official_method1_fields(
        irrep_label=irrep_label,
        opd_symbol=opd_symbol,
        opd_dir_raw=opd_dir_raw,
        space_group_number=space_group_number,
        space_group_symbol=space_group_symbol,
        basis_raw=basis_raw,
        origin_raw=origin_raw,
        size=size,
        subgroup_index=subgroup_index,
        k_coordinates=k_coordinates,
        parent_sg=parent_sg,
        k_active_raw=k_active_raw,
        basis_vectors=basis_vectors,
        origin=origin,
    )
    irrep = f"{fields['irrep']!s:<{_IRREP_WIDTH}s}"
    body = format_opd_line_body(
        irrep_label=irrep_label,
        opd_symbol=opd_symbol,
        opd_dir_raw=opd_dir_raw,
        space_group_number=space_group_number,
        space_group_symbol=space_group_symbol,
        basis_raw=basis_raw,
        origin_raw=origin_raw,
        size=size,
        subgroup_index=subgroup_index,
        k_coordinates=k_coordinates,
        parent_sg=parent_sg,
        k_active_raw=k_active_raw,
        basis_vectors=basis_vectors,
        origin=origin,
        pad_opd=True,
    )
    return f"{irrep}{body}"
