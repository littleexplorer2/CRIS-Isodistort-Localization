"""Format Method 1 order-parameter lines to match the official ISODISTORT OPD page.

Visible radio text on the website looks like::

    GM1+     P1   (a) 139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=(0,0,0), s=1, i=1, k-active= (0,0,0)
"""
from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

import numpy as np
from pymatgen.symmetry.groups import SpaceGroup

# Official radio line: irrep left-padded to 9, OPD symbol to 5.
_IRREP_WIDTH = 9
_OPD_WIDTH = 5


def format_number(value: float | int | str) -> str:
    """Render a coordinate as an integer or reduced fraction (``1/2``, not ``0.5``)."""
    if isinstance(value, str):
        text = value.strip()
        if text:
            return text
        return "0"
    frac = Fraction(value).limit_denominator(48)
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def format_tuple(values: Sequence[float | int | str], raw: str | None = None) -> str:
    """``(x,y,z)`` using the iso raw token when present."""
    if raw:
        text = raw.strip()
        if text.startswith("(") and text.endswith(")"):
            return text
        if text:
            return f"({text})"
    inner = ",".join(format_number(v) for v in values)
    return f"({inner})"


def format_basis(vectors: Sequence[Sequence[float | int | str]],
                 raw: str | None = None) -> str:
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
    """Lower is better: stay in [0,1], avoid negatives, then smaller |k|."""
    out_of_unit = int(sum(1 for x in k if x < -1e-8 or x > 1.0 + 1e-8))
    negative = int(sum(1 for x in k if x < -1e-8))
    abs_sum = float(np.sum(np.abs(k)))
    return (out_of_unit, negative, abs_sum)


def _canonical_k(k: np.ndarray, centering: str, max_g: int = 2) -> np.ndarray:
    """Pick the ISODISTORT-like representative of k modulo the reciprocal lattice."""
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


def k_star_tuples(k_coordinates: Sequence[str | float],
                  parent_sg: int) -> list[str]:
    """Star of ``k`` as ``(x,y,z)`` strings in a stable, website-like order."""
    if not k_coordinates:
        return ["(0,0,0)"]
    centering = _centering_letter(parent_sg)
    k0 = _canonical_k(_k_to_array(k_coordinates), centering)
    sg = SpaceGroup.from_int_number(int(parent_sg))
    arms: list[np.ndarray] = []
    for op in sg.symmetry_ops:
        kp = _canonical_k(np.asarray(op.rotation_matrix, dtype=float) @ k0, centering)
        if any(_k_equivalent(kp, seen, centering) for seen in arms):
            continue
        arms.append(kp)
    ordered: list[np.ndarray] = []
    for cand in [k0, *arms]:
        reduced = _canonical_k(cand, centering)
        if any(_k_equivalent(reduced, seen, centering) for seen in ordered):
            continue
        ordered.append(reduced)
    result = [format_tuple(arm.tolist()) for arm in ordered]
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
    irrep = f"{(irrep_label or '').strip():<{_IRREP_WIDTH}s}"
    opd = f"{(opd_symbol or '').strip():<{_OPD_WIDTH}s}"
    direction = opd_dir_raw.strip() if opd_dir_raw else "(a)"
    basis = format_basis(basis_vectors or [], basis_raw)
    orig = format_tuple(origin or [0, 0, 0], origin_raw)
    kact = format_k_active(
        direction,
        k_coordinates or ["0", "0", "0"],
        parent_sg,
        k_active_raw,
    )
    return (
        f"{irrep}{opd}{direction} {int(space_group_number):>3d} "
        f"{space_group_symbol}, basis={{{basis}}}, origin={orig}, "
        f"s={int(size)}, i={int(subgroup_index)}, k-active={kact}"
    )
