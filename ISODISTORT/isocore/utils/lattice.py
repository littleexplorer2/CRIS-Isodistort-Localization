"""Exact lattice arithmetic shared by the ISODISTORT search workflows.

Matrices use the convention adopted by the ISOTROPY text interface: every
row is one direct-lattice basis vector expressed in parent conventional
coordinates.  Fractions are kept exact so that centering translations such
as 1/2 and 1/3 are never classified by a floating-point tolerance.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from fractions import Fraction
from numbers import Integral, Real

Rational = Fraction
RationalMatrix = tuple[tuple[Fraction, Fraction, Fraction], ...]


def as_fraction(value: str | int | float | Fraction) -> Fraction:
    """Parse one finite scalar without changing its stated exact value.

    Text input is user/API data, so ``"1/2000001"`` and
    ``"0.5000000001"`` must not be rounded to a nearby small-denominator
    fraction.  Legacy numerical callers already holding binary floats get a
    high-denominator rational reconstruction (needed for values such as
    ``1/3``); callers that require an exact stated value should pass text or a
    :class:`Fraction` before any float conversion occurs.
    """
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    if isinstance(value, Real):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Lattice value must be finite, got {value!r}")
        return Fraction(str(number)).limit_denominator(1_000_000_000_000)
    text = str(value).strip()
    if not text:
        raise ValueError("Lattice value cannot be empty")
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"Invalid rational lattice value: {value!r}") from exc


def rational_matrix(
    rows: Sequence[Sequence[str | int | float | Fraction]],
) -> RationalMatrix:
    """Return a validated nonsingular 3x3 rational matrix."""
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        raise ValueError("Basis matrix must be 3x3")
    matrix: RationalMatrix = tuple(
        tuple(as_fraction(value) for value in row)  # type: ignore[misc]
        for row in rows
    )
    if determinant(matrix) == 0:
        raise ValueError("Basis vectors must be linearly independent")
    return matrix


def identity_matrix() -> RationalMatrix:
    return (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )


def transpose(matrix: RationalMatrix) -> RationalMatrix:
    return tuple(tuple(matrix[j][i] for j in range(3)) for i in range(3))


def multiply(left: RationalMatrix, right: RationalMatrix) -> RationalMatrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )


def determinant(matrix: RationalMatrix) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def inverse(matrix: RationalMatrix) -> RationalMatrix:
    det = determinant(matrix)
    if det == 0:
        raise ValueError("Basis vectors must be linearly independent")
    a, b, c = matrix
    cofactors = (
        (
            b[1] * c[2] - b[2] * c[1],
            -(b[0] * c[2] - b[2] * c[0]),
            b[0] * c[1] - b[1] * c[0],
        ),
        (
            -(a[1] * c[2] - a[2] * c[1]),
            a[0] * c[2] - a[2] * c[0],
            -(a[0] * c[1] - a[1] * c[0]),
        ),
        (
            a[1] * b[2] - a[2] * b[1],
            -(a[0] * b[2] - a[2] * b[0]),
            a[0] * b[1] - a[1] * b[0],
        ),
    )
    adjugate = transpose(cofactors)
    return tuple(tuple(value / det for value in row) for row in adjugate)


def is_integral(matrix: RationalMatrix) -> bool:
    return all(value.denominator == 1 for row in matrix for value in row)


def is_sublattice(child: RationalMatrix, parent: RationalMatrix) -> bool:
    """Whether every row translation of ``child`` belongs to ``parent``."""
    return is_integral(multiply(child, inverse(parent)))


def same_lattice(left: RationalMatrix, right: RationalMatrix) -> bool:
    """Whether two bases generate the same lattice (GL(3,Z) equivalence)."""
    transform = multiply(left, inverse(right))
    return is_integral(transform) and abs(determinant(transform)) == 1


def same_lattice_in_point_group_orbit(
    candidate: RationalMatrix,
    selected: RationalMatrix,
    parent_rotations: Sequence[Sequence[Sequence[str | int | float | Fraction]]] | None,
) -> bool:
    """Compare lattice classes after all parent point-group rotations.

    Fractional coordinates transform as column vectors ``x' = R x``.  With
    row-vector bases this is ``B' = B R.T``.  The identity is always tested,
    even when a caller supplies an incomplete rotation list.
    """
    orientations = [selected]
    for rotation in parent_rotations or []:
        rot = rational_matrix(rotation)
        orientations.append(multiply(selected, transpose(rot)))
    return any(same_lattice(candidate, oriented) for oriented in orientations)


def centering_primitive_matrix(letter: str) -> RationalMatrix:
    """Conventional-to-primitive row basis for a Bravais centering symbol."""
    half = Fraction(1, 2)
    third = Fraction(1, 3)
    tables: dict[str, RationalMatrix] = {
        "P": identity_matrix(),
        "I": ((-half, half, half), (half, -half, half), (half, half, -half)),
        "F": ((0, half, half), (half, 0, half), (half, half, 0)),
        "A": ((1, 0, 0), (0, half, half), (0, -half, half)),
        "B": ((half, 0, half), (0, 1, 0), (-half, 0, half)),
        "C": ((half, half, 0), (-half, half, 0), (0, 0, 1)),
        "R": (
            (2 * third, third, third),
            (-third, third, third),
            (-third, -2 * third, third),
        ),
    }
    key = (letter or "P").strip().upper()
    if key not in tables:
        raise ValueError(f"Unknown centering symbol: {letter!r}")
    return tables[key]


def to_float_rows(matrix: RationalMatrix) -> list[list[float]]:
    return [[float(value) for value in row] for row in matrix]


def matrix_key(matrix: RationalMatrix) -> tuple[tuple[tuple[int, int], ...], ...]:
    """Stable hashable representation preserving exact signs and fractions."""
    return tuple(
        tuple((value.numerator, value.denominator) for value in row)
        for row in matrix
    )
