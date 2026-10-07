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
from itertools import combinations, product
from numbers import Integral, Real

Rational = Fraction
RationalMatrix = tuple[tuple[Fraction, Fraction, Fraction], ...]
RationalVector = tuple[Fraction, Fraction, Fraction]


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


def rational_vector(
    values: Sequence[str | int | float | Fraction],
) -> RationalVector:
    """Validate a three-component vector without a denominator/period cap."""
    if len(values) != 3:
        raise ValueError("wavevector must contain exactly three components")
    return tuple(as_wavevector_fraction(value) for value in values)  # type: ignore[return-value]


def as_wavevector_fraction(value: str | int | float | Fraction) -> Fraction:
    """Recover rational k without rounding a float into a real phase class.

    Ordinary fractions such as a float representation of 1/3 retain rational
    recovery. A reconstruction onto an integer or half-integer can turn a
    genuinely complex character into Gamma or a sign-only character. Preserve
    the stated decimal float in that case; exact text/Fraction values are
    unchanged. Matrix rationalization keeps its separate numerical contract.
    """
    recovered = as_fraction(value)
    if isinstance(value, Real) and not isinstance(value, Integral):
        stated = Fraction(str(float(value)))
        if (2 * recovered).denominator == 1 and stated != recovered:
            return stated
    return recovered


def crystallographic_fraction(value: object) -> Fraction:
    """Recover a spglib rational only within an IEEE-754 roundoff budget.

    This boundary applies to computed symmetry operations, not exact user
    basis/wavevector text. The denominator bound is a reconstruction proposal;
    the residual test, rather than proximity to a desired fraction, certifies
    it. Values outside that machine-error budget fail closed.
    """
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite spglib symmetry value: {value!r}")
    candidate = Fraction(number).limit_denominator(1_000_000)
    tolerance = 128 * math.ulp(1.0) * max(1.0, abs(number))
    if abs(float(candidate) - number) > tolerance:
        raise ValueError(
            f"spglib value {number!r} is not a verified crystallographic rational"
        )
    return candidate


def wavevector_translation_phases(
    wavevector: Sequence[str | int | float | Fraction],
    translations: Sequence[Sequence[str | int | float | Fraction]],
) -> tuple[Fraction, ...]:
    """Exact ``k·t`` phases in one common fractional direct/reciprocal frame.

    The character is ``exp(2πi k·t)``. It is real/self-conjugate on a
    translation lattice iff twice every primitive-generator phase is an
    integer; Gamma requires the phases themselves to be integers. Using the
    actual primitive translations includes centered-cell characters.
    """
    k = rational_vector(wavevector)
    return tuple(
        sum(component * value for component, value in zip(
            rational_vector(translation), k, strict=True,
        ))
        for translation in translations
    )


def real_translation_character_sign(
    wavevector: Sequence[str | int | float | Fraction],
    translation: Sequence[str | int | float | Fraction],
) -> int:
    """Return an exact ±1 character, rejecting a complex translation phase."""
    phase, = wavevector_translation_phases(wavevector, [translation])
    if (2 * phase).denominator != 1:
        raise ValueError("translation character requires paired real/complex modes")
    return 1 if phase.denominator == 1 else -1


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


def translation_lattice_from_cosets(
    translations: Sequence[Sequence[str | int | float | Fraction]],
) -> RationalMatrix:
    """Recover the primitive row lattice from all pure-translation cosets.

    The actual uploaded cell contributes ``Z^3``. Identity-rotation dataset
    operations provide the finite cosets modulo that cell; a conventional HM
    centering letter is insufficient for a nonstandard or supercell input.
    The quotient order fixes ``abs(det(P)) = 1 / number_of_cosets``. Every
    generator must belong to the recovered lattice, so determinant and exact
    generator membership together certify the basis. Input must contain the
    complete quotient, not an arbitrary subset of translation generators.
    """
    zero: RationalVector = (Fraction(0), Fraction(0), Fraction(0))
    cosets: list[RationalVector] = [zero]
    for translation in translations:
        if len(translation) != 3:
            raise ValueError("translation cosets must contain three components")
        reduced = tuple(as_fraction(value) % 1 for value in translation)
        if reduced not in cosets:
            cosets.append(reduced)  # type: ignore[arg-type]
    expected_determinant = Fraction(1, len(cosets))
    candidates: set[RationalVector] = set()
    for coset in cosets:
        for shift in product((-1, 0, 1), repeat=3):
            vector = tuple(coset[index] + shift[index] for index in range(3))
            if vector != zero:
                candidates.add(vector)  # type: ignore[arg-type]
    ordered = sorted(
        candidates,
        key=lambda vector: (
            sum(value * value for value in vector),
            tuple((value.numerator, value.denominator) for value in vector),
        ),
    )
    generators = [*identity_matrix(), *cosets[1:]]
    for rows in combinations(ordered, 3):
        matrix: RationalMatrix = rows
        if abs(determinant(matrix)) != expected_determinant:
            continue
        matrix_inverse = inverse(matrix)
        if all(
            all(
                sum(item[index] * matrix_inverse[index][column] for index in range(3)).denominator == 1
                for column in range(3)
            )
            for item in generators
        ):
            return matrix
    raise ValueError("failed to reconstruct the uploaded cell's translation lattice")


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
