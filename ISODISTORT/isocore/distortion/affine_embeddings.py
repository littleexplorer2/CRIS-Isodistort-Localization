"""Exact affine-group primitives for the first stage of Method 3.

The first Method-3 page asks for affine subgroup *embeddings*.  An embedding
is a subgroup of the parent space group, not an IR/OPD label.  This module
therefore works only with Seitz operators and translation lattices.  It uses
the finite quotient construction of Stokes & Campbell, Acta Cryst. A73
(2017), Appendix B: for a selected child translation lattice ``T_s``, the
otherwise infinite search is represented by the finite group
``N_G(T_s)/T_s``.  The normalizer restriction is essential: an arbitrary
selected sublattice need not be normal in the full parent group.  Lifted
point-subgroup operations are accepted only when they are
closed under Seitz multiplication, which is the cocycle/factor-set condition

``q_RS = q_R + R q_S + (tau_R + R tau_S - tau_RS)  (mod T_s)``.

Matrices act on fractional-coordinate columns, ``x' = R x + t``.  Lattice
bases use the project's row-vector convention.  All internal arithmetic is
exact :class:`fractions.Fraction`; floats occur only at the spglib database
boundary and are accepted only after a machine-roundoff check.

This module deliberately does not attach IR labels.  The finite-quotient
enumerator is nevertheless a complete exact core *for a given point
subgroup*.  Its standardisation bridge treats spglib's Hall setting only as a
proposal: every applicable Hall choice is reconstructed with exact Seitz
operators, then converted to the project's preferred ISODISTORT setting and
verified again before a basis/origin is returned.
"""

from __future__ import annotations

import math
import re
import sys
import warnings
from collections import deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from functools import cache, lru_cache
from itertools import combinations, permutations, product
from typing import Any

import numpy as np
import spglib
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ..io.isodistort_cif import _hall_for_number
from ..utils.config_loader import get_config
from ..utils.lattice import (
    RationalMatrix,
    as_fraction,
    centering_primitive_matrix,
    determinant,
    identity_matrix,
    inverse,
    is_integral,
    multiply,
    rational_matrix,
    same_lattice,
    transpose,
)

RationalVector = tuple[Fraction, Fraction, Fraction]


def _positive_tolerance(value: float, description: str) -> float:
    """Validate an explicitly unit-labelled numerical tolerance."""
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise ValueError(f"{description} must be a finite positive number")
    return parsed


def _effective_exact_cartesian_tolerance(
    lattice: Sequence[Sequence[float]], configured_floor_angstrom: float
) -> float:
    """Return the Å error budget for a rational-to-float spglib boundary.

    The configured value is an absolute floor.  IEEE-754 coordinate
    conversion error also grows with the largest singular value of the cell
    matrix, so retain a conservative ``256 * eps * ||A||_2`` guard.  This is
    numerical roundoff control for synthetic exact data, not permission to
    absorb physical atomic displacements.
    """
    floor = _positive_tolerance(
        configured_floor_angstrom, "exact affine Cartesian tolerance (Å)"
    )
    matrix = np.asarray(lattice, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("Cartesian lattice must be a finite 3x3 matrix")
    scale_angstrom = max(1.0, float(np.linalg.norm(matrix, ord=2)))
    roundoff_guard = 256.0 * np.finfo(float).eps * scale_angstrom
    return max(floor, roundoff_guard)


@dataclass(frozen=True)
class AffineOperation:
    """One Seitz operation ``(rotation, translation)``."""

    rotation: RationalMatrix
    translation: RationalVector


@dataclass(frozen=True)
class AffineEmbedding:
    """A subgroup represented modulo its primitive translation lattice."""

    lattice: RationalMatrix
    operations: tuple[AffineOperation, ...]
    hall_number: int | None = None


@dataclass(frozen=True)
class ParentAffineGroup:
    """Parent affine group in the uploaded structure's fractional axes."""

    lattice: RationalMatrix
    operations: tuple[AffineOperation, ...]
    space_group_number: int
    hall_number: int
    symbol: str
    # Rows of the spglib-standard conventional parent basis expressed in the
    # uploaded fractional axes.  Method-3 form values and the bundled ISO
    # tables use that standard frame, whereas all exact affine operations in
    # this object use the uploaded frame.
    standard_basis: RationalMatrix | None = None
    # Uploaded fractional coordinates of the standard-setting origin.
    standard_origin: RationalVector = (
        Fraction(0),
        Fraction(0),
        Fraction(0),
    )


@dataclass(frozen=True)
class AffineQuotient:
    """Exact finite quotient ``N_G(T_s) / T_s`` for one child lattice.

    ``operations`` are immutable Seitz-coset representatives in the uploaded
    parent axes.  ``multiplication_table`` stores indices into that tuple and
    was constructed by exact Seitz multiplication modulo ``lattice``.
    """

    lattice: RationalMatrix
    operations: tuple[AffineOperation, ...]
    multiplication_table: tuple[tuple[int, ...], ...]
    identity_index: int
    translation_index: int

    def __post_init__(self) -> None:
        count = len(self.operations)
        if count == 0:
            raise ValueError("an affine quotient must contain operations")
        if determinant(self.lattice) == 0:
            raise ValueError("affine quotient lattice must be nonsingular")
        if len(self.multiplication_table) != count or any(
            len(row) != count for row in self.multiplication_table
        ):
            raise ValueError("affine quotient multiplication table has the wrong shape")
        if any(
            not isinstance(value, int) or isinstance(value, bool) or not 0 <= value < count
            for row in self.multiplication_table
            for value in row
        ):
            raise ValueError("affine quotient table contains an invalid index")
        if not 0 <= self.identity_index < count:
            raise ValueError("affine quotient identity index is invalid")
        if self.translation_index < 1:
            raise ValueError("affine quotient translation index must be positive")

    @property
    def identity(self) -> AffineOperation:
        return self.operations[self.identity_index]

    def locate(self, operation: AffineOperation) -> int:
        """Locate an operation coset exactly modulo the child lattice."""
        lattice_inverse = inverse(self.lattice)
        operation_key = _operation_coset_key(operation, lattice_inverse)
        for index, candidate in enumerate(self.operations):
            if _operation_coset_key(candidate, lattice_inverse) == operation_key:
                return index
        raise ValueError("operation is not contained in the affine quotient")


@dataclass(frozen=True)
class _AffineQuotientPreflight:
    lattice: RationalMatrix
    normalizer_operations: tuple[AffineOperation, ...]
    translation_index: int
    quotient_order: int


def _field(identity: object, *names: str) -> Any:
    for name in names:
        if isinstance(identity, Mapping) and name in identity:
            return identity[name]
        if hasattr(identity, name):
            return getattr(identity, name)
    raise KeyError(names[0])


def _optional_field(identity: object, *names: str) -> Any | None:
    """Return the first present field without conflating absence with errors."""

    for name in names:
        if isinstance(identity, Mapping) and name in identity:
            return identity[name]
        if hasattr(identity, name):
            return getattr(identity, name)
    return None


def _database_fraction(value: object) -> Fraction:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite spglib symmetry value: {value!r}")
    candidate = Fraction(number).limit_denominator(1_000_000)
    tolerance = 128 * sys.float_info.epsilon * max(1.0, abs(number))
    if abs(float(candidate) - number) > tolerance:
        raise ValueError(
            f"spglib value {number!r} is not a verified crystallographic rational"
        )
    return candidate


def _centering(symbol: object) -> str:
    match = re.search(r"[PABCIFR]", str(symbol or ""), flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"cannot determine Bravais centering from {symbol!r}")
    return match.group(0).upper()


def _add(left: RationalVector, right: RationalVector) -> RationalVector:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def _subtract(left: RationalVector, right: RationalVector) -> RationalVector:
    return tuple(left[index] - right[index] for index in range(3))  # type: ignore[return-value]


def _matrix_vector(matrix: RationalMatrix, vector: RationalVector) -> RationalVector:
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _row_matrix(vector: RationalVector, matrix: RationalMatrix) -> RationalVector:
    return tuple(
        sum(vector[row] * matrix[row][column] for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def vector_in_lattice(vector: RationalVector, lattice: RationalMatrix) -> bool:
    coordinates = _row_matrix(vector, inverse(lattice))
    return all(value.denominator == 1 for value in coordinates)


def operations_equivalent(
    left: AffineOperation,
    right: AffineOperation,
    lattice: RationalMatrix,
) -> bool:
    """Whether two Seitz operators are equal modulo ``lattice``."""

    return left.rotation == right.rotation and vector_in_lattice(
        _subtract(left.translation, right.translation), lattice
    )


def _operation_coset_key(
    operation: AffineOperation,
    lattice_inverse: RationalMatrix,
) -> tuple[RationalMatrix, RationalVector]:
    """Return the exact canonical key of a Seitz coset modulo a lattice.

    With row lattice bases, ``t @ lattice_inverse`` are the lattice
    coordinates of a translation.  Componentwise reduction modulo one is
    therefore exactly the quotient by the selected translation lattice; no
    floating tolerance or choice of coset representative enters the key.
    """

    coordinates = _row_matrix(operation.translation, lattice_inverse)
    reduced = tuple(value % 1 for value in coordinates)
    return operation.rotation, reduced  # type: ignore[return-value]


def _deduplicate_operations(
    operations: Sequence[AffineOperation], lattice: RationalMatrix
) -> tuple[AffineOperation, ...]:
    lattice_inverse = inverse(lattice)
    result: list[AffineOperation] = []
    seen: set[tuple[RationalMatrix, RationalVector]] = set()
    for operation in operations:
        key = _operation_coset_key(operation, lattice_inverse)
        if key not in seen:
            seen.add(key)
            result.append(operation)
    return tuple(result)


def _operation_sets_equal(
    left: Sequence[AffineOperation],
    right: Sequence[AffineOperation],
    lattice: RationalMatrix,
) -> bool:
    if len(left) != len(right):
        return False
    lattice_inverse = inverse(lattice)
    left_keys = {
        _operation_coset_key(operation, lattice_inverse) for operation in left
    }
    right_keys = {
        _operation_coset_key(operation, lattice_inverse) for operation in right
    }
    return len(left_keys) == len(left) and left_keys == right_keys


@cache
def _hall_numbers_for_space_group(space_group_number: int) -> tuple[int, ...]:
    """Return every spglib Hall setting for one international SG type."""

    result: list[int] = []
    for hall_number in range(1, 531):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            hall_type = spglib.get_spacegroup_type(hall_number)
        if hall_type is not None and int(hall_type.number) == int(space_group_number):
            result.append(hall_number)
    return tuple(result)


@cache
def _hall_database_operations(hall_number: int) -> tuple[AffineOperation, ...]:
    """Read one Hall setting from spglib and recover all values exactly."""

    database = spglib.get_symmetry_from_database(int(hall_number))
    if database is None:
        raise ValueError(f"unknown Hall number {hall_number}")
    operations = tuple(
        AffineOperation(
            rational_matrix(rotation.tolist()),
            tuple(_database_fraction(value) for value in translation),  # type: ignore[arg-type]
        )
        for rotation, translation in zip(
            database["rotations"], database["translations"], strict=True
        )
    )
    return _deduplicate_operations(operations, identity_matrix())


@lru_cache(maxsize=1)
def _signed_permutation_matrices() -> tuple[RationalMatrix, ...]:
    """Crystallographic conventional-axis relabellings, identity first."""

    matrices: list[RationalMatrix] = []
    for permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            rows = []
            for row in range(3):
                rows.append(
                    tuple(
                        Fraction(signs[row] if column == permutation[row] else 0)
                        for column in range(3)
                    )
                )
            matrices.append(tuple(rows))  # type: ignore[arg-type]
    identity = identity_matrix()
    return tuple(
        sorted(
            set(matrices),
            key=lambda matrix: (
                matrix != identity,
                tuple(value for row in matrix for value in row),
            ),
        )
    )


def _transform_hall_operations(
    operations: Sequence[AffineOperation],
    source_from_target: RationalMatrix,
    source_origin_of_target: RationalVector,
) -> tuple[AffineOperation, ...]:
    """Express source-Hall operators in target-Hall coordinates exactly.

    The affine coordinate relation is ``x_source = C x_target + d``.  Thus
    ``R_target = C^-1 R_source C`` and
    ``t_target = C^-1 (t_source + (R_source-I)d)``.
    """

    target_from_source = inverse(source_from_target)
    identity = identity_matrix()
    transformed = []
    for operation in operations:
        rotation = multiply(
            multiply(target_from_source, operation.rotation),
            source_from_target,
        )
        rotation_minus_identity: RationalMatrix = tuple(
            tuple(operation.rotation[row][column] - identity[row][column]
                  for column in range(3))
            for row in range(3)
        )
        translation = _matrix_vector(
            target_from_source,
            _add(
                operation.translation,
                _matrix_vector(rotation_minus_identity, source_origin_of_target),
            ),
        )
        transformed.append(AffineOperation(rotation, translation))
    return _deduplicate_operations(transformed, identity)


@cache
def _hall_setting_transforms(
    source_hall_number: int,
    target_hall_number: int,
) -> tuple[tuple[RationalMatrix, RationalVector], ...]:
    """Find exact conventional-setting transforms between two Hall choices.

    Hall choices of one IT space-group type commonly differ by an origin
    choice and/or a signed conventional-axis permutation.  Candidate axis
    maps are finite.  Origin shifts are solved from

    ``(R-I)d = C t_target - t_source  (mod Z^3)``.

    When ``R-I`` is nonsingular, its determinant ``D`` makes
    ``D Z^3`` a sublattice of ``(R-I) Z^3`` (adjugate identity), so enumerating
    integer right-hand sides modulo ``D`` is complete.  Every proposal is
    finally checked against the complete exact Seitz-operation sets; no
    numerical or symbol-only match is accepted.
    """

    source_type = spglib.get_spacegroup_type(int(source_hall_number))
    target_type = spglib.get_spacegroup_type(int(target_hall_number))
    if (
        source_type is None
        or target_type is None
        or int(source_type.number) != int(target_type.number)
    ):
        return ()
    source = _hall_database_operations(int(source_hall_number))
    target = _hall_database_operations(int(target_hall_number))
    target_rotations = {operation.rotation for operation in target}
    identity = identity_matrix()
    zero: RationalVector = (Fraction(0), Fraction(0), Fraction(0))

    results: list[tuple[RationalMatrix, RationalVector]] = []
    for source_from_target in _signed_permutation_matrices():
        target_from_source = inverse(source_from_target)
        transformed_rotations = {
            multiply(
                multiply(target_from_source, operation.rotation),
                source_from_target,
            )
            for operation in source
        }
        if transformed_rotations != target_rotations:
            continue

        shifts: list[RationalVector] = [zero]
        for source_operation in source:
            rotation_minus_identity: RationalMatrix = tuple(
                tuple(
                    source_operation.rotation[row][column] - identity[row][column]
                    for column in range(3)
                )
                for row in range(3)
            )
            determinant_value = abs(determinant(rotation_minus_identity))
            if determinant_value == 0:
                continue
            if determinant_value.denominator != 1:
                raise ValueError("Hall rotation has a non-integral determinant")
            transformed_rotation = multiply(
                multiply(target_from_source, source_operation.rotation),
                source_from_target,
            )
            matching_targets = [
                operation
                for operation in target
                if operation.rotation == transformed_rotation
            ]
            modulus = int(determinant_value)
            inverse_equation = inverse(rotation_minus_identity)
            for target_operation in matching_targets:
                right_hand_side = _subtract(
                    _matrix_vector(
                        source_from_target, target_operation.translation
                    ),
                    source_operation.translation,
                )
                for integer_shift in product(range(modulus), repeat=3):
                    shifted = tuple(
                        right_hand_side[index] + integer_shift[index]
                        for index in range(3)
                    )
                    candidate = tuple(
                        value % 1
                        for value in _matrix_vector(
                            inverse_equation, shifted  # type: ignore[arg-type]
                        )
                    )
                    if candidate not in shifts:
                        shifts.append(candidate)  # type: ignore[arg-type]
            # One nonsingular equation generates every shift class.  The
            # complete-operation verification below rejects incompatible ones.
            break

        valid_for_axes = []
        for shift in shifts:
            transformed = _transform_hall_operations(
                source, source_from_target, shift
            )
            if _operation_sets_equal(transformed, target, identity):
                valid_for_axes.append((source_from_target, shift))
        if valid_for_axes:
            results.extend(valid_for_axes)
            # Preserve a standardizer's axes whenever a pure origin-choice
            # conversion exists.  Axis permutations are a fallback for Hall
            # choices that genuinely use different conventional axes.
            if source_from_target == identity:
                break

    return tuple(results)


def _translation_lattice_from_operations(
    operations: Sequence[AffineOperation],
) -> RationalMatrix:
    """Derive the pure-translation lattice in the uploaded fractional axes.

    The international HM centering letter describes a conventional standard
    setting, not necessarily the uploaded cell.  Here ``Z^3`` translations of
    the actual input cell are augmented by every dataset operation whose
    rotation is identity.  The finite quotient size fixes the primitive-cell
    determinant exactly.  A short crystallographic coset search then selects
    a three-vector basis that contains every generator.
    """

    identity = identity_matrix()
    zero: RationalVector = (Fraction(0), Fraction(0), Fraction(0))
    cosets: list[RationalVector] = [zero]
    for operation in operations:
        if operation.rotation != identity:
            continue
        reduced = tuple(value % 1 for value in operation.translation)
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
    generators = [
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
        *cosets[1:],
    ]
    for rows in combinations(ordered, 3):
        matrix: RationalMatrix = rows
        if abs(determinant(matrix)) != expected_determinant:
            continue
        matrix_inverse = inverse(matrix)
        if all(
            all(value.denominator == 1 for value in _row_matrix(item, matrix_inverse))
            for item in generators
        ):
            return matrix
    raise ValueError("failed to reconstruct the uploaded cell's translation lattice")


def parent_affine_group(
    structure: Structure,
    symprec: float | None = None,
    angle_tolerance_degrees: float | None = None,
) -> ParentAffineGroup:
    """Detect and exactly reconstruct the parent group in uploaded axes.

    ``symprec`` is an absolute Cartesian distance in Å.  Fractional
    coordinates are deliberately not compared against it directly: spglib
    converts them with the supplied lattice metric.  Once detected, all group
    reconstruction below is exact rational arithmetic.
    """

    cfg = get_config()
    symprec_angstrom = _positive_tolerance(
        cfg.symmetry_cartesian_tolerance_angstrom
        if symprec is None
        else symprec,
        "parent symmetry Cartesian tolerance (Å)",
    )
    angle_degrees = _positive_tolerance(
        cfg.symmetry_angle_tolerance_degrees
        if angle_tolerance_degrees is None
        else angle_tolerance_degrees,
        "parent symmetry angle tolerance (degrees)",
    )
    dataset = SpacegroupAnalyzer(
        structure,
        symprec=symprec_angstrom,
        angle_tolerance=angle_degrees,
    ).get_symmetry_dataset()
    if dataset is None:
        raise ValueError("spglib could not determine the parent space group")
    symbol = str(dataset.international)
    raw_operations = [
        AffineOperation(
            rational_matrix(rotation.tolist()),
            tuple(_database_fraction(value) for value in translation),  # type: ignore[arg-type]
        )
        for rotation, translation in zip(
            dataset.rotations, dataset.translations, strict=True
        )
    ]
    lattice = _translation_lattice_from_operations(raw_operations)
    transformation = rational_matrix(
        [
            [_database_fraction(value) for value in row]
            for row in dataset.transformation_matrix
        ]
    )
    # spglib uses x_std = P x_input + p and column lattice vectors with
    # L_input = L_std P.  In this project's row-basis convention, one standard
    # conventional direct-lattice row is therefore expressed in uploaded axes
    # by S = P^{-T}.
    standard_basis = inverse(transpose(transformation))
    origin_shift: RationalVector = tuple(  # type: ignore[assignment]
        _database_fraction(value) for value in dataset.origin_shift
    )
    standard_origin = _matrix_vector(
        inverse(transformation),
        tuple(-value for value in origin_shift),  # type: ignore[arg-type]
    )
    standard_primitive = multiply(
        centering_primitive_matrix(_centering(symbol)),
        standard_basis,
    )
    if not same_lattice(standard_primitive, lattice):
        raise ValueError(
            "spglib standard-setting transform is inconsistent with the "
            "uploaded cell's exact translation lattice"
        )
    return ParentAffineGroup(
        lattice=lattice,
        operations=_deduplicate_operations(raw_operations, lattice),
        space_group_number=int(dataset.number),
        hall_number=int(dataset.hall_number),
        symbol=symbol,
        standard_basis=standard_basis,
        standard_origin=standard_origin,
    )


def embedding_from_identity(
    identity: object, parent: ParentAffineGroup
) -> AffineEmbedding:
    """Reconstruct a candidate's exact child operators in parent axes.

    Ordinary UI/ISO identities omit ``hall_number`` and therefore use the
    project's preferred Hall setting.  The standardisation bridge may supply
    an explicit Hall number while testing spglib's alternative origin/axis
    choices; it is accepted only when it belongs to the requested IT type.
    """

    basis = rational_matrix(_field(identity, "basis", "basis_vectors"))
    raw_origin = _field(identity, "origin") or (0, 0, 0)
    origin_values = tuple(as_fraction(value) for value in raw_origin)
    if len(origin_values) != 3:
        raise ValueError("embedding origin must contain three coordinates")
    origin: RationalVector = origin_values  # type: ignore[assignment]
    number = int(_field(identity, "space_group_number"))
    requested_hall_number = _optional_field(identity, "hall_number")
    if requested_hall_number is None:
        hall_number, _hall_symbol = _hall_for_number(number)
    else:
        hall_number = int(requested_hall_number)
    hall_type = spglib.get_spacegroup_type(hall_number)
    if hall_number <= 0 or hall_type is None:
        raise ValueError(f"no valid Hall setting for space group {number}")
    if int(hall_type.number) != number:
        raise ValueError(
            f"Hall number {hall_number} belongs to SG {hall_type.number}, not SG {number}"
        )
    # The Hall symbol, unlike an abbreviated/canonical HM symbol, records the
    # actual A/B/C/I/F/R conventional centering of this particular setting.
    child_lattice = multiply(
        centering_primitive_matrix(_centering(hall_type.hall_symbol)), basis
    )
    relative = multiply(child_lattice, inverse(parent.lattice))
    if not is_integral(relative):
        raise ValueError("embedding translation lattice is not a parent sublattice")

    database = spglib.get_symmetry_from_database(hall_number)
    parent_from_child = transpose(basis)
    child_from_parent = inverse(parent_from_child)
    transformed: list[AffineOperation] = []
    for raw_rotation, raw_translation in zip(
        database["rotations"], database["translations"], strict=True
    ):
        child_rotation = rational_matrix(raw_rotation.tolist())
        parent_rotation = multiply(
            multiply(parent_from_child, child_rotation), child_from_parent
        )
        if not is_integral(parent_rotation):
            raise ValueError(
                "preferred Hall rotation does not become integral in parent coordinates"
            )
        child_translation: RationalVector = tuple(  # type: ignore[assignment]
            _database_fraction(value) for value in raw_translation
        )
        parent_translation = _add(
            _matrix_vector(parent_from_child, child_translation),
            _subtract(origin, _matrix_vector(parent_rotation, origin)),
        )
        operation = AffineOperation(parent_rotation, parent_translation)
        if not any(
            operations_equivalent(operation, candidate, parent.lattice)
            for candidate in parent.operations
        ):
            raise ValueError(
                "reconstructed child operation is not contained in the parent space group"
            )
        transformed.append(operation)
    operations = _deduplicate_operations(transformed, child_lattice)
    if not operations:
        raise ValueError("embedding produced no affine operations")
    return AffineEmbedding(child_lattice, operations, hall_number)


def translation_cosets(
    parent_lattice: RationalMatrix, child_lattice: RationalMatrix
) -> tuple[RationalVector, ...]:
    """Enumerate the complete finite quotient ``T_parent/T_child``."""

    relative = multiply(child_lattice, inverse(parent_lattice))
    if not is_integral(relative):
        raise ValueError("child lattice is not a sublattice of the parent")
    index_value = abs(determinant(relative))
    if index_value.denominator != 1 or index_value < 1:
        raise ValueError("parent/child translation quotient has invalid index")
    index = int(index_value)
    zero: RationalVector = (Fraction(0), Fraction(0), Fraction(0))
    generators: tuple[RationalVector, ...] = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    representatives = [zero]
    queue: deque[RationalVector] = deque([zero])
    while queue and len(representatives) < index:
        current = queue.popleft()
        for generator in generators:
            candidate = _add(current, generator)
            if any(
                vector_in_lattice(_subtract(candidate, existing), relative)
                for existing in representatives
            ):
                continue
            representatives.append(candidate)
            queue.append(candidate)
            if len(representatives) == index:
                break
    if len(representatives) != index:
        raise ValueError("failed to enumerate the complete translation quotient")
    return tuple(_row_matrix(item, parent_lattice) for item in representatives)


def conjugate_embedding(
    source: AffineEmbedding,
    rotation: RationalMatrix,
    translation: RationalVector,
    target_lattice: RationalMatrix,
) -> AffineEmbedding:
    inverse_rotation = inverse(rotation)
    operations: list[AffineOperation] = []
    for operation in source.operations:
        transformed_rotation = multiply(
            multiply(rotation, operation.rotation), inverse_rotation
        )
        transformed_translation = _add(
            _matrix_vector(rotation, operation.translation),
            _subtract(
                translation,
                _matrix_vector(transformed_rotation, translation),
            ),
        )
        operations.append(AffineOperation(transformed_rotation, transformed_translation))
    rotated_lattice = multiply(source.lattice, transpose(rotation))
    if not same_lattice(rotated_lattice, target_lattice):
        raise ValueError("conjugating operation does not map onto target lattice")
    return AffineEmbedding(
        target_lattice,
        _deduplicate_operations(operations, target_lattice),
        source.hall_number,
    )


def affine_equivalence(
    source: AffineEmbedding,
    target: AffineEmbedding,
    parent: ParentAffineGroup,
    *,
    allow_parent_conjugacy: bool = True,
) -> str | None:
    """Return ``equal``/``parent_conjugate`` or ``None`` with exact proofs."""

    if same_lattice(source.lattice, target.lattice) and _operation_sets_equal(
        source.operations, target.operations, target.lattice
    ):
        return "equal"
    if not allow_parent_conjugacy:
        return None
    cosets = translation_cosets(parent.lattice, target.lattice)
    for parent_operation in parent.operations:
        rotated_lattice = multiply(source.lattice, transpose(parent_operation.rotation))
        if not same_lattice(rotated_lattice, target.lattice):
            continue
        for coset in cosets:
            translated = _add(parent_operation.translation, coset)
            conjugated = conjugate_embedding(
                source, parent_operation.rotation, translated, target.lattice
            )
            if _operation_sets_equal(
                conjugated.operations, target.operations, target.lattice
            ):
                return "parent_conjugate"
    return None


def _compose(left: AffineOperation, right: AffineOperation) -> AffineOperation:
    return AffineOperation(
        multiply(left.rotation, right.rotation),
        _add(left.translation, _matrix_vector(left.rotation, right.translation)),
    )


def _prepare_affine_quotient(
    parent: ParentAffineGroup,
    child_lattice: RationalMatrix,
    *,
    max_quotient_order: int,
) -> _AffineQuotientPreflight:
    """Validate and size ``N_G(T_s)/T_s`` without allocating its table."""
    lattice = rational_matrix(child_lattice)
    relative = multiply(lattice, inverse(parent.lattice))
    if not is_integral(relative):
        raise ValueError("child lattice is not a parent translation sublattice")
    index_value = abs(determinant(relative))
    if index_value.denominator != 1 or index_value < 1:
        raise ValueError("parent/child translation quotient has invalid index")
    translation_index = int(index_value)
    normalizer_operations = _deduplicate_operations(
        tuple(
            operation
            for operation in parent.operations
            if same_lattice(
                multiply(lattice, transpose(operation.rotation)),
                lattice,
            )
        ),
        parent.lattice,
    )
    quotient_order = len(normalizer_operations) * translation_index
    if max_quotient_order < 1:
        raise ValueError("max_quotient_order must be positive")
    if quotient_order > max_quotient_order:
        raise RuntimeError(
            "Method 3 affine quotient order "
            f"{quotient_order} exceeds the exact diagnostic boundary "
            f"{max_quotient_order}; no translation cosets or multiplication "
            "table were built and the result was not truncated"
        )
    return _AffineQuotientPreflight(
        lattice=lattice,
        normalizer_operations=normalizer_operations,
        translation_index=translation_index,
        quotient_order=quotient_order,
    )


def _materialize_affine_quotient(
    parent: ParentAffineGroup,
    prepared: _AffineQuotientPreflight,
) -> AffineQuotient:
    """Materialize the exact quotient after all configured size preflights."""
    cosets = translation_cosets(parent.lattice, prepared.lattice)
    lattice_inverse = inverse(prepared.lattice)
    operations: list[AffineOperation] = []
    operation_keys: set[tuple[RationalMatrix, RationalVector]] = set()
    for operation in prepared.normalizer_operations:
        for coset in cosets:
            candidate = AffineOperation(
                operation.rotation,
                _add(operation.translation, coset),
            )
            key = _operation_coset_key(candidate, lattice_inverse)
            if key not in operation_keys:
                operation_keys.add(key)
                operations.append(candidate)
    if len(operations) != prepared.quotient_order:
        raise ValueError(
            "finite affine quotient order disagrees with its exact lattice index"
        )

    operation_indices = {
        _operation_coset_key(operation, lattice_inverse): index
        for index, operation in enumerate(operations)
    }
    if len(operation_indices) != len(operations):
        raise ValueError("finite affine quotient contains duplicate operation cosets")

    def locate(operation: AffineOperation) -> int:
        key = _operation_coset_key(operation, lattice_inverse)
        try:
            return operation_indices[key]
        except KeyError as error:
            raise ValueError(
                "finite quotient is not closed under Seitz multiplication"
            ) from error

    products = tuple(
        tuple(locate(_compose(left, right)) for right in operations)
        for left in operations
    )
    identity = locate(
        AffineOperation(
            identity_matrix(),
            (Fraction(0), Fraction(0), Fraction(0)),
        )
    )
    return AffineQuotient(
        lattice=prepared.lattice,
        operations=tuple(operations),
        multiplication_table=products,
        identity_index=identity,
        translation_index=prepared.translation_index,
    )


def build_affine_quotient(
    parent: ParentAffineGroup,
    child_lattice: RationalMatrix,
    *,
    max_quotient_order: int = 512,
) -> AffineQuotient:
    """Build the exact finite normalizer quotient ``N_G(T_s)/T_s``.

    A selected child lattice need not be normal in the full parent space
    group, so ``G/T_s`` is generally undefined.  The normalizer is the unique
    parent subgroup in which ``T_s`` is normal.  This is also the finite group
    on which a fixed-translation-domain inverse-Landau calculation acts.
    """
    prepared = _prepare_affine_quotient(
        parent,
        child_lattice,
        max_quotient_order=max_quotient_order,
    )
    return _materialize_affine_quotient(parent, prepared)


def locate_embedding_subgroup(
    embedding: AffineEmbedding,
    quotient: AffineQuotient,
) -> tuple[int, ...]:
    """Locate and validate ``H/T_s`` inside an exact affine quotient.

    A Stage-A embedding is a section over its point group: exactly one
    operation coset occurs for each rotation.  Multiple lifts of one rotation
    would introduce an undeclared translation outside ``T_s`` and are
    rejected rather than silently changing the embedding's translation group.
    """
    if not same_lattice(embedding.lattice, quotient.lattice):
        raise ValueError("embedding and affine quotient use different lattices")
    indices = tuple(quotient.locate(operation) for operation in embedding.operations)
    subgroup = frozenset(indices)
    if len(subgroup) != len(indices):
        raise ValueError("embedding repeats an affine-operation coset")
    rotations = {operation.rotation for operation in embedding.operations}
    if len(rotations) != len(indices):
        raise ValueError(
            "embedding contains multiple operation cosets for one rotation; "
            "its declared translation lattice is incomplete"
        )
    if quotient.identity_index not in subgroup:
        raise ValueError("embedding subgroup does not contain the quotient identity")
    for left in subgroup:
        has_inverse = False
        for right in subgroup:
            product = quotient.multiplication_table[left][right]
            if product not in subgroup:
                raise ValueError("embedding operation cosets are not closed")
            if (
                product == quotient.identity_index
                and quotient.multiplication_table[right][left]
                == quotient.identity_index
            ):
                has_inverse = True
        if not has_inverse:
            raise ValueError("embedding operation coset has no subgroup inverse")
    return tuple(index for index in range(len(quotient.operations)) if index in subgroup)


def _rotation_subgroup_generated(
    generators: Sequence[RationalMatrix],
) -> frozenset[RationalMatrix]:
    """Return the exact finite subgroup generated by ``generators``.

    Crystallographic point operations form a finite group.  Consequently the
    positive-word monoid generated here already contains every inverse; no
    floating tolerance or separately guessed inverse operation is required.
    """

    identity = identity_matrix()
    result = {identity}
    queue: deque[RationalMatrix] = deque([identity])
    ordered_generators = tuple(generators)
    while queue:
        current = queue.popleft()
        for generator in ordered_generators:
            candidate = multiply(current, generator)
            if candidate not in result:
                result.add(candidate)
                queue.append(candidate)
    return frozenset(result)


def _point_group_generators(
    rotations: frozenset[RationalMatrix],
) -> tuple[RationalMatrix, ...]:
    """Choose a deterministic small generating set for an exact point group.

    At each step the operation producing the largest exact subgroup extension
    is selected.  Minimal cardinality is not needed for completeness, but this
    greedy construction gives the familiar one-to-three generators for the
    crystallographic point groups and avoids treating every non-identity
    rotation as an independent affine-lift variable.
    """

    identity = identity_matrix()
    if identity not in rotations:
        raise ValueError("point subgroup must contain the identity rotation")
    if any(
        multiply(left, right) not in rotations
        for left in rotations
        for right in rotations
    ):
        raise ValueError("point subgroup rotations are not closed under multiplication")

    def matrix_key(matrix: RationalMatrix) -> tuple[tuple[int, int], ...]:
        return tuple(
            (value.numerator, value.denominator)
            for row in matrix
            for value in row
        )

    generators: list[RationalMatrix] = []
    generated = frozenset({identity})
    while generated != rotations:
        candidates = sorted(rotations - generated, key=matrix_key)
        best = max(
            candidates,
            key=lambda candidate: len(
                _rotation_subgroup_generated((*generators, candidate))
            ),
        )
        generators.append(best)
        generated = _rotation_subgroup_generated(generators)
        if not generated.issubset(rotations):
            # The explicit closure check above should make this unreachable,
            # but keep the invariant local to the generator construction.
            raise ValueError("point subgroup generator escaped the supplied rotations")
    return tuple(generators)


def enumerate_lifted_point_subgroups(
    parent: ParentAffineGroup,
    child_lattice: RationalMatrix,
    point_subgroups: Sequence[Sequence[RationalMatrix]],
    *,
    max_quotient_order: int = 512,
    max_lift_assignments: int = 1_000_000,
    max_intermediate_subgroups: int = 16384,
) -> tuple[AffineEmbedding, ...]:
    """Enumerate every closed affine lift of the supplied point subgroups.

    ``point_subgroups`` must contain exact rotation subgroups of the parent
    point group.  For each point group ``H``, this function solves for every
    section of the finite extension

    ``1 -> T_parent/T_s -> pi^-1(H) -> H -> 1``.

    A section is uniquely determined by the lifts of any generating set of
    ``H``.  We therefore enumerate generator lifts, close them under the exact
    quotient multiplication table, and reject a branch immediately if the
    projection ceases to be injective.  This is the finite-group form of the
    Seitz factor-set/cocycle equations and is complete without the exponential
    Cartesian assignment over every rotation.  Projection injectivity excludes
    embeddings with a larger translation lattice than ``T_s``.
    """

    quotient_preflight = _prepare_affine_quotient(
        parent,
        child_lattice,
        max_quotient_order=max_quotient_order,
    )
    translation_index = quotient_preflight.translation_index
    normalizer_operations = quotient_preflight.normalizer_operations

    identity_rotation = identity_matrix()
    operations_per_rotation: dict[RationalMatrix, int] = {}
    for operation in normalizer_operations:
        operations_per_rotation[operation.rotation] = (
            operations_per_rotation.get(operation.rotation, 0) + 1
        )

    prepared_point_subgroups: list[
        tuple[frozenset[RationalMatrix], tuple[RationalMatrix, ...]]
    ] = []
    seen_point_subgroups: set[frozenset[RationalMatrix]] = set()
    total_lift_assignments = 0
    for raw_rotations in point_subgroups:
        rotations = frozenset(rational_matrix(rotation) for rotation in raw_rotations)
        if not rotations.issubset(operations_per_rotation):
            raise ValueError("point subgroup contains a rotation outside the parent group")
        generators = _point_group_generators(rotations)
        if rotations in seen_point_subgroups:
            continue
        seen_point_subgroups.add(rotations)

        # Only generator lifts are independent.  All remaining operation
        # lifts are fixed by quotient multiplication (equivalently, by the
        # factor-set equations).  This exact preflight retains the useful
        # no-allocation safety guarantee without the former |A|^(|H|-1)
        # overestimate, where A = T_parent/T_s.
        assignment_count = 1
        for rotation in generators:
            assignment_count *= operations_per_rotation[rotation] * translation_index
        total_lift_assignments += assignment_count
        if total_lift_assignments > max_lift_assignments:
            raise RuntimeError(
                "Method 3 affine-lift search requires "
                f"{total_lift_assignments} generator-lift assignments, exceeding "
                f"the exact diagnostic boundary {max_lift_assignments}; no "
                "translation cosets or lifts were built and the result was "
                "not truncated"
            )
        prepared_point_subgroups.append((rotations, generators))

    affine_quotient = _materialize_affine_quotient(parent, quotient_preflight)
    quotient = affine_quotient.operations
    products = affine_quotient.multiplication_table
    identity = affine_quotient.identity_index

    lifts_by_rotation: dict[RationalMatrix, tuple[int, ...]] = {
        rotation: tuple(
            index
            for index, operation in enumerate(quotient)
            if operation.rotation == rotation
        )
        for rotation in operations_per_rotation
    }

    def extend_section(
        subgroup: frozenset[int],
        lift: int,
        allowed_rotations: frozenset[RationalMatrix],
    ) -> frozenset[int] | None:
        """Close one generator lift, pruning on the section invariant.

        A repeated projected rotation represented by two distinct quotient
        elements is a non-trivial element of the translation kernel.  Such a
        branch can never become a section and is rejected as soon as it
        appears, before its full affine subgroup is materialised.
        """

        generators = (*sorted(subgroup), lift)
        result = {identity}
        by_rotation = {identity_rotation: identity}
        queue = deque([identity])
        while queue:
            current = queue.popleft()
            for generator in generators:
                candidate = products[current][generator]
                if candidate in result:
                    continue
                rotation = quotient[candidate].rotation
                if rotation not in allowed_rotations:
                    return None
                existing = by_rotation.get(rotation)
                if existing is not None and existing != candidate:
                    return None
                by_rotation[rotation] = candidate
                result.add(candidate)
                queue.append(candidate)
        return frozenset(result)

    embeddings: list[AffineEmbedding] = []
    seen_subgroups: set[frozenset[int]] = set()
    intermediate_count = 0
    for rotations, generators in prepared_point_subgroups:
        frontier = {frozenset({identity})}
        generated_rotations: tuple[RationalMatrix, ...] = ()
        for generator_rotation in generators:
            generated_rotations = (*generated_rotations, generator_rotation)
            expected_rotations = _rotation_subgroup_generated(generated_rotations)
            next_frontier: set[frozenset[int]] = set()
            for subgroup in frontier:
                for lift in lifts_by_rotation[generator_rotation]:
                    candidate = extend_section(subgroup, lift, rotations)
                    if candidate is None:
                        continue
                    candidate_rotations = {
                        quotient[index].rotation for index in candidate
                    }
                    if candidate_rotations != expected_rotations:
                        raise ValueError(
                            "affine quotient projection disagrees with point-group "
                            "generator closure"
                        )
                    if candidate not in next_frontier:
                        next_frontier.add(candidate)
                        intermediate_count += 1
                        if intermediate_count > max_intermediate_subgroups:
                            raise RuntimeError(
                                "Method 3 affine-lift search exceeded its exact "
                                f"intermediate-subgroup boundary "
                                f"({max_intermediate_subgroups}); the result was not truncated"
                            )
            frontier = next_frontier
            if not frontier:
                break

        for subgroup in frontier:
            subgroup_rotations = {quotient[index].rotation for index in subgroup}
            if (
                len(subgroup) != len(rotations)
                or subgroup_rotations != rotations
                or subgroup in seen_subgroups
            ):
                continue
            seen_subgroups.add(subgroup)
            embeddings.append(
                AffineEmbedding(
                    child_lattice,
                    tuple(quotient[index] for index in sorted(subgroup)),
                )
            )
    return tuple(embeddings)


def normalizer_point_subgroups(
    parent: ParentAffineGroup,
    child_lattice: RationalMatrix,
    target_space_group_number: int,
    parent_cartesian_lattice: Sequence[Sequence[float]],
    *,
    max_subgroups: int = 4096,
    symprec: float | None = None,
) -> tuple[tuple[RationalMatrix, ...], ...]:
    """Enumerate target-compatible point subgroups of ``N_G(T_s)``.

    Compatibility is checked in two independent ways: the exact rotations
    must form a subgroup of the lattice normalizer, and spglib must identify
    its crystallographic point-group symbol as that of the target Hall type.
    Spglib is used only for the conventional symbol classification; closure,
    membership and lattice stabilization are exact.
    """

    configured_symprec_angstrom = _positive_tolerance(
        get_config().affine_exact_cartesian_tolerance_angstrom
        if symprec is None
        else symprec,
        "exact affine Cartesian tolerance (Å)",
    )
    symprec_angstrom = _effective_exact_cartesian_tolerance(
        parent_cartesian_lattice, configured_symprec_angstrom
    )
    target_hall, _target_hall_symbol = _hall_for_number(target_space_group_number)
    target_type = spglib.get_spacegroup_type(target_hall)
    if target_type is None:
        raise ValueError(
            f"cannot determine target point group for SG {target_space_group_number}"
        )
    target_symbol = str(target_type.pointgroup_international).replace(" ", "")
    target_database = spglib.get_symmetry_from_database(target_hall)
    target_order = len(
        {
            rational_matrix(rotation.tolist())
            for rotation in target_database["rotations"]
        }
    )

    rotations: list[RationalMatrix] = []
    for operation in parent.operations:
        rotated_lattice = multiply(child_lattice, transpose(operation.rotation))
        if not same_lattice(rotated_lattice, child_lattice):
            continue
        if operation.rotation not in rotations:
            rotations.append(operation.rotation)
    identity = identity_matrix()
    if identity not in rotations:
        raise ValueError("lattice normalizer does not contain the identity")

    def rotation_index(rotation: RationalMatrix) -> int:
        try:
            return rotations.index(rotation)
        except ValueError as exc:
            raise ValueError("normalizer rotations are not closed") from exc

    products = [
        [rotation_index(multiply(left, right)) for right in rotations]
        for left in rotations
    ]
    identity_index = rotation_index(identity)

    def closure(seed: frozenset[int]) -> frozenset[int]:
        result = set(seed)
        result.add(identity_index)
        changed = True
        while changed:
            changed = False
            current = tuple(result)
            for left in current:
                for right in current:
                    product = products[left][right]
                    if product not in result:
                        result.add(product)
                        changed = True
        return frozenset(result)

    frontier = [frozenset({identity_index})]
    visited = set(frontier)
    complete: list[frozenset[int]] = []
    while frontier:
        subgroup = frontier.pop()
        if len(subgroup) == target_order:
            complete.append(subgroup)
            continue
        for element in range(len(rotations)):
            if element in subgroup:
                continue
            candidate = closure(subgroup | {element})
            if len(candidate) > target_order or candidate in visited:
                continue
            visited.add(candidate)
            if len(visited) > max_subgroups:
                raise RuntimeError(
                    "Method 3 point-subgroup enumeration exceeded its exact "
                    f"cost boundary ({max_subgroups} intermediate subgroups); "
                    "the result was not truncated"
                )
            frontier.append(candidate)

    lattice = np.asarray(parent_cartesian_lattice, dtype=float)
    result: list[tuple[RationalMatrix, ...]] = []
    for subgroup in complete:
        subgroup_rotations = tuple(rotations[index] for index in sorted(subgroup))
        identified = spglib.get_spacegroup_type_from_symmetry(
            np.asarray(
                [
                    [[int(value) for value in row] for row in rotation]
                    for rotation in subgroup_rotations
                ],
                dtype=int,
            ),
            np.zeros((len(subgroup_rotations), 3), dtype=float),
            lattice=lattice,
            symprec=symprec_angstrom,
        )
        if identified is None:
            continue
        symbol = str(identified.pointgroup_international).replace(" ", "")
        if symbol == target_symbol:
            result.append(subgroup_rotations)
    return tuple(result)


def identify_embedding_space_group(
    embedding: AffineEmbedding,
    parent_cartesian_lattice: Sequence[Sequence[float]],
    *,
    symprec: float | None = None,
) -> int | None:
    """Identify one exact embedding's international space-group type.

    Only this spglib classification boundary is floating point.  ``symprec``
    is an Å-valued Cartesian distance; the embedding itself remains exact.
    """

    configured_symprec_angstrom = _positive_tolerance(
        get_config().affine_exact_cartesian_tolerance_angstrom
        if symprec is None
        else symprec,
        "exact affine Cartesian tolerance (Å)",
    )

    parent_from_child = transpose(embedding.lattice)
    child_from_parent = inverse(parent_from_child)
    rotations: list[list[list[int]]] = []
    translations: list[list[float]] = []
    for operation in embedding.operations:
        child_rotation = multiply(
            multiply(child_from_parent, operation.rotation), parent_from_child
        )
        if not is_integral(child_rotation):
            raise ValueError("embedding rotation is not integral in child primitive axes")
        child_translation = _matrix_vector(child_from_parent, operation.translation)
        rotations.append(
            [[int(value) for value in row] for row in child_rotation]
        )
        translations.append([float(value % 1) for value in child_translation])
    child_cartesian = (
        np.asarray(embedding.lattice, dtype=float)
        @ np.asarray(parent_cartesian_lattice, dtype=float)
    )
    symprec_angstrom = _effective_exact_cartesian_tolerance(
        child_cartesian, configured_symprec_angstrom
    )
    identified = spglib.get_spacegroup_type_from_symmetry(
        np.asarray(rotations, dtype=int),
        np.asarray(translations, dtype=float),
        lattice=child_cartesian,
        symprec=symprec_angstrom,
    )
    return None if identified is None else int(identified.number)


def enumerate_target_affine_embeddings(
    parent: ParentAffineGroup,
    child_lattice: RationalMatrix,
    target_space_group_number: int,
    parent_cartesian_lattice: Sequence[Sequence[float]],
    *,
    max_point_subgroups: int = 4096,
    max_lifts: int = 16384,
    max_lift_assignments: int = 1_000_000,
    symprec: float | None = None,
) -> tuple[AffineEmbedding, ...]:
    """Enumerate exact target-SG embeddings for one selected translation lattice."""

    point_subgroups = normalizer_point_subgroups(
        parent,
        child_lattice,
        target_space_group_number,
        parent_cartesian_lattice,
        max_subgroups=max_point_subgroups,
        symprec=symprec,
    )
    lifts = enumerate_lifted_point_subgroups(
        parent,
        child_lattice,
        point_subgroups,
        max_lift_assignments=max_lift_assignments,
    )
    if len(lifts) > max_lifts:
        raise RuntimeError(
            "Method 3 affine-lift enumeration exceeded its exact cost boundary "
            f"({max_lifts} lifts); the result was not truncated"
        )
    typed = tuple(
        embedding
        for embedding in lifts
        if identify_embedding_space_group(
            embedding, parent_cartesian_lattice, symprec=symprec
        )
        == int(target_space_group_number)
    )
    # Preserve multiplicity of distinct operator subgroups.  Parent-conjugate
    # orientation/domain representatives are a proven relationship, but are
    # not equality and can be separate rows on the official first page.
    representatives: list[AffineEmbedding] = []
    for embedding in typed:
        if any(
            affine_equivalence(
                embedding,
                existing,
                parent,
                allow_parent_conjugacy=False,
            )
            == "equal"
            for existing in representatives
        ):
            continue
        representatives.append(embedding)
    return tuple(representatives)


def _verified_float_fraction(
    value: object, tolerance: float | None = None
) -> Fraction:
    number = float(value)
    tolerance_value = _positive_tolerance(
        get_config().fractional_coordinate_tolerance
        if tolerance is None
        else tolerance,
        "fractional-coordinate rational-recovery tolerance",
    )
    candidate = Fraction(number).limit_denominator(1_000_000)
    if abs(float(candidate) - number) > tolerance_value:
        raise ValueError(f"standardization value {number!r} is not rationally stable")
    return candidate


def standardized_embedding_identity(
    embedding: AffineEmbedding,
    parent: ParentAffineGroup,
    parent_cartesian_lattice: Sequence[Sequence[float]],
    target_space_group_number: int,
    *,
    symprec: float | None = None,
    fractional_tolerance: float | None = None,
) -> dict[str, object] | None:
    """Recover a conventional basis/origin and verify it exactly.

    Spglib supplies a candidate standard transformation from a two-orbit
    generic probe.  Its documented coordinate relation is
    ``x_std = P x_input + p`` and, for column lattice matrices,
    ``L_input = L_std P``.  This is converted to the project's row-basis
    convention.  Spglib commonly reports origin choice 1 even when
    ISODISTORT uses origin choice 2, so all Hall choices of the identified IT
    type are tried with exact reconstruction.  The successful source choice
    is then converted to the preferred Hall setting by an exact affine
    setting transform and reconstructed once more.  A failed bridge returns
    ``None``; no unverified setting is exposed to the UI.
    """

    configured_symprec_angstrom = _positive_tolerance(
        get_config().affine_exact_cartesian_tolerance_angstrom
        if symprec is None
        else symprec,
        "exact affine Cartesian tolerance (Å)",
    )
    fractional_tolerance_value = _positive_tolerance(
        get_config().fractional_coordinate_tolerance
        if fractional_tolerance is None
        else fractional_tolerance,
        "fractional-coordinate rational-recovery tolerance",
    )
    parent_from_primitive = transpose(embedding.lattice)
    primitive_from_parent = inverse(parent_from_primitive)
    rotations: list[list[list[int]]] = []
    translations: list[list[float]] = []
    exact_child_operations: list[tuple[RationalMatrix, RationalVector]] = []
    for operation in embedding.operations:
        child_rotation = multiply(
            multiply(primitive_from_parent, operation.rotation),
            parent_from_primitive,
        )
        if not is_integral(child_rotation):
            return None
        child_translation = _matrix_vector(
            primitive_from_parent, operation.translation
        )
        exact_child_operations.append((child_rotation, child_translation))
        rotations.append([[int(value) for value in row] for row in child_rotation])
        translations.append([float(value % 1) for value in child_translation])

    positions: list[list[float]] = []
    numbers: list[int] = []
    seeds = (
        ((0.1234567, 0.2345678, 0.3456789), 1),
        ((0.1732051, 0.2718281, 0.4142136), 2),
    )
    for seed_values, species in seeds:
        seed = tuple(Fraction(str(value)) for value in seed_values)
        orbit: list[RationalVector] = []
        for rotation, translation in exact_child_operations:
            point = tuple(
                value % 1
                for value in _add(_matrix_vector(rotation, seed), translation)
            )
            if point not in orbit:
                orbit.append(point)
        positions.extend([[float(value) for value in point] for point in orbit])
        numbers.extend([species] * len(orbit))

    child_cartesian = (
        np.asarray(embedding.lattice, dtype=float)
        @ np.asarray(parent_cartesian_lattice, dtype=float)
    )
    symprec_angstrom = _effective_exact_cartesian_tolerance(
        child_cartesian, configured_symprec_angstrom
    )
    dataset = spglib.get_symmetry_dataset(
        (child_cartesian, np.asarray(positions), numbers),
        symprec=symprec_angstrom,
    )
    if dataset is None or int(dataset.number) != int(target_space_group_number):
        return None
    try:
        transformation = rational_matrix(
            [
                [
                    _verified_float_fraction(value, fractional_tolerance_value)
                    for value in row
                ]
                for row in dataset.transformation_matrix
            ]
        )
        origin_shift: RationalVector = tuple(  # type: ignore[assignment]
            _verified_float_fraction(value, fractional_tolerance_value)
            for value in dataset.origin_shift
        )
        # Row-basis form of L_input = L_std P is
        # A_std = P^{-T} A_input.
        basis = multiply(inverse(transpose(transformation)), embedding.lattice)
        input_origin = _matrix_vector(
            inverse(transformation),
            tuple(-value for value in origin_shift),  # type: ignore[arg-type]
        )
        parent_origin = _matrix_vector(parent_from_primitive, input_origin)
    except (ArithmeticError, KeyError, TypeError, ValueError):
        return None

    space_group_number = int(dataset.number)
    reported_hall_number = int(dataset.hall_number)
    hall_numbers = _hall_numbers_for_space_group(space_group_number)
    if reported_hall_number not in hall_numbers:
        return None
    ordered_hall_numbers = (
        reported_hall_number,
        *(hall_number for hall_number in hall_numbers if hall_number != reported_hall_number),
    )
    source_identities: list[dict[str, object]] = []
    for hall_number in ordered_hall_numbers:
        hall_type = spglib.get_spacegroup_type(hall_number)
        if hall_type is None:
            continue
        source_identity: dict[str, object] = {
            "space_group_number": space_group_number,
            "symbol": str(hall_type.international_short),
            "basis": basis,
            "origin": parent_origin,
            "hall_number": hall_number,
        }
        try:
            reconstructed = embedding_from_identity(source_identity, parent)
        except (ArithmeticError, KeyError, TypeError, ValueError):
            continue
        if affine_equivalence(
            reconstructed,
            embedding,
            parent,
            allow_parent_conjugacy=False,
        ) == "equal":
            source_identities.append(source_identity)

    preferred_hall_number, _preferred_hall_symbol = _hall_for_number(
        space_group_number
    )
    if preferred_hall_number <= 0:
        return None
    preferred_type = spglib.get_spacegroup_type(preferred_hall_number)
    if preferred_type is None:
        return None

    for source_identity in source_identities:
        source_hall_number = int(source_identity["hall_number"])
        source_basis = rational_matrix(source_identity["basis"])  # type: ignore[arg-type]
        source_origin_values = tuple(
            as_fraction(value)
            for value in source_identity["origin"]  # type: ignore[union-attr]
        )
        if len(source_origin_values) != 3:
            continue
        source_origin: RationalVector = source_origin_values  # type: ignore[assignment]
        if source_hall_number == preferred_hall_number:
            setting_transforms = ((identity_matrix(), (Fraction(0),) * 3),)
        else:
            setting_transforms = _hall_setting_transforms(
                source_hall_number, preferred_hall_number
            )
        for source_from_preferred, source_origin_of_preferred in setting_transforms:
            preferred_basis = multiply(
                transpose(source_from_preferred), source_basis
            )
            preferred_origin = _add(
                source_origin,
                _matrix_vector(
                    transpose(source_basis), source_origin_of_preferred
                ),
            )
            preferred_identity: dict[str, object] = {
                "space_group_number": space_group_number,
                "symbol": str(preferred_type.international_short),
                "basis": preferred_basis,
                "origin": preferred_origin,
                "hall_number": preferred_hall_number,
            }
            try:
                reconstructed = embedding_from_identity(
                    preferred_identity, parent
                )
            except (ArithmeticError, KeyError, TypeError, ValueError):
                continue
            if affine_equivalence(
                reconstructed,
                embedding,
                parent,
                allow_parent_conjugacy=False,
            ) != "equal":
                continue
            # Callers and serialized UI identities intentionally retain the
            # historical schema: absence of hall_number means preferred Hall.
            preferred_identity.pop("hall_number")
            return preferred_identity
    return None
