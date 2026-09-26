"""Exact affine-group equivalence checks for Method 3 validation.

The website may print two bases/origins for the same subgroup embedding.  A
string comparison cannot decide whether they are merely different settings.
This module reconstructs the subgroup operations in parent fractional
coordinates and compares them modulo the *child primitive translation
lattice*.  It can additionally test conjugacy by every parent space-group
operation and every finite parent/child translation coset.

The implementation follows the crystallographic affine action

``(Q, q) (R, t) (Q, q)^-1 = (Q R Q^-1, Q t + q - Q R Q^-1 q)``

and uses :class:`fractions.Fraction` throughout.  Floating values occur only
at the spglib database boundary, where the tabulated crystallographic
fractions are reconstructed and checked against machine-roundoff error.
"""

from __future__ import annotations

import math
import re
import sys
from collections import deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from typing import Any

import spglib
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.io.isodistort_cif import _hall_for_number
from isocore.utils import get_config
from isocore.utils.lattice import (
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


@dataclass(frozen=True)
class AffineOperation:
    rotation: RationalMatrix
    translation: RationalVector


@dataclass(frozen=True)
class EmbeddingOperations:
    lattice: RationalMatrix
    operations: tuple[AffineOperation, ...]
    hall_number: int


@dataclass(frozen=True)
class ParentOperations:
    lattice: RationalMatrix
    operations: tuple[AffineOperation, ...]
    space_group_number: int
    hall_number: int
    symbol: str


def _field(identity: object, name: str) -> Any:
    if isinstance(identity, Mapping):
        return identity[name]
    return getattr(identity, name)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _vector_json(vector: RationalVector) -> list[str]:
    return [_fraction_text(value) for value in vector]


def _matrix_json(matrix: RationalMatrix) -> list[list[str]]:
    return [[_fraction_text(value) for value in row] for row in matrix]


def _database_fraction(value: object) -> Fraction:
    """Recover one rational spglib database value and verify roundoff only.

    spglib exposes its exact crystallographic tables through double precision.
    ``limit_denominator`` recovers the stored small rational; the subsequent
    128-ulp-scale check prevents an arbitrary measured float from silently
    being accepted as an exact symmetry translation.
    """

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


def _parent_translation_lattice(
    operations: Sequence[AffineOperation],
) -> RationalMatrix:
    """Recover the primitive parent translations in the uploaded axes.

    The HM centering letter belongs to the standardized conventional setting;
    it is not a statement about the basis of an uploaded primitive or
    nonstandard CIF.  In input fractional coordinates, integer translations
    generate ``Z^3`` and every identity-rotation Seitz operator supplies an
    additional translation coset.  If there are ``n`` such cosets, the
    primitive translation lattice has determinant ``1/n``.  We independently
    search a short exact rational basis containing all of those generators.
    """

    zero: RationalVector = (Fraction(0), Fraction(0), Fraction(0))
    identity = identity_matrix()
    cosets: list[RationalVector] = [zero]
    for operation in operations:
        if operation.rotation != identity:
            continue
        reduced: RationalVector = tuple(  # type: ignore[assignment]
            value % 1 for value in operation.translation
        )
        if reduced not in cosets:
            cosets.append(reduced)

    expected_determinant = Fraction(1, len(cosets))
    candidate_vectors: set[RationalVector] = set()
    for coset in cosets:
        for shift in product((-1, 0, 1), repeat=3):
            vector: RationalVector = tuple(  # type: ignore[assignment]
                coset[index] + shift[index] for index in range(3)
            )
            if vector != zero:
                candidate_vectors.add(vector)
    ordered = sorted(
        candidate_vectors,
        key=lambda vector: (
            sum(value * value for value in vector),
            tuple((value.numerator, value.denominator) for value in vector),
        ),
    )
    integer_generators: tuple[RationalVector, ...] = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    for rows in combinations(ordered, 3):
        matrix: RationalMatrix = rows
        if abs(determinant(matrix)) != expected_determinant:
            continue
        matrix_inverse = inverse(matrix)
        if all(
            all(value.denominator == 1 for value in _row_matrix(item, matrix_inverse))
            for item in (*integer_generators, *cosets[1:])
        ):
            return matrix
    raise ValueError("could not reconstruct the uploaded parent translation lattice")


def _vector_add(left: RationalVector, right: RationalVector) -> RationalVector:
    return tuple(left[i] + right[i] for i in range(3))  # type: ignore[return-value]


def _vector_subtract(left: RationalVector, right: RationalVector) -> RationalVector:
    return tuple(left[i] - right[i] for i in range(3))  # type: ignore[return-value]


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


def _vector_in_lattice(vector: RationalVector, lattice: RationalMatrix) -> bool:
    coordinates = _row_matrix(vector, inverse(lattice))
    return all(value.denominator == 1 for value in coordinates)


def _operation_equivalent(
    left: AffineOperation,
    right: AffineOperation,
    lattice: RationalMatrix,
) -> bool:
    return (
        left.rotation == right.rotation
        and _vector_in_lattice(
            _vector_subtract(left.translation, right.translation), lattice
        )
    )


def _deduplicate_operations(
    operations: Sequence[AffineOperation],
    lattice: RationalMatrix,
) -> tuple[AffineOperation, ...]:
    result: list[AffineOperation] = []
    for operation in operations:
        if not any(
            _operation_equivalent(operation, existing, lattice)
            for existing in result
        ):
            result.append(operation)
    return tuple(result)


def _operation_sets_equal(
    left: Sequence[AffineOperation],
    right: Sequence[AffineOperation],
    lattice: RationalMatrix,
) -> bool:
    return len(left) == len(right) and all(
        any(_operation_equivalent(operation, candidate, lattice) for candidate in right)
        for operation in left
    )


def parent_operations(
    structure: Structure,
    symprec: float,
    angle_tolerance_degrees: float | None = None,
) -> ParentOperations:
    """Return the detected parent affine group in the uploaded coordinates.

    ``symprec`` is an absolute Cartesian distance in Å.  The angle tolerance
    is likewise an explicit scientific input; callers that omit it use the
    configured pymatgen/spglib value rather than a library-version-dependent
    implicit default.
    """

    symprec_angstrom = float(symprec)
    angle_degrees = float(
        get_config().symmetry_angle_tolerance_degrees
        if angle_tolerance_degrees is None
        else angle_tolerance_degrees
    )
    if not math.isfinite(symprec_angstrom) or symprec_angstrom <= 0.0:
        raise ValueError("parent symmetry Cartesian tolerance must be finite and positive")
    if not math.isfinite(angle_degrees) or angle_degrees <= 0.0:
        raise ValueError("parent symmetry angle tolerance must be finite and positive")

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
            rotation=rational_matrix(rotation.tolist()),
            translation=tuple(_database_fraction(value) for value in translation),  # type: ignore[arg-type]
        )
        for rotation, translation in zip(
            dataset.rotations, dataset.translations, strict=True
        )
    ]
    lattice = _parent_translation_lattice(raw_operations)
    return ParentOperations(
        lattice=lattice,
        operations=_deduplicate_operations(raw_operations, lattice),
        space_group_number=int(dataset.number),
        hall_number=int(dataset.hall_number),
        symbol=symbol,
    )


def embedding_operations(
    identity: object,
    parent: ParentOperations,
) -> EmbeddingOperations:
    """Reconstruct one child affine subgroup in parent coordinates."""

    basis = rational_matrix(_field(identity, "basis"))
    origin_values = tuple(as_fraction(value) for value in _field(identity, "origin"))
    if len(origin_values) != 3:
        raise ValueError("embedding origin must contain three coordinates")
    origin: RationalVector = origin_values  # type: ignore[assignment]
    child_lattice = multiply(
        centering_primitive_matrix(_centering(_field(identity, "symbol"))),
        basis,
    )
    relative = multiply(child_lattice, inverse(parent.lattice))
    if not is_integral(relative):
        raise ValueError("embedding translation lattice is not a parent sublattice")

    space_group_number = int(_field(identity, "space_group_number"))
    hall_number, _hall_symbol = _hall_for_number(space_group_number)
    if hall_number <= 0:
        raise ValueError(f"no preferred Hall setting for space group {space_group_number}")
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
        child_translation: RationalVector = tuple(
            _database_fraction(value) for value in raw_translation
        )  # type: ignore[assignment]
        parent_translation = _vector_add(
            _matrix_vector(parent_from_child, child_translation),
            _vector_subtract(origin, _matrix_vector(parent_rotation, origin)),
        )
        operation = AffineOperation(parent_rotation, parent_translation)
        if not any(
            _operation_equivalent(operation, candidate, parent.lattice)
            for candidate in parent.operations
        ):
            raise ValueError(
                "reconstructed child operation is not contained in the parent space group"
            )
        transformed.append(operation)

    operations = _deduplicate_operations(transformed, child_lattice)
    if not operations:
        raise ValueError("embedding produced no affine operations")
    return EmbeddingOperations(
        lattice=child_lattice,
        operations=operations,
        hall_number=hall_number,
    )


def _translation_coset_representatives(
    parent_lattice: RationalMatrix,
    child_lattice: RationalMatrix,
) -> tuple[RationalVector, ...]:
    """Enumerate ``L_parent/L_child`` exactly by its three generators."""

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
    coordinate_representatives = [zero]
    queue: deque[RationalVector] = deque([zero])
    while queue and len(coordinate_representatives) < index:
        current = queue.popleft()
        for generator in generators:
            candidate = _vector_add(current, generator)
            if any(
                _vector_in_lattice(
                    _vector_subtract(candidate, existing), relative
                )
                for existing in coordinate_representatives
            ):
                continue
            coordinate_representatives.append(candidate)
            queue.append(candidate)
            if len(coordinate_representatives) == index:
                break
    if len(coordinate_representatives) != index:
        raise ValueError(
            "failed to enumerate the complete finite parent/child translation quotient"
        )
    return tuple(
        _row_matrix(representative, parent_lattice)
        for representative in coordinate_representatives
    )


def _conjugate_operations(
    operations: Sequence[AffineOperation],
    rotation: RationalMatrix,
    translation: RationalVector,
    target_lattice: RationalMatrix,
) -> tuple[AffineOperation, ...]:
    inverse_rotation = inverse(rotation)
    result: list[AffineOperation] = []
    for operation in operations:
        transformed_rotation = multiply(
            multiply(rotation, operation.rotation), inverse_rotation
        )
        transformed_translation = _vector_add(
            _matrix_vector(rotation, operation.translation),
            _vector_subtract(
                translation,
                _matrix_vector(transformed_rotation, translation),
            ),
        )
        result.append(
            AffineOperation(transformed_rotation, transformed_translation)
        )
    return _deduplicate_operations(result, target_lattice)


def affine_equivalence_witness(
    source: EmbeddingOperations,
    target: EmbeddingOperations,
    parent: ParentOperations,
) -> dict[str, Any] | None:
    """Return an exact equality/conjugacy witness, or ``None``."""

    if same_lattice(source.lattice, target.lattice) and _operation_sets_equal(
        source.operations, target.operations, target.lattice
    ):
        return {"relationship": "same_affine_subgroup"}

    cosets = _translation_coset_representatives(parent.lattice, target.lattice)
    for parent_operation in parent.operations:
        rotated_lattice = multiply(
            source.lattice, transpose(parent_operation.rotation)
        )
        if not same_lattice(rotated_lattice, target.lattice):
            continue
        for coset in cosets:
            effective_translation = _vector_add(
                parent_operation.translation, coset
            )
            conjugated = _conjugate_operations(
                source.operations,
                parent_operation.rotation,
                effective_translation,
                target.lattice,
            )
            if _operation_sets_equal(
                conjugated, target.operations, target.lattice
            ):
                return {
                    "relationship": "parent_affine_conjugate",
                    "parent_rotation": _matrix_json(parent_operation.rotation),
                    "parent_operation_translation": _vector_json(
                        parent_operation.translation
                    ),
                    "translation_coset": _vector_json(coset),
                    "effective_translation": _vector_json(
                        effective_translation
                    ),
                }
    return None


def compare_affine_embedding_sets(
    official: Sequence[object],
    local: Sequence[object],
    structure: Structure,
    *,
    symprec: float,
    angle_tolerance_degrees: float | None = None,
) -> dict[str, Any]:
    """Maximum-match official/local embeddings under parent affine conjugacy."""

    configured_angle_tolerance = float(
        get_config().symmetry_angle_tolerance_degrees
        if angle_tolerance_degrees is None
        else angle_tolerance_degrees
    )
    parent = parent_operations(
        structure,
        symprec,
        configured_angle_tolerance,
    )
    official_models: list[EmbeddingOperations | None] = []
    local_models: list[EmbeddingOperations | None] = []
    errors: list[dict[str, Any]] = []
    for side, identities, destination in (
        ("official", official, official_models),
        ("local", local, local_models),
    ):
        for index, identity in enumerate(identities):
            try:
                destination.append(embedding_operations(identity, parent))
            except (ArithmeticError, KeyError, TypeError, ValueError) as exc:
                destination.append(None)
                errors.append(
                    {
                        "side": side,
                        "index": index,
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                )

    edges: dict[int, list[tuple[int, dict[str, Any]]]] = {}
    for official_index, official_model in enumerate(official_models):
        if official_model is None:
            continue
        for local_index, local_model in enumerate(local_models):
            if local_model is None:
                continue
            witness = affine_equivalence_witness(
                local_model, official_model, parent
            )
            if witness is not None:
                edges.setdefault(official_index, []).append(
                    (local_index, witness)
                )
    for candidates in edges.values():
        candidates.sort(
            key=lambda item: (
                item[1]["relationship"] != "same_affine_subgroup",
                item[0],
            )
        )

    local_matches: dict[int, tuple[int, dict[str, Any]]] = {}

    def assign(official_index: int, visited: set[int]) -> bool:
        for local_index, witness in edges.get(official_index, []):
            if local_index in visited:
                continue
            visited.add(local_index)
            previous = local_matches.get(local_index)
            if previous is None or assign(previous[0], visited):
                local_matches[local_index] = (official_index, witness)
                return True
        return False

    for official_index in range(len(official)):
        assign(official_index, set())

    matches = [
        {
            "official_index": official_index,
            "local_index": local_index,
            **witness,
        }
        for local_index, (official_index, witness) in sorted(
            local_matches.items(), key=lambda item: item[1][0]
        )
    ]
    matched_official = {item["official_index"] for item in matches}
    matched_local = {item["local_index"] for item in matches}
    complete = (
        not errors
        and len(matches) == len(official)
        and len(matches) == len(local)
    )
    return {
        "status": "affine_equivalent" if complete else "incomplete",
        "comparison_semantics": (
            "Exact Fraction reconstruction in parent coordinates; affine operations "
            "are compared modulo the child primitive lattice and under parent "
            "space-group conjugacy, including every parent/child translation coset."
        ),
        "parent": {
            "space_group_number": parent.space_group_number,
            "symbol": parent.symbol,
            "hall_number": parent.hall_number,
            "symprec": float(symprec),
            "symmetry_cartesian_tolerance_angstrom": float(symprec),
            "symmetry_angle_tolerance_degrees": configured_angle_tolerance,
        },
        "match_count": len(matches),
        "official_count": len(official),
        "local_count": len(local),
        "complete_match": complete,
        "matches": matches,
        "unmatched_official_indices": [
            index for index in range(len(official)) if index not in matched_official
        ],
        "unmatched_local_indices": [
            index for index in range(len(local)) if index not in matched_local
        ],
        "verification_errors": errors,
    }
