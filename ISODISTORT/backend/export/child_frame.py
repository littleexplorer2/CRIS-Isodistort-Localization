"""Exact Seitz transformation into the emitted conventional child frame."""
from __future__ import annotations
import spglib
from backend.utils.lattice import (as_fraction, centering_primitive_matrix, identity_matrix, inverse, multiply, rational_matrix, transpose)
from backend.wrappers import SubgroupInfo
from backend.utils.affine_embeddings import translation_cosets
from .cif_displacive_model import ExactSeitzOperation

def exact_child_operations(
    subgroup: SubgroupInfo,
    affine_child,
) -> tuple[ExactSeitzOperation, ...]:
    """Conjugate the proven affine subgroup into the emitted child frame.

    ``AffineEmbedding.operations`` live in parent fractional axes and are
    represented modulo the child's primitive translation lattice.  The
    CIF model instead needs every conventional-cell Seitz operation in the
    exact frame ``x_parent = B^T x_child + q``.  Recover that frame first,
    then expand the conventional centering cosets.  The Hall database is
    used only as an independent exact-setting cross-check.
    """

    if affine_child.hall_number is None:
        raise ValueError("subgroup embedding has no exact Hall setting")
    basis = rational_matrix(subgroup.basis_vectors or identity_matrix())
    parent_from_child = transpose(basis)
    child_from_parent = inverse(parent_from_child)
    origin_values = tuple(
        as_fraction(value) for value in (subgroup.origin or (0, 0, 0))
    )
    if len(origin_values) != 3:
        raise ValueError("subgroup origin must contain three exact coordinates")
    origin = origin_values

    def matrix_vector(matrix, vector):
        return tuple(
            sum(matrix[row][column] * vector[column] for column in range(3))
            for row in range(3)
        )

    def add(left, right):
        return tuple(left[index] + right[index] for index in range(3))

    def subtract(left, right):
        return tuple(left[index] - right[index] for index in range(3))

    hall_type = spglib.get_spacegroup_type(int(affine_child.hall_number))
    if hall_type is None:
        raise ValueError(
            f"no space-group type for Hall number {affine_child.hall_number}"
        )
    hall_symbol = str(hall_type.hall_symbol or "")
    centering = next(
        (character for character in hall_symbol.upper() if character in "PABCIFR"),
        None,
    )
    if centering is None:
        raise ValueError("target Hall symbol has no conventional centering")
    centering_translations = translation_cosets(
        centering_primitive_matrix(centering),
        identity_matrix(),
    )

    operations: list[ExactSeitzOperation] = []
    for operation in affine_child.operations:
        child_rotation = multiply(
            multiply(child_from_parent, operation.rotation),
            parent_from_child,
        )
        parent_translation_about_origin = subtract(
            add(matrix_vector(operation.rotation, origin), operation.translation),
            origin,
        )
        child_translation = matrix_vector(
            child_from_parent,
            parent_translation_about_origin,
        )
        for centering_translation in centering_translations:
            operations.append(
                ExactSeitzOperation.from_values(
                    child_rotation,
                    add(child_translation, centering_translation),
                )
            )

    exact_by_key = {operation.key: operation for operation in operations}
    if len(exact_by_key) != len(operations):
        raise ValueError(
            "affine subgroup expansion produced duplicate conventional operations"
        )
    database = spglib.get_symmetry_from_database(int(affine_child.hall_number))
    if database is None:
        raise ValueError(
            f"no symmetry database entry for Hall number {affine_child.hall_number}"
        )
    database_operations = tuple(
        ExactSeitzOperation.from_values(
            rotation.tolist(),
            tuple(as_fraction(float(value)) for value in translation),
        )
        for rotation, translation in zip(
            database["rotations"],
            database["translations"],
            strict=True,
        )
    )
    if set(exact_by_key) != {operation.key for operation in database_operations}:
        raise ValueError(
            "conjugated affine subgroup disagrees with its target Hall setting"
        )
    return tuple(exact_by_key[key] for key in sorted(exact_by_key))
