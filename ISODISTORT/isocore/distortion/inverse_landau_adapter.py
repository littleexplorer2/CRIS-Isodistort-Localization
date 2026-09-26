"""Diagnostic bridge from exact Method-3 embeddings to inverse Landau tests.

This module is intentionally not connected to the web, terminal, or default
Method-3 search.  It translates the exact affine quotient
``N_G(T_s)/T_s`` into finite-group representations and reports only
embedding-level symmetry feasibility.  It does **not** infer IR labels,
coupled routes, mode amplitudes, or selectable product rows.

The structure-coordinate boundary is the only floating-point step: finite
supercell sites are matched with the configured Cartesian Å tolerance.
Every resulting site permutation and every representation matrix must then
pass exact finite-group homomorphism checks before a fixed-space result is
returned.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import TypeAlias

import numpy as np
from pymatgen.core import Lattice, Structure

from ..utils.config_loader import get_config
from ..utils.lattice import (
    RationalMatrix,
    determinant,
    inverse,
    matrix_key,
    rational_matrix,
)
from .affine_embeddings import (
    AffineEmbedding,
    AffineOperation,
    AffineQuotient,
    ParentAffineGroup,
    build_affine_quotient,
    locate_embedding_subgroup,
    translation_cosets,
)
from .inverse_landau import (
    ExactCharacterRepresentation,
    FiniteGroup,
    FixedSpaceFeasibility,
    RationalRepresentation,
    analyze_fixed_space,
    character_representation,
    direct_sum_character_representations,
    direct_sum_representations,
    remove_uniform_site_vectors,
    site_vector_character_representation,
    site_vector_representation,
    symmetric_square_representation,
)

RationalVector: TypeAlias = tuple[Fraction, Fraction, Fraction]

_TYPE_ALIASES = {"displacement": "displacive", "order": "occupational"}
_SUPPORTED_TYPES = frozenset({"strain", "displacive"})
_UNSUPPORTED_REASONS = {
    "occupational": (
        "composition-constrained occupational subspaces are not implemented"
    ),
    "rotational": (
        "rigid-unit axial-vector coordinates and their site constraints are not implemented"
    ),
    "magnetic": (
        "axial time-odd and anti-unitary magnetic representations are not implemented"
    ),
}


@dataclass(frozen=True)
class DiagnosticRepresentationBundle:
    """Selected physical representation assembled for one affine quotient."""

    group: FiniteGroup
    selected_types: tuple[str, ...]
    representation: RationalRepresentation | None
    component_names: tuple[str, ...]
    unsupported_types: tuple[tuple[str, str], ...]
    site_count: int | None = None
    max_site_match_error_angstrom: float | None = None


@dataclass(frozen=True)
class DiagnosticCharacterRepresentationBundle:
    """Character-only selected representation for large read-only audits."""

    group: FiniteGroup
    selected_types: tuple[str, ...]
    representation: ExactCharacterRepresentation | None
    component_names: tuple[str, ...]
    unsupported_types: tuple[tuple[str, str], ...]
    site_count: int | None = None
    max_site_match_error_angstrom: float | None = None


@dataclass(frozen=True)
class EmbeddingFeasibilityDiagnostic:
    """Non-product diagnostic result for one Stage-A embedding."""

    embedding_id: str
    status: str
    route_status: str
    selectable: bool
    selected_types: tuple[str, ...]
    unsupported_types: tuple[tuple[str, str], ...]
    quotient_order: int
    subgroup_indices: tuple[int, ...]
    representation_dimension: int | None
    site_count: int | None
    max_site_match_error_angstrom: float | None
    fixed_space: FixedSpaceFeasibility | None

    @property
    def is_symmetry_reachable(self) -> bool | None:
        if self.fixed_space is None:
            return None
        return self.fixed_space.is_reachable


@dataclass(frozen=True)
class _SupercellSite:
    species: tuple[tuple[str, Fraction], ...]
    parent_fractional: tuple[float, float, float]


def quotient_finite_group(quotient: AffineQuotient) -> FiniteGroup:
    """Expose an affine quotient as the validated inverse-Landau group."""
    return FiniteGroup(
        elements=tuple(range(len(quotient.operations))),
        multiplication_table=quotient.multiplication_table,
        identity_index=quotient.identity_index,
        name="N_G(T_s)/T_s",
    )


def _extended_gcd(left: int, right: int) -> tuple[int, int, int]:
    old_remainder, remainder = left, right
    old_left, coefficient_left = 1, 0
    old_right, coefficient_right = 0, 1
    while remainder:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_left, coefficient_left = (
            coefficient_left,
            old_left - quotient * coefficient_left,
        )
        old_right, coefficient_right = (
            coefficient_right,
            old_right - quotient * coefficient_right,
        )
    if old_remainder < 0:
        return -old_remainder, -old_left, -old_right
    return old_remainder, old_left, old_right


def _row_hermite_normal_form(matrix: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    """Canonical row-HNF of a full-rank square integer matrix."""
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("row-HNF requires a non-empty square integer matrix")
    work = [list(row) for row in matrix]
    for column in range(size):
        pivot_source = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot_source is None:
            raise ValueError("row-HNF requires a full-rank matrix")
        work[column], work[pivot_source] = work[pivot_source], work[column]
        for row in range(column + 1, size):
            if work[row][column] == 0:
                continue
            pivot_row = work[column][:]
            other_row = work[row][:]
            pivot = pivot_row[column]
            value = other_row[column]
            divisor, pivot_factor, other_factor = _extended_gcd(pivot, value)
            work[column] = [
                pivot_factor * pivot_row[index]
                + other_factor * other_row[index]
                for index in range(size)
            ]
            work[row] = [
                (-value // divisor) * pivot_row[index]
                + (pivot // divisor) * other_row[index]
                for index in range(size)
            ]
        if work[column][column] < 0:
            work[column] = [-value for value in work[column]]
        pivot = work[column][column]
        for row in range(column):
            quotient = work[row][column] // pivot
            work[row] = [
                work[row][index] - quotient * work[column][index]
                for index in range(size)
            ]
    return tuple(tuple(row) for row in work)


def _canonical_lattice(lattice: RationalMatrix) -> RationalMatrix:
    parsed = rational_matrix(lattice)
    denominator = 1
    for row in parsed:
        for value in row:
            denominator = math.lcm(denominator, value.denominator)
    integer = tuple(
        tuple(int(value * denominator) for value in row) for row in parsed
    )
    normal = _row_hermite_normal_form(integer)
    common = denominator
    for row in normal:
        for value in row:
            common = gcd(common, abs(value))
    denominator //= common
    return tuple(
        tuple(Fraction(value // common, denominator) for value in row)
        for row in normal
    )


def _row_matrix(vector: RationalVector, matrix: RationalMatrix) -> RationalVector:
    return tuple(
        sum(vector[row] * matrix[row][column] for row in range(3))
        for column in range(3)
    )  # type: ignore[return-value]


def _canonical_translation(
    translation: RationalVector,
    canonical_lattice: RationalMatrix,
) -> RationalVector:
    coordinates = _row_matrix(translation, inverse(canonical_lattice))
    reduced = tuple(value % 1 for value in coordinates)
    return _row_matrix(reduced, canonical_lattice)  # type: ignore[arg-type]


def _fraction_payload(value: Fraction) -> tuple[int, int]:
    return value.numerator, value.denominator


def _matrix_payload(matrix: RationalMatrix) -> tuple[tuple[tuple[int, int], ...], ...]:
    return tuple(tuple(_fraction_payload(value) for value in row) for row in matrix)


def _operation_payload(
    operation: AffineOperation,
    canonical_lattice: RationalMatrix,
) -> tuple[
    tuple[tuple[tuple[int, int], ...], ...],
    tuple[tuple[int, int], ...],
]:
    translation = _canonical_translation(operation.translation, canonical_lattice)
    return matrix_key(operation.rotation), tuple(
        _fraction_payload(value) for value in translation
    )


def _operator_set_payload(
    operations: tuple[AffineOperation, ...],
    lattice: RationalMatrix,
) -> tuple[object, ...]:
    canonical_lattice = _canonical_lattice(lattice)
    return tuple(
        sorted(_operation_payload(operation, canonical_lattice) for operation in operations)
    )


def _sha256_payload(payload: object) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def parent_frame_fingerprint(parent: ParentAffineGroup) -> str:
    """Hash the exact parent group in its uploaded coordinate frame."""
    canonical_lattice = _canonical_lattice(parent.lattice)
    payload = {
        "schema": "method3-parent-frame-v1",
        "space_group_number": parent.space_group_number,
        "lattice": _matrix_payload(canonical_lattice),
        "operations": _operator_set_payload(parent.operations, canonical_lattice),
    }
    return _sha256_payload(payload)


def stable_embedding_id(
    parent: ParentAffineGroup,
    embedding: AffineEmbedding,
    *,
    quotient: AffineQuotient | None = None,
    max_quotient_order: int = 512,
) -> str:
    """Return a content-addressed ID for one exact operator subgroup.

    GL(3,Z) changes of the child-lattice basis and operation ordering are
    removed.  Parent conjugacy is deliberately *not* removed: distinct
    origin/orientation domain embeddings remain distinct IDs.
    """
    actual_quotient = quotient or build_affine_quotient(
        parent,
        embedding.lattice,
        max_quotient_order=max_quotient_order,
    )
    locate_embedding_subgroup(embedding, actual_quotient)
    canonical_lattice = _canonical_lattice(embedding.lattice)
    payload = {
        "schema": "method3-embedding-v1",
        "parent_frame": parent_frame_fingerprint(parent),
        "lattice": _matrix_payload(canonical_lattice),
        "operations": _operator_set_payload(embedding.operations, canonical_lattice),
    }
    return f"m3e1-{_sha256_payload(payload)}"


def _normalize_types(distortion_types: str | Sequence[str]) -> tuple[str, ...]:
    raw = (distortion_types,) if isinstance(distortion_types, str) else tuple(distortion_types)
    result: list[str] = []
    for item in raw:
        normalized = _TYPE_ALIASES.get(str(item).strip().lower(), str(item).strip().lower())
        if normalized not in _SUPPORTED_TYPES and normalized not in _UNSUPPORTED_REASONS:
            raise ValueError(f"unknown distortion type {item!r}")
        if normalized not in result:
            result.append(normalized)
    if not result:
        raise ValueError("at least one distortion type is required")
    return tuple(result)


def _positive_cartesian_tolerance(value: float | None) -> float:
    tolerance = float(
        get_config().symmetry_cartesian_tolerance_angstrom
        if value is None
        else value
    )
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("site matching tolerance must be finite and positive in Å")
    return tolerance


def _species_signature(site: object) -> tuple[tuple[str, Fraction], ...]:
    species = site.species
    return tuple(
        sorted(
            (str(specie), Fraction(str(float(occupancy))))
            for specie, occupancy in species.items()
        )
    )


def _periodic_distance(
    left_parent_fractional: np.ndarray,
    right_parent_fractional: np.ndarray,
    parent_to_child: np.ndarray,
    child_cartesian: Lattice,
) -> float:
    left_child = left_parent_fractional @ parent_to_child
    right_child = right_parent_fractional @ parent_to_child
    distance, _image = child_cartesian.get_distance_and_image(left_child, right_child)
    return float(distance)


def _supercell_sites_and_permutations(
    parent: ParentAffineGroup,
    quotient: AffineQuotient,
    structure: Structure,
    *,
    tolerance_angstrom: float,
) -> tuple[tuple[tuple[int, ...], ...], int, float]:
    child_lattice_float = np.asarray(quotient.lattice, dtype=float)
    parent_to_child = np.linalg.inv(child_lattice_float)
    child_cartesian = Lattice(child_lattice_float @ structure.lattice.matrix)
    cosets = translation_cosets(parent.lattice, quotient.lattice)
    sites: list[_SupercellSite] = []
    for site in structure:
        signature = _species_signature(site)
        parent_position = np.asarray(site.frac_coords, dtype=float)
        for coset in cosets:
            position = parent_position + np.asarray(coset, dtype=float)
            equivalent = [
                existing
                for existing in sites
                if existing.species == signature
                and _periodic_distance(
                    position,
                    np.asarray(existing.parent_fractional),
                    parent_to_child,
                    child_cartesian,
                )
                <= tolerance_angstrom
            ]
            if len(equivalent) > 1:
                raise ValueError("site recovery is ambiguous within the Cartesian tolerance")
            if not equivalent:
                sites.append(
                    _SupercellSite(
                        species=signature,
                        parent_fractional=tuple(float(value) for value in position),
                    )
                )

    expected_count = Fraction(len(structure)) * abs(determinant(quotient.lattice))
    if expected_count.denominator != 1 or len(sites) != int(expected_count):
        raise ValueError(
            "recovered supercell-site multiplicity disagrees with the exact "
            f"lattice-volume invariant: recovered={len(sites)}, expected={expected_count}"
        )

    permutations: list[tuple[int, ...]] = []
    maximum_error = 0.0
    parent_positions = np.asarray(
        [site.parent_fractional for site in sites],
        dtype=float,
    )
    child_positions = parent_positions @ parent_to_child
    species_signatures = tuple(site.species for site in sites)
    for operation in quotient.operations:
        rotation = np.asarray(operation.rotation, dtype=float)
        translation = np.asarray(operation.translation, dtype=float)
        transformed_parent = parent_positions @ rotation.T + translation
        transformed_child = transformed_parent @ parent_to_child
        distances = child_cartesian.get_all_distances(
            transformed_child,
            child_positions,
        )
        permutation: list[int] = []
        for source_index, source_species in enumerate(species_signatures):
            matches = [
                target_index
                for target_index, target_species in enumerate(species_signatures)
                if target_species == source_species
                and distances[source_index, target_index] <= tolerance_angstrom
            ]
            if len(matches) != 1:
                raise ValueError(
                    "affine operation does not give a unique species-preserving "
                    "supercell-site permutation within the Cartesian tolerance"
                )
            target_index = matches[0]
            error = float(distances[source_index, target_index])
            permutation.append(target_index)
            maximum_error = max(maximum_error, error)
        if len(set(permutation)) != len(sites):
            raise ValueError("affine operation site action is not bijective")
        permutations.append(tuple(permutation))
    return tuple(permutations), len(sites), maximum_error


def build_selected_representation(
    parent: ParentAffineGroup,
    quotient: AffineQuotient,
    structure: Structure,
    distortion_types: str | Sequence[str] = ("strain", "displacive"),
    *,
    site_tolerance_angstrom: float | None = None,
) -> DiagnosticRepresentationBundle:
    """Build supported physical-Type actions without assigning IR routes."""
    selected = _normalize_types(distortion_types)
    unsupported = tuple(
        (distortion_type, _UNSUPPORTED_REASONS[distortion_type])
        for distortion_type in selected
        if distortion_type in _UNSUPPORTED_REASONS
    )
    group = quotient_finite_group(quotient)
    vector = RationalRepresentation(
        group,
        tuple(operation.rotation for operation in quotient.operations),
        name="fractional polar vector",
    )
    components: list[RationalRepresentation] = []
    names: list[str] = []
    site_count: int | None = None
    maximum_error: float | None = None

    if "strain" in selected:
        components.append(
            symmetric_square_representation(vector, name="homogeneous strain")
        )
        names.append("strain")

    if "displacive" in selected:
        tolerance = _positive_cartesian_tolerance(site_tolerance_angstrom)
        permutations, site_count, maximum_error = _supercell_sites_and_permutations(
            parent,
            quotient,
            structure,
            tolerance_angstrom=tolerance,
        )
        if site_count >= 2:
            full = site_vector_representation(
                group,
                permutations,
                vector,
                name="finite-supercell atomic displacements",
            )
            components.append(
                remove_uniform_site_vectors(
                    full,
                    site_count=site_count,
                    vector_dimension=3,
                    name="internal finite-supercell atomic displacements",
                )
            )
            names.append("displacive")

    if not components:
        representation = None
    elif len(components) == 1:
        representation = components[0]
    else:
        representation = direct_sum_representations(
            components,
            name="selected physical distortion types",
        )
    return DiagnosticRepresentationBundle(
        group=group,
        selected_types=selected,
        representation=representation,
        component_names=tuple(names),
        unsupported_types=unsupported,
        site_count=site_count,
        max_site_match_error_angstrom=maximum_error,
    )


def build_selected_character_representation(
    parent: ParentAffineGroup,
    quotient: AffineQuotient,
    structure: Structure,
    distortion_types: str | Sequence[str] = ("strain", "displacive"),
    *,
    site_tolerance_angstrom: float | None = None,
) -> DiagnosticCharacterRepresentationBundle:
    """Build the same selected physical action as exact character data.

    The site permutation is recovered and verified exactly as in
    :func:`build_selected_representation`.  For displacements, the tensor
    product character is ``fixed_sites(g) * trace(R_g)`` and quotienting the
    uniform rigid translation subtracts ``trace(R_g)``.  This avoids dense
    ``3N x 3N`` matrices while preserving every fixed-space dimension and
    pointwise-stabilizer decision used by Stage-B feasibility audits.
    """

    selected = _normalize_types(distortion_types)
    unsupported = tuple(
        (distortion_type, _UNSUPPORTED_REASONS[distortion_type])
        for distortion_type in selected
        if distortion_type in _UNSUPPORTED_REASONS
    )
    group = quotient_finite_group(quotient)
    vector = RationalRepresentation(
        group,
        tuple(operation.rotation for operation in quotient.operations),
        name="fractional polar vector",
    )
    components: list[ExactCharacterRepresentation] = []
    names: list[str] = []
    site_count: int | None = None
    maximum_error: float | None = None

    if "strain" in selected:
        components.append(
            character_representation(
                symmetric_square_representation(vector, name="homogeneous strain")
            )
        )
        names.append("strain")

    if "displacive" in selected:
        tolerance = _positive_cartesian_tolerance(site_tolerance_angstrom)
        permutations, site_count, maximum_error = _supercell_sites_and_permutations(
            parent,
            quotient,
            structure,
            tolerance_angstrom=tolerance,
        )
        if site_count >= 2:
            components.append(
                site_vector_character_representation(
                    group,
                    permutations,
                    vector,
                    remove_uniform_vectors=True,
                    name="internal finite-supercell atomic displacements",
                )
            )
            names.append("displacive")

    if not components:
        representation = None
    elif len(components) == 1:
        representation = components[0]
    else:
        representation = direct_sum_character_representations(
            components,
            name="selected physical distortion types",
        )
    return DiagnosticCharacterRepresentationBundle(
        group=group,
        selected_types=selected,
        representation=representation,
        component_names=tuple(names),
        unsupported_types=unsupported,
        site_count=site_count,
        max_site_match_error_angstrom=maximum_error,
    )


def diagnose_embedding_feasibility(
    parent: ParentAffineGroup,
    embedding: AffineEmbedding,
    structure: Structure,
    distortion_types: str | Sequence[str] = ("strain", "displacive"),
    *,
    max_quotient_order: int = 512,
    site_tolerance_angstrom: float | None = None,
) -> EmbeddingFeasibilityDiagnostic:
    """Run an embedding-level Stage-B diagnostic, never a selectable route."""
    quotient = build_affine_quotient(
        parent,
        embedding.lattice,
        max_quotient_order=max_quotient_order,
    )
    subgroup = locate_embedding_subgroup(embedding, quotient)
    identifier = stable_embedding_id(parent, embedding, quotient=quotient)
    bundle = build_selected_representation(
        parent,
        quotient,
        structure,
        distortion_types,
        site_tolerance_angstrom=site_tolerance_angstrom,
    )
    if bundle.unsupported_types:
        status = "unsupported_distortion_types"
        route_status = "not_analyzed_unsupported_types"
        analysis = None
    elif bundle.representation is None:
        status = "no_physical_degrees_of_freedom"
        route_status = "not_analyzed_no_physical_degrees_of_freedom"
        analysis = None
    else:
        analysis = analyze_fixed_space(bundle.representation, subgroup)
        status = "embedding_feasible" if analysis.is_reachable else "embedding_infeasible"
        route_status = "unresolved_irrep_decomposition"
    return EmbeddingFeasibilityDiagnostic(
        embedding_id=identifier,
        status=status,
        route_status=route_status,
        selectable=False,
        selected_types=bundle.selected_types,
        unsupported_types=bundle.unsupported_types,
        quotient_order=len(quotient.operations),
        subgroup_indices=subgroup,
        representation_dimension=(
            None if bundle.representation is None else bundle.representation.dimension
        ),
        site_count=bundle.site_count,
        max_site_match_error_angstrom=bundle.max_site_match_error_angstrom,
        fixed_space=analysis,
    )


__all__ = [
    "DiagnosticCharacterRepresentationBundle",
    "DiagnosticRepresentationBundle",
    "EmbeddingFeasibilityDiagnostic",
    "build_selected_character_representation",
    "build_selected_representation",
    "diagnose_embedding_feasibility",
    "parent_frame_fingerprint",
    "quotient_finite_group",
    "stable_embedding_id",
]
