"""Exact coupled-order-parameter witnesses for Method 3 embeddings.

A single-IR isotropy subgroup is a stabilizer ``K`` of one order parameter.
An embedding that requires coupled order parameters has stabilizer equal to an
intersection ``H = K1 intersection K2 ...``.  The subgroups may have different
translation lattices, so every stabilizer is first lifted into the same exact
finite quotient ``N_G(T_H) / T_H`` before intersections are compared.

This module deliberately does not infer activity from space-group labels.
Activity in the user-selected physical representation is an independent
fixed-space test.  The witness here proves only that an exact affine embedding
can be obtained as an intersection of enumerated single-IR stabilizers.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

from ..backend.iso_wrapper import SubgroupInfo
from .affine_embeddings import (
    AffineEmbedding,
    AffineOperation,
    AffineQuotient,
    ParentAffineGroup,
    build_affine_quotient,
    locate_embedding_subgroup,
    translation_cosets,
)


@dataclass(frozen=True)
class CoupledRouteCandidate:
    """One enumerated single-IR route and its exact affine stabilizer."""

    route: SubgroupInfo
    embedding: AffineEmbedding


@dataclass(frozen=True)
class CoupledRouteWitness:
    """A minimal exact stabilizer-intersection witness for one embedding."""

    components: tuple[SubgroupInfo, ...]
    target_subgroup_order: int
    quotient_order: int
    explored_intersection_count: int
    verification_method: str = "exact_affine_stabilizer_intersection"


def _add_translation(
    left: tuple[Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def lift_embedding_subgroup(
    embedding: AffineEmbedding,
    quotient: AffineQuotient,
) -> frozenset[int]:
    """Return ``K/T_H`` for an embedding whose translations contain ``T_H``.

    ``AffineEmbedding.operations`` contains one coset representative for each
    point operation modulo its own translation lattice ``T_K``.  When
    ``T_H`` is a proper sublattice, all cosets in ``T_K/T_H`` must be added;
    otherwise intersections from different propagation vectors are wrong.
    """

    translations = translation_cosets(embedding.lattice, quotient.lattice)
    indices = frozenset(
        quotient.locate(
            AffineOperation(
                operation.rotation,
                _add_translation(operation.translation, translation),
            )
        )
        for operation in embedding.operations
        for translation in translations
    )
    if quotient.identity_index not in indices:
        raise ValueError("lifted embedding does not contain the quotient identity")
    table = quotient.multiplication_table
    if any(table[left][right] not in indices for left in indices for right in indices):
        raise ValueError("lifted embedding operation cosets are not closed")
    return indices


def _route_key(route: SubgroupInfo) -> tuple[object, ...]:
    """Deterministic ordering only; scientific equality uses affine groups."""

    return (
        str(route.k_point_label or ""),
        tuple(str(value) for value in (route.k_parameters or ())),
        str(route.irrep_label or ""),
        str(route.opd_symbol or ""),
        str(route.opd_dir_raw or ""),
        int(route.space_group_number),
        tuple(tuple(str(value) for value in row) for row in route.basis_vectors),
        tuple(str(value) for value in route.origin),
    )


def resolve_coupled_route_witness(
    parent: ParentAffineGroup,
    target: AffineEmbedding,
    candidates: Sequence[CoupledRouteCandidate],
    *,
    max_quotient_order: int = 512,
    max_intersection_states: int = 100_000,
) -> CoupledRouteWitness | None:
    """Find a minimum-cardinality exact single-stabilizer intersection.

    The dynamic program retains the shortest deterministic component tuple for
    each distinct exact intersection.  A positive ``max_intersection_states``
    is a fail-loud resource boundary, never a silent truncation.  Zero means
    unlimited.
    """

    if max_intersection_states < 0:
        raise ValueError("max_intersection_states must be non-negative")
    quotient = build_affine_quotient(
        parent,
        target.lattice,
        max_quotient_order=max_quotient_order,
    )
    target_indices = frozenset(locate_embedding_subgroup(target, quotient))
    full_group = frozenset(range(len(quotient.operations)))

    compatible: list[tuple[SubgroupInfo, frozenset[int]]] = []
    seen: set[tuple[tuple[object, ...], frozenset[int]]] = set()
    for candidate in candidates:
        try:
            stabilizer = lift_embedding_subgroup(candidate.embedding, quotient)
        except ValueError:
            continue
        # Equal stabilizers are single-IR routes, not coupled witnesses.  The
        # full group is a trivial component and cannot reduce an intersection.
        if not target_indices < stabilizer or stabilizer == full_group:
            continue
        key = (_route_key(candidate.route), stabilizer)
        if key in seen:
            continue
        seen.add(key)
        compatible.append((candidate.route, stabilizer))
    compatible.sort(key=lambda item: _route_key(item[0]))

    states: dict[frozenset[int], tuple[SubgroupInfo, ...]] = {full_group: ()}
    for route, stabilizer in compatible:
        updated = dict(states)
        for intersection, components in states.items():
            narrowed = intersection & stabilizer
            if narrowed == intersection or not target_indices <= narrowed:
                continue
            proposed = (*components, route)
            existing = updated.get(narrowed)
            if existing is None or len(proposed) < len(existing):
                updated[narrowed] = proposed
        if max_intersection_states and len(updated) > max_intersection_states:
            raise RuntimeError(
                "Method 3 exact coupled-route search requires more than "
                f"{max_intersection_states} distinct stabilizer intersections; "
                "increase runtime.method3_max_coupled_states or set it to 0. "
                "No candidate result was truncated."
            )
        states = updated

    components = states.get(target_indices)
    if not components or len(components) < 2:
        return None
    return CoupledRouteWitness(
        components=components,
        target_subgroup_order=len(target_indices),
        quotient_order=len(quotient.operations),
        explored_intersection_count=len(states),
    )


__all__ = [
    "CoupledRouteCandidate",
    "CoupledRouteWitness",
    "lift_embedding_subgroup",
    "resolve_coupled_route_witness",
]
