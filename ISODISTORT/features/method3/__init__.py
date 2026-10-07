"""Public API for Method 3 fixed-space and coupled-route analysis."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "CoupledRouteCandidate": ("features.method3.coupled_routes", "CoupledRouteCandidate"),
    "CoupledRouteWitness": ("features.method3.coupled_routes", "CoupledRouteWitness"),
    "lift_embedding_subgroup": ("features.method3.coupled_routes", "lift_embedding_subgroup"),
    "resolve_coupled_route_witness": (
        "features.method3.coupled_routes",
        "resolve_coupled_route_witness",
    ),
    "ExactCharacterRepresentation": (
        "features.method3.inverse_landau",
        "ExactCharacterRepresentation",
    ),
    "FiniteGroup": ("features.method3.inverse_landau", "FiniteGroup"),
    "FixedSpaceFeasibility": ("features.method3.inverse_landau", "FixedSpaceFeasibility"),
    "FixedSpaceFeasibilityAnalyzer": (
        "features.method3.inverse_landau",
        "FixedSpaceFeasibilityAnalyzer",
    ),
    "FixedSpaceFeasibilitySummary": (
        "features.method3.inverse_landau",
        "FixedSpaceFeasibilitySummary",
    ),
    "RationalRepresentation": ("features.method3.inverse_landau", "RationalRepresentation"),
    "analyze_fixed_space": ("features.method3.inverse_landau", "analyze_fixed_space"),
    "character_representation": (
        "features.method3.inverse_landau",
        "character_representation",
    ),
    "direct_sum_character_representations": (
        "features.method3.inverse_landau",
        "direct_sum_character_representations",
    ),
    "direct_sum_representations": (
        "features.method3.inverse_landau",
        "direct_sum_representations",
    ),
    "fixed_space_basis": ("features.method3.inverse_landau", "fixed_space_basis"),
    "fixed_space_projector": ("features.method3.inverse_landau", "fixed_space_projector"),
    "homogeneous_strain_representation": (
        "features.method3.inverse_landau",
        "homogeneous_strain_representation",
    ),
    "pointwise_stabilizer": ("features.method3.inverse_landau", "pointwise_stabilizer"),
    "polar_vector_representation": (
        "features.method3.inverse_landau",
        "polar_vector_representation",
    ),
    "remove_uniform_site_vectors": (
        "features.method3.inverse_landau",
        "remove_uniform_site_vectors",
    ),
    "site_displacement_representation": (
        "features.method3.inverse_landau",
        "site_displacement_representation",
    ),
    "site_scalar_representation": (
        "features.method3.inverse_landau",
        "site_scalar_representation",
    ),
    "site_vector_character_representation": (
        "features.method3.inverse_landau",
        "site_vector_character_representation",
    ),
    "site_vector_representation": (
        "features.method3.inverse_landau",
        "site_vector_representation",
    ),
    "symmetric_square_representation": (
        "features.method3.inverse_landau",
        "symmetric_square_representation",
    ),
    "DiagnosticCharacterRepresentationBundle": (
        "features.method3.inverse_landau_adapter",
        "DiagnosticCharacterRepresentationBundle",
    ),
    "DiagnosticRepresentationBundle": (
        "features.method3.inverse_landau_adapter",
        "DiagnosticRepresentationBundle",
    ),
    "EmbeddingFeasibilityDiagnostic": (
        "features.method3.inverse_landau_adapter",
        "EmbeddingFeasibilityDiagnostic",
    ),
    "build_selected_character_representation": (
        "features.method3.inverse_landau_adapter",
        "build_selected_character_representation",
    ),
    "build_selected_representation": (
        "features.method3.inverse_landau_adapter",
        "build_selected_representation",
    ),
    "diagnose_embedding_feasibility": (
        "features.method3.inverse_landau_adapter",
        "diagnose_embedding_feasibility",
    ),
    "diagnose_embeddings_feasibility": (
        "features.method3.inverse_landau_adapter",
        "diagnose_embeddings_feasibility",
    ),
    "parent_frame_fingerprint": (
        "features.method3.inverse_landau_adapter",
        "parent_frame_fingerprint",
    ),
    "quotient_finite_group": (
        "features.method3.inverse_landau_adapter",
        "quotient_finite_group",
    ),
    "stable_embedding_id": (
        "features.method3.inverse_landau_adapter",
        "stable_embedding_id",
    ),
}

__all__ = sorted(_EXPORTS)  # noqa: PLE0605 - lazy export map is the source of truth


def __getattr__(name: str) -> Any:
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    value = getattr(import_module(module_name), attribute)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted((*globals(), *__all__))
