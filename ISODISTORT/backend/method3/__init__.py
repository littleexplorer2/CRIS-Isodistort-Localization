"""Public API for Method 3 fixed-space and coupled-route analysis."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "CoupledRouteCandidate": ("backend.method3.coupled_routes", "CoupledRouteCandidate"),
    "CoupledRouteWitness": ("backend.method3.coupled_routes", "CoupledRouteWitness"),
    "lift_embedding_subgroup": ("backend.method3.coupled_routes", "lift_embedding_subgroup"),
    "resolve_coupled_route_witness": ("backend.method3.coupled_routes", "resolve_coupled_route_witness"),
    "ExactCharacterRepresentation": ("backend.method3.inverse_landau", "ExactCharacterRepresentation"),
    "FiniteGroup": ("backend.method3.inverse_landau", "FiniteGroup"),
    "FixedSpaceFeasibility": ("backend.method3.inverse_landau", "FixedSpaceFeasibility"),
    "FixedSpaceFeasibilityAnalyzer": ("backend.method3.inverse_landau", "FixedSpaceFeasibilityAnalyzer"),
    "FixedSpaceFeasibilitySummary": ("backend.method3.inverse_landau", "FixedSpaceFeasibilitySummary"),
    "RationalRepresentation": ("backend.method3.inverse_landau", "RationalRepresentation"),
    "analyze_fixed_space": ("backend.method3.inverse_landau", "analyze_fixed_space"),
    "character_representation": ("backend.method3.inverse_landau", "character_representation"),
    "direct_sum_character_representations": ("backend.method3.inverse_landau", "direct_sum_character_representations"),
    "direct_sum_representations": ("backend.method3.inverse_landau", "direct_sum_representations"),
    "fixed_space_basis": ("backend.method3.inverse_landau", "fixed_space_basis"),
    "fixed_space_projector": ("backend.method3.inverse_landau", "fixed_space_projector"),
    "homogeneous_strain_representation": ("backend.method3.inverse_landau", "homogeneous_strain_representation"),
    "pointwise_stabilizer": ("backend.method3.inverse_landau", "pointwise_stabilizer"),
    "polar_vector_representation": ("backend.method3.inverse_landau", "polar_vector_representation"),
    "remove_uniform_site_vectors": ("backend.method3.inverse_landau", "remove_uniform_site_vectors"),
    "site_displacement_representation": ("backend.method3.inverse_landau", "site_displacement_representation"),
    "site_scalar_representation": ("backend.method3.inverse_landau", "site_scalar_representation"),
    "site_vector_character_representation": ("backend.method3.inverse_landau", "site_vector_character_representation"),
    "site_vector_representation": ("backend.method3.inverse_landau", "site_vector_representation"),
    "symmetric_square_representation": ("backend.method3.inverse_landau", "symmetric_square_representation"),
    "DiagnosticCharacterRepresentationBundle": (
        "backend.method3.inverse_landau_adapter",
        "DiagnosticCharacterRepresentationBundle",
    ),
    "DiagnosticRepresentationBundle": ("backend.method3.inverse_landau_adapter", "DiagnosticRepresentationBundle"),
    "EmbeddingFeasibilityDiagnostic": ("backend.method3.inverse_landau_adapter", "EmbeddingFeasibilityDiagnostic"),
    "build_selected_character_representation": (
        "backend.method3.inverse_landau_adapter",
        "build_selected_character_representation",
    ),
    "build_selected_representation": ("backend.method3.inverse_landau_adapter", "build_selected_representation"),
    "diagnose_embedding_feasibility": ("backend.method3.inverse_landau_adapter", "diagnose_embedding_feasibility"),
    "diagnose_embeddings_feasibility": ("backend.method3.inverse_landau_adapter", "diagnose_embeddings_feasibility"),
    "parent_frame_fingerprint": ("backend.method3.inverse_landau_adapter", "parent_frame_fingerprint"),
    "quotient_finite_group": ("backend.method3.inverse_landau_adapter", "quotient_finite_group"),
    "stable_embedding_id": ("backend.method3.inverse_landau_adapter", "stable_embedding_id"),
    "Method3Query": ("backend.method3.search", "Method3Query"),
    "Method3ResultItem": ("backend.method3.search", "Method3ResultItem"),
    "Method3Solver": ("backend.method3.search", "Method3Solver"),
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
