"""Immutable session result and failure contracts."""
from __future__ import annotations
from dataclasses import dataclass
from backend.export import SubgroupExportSpec
from backend.models.iso_mode_models import ModeIdentityFailure
from backend.utils import IsodistortError


@dataclass(frozen=True)
class DisplacementColumn:
    """One raw column in the same row order and fractional axes as its frame."""

    key: str
    label: str
    fractional_vectors: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class DisplacementSnapshot:
    """Versioned, immutable public contract for displacement consumers.

    Row lattice vectors are in angstroms; fractional rows obey r = x @ A.
    This frame is the basis supercell BEFORE shifting to subgroup origin.
    Raw columns retain their calculated scale/sign/order. Consumers wanting
    Cartesian normalization use the existing cart_normalized_mode_matrix.
    """

    schema: str
    candidate_identity_json: str
    calculation_context: tuple
    nmod: int
    parent_space_group_number: int
    parent_space_group_symbol: str
    child_space_group_symbol: str
    lattice_angstrom: tuple[tuple[float, float, float], ...]
    fractional_coordinates: tuple[tuple[float, float, float], ...]
    species: tuple[str, ...]
    columns: tuple[DisplacementColumn, ...]
    coordinate_convention: str = "row_fractional_basis_supercell_before_origin_shift"
    normalization: str = "raw_calculated_supercell_columns"

@dataclass(frozen=True)
class ExportCandidateFailure:
    """One candidate that failed calculation, preparation, or writer preflight."""

    position: int
    candidate_index: int
    subgroup_index: int
    parent_space_group_number: int
    k_point_label: str
    k_coordinates: tuple[str, ...]
    k_parameters: tuple[str, ...]
    irrep_label: str
    opd_symbol: str
    space_group_number: int
    basis_vectors: tuple[tuple[str, ...], ...]
    origin: tuple[str, ...]
    embedding_id: str
    error_type: str
    message: str


class UnresolvedModeIdentityError(IsodistortError):
    """A candidate has numerical columns but cannot enter authoritative writers."""

    def __init__(self, failures: list[ModeIdentityFailure]) -> None:
        if not failures:
            raise ValueError("unresolved mode identity error requires failures")
        self.failures = tuple(failures)
        details = "; ".join(
            (
                f"mode={failure.amplitude_key} orbit={failure.orbit_id or '<missing>'} "
                f"IR={failure.global_irrep or '<missing>'} "
                f"k=({','.join(failure.k_coordinates)}) "
                f"criterion={failure.reason}"
            )
            for failure in self.failures
        )
        super().__init__(
            "displacement modes lack verified ISO microscopic identities: " + details
        )


class ExportCancelledError(IsodistortError):
    """A batch stopped at a candidate boundary after a user cancellation."""


@dataclass(frozen=True)
class ExportBatchPlan:
    """Prepared candidates plus explicitly ineligible identity failures."""

    specs: tuple[SubgroupExportSpec, ...]
    skipped: tuple[ExportCandidateFailure, ...]


class ExportPreparationError(IsodistortError):
    """Batch export failed before publication because candidates were not ready."""

    def __init__(
        self,
        method: int,
        failures: list[ExportCandidateFailure],
    ) -> None:
        self.method = int(method)
        self.failures = tuple(failures)
        details = "; ".join(
            (
                f"#{failure.position} candidate_index={failure.candidate_index} "
                f"subgroup_index={failure.subgroup_index} "
                f"parent_SG={failure.parent_space_group_number} "
                f"{failure.k_point_label}{failure.k_coordinates}/"
                f"{failure.k_parameters}/{failure.irrep_label}/"
                f"{failure.opd_symbol}/SG{failure.space_group_number} "
                f"basis={failure.basis_vectors} origin={failure.origin} "
                f"embedding={failure.embedding_id}: "
                f"{failure.error_type}: {failure.message}"
            )
            for failure in self.failures
        )
        super().__init__(
            f"Method {self.method} batch export preparation failed for "
            f"{len(self.failures)} candidate(s): {details}"
        )
