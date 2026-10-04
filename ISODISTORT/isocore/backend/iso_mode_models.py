"""Exact models for ISO microscopic displacement and site-irrep data.

The interactive ISO program reports physically irreducible representations
in its standard conventional setting.  This module deliberately keeps every
reported rational number as :class:`fractions.Fraction`; conversion to a
floating Cartesian displacement happens only at the structure boundary.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Literal, TypeAlias

import numpy as np

ExactVector: TypeAlias = tuple[Fraction, Fraction, Fraction]
ExactMatrix3: TypeAlias = tuple[ExactVector, ExactVector, ExactVector]
PrintedCoefficient: TypeAlias = Fraction | Decimal


@dataclass(frozen=True)
class InvariantDirection:
    """One parent irrep with a nonzero fixed space for an exact child embedding."""

    irrep_label: str
    direction_raw: str
    subgroup_number: int
    subgroup_symbol: str
    size: int


TensorCoefficientQuality: TypeAlias = Literal[
    "exact_rational",
    "printed_decimal",
    "unresolved",
]


@dataclass(frozen=True)
class SymbolicTensorComponent:
    """One printed tensor basis function resolved against a named basis.

    Decimal coefficients retain exactly the finite precision printed by ISO;
    they are not promoted to an algebraic number such as ``sqrt(3)``.  A
    caller may use them for a tolerance-based whole-space comparison, but may
    only treat :attr:`exact_coefficients` as exact arithmetic.
    """

    expression_raw: str
    coefficient_basis: tuple[str, ...]
    coefficients: tuple[PrintedCoefficient, ...] | None
    quality: TensorCoefficientQuality
    reason: str | None = None

    def __post_init__(self) -> None:
        if self.coefficients is not None and (
            len(self.coefficients) != len(self.coefficient_basis)
        ):
            raise ValueError("tensor coefficient vector has the wrong dimension")
        if self.quality == "unresolved" and self.coefficients is not None:
            raise ValueError("an unresolved tensor component cannot expose coefficients")
        if self.quality != "unresolved" and self.coefficients is None:
            raise ValueError("a resolved tensor component requires coefficients")
        if self.quality == "exact_rational" and any(
            not isinstance(value, Fraction) for value in (self.coefficients or ())
        ):
            raise ValueError("exact tensor coefficients must be Fraction values")
        if self.quality == "printed_decimal" and not any(
            isinstance(value, Decimal) for value in (self.coefficients or ())
        ):
            raise ValueError("printed-decimal tensor coefficients require a Decimal")

    @property
    def exact_coefficients(self) -> tuple[Fraction, ...] | None:
        if self.quality != "exact_rational" or self.coefficients is None:
            return None
        return tuple(self.coefficients)  # type: ignore[return-value]

    def as_float_vector(self) -> np.ndarray | None:
        if self.coefficients is None:
            return None
        return np.asarray([float(value) for value in self.coefficients], dtype=float)


@dataclass(frozen=True)
class MacroscopicTensorBlock:
    """One ISO macroscopic-irrep copy and its ordered components."""

    global_irrep: str
    components: tuple[SymbolicTensorComponent, ...]
    source_order: int
    copy_index: int | None = None

    def __post_init__(self) -> None:
        if not self.global_irrep:
            raise ValueError("macroscopic tensor block requires an irrep label")
        if not self.components:
            raise ValueError("macroscopic tensor block requires at least one component")
        bases = {component.coefficient_basis for component in self.components}
        if len(bases) != 1:
            raise ValueError("macroscopic tensor components use inconsistent bases")


MacroscopicTensorStatus: TypeAlias = Literal["verified", "partial", "unresolved"]


@dataclass(frozen=True)
class MacroscopicTensorBasis:
    """Auditable result of a ``SHOW MACROSCOPIC`` distortion query."""

    rank_signature: str
    coefficient_basis: tuple[str, ...]
    blocks: tuple[MacroscopicTensorBlock, ...]
    status: MacroscopicTensorStatus
    reason: str | None = None
    invariant_directions: tuple[InvariantDirection, ...] = ()

    def __post_init__(self) -> None:
        if self.status == "unresolved" and self.blocks:
            raise ValueError("an unresolved tensor basis cannot expose trusted blocks")
        if self.status != "unresolved" and not self.blocks:
            raise ValueError("a resolved tensor basis requires blocks")
        if any(
            component.coefficient_basis != self.coefficient_basis
            for block in self.blocks
            for component in block.components
        ):
            raise ValueError("tensor basis metadata does not match its components")

    @classmethod
    def unresolved(
        cls,
        *,
        rank_signature: str,
        coefficient_basis: tuple[str, ...],
        reason: str,
        invariant_directions: tuple[InvariantDirection, ...] = (),
    ) -> MacroscopicTensorBasis:
        return cls(
            rank_signature=rank_signature,
            coefficient_basis=coefficient_basis,
            blocks=(),
            status="unresolved",
            reason=reason,
            invariant_directions=invariant_directions,
        )


@dataclass(frozen=True)
class MicroscopicVectorRow:
    """One representative point and all ISO standard displacement columns."""

    point_raw: tuple[str, str, str]
    displacements: tuple[ExactVector, ...]


@dataclass(frozen=True)
class MicroscopicVectorBlock:
    """One ordered local-site-irrep block from ``DISPLAY DISTORTION``."""

    global_irrep: str
    wyckoff_letter: str
    site_irrep: str
    rows: tuple[MicroscopicVectorRow, ...]
    source_order: int
    direction_symbol: str = ""

    @property
    def column_count(self) -> int:
        return len(self.rows[0].displacements) if self.rows else 0

    def __post_init__(self) -> None:
        if not self.rows:
            raise ValueError("a microscopic-vector block must contain at least one row")
        width = self.column_count
        if width <= 0:
            raise ValueError("a microscopic-vector block must contain displacement columns")
        if any(len(row.displacements) != width for row in self.rows):
            raise ValueError("microscopic-vector continuation rows have inconsistent widths")


_PROVENANCE_DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")


def _payload_digest(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def microscopic_atom_order_id(atom_ids: Sequence[str]) -> str:
    """Return the deterministic identity of one ordered physical-atom list."""

    identifiers = tuple(str(value).strip() for value in atom_ids)
    if not identifiers or any(not value for value in identifiers):
        raise ValueError("microscopic atom IDs must be nonempty")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("microscopic atom IDs must be unique")
    return f"microscopic-atom-order-{_payload_digest(identifiers)[:20]}"


def microscopic_float_column_digest(
    rows: Sequence[tuple[Sequence[str], Sequence[object]]],
) -> str:
    """Digest the deterministic float realization of one exact ISO column."""

    payload: list[dict[str, object]] = []
    for point_raw, displacement in rows:
        point = tuple(str(value) for value in point_raw)
        vector = tuple(float(value) for value in displacement)
        if len(point) != 3 or len(vector) != 3 or any(
            not np.isfinite(value) for value in vector
        ):
            raise ValueError("microscopic float-column rows must be finite 3-vectors")
        payload.append({
            "point_raw": point,
            "displacement_float_hex": tuple(value.hex() for value in vector),
        })
    if not payload:
        raise ValueError("microscopic float-column digest requires at least one row")
    return _payload_digest(payload)


@dataclass(frozen=True)
class MicroscopicColumnProvenance:
    """Immutable source evidence for one ISO microscopic displacement column.

    ``DISPLAY DISTORTION`` proves the ordered symbolic rows and their exact
    rational vectors.  It does not, by itself, identify the downstream child
    coordinate frame or a full physical-atom ordering.  Those facts therefore
    remain explicit ``None`` values instead of being inferred from the current
    structure or from list position.
    """

    query_digest: str
    query_order: int | None
    query_irrep_label: str | None
    source_parent_sg: int
    source_k_coordinates: tuple[str, ...]
    source_subgroup_space_group_number: int
    source_subgroup_basis: tuple[tuple[str, ...], ...]
    source_subgroup_origin: tuple[str, ...]
    source_subgroup_irrep_label: str
    source_subgroup_opd_symbol: str
    source_subgroup_primary_direction_selector: str
    source_block_order: int
    source_global_irrep: str
    source_wyckoff_letter: str
    source_site_irrep: str
    source_direction_symbol: str
    source_column_index: int
    source_block_column_count: int
    source_row_count: int
    exact_row_digest: str
    exact_vector_digest: str
    source_float_column_digest: str
    source_frame_id: str | None
    source_atom_order_id: str | None
    unresolved_fields: tuple[str, ...]
    unresolved_diagnostics: tuple[str, ...]

    def __post_init__(self) -> None:
        for name, value in (
            ("query_digest", self.query_digest),
            ("exact_row_digest", self.exact_row_digest),
            ("exact_vector_digest", self.exact_vector_digest),
            ("source_float_column_digest", self.source_float_column_digest),
        ):
            if _PROVENANCE_DIGEST_RE.fullmatch(str(value)) is None:
                raise ValueError(f"{name} must be one lowercase SHA-256 digest")
        if self.query_order is not None and self.query_order < 0:
            raise ValueError("microscopic query order must be non-negative")
        if self.query_irrep_label is not None and not self.query_irrep_label.strip():
            raise ValueError("resolved microscopic query irrep cannot be empty")
        if (self.query_order is None) != (self.query_irrep_label is None):
            raise ValueError(
                "microscopic query order and irrep label must resolve together"
            )
        if self.source_parent_sg <= 0:
            raise ValueError("microscopic source context requires a parent space group")
        if self.source_k_coordinates:
            if len(self.source_k_coordinates) != 3:
                raise ValueError(
                    "microscopic source context requires one exact three-component k point"
                )
            try:
                normalized_k = tuple(
                    str(Fraction(str(value).strip()))
                    for value in self.source_k_coordinates
                )
            except (ValueError, ZeroDivisionError) as exc:
                raise ValueError(
                    "microscopic source context requires exact rational k coordinates"
                ) from exc
            if tuple(self.source_k_coordinates) != normalized_k:
                raise ValueError(
                    "microscopic source k coordinates must use canonical rational tokens"
                )
        if not 1 <= self.source_subgroup_space_group_number <= 230:
            raise ValueError(
                "microscopic source context requires a child space-group number"
            )
        if (
            len(self.source_subgroup_basis) != 3
            or any(len(row) != 3 for row in self.source_subgroup_basis)
        ):
            raise ValueError("microscopic source subgroup basis must be 3x3")
        try:
            exact_basis = tuple(
                tuple(Fraction(value) for value in row)
                for row in self.source_subgroup_basis
            )
            exact_origin = tuple(
                Fraction(value) for value in self.source_subgroup_origin
            )
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(
                "microscopic source subgroup embedding must be exact rational data"
            ) from exc
        if tuple(
            tuple(str(value) for value in row) for row in exact_basis
        ) != self.source_subgroup_basis:
            raise ValueError(
                "microscopic source subgroup basis must use canonical rational tokens"
            )
        if len(exact_origin) != 3 or tuple(
            str(value) for value in exact_origin
        ) != self.source_subgroup_origin:
            raise ValueError(
                "microscopic source subgroup origin must use canonical rational tokens"
            )
        determinant = (
            exact_basis[0][0]
            * (exact_basis[1][1] * exact_basis[2][2]
               - exact_basis[1][2] * exact_basis[2][1])
            - exact_basis[0][1]
            * (exact_basis[1][0] * exact_basis[2][2]
               - exact_basis[1][2] * exact_basis[2][0])
            + exact_basis[0][2]
            * (exact_basis[1][0] * exact_basis[2][1]
               - exact_basis[1][1] * exact_basis[2][0])
        )
        if determinant == 0:
            raise ValueError("microscopic source subgroup basis must be nonsingular")
        if not (
            self.source_subgroup_irrep_label.strip()
            and self.source_subgroup_opd_symbol.strip()
        ):
            raise ValueError(
                "microscopic source context requires a subgroup irrep and OPD"
            )
        if not self.source_subgroup_primary_direction_selector.startswith(
            "VECTOR,"
        ):
            raise ValueError(
                "microscopic source context requires the primary VECTOR selector"
            )
        if not self.source_direction_symbol.startswith("VECTOR,"):
            raise ValueError(
                "microscopic source direction must retain its exact VECTOR selector"
            )
        if self.source_block_order < 0 or self.source_column_index < 0:
            raise ValueError("microscopic source block and column must be non-negative")
        if not (
            self.source_global_irrep.strip()
            and self.source_wyckoff_letter.strip()
            and self.source_site_irrep.strip()
        ):
            raise ValueError("microscopic source block identity is incomplete")
        if self.source_block_column_count <= 0:
            raise ValueError("microscopic source block must contain at least one column")
        if self.source_column_index >= self.source_block_column_count:
            raise ValueError("microscopic source column is outside its source block")
        if self.source_row_count <= 0:
            raise ValueError("microscopic source column must contain at least one row")
        if self.source_frame_id is not None and not self.source_frame_id.strip():
            raise ValueError("resolved microscopic source frame cannot be empty")
        if self.source_atom_order_id is not None and not self.source_atom_order_id.strip():
            raise ValueError("resolved microscopic source atom order cannot be empty")

        expected_fields: list[str] = []
        expected_diagnostics: list[str] = []
        if self.query_order is None:
            expected_fields.append("query_order")
            expected_diagnostics.append("query_assignment_unresolved")
        if not self.source_k_coordinates:
            expected_fields.append("source_k_coordinates")
            expected_diagnostics.append("source_k_coordinates_unresolved")
        if self.source_frame_id is None:
            expected_fields.append("source_frame_id")
            expected_diagnostics.append("source_frame_unresolved")
        if self.source_atom_order_id is None:
            expected_fields.append("source_atom_order_id")
            expected_diagnostics.append("source_atom_order_unresolved")
        if self.unresolved_fields != tuple(expected_fields):
            raise ValueError(
                "microscopic unresolved fields do not match the absent source facts"
            )
        if self.unresolved_diagnostics != tuple(expected_diagnostics):
            raise ValueError(
                "microscopic unresolved diagnostics do not match the absent source facts"
            )

    @classmethod
    def from_exact_column(
        cls,
        *,
        query_digest: str,
        query_order: int | None,
        query_irrep_label: str | None,
        source_parent_sg: int,
        source_k_coordinates: Sequence[str],
        source_subgroup_space_group_number: int,
        source_subgroup_basis: Sequence[Sequence[object]],
        source_subgroup_origin: Sequence[object],
        source_subgroup_irrep_label: str,
        source_subgroup_opd_symbol: str,
        source_subgroup_primary_direction_selector: str,
        block: MicroscopicVectorBlock,
        column_index: int,
        source_frame_id: str | None = None,
        source_atom_order_id: str | None = None,
    ) -> MicroscopicColumnProvenance:
        """Build provenance before any exact vector is converted to ``float``."""

        column = int(column_index)
        if column < 0 or column >= block.column_count:
            raise ValueError("microscopic source column is outside its source block")
        row_payload = [list(row.point_raw) for row in block.rows]
        vector_payload = [
            [
                [value.numerator, value.denominator]
                for value in row.displacements[column]
            ]
            for row in block.rows
        ]
        unresolved_fields: list[str] = []
        unresolved_diagnostics: list[str] = []
        if query_order is None:
            unresolved_fields.append("query_order")
            unresolved_diagnostics.append("query_assignment_unresolved")
        exact_k = tuple(str(value).strip() for value in source_k_coordinates)
        exact_subgroup_basis = tuple(
            tuple(str(Fraction(str(value).strip())) for value in row)
            for row in source_subgroup_basis
        )
        exact_subgroup_origin = tuple(
            str(Fraction(str(value).strip()))
            for value in source_subgroup_origin
        )
        if not exact_k:
            unresolved_fields.append("source_k_coordinates")
            unresolved_diagnostics.append("source_k_coordinates_unresolved")
        if source_frame_id is None:
            unresolved_fields.append("source_frame_id")
            unresolved_diagnostics.append("source_frame_unresolved")
        if source_atom_order_id is None:
            unresolved_fields.append("source_atom_order_id")
            unresolved_diagnostics.append("source_atom_order_unresolved")
        return cls(
            query_digest=str(query_digest),
            query_order=query_order,
            query_irrep_label=query_irrep_label,
            source_parent_sg=int(source_parent_sg),
            source_k_coordinates=exact_k,
            source_subgroup_space_group_number=int(
                source_subgroup_space_group_number
            ),
            source_subgroup_basis=exact_subgroup_basis,
            source_subgroup_origin=exact_subgroup_origin,
            source_subgroup_irrep_label=str(source_subgroup_irrep_label),
            source_subgroup_opd_symbol=str(source_subgroup_opd_symbol),
            source_subgroup_primary_direction_selector=str(
                source_subgroup_primary_direction_selector
            ),
            source_block_order=int(block.source_order),
            source_global_irrep=str(block.global_irrep),
            source_wyckoff_letter=str(block.wyckoff_letter),
            source_site_irrep=str(block.site_irrep),
            source_direction_symbol=str(block.direction_symbol),
            source_column_index=column,
            source_block_column_count=block.column_count,
            source_row_count=len(block.rows),
            exact_row_digest=_payload_digest(row_payload),
            exact_vector_digest=_payload_digest(vector_payload),
            source_float_column_digest=microscopic_float_column_digest([
                (row.point_raw, row.displacements[column]) for row in block.rows
            ]),
            source_frame_id=source_frame_id,
            source_atom_order_id=source_atom_order_id,
            unresolved_fields=tuple(unresolved_fields),
            unresolved_diagnostics=tuple(unresolved_diagnostics),
        )

    @property
    def exact_source_token(self) -> str:
        """Stable identity of the direct ISO column, independent of later mapping."""

        payload = {
            "query_digest": self.query_digest,
            "query_order": self.query_order,
            "query_irrep_label": self.query_irrep_label,
            "source_parent_sg": self.source_parent_sg,
            "source_k_coordinates": self.source_k_coordinates,
            "source_subgroup_space_group_number": (
                self.source_subgroup_space_group_number
            ),
            "source_subgroup_basis": self.source_subgroup_basis,
            "source_subgroup_origin": self.source_subgroup_origin,
            "source_subgroup_irrep_label": self.source_subgroup_irrep_label,
            "source_subgroup_opd_symbol": self.source_subgroup_opd_symbol,
            "source_subgroup_primary_direction_selector": (
                self.source_subgroup_primary_direction_selector
            ),
            "source_block_order": self.source_block_order,
            "source_global_irrep": self.source_global_irrep,
            "source_wyckoff_letter": self.source_wyckoff_letter,
            "source_site_irrep": self.source_site_irrep,
            "source_direction_symbol": self.source_direction_symbol,
            "source_column_index": self.source_column_index,
            "source_block_column_count": self.source_block_column_count,
            "source_row_count": self.source_row_count,
            "exact_row_digest": self.exact_row_digest,
            "exact_vector_digest": self.exact_vector_digest,
            "source_float_column_digest": self.source_float_column_digest,
        }
        return f"microscopic-source-{_payload_digest(payload)[:20]}"

    @property
    def stable_token(self) -> str:
        payload = {
            "exact_source_token": self.exact_source_token,
            "source_frame_id": self.source_frame_id,
            "source_atom_order_id": self.source_atom_order_id,
        }
        return f"microscopic-column-{_payload_digest(payload)[:20]}"

    @property
    def source_order_key(self) -> tuple[int, int, int]:
        if self.query_order is None:
            raise ValueError("microscopic query assignment is unresolved")
        return (
            self.query_order,
            self.source_block_order,
            self.source_column_index,
        )

    def bind_mapping(
        self,
        *,
        frame_id: str,
        atom_ids: Sequence[str],
    ) -> MicroscopicColumnProvenance:
        """Bind preserved ISO evidence to an explicit mapped frame and row order.

        This operation never invents a query assignment or changes source
        digests. It records only mapping facts supplied by the caller and
        recomputes the unresolved-field diagnostics from those facts.
        """

        resolved_frame = str(frame_id).strip()
        if not resolved_frame:
            raise ValueError("microscopic mapped frame ID cannot be empty")
        atom_order_id = microscopic_atom_order_id(atom_ids)
        unresolved_fields: list[str] = []
        unresolved_diagnostics: list[str] = []
        if self.query_order is None:
            unresolved_fields.append("query_order")
            unresolved_diagnostics.append("query_assignment_unresolved")
        if not self.source_k_coordinates:
            unresolved_fields.append("source_k_coordinates")
            unresolved_diagnostics.append("source_k_coordinates_unresolved")
        return replace(
            self,
            source_frame_id=resolved_frame,
            source_atom_order_id=atom_order_id,
            unresolved_fields=tuple(unresolved_fields),
            unresolved_diagnostics=tuple(unresolved_diagnostics),
        )


@dataclass(frozen=True)
class SiteOperation:
    """One oriented site-stabilizer element printed by ISO."""

    point_operation: str
    translation: ExactVector


@dataclass(frozen=True)
class SiteIrrepCharacter:
    """Characters of one physically irreducible real site representation."""

    label: str
    characters: tuple[Fraction, ...]


@dataclass(frozen=True)
class OrientedSiteCharacterTable:
    """Orientation-resolved character table for one parent Wyckoff position."""

    wyckoff_letter: str
    point_group: str
    operations: tuple[SiteOperation, ...]
    irreps: tuple[SiteIrrepCharacter, ...]

    def __post_init__(self) -> None:
        if not self.operations:
            raise ValueError("site character table contains no operations")
        count = len(self.operations)
        if not self.irreps:
            raise ValueError("site character table contains no irreps")
        if any(len(ir.characters) != count for ir in self.irreps):
            raise ValueError("site character row length does not match operation count")


ModeIdentityStatus: TypeAlias = Literal["verified", "partial", "unresolved"]
ModeIdentitySource: TypeAlias = Literal[
    "iso_microscopic",
    "iso_macroscopic",
    "character_projector",
    "unresolved",
]


@dataclass(frozen=True)
class ModeIdentity:
    """Scientific identity of one canonical displacement-mode column.

    ``site_irrep`` and copy/component metadata remain ``None`` when no
    auditable source determines them.  Callers must never fill those fields
    from vector direction, mode count, material name, or an expected label.
    """

    parent_sg: int
    global_irrep: str
    k_coordinates: tuple[str, ...]
    wyckoff_letter: str
    orbit_id: str = ""
    site_irrep: str | None = None
    copy_index: int | None = None
    component_index: int | None = None
    component_label: str | None = None
    source: ModeIdentitySource = "unresolved"
    status: ModeIdentityStatus = "unresolved"
    reason: str | None = None

    @property
    def display_site_irrep(self) -> str | None:
        if self.site_irrep is None:
            return None
        if self.copy_index is None:
            return self.site_irrep
        return f"{self.site_irrep}_{self.copy_index}"

    @property
    def stable_token(self) -> str:
        payload = {
            "parent_sg": self.parent_sg,
            "global_irrep": self.global_irrep,
            "k_coordinates": self.k_coordinates,
            "wyckoff_letter": self.wyckoff_letter,
            "orbit_id": self.orbit_id,
            "site_irrep": self.site_irrep,
            "copy_index": self.copy_index,
            "component_index": self.component_index,
        }
        digest = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        return f"mode-{digest[:20]}"

    @classmethod
    def unresolved(
        cls,
        *,
        parent_sg: int,
        global_irrep: str,
        k_coordinates: tuple[str, ...],
        wyckoff_letter: str,
        orbit_id: str = "",
        reason: str,
    ) -> ModeIdentity:
        return cls(
            parent_sg=parent_sg,
            global_irrep=global_irrep,
            k_coordinates=k_coordinates,
            wyckoff_letter=wyckoff_letter,
            orbit_id=orbit_id,
            source="unresolved",
            status="unresolved",
            reason=reason,
        )


def validate_microscopic_mode_source(
    identity: ModeIdentity,
    provenance: MicroscopicColumnProvenance,
    *,
    mode_id: str | None = None,
    require_resolved_provenance: bool = False,
) -> None:
    """Validate one canonical ISO microscopic column's scientific identity.

    The exact vector digest and the mode identity originate in separate ISO
    records. Keeping both objects is only useful when their shared fields are
    checked together. Frame and physical-atom-order facts may remain
    unresolved while a column is being mapped, but an export-ready caller can
    require the complete provenance explicitly.
    """

    if not isinstance(identity, ModeIdentity):
        raise ValueError("microscopic column requires a ModeIdentity")
    if not isinstance(provenance, MicroscopicColumnProvenance):
        raise ValueError("microscopic column requires exact source provenance")
    if identity.source != "iso_microscopic" or identity.status != "verified":
        raise ValueError("microscopic mode identity is not verified")
    if identity.reason is not None:
        raise ValueError("verified microscopic mode identity cannot retain an unresolved reason")
    if identity.parent_sg <= 0:
        raise ValueError("verified microscopic mode identity requires a parent space group")
    if len(identity.k_coordinates) != 3 or any(
        not str(value).strip() for value in identity.k_coordinates
    ):
        raise ValueError("verified microscopic mode identity requires one exact three-component k point")
    try:
        tuple(Fraction(str(value).strip()) for value in identity.k_coordinates)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(
            "verified microscopic mode identity requires exact rational k coordinates"
        ) from exc
    if not str(identity.global_irrep).strip():
        raise ValueError("verified microscopic mode identity requires a global irrep")
    if not str(identity.wyckoff_letter).strip():
        raise ValueError("verified microscopic mode identity requires a Wyckoff letter")
    if not str(identity.orbit_id).strip():
        raise ValueError("verified microscopic mode identity requires a physical orbit ID")
    if identity.site_irrep is None or not str(identity.site_irrep).strip():
        raise ValueError("verified microscopic mode identity requires a site irrep")
    if identity.component_index is None or identity.component_index < 0:
        raise ValueError("verified microscopic mode identity requires a component index")
    if identity.component_label is None or not str(identity.component_label).strip():
        raise ValueError("verified microscopic mode identity requires a component label")
    if mode_id is not None and str(mode_id) != identity.stable_token:
        raise ValueError("microscopic mode ID does not match its scientific identity")
    if (
        identity.parent_sg != provenance.source_parent_sg
        or identity.k_coordinates != provenance.source_k_coordinates
        or identity.component_index != provenance.source_column_index
        or identity.global_irrep != provenance.query_irrep_label
        or identity.global_irrep != provenance.source_global_irrep
        or identity.wyckoff_letter != provenance.source_wyckoff_letter
        or identity.site_irrep != provenance.source_site_irrep
    ):
        raise ValueError("microscopic mode identity and source provenance disagree")
    if require_resolved_provenance and provenance.unresolved_fields:
        raise ValueError(
            "microscopic source provenance remains unresolved: "
            + ", ".join(provenance.unresolved_fields)
        )


@dataclass(frozen=True)
class ModeSubspaceValidation:
    """Rank and principal-angle comparison of two complete mode spaces."""

    matched: bool
    reference_rank: int
    candidate_rank: int
    minimum_principal_cosine: float
    reason: str = ""


@dataclass(frozen=True)
class SignedColumnAlignment:
    """One-to-one relation ``current[j] = scale * canonical[i]``."""

    matched: bool
    canonical_rank: int
    current_rank: int
    canonical_to_current: tuple[int, ...] = ()
    current_over_canonical_scales: tuple[float, ...] = ()
    cartesian_residuals: tuple[float, ...] = ()
    maximum_cartesian_residual: float = 0.0
    reason: str = ""


@dataclass(frozen=True)
class SiteIrrepProjector:
    """Exact central projector onto one polar-vector isotypic component."""

    label: str
    matrix: ExactMatrix3
    rank: int
    character_norm: Fraction


@dataclass(frozen=True)
class IsoPointOperationCatalog:
    """The 72 ISO point-operation labels and exact matrices in ``data_space``."""

    labels: tuple[str, ...]
    matrices: tuple[ExactMatrix3, ...]
    source_sha256: str

    @classmethod
    def from_data_space(cls, path: str | Path) -> IsoPointOperationCatalog:
        source = Path(path)
        raw = source.read_bytes()
        text = raw.decode("utf-8")

        def _section(start: str, end: str) -> str:
            match = re.search(
                rf"(?ms)^{re.escape(start)}\s*$\n(.*?)^{re.escape(end)}\s*$",
                text,
            )
            if match is None:
                raise ValueError(f"missing data_space section {start!r}")
            return match.group(1)

        labels = tuple(
            value.strip()
            for value in re.findall(
                r'"([^"]*)"', _section("point_op_label", "point_op_label_stokes")
            )
        )
        values = tuple(
            int(value)
            for value in re.findall(r"[-+]?\d+", _section("ipoint_op", "ipoint_op_mlt"))
        )
        if len(labels) != 72:
            raise ValueError(f"expected 72 point-operation labels, found {len(labels)}")
        if len(values) != 72 * 3 * 3:
            raise ValueError(f"expected 648 point-operation entries, found {len(values)}")
        matrices: list[ExactMatrix3] = []
        for offset in range(0, len(values), 9):
            chunk = values[offset:offset + 9]
            matrix: ExactMatrix3 = (
                tuple(Fraction(value) for value in chunk[0:3]),  # type: ignore[assignment]
                tuple(Fraction(value) for value in chunk[3:6]),  # type: ignore[assignment]
                tuple(Fraction(value) for value in chunk[6:9]),  # type: ignore[assignment]
            )
            det = _determinant3(matrix)
            if abs(det) != 1:
                raise ValueError("data_space point operation is not unimodular")
            matrices.append(matrix)
        return cls(labels, tuple(matrices), hashlib.sha256(raw).hexdigest())

    def matrix_for_label(self, label: str) -> ExactMatrix3:
        matches = [
            matrix
            for candidate, matrix in zip(self.labels, self.matrices, strict=True)
            if candidate == label
        ]
        if not matches:
            raise KeyError(f"unknown ISO point-operation label {label!r}")
        if any(matrix != matches[0] for matrix in matches[1:]):
            raise ValueError(f"ambiguous ISO point-operation label {label!r}")
        return matches[0]


def _zero3() -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(3)] for _ in range(3)]


def _identity3() -> ExactMatrix3:
    return tuple(
        tuple(Fraction(row == column) for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _matmul(left: ExactMatrix3, right: ExactMatrix3) -> ExactMatrix3:
    return tuple(
        tuple(
            sum((left[row][k] * right[k][column] for k in range(3)), Fraction(0))
            for column in range(3)
        )
        for row in range(3)
    )  # type: ignore[return-value]


def _determinant3(matrix: ExactMatrix3) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def _rank(matrix: ExactMatrix3) -> int:
    work = [list(row) for row in matrix]
    rank = 0
    for column in range(3):
        pivot = next((row for row in range(rank, 3) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [value / scale for value in work[rank]]
        for row in range(3):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index]
                for index in range(3)
            ]
        rank += 1
    return rank


def build_site_irrep_projectors(
    table: OrientedSiteCharacterTable,
    catalog: IsoPointOperationCatalog,
) -> tuple[SiteIrrepProjector, ...]:
    """Build exact polar-vector projectors for ISO physical irreps.

    ISO reports physically irreducible *real* representations.  Their
    character norm is 1, 2, or 4 for real, complex, or quaternionic type.
    The usual complex-irrep coefficient ``d/|G|`` is therefore corrected by
    that norm; this is essential for labels such as the C4 ``E*`` irrep.
    """

    matrices = tuple(
        catalog.matrix_for_label(operation.point_operation)
        for operation in table.operations
    )
    matrix_to_index = {matrix: index for index, matrix in enumerate(matrices)}
    if len(matrix_to_index) != len(matrices):
        raise ValueError("site point group contains duplicate rotation matrices")
    inverses: list[int] = []
    identity = _identity3()
    for left in matrices:
        inverse_index = next(
            (
                index
                for index, right in enumerate(matrices)
                if _matmul(left, right) == identity and _matmul(right, left) == identity
            ),
            None,
        )
        if inverse_index is None:
            raise ValueError("site point group is not closed under inverses")
        inverses.append(inverse_index)
        for right in matrices:
            if _matmul(left, right) not in matrix_to_index:
                raise ValueError("site point group is not closed under multiplication")

    order = Fraction(len(matrices))
    projectors: list[SiteIrrepProjector] = []
    for irrep in table.irreps:
        dimension = irrep.characters[matrix_to_index[identity]]
        if dimension <= 0:
            raise ValueError(f"site irrep {irrep.label!r} has invalid dimension")
        character_norm = sum(
            (value * value for value in irrep.characters), Fraction(0)
        ) / order
        if character_norm <= 0 or character_norm.denominator != 1:
            raise ValueError(f"site irrep {irrep.label!r} has invalid character norm")
        coefficient = dimension / (character_norm * order)
        accumulator = _zero3()
        for index, matrix in enumerate(matrices):
            character = irrep.characters[inverses[index]]
            for row in range(3):
                for column in range(3):
                    accumulator[row][column] += coefficient * character * matrix[row][column]
        projector: ExactMatrix3 = tuple(
            tuple(row) for row in accumulator
        )  # type: ignore[assignment]
        if _matmul(projector, projector) != projector:
            raise ValueError(f"site-irrep projector {irrep.label!r} is not idempotent")
        projectors.append(
            SiteIrrepProjector(irrep.label, projector, _rank(projector), character_norm)
        )

    nonzero = [projector for projector in projectors if projector.rank]
    total = _zero3()
    for index, left in enumerate(nonzero):
        for right in nonzero[index + 1:]:
            if _matmul(left.matrix, right.matrix) != tuple(tuple(row) for row in _zero3()):
                raise ValueError("site-irrep projectors are not pairwise disjoint")
        for row in range(3):
            for column in range(3):
                total[row][column] += left.matrix[row][column]
    if tuple(tuple(row) for row in total) != identity:
        raise ValueError("polar-vector site-irrep projectors do not resolve the identity")
    return tuple(projectors)


def validate_signed_cartesian_columns(
    canonical: Sequence[np.ndarray],
    current: Sequence[np.ndarray],
    lattice: np.ndarray,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> SignedColumnAlignment:
    """Require a unique signed scalar and permutation for every mode column.

    Both inputs use fractional row-vector displacements.  Residuals and ranks
    are evaluated after multiplication by the row-vector Cartesian lattice.
    A whole-space basis rotation therefore fails even when the two spans are
    identical: this predicate establishes column provenance, not merely span
    equivalence.
    """

    lattice_array = np.asarray(lattice, dtype=float)
    if (
        lattice_array.shape != (3, 3)
        or not np.all(np.isfinite(lattice_array))
        or abs(float(np.linalg.det(lattice_array))) <= np.finfo(float).eps
    ):
        raise ValueError("lattice must be one finite nonsingular 3x3 matrix")
    if not np.isfinite(rtol) or rtol <= 0 or not np.isfinite(atol) or atol <= 0:
        raise ValueError("column-alignment tolerances must be finite and positive")
    if len(canonical) == 0 or len(current) == 0:
        return SignedColumnAlignment(
            False, 0, 0, reason="empty_column_set"
        )
    if len(canonical) != len(current):
        return SignedColumnAlignment(
            False, 0, 0, reason="column_count_mismatch"
        )

    atom_count: int | None = None

    def _flatten(items: Sequence[np.ndarray]) -> list[np.ndarray]:
        nonlocal atom_count
        result: list[np.ndarray] = []
        for item in items:
            array = np.asarray(item, dtype=float)
            if array.ndim != 2 or array.shape[1] != 3 or not np.all(np.isfinite(array)):
                raise ValueError("each mode column must be one finite (n_atoms, 3) array")
            if atom_count is None:
                atom_count = array.shape[0]
            elif array.shape[0] != atom_count:
                raise ValueError("all mode columns must use one atom ordering")
            result.append((array @ lattice_array).reshape(-1))
        return result

    canonical_vectors = _flatten(canonical)
    current_vectors = _flatten(current)
    canonical_matrix = np.column_stack(canonical_vectors)
    current_matrix = np.column_stack(current_vectors)

    def _rank(matrix: np.ndarray) -> int:
        singular = np.linalg.svd(matrix, compute_uv=False)
        threshold = max(
            atol,
            rtol * (float(singular[0]) if singular.size else 0.0),
        )
        return int(np.count_nonzero(singular > threshold))

    canonical_rank = _rank(canonical_matrix)
    current_rank = _rank(current_matrix)
    if canonical_rank != len(canonical_vectors):
        return SignedColumnAlignment(
            False,
            canonical_rank,
            current_rank,
            reason="canonical_columns_not_independent",
        )
    if current_rank != len(current_vectors):
        return SignedColumnAlignment(
            False,
            canonical_rank,
            current_rank,
            reason="current_columns_not_independent",
        )

    matches: list[list[tuple[int, float, float]]] = []
    for canonical_vector in canonical_vectors:
        norm_squared = float(canonical_vector @ canonical_vector)
        if norm_squared <= atol * atol:
            return SignedColumnAlignment(
                False,
                canonical_rank,
                current_rank,
                reason="canonical_zero_column",
            )
        candidates: list[tuple[int, float, float]] = []
        canonical_norm = float(np.sqrt(norm_squared))
        for index, current_vector in enumerate(current_vectors):
            scale = float((canonical_vector @ current_vector) / norm_squared)
            if (
                not np.isfinite(scale)
                or abs(scale) * canonical_norm <= atol
            ):
                continue
            residual = float(np.linalg.norm(current_vector - scale * canonical_vector))
            current_norm = float(np.linalg.norm(current_vector))
            threshold = atol + rtol * max(current_norm, abs(scale) * canonical_norm)
            if residual <= threshold:
                candidates.append((index, scale, residual))
        matches.append(candidates)

    if any(not candidates for candidates in matches):
        return SignedColumnAlignment(
            False,
            canonical_rank,
            current_rank,
            reason="non_diagonal_column_mixing",
        )
    if any(len(candidates) != 1 for candidates in matches):
        return SignedColumnAlignment(
            False,
            canonical_rank,
            current_rank,
            reason="non_unique_scalar_column_mapping",
        )
    mapping = tuple(candidates[0][0] for candidates in matches)
    if len(set(mapping)) != len(mapping):
        return SignedColumnAlignment(
            False,
            canonical_rank,
            current_rank,
            reason="non_bijective_scalar_column_mapping",
        )
    scales = tuple(candidates[0][1] for candidates in matches)
    residuals = tuple(candidates[0][2] for candidates in matches)
    return SignedColumnAlignment(
        True,
        canonical_rank,
        current_rank,
        canonical_to_current=mapping,
        current_over_canonical_scales=scales,
        cartesian_residuals=residuals,
        maximum_cartesian_residual=max(residuals, default=0.0),
    )


def validate_cartesian_subspaces(
    reference: list[np.ndarray] | tuple[np.ndarray, ...],
    candidate: list[np.ndarray] | tuple[np.ndarray, ...],
    lattice: np.ndarray,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> ModeSubspaceValidation:
    """Compare complete displacement spaces by rank and principal angles."""

    lattice_array = np.asarray(lattice, dtype=float)
    if lattice_array.shape != (3, 3) or not np.all(np.isfinite(lattice_array)):
        raise ValueError("lattice must be a finite 3x3 matrix")

    def _columns(items: list[np.ndarray] | tuple[np.ndarray, ...]) -> np.ndarray:
        vectors: list[np.ndarray] = []
        atom_count: int | None = None
        for item in items:
            array = np.asarray(item, dtype=float)
            if array.ndim != 2 or array.shape[1] != 3:
                raise ValueError("each mode must have shape (n_atoms, 3)")
            if atom_count is None:
                atom_count = array.shape[0]
            elif array.shape[0] != atom_count:
                raise ValueError("all modes must use the same atom ordering")
            vectors.append((array @ lattice_array).reshape(-1))
        if not vectors:
            return np.zeros((0, 0), dtype=float)
        return np.column_stack(vectors)

    left = _columns(reference)
    right = _columns(candidate)
    if left.size == 0 and right.size == 0:
        return ModeSubspaceValidation(True, 0, 0, 1.0)
    if left.shape[0] != right.shape[0]:
        return ModeSubspaceValidation(False, 0, 0, 0.0, "atom_dimension_mismatch")

    def _basis(matrix: np.ndarray) -> tuple[np.ndarray, int]:
        if matrix.size == 0:
            return np.zeros((matrix.shape[0], 0)), 0
        u, singular, _vh = np.linalg.svd(matrix, full_matrices=False)
        threshold = max(atol, rtol * (float(singular[0]) if singular.size else 0.0))
        rank = int(np.count_nonzero(singular > threshold))
        return u[:, :rank], rank

    q_left, rank_left = _basis(left)
    q_right, rank_right = _basis(right)
    if rank_left != rank_right:
        return ModeSubspaceValidation(
            False, rank_left, rank_right, 0.0, "rank_mismatch"
        )
    if rank_left == 0:
        return ModeSubspaceValidation(True, 0, 0, 1.0)
    cosines = np.linalg.svd(q_left.T @ q_right, compute_uv=False)
    minimum = float(np.min(cosines)) if cosines.size else 0.0
    matched = minimum >= 1.0 - max(rtol * 10.0, 1e-10)
    return ModeSubspaceValidation(
        matched,
        rank_left,
        rank_right,
        minimum,
        "" if matched else "principal_angle_mismatch",
    )


__all__ = [
    "ExactMatrix3",
    "ExactVector",
    "InvariantDirection",
    "IsoPointOperationCatalog",
    "MacroscopicTensorBasis",
    "MacroscopicTensorBlock",
    "MicroscopicColumnProvenance",
    "MicroscopicVectorBlock",
    "MicroscopicVectorRow",
    "ModeIdentity",
    "ModeSubspaceValidation",
    "OrientedSiteCharacterTable",
    "SignedColumnAlignment",
    "SiteIrrepCharacter",
    "SiteIrrepProjector",
    "SiteOperation",
    "SymbolicTensorComponent",
    "build_site_irrep_projectors",
    "microscopic_atom_order_id",
    "microscopic_float_column_digest",
    "validate_cartesian_subspaces",
    "validate_microscopic_mode_source",
    "validate_signed_cartesian_columns",
]
