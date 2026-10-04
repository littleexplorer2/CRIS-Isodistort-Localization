"""Pure scientific model for displacive-mode CIF loops.

The CIF writer needs three related objects in one and the same child setting:

* exact Seitz operations, acting on fractional-coordinate columns as
  ``x' = R x + t``;
* a complete conventional-cell atom list and its asymmetric representatives;
* authoritative, unmixed mode directions and signs in that atom order/frame.

This module derives site displacement constraints and orbit transporters,
canonicalizes every free-coordinate column to max-component-one, and computes
the matching primitive-cell Cartesian normalization.
It performs no I/O and knows no space-group names or example structures.

Fractional positions and Seitz data use :class:`fractions.Fraction`.  Raw mode
coefficients may be irrational, so projection and Cartesian metric checks use
finite floats only at that boundary.  Every ambiguous or incomplete mapping
raises :class:`DisplaciveModelError`; callers must not emit official-looking
mode loops from a partial result.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral, Real

import numpy as np

from ..backend.iso_mode_models import (
    MicroscopicColumnProvenance,
    ModeIdentity,
    microscopic_atom_order_id,
    validate_microscopic_mode_source,
)
from ..backend.iso_wrapper import (
    MicroscopicDomainExtensionEvidence,
    validate_microscopic_domain_extension_bindings,
)
from ..utils.lattice import (
    RationalMatrix,
    as_fraction,
    determinant,
    identity_matrix,
    inverse,
    multiply,
    rational_matrix,
    transpose,
)

RationalVector = tuple[Fraction, Fraction, Fraction]
FloatVector = tuple[float, float, float]

_AXIS_NAMES = ("x", "y", "z")
_PARAMETER_NAMES = ("Dx", "Dy", "Dz")
_AS_AMPLITUDE_CONVENTION = "As_unit_mode_angstrom"


class DisplaciveModelError(ValueError):
    """A fail-closed scientific-contract violation.

    ``code`` is stable enough for an API layer to translate into a diagnostic
    without inspecting the human-readable message.
    """

    def __init__(self, code: str, message: str) -> None:
        self.code = str(code)
        super().__init__(f"{self.code}: {message}")


def _fail(code: str, message: str) -> None:
    raise DisplaciveModelError(code, message)


def _identifier(value: object, description: str, code: str) -> str:
    identifier = str(value).strip()
    if not identifier:
        _fail(code, f"{description} cannot be empty or whitespace")
    return identifier


def _nonnegative_integer(value: object, description: str, code: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        _fail(code, f"{description} must be an integer")
    parsed = int(value)
    if parsed < 0:
        _fail(code, f"{description} must be nonnegative")
    return parsed


def _vector3(values: Sequence[object], description: str) -> RationalVector:
    if len(values) != 3:
        _fail("invalid_exact_vector", f"{description} must contain three values")
    try:
        result = tuple(as_fraction(value) for value in values)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise DisplaciveModelError("invalid_exact_vector", f"{description} is not an exact finite vector") from exc
    return result  # type: ignore[return-value]


def _float_vector3(values: Sequence[object], description: str) -> FloatVector:
    if len(values) != 3:
        _fail("invalid_float_vector", f"{description} must contain three values")
    try:
        result = tuple(float(value) for value in values)
    except (TypeError, ValueError, OverflowError) as exc:
        raise DisplaciveModelError("invalid_float_vector", f"{description} is not a finite vector") from exc
    if not all(math.isfinite(value) for value in result):
        _fail("invalid_float_vector", f"{description} is not a finite vector")
    return result  # type: ignore[return-value]


def _mod_one(vector: RationalVector) -> RationalVector:
    return tuple(value % 1 for value in vector)  # type: ignore[return-value]


def _add(left: RationalVector, right: RationalVector) -> RationalVector:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def _subtract(left: RationalVector, right: RationalVector) -> RationalVector:
    return tuple(left[index] - right[index] for index in range(3))  # type: ignore[return-value]


def _matrix_vector(matrix: RationalMatrix, vector: RationalVector) -> RationalVector:
    return tuple(sum(matrix[row][column] * vector[column] for column in range(3)) for row in range(3))  # type: ignore[return-value]


def _rotation_key(rotation: RationalMatrix) -> tuple[tuple[tuple[int, int], ...], ...]:
    return tuple(tuple((value.numerator, value.denominator) for value in row) for row in rotation)


def _vector_key(vector: RationalVector) -> tuple[tuple[int, int], ...]:
    return tuple((value.numerator, value.denominator) for value in vector)


@dataclass(frozen=True)
class ExactSeitzOperation:
    """One exact emitted-setting Seitz operation modulo integer translations."""

    rotation: RationalMatrix
    translation: RationalVector
    label: str = ""

    @classmethod
    def from_values(
        cls,
        rotation: Sequence[Sequence[object]],
        translation: Sequence[object],
        *,
        label: str = "",
    ) -> ExactSeitzOperation:
        try:
            exact_rotation = rational_matrix(rotation)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise DisplaciveModelError(
                "invalid_seitz_rotation", "Seitz rotation must be one exact nonsingular 3x3 matrix"
            ) from exc
        if any(value.denominator != 1 for row in exact_rotation for value in row):
            _fail("invalid_seitz_rotation", "emitted fractional rotation must be integral")
        if abs(determinant(exact_rotation)) != 1:
            _fail("invalid_seitz_rotation", "emitted fractional rotation must be unimodular")
        exact_translation = _mod_one(_vector3(translation, "Seitz translation"))
        return cls(exact_rotation, exact_translation, str(label))

    @property
    def key(
        self,
    ) -> tuple[tuple[tuple[tuple[int, int], ...], ...], tuple[tuple[int, int], ...]]:
        return _rotation_key(self.rotation), _vector_key(self.translation)

    def operate(self, coordinate: RationalVector) -> RationalVector:
        return _mod_one(_add(_matrix_vector(self.rotation, coordinate), self.translation))

    def transform_displacement(self, displacement: RationalVector) -> RationalVector:
        return _matrix_vector(self.rotation, displacement)


def _compose(left: ExactSeitzOperation, right: ExactSeitzOperation) -> ExactSeitzOperation:
    """Return ``left * right`` for column-coordinate affine actions."""

    rotation = multiply(left.rotation, right.rotation)
    translation = _mod_one(_add(_matrix_vector(left.rotation, right.translation), left.translation))
    return ExactSeitzOperation(rotation, translation)


def canonicalize_emitted_seitz(
    operations: Sequence[ExactSeitzOperation], *, reject_duplicates: bool = True
) -> tuple[ExactSeitzOperation, ...]:
    """Canonicalize and exactly validate a finite conventional-cell Seitz group."""

    if not operations:
        _fail("missing_seitz_operations", "the emitted child Seitz set is empty")
    by_key: dict[object, ExactSeitzOperation] = {}
    for operation in operations:
        if not isinstance(operation, ExactSeitzOperation):
            _fail("invalid_seitz_operation", "all operations must be ExactSeitzOperation values")
        canonical = ExactSeitzOperation.from_values(operation.rotation, operation.translation, label=operation.label)
        if canonical.key in by_key and reject_duplicates:
            _fail("duplicate_seitz_operation", "the emitted Seitz set contains a duplicate modulo 1")
        by_key.setdefault(canonical.key, canonical)
    ordered = tuple(by_key[key] for key in sorted(by_key))
    keys = {operation.key for operation in ordered}
    identity = ExactSeitzOperation.from_values(identity_matrix(), (0, 0, 0)).key
    if identity not in keys:
        _fail("missing_seitz_identity", "the emitted Seitz set has no identity operation")
    for left in ordered:
        for right in ordered:
            if _compose(left, right).key not in keys:
                _fail("seitz_group_not_closed", "the emitted Seitz set is not closed modulo 1")
    return ordered


def conjugate_seitz_operations(
    operations: Sequence[ExactSeitzOperation],
    coordinate_transform: Sequence[Sequence[object]],
    origin_shift: Sequence[object],
) -> tuple[ExactSeitzOperation, ...]:
    """Conjugate operations under ``x_new = A x_old + q`` exactly."""

    source = canonicalize_emitted_seitz(operations)
    try:
        transform = rational_matrix(coordinate_transform)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise DisplaciveModelError(
            "invalid_setting_transform", "setting transform must be one exact nonsingular 3x3 matrix"
        ) from exc
    if abs(determinant(transform)) != 1 or any(value.denominator != 1 for row in transform for value in row):
        _fail("invalid_setting_transform", "setting transform must belong to GL(3,Z)")
    transform_inverse = inverse(transform)
    shift = _vector3(origin_shift, "setting origin shift")
    result: list[ExactSeitzOperation] = []
    for operation in source:
        rotation = multiply(multiply(transform, operation.rotation), transform_inverse)
        translation = _add(
            _matrix_vector(transform, operation.translation),
            _subtract(shift, _matrix_vector(rotation, shift)),
        )
        result.append(ExactSeitzOperation.from_values(rotation, translation))
    return canonicalize_emitted_seitz(result)


def transform_fractional_position(
    coordinate: Sequence[object],
    coordinate_transform: Sequence[Sequence[object]],
    origin_shift: Sequence[object],
) -> RationalVector:
    """Apply ``x_new = A x_old + q`` to one exact fractional position."""

    transform = rational_matrix(coordinate_transform)
    old = _vector3(coordinate, "fractional position")
    shift = _vector3(origin_shift, "setting origin shift")
    return _mod_one(_add(_matrix_vector(transform, old), shift))


def transform_fractional_displacement(
    displacement: Sequence[object], coordinate_transform: Sequence[Sequence[object]]
) -> RationalVector:
    """Apply ``u_new = A u_old``; an origin shift never enters a displacement."""

    transform = rational_matrix(coordinate_transform)
    return _matrix_vector(transform, _vector3(displacement, "fractional displacement"))


def transformed_reference_lattice(
    reference_lattice: Sequence[Sequence[object]],
    coordinate_transform: Sequence[Sequence[object]],
) -> tuple[FloatVector, FloatVector, FloatVector]:
    """Return ``L_new=A^-T L_old`` for row-vector Cartesian lattices."""

    lattice = _finite_lattice(reference_lattice)
    transform = rational_matrix(coordinate_transform)
    inverse_transpose = np.asarray(transpose(inverse(transform)), dtype=float)
    changed = inverse_transpose @ np.asarray(lattice, dtype=float)
    return tuple(tuple(float(value) for value in row) for row in changed)  # type: ignore[return-value]


def _rref(rows: Sequence[Sequence[Fraction]], *, column_count: int) -> tuple[list[list[Fraction]], tuple[int, ...]]:
    matrix = [list(row) for row in rows]
    if any(len(row) != column_count for row in matrix):
        _fail("invalid_rational_matrix", "RREF rows have inconsistent lengths")
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        source = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column] != 0),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [matrix[row][index] - factor * matrix[pivot_row][index] for index in range(column_count)]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    nonzero = [row for row in matrix if any(value != 0 for value in row)]
    return nonzero, tuple(pivot_columns)


def _nullspace3(rows: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    reduced, pivots = _rref(rows, column_count=3)
    free_columns = [column for column in range(3) if column not in pivots]
    columns: list[list[Fraction]] = []
    for free in free_columns:
        vector = [Fraction(0), Fraction(0), Fraction(0)]
        vector[free] = Fraction(1)
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row_index][free]
        columns.append(vector)
    return columns


def _canonical_fixed_basis(
    constraint_rows: Sequence[Sequence[Fraction]],
) -> tuple[tuple[tuple[Fraction, ...], ...], tuple[int, ...]]:
    columns = _nullspace3(constraint_rows)
    dimension = len(columns)
    if dimension == 0:
        return ((), (), ()), ()
    basis_rows = [[columns[column][axis] for axis in range(3)] for column in range(dimension)]
    reduced_rows, pivots = _rref(basis_rows, column_count=3)
    if len(reduced_rows) != dimension or len(pivots) != dimension:
        _fail("fixed_space_rank_failure", "could not canonicalize the exact site fixed space")
    basis = tuple(tuple(reduced_rows[column][axis] for column in range(dimension)) for axis in range(3))
    return basis, pivots


def _format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _format_linear_expression(coefficients: Sequence[Fraction], parameters: Sequence[str]) -> str:
    terms: list[str] = []
    for coefficient, parameter in zip(coefficients, parameters, strict=True):
        if coefficient == 0:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = abs(coefficient)
        body = parameter if magnitude == 1 else f"{_format_fraction(magnitude)}*{parameter}"
        if not terms:
            terms.append(body if sign == "+" else f"-{body}")
        else:
            terms.append(f"{sign}{body}")
    return "".join(terms) or "0"


@dataclass(frozen=True)
class DeltaParameter:
    axis: int
    symbol: str
    suffix: str


@dataclass(frozen=True)
class SiteCoordinateModel:
    representative_id: str
    representative_frac: RationalVector
    stabilizer_indices: tuple[int, ...]
    fixed_basis: tuple[tuple[Fraction, ...], ...]
    pivot_axes: tuple[int, ...]
    parameters: tuple[DeltaParameter, ...]
    component_expressions: tuple[str, str, str]
    symmform: str

    @property
    def dimension(self) -> int:
        return len(self.pivot_axes)


def site_fixed_coordinate_model(
    representative_id: str,
    representative_frac: Sequence[object],
    operations: Sequence[ExactSeitzOperation],
) -> SiteCoordinateModel:
    """Derive one site's exact polar-displacement fixed space and CIF form."""

    emitted = canonicalize_emitted_seitz(operations)
    coordinate = _mod_one(_vector3(representative_frac, "site representative"))
    stabilizer = tuple(index for index, operation in enumerate(emitted) if operation.operate(coordinate) == coordinate)
    if not stabilizer:
        _fail("missing_site_stabilizer", f"site {representative_id!r} has no stabilizer identity")
    constraints: list[list[Fraction]] = []
    for index in stabilizer:
        rotation = emitted[index].rotation
        for row in range(3):
            constraints.append([rotation[row][column] - (1 if row == column else 0) for column in range(3)])
    basis, pivots = _canonical_fixed_basis(constraints)
    parameters = tuple(DeltaParameter(axis, _PARAMETER_NAMES[axis], f"d{_AXIS_NAMES[axis]}") for axis in pivots)
    names = tuple(parameter.symbol for parameter in parameters)
    expressions = tuple(_format_linear_expression(basis[axis], names) for axis in range(3))
    return SiteCoordinateModel(
        representative_id=str(representative_id),
        representative_frac=coordinate,
        stabilizer_indices=stabilizer,
        fixed_basis=basis,
        pivot_axes=pivots,
        parameters=parameters,
        component_expressions=expressions,  # type: ignore[arg-type]
        symmform=",".join(expressions),
    )


ChemicalSignature = tuple[tuple[str, Fraction], ...]


def _chemical_signature(
    species: str | Mapping[str, object] | Sequence[tuple[str, object]],
) -> ChemicalSignature:
    if isinstance(species, str):
        items = [(species, Fraction(1))]
    elif isinstance(species, Mapping):
        items = [(str(name), as_fraction(value)) for name, value in species.items()]
    else:
        try:
            items = [(str(name), as_fraction(value)) for name, value in species]
        except (TypeError, ValueError) as exc:
            raise DisplaciveModelError(
                "invalid_chemical_signature", "chemical signature must contain species/occupancy pairs"
            ) from exc
    if not items or any(not name or value <= 0 or value > 1 for name, value in items):
        _fail("invalid_chemical_signature", "species occupancies must be in the interval (0, 1]")
    combined: dict[str, Fraction] = {}
    for name, value in items:
        combined[name] = combined.get(name, Fraction(0)) + value
    if sum(combined.values(), Fraction(0)) > 1:
        _fail("invalid_chemical_signature", "total site occupancy exceeds 1")
    return tuple(sorted(combined.items()))


@dataclass(frozen=True)
class ChildAtom:
    atom_id: str
    frac: RationalVector
    chemical_signature: ChemicalSignature

    @classmethod
    def from_values(
        cls,
        atom_id: str,
        frac: Sequence[object],
        species: str | Mapping[str, object] | Sequence[tuple[str, object]] = "X",
    ) -> ChildAtom:
        identifier = _identifier(atom_id, "child atom ID", "invalid_atom_id")
        return cls(identifier, _mod_one(_vector3(frac, "child atom coordinate")), _chemical_signature(species))


@dataclass(frozen=True)
class AtomTransporter:
    member_atom_id: str
    candidate_operation_indices: tuple[int, ...]
    canonical_operation_index: int
    induced_fixed_map: tuple[tuple[Fraction, ...], ...]


@dataclass(frozen=True)
class ChildOrbitCoordinate:
    site: SiteCoordinateModel
    representative_atom_id: str
    member_atom_ids: tuple[str, ...]
    transporters: tuple[AtomTransporter, ...]


def _rotation_fixed_map(
    rotation: RationalMatrix, basis: tuple[tuple[Fraction, ...], ...]
) -> tuple[tuple[Fraction, ...], ...]:
    dimension = len(basis[0]) if basis else 0
    return tuple(
        tuple(sum(rotation[row][axis] * basis[axis][column] for axis in range(3)) for column in range(dimension))
        for row in range(3)
    )


def build_child_orbits(
    atoms: Sequence[ChildAtom],
    representative_atom_ids: Sequence[str],
    operations: Sequence[ExactSeitzOperation],
) -> tuple[ChildOrbitCoordinate, ...]:
    """Partition atoms and prove a unique induced transporter for every member."""

    emitted = canonicalize_emitted_seitz(operations)
    atom_by_id: dict[str, ChildAtom] = {}
    for atom in atoms:
        if atom.atom_id in atom_by_id:
            _fail("duplicate_atom_id", f"child atom ID {atom.atom_id!r} is duplicated")
        atom_by_id[atom.atom_id] = atom
    if not atom_by_id:
        _fail("missing_child_atoms", "the conventional child atom list is empty")
    representative_ids = tuple(str(value) for value in representative_atom_ids)
    if not representative_ids or len(set(representative_ids)) != len(representative_ids):
        _fail("invalid_representatives", "representative atom IDs must be nonempty and unique")
    if any(identifier not in atom_by_id for identifier in representative_ids):
        _fail("invalid_representatives", "a representative atom ID is absent from the child atom list")

    owners: dict[str, str] = {}
    orbits: list[ChildOrbitCoordinate] = []
    for representative_id in representative_ids:
        representative = atom_by_id[representative_id]
        site = site_fixed_coordinate_model(representative_id, representative.frac, emitted)
        transports: list[AtomTransporter] = []
        generated_positions = {operation.operate(representative.frac) for operation in emitted}
        for atom in atoms:
            if atom.chemical_signature != representative.chemical_signature:
                continue
            candidates = tuple(
                index for index, operation in enumerate(emitted) if operation.operate(representative.frac) == atom.frac
            )
            if not candidates:
                continue
            previous_owner = owners.get(atom.atom_id)
            if previous_owner is not None:
                _fail(
                    "overlapping_child_orbits",
                    f"atom {atom.atom_id!r} belongs to representatives {previous_owner!r} and {representative_id!r}",
                )
            induced = tuple(_rotation_fixed_map(emitted[index].rotation, site.fixed_basis) for index in candidates)
            if any(value != induced[0] for value in induced[1:]):
                _fail(
                    "ambiguous_transporter_action",
                    f"candidate transporters to atom {atom.atom_id!r} disagree on the site fixed space",
                )
            owners[atom.atom_id] = representative_id
            transports.append(
                AtomTransporter(
                    member_atom_id=atom.atom_id,
                    candidate_operation_indices=candidates,
                    canonical_operation_index=candidates[0],
                    induced_fixed_map=induced[0],
                )
            )
        if len(transports) != len(generated_positions):
            _fail(
                "incomplete_child_orbit",
                f"representative {representative_id!r} generates "
                f"{len(generated_positions)} positions but maps to "
                f"{len(transports)} atoms",
            )
        transports.sort(key=lambda item: item.member_atom_id)
        orbits.append(
            ChildOrbitCoordinate(
                site=site,
                representative_atom_id=representative_id,
                member_atom_ids=tuple(item.member_atom_id for item in transports),
                transporters=tuple(transports),
            )
        )
    missing = sorted(set(atom_by_id).difference(owners))
    if missing:
        _fail("unassigned_child_atoms", f"child atoms have no asymmetric representative: {', '.join(missing)}")
    return tuple(orbits)


@dataclass(frozen=True)
class ModeScaleProvenance:
    source: str
    convention: str
    direction_resolved: bool
    sign_resolved: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.source, str) or not isinstance(self.convention, str):
            _fail("missing_mode_scale_provenance", "mode source and scale convention must be text")
        if not self.source.strip() or not self.convention.strip():
            _fail("missing_mode_scale_provenance", "mode source and scale convention are required")


@dataclass(frozen=True)
class RawModeColumn:
    mode_id: str
    column_index: int
    frame_id: str
    atom_ids: tuple[str, ...]
    displacements: tuple[FloatVector, ...]
    scale_provenance: ModeScaleProvenance
    label: str = ""
    mode_identity: ModeIdentity | None = None
    microscopic_provenance: MicroscopicColumnProvenance | None = None
    microscopic_domain_extension: MicroscopicDomainExtensionEvidence | None = None

    @classmethod
    def from_values(
        cls,
        mode_id: str,
        column_index: int,
        frame_id: str,
        atom_ids: Sequence[str],
        displacements: Sequence[Sequence[object]],
        scale_provenance: ModeScaleProvenance,
        *,
        mode_identity: ModeIdentity | None = None,
        microscopic_provenance: MicroscopicColumnProvenance | None = None,
        microscopic_domain_extension: MicroscopicDomainExtensionEvidence | None = None,
        label: str = "",
    ) -> RawModeColumn:
        identifier = _identifier(mode_id, "mode ID", "invalid_mode_identity")
        column = _nonnegative_integer(column_index, "mode column index", "invalid_mode_identity")
        frame_identifier = _identifier(frame_id, "mode frame ID", "invalid_mode_frame")
        identifiers = tuple(_identifier(value, "mode atom ID", "invalid_mode_atom_order") for value in atom_ids)
        vectors = tuple(_float_vector3(row, f"mode {mode_id!r} displacement") for row in displacements)
        if len(identifiers) != len(vectors) or len(set(identifiers)) != len(identifiers):
            _fail("invalid_mode_atom_order", "mode atom IDs must be unique and match displacement rows")
        if not isinstance(scale_provenance, ModeScaleProvenance):
            _fail("missing_mode_scale_provenance", "mode scale provenance has the wrong type")
        return cls(
            mode_id=identifier,
            column_index=column,
            frame_id=frame_identifier,
            atom_ids=identifiers,
            displacements=vectors,
            scale_provenance=scale_provenance,
            mode_identity=mode_identity,
            microscopic_provenance=microscopic_provenance,
            microscopic_domain_extension=microscopic_domain_extension,
            label=str(label),
        )

    @property
    def global_irrep(self) -> str:
        """Return the source-bound global irrep, or empty text for diagnostics."""

        if self.mode_identity is None:
            return ""
        return str(self.mode_identity.global_irrep).strip()

    @property
    def display_label(self) -> str:
        """Return the canonical source display label used by export adapters."""

        return str(self.label).strip()

    def source_binding_issue(self) -> str | None:
        """Return why identity/provenance do not bind this mapped column."""

        if self.mode_identity is None:
            return "missing verified microscopic ModeIdentity"
        if self.microscopic_provenance is None:
            return "missing microscopic source-column provenance"
        extension = self.microscopic_domain_extension
        expected_extension_source = (
            "iso_microscopic_with_verified_bush_domain_extension"
        )
        expected_direct_source = "iso_microscopic_display_distortion"
        scale_source = self.scale_provenance.source.strip()
        if extension is not None and scale_source != expected_extension_source:
            return "microscopic scale source omits verified domain extension"
        if extension is None and scale_source != expected_direct_source:
            if scale_source == expected_extension_source:
                return "microscopic scale source claims missing domain extension"
            return "unknown direct microscopic scale source"
        if extension is not None and (
            extension.canonical_source_column_token
            != self.microscopic_provenance.exact_source_token
        ):
            return "microscopic domain extension is bound to another source column"
        if extension is not None and (
            extension.canonical_source_order_key
            != self.microscopic_provenance.source_order_key
        ):
            return "microscopic domain extension is bound to another source index"
        try:
            validate_microscopic_mode_source(
                self.mode_identity,
                self.microscopic_provenance,
                mode_id=self.mode_id,
                require_resolved_provenance=True,
            )
        except ValueError as exc:
            return str(exc)
        if self.microscopic_provenance.source_frame_id != self.frame_id:
            return "microscopic source frame differs from the raw mode frame"
        if (
            self.microscopic_provenance.source_atom_order_id
            != microscopic_atom_order_id(self.atom_ids)
        ):
            return "microscopic source atom order differs from the raw mode rows"
        return None

    def export_validation_issue(self) -> str | None:
        """Return why this source column cannot enter a validated export model."""

        binding_issue = self.source_binding_issue()
        if binding_issue is not None:
            return binding_issue
        if not self.scale_provenance.direction_resolved:
            return "mode direction is unresolved"
        if not self.scale_provenance.sign_resolved:
            return "mode sign is unresolved"
        if not self.display_label:
            return "canonical mode display label is empty"
        return None

    @property
    def has_resolved_source(self) -> bool:
        """Whether identity, source evidence, scale and label are export-ready."""

        return self.export_validation_issue() is None


@dataclass(frozen=True)
class FreeCoordinateRow:
    representative_atom_id: str
    local_parameter_index: int
    axis: int
    symbol: str
    delta_label: str


@dataclass(frozen=True)
class ModeProjectionDiagnostics:
    free_coordinate_count: int
    mode_count: int
    matrix_rank: int
    minimum_singular_value: float
    max_transporter_residual_angstrom: float
    max_fixed_space_residual_angstrom: float
    max_reconstruction_residual_angstrom: float


@dataclass(frozen=True)
class ModeProjection:
    free_rows: tuple[FreeCoordinateRow, ...]
    mode_ids: tuple[str, ...]
    matrix: tuple[tuple[float, ...], ...]
    input_to_canonical_scale_factors: tuple[float, ...]
    diagnostics: ModeProjectionDiagnostics

    def matrix_array(self) -> np.ndarray:
        return np.asarray(self.matrix, dtype=float).reshape(
            self.diagnostics.free_coordinate_count, self.diagnostics.mode_count
        )


def _finite_lattice(
    rows: Sequence[Sequence[object]],
) -> tuple[FloatVector, FloatVector, FloatVector]:
    if len(rows) != 3:
        _fail("invalid_reference_lattice", "reference lattice must be 3x3")
    lattice = tuple(_float_vector3(row, "reference lattice row") for row in rows)
    matrix = np.asarray(lattice, dtype=float)
    if abs(float(np.linalg.det(matrix))) <= np.finfo(float).eps:
        _fail("invalid_reference_lattice", "reference lattice must be nonsingular")
    return lattice  # type: ignore[return-value]


def _positive_finite(value: float, description: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        _fail("invalid_numeric_tolerance", f"{description} must be finite and positive")
    return parsed


def _cartesian_norm(vector: np.ndarray, lattice: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(vector, dtype=float) @ lattice))


def _validate_metric_compatibility(
    operations: Sequence[ExactSeitzOperation],
    reference_lattice: Sequence[Sequence[object]],
    relative_tolerance: float,
) -> None:
    """Require every fractional rotation to preserve the row-lattice metric.

    Fractional coordinates are columns while lattice vectors are rows, so
    ``G = L L.T`` and a valid child rotation obeys ``R.T G R = G``.  The
    tolerance is dimensionless and relative to the largest entry of ``G``;
    the decision is unchanged when all lengths use another consistent unit.
    """

    tolerance = _positive_finite(relative_tolerance, "reference metric relative tolerance")
    lattice = np.asarray(_finite_lattice(reference_lattice), dtype=float)
    metric = lattice @ lattice.T
    metric_scale = float(np.max(np.abs(metric)))
    if not math.isfinite(metric_scale) or metric_scale <= 0:
        _fail("invalid_reference_lattice", "reference lattice has an invalid Cartesian metric")
    threshold = tolerance * metric_scale
    seen: set[object] = set()
    for operation in operations:
        rotation_key = _rotation_key(operation.rotation)
        if rotation_key in seen:
            continue
        seen.add(rotation_key)
        rotation = np.asarray(operation.rotation, dtype=float)
        residual = float(np.max(np.abs(rotation.T @ metric @ rotation - metric)))
        if not math.isfinite(residual) or residual > threshold:
            _fail(
                "reference_metric_symmetry_mismatch",
                "a child rotation does not preserve the unstrained reference metric "
                f"within relative tolerance {tolerance:.3g}",
            )


def project_modes_to_free_rows(
    orbits: Sequence[ChildOrbitCoordinate],
    atoms: Sequence[ChildAtom],
    operations: Sequence[ExactSeitzOperation],
    modes: Sequence[RawModeColumn],
    *,
    frame_id: str,
    reference_lattice: Sequence[Sequence[object]],
    cartesian_tolerance_angstrom: float = 1.0e-8,
    rank_relative_tolerance: float = 1.0e-10,
    require_resolved_provenance: bool = True,
    require_complete: bool = True,
) -> ModeProjection:
    """Project resolved modes and canonicalize each free-coordinate column."""

    emitted = canonicalize_emitted_seitz(operations)
    lattice = np.asarray(_finite_lattice(reference_lattice), dtype=float)
    cart_tolerance = _positive_finite(cartesian_tolerance_angstrom, "Cartesian projection tolerance (Angstrom)")
    rank_tolerance = _positive_finite(rank_relative_tolerance, "relative rank tolerance")
    atom_ids = tuple(atom.atom_id for atom in atoms)
    if len(set(atom_ids)) != len(atom_ids):
        _fail("duplicate_atom_id", "child atom IDs must be unique")
    ordered_modes = tuple(sorted(modes, key=lambda mode: mode.column_index))
    if tuple(mode.column_index for mode in ordered_modes) != tuple(range(len(ordered_modes))):
        _fail("noncontiguous_mode_columns", "mode columns must be contiguous and zero-based")
    if len({mode.mode_id for mode in ordered_modes}) != len(ordered_modes):
        _fail("duplicate_mode_identity", "mode IDs must be unique")
    mode_maps: list[dict[str, np.ndarray]] = []
    for mode in ordered_modes:
        if mode.frame_id != frame_id:
            _fail("mode_frame_mismatch", f"mode {mode.mode_id!r} is not in frame {frame_id!r}")
        if set(mode.atom_ids) != set(atom_ids) or len(mode.atom_ids) != len(atom_ids):
            _fail("mode_atom_order_mismatch", f"mode {mode.mode_id!r} does not cover the child atom set")
        source_issue = mode.export_validation_issue()
        if require_resolved_provenance and source_issue is not None:
            _fail(
                "unresolved_mode_provenance",
                f"mode {mode.mode_id!r} is not export-ready: {source_issue}",
            )
        mode_maps.append(
            {
                atom_id: np.asarray(vector, dtype=float)
                for atom_id, vector in zip(mode.atom_ids, mode.displacements, strict=True)
            }
        )
    if require_resolved_provenance:
        try:
            validate_microscopic_domain_extension_bindings([
                (
                    mode.mode_identity,  # type: ignore[arg-type]
                    mode.microscopic_provenance,  # type: ignore[arg-type]
                    mode.microscopic_domain_extension,
                )
                for mode in ordered_modes
            ])
        except ValueError as exc:
            _fail(
                "invalid_domain_extension_evidence",
                str(exc),
            )

    free_rows = tuple(
        FreeCoordinateRow(
            orbit.representative_atom_id,
            local_index,
            parameter.axis,
            parameter.symbol,
            f"{orbit.representative_atom_id}_{parameter.suffix}",
        )
        for orbit in orbits
        for local_index, parameter in enumerate(orbit.site.parameters)
    )
    free_count = len(free_rows)
    mode_count = len(ordered_modes)
    if require_complete and mode_count != free_count:
        _fail(
            "incomplete_mode_basis",
            f"complete CIF matrix needs {free_count} modes, received {mode_count}",
        )
    matrix = np.zeros((free_count, mode_count), dtype=float)
    input_to_canonical: list[float] = []
    max_transporter = 0.0
    max_fixed = 0.0
    max_reconstruction = 0.0
    row_offsets: list[int] = []
    row_offset = 0
    for orbit in orbits:
        row_offsets.append(row_offset)
        row_offset += orbit.site.dimension

    for mode_index, (mode, displacement_by_atom) in enumerate(zip(ordered_modes, mode_maps, strict=True)):
        free_coefficients: list[np.ndarray] = []
        for orbit in orbits:
            site = orbit.site
            representative = displacement_by_atom[orbit.representative_atom_id]
            coefficients = (
                representative[np.asarray(site.pivot_axes, dtype=int)] if site.dimension else np.empty(0, dtype=float)
            )
            free_coefficients.append(coefficients)
        maximum = max((float(np.max(np.abs(values))) for values in free_coefficients if values.size), default=0.0)
        if not math.isfinite(maximum) or maximum <= 0:
            _fail(
                "zero_free_coordinate_mode",
                f"mode {mode.mode_id!r} has no nonzero child free-coordinate component",
            )
        scale_factor = 1.0 / maximum
        input_to_canonical.append(scale_factor)
        canonical_by_atom = {
            atom_id: vector * scale_factor for atom_id, vector in displacement_by_atom.items()
        }

        for orbit_index, orbit in enumerate(orbits):
            site = orbit.site
            basis = np.asarray(site.fixed_basis, dtype=float).reshape(3, site.dimension)
            representative = canonical_by_atom[orbit.representative_atom_id]
            coefficients = free_coefficients[orbit_index] * scale_factor
            if not orbit.transporters:
                _fail("missing_mode_transporter", f"mode {mode.mode_id!r} has an empty child orbit")
            for transporter in orbit.transporters:
                actual = canonical_by_atom[transporter.member_atom_id]
                for operation_index in transporter.candidate_operation_indices:
                    rotation = np.asarray(emitted[operation_index].rotation, dtype=float)
                    expected = representative @ rotation.T
                    residual = _cartesian_norm(actual - expected, lattice)
                    max_transporter = max(max_transporter, residual)
                    if residual > cart_tolerance:
                        _fail(
                            "conflicting_mode_transporters",
                            f"mode {mode.mode_id!r} gives conflicting symmetry arrows "
                            f"for orbit {orbit.representative_atom_id!r}",
                        )
            projected = basis @ coefficients
            fixed_residual = _cartesian_norm(representative - projected, lattice)
            max_fixed = max(max_fixed, fixed_residual)
            if fixed_residual > cart_tolerance:
                _fail(
                    "mode_outside_site_fixed_space",
                    f"mode {mode.mode_id!r} violates site constraints at {orbit.representative_atom_id!r}",
                )
            start = row_offsets[orbit_index]
            matrix[start : start + site.dimension, mode_index] = coefficients
            for transporter in orbit.transporters:
                actual = canonical_by_atom[transporter.member_atom_id]
                operation = emitted[transporter.canonical_operation_index]
                expected = projected @ np.asarray(operation.rotation, dtype=float).T
                residual = _cartesian_norm(actual - expected, lattice)
                max_reconstruction = max(max_reconstruction, residual)
                if residual > cart_tolerance:
                    _fail(
                        "mode_reconstruction_failure",
                        f"mode {mode.mode_id!r} does not reconstruct atom {transporter.member_atom_id!r}",
                    )

    if matrix.size:
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        largest = float(singular_values[0]) if singular_values.size else 0.0
        threshold = rank_tolerance * largest
        rank = int(np.count_nonzero(singular_values > threshold))
        minimum = float(singular_values[-1]) if singular_values.size else 0.0
    else:
        rank = 0
        minimum = 0.0
    if require_complete and rank != free_count:
        _fail(
            "rank_deficient_mode_matrix",
            f"complete {free_count}x{mode_count} mode matrix has rank {rank}",
        )
    return ModeProjection(
        free_rows=free_rows,
        mode_ids=tuple(mode.mode_id for mode in ordered_modes),
        matrix=tuple(tuple(float(value) for value in row) for row in matrix),
        input_to_canonical_scale_factors=tuple(input_to_canonical),
        diagnostics=ModeProjectionDiagnostics(
            free_coordinate_count=free_count,
            mode_count=mode_count,
            matrix_rank=rank,
            minimum_singular_value=minimum,
            max_transporter_residual_angstrom=max_transporter,
            max_fixed_space_residual_angstrom=max_fixed,
            max_reconstruction_residual_angstrom=max_reconstruction,
        ),
    )


@dataclass(frozen=True)
class CenteringInfo:
    translations: tuple[RationalVector, ...]
    index: int
    rotation_coset_count: int
    determinant_index: int | None


def derive_centering_index(
    operations: Sequence[ExactSeitzOperation],
    *,
    conventional_basis: Sequence[Sequence[object]] | None = None,
    primitive_translation_basis: Sequence[Sequence[object]] | None = None,
) -> CenteringInfo:
    """Derive ``C`` from pure translations and optionally verify exact volumes."""

    emitted = canonicalize_emitted_seitz(operations)
    identity = identity_matrix()
    translations = tuple(
        sorted(
            (operation.translation for operation in emitted if operation.rotation == identity),
            key=_vector_key,
        )
    )
    if not translations or (Fraction(0), Fraction(0), Fraction(0)) not in translations:
        _fail("invalid_centering_subgroup", "pure-translation subgroup has no identity")
    translation_set = set(translations)
    for left in translations:
        for right in translations:
            if _mod_one(_add(left, right)) not in translation_set:
                _fail("invalid_centering_subgroup", "pure translations are not closed modulo 1")
    index = len(translations)
    counts: dict[object, int] = {}
    for operation in emitted:
        key = _rotation_key(operation.rotation)
        counts[key] = counts.get(key, 0) + 1
    if any(count != index for count in counts.values()):
        _fail("inconsistent_centering_cosets", "rotation cosets do not all contain C translations")

    determinant_index: int | None = None
    if (conventional_basis is None) != (primitive_translation_basis is None):
        _fail(
            "incomplete_centering_lattices",
            "conventional and primitive bases must be supplied together",
        )
    if conventional_basis is not None and primitive_translation_basis is not None:
        conventional = rational_matrix(conventional_basis)
        primitive = rational_matrix(primitive_translation_basis)
        ratio = abs(determinant(conventional) / determinant(primitive))
        if ratio.denominator != 1 or ratio <= 0:
            _fail("invalid_centering_determinant", "conventional/primitive volume ratio is not a positive integer")
        determinant_index = int(ratio)
        if determinant_index != index:
            _fail(
                "centering_index_mismatch",
                f"Seitz translations give C={index}, lattice determinants give C={determinant_index}",
            )
    return CenteringInfo(translations, index, len(counts), determinant_index)


@dataclass(frozen=True)
class NormDiagnostics:
    conventional_rss: tuple[float, ...]
    primitive_rss: tuple[float, ...]
    max_unit_rss_error: float


@dataclass(frozen=True)
class DisplaciveNorms:
    mode_ids: tuple[str, ...]
    normfactors: tuple[float, ...]
    diagnostics: NormDiagnostics


def primitive_rss_normfactors(
    modes: Sequence[RawModeColumn],
    reference_lattice: Sequence[Sequence[object]],
    centering_index: int,
    *,
    rss_tolerance: float = 1.0e-12,
) -> DisplaciveNorms:
    """Normalize raw conventional-cell modes using the unstrained child metric."""

    lattice = np.asarray(_finite_lattice(reference_lattice), dtype=float)
    if isinstance(centering_index, bool) or not isinstance(centering_index, Integral):
        _fail("invalid_centering_index", "centering index must be a positive integer")
    if int(centering_index) <= 0:
        _fail("invalid_centering_index", "centering index must be positive")
    tolerance = _positive_finite(rss_tolerance, "RSS normalization tolerance")
    ordered = tuple(sorted(modes, key=lambda mode: mode.column_index))
    conventional_values: list[float] = []
    primitive_values: list[float] = []
    normfactors: list[float] = []
    max_error = 0.0
    for mode in ordered:
        array = np.asarray(mode.displacements, dtype=float)
        if array.ndim != 2 or array.shape[1] != 3 or not np.all(np.isfinite(array)):
            _fail("invalid_mode_displacements", f"mode {mode.mode_id!r} is not one finite Nx3 array")
        cartesian = array @ lattice
        conventional = float(np.sum(cartesian * cartesian))
        primitive = conventional / int(centering_index)
        if not math.isfinite(primitive) or primitive <= 0:
            _fail("zero_mode_norm", f"mode {mode.mode_id!r} has zero or invalid primitive RSS")
        normfactor = 1.0 / math.sqrt(primitive)
        unit_rss = primitive * normfactor * normfactor
        error = abs(unit_rss - 1.0)
        if error > tolerance:
            _fail("mode_norm_roundtrip_failure", f"mode {mode.mode_id!r} failed primitive RSS normalization")
        conventional_values.append(conventional)
        primitive_values.append(primitive)
        normfactors.append(normfactor)
        max_error = max(max_error, error)
    return DisplaciveNorms(
        mode_ids=tuple(mode.mode_id for mode in ordered),
        normfactors=tuple(normfactors),
        diagnostics=NormDiagnostics(
            conventional_rss=tuple(conventional_values),
            primitive_rss=tuple(primitive_values),
            max_unit_rss_error=max_error,
        ),
    )


@dataclass(frozen=True)
class ExactChildFrame:
    frame_id: str
    operations: tuple[ExactSeitzOperation, ...]
    reference_lattice: tuple[FloatVector, FloatVector, FloatVector]
    conventional_basis: RationalMatrix | None
    primitive_translation_basis: RationalMatrix | None
    metric_relative_tolerance: float

    @classmethod
    def from_values(
        cls,
        frame_id: str,
        operations: Sequence[ExactSeitzOperation],
        reference_lattice: Sequence[Sequence[object]],
        *,
        conventional_basis: Sequence[Sequence[object]] | None = None,
        primitive_translation_basis: Sequence[Sequence[object]] | None = None,
        metric_relative_tolerance: float = 1.0e-8,
    ) -> ExactChildFrame:
        identifier = _identifier(frame_id, "child frame ID", "invalid_frame_id")
        if (conventional_basis is None) != (primitive_translation_basis is None):
            _fail(
                "incomplete_centering_lattices",
                "conventional and primitive bases must be supplied together",
            )
        conventional = None if conventional_basis is None else rational_matrix(conventional_basis)
        primitive = None if primitive_translation_basis is None else rational_matrix(primitive_translation_basis)
        emitted = canonicalize_emitted_seitz(operations)
        lattice = _finite_lattice(reference_lattice)
        metric_tolerance = _positive_finite(
            metric_relative_tolerance, "reference metric relative tolerance"
        )
        _validate_metric_compatibility(emitted, lattice, metric_tolerance)
        return cls(
            identifier,
            emitted,
            lattice,
            conventional,
            primitive,
            metric_tolerance,
        )


@dataclass(frozen=True)
class DisplaciveModelDiagnostics:
    status: str
    operation_count: int
    atom_count: int
    orbit_count: int
    free_coordinate_count: int
    mode_count: int
    centering_index: int
    matrix_rank: int
    max_cartesian_residual_angstrom: float
    max_unit_rss_error: float
    validation_issues: tuple[str, ...] = ()


@dataclass(frozen=True)
class DisplaciveCifModel:
    frame: ExactChildFrame
    atoms: tuple[ChildAtom, ...]
    orbits: tuple[ChildOrbitCoordinate, ...]
    projection: ModeProjection
    canonical_modes: tuple[RawModeColumn, ...]
    centering: CenteringInfo
    norms: DisplaciveNorms
    amplitudes: tuple[float, ...]
    amplitude_convention: str
    diagnostics: DisplaciveModelDiagnostics


def build_displacive_cif_model(
    frame: ExactChildFrame,
    atoms: Sequence[ChildAtom],
    representative_atom_ids: Sequence[str],
    modes: Sequence[RawModeColumn],
    *,
    amplitudes: Sequence[Real] | None = None,
    amplitude_convention: str = _AS_AMPLITUDE_CONVENTION,
    cartesian_tolerance_angstrom: float = 1.0e-8,
    rank_relative_tolerance: float = 1.0e-10,
    require_resolved_provenance: bool = True,
    require_complete: bool = True,
) -> DisplaciveCifModel:
    """Build the complete validated source for CIF displacive-mode loops."""

    emitted = canonicalize_emitted_seitz(frame.operations)
    parsed_amplitude_convention = str(amplitude_convention).strip()
    if parsed_amplitude_convention != _AS_AMPLITUDE_CONVENTION:
        _fail(
            "unsupported_mode_amplitude_convention",
            "mode amplitudes must be As values in Angstrom on primitive-cell unit modes; "
            "raw mode coefficients are not accepted",
        )
    atom_values = tuple(atoms)
    input_modes = tuple(modes)
    if amplitudes is None:
        input_amplitudes = (0.0,) * len(input_modes)
    else:
        try:
            input_amplitudes = tuple(float(value) for value in amplitudes)
        except (TypeError, ValueError, OverflowError) as exc:
            raise DisplaciveModelError("invalid_mode_amplitudes", "mode amplitudes must be finite numbers") from exc
        if len(input_amplitudes) != len(input_modes) or not all(math.isfinite(value) for value in input_amplitudes):
            _fail(
                "invalid_mode_amplitudes",
                "mode amplitudes must be finite and match the mode count",
            )
    ordered_pairs = tuple(
        sorted(
            zip(input_modes, input_amplitudes, strict=True),
            key=lambda item: item[0].column_index,
        )
    )
    mode_values = tuple(item[0] for item in ordered_pairs)
    amplitude_values = tuple(item[1] for item in ordered_pairs)
    orbits = build_child_orbits(atom_values, representative_atom_ids, emitted)
    centering = derive_centering_index(
        emitted,
        conventional_basis=frame.conventional_basis,
        primitive_translation_basis=frame.primitive_translation_basis,
    )
    if any(len(orbit.member_atom_ids) % centering.index != 0 for orbit in orbits):
        _fail(
            "orbit_centering_mismatch",
            "a conventional-child orbit multiplicity is not divisible by C",
        )
    projection = project_modes_to_free_rows(
        orbits,
        atom_values,
        emitted,
        mode_values,
        frame_id=frame.frame_id,
        reference_lattice=frame.reference_lattice,
        cartesian_tolerance_angstrom=cartesian_tolerance_angstrom,
        rank_relative_tolerance=rank_relative_tolerance,
        require_resolved_provenance=require_resolved_provenance,
        require_complete=require_complete,
    )
    canonical_atom_ids = tuple(atom.atom_id for atom in atom_values)
    canonical_modes_list: list[RawModeColumn] = []
    for mode, scale_factor in zip(
        mode_values, projection.input_to_canonical_scale_factors, strict=True
    ):
        displacement_by_atom = dict(
            zip(mode.atom_ids, mode.displacements, strict=True)
        )
        canonical_displacements = np.asarray(
            [displacement_by_atom[atom_id] for atom_id in canonical_atom_ids],
            dtype=float,
        ) * scale_factor
        canonical_provenance = mode.microscopic_provenance
        if canonical_provenance is not None and mode.source_binding_issue() is None:
            canonical_provenance = canonical_provenance.bind_mapping(
                frame_id=mode.frame_id,
                atom_ids=canonical_atom_ids,
            )
        canonical_modes_list.append(
            RawModeColumn.from_values(
                mode.mode_id,
                mode.column_index,
                mode.frame_id,
                canonical_atom_ids,
                canonical_displacements,
                ModeScaleProvenance(
                    source=mode.scale_provenance.source,
                    convention="free_coordinate_max_abs_one",
                    direction_resolved=mode.scale_provenance.direction_resolved,
                    sign_resolved=mode.scale_provenance.sign_resolved,
                ),
                mode_identity=mode.mode_identity,
                microscopic_provenance=canonical_provenance,
                microscopic_domain_extension=mode.microscopic_domain_extension,
                label=mode.label,
            )
        )
    canonical_modes = tuple(canonical_modes_list)
    norms = primitive_rss_normfactors(canonical_modes, frame.reference_lattice, centering.index)
    if projection.mode_ids != norms.mode_ids:
        _fail("mode_order_mismatch", "projection and normalization use different mode orders")
    projection_diag = projection.diagnostics
    max_residual = max(
        projection_diag.max_transporter_residual_angstrom,
        projection_diag.max_fixed_space_residual_angstrom,
        projection_diag.max_reconstruction_residual_angstrom,
    )
    validation_issues = tuple(
        [
            issue
            for mode in canonical_modes
            if (source_issue := mode.export_validation_issue()) is not None
            for issue in (f"mode {mode.mode_id!r}: {source_issue}",)
        ]
        + (
            []
            if projection_diag.mode_count == projection_diag.free_coordinate_count
            else [
                "mode basis is incomplete: "
                f"{projection_diag.mode_count} columns for "
                f"{projection_diag.free_coordinate_count} free coordinates"
            ]
        )
        + (
            []
            if projection_diag.matrix_rank == projection_diag.free_coordinate_count
            else [
                "mode matrix is not full rank: "
                f"rank {projection_diag.matrix_rank} for "
                f"{projection_diag.free_coordinate_count} free coordinates"
            ]
        )
    )
    return DisplaciveCifModel(
        frame=ExactChildFrame(
            frame.frame_id,
            emitted,
            frame.reference_lattice,
            frame.conventional_basis,
            frame.primitive_translation_basis,
            frame.metric_relative_tolerance,
        ),
        atoms=atom_values,
        orbits=orbits,
        projection=projection,
        canonical_modes=canonical_modes,
        centering=centering,
        norms=norms,
        amplitudes=amplitude_values,
        amplitude_convention=parsed_amplitude_convention,
        diagnostics=DisplaciveModelDiagnostics(
            status="validated" if not validation_issues else "diagnostic",
            operation_count=len(emitted),
            atom_count=len(atom_values),
            orbit_count=len(orbits),
            free_coordinate_count=projection_diag.free_coordinate_count,
            mode_count=projection_diag.mode_count,
            centering_index=centering.index,
            matrix_rank=projection_diag.matrix_rank,
            max_cartesian_residual_angstrom=max_residual,
            max_unit_rss_error=norms.diagnostics.max_unit_rss_error,
            validation_issues=validation_issues,
        ),
    )


__all__ = [
    "AtomTransporter",
    "CenteringInfo",
    "ChildAtom",
    "ChildOrbitCoordinate",
    "DeltaParameter",
    "DisplaciveCifModel",
    "DisplaciveModelDiagnostics",
    "DisplaciveModelError",
    "DisplaciveNorms",
    "ExactChildFrame",
    "ExactSeitzOperation",
    "FreeCoordinateRow",
    "ModeProjection",
    "ModeProjectionDiagnostics",
    "ModeScaleProvenance",
    "NormDiagnostics",
    "RawModeColumn",
    "SiteCoordinateModel",
    "build_child_orbits",
    "build_displacive_cif_model",
    "canonicalize_emitted_seitz",
    "conjugate_seitz_operations",
    "derive_centering_index",
    "primitive_rss_normfactors",
    "project_modes_to_free_rows",
    "site_fixed_coordinate_model",
    "transform_fractional_displacement",
    "transform_fractional_position",
    "transformed_reference_lattice",
]
