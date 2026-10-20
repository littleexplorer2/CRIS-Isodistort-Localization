"""Ordinary 3D structures. Mathematical contract: docs/reference_library/STRUCTURE_CONTRACT.md."""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
import math
from typing import Iterable

Vector = tuple[Fraction, Fraction, Fraction]
IntegerMatrix = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
IDENTITY: IntegerMatrix = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
METRIC_RELATIVE_TOLERANCE = 1e-8
MAX_OPERATIONS = 384


class StructureError(ValueError):
    """A structured rejection, never a silently repaired scientific input."""

    def __init__(self, code: str, message: str, **details: object):
        super().__init__(message)
        self.code = code
        self.details = details


def _fraction(value: object) -> Fraction:
    try:
        return value if isinstance(value, Fraction) else Fraction(str(value))
    except (ValueError, ZeroDivisionError, TypeError) as exc:
        raise StructureError("invalid_number", "A finite rational value is required.") from exc


def _vector(values: Iterable[object]) -> Vector:
    result = tuple(_fraction(value) for value in values)
    if len(result) != 3:
        raise StructureError("invalid_dimension", "Exactly three coordinates are required.")
    return result  # type: ignore[return-value]


def _determinant(r: IntegerMatrix) -> int:
    return (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))


@dataclass(frozen=True)
class Lattice:
    """Lengths in angstrom, angles in degrees; matrix columns are a, b, c."""

    a: Fraction
    b: Fraction
    c: Fraction
    alpha: Fraction
    beta: Fraction
    gamma: Fraction
    standard_uncertainties: tuple[Fraction | None, ...] = (None,) * 6

    def __post_init__(self) -> None:
        for name in ("a", "b", "c", "alpha", "beta", "gamma"):
            value = _fraction(getattr(self, name))
            try:
                finite = math.isfinite(float(value))
            except OverflowError:
                finite = False
            if not finite:
                raise StructureError("invalid_lattice", "Lattice values must be finite.")
            object.__setattr__(self, name, value)
        if min(self.a, self.b, self.c) <= 0:
            raise StructureError("invalid_lattice", "All lattice lengths must be positive.")
        if any(not 0 < v < 180 for v in (self.alpha, self.beta, self.gamma)):
            raise StructureError("invalid_lattice", "Angles must lie strictly between 0 and 180 degrees.")
        if len(self.standard_uncertainties) != 6:
            raise StructureError("invalid_dimension", "Six lattice uncertainties are required.")
        # Reject non-positive and numerically singular Gram matrices before use.
        matrix = self.matrix
        volume = matrix[0][0] * matrix[1][1] * matrix[2][2]
        if not math.isfinite(volume) or volume <= 0:
            raise StructureError("invalid_lattice", "Cell volume overflows or underflows binary64.")
        object.__setattr__(self, "standard_uncertainties", tuple(self.standard_uncertainties))

    @property
    def matrix(self) -> tuple[tuple[float, float, float], ...]:
        a, b, c = map(float, (self.a, self.b, self.c))
        ca, cb, cg = (math.cos(math.radians(float(x)))
                      for x in (self.alpha, self.beta, self.gamma))
        sg = math.sin(math.radians(float(self.gamma)))
        if abs(sg) <= 1e-12:
            raise StructureError("singular_lattice", "The lattice is numerically singular.")
        cy = (ca - cb * cg) / sg
        cz2 = 1 - cb * cb - cy * cy
        if cz2 <= 1e-12:
            raise StructureError("singular_lattice", "The cell Gram matrix is non-positive or numerically singular.")
        matrix = ((a, b * cg, c * cb), (0.0, b * sg, c * cy),
                  (0.0, 0.0, c * math.sqrt(cz2)))
        if not all(math.isfinite(v) for row in matrix for v in row):
            raise StructureError("invalid_lattice", "Cartesian matrix overflow.")
        return matrix

    @property
    def metric(self) -> tuple[tuple[float, float, float], ...]:
        a = self.matrix
        return tuple(tuple(sum(a[k][i] * a[k][j] for k in range(3))
                           for j in range(3)) for i in range(3))

    @property
    def volume(self) -> float:
        a = self.matrix
        return a[0][0] * a[1][1] * a[2][2]

    def fractional_to_cartesian(self, coordinates: Iterable[object]) -> tuple[float, float, float]:
        x = tuple(map(float, _vector(coordinates)))
        a = self.matrix
        return tuple(sum(a[i][j] * x[j] for j in range(3)) for i in range(3))

    def cartesian_to_fractional(self, coordinates: Iterable[object]) -> tuple[float, float, float]:
        x = tuple(map(float, _vector(coordinates)))
        a = self.matrix
        z = x[2] / a[2][2]
        y = (x[1] - a[1][2] * z) / a[1][1]
        return ((x[0] - a[0][1] * y - a[0][2] * z) / a[0][0], y, z)


@dataclass(frozen=True)
class SymmetryOperation:
    """Active fractional action x' = R x + t, with t reduced modulo integers."""

    rotation: IntegerMatrix
    translation: Vector
    expression: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        if len(self.rotation) != 3 or any(len(row) != 3 for row in self.rotation):
            raise StructureError("invalid_symmetry_operation", "A 3 by 3 rotation is required.")
        entries = tuple(tuple(_fraction(v) for v in row) for row in self.rotation)
        if any(v.denominator != 1 for row in entries for v in row):
            raise StructureError("nonintegral_rotation", "Fractional-coordinate lattice rotations must be integral.")
        r = tuple(tuple(int(v) for v in row) for row in entries)
        if _determinant(r) not in (-1, 1):
            raise StructureError("nonunimodular_rotation", "Rotation determinant must be +1 or -1.")
        object.__setattr__(self, "rotation", r)
        object.__setattr__(self, "translation", tuple(v % 1 for v in _vector(self.translation)))

    def apply(self, coordinates: Iterable[object], *, wrap: bool = True) -> Vector:
        x = _vector(coordinates)
        result = tuple(sum(self.rotation[i][j] * x[j] for j in range(3))
                       + self.translation[i] for i in range(3))
        return tuple(v % 1 for v in result) if wrap else result  # type: ignore[return-value]

    def compose(self, other: SymmetryOperation) -> SymmetryOperation:
        """self after other; translations are equivalent modulo the given cell."""
        r = tuple(tuple(sum(self.rotation[i][k] * other.rotation[k][j] for k in range(3))
                        for j in range(3)) for i in range(3))
        t = self.apply(other.translation, wrap=True)
        return SymmetryOperation(r, t)


@dataclass(frozen=True)
class AtomicSite:
    label: str
    type_symbol: str
    fractional: Vector
    occupancy: Fraction = Fraction(1)
    wyckoff_symbol: str | None = None
    site_symmetry_multiplicity: int | None = None
    coordinate_uncertainties: tuple[Fraction | None, ...] = (None,) * 3
    occupancy_uncertainty: Fraction | None = None
    multiplicity_source_tag: str | None = None

    def __post_init__(self) -> None:
        if not self.label or not self.type_symbol or any(
                c.isspace() for c in self.label + self.type_symbol):
            raise StructureError("invalid_atom_identity", "Atom label and explicit type must be non-empty single tokens.")
        if any(c in self.type_symbol for c in ("/", ",", ";")):
            raise StructureError("mixed_occupancy_unsupported", "Composite atom-type expressions are unsupported.")
        object.__setattr__(self, "fractional", _vector(self.fractional))
        object.__setattr__(self, "occupancy", _fraction(self.occupancy))
        if not 0 < self.occupancy <= 1:
            raise StructureError("invalid_occupancy", "Occupancy must be in (0, 1].")
        if self.site_symmetry_multiplicity is not None and (
                not isinstance(self.site_symmetry_multiplicity, int)
                or self.site_symmetry_multiplicity < 1):
            raise StructureError("invalid_multiplicity", "Declared multiplicity must be a positive integer.")
        uncertainties = tuple(self.coordinate_uncertainties)
        if len(uncertainties) != 3:
            raise StructureError("invalid_dimension", "Three coordinate uncertainties are required.")
        object.__setattr__(self, "coordinate_uncertainties", uncertainties)

    @property
    def symbol(self) -> str:
        """Explicit species/type identity; it is not inferred from the site label."""
        return self.type_symbol


@dataclass(frozen=True)
class CrystalStructure:
    data_name: str
    lattice: Lattice
    sites: tuple[AtomicSite, ...]
    symmetry_operations: tuple[SymmetryOperation, ...]
    space_group_number: int | None = None
    space_group_name: str | None = None
    source_sha256: str | None = None
    notices: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "sites", tuple(self.sites))
        object.__setattr__(self, "symmetry_operations", tuple(self.symmetry_operations))
        object.__setattr__(self, "notices", tuple(self.notices))


@dataclass(frozen=True)
class ExpandedSite:
    """One distinct cell site; source_index and operation_indices are zero based."""

    label: str
    type_symbol: str
    fractional: Vector
    occupancy: Fraction
    source_label: str
    source_index: int
    operation_indices: tuple[int, ...]

    @property
    def symbol(self) -> str:
        return self.type_symbol

    @property
    def asymmetric_index(self) -> int:
        return self.source_index


def validate_symmetry_operations(
    operations: Iterable[SymmetryOperation], lattice: Lattice | None = None,
    *, metric_relative_tolerance: float = METRIC_RELATIVE_TOLERANCE,
) -> tuple[SymmetryOperation, ...]:
    ops = tuple(operations)
    if not ops:
        raise StructureError("missing_symmetry_operations", "Explicit space-group operations are required.")
    if len(ops) > MAX_OPERATIONS:
        raise StructureError("operation_limit", f"At most {MAX_OPERATIONS} explicit operations are supported.")
    if len(set(ops)) != len(ops):
        raise StructureError("duplicate_symmetry_operation", "Operations repeat modulo integer translations.")
    group = set(ops)
    if SymmetryOperation(IDENTITY, (Fraction(0),) * 3) not in group:
        raise StructureError("missing_identity", "The explicit operation list lacks the identity.")
    for i, g in enumerate(ops):
        for j, h in enumerate(ops):
            if g.compose(h) not in group:
                raise StructureError(
                    "nonclosed_symmetry_operations",
                    "The explicit operations are not exactly closed modulo integer translations. "
                    "Rounded decimal translations such as 0.333333 are not silently changed to 1/3.",
                    left_index=i, right_index=j)
    if not math.isfinite(metric_relative_tolerance) or not 0 <= metric_relative_tolerance <= 1e-4:
        raise StructureError("invalid_tolerance", "Metric relative tolerance must be finite and in [0, 1e-4].")
    if lattice is not None:
        metric = lattice.metric
        if (not all(math.isfinite(v) for row in metric for v in row)
                or any(metric[i][i] <= 0 for i in range(3))):
            raise StructureError("invalid_lattice", "The cell metric cannot be represented in binary64.")
        lengths = tuple(math.sqrt(metric[i][i]) for i in range(3))
        for index, op in enumerate(ops):
            transformed = tuple(tuple(sum(op.rotation[k][i] * metric[k][l] * op.rotation[l][j]
                                          for k in range(3) for l in range(3))
                                      for j in range(3)) for i in range(3))
            # Normalize each Gram entry by its own pair of axis lengths: a long c
            # must not hide a broken a/b isometry. See STRUCTURE_CONTRACT.md.
            residual = max(abs(transformed[i][j] - metric[i][j]) / (lengths[i] * lengths[j])
                           for i in range(3) for j in range(3))
            if residual > metric_relative_tolerance:
                raise StructureError("symmetry_metric_mismatch", "An operation does not preserve the supplied cell metric.",
                                     operation_index=index, relative_residual=residual,
                                     tolerance=metric_relative_tolerance)
    return ops


def expand_structure(
    structure: CrystalStructure, *, metric_relative_tolerance: float = METRIC_RELATIVE_TOLERANCE,
) -> tuple[ExpandedSite, ...]:
    """Expand independent representatives exactly. No distance snapping or averaging."""
    ops = validate_symmetry_operations(structure.symmetry_operations, structure.lattice,
                                       metric_relative_tolerance=metric_relative_tolerance)
    if not structure.sites:
        raise StructureError("missing_atoms", "At least one atomic site is required.")
    if len({site.label for site in structure.sites}) != len(structure.sites):
        raise StructureError("duplicate_atom_label", "Independent atom labels must be unique.")
    occupied: dict[Vector, AtomicSite] = {}
    expanded: list[ExpandedSite] = []
    for source_index, site in enumerate(structure.sites):
        orbit: dict[Vector, list[int]] = {}
        for index, operation in enumerate(ops):
            orbit.setdefault(operation.apply(site.fractional), []).append(index)
        multiplicity = site.site_symmetry_multiplicity
        if multiplicity is not None and multiplicity != len(orbit):
            message = "Declared multiplicity differs from the exact orbit size."
            if site.multiplicity_source_tag == "_atom_site_symmetry_multiplicity":
                message += " The legacy CIF tag is historically ambiguous; confirm whether it meant orbit size or site-group order."
            raise StructureError("multiplicity_mismatch", message, label=site.label,
                                 declared=multiplicity, computed=len(orbit),
                                 source_tag=site.multiplicity_source_tag)
        for orbit_index, (position, generating_ops) in enumerate(orbit.items()):
            previous = occupied.get(position)
            if previous is not None:
                code = ("mixed_occupancy_unsupported" if previous.type_symbol != site.type_symbol
                        else "overlapping_independent_sites")
                raise StructureError(code, "Independent sites generate an identical fractional position.",
                                     first_label=previous.label, second_label=site.label)
            occupied[position] = site
            expanded.append(ExpandedSite(
                f"{site.label}:{orbit_index + 1}", site.type_symbol, position, site.occupancy,
                site.label, source_index, tuple(generating_ops)))
    return tuple(expanded)
