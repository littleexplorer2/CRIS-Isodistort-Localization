"""Validated, shared contract for all displacive-mode export writers.

The low-level CIF model proves the child setting, atom order, orbit
transporters, mode directions and normalization.  This module binds that
scientific model to two distinct structures:

* ``reference_structure`` is the undistorted child used to define every mode;
* ``final_structure`` is reconstructed from the reference plus the declared
  ``As`` amplitudes.

Writers must consume this object rather than interpreting an unlabelled
``dict[str, ndarray]`` independently.  In particular, a final structure is
never reused as the coordinate-formula reference, which would apply the same
displacement twice.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from numbers import Integral
from types import MappingProxyType

import numpy as np
import spglib
from pymatgen.core import IStructure, Lattice, Structure

from ..utils.lattice import (
    RationalMatrix,
    as_fraction,
    centering_primitive_matrix,
    determinant,
    inverse,
    is_integral,
    matrix_key,
    multiply,
    rational_matrix,
)
from ..utils.opd_format import _centering_letter
from .cif_displacive_model import (
    DisplaciveCifModel,
    DisplaciveModelError,
    ExactSeitzOperation,
)


def _fail(code: str, message: str) -> None:
    raise DisplaciveModelError(code, message)


def _periodic_cartesian_residual(
    first: np.ndarray,
    second: np.ndarray,
    lattice: np.ndarray,
) -> float:
    first_array = np.asarray(first, dtype=float)
    second_array = np.asarray(second, dtype=float)
    if first_array.shape != second_array.shape:
        raise ValueError("periodic coordinate arrays must have matching shapes")
    if first_array.size == 0:
        return 0.0
    periodic_lattice = Lattice(np.asarray(lattice, dtype=float))
    return max(
        float(
            periodic_lattice.get_distance_and_image(
                first_coordinate,
                second_coordinate,
            )[0]
        )
        for first_coordinate, second_coordinate in zip(
            first_array,
            second_array,
            strict=True,
        )
    )


def _chemical_signature(site) -> tuple[tuple[str, Fraction], ...]:
    return tuple(
        sorted(
            (str(species), Fraction(str(float(occupancy))))
            for species, occupancy in site.species.items()
        )
    )


def _snapshot_structure(structure: Structure | IStructure) -> IStructure:
    """Detach an immutable crystallographic snapshot from caller-owned data."""

    if not isinstance(structure, (Structure, IStructure)):
        _fail("invalid_structure_snapshot", "export structures must be pymatgen structures")
    return IStructure(
        lattice=structure.lattice,
        species=[site.species for site in structure],
        coords=np.asarray(structure.frac_coords, dtype=float),
        charge=structure.charge,
        coords_are_cartesian=False,
        site_properties=deepcopy(structure.site_properties),
        labels=[site.label for site in structure],
        properties=deepcopy(structure.properties),
    )


def _require_full_single_species(structure: IStructure, description: str) -> None:
    for index, site in enumerate(structure):
        species = tuple(site.species.items())
        if len(species) != 1 or Fraction(str(float(species[0][1]))) != 1:
            _fail(
                "unsupported_site_occupancy",
                f"{description} site {index} must contain one fully occupied species",
            )


def _exact_origin(values) -> tuple[Fraction, Fraction, Fraction]:
    if len(values) != 3:
        _fail("invalid_parent_child_origin", "parent-to-child origin must contain three values")
    try:
        parsed = tuple(as_fraction(value) for value in values)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise DisplaciveModelError(
            "invalid_parent_child_origin",
            "parent-to-child origin must be an exact finite rational vector",
        ) from exc
    return parsed  # type: ignore[return-value]


def _exact_integer_translation(values) -> tuple[int, int, int]:
    if len(values) != 3:
        _fail(
            "invalid_parent_child_atom_mapping",
            "parent-cell translation must contain three integer values",
        )
    try:
        parsed = tuple(as_fraction(value) for value in values)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise DisplaciveModelError(
            "invalid_parent_child_atom_mapping",
            "parent-cell translation must contain three exact integer values",
        ) from exc
    if any(value.denominator != 1 for value in parsed):
        _fail(
            "invalid_parent_child_atom_mapping",
            "parent-cell translation must contain three exact integer values",
        )
    return tuple(int(value) for value in parsed)  # type: ignore[return-value]


def _affine_embedding_id(
    parent_space_group_number: int,
    basis: RationalMatrix,
    origin: tuple[Fraction, Fraction, Fraction],
) -> str:
    payload = (
        int(parent_space_group_number),
        matrix_key(basis),
        tuple((value.numerator, value.denominator) for value in origin),
    )
    digest = sha256(repr(payload).encode("ascii")).hexdigest()
    return f"parent-child-affine-v1:{digest}"


def _centering_multiplicity_from_letter(letter: str) -> int:
    letter = str(letter).strip().upper()
    if letter == "F":
        return 4
    if letter == "R":
        return 3
    if letter in {"A", "B", "C", "I"}:
        return 2
    return 1


def _centering_multiplicity(space_group_number: int) -> int:
    return _centering_multiplicity_from_letter(
        _centering_letter(int(space_group_number))
    )


def _database_seitz_operations(
    hall_number: int,
) -> tuple[ExactSeitzOperation, ...]:
    database = spglib.get_symmetry_from_database(int(hall_number))
    if database is None:
        _fail(
            "invalid_target_hall_setting",
            f"Hall number {hall_number} has no spglib symmetry database entry",
        )
    return tuple(
        ExactSeitzOperation.from_values(rotation, translation)
        for rotation, translation in zip(
            database["rotations"],
            database["translations"],
            strict=True,
        )
    )


def _hall_centering_multiplicity(hall_number: int) -> int:
    operations = _database_seitz_operations(hall_number)
    identity = rational_matrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    return sum(operation.rotation == identity for operation in operations)


def _validate_emitted_space_group(
    operations: tuple[ExactSeitzOperation, ...],
    *,
    model_centering_index: int,
    target_space_group_number: int,
    target_hall_number: int | None,
) -> None:
    """Bind the exact emitted Seitz group to the declared target identity."""

    expected_centering = (
        _centering_multiplicity(target_space_group_number)
        if target_hall_number is None
        else _hall_centering_multiplicity(target_hall_number)
    )
    if int(model_centering_index) != expected_centering:
        _fail(
            "target_centering_mismatch",
            "emitted Seitz translations give centering multiplicity "
            f"{model_centering_index}, but the target space-group setting requires "
            f"{expected_centering}",
        )

    rotations = np.asarray(
        [operation.rotation for operation in operations],
        dtype=np.intc,
    )
    translations = np.asarray(
        [operation.translation for operation in operations],
        dtype=float,
    )
    identified = spglib.get_spacegroup_type_from_symmetry(
        rotations,
        translations,
        symprec=1.0e-8,
    )
    if identified is None:
        _fail(
            "unrecognized_emitted_space_group",
            "spglib could not identify the complete emitted exact Seitz group",
        )
    if int(identified.number) != int(target_space_group_number):
        _fail(
            "target_space_group_type_mismatch",
            "emitted exact Seitz operations identify space group "
            f"{identified.number}, not declared target {target_space_group_number}",
        )

    if target_hall_number is not None:
        expected_keys = {
            operation.key
            for operation in _database_seitz_operations(target_hall_number)
        }
        emitted_keys = {operation.key for operation in operations}
        if emitted_keys != expected_keys:
            _fail(
                "target_space_group_setting_mismatch",
                "emitted exact Seitz operations do not match the declared Hall setting "
                f"{target_hall_number}",
            )


@dataclass(frozen=True)
class ParentChildSiteMapping:
    """One explicit parent-site image in the emitted reference child cell.

    ``parent_cell_translation`` is an integer translation of the conventional
    parent cell.  Together with ``parent_site_index`` it fixes the unwrapped
    parent coordinate corresponding to ``child_atom_id``; no nearest-site or
    same-element guess is permitted at export time.
    """

    child_atom_id: str
    parent_site_index: int
    parent_orbit_id: str
    parent_cell_translation: tuple[int, int, int]

    def __post_init__(self) -> None:
        child_atom_id = str(self.child_atom_id).strip()
        parent_orbit_id = str(self.parent_orbit_id).strip()
        parent_site_index = self.parent_site_index
        if not child_atom_id or not parent_orbit_id:
            _fail(
                "invalid_parent_child_atom_mapping",
                "parent-child mapping requires nonempty child-atom and parent-orbit IDs",
            )
        if (
            isinstance(parent_site_index, bool)
            or not isinstance(parent_site_index, Integral)
            or int(parent_site_index) < 0
        ):
            _fail(
                "invalid_parent_child_atom_mapping",
                "parent-child mapping parent-site index must be a nonnegative integer",
            )
        translation = _exact_integer_translation(self.parent_cell_translation)
        object.__setattr__(self, "child_atom_id", child_atom_id)
        object.__setattr__(self, "parent_site_index", int(parent_site_index))
        object.__setattr__(self, "parent_orbit_id", parent_orbit_id)
        object.__setattr__(self, "parent_cell_translation", translation)


@dataclass(frozen=True)
class ExactParentChildEmbedding:
    """Exact lattice, origin and atom-image embedding shared by every writer."""

    parent_structure: Structure | IStructure
    parent_space_group_number: int
    basis: RationalMatrix
    origin: tuple[Fraction, Fraction, Fraction]
    atom_mappings: tuple[ParentChildSiteMapping, ...]
    embedding_id: str = ""

    def __post_init__(self) -> None:
        parent = _snapshot_structure(self.parent_structure)
        _require_full_single_species(parent, "parent")
        number = self.parent_space_group_number
        if (
            isinstance(number, bool)
            or not isinstance(number, Integral)
            or not 1 <= int(number) <= 230
        ):
            _fail(
                "invalid_parent_space_group",
                "parent space-group number must be an integer from 1 through 230",
            )
        try:
            basis = rational_matrix(self.basis)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise DisplaciveModelError(
                "invalid_parent_child_basis",
                "parent-to-child basis must be an exact nonsingular 3x3 matrix",
            ) from exc
        origin = _exact_origin(self.origin)
        try:
            mappings = tuple(self.atom_mappings)
        except TypeError as exc:
            raise DisplaciveModelError(
                "invalid_parent_child_atom_mapping",
                "parent-child atom mappings must be a finite immutable sequence",
            ) from exc
        if not mappings or any(
            not isinstance(mapping, ParentChildSiteMapping) for mapping in mappings
        ):
            _fail(
                "invalid_parent_child_atom_mapping",
                "parent-child atom mappings must contain explicit mapping records",
            )
        child_ids = [mapping.child_atom_id for mapping in mappings]
        if len(set(child_ids)) != len(child_ids):
            _fail(
                "invalid_parent_child_atom_mapping",
                "every emitted child atom must occur in the parent-child mapping exactly once",
            )
        if any(mapping.parent_site_index >= len(parent) for mapping in mappings):
            _fail(
                "invalid_parent_child_atom_mapping",
                "parent-child mapping references a parent-site index outside the parent structure",
            )
        parent_primitive = centering_primitive_matrix(
            _centering_letter(int(number))
        )
        relative = multiply(basis, inverse(parent_primitive))
        if not is_integral(relative):
            _fail(
                "invalid_parent_child_translation_lattice",
                "child conventional lattice is not a sublattice of the parent primitive translation lattice",
            )
        translation_index = abs(determinant(relative))
        if translation_index.denominator != 1 or translation_index < 1:
            _fail(
                "invalid_parent_child_translation_lattice",
                "parent/child conventional translation quotient has no positive finite integer index",
            )
        stated_embedding_id = str(self.embedding_id).strip()
        embedding_id = stated_embedding_id or _affine_embedding_id(
            int(number), basis, origin
        )
        object.__setattr__(self, "parent_structure", parent)
        object.__setattr__(self, "parent_space_group_number", int(number))
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "origin", origin)
        object.__setattr__(self, "atom_mappings", mappings)
        object.__setattr__(self, "embedding_id", embedding_id)
        object.__setattr__(self, "_conventional_translation_index", int(translation_index))
        object.__setattr__(self, "_parent_primitive_lattice", parent_primitive)

    @property
    def conventional_translation_index(self) -> int:
        """Order of ``T_parent,primitive / T_child,conventional``."""

        return int(self._conventional_translation_index)

    @property
    def parent_primitive_lattice(self) -> RationalMatrix:
        return self._parent_primitive_lattice

    def expected_reference_lattice(self) -> np.ndarray:
        return np.asarray(self.basis, dtype=float) @ np.asarray(
            self.parent_structure.lattice.matrix, dtype=float
        )

    def validate_reference_lattice(
        self,
        reference: Structure | IStructure,
        *,
        tolerance_angstrom: float,
    ) -> None:
        actual = np.asarray(reference.lattice.matrix, dtype=float)
        expected = self.expected_reference_lattice()
        if not np.allclose(
            actual,
            expected,
            rtol=2.0e-12,
            atol=float(tolerance_angstrom),
        ):
            _fail(
                "parent_child_lattice_mismatch",
                "reference child lattice is not the exact parent-to-child basis times the parent lattice",
            )

    def validate_subgroup(self, subgroup) -> None:
        try:
            subgroup_basis = rational_matrix(subgroup.basis_vectors)
            subgroup_origin = _exact_origin(subgroup.origin)
        except DisplaciveModelError:
            raise
        except (AttributeError, TypeError, ValueError, ZeroDivisionError) as exc:
            raise DisplaciveModelError(
                "subgroup_embedding_mismatch",
                "subgroup does not provide the validated exact basis and origin",
            ) from exc
        if subgroup_basis != self.basis or subgroup_origin != self.origin:
            _fail(
                "subgroup_embedding_mismatch",
                "subgroup basis/origin differs from the validated parent-to-child embedding",
            )
        try:
            subgroup_parent = int(getattr(subgroup, "parent_sg", 0) or 0)
        except (TypeError, ValueError) as exc:
            raise DisplaciveModelError(
                "subgroup_embedding_mismatch",
                "subgroup parent space group is not a valid integer",
            ) from exc
        if subgroup_parent != self.parent_space_group_number:
            _fail(
                "subgroup_embedding_mismatch",
                "subgroup parent space group differs from the validated embedding",
            )
        supplied_embedding_id = str(
            getattr(subgroup, "_displacive_embedding_id", "")
            or getattr(subgroup, "_method3_embedding_id", "")
            or getattr(subgroup, "embedding_id", "")
            or ""
        ).strip()
        if not supplied_embedding_id:
            supplied_embedding_id = _affine_embedding_id(
                subgroup_parent,
                subgroup_basis,
                subgroup_origin,
            )
        if supplied_embedding_id != self.embedding_id:
            _fail(
                "subgroup_embedding_mismatch",
                "subgroup embedding identity differs from the validated parent-child embedding",
            )


@dataclass(frozen=True, init=False)
class DisplaciveSubgroupIdentity:
    """Complete immutable scientific identity bound to one export contract.

    The primitive-cell size ratio is derived from exact translation lattices
    and the target centering.  Callers cannot inject the mutable ``sg.size``
    value that later controls the official ``As`` to ``Ap`` conversion.
    """

    irrep_label: str
    opd_symbol: str
    primary_direction_selector: str
    size: int
    parent_space_group_number: int
    target_space_group_number: int
    target_hall_number: int | None
    basis: RationalMatrix
    origin: tuple[Fraction, Fraction, Fraction]
    embedding_id: str

    def __init__(
        self,
        embedding: ExactParentChildEmbedding,
        *,
        irrep_label: str,
        opd_symbol: str,
        primary_direction_selector: str,
        target_space_group_number: int,
        target_hall_number: int | None = None,
    ) -> None:
        if not isinstance(embedding, ExactParentChildEmbedding):
            _fail(
                "invalid_subgroup_identity",
                "subgroup identity must be bound to an exact parent-child embedding",
            )
        irrep = str(irrep_label).strip()
        opd = str(opd_symbol).strip()
        direction_selector = str(primary_direction_selector).strip()
        target = target_space_group_number
        if (
            not irrep
            or not opd
            or not direction_selector.startswith("VECTOR,")
        ):
            _fail(
                "invalid_subgroup_identity",
                "subgroup identity requires explicit irrep, OPD and exact VECTOR direction",
            )
        if (
            isinstance(target, bool)
            or not isinstance(target, Integral)
            or not 1 <= int(target) <= 230
        ):
            _fail(
                "invalid_subgroup_identity",
                "target space-group number must be an integer from 1 through 230",
            )
        hall: int | None = None
        if target_hall_number is not None:
            if (
                isinstance(target_hall_number, bool)
                or not isinstance(target_hall_number, Integral)
                or not 1 <= int(target_hall_number) <= 530
            ):
                _fail(
                    "invalid_target_hall_setting",
                    "target Hall number must be an integer from 1 through 530",
                )
            hall = int(target_hall_number)
            hall_type = spglib.get_spacegroup_type(hall)
            if hall_type is None or int(hall_type.number) != int(target):
                _fail(
                    "invalid_target_hall_setting",
                    f"Hall number {hall} does not belong to target space group {target}",
                )
        child_centering = (
            _centering_multiplicity(int(target))
            if hall is None
            else _hall_centering_multiplicity(hall)
        )
        translation_index = embedding.conventional_translation_index
        if translation_index % child_centering:
            _fail(
                "invalid_subgroup_size",
                "parent/child translation index is incompatible with target-cell centering",
            )
        size = translation_index // child_centering
        if size < 1:
            _fail("invalid_subgroup_size", "derived subgroup primitive-cell size must be positive")
        object.__setattr__(self, "irrep_label", irrep)
        object.__setattr__(self, "opd_symbol", opd)
        object.__setattr__(
            self, "primary_direction_selector", direction_selector,
        )
        object.__setattr__(self, "size", size)
        object.__setattr__(
            self,
            "parent_space_group_number",
            embedding.parent_space_group_number,
        )
        object.__setattr__(self, "target_space_group_number", int(target))
        object.__setattr__(self, "target_hall_number", hall)
        object.__setattr__(self, "basis", embedding.basis)
        object.__setattr__(self, "origin", embedding.origin)
        object.__setattr__(self, "embedding_id", embedding.embedding_id)

    def validate_subgroup(self, subgroup, embedding: ExactParentChildEmbedding) -> None:
        embedding.validate_subgroup(subgroup)
        try:
            actual = (
                str(getattr(subgroup, "irrep_label", "") or "").strip(),
                str(getattr(subgroup, "opd_symbol", "") or "").strip(),
                int(getattr(subgroup, "size", 0) or 0),
                int(getattr(subgroup, "parent_sg", 0) or 0),
                int(getattr(subgroup, "space_group_number", 0) or 0),
            )
        except (TypeError, ValueError) as exc:
            raise DisplaciveModelError(
                "subgroup_identity_mismatch",
                "subgroup scientific identity contains a non-integer SG or size field",
            ) from exc
        expected = (
            self.irrep_label,
            self.opd_symbol,
            self.size,
            self.parent_space_group_number,
            self.target_space_group_number,
        )
        if actual != expected:
            _fail(
                "subgroup_identity_mismatch",
                "subgroup irrep/OPD/size/parent SG/target SG differs from the validated export identity",
            )
        try:
            from ..backend.iso_wrapper import IsoWrapper  # noqa: PLC0415

            actual_direction = IsoWrapper._direction_vector_selector(
                str(getattr(subgroup, "opd_dir_raw", "") or "")
            )
        except ValueError as exc:
            raise DisplaciveModelError(
                "subgroup_identity_mismatch",
                "subgroup exact OPD direction is missing or invalid",
            ) from exc
        if actual_direction != self.primary_direction_selector:
            _fail(
                "subgroup_identity_mismatch",
                "subgroup exact OPD direction differs from the validated primary VECTOR selector",
            )


def _translation_lattice_residual(
    vector_parent_fractional: np.ndarray,
    translation_lattice: RationalMatrix,
    parent_cartesian_lattice: np.ndarray,
) -> float:
    lattice = np.asarray(translation_lattice, dtype=float)
    coefficients = np.asarray(vector_parent_fractional, dtype=float) @ np.linalg.inv(
        lattice
    )
    nearest = np.rint(coefficients)
    remainder = np.asarray(vector_parent_fractional, dtype=float) - nearest @ lattice
    return float(np.linalg.norm(remainder @ parent_cartesian_lattice))


def _parent_translation_classes(
    parent: IStructure,
    parent_primitive_lattice: RationalMatrix,
    *,
    tolerance_angstrom: float,
) -> tuple[tuple[int, ...], ...]:
    """Partition conventional parent sites into primitive-translation classes."""

    coordinates = np.asarray(parent.frac_coords, dtype=float)
    cartesian_lattice = np.asarray(parent.lattice.matrix, dtype=float)
    remaining = set(range(len(parent)))
    classes: list[tuple[int, ...]] = []
    while remaining:
        representative = min(remaining)
        signature = _chemical_signature(parent[representative])
        members = tuple(
            index
            for index in sorted(remaining)
            if _chemical_signature(parent[index]) == signature
            and _translation_lattice_residual(
                coordinates[index] - coordinates[representative],
                parent_primitive_lattice,
                cartesian_lattice,
            )
            <= tolerance_angstrom
        )
        for index in members:
            remaining.remove(index)
        classes.append(members)
    return tuple(classes)


def _validate_parent_child_site_mapping(
    *,
    embedding: ExactParentChildEmbedding,
    atoms: tuple,
    atom_owner: Mapping[str, str],
    tolerance_angstrom: float,
) -> None:
    """Prove atom origin, chemistry, orbit ownership and all translation cosets."""

    atom_by_id = {atom.atom_id: atom for atom in atoms}
    mappings = embedding.atom_mappings
    mapping_ids = {mapping.child_atom_id for mapping in mappings}
    if mapping_ids != set(atom_by_id):
        missing = sorted(set(atom_by_id).difference(mapping_ids))
        extra = sorted(mapping_ids.difference(atom_by_id))
        _fail(
            "parent_child_atom_mapping_mismatch",
            "parent-child mapping must cover the validated child atom order exactly; "
            f"missing={missing}, extra={extra}",
        )

    parent = embedding.parent_structure
    parent_coordinates = np.asarray(parent.frac_coords, dtype=float)
    parent_lattice = np.asarray(parent.lattice.matrix, dtype=float)
    basis = np.asarray(embedding.basis, dtype=float)
    origin = np.asarray(embedding.origin, dtype=float)
    mappings_by_parent: dict[int, list[ParentChildSiteMapping]] = {}
    for mapping in mappings:
        atom = atom_by_id[mapping.child_atom_id]
        expected_orbit = atom_owner[mapping.child_atom_id]
        if mapping.parent_orbit_id != expected_orbit:
            _fail(
                "parent_child_orbit_mapping_mismatch",
                f"child atom {mapping.child_atom_id!r} is assigned to parent orbit "
                f"{mapping.parent_orbit_id!r}, not {expected_orbit!r}",
            )
        parent_site = parent[mapping.parent_site_index]
        if atom.chemical_signature != _chemical_signature(parent_site):
            _fail(
                "parent_child_chemistry_mismatch",
                f"child atom {mapping.child_atom_id!r} chemistry differs from mapped parent site",
            )
        mappings_by_parent.setdefault(mapping.parent_site_index, []).append(mapping)

    parent_classes = _parent_translation_classes(
        parent,
        embedding.parent_primitive_lattice,
        tolerance_angstrom=tolerance_angstrom,
    )
    expected_cosets = embedding.conventional_translation_index
    for parent_class in parent_classes:
        class_mappings = [
            mapping
            for parent_index in parent_class
            for mapping in mappings_by_parent.get(parent_index, ())
        ]
        if len(class_mappings) != expected_cosets:
            _fail(
                "incomplete_parent_child_translation_cosets",
                "each primitive parent-site class must map to every parent/child "
                f"translation coset exactly once; expected={expected_cosets}, "
                f"found={len(class_mappings)}",
            )
        orbit_ids = {mapping.parent_orbit_id for mapping in class_mappings}
        if len(orbit_ids) != 1:
            _fail(
                "parent_child_orbit_mapping_mismatch",
                "one primitive parent-site class cannot be split across physical parent orbits",
            )
        representative = parent_class[0]
        offsets: list[np.ndarray] = []
        for mapping in class_mappings:
            offset = (
                parent_coordinates[mapping.parent_site_index]
                + np.asarray(mapping.parent_cell_translation, dtype=float)
                - parent_coordinates[representative]
            )
            if (
                _translation_lattice_residual(
                    offset,
                    embedding.parent_primitive_lattice,
                    parent_lattice,
                )
                > tolerance_angstrom
            ):
                _fail(
                    "invalid_parent_child_translation_cosets",
                    "mapped parent images are not translations of one primitive parent site",
                )
            if any(
                _translation_lattice_residual(
                    offset - existing,
                    embedding.basis,
                    parent_lattice,
                )
                <= tolerance_angstrom
                for existing in offsets
            ):
                _fail(
                    "duplicate_parent_child_translation_coset",
                    "parent-child mapping repeats a translation coset in the emitted child cell",
                )
            offsets.append(offset)

    for mapping in mappings:
        atom = atom_by_id[mapping.child_atom_id]
        child_fractional = np.asarray(
            [float(value) for value in atom.frac], dtype=float
        )
        mapped_parent = (
            parent_coordinates[mapping.parent_site_index]
            + np.asarray(mapping.parent_cell_translation, dtype=float)
        )
        transformed_child = child_fractional @ basis + origin
        residual = float(
            np.linalg.norm((transformed_child - mapped_parent) @ parent_lattice)
        )
        if residual > tolerance_angstrom:
            _fail(
                "parent_child_coordinate_mismatch",
                f"child atom {mapping.child_atom_id!r} does not satisfy "
                "x_parent = x_child @ B + q for its explicit parent image",
            )


@dataclass(frozen=True)
class ParentOrbitType:
    """Exact relation between one parent orbit and one IsoVIZ atom type.

    ``child_atom_ids`` are expressed in the validated emitted child frame.
    Keeping this partition beside the parent ``orbit_id`` prevents a writer
    from reconstructing physical orbit identity from a display label, element,
    or displacement magnitude.
    """

    type_index: int
    type_label: str
    child_atom_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        index = self.type_index
        if isinstance(index, bool) or not isinstance(index, Integral) or int(index) <= 0:
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit type index must be a positive integer",
            )
        label = str(self.type_label).strip()
        atom_ids = tuple(str(value).strip() for value in self.child_atom_ids)
        if not label:
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit type label cannot be empty",
            )
        if not atom_ids or any(not value for value in atom_ids):
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit child atom IDs must be nonempty",
            )
        if len(set(atom_ids)) != len(atom_ids):
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit child atom IDs must be unique",
            )
        object.__setattr__(self, "type_index", int(index))
        object.__setattr__(self, "type_label", label)
        object.__setattr__(self, "child_atom_ids", atom_ids)


@dataclass(frozen=True)
class DisplaciveExportMode:
    """One writer-ready mode in the shared atom order and child frame."""

    key: str
    label: str
    global_irrep_label: str
    mode_id: str
    column_index: int
    parent_orbit_id: str
    parent_type_index: int
    parent_type_label: str
    raw_fractional: np.ndarray
    normalized_fractional_per_angstrom: np.ndarray
    normfactor_per_angstrom: float
    amplitude_as_angstrom: float
    max_amplitude_angstrom: float

    def __post_init__(self) -> None:
        global_irrep = str(self.global_irrep_label).strip()
        if not global_irrep:
            _fail(
                "mode_irrep_mismatch",
                "every displacive export mode requires an explicit global irrep label",
            )
        parent_orbit_id = str(self.parent_orbit_id).strip()
        parent_type_label = str(self.parent_type_label).strip()
        parent_type_index = self.parent_type_index
        if not parent_orbit_id or not parent_type_label:
            _fail(
                "mode_parent_orbit_type_mismatch",
                "every displacive export mode requires an explicit parent orbit and type label",
            )
        if (
            isinstance(parent_type_index, bool)
            or not isinstance(parent_type_index, Integral)
            or int(parent_type_index) <= 0
        ):
            _fail(
                "mode_parent_orbit_type_mismatch",
                "every displacive export mode requires a positive parent type index",
            )
        raw = np.asarray(self.raw_fractional, dtype=float)
        normalized = np.asarray(self.normalized_fractional_per_angstrom, dtype=float)
        if raw.ndim != 2 or raw.shape[1] != 3 or normalized.shape != raw.shape:
            _fail("invalid_export_mode_shape", "displacive export modes must be matching finite Nx3 arrays")
        if not np.all(np.isfinite(raw)) or not np.all(np.isfinite(normalized)):
            _fail("invalid_export_mode_shape", "displacive export modes must be finite")
        for name, value in (
            ("normfactor", self.normfactor_per_angstrom),
            ("amplitude", self.amplitude_as_angstrom),
            ("maximum amplitude", self.max_amplitude_angstrom),
        ):
            if not math.isfinite(float(value)):
                _fail("invalid_export_mode_value", f"{name} must be finite")
        if self.normfactor_per_angstrom <= 0 or self.max_amplitude_angstrom <= 0:
            _fail("invalid_export_mode_value", "mode normalization and maximum amplitude must be positive")
        if not np.allclose(
            normalized,
            raw * float(self.normfactor_per_angstrom),
            rtol=2.0e-12,
            atol=2.0e-12,
        ):
            _fail("mode_normalization_mismatch", "normalized mode is not normfactor times the raw mode")
        raw_copy = np.array(raw, dtype=float, copy=True)
        normalized_copy = np.array(normalized, dtype=float, copy=True)
        raw_copy.setflags(write=False)
        normalized_copy.setflags(write=False)
        object.__setattr__(self, "global_irrep_label", global_irrep)
        object.__setattr__(self, "parent_orbit_id", parent_orbit_id)
        object.__setattr__(self, "parent_type_index", int(parent_type_index))
        object.__setattr__(self, "parent_type_label", parent_type_label)
        object.__setattr__(self, "raw_fractional", raw_copy)
        object.__setattr__(self, "normalized_fractional_per_angstrom", normalized_copy)


@dataclass(frozen=True)
class DisplaciveExportData:
    """Single validated source for CIF, IsoVIZ, Complete and TOPAS writers."""

    model: DisplaciveCifModel
    embedding: ExactParentChildEmbedding
    subgroup_identity: DisplaciveSubgroupIdentity
    reference_structure: Structure | IStructure
    representative_labels: Mapping[str, str]
    parent_orbit_types: Mapping[str, ParentOrbitType]
    final_structure: Structure | IStructure | None = None
    coordinate_tolerance_angstrom: float = 1.0e-8

    def __post_init__(self) -> None:
        if not isinstance(self.model, DisplaciveCifModel) or self.model.diagnostics.status != "validated":
            _fail("unvalidated_displacive_model", "export data requires a validated displacive CIF model")
        tolerance = float(self.coordinate_tolerance_angstrom)
        if not math.isfinite(tolerance) or tolerance <= 0:
            _fail("invalid_export_tolerance", "coordinate tolerance must be finite and positive")

        if not isinstance(self.embedding, ExactParentChildEmbedding):
            _fail(
                "missing_parent_child_embedding",
                "displacive export requires an exact parent-to-child embedding",
            )
        if not isinstance(self.subgroup_identity, DisplaciveSubgroupIdentity):
            _fail(
                "missing_subgroup_identity",
                "displacive export requires a complete immutable subgroup identity",
            )
        identity = self.subgroup_identity
        if (
            identity.parent_space_group_number
            != self.embedding.parent_space_group_number
            or identity.basis != self.embedding.basis
            or identity.origin != self.embedding.origin
            or identity.embedding_id != self.embedding.embedding_id
        ):
            _fail(
                "subgroup_identity_mismatch",
                "subgroup identity was not derived from the export parent-child embedding",
            )
        model = self.model
        _validate_emitted_space_group(
            model.frame.operations,
            model_centering_index=model.centering.index,
            target_space_group_number=identity.target_space_group_number,
            target_hall_number=identity.target_hall_number,
        )
        atoms = tuple(model.atoms)
        atom_ids = tuple(atom.atom_id for atom in atoms)
        if any(tuple(mode.atom_ids) != atom_ids for mode in model.canonical_modes):
            _fail(
                "atom_order_mismatch",
                "canonical mode columns do not use the validated model atom order",
            )
        try:
            raw_parent_orbit_types = dict(self.parent_orbit_types)
            parent_orbit_types = {
                str(orbit_id).strip(): parent_type
                for orbit_id, parent_type in raw_parent_orbit_types.items()
            }
        except (TypeError, ValueError) as exc:
            raise DisplaciveModelError(
                "invalid_parent_orbit_type",
                "parent-orbit type mapping must be a finite mapping",
            ) from exc
        if (
            not parent_orbit_types
            or "" in parent_orbit_types
            or len(parent_orbit_types) != len(raw_parent_orbit_types)
        ):
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit type mapping requires unique nonempty physical orbit IDs",
            )
        if any(
            not isinstance(parent_type, ParentOrbitType)
            for parent_type in parent_orbit_types.values()
        ):
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit type mapping values must be ParentOrbitType records",
            )
        type_indices = [parent_type.type_index for parent_type in parent_orbit_types.values()]
        type_labels = [parent_type.type_label for parent_type in parent_orbit_types.values()]
        if len(set(type_indices)) != len(type_indices) or set(type_indices) != set(
            range(1, len(type_indices) + 1)
        ):
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit type indices must be unique, contiguous and one-based",
            )
        if len(set(type_labels)) != len(type_labels):
            _fail(
                "invalid_parent_orbit_type",
                "parent-orbit type labels must be unique",
            )
        atom_owner: dict[str, str] = {}
        atom_by_id = {atom.atom_id: atom for atom in atoms}
        for orbit_id, parent_type in parent_orbit_types.items():
            signatures = set()
            for atom_id in parent_type.child_atom_ids:
                if atom_id not in atom_by_id:
                    _fail(
                        "parent_orbit_atom_mismatch",
                        f"parent orbit {orbit_id!r} references unknown child atom {atom_id!r}",
                    )
                if atom_id in atom_owner:
                    _fail(
                        "parent_orbit_atom_mismatch",
                        f"child atom {atom_id!r} belongs to more than one parent orbit type",
                    )
                atom_owner[atom_id] = orbit_id
                signatures.add(atom_by_id[atom_id].chemical_signature)
            if len(signatures) != 1:
                _fail(
                    "parent_orbit_atom_mismatch",
                    f"parent orbit {orbit_id!r} mixes distinct child-site chemistries",
                )
        if set(atom_owner) != set(atom_ids):
            missing = sorted(set(atom_ids).difference(atom_owner))
            _fail(
                "parent_orbit_atom_mismatch",
                "parent-orbit type mapping does not cover every child atom: "
                + ", ".join(missing),
            )
        for mode in model.canonical_modes:
            mode_identity = getattr(mode, "mode_identity", None)
            orbit_id = str(
                getattr(mode_identity, "orbit_id", "") or ""
            ).strip()
            if orbit_id and orbit_id not in parent_orbit_types:
                _fail(
                    "missing_parent_orbit_type",
                    f"canonical mode {mode.mode_id!r} refers to unmapped parent orbit {orbit_id!r}",
                )
        reference = _snapshot_structure(self.reference_structure)
        _require_full_single_species(reference, "reference child")
        if len(reference) != len(atoms):
            _fail("atom_order_mismatch", "reference child atom count does not match the validated model")
        self.embedding.validate_reference_lattice(
            reference,
            tolerance_angstrom=tolerance,
        )
        frame_lattice = np.asarray(model.frame.reference_lattice, dtype=float)
        reference_lattice = np.asarray(reference.lattice.matrix, dtype=float)
        if not np.allclose(frame_lattice, reference_lattice, rtol=2.0e-12, atol=tolerance):
            _fail("reference_frame_mismatch", "reference child lattice differs from the validated mode frame")
        exact_coords = np.asarray([[float(value) for value in atom.frac] for atom in atoms], dtype=float)
        residual = _periodic_cartesian_residual(
            np.asarray(reference.frac_coords, dtype=float), exact_coords, reference_lattice
        )
        if residual > tolerance:
            _fail("reference_atom_mapping_mismatch", "reference child coordinates do not match the exact atom order")
        for atom, site in zip(atoms, reference, strict=True):
            if atom.chemical_signature != _chemical_signature(site):
                _fail("reference_atom_mapping_mismatch", "reference child chemistry differs from the exact atom order")
        _validate_parent_child_site_mapping(
            embedding=self.embedding,
            atoms=atoms,
            atom_owner=atom_owner,
            tolerance_angstrom=tolerance,
        )

        amplitudes = np.asarray(model.amplitudes, dtype=float)
        norms = np.asarray(model.norms.normfactors, dtype=float)
        displacement = np.zeros((len(reference), 3), dtype=float)
        for amplitude, normfactor, mode in zip(
            amplitudes, norms, model.canonical_modes, strict=True
        ):
            displacement += amplitude * normfactor * np.asarray(mode.displacements, dtype=float)
        expected_coords = np.mod(np.asarray(reference.frac_coords, dtype=float) + displacement, 1.0)
        final = self.final_structure
        if final is None:
            final = IStructure(
                lattice=reference.lattice,
                species=[site.species for site in reference],
                coords=expected_coords,
                coords_are_cartesian=False,
                site_properties=deepcopy(reference.site_properties),
                labels=[site.label for site in reference],
                properties=deepcopy(reference.properties),
            )
        else:
            final = _snapshot_structure(final)
        _require_full_single_species(final, "final child")
        if len(final) != len(reference):
            _fail("final_atom_mapping_mismatch", "final child atom count differs from the reference child")
        final_lattice = np.asarray(final.lattice.matrix, dtype=float)
        if not np.allclose(
            final_lattice,
            reference_lattice,
            rtol=2.0e-12,
            atol=tolerance,
        ):
            _fail(
                "final_frame_mismatch",
                "pure displacive export requires the final child lattice to match the validated reference frame",
            )
        for before, after in zip(reference, final, strict=True):
            if before.species != after.species:
                _fail("final_atom_mapping_mismatch", "final child atom order or chemistry differs from the reference")
        final_residual = _periodic_cartesian_residual(
            np.asarray(final.frac_coords, dtype=float), expected_coords, reference_lattice
        )
        if final_residual > tolerance:
            _fail(
                "final_displacement_mismatch",
                "final child coordinates are not reference plus the declared normalized mode amplitudes",
            )

        representative_ids = {orbit.representative_atom_id for orbit in model.orbits}
        representative_labels = {
            str(key): str(value).strip()
            for key, value in dict(self.representative_labels).items()
        }
        if set(representative_labels) != representative_ids or any(
            not value for value in representative_labels.values()
        ):
            _fail(
                "representative_label_mismatch",
                "representative labels must cover every validated child orbit exactly once",
            )
        if len(set(representative_labels.values())) != len(representative_labels):
            _fail("representative_label_mismatch", "representative labels must be unique")

        object.__setattr__(self, "reference_structure", reference)
        object.__setattr__(
            self,
            "parent_orbit_types",
            MappingProxyType(parent_orbit_types),
        )
        object.__setattr__(
            self,
            "representative_labels",
            MappingProxyType(representative_labels),
        )
        object.__setattr__(self, "final_structure", final)
        # Preserve the more specific canonical source-identity diagnostics
        # before validating the encompassing ISO query context.
        self.modes()
        primary_direction_selectors: set[str] = set()
        for mode in model.canonical_modes:
            provenance = mode.microscopic_provenance
            if provenance is None:
                _fail(
                    "microscopic_query_context_mismatch",
                    "canonical mode lacks microscopic query provenance",
                )
            source_basis = tuple(
                tuple(Fraction(value) for value in row)
                for row in provenance.source_subgroup_basis
            )
            source_origin = tuple(
                Fraction(value) for value in provenance.source_subgroup_origin
            )
            if (
                provenance.source_parent_sg
                != self.embedding.parent_space_group_number
                or provenance.source_subgroup_space_group_number
                != identity.target_space_group_number
                or source_basis != self.embedding.basis
                or source_origin != self.embedding.origin
                or provenance.source_subgroup_irrep_label != identity.irrep_label
                or provenance.source_subgroup_opd_symbol != identity.opd_symbol
            ):
                _fail(
                    "microscopic_query_context_mismatch",
                    "canonical mode query context differs from the exact export embedding "
                    "or subgroup identity",
                )
            primary_direction_selectors.add(
                provenance.source_subgroup_primary_direction_selector
            )
        if primary_direction_selectors != {identity.primary_direction_selector}:
            _fail(
                "microscopic_query_context_mismatch",
                "canonical primary modes do not identify one exact subgroup VECTOR direction",
            )

    @property
    def parent_structure(self) -> IStructure:
        return self.embedding.parent_structure

    @property
    def parent_space_group_number(self) -> int:
        return self.embedding.parent_space_group_number

    @property
    def atom_ids(self) -> tuple[str, ...]:
        return tuple(atom.atom_id for atom in self.model.atoms)

    @property
    def mode_keys(self) -> tuple[str, ...]:
        return tuple(row.key for row in self.modes())

    @property
    def mode_labels(self) -> tuple[str, ...]:
        return tuple(row.label for row in self.modes())

    @property
    def mode_global_irreps(self) -> tuple[str, ...]:
        return tuple(row.global_irrep_label for row in self.modes())

    def ordered_parent_orbit_types(self) -> tuple[tuple[str, ParentOrbitType], ...]:
        """Return parent orbit types in their explicit IsoVIZ order."""

        return tuple(
            sorted(
                self.parent_orbit_types.items(),
                key=lambda item: item[1].type_index,
            )
        )

    def parent_type_species(self, orbit_id: str) -> str:
        """Return the single validated species token for one parent type."""

        try:
            parent_type = self.parent_orbit_types[str(orbit_id)]
        except KeyError as exc:
            raise DisplaciveModelError(
                "missing_parent_orbit_type",
                f"parent orbit {orbit_id!r} has no export type",
            ) from exc
        atom_by_id = {atom.atom_id: atom for atom in self.model.atoms}
        signature = atom_by_id[parent_type.child_atom_ids[0]].chemical_signature
        if len(signature) != 1 or signature[0][1] != 1:
            _fail(
                "unsupported_site_occupancy",
                f"parent orbit {orbit_id!r} is not a fully occupied single-species type",
            )
        return signature[0][0]

    def validate_subgroup(self, subgroup) -> None:
        self.subgroup_identity.validate_subgroup(subgroup, self.embedding)

    def snapshot_parent(self) -> IStructure:
        return _snapshot_structure(self.parent_structure)

    def snapshot_reference(self) -> IStructure:
        return _snapshot_structure(self.reference_structure)

    def snapshot_final(self) -> IStructure:
        if self.final_structure is None:  # pragma: no cover - established in __post_init__
            _fail("missing_final_structure", "validated export data has no final child")
        return _snapshot_structure(self.final_structure)

    def modes(self) -> tuple[DisplaciveExportMode, ...]:
        """Return writer rows in the validated source-column order."""

        lattice = np.asarray(self.reference_structure.lattice.matrix, dtype=float)
        result: list[DisplaciveExportMode] = []
        for mode, normfactor, amplitude in zip(
            self.model.canonical_modes,
            self.model.norms.normfactors,
            self.model.amplitudes,
            strict=True,
        ):
            if not bool(getattr(mode, "has_resolved_source", False)):
                _fail(
                    "unresolved_mode_identity",
                    f"mode {mode.mode_id!r} has no verified canonical microscopic identity",
                )
            key = str(mode.mode_id).strip()
            label = str(getattr(mode, "display_label", "")).strip()
            global_irrep = str(getattr(mode, "global_irrep", "")).strip()
            parent_orbit_id = str(mode.mode_identity.orbit_id).strip()
            try:
                parent_type = self.parent_orbit_types[parent_orbit_id]
            except KeyError as exc:
                raise DisplaciveModelError(
                    "missing_parent_orbit_type",
                    f"mode {mode.mode_id!r} has no parent-orbit type",
                ) from exc
            identity_parent_sg = int(mode.mode_identity.parent_sg)
            if identity_parent_sg != self.parent_space_group_number:
                _fail(
                    "mode_parent_space_group_mismatch",
                    f"mode {mode.mode_id!r} belongs to parent space group "
                    f"{identity_parent_sg}, not {self.parent_space_group_number}",
                )
            if not key or not label or not global_irrep:
                _fail(
                    "unresolved_mode_identity",
                    "canonical mode key, display label and global irrep must be explicit",
                )
            raw = np.asarray(mode.displacements, dtype=float)
            normalized = raw * float(normfactor)
            cartesian_norms = np.linalg.norm(normalized @ lattice, axis=1)
            dmax = float(np.max(cartesian_norms))
            if not math.isfinite(dmax) or dmax <= 0:
                _fail("invalid_mode_maximum", f"mode {mode.mode_id!r} has no finite Cartesian displacement")
            active_tolerance = max(
                1.0e-12,
                256.0 * np.finfo(float).eps * max(1.0, dmax),
            )
            active_atom_ids = {
                atom_id
                for atom_id, norm in zip(self.atom_ids, cartesian_norms, strict=True)
                if float(norm) > active_tolerance
            }
            unexpected = sorted(
                active_atom_ids.difference(parent_type.child_atom_ids)
            )
            if unexpected:
                _fail(
                    "mode_parent_orbit_type_mismatch",
                    f"mode {mode.mode_id!r} for parent orbit {parent_orbit_id!r} "
                    "has nonzero rows on child atoms outside that orbit: "
                    + ", ".join(unexpected),
                )
            result.append(
                DisplaciveExportMode(
                    key=key,
                    label=label,
                    global_irrep_label=global_irrep,
                    mode_id=mode.mode_id,
                    column_index=mode.column_index,
                    parent_orbit_id=parent_orbit_id,
                    parent_type_index=parent_type.type_index,
                    parent_type_label=parent_type.type_label,
                    raw_fractional=raw,
                    normalized_fractional_per_angstrom=normalized,
                    normfactor_per_angstrom=float(normfactor),
                    amplitude_as_angstrom=float(amplitude),
                    max_amplitude_angstrom=1.0 / dmax,
                )
            )
        return tuple(result)

    def representative_label(self, atom_id: str) -> str:
        try:
            return self.representative_labels[str(atom_id)]
        except KeyError as exc:
            raise DisplaciveModelError(
                "representative_label_mismatch",
                f"validated representative {atom_id!r} has no export label",
            ) from exc


__all__ = [
    "DisplaciveExportData",
    "DisplaciveExportMode",
    "DisplaciveSubgroupIdentity",
    "ExactParentChildEmbedding",
    "ParentChildSiteMapping",
    "ParentOrbitType",
]
