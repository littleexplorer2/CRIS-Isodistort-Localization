"""CIF structure and metadata comparison, independent of UI and file pairing."""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import spglib
from pymatgen.core import Structure

SPACE_GROUP_NUMBER_TAGS = ("_symmetry_Int_Tables_number", "_space_group_IT_number")
SPACE_GROUP_SYMBOL_TAGS = ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt")
LABEL_TAG = "_atom_site_label"
MAGMOM_TAGS = (
    "_atom_site_moment_crystalaxis_x",
    "_atom_site_moment_crystalaxis_y",
    "_atom_site_moment_crystalaxis_z",
)


def _first_value(block: dict[str, Any], tags: tuple[str, ...]) -> Any | None:
    for tag in tags:
        if tag in block:
            value = block[tag]
            return value[0] if isinstance(value, list) else value
    return None


def _numeric(value: Any) -> float:
    text = re.sub(r"\([^)]*\)$", "", str(value).strip().strip("'\""))
    return float(text)


def _numeric_list(block: dict[str, Any], tag: str) -> list[float] | None:
    if tag not in block:
        return None
    values = block[tag]
    if not isinstance(values, list):
        values = [values]
    return [_numeric(value) for value in values]


def _periodic_delta(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    return (second - first) - np.round(second - first)


def _match_atoms(
    local: Structure,
    reference: Structure,
    coordinate_tolerance: float,
    ignore_atom_order: bool,
) -> tuple[list[int] | None, float | None]:
    if len(local) != len(reference):
        return None, None
    if not ignore_atom_order:
        delta = _periodic_delta(local.frac_coords, reference.frac_coords)
        maximum = float(np.max(np.abs(delta))) if len(local) else 0.0
        if maximum > coordinate_tolerance:
            return None, maximum
        return list(range(len(local))), maximum

    candidates: dict[int, list[tuple[float, int]]] = {index: [] for index in range(len(local))}
    for local_index, local_site in enumerate(local):
        for reference_index, reference_site in enumerate(reference):
            if local_site.species_string != reference_site.species_string:
                continue
            delta = _periodic_delta(
                np.asarray(local_site.frac_coords),
                np.asarray(reference_site.frac_coords),
            )
            distance = float(np.max(np.abs(delta)))
            if distance <= coordinate_tolerance:
                candidates[local_index].append((distance, reference_index))

    # A greedy edge choice can reject a valid mapping when one site has several
    # neighbours but another site has only one. Augmenting paths find a perfect
    # species-constrained bipartite matching whenever one exists.
    for values in candidates.values():
        values.sort()
    reference_to_local: dict[int, int] = {}

    def assign(local_index: int, visited: set[int]) -> bool:
        for _, reference_index in candidates[local_index]:
            if reference_index in visited:
                continue
            visited.add(reference_index)
            previous = reference_to_local.get(reference_index)
            if previous is None or assign(previous, visited):
                reference_to_local[reference_index] = local_index
                return True
        return False

    if any(not assign(local_index, set()) for local_index in range(len(local))):
        return None, None
    assignments = {local_index: reference_index for reference_index, local_index in reference_to_local.items()}
    ordered = [assignments[index] for index in range(len(local))]
    distances = [
        next(distance for distance, candidate in candidates[local_index] if candidate == reference_index)
        for local_index, reference_index in assignments.items()
    ]
    return ordered, max(distances, default=0.0)


def _spglib_space_group(structure: Structure, symprec: float) -> dict[str, Any]:
    numbers = [site.specie.Z for site in structure]
    dataset = spglib.get_symmetry_dataset((structure.lattice.matrix, structure.frac_coords, numbers), symprec=symprec)
    if dataset is None:
        return {"number": None, "symbol": None}
    return {"number": int(dataset.number), "symbol": str(dataset.international)}


def _compare_mapped_scalars(
    local: list[float] | None,
    reference: list[float] | None,
    assignments: list[int] | None,
    tolerance: float,
) -> bool | None:
    if local is None and reference is None:
        return None
    if local is None or reference is None or assignments is None:
        return False
    if len(local) != len(reference) or len(local) != len(assignments):
        return False
    return all(abs(local[i] - reference[assignments[i]]) <= tolerance for i in range(len(local)))


def _site_occupancies(structure: Structure) -> list[float]:
    """Occupancy totals for the fully expanded sites returned by pymatgen."""
    return [float(sum(site.species.values())) for site in structure]


def _compare_mapped_vectors(
    local: list[list[float] | None] | None,
    reference: list[list[float] | None] | None,
    assignments: list[int] | None,
    tolerance: float,
) -> bool | None:
    if local is None and reference is None:
        return None
    if local is None or reference is None or assignments is None:
        return False
    if any(value is None for value in local + reference):
        return False
    local_columns = [value for value in local if value is not None]
    reference_columns = [value for value in reference if value is not None]
    if len(local_columns) != 3 or len(reference_columns) != 3:
        return False
    if any(len(column) != len(assignments) for column in local_columns + reference_columns):
        return False
    return all(
        abs(local_columns[axis][index] - reference_columns[axis][reference_index]) <= tolerance
        for index, reference_index in enumerate(assignments)
        for axis in range(3)
    )


def _vector_row_count(values: list[list[float] | None]) -> int | None:
    if len(values) != 3 or any(column is None for column in values):
        return None
    lengths = {len(column) for column in values if column is not None}
    return lengths.pop() if len(lengths) == 1 else None


def _declared_space_group(block: dict[str, Any]) -> dict[str, Any]:
    number = _first_value(block, SPACE_GROUP_NUMBER_TAGS)
    symbol = _first_value(block, SPACE_GROUP_SYMBOL_TAGS)
    normalized_symbol = None
    if symbol is not None:
        normalized_symbol = " ".join(str(symbol).strip("'\"").split())
    return {
        "number": int(_numeric(number)) if number is not None else None,
        "symbol": normalized_symbol,
    }


def _compare_optional_metadata(first: Any, second: Any) -> bool | None:
    if first is None and second is None:
        return None
    return first == second


def _validate_declared_space_group(declared: dict[str, Any], inferred: dict[str, Any]) -> bool:
    return declared["number"] is None or declared["number"] == inferred["number"]


def compare_structures(
    local: Structure,
    reference: Structure,
    local_block: dict[str, Any],
    reference_block: dict[str, Any],
    lattice_tolerance: float,
    coordinate_tolerance: float,
    scalar_tolerance: float,
    ignore_atom_order: bool,
) -> tuple[bool, dict[str, Any]]:
    details: dict[str, Any] = {
        "lattice_tolerance": lattice_tolerance,
        "coordinate_tolerance": coordinate_tolerance,
        "scalar_tolerance": scalar_tolerance,
        "ignore_atom_order": ignore_atom_order,
        "local_atom_count": len(local),
        "reference_atom_count": len(reference),
        "lattice_max_abs_difference": None,
        "coordinate_max_periodic_difference": None,
        "species_order_equal": None,
        "atom_order_equal": None,
        "lattice_equal": None,
        "coordinates_equal": None,
        "occupancies_equal": None,
        "magnetic_moments_equal": None,
        "labels_equal": None,
        "declared_space_group_equal": None,
        "inferred_space_group_equal": None,
        "space_group": {},
        "issues": [],
    }
    if len(local) != len(reference):
        details["issues"].append("atom count differs")
        return False, details

    lattice_max = float(np.max(np.abs(local.lattice.matrix - reference.lattice.matrix)))
    species_local = [site.species_string for site in local]
    species_reference = [site.species_string for site in reference]
    species_order_equal = species_local == species_reference
    assignments, coordinate_max = _match_atoms(local, reference, coordinate_tolerance, ignore_atom_order)
    coordinates_equal = assignments is not None
    if not coordinates_equal and coordinate_max is None:
        coordinate_max = float("inf")
    lattice_equal = lattice_max <= lattice_tolerance
    atom_order_equal = species_order_equal and not ignore_atom_order

    occupancies_equal = _compare_mapped_scalars(
        _site_occupancies(local),
        _site_occupancies(reference),
        assignments,
        scalar_tolerance,
    )
    local_moments = [_numeric_list(local_block, tag) for tag in MAGMOM_TAGS]
    reference_moments = [_numeric_list(reference_block, tag) for tag in MAGMOM_TAGS]
    magnetic_present = any(value is not None for value in local_moments + reference_moments)
    if not magnetic_present:
        magnetic_equal = None
    elif ignore_atom_order:
        # Raw magnetic loops describe ASU rows.  They can use the expanded-site
        # mapping only when the CIF explicitly lists every expanded site.
        expanded = len(local)
        can_map = _vector_row_count(local_moments) == expanded and _vector_row_count(reference_moments) == len(
            reference
        )
        magnetic_equal = (
            _compare_mapped_vectors(local_moments, reference_moments, assignments, scalar_tolerance)
            if can_map
            else None
        )
    else:
        local_rows = _vector_row_count(local_moments)
        reference_rows = _vector_row_count(reference_moments)
        identity = list(range(local_rows)) if local_rows == reference_rows and local_rows is not None else None
        magnetic_equal = _compare_mapped_vectors(local_moments, reference_moments, identity, scalar_tolerance)
    # CIF labels identify rows in the asymmetric-unit loop; they are not site
    # properties on pymatgen's symmetry-expanded Structure.  Consequently the
    # expanded-site assignment cannot index those raw lists.  When row order is
    # significant compare labels in their declared order.  When row order is
    # intentionally ignored, labels are display metadata and are not a semantic
    # failure condition; species, coordinates and expanded occupancies still are.
    labels_equal = (
        None
        if ignore_atom_order
        else _compare_optional_metadata(local_block.get(LABEL_TAG), reference_block.get(LABEL_TAG))
    )

    declared_local = _declared_space_group(local_block)
    declared_reference = _declared_space_group(reference_block)
    declared_equal = _compare_optional_metadata(declared_local, declared_reference)
    symprec = max(lattice_tolerance, coordinate_tolerance * max(local.lattice.abc))
    inferred_local = _spglib_space_group(local, symprec)
    inferred_reference = _spglib_space_group(reference, symprec)
    inferred_equal = inferred_local == inferred_reference
    declared_valid = _validate_declared_space_group(declared_local, inferred_local)
    declared_valid = declared_valid and _validate_declared_space_group(declared_reference, inferred_reference)

    issues: list[str] = []
    if not lattice_equal:
        issues.append("lattice differs beyond lattice tolerance")
    if not species_order_equal and not ignore_atom_order:
        issues.append("species or atom order differs")
    if not coordinates_equal:
        issues.append("fractional coordinates differ beyond tolerance")
    if occupancies_equal is False:
        issues.append("occupancies differ beyond scalar tolerance")
    if magnetic_equal is False:
        issues.append("magnetic moments differ beyond scalar tolerance")
    if labels_equal is False:
        issues.append("atom labels differ")
    if declared_equal is False:
        issues.append("declared space group differs")
    if not inferred_equal:
        issues.append("spglib-inferred space group differs")
    if not declared_valid:
        issues.append("declared space group does not match spglib inference")

    structure_equal = (
        lattice_equal
        and coordinates_equal
        and occupancies_equal is not False
        and magnetic_equal is not False
        and labels_equal is not False
        and declared_equal is not False
        and inferred_equal
    )
    details.update(
        {
            "lattice_max_abs_difference": lattice_max,
            "coordinate_max_periodic_difference": coordinate_max,
            "species_order_equal": species_order_equal,
            "atom_order_equal": atom_order_equal,
            "lattice_equal": lattice_equal,
            "coordinates_equal": coordinates_equal,
            "occupancies_equal": occupancies_equal,
            "magnetic_moments_equal": magnetic_equal,
            "labels_equal": labels_equal,
            "declared_space_group_equal": declared_equal,
            "inferred_space_group_equal": inferred_equal,
            "space_group": {
                "local_declared": declared_local,
                "reference_declared": declared_reference,
                "local_inferred": inferred_local,
                "reference_inferred": inferred_reference,
            },
            "issues": issues,
        }
    )
    return structure_equal, details
