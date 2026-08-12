"""Search workflows aligned with ISODISTORT Method 1-4."""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Optional, Sequence

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.groups import SpaceGroup

from ..backend import DistortionMode, IsoWrapper, SubgroupInfo
from .phase_path import normalize_distortion_types


CRYSTAL_SYSTEMS = {
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic",
}


def _sg_to_crystal_system(space_group_number: int) -> str:
    if 1 <= space_group_number <= 2:
        return "triclinic"
    if 3 <= space_group_number <= 15:
        return "monoclinic"
    if 16 <= space_group_number <= 74:
        return "orthorhombic"
    if 75 <= space_group_number <= 142:
        return "tetragonal"
    if 143 <= space_group_number <= 167:
        return "trigonal"
    if 168 <= space_group_number <= 194:
        return "hexagonal"
    if 195 <= space_group_number <= 230:
        return "cubic"
    return "unknown"


def _to_float(value: str | int | float) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if "/" in text:
        return float(Fraction(text))
    return float(text)


def _parse_basis_rows(rows: Sequence[Sequence[str | int | float]]) -> List[List[float]]:
    basis = []
    for row in rows:
        if len(row) != 3:
            raise ValueError("Each basis row must contain 3 values")
        basis.append([_to_float(v) for v in row])
    if len(basis) != 3:
        raise ValueError("Basis matrix must be 3x3")
    return basis


def _matrix_matches_diagonal_supercell(matrix: Sequence[Sequence[float]],
                                       supercell: Sequence[int]) -> bool:
    if len(supercell) != 3:
        return False
    target = np.diag(np.array(supercell, dtype=float))
    got = np.array(matrix, dtype=float)
    return np.allclose(got, target, atol=1e-8)


def _space_group_to_point_group(space_group_number: int) -> str:
    try:
        return SpaceGroup.from_int_number(space_group_number).point_group
    except Exception:
        return ""


@dataclass
class Method1Query:
    """Method 1: search over all special k points."""

    distortion_types: str | Sequence[str] | None = None
    crystal_system: Optional[str] = None
    subgroup_space_group: Optional[int] = None
    direct_sublattice: Optional[Sequence[int]] = None
    maximal_subgroup_only: bool = False


@dataclass
class Method1ResultItem:
    subgroup: SubgroupInfo
    crystal_system: str
    is_maximal: bool
    direct_sublattice: List[int] = field(default_factory=lambda: [1, 1, 1])


@dataclass
class Method2Query:
    """Method 2: general method over specified k point(s)."""

    subgroup_idx: int
    distortion_type: str = "displacement"
    k_point_label: Optional[str] = None
    k_point_coordinates: Optional[Sequence[str | int | float]] = None
    k_parameters: Dict[str, str | int | float] = field(default_factory=dict)
    number_of_independent_modulations: int = 0
    number_of_superposed_irs: int = 1
    specified_opd: Optional[str] = None


@dataclass
class Method2Result:
    subgroup: SubgroupInfo
    modes: List[DistortionMode]
    metadata: Dict[str, object]


@dataclass
class Method3Query:
    """Method 3: search over arbitrary k for point/space group + supercell."""

    distortion_types: str | Sequence[str] | None = None
    point_group: Optional[str] = None
    space_group_type: Optional[int] = None
    supercell_basis: Optional[Sequence[Sequence[str | int | float]]] = None
    direct_sublattice_centering: Optional[str] = None


@dataclass
class Method3ResultItem:
    subgroup: SubgroupInfo
    point_group: str
    basis: List[List[float]]


@dataclass
class Method4Query:
    """Method 4: mode decomposition of distorted structure."""

    atom_matching_method: str = "nearest-site"
    robust_distance_threshold: float = 0.25
    provided_origin_shift: Optional[Sequence[float]] = None


@dataclass
class Method4Result:
    amplitudes: Dict[str, float]
    rms_residual: float
    max_abs_residual: float
    assignments: List[int]
    metadata: Dict[str, object]


class IsoSearchEngine:
    """Implements local Method 1-4 search workflows."""

    def __init__(self, iso_wrapper: IsoWrapper):
        self._iso = iso_wrapper

    def method_1_search(self, parent_sg: int, query: Method1Query) -> List[Method1ResultItem]:
        distortion_types = normalize_distortion_types(query.distortion_types)
        subgroups = self._iso.list_subgroups(parent_sg, distortion_types)

        result: List[Method1ResultItem] = []
        for sg in subgroups:
            crystal_system = _sg_to_crystal_system(sg.space_group_number)
            item = Method1ResultItem(
                subgroup=sg,
                crystal_system=crystal_system,
                is_maximal=abs(parent_sg - sg.space_group_number) <= 2,
            )
            if query.crystal_system and crystal_system != query.crystal_system.lower().strip():
                continue
            if query.subgroup_space_group and sg.space_group_number != query.subgroup_space_group:
                continue
            if query.maximal_subgroup_only and not item.is_maximal:
                continue
            if query.direct_sublattice and item.direct_sublattice != list(query.direct_sublattice):
                continue
            result.append(item)
        return result

    def method_2_search(self, parent_sg: int, subgroups: Sequence[SubgroupInfo],
                        query: Method2Query) -> Method2Result:
        target = next((s for s in subgroups if s.index == query.subgroup_idx), None)
        if target is None:
            raise ValueError(f"Subgroup index {query.subgroup_idx} not found")

        modes = self._iso.calc_distortion_modes(
            parent_sg,
            query.subgroup_idx,
            distortion_type=query.distortion_type,
        )

        metadata = {
            "k_point_label": query.k_point_label,
            "k_point_coordinates": [
                _to_float(v) for v in query.k_point_coordinates
            ] if query.k_point_coordinates else None,
            "k_parameters": {k: _to_float(v) for k, v in query.k_parameters.items()},
            "number_of_independent_modulations": query.number_of_independent_modulations,
            "number_of_superposed_irs": query.number_of_superposed_irs,
            "specified_opd": query.specified_opd,
        }

        return Method2Result(subgroup=target, modes=modes, metadata=metadata)

    def method_3_search(self, parent_sg: int, query: Method3Query) -> List[Method3ResultItem]:
        distortion_types = normalize_distortion_types(query.distortion_types)
        subgroups = self._iso.list_subgroups(parent_sg, distortion_types)

        if query.point_group and query.space_group_type:
            # Follow ISODISTORT rule: space-group selection supersedes point-group selection.
            point_group_filter = None
        else:
            point_group_filter = query.point_group

        basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        if query.supercell_basis:
            basis = _parse_basis_rows(query.supercell_basis)

        result: List[Method3ResultItem] = []
        for sg in subgroups:
            point_group = _space_group_to_point_group(sg.space_group_number)
            if query.space_group_type and sg.space_group_number != query.space_group_type:
                continue
            if point_group_filter and point_group != point_group_filter:
                continue
            if query.supercell_basis and not _matrix_matches_diagonal_supercell(basis, [1, 1, 1]):
                # Current local backend does not expose full arbitrary-basis subgroup regeneration.
                # Keep entries but mark the requested basis in the result for downstream checks.
                pass

            result.append(Method3ResultItem(subgroup=sg, point_group=point_group, basis=basis))

        return result

    def method_4_decompose(self,
                           parent_structure: Structure,
                           distorted_structure: Structure,
                           mode_displacements: Dict[str, np.ndarray],
                           query: Method4Query) -> Method4Result:
        if len(parent_structure) != len(distorted_structure):
            raise ValueError(
                "Mode decomposition currently requires parent and distorted structures "
                "to have the same atom count"
            )

        assignments = self._match_atoms(parent_structure, distorted_structure, query)
        delta = self._build_delta_vector(parent_structure, distorted_structure, assignments)

        mode_labels = list(mode_displacements.keys())
        if not mode_labels:
            raise ValueError("No mode displacements are available for decomposition")

        columns = []
        for label in mode_labels:
            vec = np.asarray(mode_displacements[label], dtype=float).reshape(-1)
            if vec.size != delta.size:
                raise ValueError(
                    f"Mode {label} has incompatible size {vec.size}, expected {delta.size}"
                )
            columns.append(vec)
        a_matrix = np.column_stack(columns)

        coeffs, _, _, _ = np.linalg.lstsq(a_matrix, delta, rcond=None)
        reconstructed = a_matrix @ coeffs
        residual = delta - reconstructed

        amplitudes = {mode_labels[i]: float(coeffs[i]) for i in range(len(mode_labels))}
        rms = float(np.sqrt(np.mean(residual ** 2)))
        max_abs = float(np.max(np.abs(residual)))

        return Method4Result(
            amplitudes=amplitudes,
            rms_residual=rms,
            max_abs_residual=max_abs,
            assignments=assignments,
            metadata={
                "atom_matching_method": query.atom_matching_method,
                "provided_origin_shift": list(query.provided_origin_shift)
                if query.provided_origin_shift is not None else None,
            },
        )

    @staticmethod
    def _fractional_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        delta = b - a
        delta -= np.round(delta)
        return delta

    def _match_atoms(self,
                     parent_structure: Structure,
                     distorted_structure: Structure,
                     query: Method4Query) -> List[int]:
        if query.atom_matching_method not in {"nearest-site", "robust"}:
            raise ValueError("atom_matching_method must be 'nearest-site' or 'robust'")

        assignments: List[int] = []
        used: set[int] = set()

        for i, site in enumerate(parent_structure):
            species = site.species_string
            parent_coord = np.asarray(site.frac_coords)
            candidates = []
            for j, dst in enumerate(distorted_structure):
                if j in used:
                    continue
                if dst.species_string != species:
                    continue
                delta = self._fractional_delta(parent_coord, np.asarray(dst.frac_coords))
                dist = float(np.linalg.norm(delta))
                if query.atom_matching_method == "nearest-site":
                    candidates.append((dist, j))
                else:
                    if dist <= query.robust_distance_threshold:
                        candidates.append((dist, j))

            if not candidates:
                raise ValueError(f"Cannot match parent atom index {i} ({species})")

            candidates.sort(key=lambda t: t[0])
            chosen = candidates[0][1]
            assignments.append(chosen)
            used.add(chosen)

        return assignments

    def _build_delta_vector(self,
                            parent_structure: Structure,
                            distorted_structure: Structure,
                            assignments: Sequence[int]) -> np.ndarray:
        delta_rows = []
        for i, j in enumerate(assignments):
            p = np.asarray(parent_structure[i].frac_coords)
            d = np.asarray(distorted_structure[j].frac_coords)
            delta_rows.append(self._fractional_delta(p, d))
        return np.asarray(delta_rows, dtype=float).reshape(-1)