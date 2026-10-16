"""Method 4 decomposition: Cartesian least squares, atom matching and strain.

Scientific implementation moved without changing formulas or conventions.
The Method 1 search facade inherits this implementation for existing callers.
"""
from __future__ import annotations

import math
import warnings
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from pymatgen.core import Lattice, Structure
from scipy.optimize import linear_sum_assignment

from backend.utils.config_loader import get_config
from .strain import decompose_homogeneous_strain

@dataclass
class Method4Query:
    """Method 4: mode decomposition of distorted structure."""

    atom_matching_method: str = "nearest-site"
    # Cartesian distance in angstrom.  A fractional-coordinate cutoff would
    # have a different physical meaning for every lattice and crystal axis.
    robust_distance_threshold: float = 0.25
    provided_origin_shift: Sequence[float] | None = None
    # Number of primitive daughter cells represented by the comparison cell.
    # Official As is normalized within one primitive supercell.
    primitive_cell_multiplicity: int = 1
    # Primitive-supercell volume divided by primitive-parent volume.
    supercell_size: float = 1.0
    # Conventional parent lattice and the exact Method-4 basis. Supplying both
    # expresses the official symmetric multiplier in parent-lattice-basis
    # coordinates. Direct low-level callers may omit them; the undistorted
    # comparison cell and an identity basis are then used.
    reference_parent_lattice: Sequence[Sequence[float]] | None = None
    parent_to_child_basis: Sequence[Sequence[float]] | None = None
    strain_mode_labels: Sequence[str] | None = None
    strain_mode_irrep_labels: Sequence[str] | None = None
    strain_mode_q_raw: Sequence[Sequence[float]] | None = None
    strain_mode_q_unit: Sequence[Sequence[float]] | None = None
    strain_mode_normfactors: Sequence[float] | None = None


@dataclass
class Method4Result:
    # Official-style primitive-supercell-normalized As amplitudes (angstrom).
    amplitudes: dict[str, float]
    parent_cell_amplitudes: dict[str, float]
    raw_coefficients: dict[str, float]
    mode_normfactors: dict[str, float]
    rms_residual: float
    max_abs_residual: float
    assignments: list[int]
    strain_mode_amplitudes: dict[str, float]
    strain_modes: list[dict[str, object]]
    strain_raw_coordinate_sum_parent_basis: dict[str, float] | None
    strain_applied_engineering_q_parent_basis: dict[str, float]
    strain_tensor_parent_basis: list[list[float]]
    strain_multiplier_parent_basis: list[list[float]]
    metadata: dict[str, object]

    @property
    def strain_voigt_engineering(self) -> dict[str, float]:
        """Deprecated alias for the same applied parent-basis q object."""

        warnings.warn(
            "strain_voigt_engineering is deprecated; use "
            "strain_applied_engineering_q_parent_basis",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.strain_applied_engineering_q_parent_basis

    @property
    def strain_tensor(self) -> list[list[float]]:
        """Deprecated alias for the same parent-basis tensor object."""

        warnings.warn(
            "strain_tensor is deprecated; use strain_tensor_parent_basis",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.strain_tensor_parent_basis



class ModeDecomposer:
    """Single owner of Method 4 decomposition and its coordinate conventions."""

    def method_4_decompose(self,
                           parent_structure: Structure,
                           distorted_structure: Structure,
                           mode_displacements: dict[str, np.ndarray],
                           query: Method4Query) -> Method4Result:
        if len(parent_structure) != len(distorted_structure):
            raise ValueError(
                "Mode decomposition currently requires parent and distorted structures "
                "to have the same atom count"
            )

        parent_metric = np.asarray(parent_structure.lattice.metric_tensor, dtype=float)
        distorted_metric = np.asarray(
            distorted_structure.lattice.metric_tensor, dtype=float,
        )
        metric_scale = max(float(np.linalg.norm(parent_metric)), 1.0)
        lattice_residual = float(
            np.linalg.norm(distorted_metric - parent_metric) / metric_scale
        )
        lattice_tolerance = float(get_config().lattice_tolerance)
        reference_parent_lattice = np.asarray(
            query.reference_parent_lattice
            if query.reference_parent_lattice is not None
            else parent_structure.lattice.matrix,
            dtype=float,
        )
        parent_to_child_basis = np.asarray(
            query.parent_to_child_basis
            if query.parent_to_child_basis is not None
            else np.eye(3),
            dtype=float,
        )
        if (query.reference_parent_lattice is None) != (
            query.parent_to_child_basis is None
        ):
            raise ValueError(
                "reference_parent_lattice and parent_to_child_basis must be "
                "provided together"
            )
        expected_reference_lattice = parent_to_child_basis @ reference_parent_lattice
        expected_reference_metric = (
            expected_reference_lattice @ expected_reference_lattice.T
        )
        reference_metric_scale = max(
            float(np.linalg.norm(expected_reference_metric)), 1.0
        )
        reference_metric_residual = float(
            np.linalg.norm(parent_metric - expected_reference_metric)
            / reference_metric_scale
        )
        if reference_metric_residual > lattice_tolerance:
            raise ValueError(
                "Selected Method 4 basis is incompatible with the undistorted "
                f"reference cell (relative metric residual "
                f"{reference_metric_residual:.6g} > {lattice_tolerance:.6g})"
            )
        strain = decompose_homogeneous_strain(
            reference_parent_lattice,
            parent_to_child_basis,
            distorted_structure.lattice.matrix,
        )

        assignments = self._match_atoms(parent_structure, distorted_structure, query)
        delta_fractional = self._build_delta_vector(
            parent_structure,
            distorted_structure,
            assignments,
            query,
        ).reshape((-1, 3))
        lattice = np.asarray(parent_structure.lattice.matrix, dtype=float)
        # Fit and report residuals in Cartesian angstroms.  This makes the
        # least-squares metric invariant to the choice and scale of fractional
        # axes, while retaining the exact generation coefficients.
        delta = (delta_fractional @ lattice).reshape(-1)

        mode_labels = list(mode_displacements.keys())
        if not mode_labels:
            raise ValueError("No mode displacements are available for decomposition")

        columns = []
        for label in mode_labels:
            fractional = np.asarray(mode_displacements[label], dtype=float)
            expected_shape = (len(parent_structure), 3)
            if fractional.shape != expected_shape:
                raise ValueError(
                    f"Mode {label} has incompatible shape {fractional.shape}, "
                    f"expected {expected_shape}"
                )
            if not bool(np.all(np.isfinite(fractional))):
                raise ValueError(f"Mode {label} contains non-finite displacements")
            vec = (fractional @ lattice).reshape(-1)
            columns.append(vec)
        a_matrix = np.column_stack(columns)

        coeffs, _, rank, singular_values = np.linalg.lstsq(
            a_matrix,
            delta,
            rcond=None,
        )
        if int(rank) < len(mode_labels):
            raise ValueError(
                "Method 4 mode basis is linearly dependent: "
                f"rank {rank} for {len(mode_labels)} modes; amplitudes are not unique"
            )
        reconstructed = a_matrix @ coeffs
        residual = delta - reconstructed

        primitive_multiplicity = int(query.primitive_cell_multiplicity)
        if primitive_multiplicity <= 0:
            raise ValueError("primitive_cell_multiplicity must be a positive integer")
        supercell_size = float(query.supercell_size)
        if not math.isfinite(supercell_size) or supercell_size <= 0:
            raise ValueError("supercell_size must be a positive finite value")
        # The fitted columns retain the mapper's internal fractional scale so
        # ``raw_coefficients`` reconstruct the input exactly.  Official
        # complete-mode tables use a different presentation coordinate: every
        # child-fractional column is first divided by its largest absolute
        # component.  Keep those two coordinates explicit.  If
        # ``q_canonical = q_internal / scale``, then
        #
        #   c_internal q_internal = (c_internal scale) q_canonical
        #                         = As normfactor q_canonical.
        #
        # A conventional centered cell contains ``primitive_multiplicity``
        # identical primitive cells, so its Cartesian RSS norm is divided by
        # sqrt(multiplicity) before constructing the official normfactor.
        conventional_norms = np.linalg.norm(a_matrix, axis=0)
        primitive_norms = conventional_norms / math.sqrt(primitive_multiplicity)
        if bool(np.any(primitive_norms <= 0)):
            raise ValueError("Method 4 mode basis contains a zero-norm vector")
        canonical_scales = np.asarray([
            float(np.max(np.abs(np.asarray(mode_displacements[label], dtype=float))))
            for label in mode_labels
        ])
        if bool(np.any(canonical_scales <= 0)):
            raise ValueError("Method 4 mode basis contains a zero fractional column")
        internal_normfactors = 1.0 / primitive_norms
        normfactors_array = canonical_scales * internal_normfactors
        raw_coefficients = {
            mode_labels[i]: float(coeffs[i]) for i in range(len(mode_labels))
        }
        canonical_raw_coefficients = {
            mode_labels[i]: float(coeffs[i] * canonical_scales[i])
            for i in range(len(mode_labels))
        }
        amplitudes = {
            mode_labels[i]: float(coeffs[i] / internal_normfactors[i])
            for i in range(len(mode_labels))
        }
        parent_cell_amplitudes = {
            label: float(value / math.sqrt(supercell_size))
            for label, value in amplitudes.items()
        }
        mode_normfactors = {
            mode_labels[i]: float(normfactors_array[i])
            for i in range(len(mode_labels))
        }
        rms = float(np.sqrt(np.mean(residual ** 2)))
        max_abs = float(np.max(np.abs(residual)))
        condition_number = (
            float(singular_values[0] / singular_values[-1])
            if len(singular_values) and singular_values[-1] > 0
            else float("inf")
        )

        strain_mode_labels = list(query.strain_mode_labels or [])
        strain_irreps = list(query.strain_mode_irrep_labels or [])
        strain_raw_rows = list(query.strain_mode_q_raw or [])
        strain_unit_rows = list(query.strain_mode_q_unit or [])
        strain_normfactors = list(query.strain_mode_normfactors or [])
        strain_mode_amplitudes: dict[str, float] = {}
        strain_mode_definitions: list[dict[str, object]] = []
        strain_raw_coordinate_sum: dict[str, float] | None = None
        strain_mode_status = "unresolved_missing_canonical_basis"
        supplied_lengths = {
            len(strain_mode_labels),
            len(strain_irreps),
            len(strain_raw_rows),
            len(strain_unit_rows),
            len(strain_normfactors),
        }
        if supplied_lengths != {0}:
            if len(supplied_lengths) != 1 or not strain_mode_labels:
                raise ValueError(
                    "Method 4 canonical strain mode fields must be supplied "
                    "together with equal nonzero lengths"
                )
            raw_matrix = np.column_stack(
                [np.asarray(row, dtype=float) for row in strain_raw_rows]
            )
            unit_matrix = np.column_stack(
                [np.asarray(row, dtype=float) for row in strain_unit_rows]
            )
            if raw_matrix.shape[0] != 6 or unit_matrix.shape != raw_matrix.shape:
                raise ValueError(
                    "Method 4 canonical strain modes must contain six-component rows"
                )
            if not np.all(np.isfinite(raw_matrix)) or not np.all(
                np.isfinite(unit_matrix)
            ):
                raise ValueError("Method 4 canonical strain modes must be finite")
            for column, normfactor in enumerate(strain_normfactors):
                if not math.isfinite(float(normfactor)) or float(normfactor) <= 0:
                    raise ValueError("Method 4 strain normfactors must be positive")
                if not np.allclose(
                    unit_matrix[:, column],
                    float(normfactor) * raw_matrix[:, column],
                    rtol=5.0e-5,
                    atol=5.0e-5,
                ):
                    raise ValueError(
                        "Method 4 q_unit must equal normfactor*q_raw"
                    )
            strain_coefficients, _, strain_rank, _ = np.linalg.lstsq(
                unit_matrix,
                strain.applied_engineering_q_parent_basis,
                rcond=None,
            )
            if int(strain_rank) < len(strain_mode_labels):
                raise ValueError("Method 4 canonical strain basis is dependent")
            strain_fit_residual = float(
                np.linalg.norm(
                    unit_matrix @ strain_coefficients
                    - strain.applied_engineering_q_parent_basis
                )
            )
            if strain_fit_residual > max(5.0e-7, lattice_tolerance * 10.0):
                raise ValueError(
                    "Method 4 applied strain is outside the selected subgroup's "
                    "canonical strain fixed space"
                )
            strain_mode_amplitudes = {
                label: float(strain_coefficients[index])
                for index, label in enumerate(strain_mode_labels)
            }
            raw_sum = raw_matrix @ strain_coefficients
            strain_raw_coordinate_sum = {
                label: float(value)
                for label, value in zip(strain.voigt_order, raw_sum, strict=True)
            }
            strain_mode_definitions = [
                {
                    "label": label,
                    "irrep_label": strain_irreps[index],
                    "q_raw": [float(value) for value in raw_matrix[:, index]],
                    "q_unit": [float(value) for value in unit_matrix[:, index]],
                    "normfactor": float(strain_normfactors[index]),
                }
                for index, label in enumerate(strain_mode_labels)
            ]
            strain_mode_status = "canonical_iso_rank_12"

        return Method4Result(
            amplitudes=amplitudes,
            parent_cell_amplitudes=parent_cell_amplitudes,
            raw_coefficients=raw_coefficients,
            mode_normfactors=mode_normfactors,
            rms_residual=rms,
            max_abs_residual=max_abs,
            assignments=assignments,
            strain_mode_amplitudes=strain_mode_amplitudes,
            strain_modes=strain_mode_definitions,
            strain_raw_coordinate_sum_parent_basis=strain_raw_coordinate_sum,
            strain_applied_engineering_q_parent_basis={
                label: float(value)
                for label, value in zip(
                    strain.voigt_order,
                    strain.applied_engineering_q_parent_basis,
                    strict=True,
                )
            },
            strain_tensor_parent_basis=[
                [float(value) for value in row]
                for row in strain.strain_tensor_parent_basis
            ],
            strain_multiplier_parent_basis=[
                [float(value) for value in row]
                for row in strain.multiplier_parent_basis
            ],
            metadata={
                "atom_matching_method": query.atom_matching_method,
                "provided_origin_shift": list(query.provided_origin_shift)
                if query.provided_origin_shift is not None else None,
                "residual_unit": "angstrom",
                "lattice_metric_relative_residual": lattice_residual,
                "reference_lattice_metric_relative_residual": (
                    reference_metric_residual
                ),
                "strain_reconstruction_relative_metric_residual": (
                    strain.relative_metric_residual
                ),
                "strain_coordinate_frame": strain.coordinate_frame,
                "strain_voigt_order": list(strain.voigt_order),
                "strain_column_lattice_action": strain.column_lattice_action,
                "strain_row_lattice_action": strain.row_lattice_action,
                "strain_child_lattice_action": strain.child_lattice_action,
                "strain_metric_equation": strain.metric_equation,
                "strain_convention": (
                    "parent-lattice-basis engineering q order "
                    "(E11,E22,E33,2E23,2E13,2E12); M=I+E; "
                    "child row lattice B@M@P"
                ),
                "strain_mode_status": strain_mode_status,
                "deprecated_result_aliases": {
                    "strain_voigt_engineering": (
                        "strain_applied_engineering_q_parent_basis"
                    ),
                    "strain_tensor": "strain_tensor_parent_basis",
                },
                "mode_matrix_rank": int(rank),
                "mode_matrix_columns": len(mode_labels),
                "mode_matrix_condition_number": condition_number,
                "raw_coefficient_convention": (
                    "coefficient of the installed internal child-fractional column"
                ),
                "canonical_mode_fractional_scale": {
                    mode_labels[i]: float(canonical_scales[i])
                    for i in range(len(mode_labels))
                },
                "canonical_raw_coefficients": canonical_raw_coefficients,
                "canonical_mode_convention": (
                    "q_canonical = q_internal / max(abs(q_internal))"
                ),
                "mode_normfactor_convention": (
                    "inverse Cartesian RSS of q_canonical over one primitive child cell"
                ),
                "amplitude_unit": "angstrom",
                "amplitude_convention": "As primitive-supercell normalized",
                "parent_cell_amplitude_convention": "Ap = As / sqrt(s)",
                "primitive_cell_multiplicity": primitive_multiplicity,
                "supercell_size": supercell_size,
            },
        )


    @staticmethod
    def _fractional_delta(a: np.ndarray, b: np.ndarray,
                          lattice: Lattice) -> np.ndarray:
        """Return the shortest displacement in the reference cell's metric.

        For row lattice L this minimizes ||(b + n - a) @ L|| over integer
        lattice images n. Componentwise rounding only solves that problem
        for orthogonal axes. The image translates b, as in pymatgen's atom
        assignment distances.
        """
        distance, image = lattice.get_distance_and_image(a, b)
        delta = b - a
        rounded = delta - np.round(delta)
        # At equal-distance boundaries, retain the previous signed branch
        # (including round-to-even half-cell ties). Only permit differences
        # at the scale of floating-point error in a three-component norm.
        rounded_distance = float(np.linalg.norm(rounded @ lattice.matrix))
        if math.isclose(
            rounded_distance, float(distance),
            rel_tol=16.0 * np.finfo(float).eps, abs_tol=0.0,
        ):
            return rounded
        return delta + image


    def _match_atoms(self,
                     parent_structure: Structure,
                     distorted_structure: Structure,
                     query: Method4Query) -> list[int]:
        if query.atom_matching_method not in {"nearest-site", "robust"}:
            raise ValueError("atom_matching_method must be 'nearest-site' or 'robust'")

        threshold = float(query.robust_distance_threshold)
        if not math.isfinite(threshold) or threshold <= 0:
            raise ValueError("robust_distance_threshold must be a positive finite angstrom value")
        shift = self._origin_shift(query)
        distorted_coords = np.asarray(distorted_structure.frac_coords, dtype=float) - shift
        parent_coords = np.asarray(parent_structure.frac_coords, dtype=float)
        assignments = [-1] * len(parent_structure)
        parent_species = [site.species_string for site in parent_structure]
        distorted_species = [site.species_string for site in distorted_structure]

        for species in sorted(set(parent_species) | set(distorted_species)):
            parent_indices = [
                index for index, value in enumerate(parent_species) if value == species
            ]
            distorted_indices = [
                index for index, value in enumerate(distorted_species) if value == species
            ]
            if len(parent_indices) != len(distorted_indices):
                raise ValueError(
                    f"Cannot match species {species}: parent has {len(parent_indices)} "
                    f"site(s), distorted structure has {len(distorted_indices)}"
                )
            if not parent_indices:
                continue
            distances = parent_structure.lattice.get_all_distances(
                parent_coords[parent_indices],
                distorted_coords[distorted_indices],
            )
            rows, columns = linear_sum_assignment(distances)
            for row, column in zip(rows, columns, strict=True):
                distance = float(distances[row, column])
                parent_index = parent_indices[int(row)]
                if query.atom_matching_method == "robust" and distance > threshold:
                    raise ValueError(
                        f"Cannot match parent atom index {parent_index} ({species}): "
                        f"best global assignment distance {distance:.6g} angstrom exceeds "
                        f"threshold {threshold:.6g} angstrom"
                    )
                assignments[parent_index] = distorted_indices[int(column)]

        if any(index < 0 for index in assignments):
            raise ValueError("Cannot construct a complete species-preserving atom assignment")
        return assignments


    @staticmethod
    def _origin_shift(query: Method4Query) -> np.ndarray:
        if query.provided_origin_shift is None:
            return np.zeros(3, dtype=float)
        shift = np.asarray(query.provided_origin_shift, dtype=float)
        if shift.shape != (3,) or not bool(np.all(np.isfinite(shift))):
            raise ValueError("provided_origin_shift must contain three finite fractional values")
        return shift


    def _build_delta_vector(self,
                            parent_structure: Structure,
                            distorted_structure: Structure,
                            assignments: Sequence[int],
                            query: Method4Query) -> np.ndarray:
        shift = self._origin_shift(query)
        delta_rows = []
        for i, j in enumerate(assignments):
            p = np.asarray(parent_structure[i].frac_coords)
            d = np.asarray(distorted_structure[j].frac_coords) - shift
            # The daughter fractional axes refer to the undistorted comparison
            # cell. Homogeneous strain is decomposed separately; its final
            # metric must not change the displacement branch or mode norm.
            delta_rows.append(self._fractional_delta(p, d, parent_structure.lattice))
        return np.asarray(delta_rows, dtype=float).reshape(-1)
