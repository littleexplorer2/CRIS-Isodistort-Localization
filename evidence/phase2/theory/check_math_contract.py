"""Independent finite algebra examples; no production imports or website answers."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BOUND = 512 * np.finfo(float).eps
CHECKS = []


def equal(name, actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    residual = float(np.max(np.abs(actual - expected)))
    if residual > BOUND:
        raise AssertionError((name, residual, BOUND))
    CHECKS.append({"name": name, "max_absolute_residual": residual,
                   "passed": True})


def nonzero(name, residual, lower_bound):
    residual = float(residual)
    if residual <= lower_bound:
        raise AssertionError((name, residual, lower_bound))
    CHECKS.append({"name": name, "counterexample_residual": residual,
                   "required_lower_bound": lower_bound, "passed": True})


def kelvin(matrix):
    return np.array([matrix[0, 0], matrix[1, 1], matrix[2, 2],
                     np.sqrt(2) * matrix[1, 2], np.sqrt(2) * matrix[0, 2],
                     np.sqrt(2) * matrix[0, 1]])


def main():
    # H={e}: every vector is fixed, but a rotated basis mixes two parent IRs.
    parent_op = np.diag([1.0, -1.0])
    mixed = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2)
    equal("mixed_basis_is_complete_orthonormal_H_fixed", mixed.T @ mixed, np.eye(2))
    nonzero("H_fixed_does_not_establish_parent_IR_eigenvector",
            np.linalg.norm(parent_op @ mixed[:, 0] - mixed[:, 0]), 1)

    # Two repeated conventional-cell representatives, one primitive DOF.
    swap = np.array([[0., 1.], [1., 0.]])
    projector = (np.eye(2) + swap) / 2
    equal("centering_projector_idempotence", projector @ projector, projector)
    equal("centering_projector_dimension", np.trace(projector), 1)
    repeated = np.ones(2)
    equal("primitive_norm_uses_centering_weight", repeated @ (.5*np.eye(2)) @ repeated, 1)
    nonzero("unweighted_conventional_norm_is_wrong", repeated @ repeated - 1, .9)

    # A skew column basis: fractional Euclidean normalization is not Cartesian.
    cell = np.array([[2., 1., 0.], [0., 3., 0.], [0., 0., 4.]])
    metric = cell.T @ cell
    nonzero("fractional_axes_are_not_metric_orthogonal", abs(metric[0, 1]), 1)
    modes = np.linalg.inv(cell)
    equal("whitened_fractional_basis_metric_orthogonality", modes.T @ metric @ modes, np.eye(3))
    amplitudes = np.array([.12, -.08, .03])
    displacement = modes @ amplitudes
    equal("displacement_amplitude_roundtrip", modes.T @ metric @ displacement, amplitudes)
    equal("amplitude_RSS_equals_reference_norm",
          amplitudes @ amplitudes, displacement @ metric @ displacement)
    change = np.array([[1., 1., 0.], [0., 1., 0.], [0., 0., 1.]])
    cell_changed = cell @ change
    modes_changed = np.linalg.solve(change, modes)
    equal("non_diagonal_basis_covariance",
          cell_changed @ modes_changed, cell @ modes)
    equal("non_diagonal_basis_norm_covariance",
          modes_changed.T @ (cell_changed.T @ cell_changed) @ modes_changed, np.eye(3))

    # A child translation lattice need not support the full parent point group.
    rotate2 = np.array([[0., -1.], [1., 0.]])
    child = np.diag([2., 1.])
    induced = np.linalg.solve(child, rotate2 @ child)
    nonzero("non_normal_supercell_full_parent_action_fails",
            np.max(np.abs(induced - np.rint(induced))), .4)
    core = 2*np.eye(2)
    equal("normal_core_allows_parent_action", np.linalg.solve(core, rotate2 @ core), rotate2)

    # H=C4 in fixed Cartesian axes: one uniform translation, two strains.
    rotate = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    group = [np.linalg.matrix_power(rotate, i) for i in range(4)]
    acoustic_projector = sum(group) / 4
    equal("fixed_origin_acoustic_dimension_is_one_not_three",
          np.trace(acoustic_projector), 1)
    equal("acoustic_axis_is_z", acoustic_projector, np.diag([0., 0., 1.]))
    tensor_basis = [np.diag(row) for row in np.eye(3)]
    for i, j in [(1, 2), (0, 2), (0, 1)]:
        tensor = np.zeros((3, 3))
        tensor[i, j] = tensor[j, i] = 1/np.sqrt(2)
        tensor_basis.append(tensor)
    operators = [np.column_stack([kelvin(q @ e @ q.T) for e in tensor_basis])
                 for q in group]
    strain_projector = sum(operators) / 4
    expected = np.zeros((6, 6))
    expected[:2, :2] = .5
    expected[2, 2] = 1
    equal("C4_strain_projector_exact_form", strain_projector, expected)
    equal("C4_strain_dimension_two", np.trace(strain_projector), 2)
    equal("C4_strain_projector_idempotence", strain_projector @ strain_projector, strain_projector)
    for i in range(4):
        for j in range(4):
            equal(f"strain_representation_product_{i}_{j}",
                  operators[i] @ operators[j], operators[(i+j) % 4])
    e1 = np.diag([1., 1., 0.]) / np.sqrt(2)
    e2 = np.diag([0., 0., 1.])
    e = .03 * e1 - .01 * e2
    equal("strain_amplitude_Frobenius_norm", np.sum(e*e), .03**2 + .01**2)
    equal("Kelvin_equals_Frobenius_shear_metric",
          kelvin(tensor_basis[3]) @ kelvin(tensor_basis[3]),
          np.sum(tensor_basis[3]**2))
    nonzero("engineering_shear_cannot_use_unweighted_tensor_metric",
            abs(1 - (tensor_basis[3][1, 2]**2)), .4)
    stretch = np.eye(3) + e
    nonzero("finite_symmetric_stretch_positive", np.min(np.linalg.eigvalsh(stretch)), .9)
    strained_cell = stretch @ cell
    equal("strain_cell_roundtrip", strained_cell @ np.linalg.inv(cell) - np.eye(3), e)
    for i, q in enumerate(group):
        equal(f"strain_H_invariance_{i}", q @ e @ q.T, e)

    record = {
        "scope": "Finite mathematical examples of MODE_CONTRACT; not production/scientific end-to-end acceptance",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "numpy_version": np.__version__,
        "source_imports": ["Python standard library", "numpy"],
        "website_case_results_read": False,
        "production_modules_imported": False,
        "tolerance": {
            "absolute_bound": float(BOUND),
            "derivation": "512*binary64 epsilon for these <=6-dimensional, well-conditioned explicit toy matrices only; not a production rank threshold",
        },
        "checks": CHECKS,
    }
    (HERE / "math_checks.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"check_count": len(CHECKS), "all_passed": True,
                      "scope": record["scope"], "absolute_bound": float(BOUND)}))


if __name__ == "__main__":
    main()
