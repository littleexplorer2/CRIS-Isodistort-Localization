"""Homogeneous-strain modes in ISODISTORT parent-lattice coordinates.

The parent direct lattice is stored as a row matrix ``P``. ISODISTORT uses a
symmetric multiplier ``M = I + E`` in the *parent lattice basis*: the same
operation is ``C -> C @ M`` for the column-basis convention ``C = P.T`` and
``P -> M @ P`` in the row convention used here. A child basis ``B`` therefore
reconstructs as ``B @ M @ P``.

The six coordinates are the engineering-Voigt row

``q = (E11, E22, E33, 2 E23, 2 E13, 2 E12)``.

This factor of two is fixed by official ISODISTORT files: a raw unit shear has
``normfactor=sqrt(2)`` and its IsoVIZ vector is ``sqrt(2)``. Thus
``q.T @ W @ q`` with ``W=diag(1,1,1,1/2,1/2,1/2)`` is the Frobenius norm of
``E``. The coordinates are not Cartesian components when ``P`` is oblique.

For a fractional-column rotation ``R``, ``R.T @ G @ R = G`` where
``G=P@P.T``. Linearizing ``M @ G @ M`` gives ``dG=E@G+G@E``. A strain mode
preserves the embedded subgroup to first order when
``R.T @ dG @ R = dG`` for every subgroup rotation. The implementation builds
these constraints from the actual embedded operations and actual parent
lattice; it does not dispatch on a space-group number or crystal system.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

ENGINEERING_VOIGT_ORDER = ("11", "22", "33", "2*23", "2*13", "2*12")
ENGINEERING_VOIGT_METRIC = np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])
PARENT_BASIS_COORDINATE_FRAME = "parent_lattice_basis"
COLUMN_LATTICE_ACTION = "C_deformed = C_parent @ M"
ROW_LATTICE_ACTION = "P_deformed = M @ P_parent"
CHILD_LATTICE_ACTION = "L_child = B @ M @ P_parent"
METRIC_EQUATION = "M @ G_parent @ M = H_parent_basis"
_SQRT_VOIGT_METRIC = np.diag(
    [1.0, 1.0, 1.0, 1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)]
)
_INV_SQRT_VOIGT_METRIC = np.diag(
    [1.0, 1.0, 1.0, np.sqrt(2.0), np.sqrt(2.0), np.sqrt(2.0)]
)


@dataclass(frozen=True, slots=True)
class HomogeneousStrainMode:
    """One deterministic basis vector of the embedded strain fixed space.

    ``q_raw`` is the max-component-one vector written to the CIF conversion
    matrix. ``q_unit = normfactor * q_raw`` is its unit-Frobenius form written
    to IsoVIZ. A mode amplitude multiplies ``q_unit`` when reconstructing the
    applied multiplier ``M``.
    """

    index: int
    q_unit: np.ndarray
    q_raw: np.ndarray
    normfactor: float
    irrep_label: str | None
    label_status: str
    invariance_max_relative_metric_error: float
    unit_norm_error: float
    raw_max_component_error: float
    canonical_label: str | None = None
    irrep_direction: str | None = None


@dataclass(frozen=True, slots=True)
class HomogeneousStrainModeDiagnostics:
    """Numerical proof obligations for the parent-basis fixed space."""

    input_rotation_count: int
    point_group_order: int
    fixed_dimension: int
    rotation_metric_preservation_max_relative_error: float
    rotation_determinant_max_error: float
    group_identity_error: float
    group_closure_max_error: float
    representation_identity_error: float
    representation_closure_max_error: float
    fixed_projector_idempotence_max_error: float
    fixed_projector_constraint_max_error: float
    invariance_max_relative_metric_error: float
    orthonormality_max_abs_error: float


@dataclass(frozen=True, slots=True)
class HomogeneousStrainModeResult:
    """Complete parent-basis strain fixed space and diagnostics."""

    modes: tuple[HomogeneousStrainMode, ...]
    diagnostics: HomogeneousStrainModeDiagnostics
    projector_engineering: np.ndarray
    projector_orthonormal: np.ndarray
    parent_metric: np.ndarray
    voigt_order: tuple[str, ...] = ENGINEERING_VOIGT_ORDER
    coordinate_frame: str = PARENT_BASIS_COORDINATE_FRAME
    column_lattice_action: str = COLUMN_LATTICE_ACTION
    row_lattice_action: str = ROW_LATTICE_ACTION
    child_lattice_action: str = CHILD_LATTICE_ACTION
    metric_equation: str = METRIC_EQUATION
    basis_source: str = "metric_fixed_space_fallback"
    export_ready: bool = False
    fixed_space_validation_max_error: float = 0.0
    irrep_labels_available: bool = False
    irrep_label_reason: str = (
        "unresolved: canonical parent-Gamma macroscopic basis was not supplied"
    )


@dataclass(frozen=True, slots=True)
class CanonicalStrainModeDefinition:
    """One ISO ``DISPLAY DISTORTION`` rank-[12] macroscopic basis vector."""

    label: str
    q_raw: np.ndarray
    irrep_label: str
    irrep_direction: str = ""


def _invariant_direction_parameter_matrix(
    direction_raw: str, component_count: int
) -> tuple[tuple[str, ...], np.ndarray]:
    """Resolve a linear ISO direction such as ``(a,0)`` or ``(a,a)``."""

    raw = str(direction_raw).strip()
    if not (raw.startswith("(") and raw.endswith(")")):
        raise ValueError("ISO invariant direction must be parenthesized")
    body = raw[1:-1]
    if ";" in body:
        raise ValueError("Gamma strain direction unexpectedly contains star arms")
    coordinates = [value.strip() for value in body.split(",")]
    if len(coordinates) != component_count:
        raise ValueError(
            "ISO invariant direction dimension does not match macroscopic irrep"
        )
    parameter_order: list[str] = []
    rows: list[dict[str, Fraction]] = []
    term_pattern = re.compile(
        r"(?:(?P<numerator>(?:\d+(?:\.\d*)?|\.\d+)(?:/\d+)?)\*?)?"
        r"(?P<parameter>[A-Za-z][A-Za-z0-9_]*)"
        r"(?:/(?P<denominator>\d+))?"
    )
    for coordinate in coordinates:
        text = re.sub(r"\s+", "", coordinate)
        if text in {"", "0", "+0", "-0"}:
            rows.append({})
            continue
        signed = text if text[0] in "+-" else f"+{text}"
        terms = re.findall(r"[+-][^+-]+", signed)
        if not terms or "".join(terms) != signed:
            raise ValueError("unsupported nonlinear ISO invariant direction")
        coefficients: dict[str, Fraction] = {}
        for term in terms:
            sign = -1 if term[0] == "-" else 1
            match = term_pattern.fullmatch(term[1:])
            if match is None:
                raise ValueError("unsupported ISO invariant-direction term")
            parameter = match.group("parameter")
            # Fraction parses both rational tokens and finite decimal strings
            # exactly, so an ISO-printed approximation such as ``1.732`` is
            # retained as printed.  The independent fixed-space comparison,
            # rather than an assumed radical, decides whether it is usable.
            numerator = Fraction(match.group("numerator") or 1)
            denominator = int(match.group("denominator") or 1)
            coefficient = sign * numerator / denominator
            coefficients[parameter] = coefficients.get(
                parameter, Fraction(0)
            ) + coefficient
            if parameter not in parameter_order:
                parameter_order.append(parameter)
        rows.append(coefficients)
    if not parameter_order:
        return (), np.zeros((component_count, 0), dtype=float)
    matrix = np.array(
        [
            [float(row.get(parameter, Fraction(0))) for parameter in parameter_order]
            for row in rows
        ],
        dtype=float,
    )
    return tuple(parameter_order), matrix


def canonical_strain_definitions_from_iso(
    macroscopic_basis: object,
) -> tuple[CanonicalStrainModeDefinition, ...]:
    """Apply exact-embedding invariant directions to ISO rank-[12] blocks.

    The wrapper deliberately keeps the macro basis and ``DISPLAY DIRECTION``
    evidence separate. This function combines them without guessing irrep
    labels or direction components. Printed decimal tensor coefficients remain
    floats here and are subsequently accepted only if the whole-space metric
    check in :func:`apply_canonical_strain_basis` passes.
    """

    if getattr(macroscopic_basis, "status", "unresolved") == "unresolved":
        reason = getattr(macroscopic_basis, "reason", None) or "unknown reason"
        raise ValueError(f"ISO macroscopic strain basis is unresolved: {reason}")
    if str(getattr(macroscopic_basis, "rank_signature", "")) != "[12]":
        raise ValueError("ISO macroscopic strain basis must use rank [12]")
    coefficient_basis = tuple(
        str(value) for value in getattr(macroscopic_basis, "coefficient_basis", ())
    )
    if coefficient_basis != ("xx", "yy", "zz", "yz", "xz", "xy"):
        raise ValueError("ISO macroscopic strain basis uses an unknown component order")

    directions_by_irrep: dict[str, list[object]] = {}
    for direction in getattr(macroscopic_basis, "invariant_directions", ()):
        directions_by_irrep.setdefault(str(direction.irrep_label), []).append(
            direction
        )
    blocks = tuple(getattr(macroscopic_basis, "blocks", ()))
    multiplicities: dict[str, int] = {}
    for block in blocks:
        label = str(block.global_irrep)
        multiplicities[label] = multiplicities.get(label, 0) + 1

    definitions: list[CanonicalStrainModeDefinition] = []
    for block in blocks:
        irrep = str(block.global_irrep)
        matching = directions_by_irrep.get(irrep, [])
        if not matching:
            continue
        unique_directions = {str(item.direction_raw) for item in matching}
        if len(unique_directions) != 1:
            raise ValueError(
                f"ISO returned multiple invariant directions for strain irrep {irrep}"
            )
        invariant_direction = next(iter(unique_directions))
        components = tuple(block.components)
        parameters, direction_matrix = _invariant_direction_parameter_matrix(
            invariant_direction, len(components)
        )
        component_matrix = []
        for component in components:
            vector = component.as_float_vector()
            if vector is None or np.asarray(vector).shape != (6,):
                raise ValueError(
                    f"ISO strain component for {irrep} is numerically unresolved"
                )
            component_matrix.append(np.asarray(vector, dtype=float))
        tensor_matrix = np.column_stack(component_matrix)
        invariant_vectors = tensor_matrix @ direction_matrix
        copy_suffix = (
            f"_{int(block.copy_index)}"
            if multiplicities[irrep] > 1 and block.copy_index is not None
            else ""
        )
        for column, parameter in enumerate(parameters):
            q_raw = np.array(invariant_vectors[:, column], copy=True)
            sign_threshold = max(
                1.0e-12,
                float(np.max(np.abs(q_raw), initial=0.0)) * 1.0e-10,
            )
            first_nonzero = next(
                (value for value in q_raw if abs(float(value)) > sign_threshold),
                0.0,
            )
            # ISODISTORT's final one-dimensional mode convention makes the
            # first nonzero displayed coefficient positive. For SG139 GM5+
            # this turns ISO's parent macro component ``-xz`` into the +xz
            # vector archived in Complete modes, CIF, and IsoVIZ.
            if first_nonzero < 0.0:
                q_raw *= -1.0
            maximum_component = float(np.max(np.abs(q_raw), initial=0.0))
            if maximum_component <= sign_threshold:
                raise ValueError(
                    f"ISO invariant direction for strain irrep {irrep} is zero"
                )
            # Website CIF conversion columns are max-component-one raw
            # coordinates.  This step is essential for directions such as
            # cubic ``xx+yy-2zz``: its official q_raw is (.5,.5,-1), while
            # q_unit retains the same physical direction after normalization.
            q_raw /= maximum_component
            definitions.append(
                CanonicalStrainModeDefinition(
                    label=f"{irrep}strain{copy_suffix}({parameter})",
                    q_raw=q_raw,
                    irrep_label=irrep,
                    irrep_direction=invariant_direction,
                )
            )
    return tuple(definitions)


def _matrix3(values: object, description: str) -> np.ndarray:
    raw = getattr(values, "rotation_matrix", values)
    matrix = np.asarray(raw, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"{description} must be a 3x3 matrix")
    if not bool(np.all(np.isfinite(matrix))):
        raise ValueError(f"{description} must contain only finite values")
    return matrix


def engineering_voigt_to_tensor(values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Convert official engineering coordinates to symmetric parent-basis ``E``."""

    q = np.asarray(values, dtype=float)
    if q.shape != (6,):
        raise ValueError("engineering Voigt strain must contain six values")
    if not bool(np.all(np.isfinite(q))):
        raise ValueError("engineering Voigt strain must contain only finite values")
    return np.array(
        [
            [q[0], q[5] / 2.0, q[4] / 2.0],
            [q[5] / 2.0, q[1], q[3] / 2.0],
            [q[4] / 2.0, q[3] / 2.0, q[2]],
        ],
        dtype=float,
    )


def tensor_to_engineering_voigt(
    values: object, *, symmetry_tolerance: float = 1.0e-10
) -> np.ndarray:
    """Convert symmetric parent-basis ``E`` to official engineering order."""

    tensor = _matrix3(values, "strain tensor")
    antisymmetric = tensor - tensor.T
    if float(np.max(np.abs(antisymmetric))) > symmetry_tolerance:
        raise ValueError("strain tensor must be symmetric")
    tensor = (tensor + tensor.T) / 2.0
    return np.array(
        [
            tensor[0, 0],
            tensor[1, 1],
            tensor[2, 2],
            2.0 * tensor[1, 2],
            2.0 * tensor[0, 2],
            2.0 * tensor[0, 1],
        ],
        dtype=float,
    )


def parent_basis_metric_perturbation(
    q: Sequence[float] | np.ndarray, parent_metric: object
) -> np.ndarray:
    """Return ``dG = E @ G + G @ E`` for parent-basis engineering ``q``."""

    metric = _matrix3(parent_metric, "parent metric")
    strain = engineering_voigt_to_tensor(q)
    return strain @ metric + metric @ strain


def _metric_congruence_action(rotation: np.ndarray) -> np.ndarray:
    """Return the six-dimensional action ``S -> R.T @ S @ R``."""

    columns: list[np.ndarray] = []
    for index in range(6):
        q = np.zeros(6, dtype=float)
        q[index] = 1.0
        symmetric = engineering_voigt_to_tensor(q)
        columns.append(tensor_to_engineering_voigt(rotation.T @ symmetric @ rotation))
    return np.column_stack(columns)


def _metric_perturbation_map(parent_metric: np.ndarray) -> np.ndarray:
    columns = []
    for index in range(6):
        q = np.zeros(6, dtype=float)
        q[index] = 1.0
        columns.append(
            tensor_to_engineering_voigt(
                parent_basis_metric_perturbation(q, parent_metric)
            )
        )
    return np.column_stack(columns)


def parent_basis_strain_action(
    rotation_fractional: object, parent_lattice: object
) -> np.ndarray:
    """Return the induced action on parent-basis engineering strain.

    It applies ``R.T @ dG @ R`` and then solves the invertible Sylvester map
    ``E @ G + G @ E = dG``. The action is generally not orthogonal in the
    simple engineering-Voigt metric for an oblique lattice.
    """

    rotation = _matrix3(rotation_fractional, "fractional rotation")
    parent = _matrix3(parent_lattice, "parent lattice")
    if abs(float(np.linalg.det(parent))) <= np.finfo(float).eps:
        raise ValueError("parent lattice must be nonsingular")
    metric = parent @ parent.T
    sylvester = _metric_perturbation_map(metric)
    congruence = _metric_congruence_action(rotation)
    return np.linalg.solve(sylvester, congruence @ sylvester)


def _validated_fractional_rotations(
    rotations: Sequence[object],
    parent_metric: np.ndarray,
    metric_tolerance: float,
) -> tuple[tuple[np.ndarray, ...], float, float]:
    if not rotations:
        raise ValueError("at least one fractional rotation is required")
    validated: list[np.ndarray] = []
    maximum_metric_error = 0.0
    maximum_determinant_error = 0.0
    metric_scale = max(float(np.linalg.norm(parent_metric)), 1.0)
    for index, value in enumerate(rotations):
        rotation = _matrix3(value, f"fractional rotation {index}")
        metric_error = float(
            np.linalg.norm(rotation.T @ parent_metric @ rotation - parent_metric)
            / metric_scale
        )
        determinant_error = abs(abs(float(np.linalg.det(rotation))) - 1.0)
        maximum_metric_error = max(maximum_metric_error, metric_error)
        maximum_determinant_error = max(
            maximum_determinant_error, determinant_error
        )
        if max(metric_error, determinant_error) > metric_tolerance:
            raise ValueError(
                f"fractional rotation {index} does not preserve the supplied parent metric"
            )
        if not any(
            np.allclose(rotation, old, rtol=0.0, atol=metric_tolerance)
            for old in validated
        ):
            validated.append(rotation)
    validated.sort(key=lambda matrix: tuple(np.round(matrix, 12).reshape(-1)))
    return tuple(validated), maximum_metric_error, maximum_determinant_error


def _set_identity_closure_errors(
    matrices: Sequence[np.ndarray],
) -> tuple[float, float]:
    identity = np.eye(matrices[0].shape[0])
    identity_error = min(
        float(np.max(np.abs(matrix - identity))) for matrix in matrices
    )
    closure_error = 0.0
    for left in matrices:
        for right in matrices:
            product = left @ right
            residual = min(
                float(np.max(np.abs(product - candidate)))
                for candidate in matrices
            )
            closure_error = max(closure_error, residual)
    return identity_error, closure_error


def _deterministic_projector_basis(
    projector: np.ndarray, dimension: int, tolerance: float
) -> np.ndarray:
    vectors: list[np.ndarray] = []
    for coordinate in range(projector.shape[0]):
        vector = np.array(projector[:, coordinate], dtype=float, copy=True)
        for _pass in range(2):
            for previous in vectors:
                vector -= previous * float(previous @ vector)
        norm = float(np.linalg.norm(vector))
        if norm <= tolerance:
            continue
        vector /= norm
        first = next(
            (value for value in vector if abs(float(value)) > tolerance), 0.0
        )
        if first < 0.0:
            vector *= -1.0
        vectors.append(vector)
        if len(vectors) == dimension:
            break
    if len(vectors) != dimension:
        raise RuntimeError("could not construct a deterministic fixed-space basis")
    return np.column_stack(vectors)


def compute_homogeneous_strain_modes(
    rotations_fractional: Sequence[object],
    parent_lattice: object,
    *,
    metric_tolerance: float = 1.0e-7,
    group_tolerance: float = 1.0e-7,
    rank_relative_tolerance: float = 1.0e-10,
    rank_absolute_tolerance: float = 1.0e-12,
) -> HomogeneousStrainModeResult:
    """Compute ``Fix_H`` from an actual fractional embedding and parent cell.

    ``rotations_fractional`` must contain the complete embedded point group in
    the spglib convention ``x' = R @ x``. Centering-related duplicates are
    removed before group validation. ``parent_lattice`` stores direct vectors
    as rows.
    """

    if metric_tolerance <= 0.0 or group_tolerance <= 0.0:
        raise ValueError("metric and group tolerances must be positive")
    if rank_relative_tolerance <= 0.0 or rank_absolute_tolerance <= 0.0:
        raise ValueError("rank tolerances must be positive")
    parent = _matrix3(parent_lattice, "parent lattice")
    if abs(float(np.linalg.det(parent))) <= np.finfo(float).eps:
        raise ValueError("parent lattice must be nonsingular")
    parent_metric = parent @ parent.T
    rotations, metric_error, determinant_error = _validated_fractional_rotations(
        rotations_fractional, parent_metric, metric_tolerance
    )
    identity_error, closure_error = _set_identity_closure_errors(rotations)
    if identity_error > group_tolerance:
        raise ValueError("fractional rotations do not contain the identity")
    if closure_error > group_tolerance:
        raise ValueError("fractional rotations are not closed under multiplication")

    sylvester = _metric_perturbation_map(parent_metric)
    congruence_actions = tuple(
        _metric_congruence_action(rotation) for rotation in rotations
    )
    strain_actions = tuple(
        np.linalg.solve(sylvester, action @ sylvester)
        for action in congruence_actions
    )
    representation_identity, representation_closure = (
        _set_identity_closure_errors(strain_actions)
    )
    if representation_identity > group_tolerance * 10.0:
        raise RuntimeError("parent-basis strain representation lacks an identity")
    if representation_closure > group_tolerance * 10.0:
        raise RuntimeError("parent-basis strain representation is not closed")

    # Work in y=sqrt(W)q so an ordinary Euclidean null-space basis is exactly
    # W-orthonormal in the official engineering coordinates.
    constraint_blocks = [
        (action @ sylvester - sylvester) @ _INV_SQRT_VOIGT_METRIC
        for action in congruence_actions
    ]
    constraints = np.vstack(constraint_blocks)
    _left, singular, right_t = np.linalg.svd(constraints, full_matrices=True)
    scale = float(singular[0]) if len(singular) else 0.0
    rank_tolerance = max(rank_absolute_tolerance, rank_relative_tolerance * scale)
    constraint_rank = int(np.sum(singular > rank_tolerance))
    fixed_dimension = 6 - constraint_rank
    null_basis = right_t[constraint_rank:].T
    projector_orthonormal = null_basis @ null_basis.T
    projector_orthonormal = (
        projector_orthonormal + projector_orthonormal.T
    ) / 2.0
    basis_orthonormal = _deterministic_projector_basis(
        projector_orthonormal,
        fixed_dimension,
        max(rank_absolute_tolerance, rank_relative_tolerance),
    )
    q_unit_matrix = _INV_SQRT_VOIGT_METRIC @ basis_orthonormal
    gram = q_unit_matrix.T @ ENGINEERING_VOIGT_METRIC @ q_unit_matrix
    orthonormality_error = float(
        np.max(np.abs(gram - np.eye(fixed_dimension)), initial=0.0)
    )

    modes: list[HomogeneousStrainMode] = []
    global_invariance_error = 0.0
    metric_scale = max(float(np.linalg.norm(parent_metric)), 1.0)
    for column in range(fixed_dimension):
        q_unit = q_unit_matrix[:, column]
        normfactor = float(np.max(np.abs(q_unit)))
        if normfactor <= rank_absolute_tolerance:
            raise RuntimeError("fixed-space basis contains a zero strain mode")
        q_raw = q_unit / normfactor
        delta_metric = parent_basis_metric_perturbation(q_unit, parent_metric)
        delta_scale = max(float(np.linalg.norm(delta_metric)), metric_scale)
        invariance_error = max(
            float(
                np.linalg.norm(
                    rotation.T @ delta_metric @ rotation - delta_metric
                )
                / delta_scale
            )
            for rotation in rotations
        )
        unit_norm_error = abs(
            float(q_unit.T @ ENGINEERING_VOIGT_METRIC @ q_unit) - 1.0
        )
        raw_error = abs(float(np.max(np.abs(q_raw))) - 1.0)
        global_invariance_error = max(global_invariance_error, invariance_error)
        modes.append(
            HomogeneousStrainMode(
                index=column + 1,
                q_unit=np.array(q_unit, copy=True),
                q_raw=np.array(q_raw, copy=True),
                normfactor=normfactor,
                irrep_label=None,
                label_status="unlabeled",
                invariance_max_relative_metric_error=invariance_error,
                unit_norm_error=unit_norm_error,
                raw_max_component_error=raw_error,
            )
        )

    fixed_projector_idempotence = float(
        np.max(
            np.abs(
                projector_orthonormal @ projector_orthonormal
                - projector_orthonormal
            )
        )
    )
    fixed_projector_constraint = float(
        np.max(np.abs(constraints @ projector_orthonormal), initial=0.0)
    )
    projector_engineering = (
        _INV_SQRT_VOIGT_METRIC
        @ projector_orthonormal
        @ _SQRT_VOIGT_METRIC
    )
    diagnostics = HomogeneousStrainModeDiagnostics(
        input_rotation_count=len(rotations_fractional),
        point_group_order=len(rotations),
        fixed_dimension=fixed_dimension,
        rotation_metric_preservation_max_relative_error=metric_error,
        rotation_determinant_max_error=determinant_error,
        group_identity_error=identity_error,
        group_closure_max_error=closure_error,
        representation_identity_error=representation_identity,
        representation_closure_max_error=representation_closure,
        fixed_projector_idempotence_max_error=fixed_projector_idempotence,
        fixed_projector_constraint_max_error=fixed_projector_constraint,
        invariance_max_relative_metric_error=global_invariance_error,
        orthonormality_max_abs_error=orthonormality_error,
    )
    return HomogeneousStrainModeResult(
        modes=tuple(modes),
        diagnostics=diagnostics,
        projector_engineering=projector_engineering,
        projector_orthonormal=projector_orthonormal,
        parent_metric=np.array(parent_metric, copy=True),
    )


def apply_canonical_strain_basis(
    fixed_space: HomogeneousStrainModeResult,
    definitions: Sequence[CanonicalStrainModeDefinition],
    *,
    validation_tolerance: float = 5.0e-4,
) -> HomogeneousStrainModeResult:
    """Replace a fallback basis with ISO's canonical parent-Gamma basis.

    ``definitions`` must already have been filtered by the exact embedding's
    invariant directions. This function does not trust labels as a numerical
    proof: it independently requires the W-orthonormalized canonical span to
    equal the metric-derived fixed space. The relaxed default reflects ISO's
    printed decimal coefficients (for example ``1.732``), rather than treating
    those display values as exact radicals.
    """

    expected = fixed_space.diagnostics.fixed_dimension
    if len(definitions) != expected:
        raise ValueError(
            "canonical strain mode count does not equal metric fixed-space dimension"
        )
    if validation_tolerance <= 0.0:
        raise ValueError("canonical strain validation tolerance must be positive")

    raw_columns: list[np.ndarray] = []
    unit_columns: list[np.ndarray] = []
    normfactors: list[float] = []
    for definition in definitions:
        raw = np.asarray(definition.q_raw, dtype=float)
        if raw.shape != (6,) or not np.all(np.isfinite(raw)):
            raise ValueError("canonical strain q_raw must contain six finite values")
        raw_maximum = float(np.max(np.abs(raw), initial=0.0))
        if abs(raw_maximum - 1.0) > 1.0e-10:
            raise ValueError(
                "canonical strain q_raw must use max-component-one website coordinates"
            )
        squared_norm = float(raw @ ENGINEERING_VOIGT_METRIC @ raw)
        if squared_norm <= np.finfo(float).eps:
            raise ValueError("canonical strain q_raw must be nonzero")
        normfactor = 1.0 / np.sqrt(squared_norm)
        raw_columns.append(np.array(raw, copy=True))
        unit_columns.append(normfactor * raw)
        normfactors.append(float(normfactor))

    unit_matrix = np.column_stack(unit_columns)
    y_matrix = _SQRT_VOIGT_METRIC @ unit_matrix
    singular = np.linalg.svd(y_matrix, compute_uv=False)
    if len(singular) != expected or float(np.min(singular, initial=1.0)) <= 1.0e-10:
        raise ValueError("canonical strain modes are linearly dependent")
    orthonormal_y, _triangular = np.linalg.qr(y_matrix)
    canonical_projector = orthonormal_y @ orthonormal_y.T
    projector_error = float(
        np.max(
            np.abs(canonical_projector - fixed_space.projector_orthonormal),
            initial=0.0,
        )
    )
    gram_error = float(
        np.max(
            np.abs(
                unit_matrix.T @ ENGINEERING_VOIGT_METRIC @ unit_matrix
                - np.eye(expected)
            ),
            initial=0.0,
        )
    )
    validation_error = max(projector_error, gram_error)
    if validation_error > validation_tolerance:
        raise ValueError(
            "canonical ISO strain basis does not match the metric fixed space"
        )

    modes: list[HomogeneousStrainMode] = []
    for index, (definition, raw, unit, normfactor) in enumerate(
        zip(
            definitions,
            raw_columns,
            unit_columns,
            normfactors,
            strict=True,
        ),
        start=1,
    ):
        y = _SQRT_VOIGT_METRIC @ unit
        fixed_residual = float(
            np.linalg.norm((np.eye(6) - fixed_space.projector_orthonormal) @ y)
        )
        modes.append(
            HomogeneousStrainMode(
                index=index,
                q_unit=np.array(unit, copy=True),
                q_raw=np.array(raw, copy=True),
                normfactor=normfactor,
                irrep_label=definition.irrep_label,
                label_status="canonical_iso_macro",
                invariance_max_relative_metric_error=fixed_residual,
                unit_norm_error=abs(
                    float(unit @ ENGINEERING_VOIGT_METRIC @ unit) - 1.0
                ),
                raw_max_component_error=abs(
                    float(np.max(np.abs(raw))) - 1.0
                ),
                canonical_label=definition.label,
                irrep_direction=definition.irrep_direction,
            )
        )
    return HomogeneousStrainModeResult(
        modes=tuple(modes),
        diagnostics=fixed_space.diagnostics,
        projector_engineering=fixed_space.projector_engineering,
        projector_orthonormal=fixed_space.projector_orthonormal,
        parent_metric=fixed_space.parent_metric,
        basis_source="iso_display_distortion_rank_12",
        export_ready=True,
        fixed_space_validation_max_error=validation_error,
        irrep_labels_available=True,
        irrep_label_reason="canonical labels parsed from ISO rank-[12] macroscopic blocks",
    )


__all__ = [
    "CHILD_LATTICE_ACTION",
    "COLUMN_LATTICE_ACTION",
    "ENGINEERING_VOIGT_METRIC",
    "ENGINEERING_VOIGT_ORDER",
    "METRIC_EQUATION",
    "PARENT_BASIS_COORDINATE_FRAME",
    "ROW_LATTICE_ACTION",
    "CanonicalStrainModeDefinition",
    "HomogeneousStrainMode",
    "HomogeneousStrainModeDiagnostics",
    "HomogeneousStrainModeResult",
    "apply_canonical_strain_basis",
    "canonical_strain_definitions_from_iso",
    "compute_homogeneous_strain_modes",
    "engineering_voigt_to_tensor",
    "parent_basis_metric_perturbation",
    "parent_basis_strain_action",
    "tensor_to_engineering_voigt",
]
