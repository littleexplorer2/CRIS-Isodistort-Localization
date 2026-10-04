from __future__ import annotations

import hashlib
import math
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from isocore.backend.iso_mode_models import (
    MicroscopicColumnProvenance,
    ModeIdentity,
    microscopic_atom_order_id,
)
from isocore.io.cif_displacive_model import (
    ChildAtom,
    DisplaciveModelError,
    ExactChildFrame,
    ExactSeitzOperation,
    ModeScaleProvenance,
    RawModeColumn,
    build_child_orbits,
    build_displacive_cif_model,
    conjugate_seitz_operations,
    derive_centering_index,
    primitive_rss_normfactors,
    project_modes_to_free_rows,
    site_fixed_coordinate_model,
    transform_fractional_displacement,
    transform_fractional_position,
    transformed_reference_lattice,
)

I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ZERO = (0, 0, 0)
ANTI_DIAGONAL_MIRROR = ((0, -1, 0), (-1, 0, 0), (0, 0, 1))
# This skew metric has Gxx=Gyy and is invariant under ANTI_DIAGONAL_MIRROR.
SKEW_LATTICE = ((2.0, 0.0, 0.0), (0.5, math.sqrt(15.0) / 2.0, 0.0), (0.0, 0.0, 4.0))
HEXAGONAL_SIXFOLD = ((1, -1, 0), (1, 0, 0), (0, 0, 1))
HEXAGONAL_LATTICE = ((1.0, 0.0, 0.0), (-0.5, math.sqrt(3.0) / 2.0, 0.0), (0.0, 0.0, 2.0))
F_PRIMITIVE = (
    (0, Fraction(1, 2), Fraction(1, 2)),
    (Fraction(1, 2), 0, Fraction(1, 2)),
    (Fraction(1, 2), Fraction(1, 2), 0),
)


def _op(rotation=I3, translation=ZERO) -> ExactSeitzOperation:
    return ExactSeitzOperation.from_values(rotation, translation)


def _canonical_source() -> ModeScaleProvenance:
    return ModeScaleProvenance(
        source="iso_microscopic_display_distortion",
        convention="authoritative_fractional_direction_arbitrary_positive_scale",
        direction_resolved=True,
        sign_resolved=True,
    )


def _mode(
    display_label: str,
    column: int,
    frame_id: str,
    atoms: tuple[ChildAtom, ...],
    vectors,
    *,
    provenance: ModeScaleProvenance | None = None,
) -> RawModeColumn:
    atom_ids = tuple(atom.atom_id for atom in atoms)
    global_irrep = "TEST1+"
    site_irrep = "A"
    identity = ModeIdentity(
        parent_sg=1,
        global_irrep=global_irrep,
        k_coordinates=("0", "0", "0"),
        wyckoff_letter="a",
        orbit_id=atom_ids[0] if atom_ids else "",
        site_irrep=site_irrep,
        component_index=column,
        component_label=display_label,
        source="iso_microscopic",
        status="verified",
    )
    row_digest = hashlib.sha256(
        repr((frame_id, atom_ids)).encode("utf-8")
    ).hexdigest()
    vector_digest = hashlib.sha256(
        repr((display_label, vectors)).encode("utf-8")
    ).hexdigest()
    microscopic_provenance = MicroscopicColumnProvenance(
        query_digest=hashlib.sha256(b"synthetic exact query").hexdigest(),
        query_order=0,
        query_irrep_label=global_irrep,
        source_parent_sg=1,
        source_k_coordinates=("0", "0", "0"),
        source_subgroup_space_group_number=1,
        source_subgroup_basis=tuple(
            tuple(str(value) for value in row) for row in I3
        ),
        source_subgroup_origin=tuple(str(value) for value in ZERO),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        source_block_order=0,
        source_global_irrep=global_irrep,
        source_wyckoff_letter="a",
        source_site_irrep=site_irrep,
        source_direction_symbol="VECTOR,A",
        source_column_index=column,
        source_block_column_count=64,
        source_row_count=len(atom_ids),
        exact_row_digest=row_digest,
        exact_vector_digest=vector_digest,
        source_float_column_digest=vector_digest,
        source_frame_id=frame_id,
        source_atom_order_id=microscopic_atom_order_id(atom_ids),
        unresolved_fields=(),
        unresolved_diagnostics=(),
    )
    return RawModeColumn.from_values(
        identity.stable_token,
        column,
        frame_id,
        atom_ids,
        vectors,
        provenance or _canonical_source(),
        mode_identity=identity,
        microscopic_provenance=microscopic_provenance,
        label=display_label,
    )


def _anti_coupled_fixture(frame_id: str = "anti"):
    operations = (_op(), _op(ANTI_DIAGONAL_MIRROR))
    atoms = (ChildAtom.from_values("A", ZERO, "X"),)
    modes = (
        _mode("anti_xy", 0, frame_id, atoms, ((1.0, -1.0, 0.0),)),
        _mode("z", 1, frame_id, atoms, ((0.0, 0.0, 1.0),)),
    )
    frame = ExactChildFrame.from_values(frame_id, operations, SKEW_LATTICE)
    return frame, atoms, modes


def _f_centered_fixture(frame_id: str = "F"):
    half = Fraction(1, 2)
    translations = (
        (0, 0, 0),
        (0, half, half),
        (half, 0, half),
        (half, half, 0),
    )
    operations = tuple(_op(I3, translation) for translation in translations)
    atoms = tuple(
        ChildAtom.from_values(f"A{index}", translation, "X") for index, translation in enumerate(translations)
    )
    vectors = np.eye(3)
    modes = tuple(
        _mode(
            f"mode_{axis}",
            axis,
            frame_id,
            atoms,
            tuple(tuple(vectors[axis]) for _atom in atoms),
        )
        for axis in range(3)
    )
    return operations, atoms, modes


def _rotation_group(generator, order: int) -> tuple[ExactSeitzOperation, ...]:
    current = np.eye(3, dtype=int)
    matrix = np.asarray(generator, dtype=int)
    operations = []
    for _index in range(order):
        operations.append(_op(current.tolist()))
        current = matrix @ current
    np.testing.assert_array_equal(current, np.eye(3, dtype=int))
    return tuple(operations)


def test_fraction_stabilizer_rref_emits_anti_coupled_form_and_delta_rows() -> None:
    operations = (_op(), _op(ANTI_DIAGONAL_MIRROR))

    site = site_fixed_coordinate_model("A", ZERO, operations)

    assert site.symmform == "Dx,-Dx,Dz"
    assert site.component_expressions == ("Dx", "-Dx", "Dz")
    assert site.pivot_axes == (0, 2)
    assert tuple(parameter.suffix for parameter in site.parameters) == ("dx", "dz")
    assert site.fixed_basis == (
        (Fraction(1), Fraction(0)),
        (Fraction(-1), Fraction(0)),
        (Fraction(0), Fraction(1)),
    )

    frame, atoms, modes = _anti_coupled_fixture()
    model = build_displacive_cif_model(frame, atoms, ("A",), modes)
    assert tuple(row.delta_label for row in model.projection.free_rows) == (
        "A_dx",
        "A_dz",
    )
    np.testing.assert_allclose(model.projection.matrix_array(), np.eye(2), atol=0.0)


def test_nonzero_seitz_translation_selects_stabilizer_but_not_tangent_action() -> None:
    translated_mirror = _op(
        ((-1, 0, 0), (0, 1, 0), (0, 0, 1)),
        (Fraction(1, 2), 0, 0),
    )
    operations = (_op(), translated_mirror)

    fixed = site_fixed_coordinate_model("special", (Fraction(1, 4), Fraction(1, 7), Fraction(2, 9)), operations)
    general = site_fixed_coordinate_model("general", ZERO, operations)

    assert fixed.stabilizer_indices == (0, 1)
    assert fixed.symmform == "0,Dy,Dz"
    assert len(general.stabilizer_indices) == 1
    assert general.symmform == "Dx,Dy,Dz"


def test_equivalent_transporters_are_accepted_and_conflicting_arrows_fail_closed() -> None:
    frame, atoms, modes = _anti_coupled_fixture()
    orbits = build_child_orbits(atoms, ("A",), frame.operations)

    transporter = orbits[0].transporters[0]
    assert len(transporter.candidate_operation_indices) == 2
    assert transporter.induced_fixed_map == (
        (Fraction(1), Fraction(0)),
        (Fraction(-1), Fraction(0)),
        (Fraction(0), Fraction(1)),
    )
    valid = project_modes_to_free_rows(
        orbits,
        atoms,
        frame.operations,
        modes,
        frame_id=frame.frame_id,
        reference_lattice=frame.reference_lattice,
    )
    assert valid.diagnostics.max_transporter_residual_angstrom == 0.0

    invalid = (_mode("forbidden_x", 0, frame.frame_id, atoms, ((1.0, 0.0, 0.0),)),)
    with pytest.raises(DisplaciveModelError) as caught:
        project_modes_to_free_rows(
            orbits,
            atoms,
            frame.operations,
            invalid,
            frame_id=frame.frame_id,
            reference_lattice=frame.reference_lattice,
            require_complete=False,
        )
    assert caught.value.code == "conflicting_mode_transporters"


def test_non_diagonal_gl3z_and_origin_conjugacy_preserve_cartesian_modes_and_norms() -> None:
    old_frame, old_atoms, old_modes = _anti_coupled_fixture("old")
    transform = ((1, 1, 0), (0, 1, 0), (0, 0, 1))
    origin = (Fraction(1, 3), Fraction(1, 4), Fraction(1, 5))
    new_operations = conjugate_seitz_operations(old_frame.operations, transform, origin)
    new_atoms = (ChildAtom.from_values("A", transform_fractional_position(ZERO, transform, origin), "X"),)
    new_modes = tuple(
        _mode(
            mode.display_label,
            mode.column_index,
            "new",
            new_atoms,
            (
                tuple(
                    float(value)
                    for value in transform_fractional_displacement(
                        old_modes[mode.column_index].displacements[0], transform
                    )
                ),
            ),
        )
        for mode in old_modes
    )
    new_lattice = transformed_reference_lattice(SKEW_LATTICE, transform)
    new_frame = ExactChildFrame.from_values("new", new_operations, new_lattice)

    old_model = build_displacive_cif_model(old_frame, old_atoms, ("A",), old_modes)
    new_model = build_displacive_cif_model(new_frame, new_atoms, ("A",), new_modes)

    assert new_model.orbits[0].site.symmform == "0,Dy,Dz"
    np.testing.assert_allclose(new_model.norms.normfactors, old_model.norms.normfactors, rtol=1.0e-14)
    for old_mode, new_mode in zip(old_modes, new_modes, strict=True):
        old_cart = np.asarray(old_mode.displacements) @ np.asarray(SKEW_LATTICE)
        new_cart = np.asarray(new_mode.displacements) @ np.asarray(new_lattice)
        np.testing.assert_allclose(new_cart, old_cart, rtol=1.0e-14, atol=1.0e-14)


def test_reference_metric_uses_column_coordinate_rt_g_r_convention() -> None:
    operations = _rotation_group(HEXAGONAL_SIXFOLD, 6)
    metric = np.asarray(HEXAGONAL_LATTICE) @ np.asarray(HEXAGONAL_LATTICE).T
    rotation = np.asarray(HEXAGONAL_SIXFOLD, dtype=float)

    np.testing.assert_allclose(rotation.T @ metric @ rotation, metric, atol=2.0e-15)
    assert np.max(np.abs(rotation @ metric @ rotation.T - metric)) > 1.0
    frame = ExactChildFrame.from_values("hex", operations, HEXAGONAL_LATTICE)

    transform = ((1, 1, 0), (0, 1, 0), (0, 0, 1))
    transformed = ExactChildFrame.from_values(
        "hex_sheared",
        conjugate_seitz_operations(operations, transform, (Fraction(1, 7), 0, 0)),
        transformed_reference_lattice(HEXAGONAL_LATTICE, transform),
    )
    assert frame.frame_id == "hex"
    assert transformed.frame_id == "hex_sheared"


def test_reference_metric_rejects_rotation_incompatible_lattice() -> None:
    operations = _rotation_group(((0, -1, 0), (1, 0, 0), (0, 0, 1)), 4)

    with pytest.raises(DisplaciveModelError) as caught:
        ExactChildFrame.from_values(
            "invalid_c4",
            operations,
            ((2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 4.0)),
        )
    assert caught.value.code == "reference_metric_symmetry_mismatch"


def test_f_centering_and_skew_reference_metric_give_primitive_rss_norms() -> None:
    operations, atoms, modes = _f_centered_fixture()
    frame = ExactChildFrame.from_values(
        "F",
        operations,
        SKEW_LATTICE,
        conventional_basis=I3,
        primitive_translation_basis=F_PRIMITIVE,
    )

    model = build_displacive_cif_model(frame, atoms, ("A0",), modes)

    assert model.centering.index == 4
    assert model.centering.determinant_index == 4
    assert model.diagnostics.status == "validated"
    np.testing.assert_allclose(model.projection.matrix_array(), np.eye(3), atol=0.0)
    expected = tuple(1.0 / np.linalg.norm(np.asarray(SKEW_LATTICE)[axis]) for axis in range(3))
    np.testing.assert_allclose(model.norms.normfactors, expected, rtol=1.0e-14)
    np.testing.assert_allclose(
        model.norms.diagnostics.conventional_rss,
        tuple(4.0 * np.dot(row, row) for row in np.asarray(SKEW_LATTICE)),
        rtol=1.0e-14,
    )


def test_centered_frame_without_lattice_witness_uses_only_exact_translation_index() -> None:
    operations, atoms, modes = _f_centered_fixture("F_default")
    frame = ExactChildFrame.from_values("F_default", operations, SKEW_LATTICE)

    model = build_displacive_cif_model(frame, atoms, ("A0",), modes)

    assert model.centering.index == 4
    assert model.centering.determinant_index is None


def test_centering_translation_and_lattice_indices_must_agree() -> None:
    half = Fraction(1, 2)
    operations = tuple(
        _op(I3, translation)
        for translation in (
            ZERO,
            (0, half, half),
            (half, 0, half),
            (half, half, 0),
        )
    )

    with pytest.raises(DisplaciveModelError) as caught:
        derive_centering_index(
            operations,
            conventional_basis=I3,
            primitive_translation_basis=I3,
        )
    assert caught.value.code == "centering_index_mismatch"


def test_norm_uses_reference_lattice_instead_of_final_strained_lattice() -> None:
    atoms = (ChildAtom.from_values("A", ZERO, "X"),)
    mode = _mode("x", 0, "P", atoms, ((1.0, 0.0, 0.0),))
    strained = np.asarray(SKEW_LATTICE, dtype=float).copy()
    strained[0] *= 1.25

    reference_norm = primitive_rss_normfactors((mode,), SKEW_LATTICE, 1)
    wrong_strained_norm = primitive_rss_normfactors((mode,), strained, 1)

    assert reference_norm.normfactors[0] == pytest.approx(0.5)
    assert wrong_strained_norm.normfactors[0] == pytest.approx(0.4)
    assert reference_norm.normfactors != wrong_strained_norm.normfactors


def test_unsorted_mode_input_reorders_columns_norms_and_amplitudes_together() -> None:
    _frame, atoms, modes = _anti_coupled_fixture()
    frame = ExactChildFrame.from_values(
        "anti",
        _frame.operations,
        np.asarray(SKEW_LATTICE),
        conventional_basis=np.eye(3, dtype=int),
        primitive_translation_basis=np.eye(3, dtype=int),
    )

    model = build_displacive_cif_model(
        frame,
        atoms,
        ("A",),
        (modes[1], modes[0]),
        amplitudes=(7.0, 3.0),
    )

    assert model.projection.mode_ids == tuple(mode.mode_id for mode in modes)
    assert model.norms.mode_ids == model.projection.mode_ids
    assert model.amplitudes == (3.0, 7.0)
    np.testing.assert_allclose(model.projection.matrix_array(), np.eye(2), atol=0.0)


def test_canonical_modes_reorder_displacement_rows_to_model_atom_order() -> None:
    operations = (_op(),)
    atoms = (
        ChildAtom.from_values("A", ZERO, "Fe"),
        ChildAtom.from_values("B", (Fraction(1, 4),) * 3, "Al"),
    )
    reversed_atoms = tuple(reversed(atoms))
    frame = ExactChildFrame.from_values("P", operations, SKEW_LATTICE)
    modes: list[RawModeColumn] = []
    expected: list[np.ndarray] = []
    column = 0
    for atom_index in range(2):
        for axis in range(3):
            model_order = np.zeros((2, 3), dtype=float)
            model_order[atom_index, axis] = 1.0
            expected.append(model_order)
            modes.append(
                _mode(
                    f"atom-{atom_index}-axis-{axis}",
                    column,
                    frame.frame_id,
                    reversed_atoms,
                    tuple(tuple(row) for row in model_order[::-1]),
                )
            )
            column += 1

    model = build_displacive_cif_model(
        frame,
        atoms,
        ("A", "B"),
        modes,
    )

    assert model.diagnostics.status == "validated"
    for source, canonical, expected_rows in zip(
        modes, model.canonical_modes, expected, strict=True
    ):
        assert source.atom_ids == ("B", "A")
        assert canonical.atom_ids == ("A", "B")
        np.testing.assert_allclose(canonical.displacements, expected_rows, atol=0.0)
        assert canonical.mode_identity is source.mode_identity
        assert canonical.microscopic_provenance is not source.microscopic_provenance
        assert (
            canonical.microscopic_provenance.exact_row_digest
            == source.microscopic_provenance.exact_row_digest
        )
        assert (
            canonical.microscopic_provenance.exact_vector_digest
            == source.microscopic_provenance.exact_vector_digest
        )
        assert (
            source.microscopic_provenance.source_atom_order_id
            == microscopic_atom_order_id(("B", "A"))
        )
        assert (
            canonical.microscopic_provenance.source_atom_order_id
            == microscopic_atom_order_id(("A", "B"))
        )


def test_free_coordinate_columns_are_canonicalized_with_full_modes_and_norms() -> None:
    operations = (_op(),)
    atoms = (ChildAtom.from_values("A", ZERO, "X"),)
    frame = ExactChildFrame.from_values("P", operations, SKEW_LATTICE)
    modes = (
        _mode("x", 0, "P", atoms, ((7.0, 0.0, 0.0),)),
        _mode("negative_y", 1, "P", atoms, ((0.0, -3.0, 0.0),)),
        _mode("z", 2, "P", atoms, ((0.0, 0.0, 2.0),)),
    )

    model = build_displacive_cif_model(frame, atoms, ("A",), modes, amplitudes=(1.0, -2.0, 3.0))

    np.testing.assert_allclose(
        model.projection.matrix_array(),
        np.diag((1.0, -1.0, 1.0)),
        atol=0.0,
    )
    np.testing.assert_allclose(
        model.projection.input_to_canonical_scale_factors,
        (1.0 / 7.0, 1.0 / 3.0, 1.0 / 2.0),
        rtol=1.0e-15,
    )
    np.testing.assert_allclose(
        np.asarray([mode.displacements[0] for mode in model.canonical_modes]),
        np.diag((1.0, -1.0, 1.0)),
        atol=0.0,
    )
    expected_norms = tuple(1.0 / np.linalg.norm(np.asarray(SKEW_LATTICE)[axis]) for axis in range(3))
    np.testing.assert_allclose(model.norms.normfactors, expected_norms, rtol=1.0e-14)
    assert model.amplitudes == (1.0, -2.0, 3.0)
    assert model.amplitude_convention == "As_unit_mode_angstrom"

    for source, canonical, normfactor in zip(modes, model.canonical_modes, model.norms.normfactors, strict=True):
        source_cartesian = np.asarray(source.displacements) @ np.asarray(SKEW_LATTICE)
        source_unit = source_cartesian / np.linalg.norm(source_cartesian)
        canonical_unit = np.asarray(canonical.displacements) @ np.asarray(SKEW_LATTICE) * normfactor
        np.testing.assert_allclose(canonical_unit, source_unit, rtol=1.0e-14, atol=1.0e-14)


def test_raw_coefficient_amplitudes_and_zero_or_unresolved_columns_fail_closed() -> None:
    operations = (_op(),)
    atoms = (ChildAtom.from_values("A", ZERO, "X"),)
    frame = ExactChildFrame.from_values("P", operations, SKEW_LATTICE)
    basis = (
        _mode("x", 0, "P", atoms, ((1.0, 0.0, 0.0),)),
        _mode("y", 1, "P", atoms, ((0.0, 1.0, 0.0),)),
        _mode("z", 2, "P", atoms, ((0.0, 0.0, 1.0),)),
    )

    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(
            frame,
            atoms,
            ("A",),
            basis,
            amplitudes=(1.0, 2.0, 3.0),
            amplitude_convention="raw_coefficient",
        )
    assert caught.value.code == "unsupported_mode_amplitude_convention"

    zero = (_mode("zero", 0, "P", atoms, ((0.0, 0.0, 0.0),)),)
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), zero, require_complete=False)
    assert caught.value.code == "zero_free_coordinate_mode"

    mixed_source = ModeScaleProvenance(
        source="gram_schmidt",
        convention="arbitrary_span_basis",
        direction_resolved=False,
        sign_resolved=True,
    )
    mixed = (_mode("mixed", 0, "P", atoms, ((1.0, 0.0, 0.0),), provenance=mixed_source),)
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), mixed, require_complete=False)
    assert caught.value.code == "unresolved_mode_provenance"


def test_relaxed_diagnostic_build_never_launders_unresolved_source_or_completeness() -> None:
    operations = (_op(),)
    atoms = (ChildAtom.from_values("A", ZERO, "X"),)
    frame = ExactChildFrame.from_values("P", operations, SKEW_LATTICE)
    unresolved_scale = ModeScaleProvenance(
        source="gram_schmidt",
        convention="arbitrary_span_basis",
        direction_resolved=False,
        sign_resolved=False,
    )
    full_modes = tuple(
        replace(
            _mode(
                f"axis-{axis}",
                axis,
                frame.frame_id,
                atoms,
                (tuple(np.eye(3)[axis]),),
            ),
            scale_provenance=unresolved_scale,
            mode_identity=None,
            microscopic_provenance=None,
        )
        for axis in range(3)
    )

    diagnostic = build_displacive_cif_model(
        frame,
        atoms,
        ("A",),
        full_modes,
        require_resolved_provenance=False,
    )
    assert diagnostic.diagnostics.status == "diagnostic"
    assert diagnostic.diagnostics.validation_issues
    assert all(
        not mode.scale_provenance.direction_resolved
        and not mode.scale_provenance.sign_resolved
        for mode in diagnostic.canonical_modes
    )

    incomplete = build_displacive_cif_model(
        frame,
        atoms,
        ("A",),
        (_mode("x-only", 0, frame.frame_id, atoms, ((1.0, 0.0, 0.0),)),),
        require_complete=False,
    )
    assert incomplete.diagnostics.status == "diagnostic"
    assert incomplete.diagnostics.mode_count == 1
    assert incomplete.diagnostics.free_coordinate_count == 3
    assert any("incomplete" in issue for issue in incomplete.diagnostics.validation_issues)


def test_empty_display_label_cannot_enter_validated_model() -> None:
    frame, atoms, modes = _anti_coupled_fixture()
    empty_label = (replace(modes[0], label=""), modes[1])

    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), empty_label)
    assert caught.value.code == "unresolved_mode_provenance"


@pytest.mark.parametrize(
    ("identity_changes", "message"),
    (
        ({"orbit_id": ""}, "physical orbit ID"),
        ({"k_coordinates": ("0", "g", "0")}, "exact rational k coordinates"),
    ),
)
def test_verified_mode_identity_requires_orbit_and_exact_rational_k(
    identity_changes: dict[str, object],
    message: str,
) -> None:
    frame, atoms, modes = _anti_coupled_fixture()
    identity = replace(modes[0].mode_identity, **identity_changes)
    invalid = replace(modes[0], mode_id=identity.stable_token, mode_identity=identity)

    assert message in invalid.source_binding_issue()
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), (invalid, modes[1]))
    assert caught.value.code == "unresolved_mode_provenance"


def test_raw_mode_rejects_forged_atom_order_binding_for_validated_export() -> None:
    frame, atoms, modes = _anti_coupled_fixture()
    forged = replace(
        modes[0],
        microscopic_provenance=replace(
            modes[0].microscopic_provenance,
            source_atom_order_id=microscopic_atom_order_id(("forged",)),
        ),
    )

    assert "atom order" in forged.source_binding_issue()
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), (forged, modes[1]))
    assert caught.value.code == "unresolved_mode_provenance"


def test_raw_mode_legacy_positional_fields_remain_diagnostic_compatible() -> None:
    mode = RawModeColumn(
        "legacy-mode",
        0,
        "legacy-frame",
        ("A",),
        ((1.0, 0.0, 0.0),),
        _canonical_source(),
        "legacy-label",
    )

    assert mode.label == "legacy-label"
    assert mode.mode_identity is None
    assert mode.microscopic_provenance is None
    assert not mode.has_resolved_source


def test_projection_does_not_use_floating_matrix_inverse(monkeypatch: pytest.MonkeyPatch) -> None:
    frame, atoms, modes = _anti_coupled_fixture()

    def forbidden_inverse(*_args, **_kwargs):
        raise AssertionError("floating matrix inverse must not be used for exact Seitz transport")

    monkeypatch.setattr(np.linalg, "inv", forbidden_inverse)
    model = build_displacive_cif_model(frame, atoms, ("A",), modes)
    np.testing.assert_allclose(model.projection.matrix_array(), np.eye(2), atol=0.0)


@pytest.mark.parametrize(
    ("factory", "code"),
    (
        (lambda: ChildAtom.from_values("   ", ZERO, "X"), "invalid_atom_id"),
        (
            lambda: RawModeColumn.from_values(
                "mode", 0.5, "frame", ("A",), ((1.0, 0.0, 0.0),), _canonical_source()
            ),
            "invalid_mode_identity",
        ),
        (
            lambda: RawModeColumn.from_values(
                "mode", 0, "   ", ("A",), ((1.0, 0.0, 0.0),), _canonical_source()
            ),
            "invalid_mode_frame",
        ),
        (
            lambda: RawModeColumn.from_values(
                "mode", 0, "frame", ("   ",), ((1.0, 0.0, 0.0),), _canonical_source()
            ),
            "invalid_mode_atom_order",
        ),
        (
            lambda: primitive_rss_normfactors(
                (
                    RawModeColumn.from_values(
                        "mode", 0, "frame", ("A",), ((1.0, 0.0, 0.0),), _canonical_source()
                    ),
                ),
                SKEW_LATTICE,
                1.5,
            ),
            "invalid_centering_index",
        ),
    ),
)
def test_invalid_identifiers_and_noninteger_indices_fail_closed(factory, code: str) -> None:
    with pytest.raises(DisplaciveModelError) as caught:
        factory()
    assert caught.value.code == code


def test_unresolved_or_rank_deficient_complete_modes_fail_closed() -> None:
    frame, atoms, modes = _anti_coupled_fixture()
    unresolved_source = ModeScaleProvenance(
        source="gram_schmidt",
        convention="arbitrary_span_basis",
        direction_resolved=False,
        sign_resolved=False,
    )
    unresolved = (
        _mode(
            "unresolved",
            0,
            frame.frame_id,
            atoms,
            ((1.0, -1.0, 0.0),),
            provenance=unresolved_source,
        ),
        modes[1],
    )
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), unresolved)
    assert caught.value.code == "unresolved_mode_provenance"

    duplicate_columns = (
        modes[0],
        _mode(
            "duplicate",
            1,
            frame.frame_id,
            atoms,
            ((1.0, -1.0, 0.0),),
        ),
    )
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), duplicate_columns)
    assert caught.value.code == "rank_deficient_mode_matrix"
