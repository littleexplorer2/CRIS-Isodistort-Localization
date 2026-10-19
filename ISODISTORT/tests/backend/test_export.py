"""Tests for backend.export; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.export import StructureExporter
from backend.export import exact_child_operations
from backend.export import isodistort_cif as cif_writer
from backend.export.cif_displacive_model import ChildAtom
from backend.export.cif_displacive_model import DisplaciveModelError
from backend.export.cif_displacive_model import ExactChildFrame
from backend.export.cif_displacive_model import ExactSeitzOperation
from backend.export.cif_displacive_model import ModeScaleProvenance
from backend.export.cif_displacive_model import RawModeColumn
from backend.export.cif_displacive_model import build_child_orbits
from backend.export.cif_displacive_model import build_displacive_cif_model
from backend.export.cif_displacive_model import conjugate_seitz_operations
from backend.export.cif_displacive_model import derive_centering_index
from backend.export.cif_displacive_model import primitive_rss_normfactors
from backend.export.cif_displacive_model import project_modes_to_free_rows
from backend.export.cif_displacive_model import site_fixed_coordinate_model
from backend.export.cif_displacive_model import transform_fractional_displacement
from backend.export.cif_displacive_model import transform_fractional_position
from backend.export.cif_displacive_model import transformed_reference_lattice
from backend.export.displacive_export import DisplaciveExportData
from backend.export.displacive_export import DisplaciveSubgroupIdentity
from backend.export.displacive_export import ExactParentChildEmbedding
from backend.export.displacive_export import ParentChildSiteMapping
from backend.export.displacive_export import ParentOrbitType
from backend.export.displacive_export import _periodic_cartesian_residual
from backend.export.displacive_export import _validate_emitted_space_group
from backend.export.distortion_formats import FORMAT_CIF
from backend.export.distortion_formats import FORMAT_ISOVIZ
from backend.export.distortion_formats import FORMAT_MODES
from backend.export.distortion_formats import FORMAT_TOPAS
from backend.export.distortion_formats import StrainExportData
from backend.export.distortion_formats import SubgroupExportSpec
from backend.export.distortion_formats import build_export_zip
from backend.export.distortion_formats import cart_normalized_mode_matrix
from backend.export.distortion_formats import folder_label_for_subgroup
from backend.export.distortion_formats import format_filename
from backend.export.distortion_formats import method3_case_folder
from backend.export.distortion_formats import parse_export_formats
from backend.export.distortion_formats import parse_export_method
from backend.export.distortion_formats import render_cif
from backend.export.distortion_formats import render_complete_modes
from backend.export.distortion_formats import render_isoviz
from backend.export.distortion_formats import render_topas
from backend.export.distortion_formats import resolve_strained_export_structure
from backend.export.distortion_formats import subgroup_identity_digest
from backend.export.distortion_formats import subgroup_label
from backend.export.distortion_formats import unique_folder_name
from backend.export.distortion_formats import write_subgroup_files
from backend.export.isodistort_cif import _apply_origin_choice
from backend.export.isodistort_cif import _coordinate_formula
from backend.export.isodistort_cif import _irrep_kovalev_tag
from backend.export.isodistort_cif import _k_point_comment
from backend.export.isodistort_cif import _setting
from backend.export.isodistort_cif import _space_hm_symbol
from backend.export.isodistort_cif import _validate_displacive_writer_frame
from backend.export.isodistort_cif import _validated_subgroup_sites
from backend.export.isodistort_cif import _wyckoff_letter
from backend.export.isodistort_isoviz import _AtomLayout
from backend.export.isodistort_isoviz import _authoritative_mode_parentatom
from backend.export.isodistort_isoviz import render_isodistort_isoviz
from backend.input_cif import read_cif
from backend.input_cif import read_structure
from backend.method2.distortion_mapper import DistortionMapper as _microscopic_column_provenance_DistortionMapper
from backend.method4.strain_modes import CanonicalStrainModeDefinition
from backend.method4.strain_modes import apply_canonical_strain_basis
from backend.method4.strain_modes import compute_homogeneous_strain_modes
from backend.method4.strain_modes import engineering_voigt_to_tensor
from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import MicroscopicVectorBlock
from backend.models.iso_mode_models import MicroscopicVectorRow
from backend.models.iso_mode_models import ModeIdentity
from backend.models.iso_mode_models import microscopic_atom_order_id
from backend.tables.spacegroup_meta import hall_for_number
from backend.utils.affine_embeddings import AffineOperation
from backend.utils.lattice import as_fraction
from backend.utils.lattice import inverse
from backend.utils.lattice import multiply
from backend.utils.lattice import rational_matrix
from backend.utils.lattice import transpose
from backend.wrappers import SubgroupInfo
from backend.wrappers import SubgroupInfo as _formats_SubgroupInfo
from dataclasses import FrozenInstanceError
from dataclasses import replace
from fractions import Fraction
from pymatgen.core import IStructure
from pymatgen.core import Lattice
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from tests.support_backend_export import ANTI_DIAGONAL_MIRROR
from tests.support_backend_export import SKEW_LATTICE
from tests.support_backend_export import _TetragonalMacroIso
from tests.support_backend_export import _anti_coupled_fixture
from tests.support_backend_export import _bare_api
from tests.support_backend_export import _canonical_source
from tests.support_backend_export import _cif_displacive_model_I3
from tests.support_backend_export import _cif_displacive_model_ZERO
from tests.support_backend_export import _cubic
from tests.support_backend_export import _formats_REPO
from tests.support_backend_export import _minimal_export_api
from tests.support_backend_export import _mode
from tests.support_backend_export import _op
from tests.support_backend_export import _sg
from tests.support_backend_export import _strain_export_subgroup
from tests.support_backend_method2 import _microscopic_column_provenance_subgroup
from tests.support_backend_method2 import _sparse_domain_extension_fixture
from tests.support_paths import TESTS_ROOT
from types import SimpleNamespace
import inspect
import io
import json
import math
import numpy as np
import pytest
import re
import spglib
import zipfile


# Original coverage: test_basics.py.

def test_export_cif_uses_cif_wyckoff_lines_ndnio2():
    """Distortion CIF comments must reuse the same CIF-ordered parent header."""

    from application.session import IsoDistort
    from backend.wrappers import SubgroupInfo
    from backend.export.distortion_formats import render_cif

    cris = TESTS_ROOT.parent.parent
    cif = cris / "experiment_data" / "NdNiO2 own.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("NdNiO2 own.cif not in experiment_data")
    iso = IsoDistort()
    iso.load_structure(cif)
    sg = SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        size=2,
        subgroup_index=2,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label="A1+",
        k_point_label="A",
        parent_sg=123,
        k_coordinates=["1/2", "1/2", "1/2"],
        basis_vectors=[[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
        origin=[0.0, 0.0, 0.0],
    )
    spec = iso._spec_for_subgroup(
        sg, use_current_modes=False, use_generated_structure=False
    )
    assert spec.parent_wyckoff_lines
    assert spec.parent_wyckoff_lines[0].startswith("O 2f")
    text = render_cif(spec.structure, spec)
    assert "# O 2f (0,1/2,0)" in text
    assert "# ND 1d (1/2,1/2,1/2)" in text
    assert "# NI 1a (0,0,0)" in text
    assert "# Nd1" not in text

# Original coverage: test_cif_displacive_model.py.

HEXAGONAL_SIXFOLD = ((1, -1, 0), (1, 0, 0), (0, 0, 1))

HEXAGONAL_LATTICE = ((1.0, 0.0, 0.0), (-0.5, math.sqrt(3.0) / 2.0, 0.0), (0.0, 0.0, 2.0))

F_PRIMITIVE = (
    (0, Fraction(1, 2), Fraction(1, 2)),
    (Fraction(1, 2), 0, Fraction(1, 2)),
    (Fraction(1, 2), Fraction(1, 2), 0),
)

def _f_centered_fixture(frame_id: str = "F"):
    half = Fraction(1, 2)
    translations = (
        (0, 0, 0),
        (0, half, half),
        (half, 0, half),
        (half, half, 0),
    )
    operations = tuple(_op(_cif_displacive_model_I3, translation) for translation in translations)
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

    site = site_fixed_coordinate_model("A", _cif_displacive_model_ZERO, operations)

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
    general = site_fixed_coordinate_model("general", _cif_displacive_model_ZERO, operations)

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
    new_atoms = (ChildAtom.from_values("A", transform_fractional_position(_cif_displacive_model_ZERO, transform, origin), "X"),)
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
        conventional_basis=_cif_displacive_model_I3,
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
        _op(_cif_displacive_model_I3, translation)
        for translation in (
            _cif_displacive_model_ZERO,
            (0, half, half),
            (half, 0, half),
            (half, half, 0),
        )
    )

    with pytest.raises(DisplaciveModelError) as caught:
        derive_centering_index(
            operations,
            conventional_basis=_cif_displacive_model_I3,
            primitive_translation_basis=_cif_displacive_model_I3,
        )
    assert caught.value.code == "centering_index_mismatch"

def test_norm_uses_reference_lattice_instead_of_final_strained_lattice() -> None:
    atoms = (ChildAtom.from_values("A", _cif_displacive_model_ZERO, "X"),)
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
        ChildAtom.from_values("A", _cif_displacive_model_ZERO, "Fe"),
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
    atoms = (ChildAtom.from_values("A", _cif_displacive_model_ZERO, "X"),)
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
    atoms = (ChildAtom.from_values("A", _cif_displacive_model_ZERO, "X"),)
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
    atoms = (ChildAtom.from_values("A", _cif_displacive_model_ZERO, "X"),)
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
        (lambda: ChildAtom.from_values("   ", _cif_displacive_model_ZERO, "X"), "invalid_atom_id"),
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

# Original coverage: test_cif_headers.py.

def test_space_hm_symbol_splits_mm2_and_1m1():
    assert _space_hm_symbol("Pmm2") == "P m m 2"
    assert _space_hm_symbol("P1m1") == "P 1 m 1"
    assert _space_hm_symbol("Imm2") == "I m m 2"
    assert _space_hm_symbol("I2/m2/m2/m") == "I 2/m 2/m 2/m"

def test_cif_header_helpers_wire_parent_and_k():
    sg = SubgroupInfo(
        index=0,
        space_group_number=25,
        space_group_symbol="Pmm2",
        irrep_label="A1+",
        k_point_label="A",
        parent_sg=123,
        k_coordinates=["1/2", "1/2", "1/2"],
    )
    assert _k_point_comment(sg) == "k point: A, k20 (1/2,1/2,1/2)"
    assert _irrep_kovalev_tag(sg) == "k20t1"

# Original coverage: test_displacive_export.py.

_displacive_export_I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

_displacive_export_ZERO = (0, 0, 0)

def _site_mapping(
    child_atom_id: str,
    parent_site_index: int,
    parent_orbit_id: str,
    translation=_displacive_export_ZERO,
) -> ParentChildSiteMapping:
    return ParentChildSiteMapping(
        child_atom_id=child_atom_id,
        parent_site_index=parent_site_index,
        parent_orbit_id=parent_orbit_id,
        parent_cell_translation=translation,
    )

def _subgroup_identity(
    embedding: ExactParentChildEmbedding,
    *,
    target_space_group_number: int = 1,
    target_hall_number: int | None = None,
    irrep_label: str = "GM1+",
    opd_symbol: str = "P1",
    primary_direction_selector: str | None = "VECTOR,A",
    context_kind: str = "single_irrep",
) -> DisplaciveSubgroupIdentity:
    return DisplaciveSubgroupIdentity(
        embedding,
        target_space_group_number=target_space_group_number,
        target_hall_number=target_hall_number,
        irrep_label=(None if context_kind == "exact_fixed_space" else irrep_label),
        opd_symbol=(None if context_kind == "exact_fixed_space" else opd_symbol),
        primary_direction_selector=primary_direction_selector,
        context_kind=context_kind,
    )

def _resolved_mode(
    *,
    axis: int,
    frame_id: str,
    atom_ids: tuple[str, ...],
    displacements,
    global_irrep: str,
    parent_sg: int,
    wyckoff_letter: str,
    orbit_id: str | None = None,
    column_index: int | None = None,
    source_child_sg: int = 1,
    source_basis=_displacive_export_I3,
    source_origin=_displacive_export_ZERO,
    source_subgroup_irrep: str = "GM1+",
    source_subgroup_opd: str = "P1",
    context_kind: str = "single_irrep",
) -> RawModeColumn:
    source_order = axis if column_index is None else column_index
    site_irrep = f"A{axis + 1}"
    block = MicroscopicVectorBlock(
        global_irrep=global_irrep,
        wyckoff_letter=wyckoff_letter,
        site_irrep=site_irrep,
        rows=(
            MicroscopicVectorRow(
                point_raw=("0", "0", "0"),
                displacements=((Fraction(1), Fraction(0), Fraction(0)),),
            ),
        ),
        source_order=source_order,
        direction_symbol="VECTOR,A",
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=f"{source_order:064x}",
        query_order=source_order,
        query_irrep_label=global_irrep,
        source_parent_sg=parent_sg,
        source_k_coordinates=("0", "0", "0"),
        source_subgroup_space_group_number=source_child_sg,
        source_subgroup_basis=source_basis,
        source_subgroup_origin=source_origin,
        source_subgroup_irrep_label=(
            None if context_kind == "exact_fixed_space"
            else source_subgroup_irrep
        ),
        source_subgroup_opd_symbol=(
            None if context_kind == "exact_fixed_space"
            else source_subgroup_opd
        ),
        source_subgroup_primary_direction_selector=(
            None if context_kind == "exact_fixed_space" else "VECTOR,A"
        ),
        source_subgroup_context_kind=context_kind,
        block=block,
        column_index=0,
        source_frame_id=frame_id,
        source_atom_order_id=microscopic_atom_order_id(atom_ids),
    )
    identity = ModeIdentity(
        parent_sg=parent_sg,
        global_irrep=global_irrep,
        k_coordinates=("0", "0", "0"),
        wyckoff_letter=wyckoff_letter,
        orbit_id=orbit_id or f"orbit-{wyckoff_letter}",
        site_irrep=site_irrep,
        component_index=0,
        component_label=f"c{axis + 1}",
        source="iso_microscopic",
        status="verified",
    )
    return RawModeColumn.from_values(
        identity.stable_token,
        source_order,
        frame_id,
        atom_ids,
        displacements,
        ModeScaleProvenance(
            source="iso_microscopic_display_distortion",
            convention="unmixed_fractional_source_column",
            direction_resolved=True,
            sign_resolved=True,
        ),
        mode_identity=identity,
        microscopic_provenance=provenance,
        label=f"mode-{source_order}",
    )

def _export_data(
    *,
    final_x: float | None = None,
    global_irreps: tuple[str, ...] = ("GM1+", "GM1+", "GM1+"),
    context_kind: str = "single_irrep",
) -> DisplaciveExportData:
    operation = ExactSeitzOperation.from_values(_displacive_export_I3, _displacive_export_ZERO)
    frame = ExactChildFrame.from_values(
        "p1_exact",
        (operation,),
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atom = ChildAtom.from_values(
        "a0", (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)), "Fe"
    )
    modes = tuple(
        _resolved_mode(
            axis=axis,
            frame_id=frame.frame_id,
            atom_ids=(atom.atom_id,),
            displacements=(tuple(float(value) for value in np.eye(3)[axis]),),
            global_irrep=global_irreps[axis],
            parent_sg=1,
            wyckoff_letter="a",
            context_kind=context_kind,
        )
        for axis in range(3)
    )
    model = build_displacive_cif_model(
        frame,
        (atom,),
        (atom.atom_id,),
        modes,
        amplitudes=(0.4, 0.0, 0.0),
    )
    reference = Structure(
        Lattice.cubic(4.0),
        ["Fe"],
        [[0.25, 0.5, 0.75]],
        labels=["Fe1"],
    )
    final = None
    if final_x is not None:
        final = Structure(
            reference.lattice,
            ["Fe"],
            [[final_x, 0.5, 0.75]],
            labels=["Fe1"],
        )
    parent = Structure.from_sites(reference)
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        _displacive_export_I3,
        _displacive_export_ZERO,
        (_site_mapping("a0", 0, "orbit-a"),),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(
            embedding,
            primary_direction_selector=(
                None if context_kind == "exact_fixed_space" else "VECTOR,A"
            ),
            context_kind=context_kind,
        ),
        reference_structure=reference,
        final_structure=final,
        representative_labels={"a0": "Fe1_1"},
        parent_orbit_types={
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
        },
    )

def _spec(data: DisplaciveExportData) -> SubgroupExportSpec:
    exact_fixed_space = (
        data.subgroup_identity.context_kind == "exact_fixed_space"
    )
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=1,
        space_group_symbol="P1",
        size=1,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
        irrep_label="" if exact_fixed_space else "GM1+",
        opd_symbol="" if exact_fixed_space else "P1",
        opd_dir_raw="" if exact_fixed_space else "(a)",
        parent_sg=1,
    )
    if exact_fixed_space:
        subgroup.method3_route_resolution = "exact_fixed_space"
    subgroup._displacive_embedding_id = data.embedding.embedding_id
    return SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )

def _two_same_species_orbit_export_data() -> DisplaciveExportData:
    operation = ExactSeitzOperation.from_values(_displacive_export_I3, _displacive_export_ZERO)
    frame = ExactChildFrame.from_values(
        "p1_two_orbits",
        (operation,),
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atoms = (
        ChildAtom.from_values("a0", (Fraction(1, 4),) * 3, "Fe"),
        ChildAtom.from_values("b0", (Fraction(3, 4),) * 3, "Fe"),
    )
    atom_ids = tuple(atom.atom_id for atom in atoms)
    modes = []
    for orbit_position, (orbit_id, atom_position) in enumerate(
        (("orbit-a", 0), ("orbit-b", 1))
    ):
        for axis in range(3):
            rows = np.zeros((2, 3), dtype=float)
            rows[atom_position, axis] = 1.0
            modes.append(
                _resolved_mode(
                    axis=axis,
                    frame_id=frame.frame_id,
                    atom_ids=atom_ids,
                    displacements=rows,
                    global_irrep="GM1+",
                    parent_sg=1,
                    wyckoff_letter="a",
                    orbit_id=orbit_id,
                    column_index=3 * orbit_position + axis,
                )
            )
    model = build_displacive_cif_model(
        frame,
        atoms,
        atom_ids,
        tuple(modes),
        amplitudes=(0.0,) * 6,
    )
    reference = Structure(
        Lattice.cubic(4.0),
        ["Fe", "Fe"],
        [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
        labels=["Fe1", "Fe2"],
    )
    parent = Structure.from_sites(reference)
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        _displacive_export_I3,
        _displacive_export_ZERO,
        (
            _site_mapping("a0", 0, "orbit-a"),
            _site_mapping("b0", 1, "orbit-b"),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels={"a0": "Fe1_1", "b0": "Fe2_1"},
        parent_orbit_types={
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
            "orbit-b": ParentOrbitType(2, "Fe2", ("b0",)),
        },
    )

def _index_two_export_data() -> DisplaciveExportData:
    basis = ((2, 0, 0), (0, 1, 0), (0, 0, 1))
    operation = ExactSeitzOperation.from_values(_displacive_export_I3, _displacive_export_ZERO)
    frame = ExactChildFrame.from_values(
        "p1_index_two",
        (operation,),
        ((8.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atoms = (
        ChildAtom.from_values(
            "a0", (Fraction(1, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
        ChildAtom.from_values(
            "a1", (Fraction(5, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
    )
    atom_ids = tuple(atom.atom_id for atom in atoms)
    modes = []
    for atom_position, global_irrep in enumerate(("GM1+", "GM2+")):
        for axis in range(3):
            rows = np.zeros((2, 3), dtype=float)
            rows[atom_position, axis] = 1.0
            modes.append(
                _resolved_mode(
                    axis=axis,
                    frame_id=frame.frame_id,
                    atom_ids=atom_ids,
                    displacements=rows,
                    global_irrep=global_irrep,
                    parent_sg=1,
                    wyckoff_letter="a",
                    orbit_id="orbit-a",
                    column_index=3 * atom_position + axis,
                    source_basis=basis,
                )
            )
    model = build_displacive_cif_model(
        frame,
        atoms,
        atom_ids,
        tuple(modes),
        amplitudes=(0.4, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
    parent = Structure(
        Lattice.cubic(4.0),
        ["Fe"],
        [[0.25, 0.25, 0.25]],
        labels=["Fe1"],
    )
    reference = Structure(
        Lattice.orthorhombic(8.0, 4.0, 4.0),
        ["Fe", "Fe"],
        [[0.125, 0.25, 0.25], [0.625, 0.25, 0.25]],
        labels=["Fe1", "Fe1"],
    )
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        basis,
        _displacive_export_ZERO,
        (
            _site_mapping("a0", 0, "orbit-a", (0, 0, 0)),
            _site_mapping("a1", 0, "orbit-a", (1, 0, 0)),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels={"a0": "Fe1_1", "a1": "Fe1_2"},
        parent_orbit_types={
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0", "a1")),
        },
    )

def _subgroup_for_data(data: DisplaciveExportData) -> SubgroupInfo:
    identity = data.subgroup_identity
    exact_fixed_space = identity.context_kind == "exact_fixed_space"
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=identity.target_space_group_number,
        space_group_symbol="P1",
        size=identity.size,
        basis_vectors=[
            [value for value in row] for row in identity.basis
        ],
        origin=[value for value in identity.origin],
        irrep_label=identity.irrep_label or "",
        opd_symbol=identity.opd_symbol or "",
        opd_dir_raw="" if exact_fixed_space else "(a)",
        parent_sg=identity.parent_space_group_number,
    )
    if exact_fixed_space:
        subgroup.method3_route_resolution = "exact_fixed_space"
    subgroup._displacive_embedding_id = identity.embedding_id
    return subgroup

def test_shared_export_contract_separates_reference_and_final_coordinates() -> None:
    data = _export_data()
    spec = _spec(data)

    assert spec.structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert spec.cif_structure is not None
    assert spec.cif_structure.frac_coords[0, 0] == pytest.approx(0.35)
    assert tuple(spec.mode_displacements_sc or {}) == data.mode_keys
    assert spec.amplitudes == dict(zip(data.mode_keys, (0.4, 0.0, 0.0), strict=True))
    assert spec.mode_global_irreps == {
        key: "GM1+" for key in data.mode_keys
    }
    assert spec.mode_keys == data.mode_keys
    assert spec.atom_ids == ("a0",)
    assert spec.parent_sg == 1
    assert spec.parent_symbol == "P1"
    assert spec.parent_orbit_types is data.parent_orbit_types
    rows = data.modes()
    assert {row.parent_orbit_id for row in rows} == {"orbit-a"}
    assert {row.parent_type_index for row in rows} == {1}
    assert {row.parent_type_label for row in rows} == {"Fe1"}

    with pytest.raises(DisplaciveModelError) as caught:
        _export_data(final_x=0.25)
    assert caught.value.code == "final_displacement_mismatch"

def test_parent_child_mapping_proves_complete_index_two_translation_quotient() -> None:
    data = _index_two_export_data()
    subgroup = _subgroup_for_data(data)

    assert data.embedding.conventional_translation_index == 2
    assert data.subgroup_identity.size == 2
    assert tuple(
        mapping.parent_cell_translation
        for mapping in data.embedding.atom_mappings
    ) == ((0, 0, 0), (1, 0, 0))
    with pytest.raises(FrozenInstanceError):
        data.embedding.atom_mappings[0].parent_site_index = 1  # type: ignore[misc]
    data.validate_subgroup(subgroup)

    spec = SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    complete = render_complete_modes(spec)
    assert "As=  0.400000" in complete
    assert "Ap=  0.282843" in complete

def test_parent_child_mapping_rejects_duplicate_translation_coset() -> None:
    source = _index_two_export_data()
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        _displacive_export_ZERO,
        (
            _site_mapping("a0", 0, "orbit-a", (0, 0, 0)),
            _site_mapping("a1", 0, "orbit-a", (0, 0, 0)),
        ),
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=embedding,
            subgroup_identity=_subgroup_identity(embedding),
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "duplicate_parent_child_translation_coset"

def test_parent_child_mapping_rejects_identity_basis_zero_origin_global_shift() -> None:
    source = _export_data()
    shifted_atom = ChildAtom.from_values(
        "a0", (Fraction(3, 8), Fraction(1, 2), Fraction(3, 4)), "Fe"
    )
    shifted_model = build_displacive_cif_model(
        source.model.frame,
        (shifted_atom,),
        (shifted_atom.atom_id,),
        source.model.canonical_modes,
        amplitudes=source.model.amplitudes,
    )
    shifted_reference = Structure(
        source.reference_structure.lattice,
        ["Fe"],
        [[0.375, 0.5, 0.75]],
        labels=["Fe1"],
    )
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        _displacive_export_I3,
        _displacive_export_ZERO,
        source.embedding.atom_mappings,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=shifted_model,
            embedding=embedding,
            subgroup_identity=_subgroup_identity(embedding),
            reference_structure=shifted_reference,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "parent_child_coordinate_mismatch"

def test_parent_child_mapping_rejects_chemical_or_orbit_substitution() -> None:
    source = _export_data()
    wrong_chemistry_parent = Structure(
        source.parent_structure.lattice,
        ["Co"],
        source.parent_structure.frac_coords,
        labels=["Co1"],
    )
    chemical_embedding = ExactParentChildEmbedding(
        wrong_chemistry_parent,
        1,
        _displacive_export_I3,
        _displacive_export_ZERO,
        source.embedding.atom_mappings,
    )
    with pytest.raises(DisplaciveModelError) as chemical:
        DisplaciveExportData(
            model=source.model,
            embedding=chemical_embedding,
            subgroup_identity=_subgroup_identity(chemical_embedding),
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )
    assert chemical.value.code == "parent_child_chemistry_mismatch"

    orbit_embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        _displacive_export_I3,
        _displacive_export_ZERO,
        (_site_mapping("a0", 0, "another-orbit"),),
    )
    with pytest.raises(DisplaciveModelError) as orbit:
        DisplaciveExportData(
            model=source.model,
            embedding=orbit_embedding,
            subgroup_identity=_subgroup_identity(orbit_embedding),
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )
    assert orbit.value.code == "parent_child_orbit_mapping_mismatch"

def test_subgroup_identity_derives_size_and_rejects_every_mutable_identity_field() -> None:
    data = _index_two_export_data()
    assert "size" not in inspect.signature(DisplaciveSubgroupIdentity).parameters

    mutations = (
        ("irrep_label", "GM2+", "subgroup_identity_mismatch"),
        ("opd_symbol", "P2", "subgroup_identity_mismatch"),
        ("opd_dir_raw", "(0,a)", "subgroup_identity_mismatch"),
        ("size", 1, "subgroup_identity_mismatch"),
        ("parent_sg", 2, "subgroup_embedding_mismatch"),
        ("space_group_number", 2, "subgroup_identity_mismatch"),
    )
    for field, value, expected_code in mutations:
        subgroup = _subgroup_for_data(data)
        setattr(subgroup, field, value)
        with pytest.raises(DisplaciveModelError) as caught:
            data.validate_subgroup(subgroup)
        assert caught.value.code == expected_code

    subgroup = _subgroup_for_data(data)
    subgroup.basis_vectors[0][0] = 3
    with pytest.raises(DisplaciveModelError) as basis:
        data.validate_subgroup(subgroup)
    assert basis.value.code == "subgroup_embedding_mismatch"

    subgroup = _subgroup_for_data(data)
    subgroup.origin[0] = 1
    with pytest.raises(DisplaciveModelError) as origin:
        data.validate_subgroup(subgroup)
    assert origin.value.code == "subgroup_embedding_mismatch"

    subgroup = _subgroup_for_data(data)
    subgroup._displacive_embedding_id = "forged-embedding"
    with pytest.raises(DisplaciveModelError) as embedding_id:
        data.validate_subgroup(subgroup)
    assert embedding_id.value.code == "subgroup_embedding_mismatch"

def test_export_query_context_distinguishes_exact_embedding_and_primary_direction() -> None:
    source = _export_data()

    def _model_with_provenance_change(**changes):
        modes = tuple(
            replace(
                mode,
                microscopic_provenance=replace(
                    mode.microscopic_provenance,
                    **changes,
                ),
            )
            for mode in source.model.canonical_modes
        )
        return build_displacive_cif_model(
            source.model.frame,
            source.model.atoms,
            tuple(orbit.representative_atom_id for orbit in source.model.orbits),
            modes,
            amplitudes=source.model.amplitudes,
        )

    context_changes = (
        {"source_subgroup_space_group_number": 2},
        {"source_subgroup_basis": (("2", "0", "0"), ("0", "1", "0"), ("0", "0", "1"))},
        {"source_subgroup_origin": ("1", "0", "0")},
        {"source_subgroup_irrep_label": "GM2+"},
        {"source_subgroup_opd_symbol": "P2"},
        {"source_subgroup_primary_direction_selector": "VECTOR,0,A"},
    )
    for changes in context_changes:
        with pytest.raises(DisplaciveModelError) as caught:
            DisplaciveExportData(
                model=_model_with_provenance_change(**changes),
                embedding=source.embedding,
                subgroup_identity=source.subgroup_identity,
                reference_structure=source.reference_structure,
                final_structure=source.final_structure,
                representative_labels=source.representative_labels,
                parent_orbit_types=source.parent_orbit_types,
            )
        assert caught.value.code == "microscopic_query_context_mismatch"

def test_export_query_context_does_not_require_a_primary_displacement_block() -> None:
    data = _export_data(global_irreps=("GM2+", "GM2+", "GM2+"))
    assert data.subgroup_identity.primary_direction_selector == "VECTOR,A"
    assert {row.global_irrep_label for row in data.modes()} == {"GM2+"}

def test_export_contract_has_no_parallel_mode_or_atom_identity_inputs() -> None:
    parameters = inspect.signature(DisplaciveExportData).parameters
    assert not {
        "mode_keys",
        "mode_labels",
        "mode_global_irreps",
        "atom_ids",
    }.intersection(parameters)

def test_parent_orbit_type_contract_rejects_unmapped_mode_orbit() -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types={
                "another-orbit": ParentOrbitType(1, "Fe1", ("a0",)),
            },
        )

    assert caught.value.code == "missing_parent_orbit_type"

@pytest.mark.parametrize(
    "parent_orbit_types",
    [
        {
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
            "orbit-b": ParentOrbitType(1, "Fe2", ("other",)),
        },
        {
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
            "orbit-b": ParentOrbitType(2, "Fe1", ("other",)),
        },
    ],
)
def test_parent_orbit_type_contract_requires_unique_indices_and_labels(
    parent_orbit_types: dict[str, ParentOrbitType],
) -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=parent_orbit_types,
        )

    assert caught.value.code == "invalid_parent_orbit_type"

def test_parent_orbit_type_contract_requires_exact_child_atom_partition() -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types={
                "orbit-a": ParentOrbitType(1, "Fe1", ("unknown",)),
            },
        )

    assert caught.value.code == "parent_orbit_atom_mismatch"

def test_mode_source_rejects_empty_global_irrep_before_export_contract() -> None:
    with pytest.raises(ValueError, match="query irrep cannot be empty"):
        _export_data(global_irreps=("GM1+", "", "GM1+"))

def test_production_writers_fail_closed_without_verified_displacive_data() -> None:
    reference = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
    spec = SubgroupExportSpec(
        subgroup=SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1"),
        structure=reference,
        mode_displacements_sc={"anonymous": np.array([[1.0, 0.0, 0.0]])},
        require_verified_displacive_data=True,
    )

    for writer in (
        lambda: render_cif(reference, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    ):
        with pytest.raises(ValueError, match="exact emitted-setting frame"):
            writer()

def test_four_writers_share_validated_modes_without_double_counting() -> None:
    spec = _spec(_export_data())

    cif = render_cif(spec.cif_structure, spec)
    complete = render_complete_modes(spec)
    isoviz = render_isoviz(spec)
    topas = render_topas(spec)

    assert "0.35000  0.50000  0.75000" in cif
    assert '"1/4 + Fe1_1_dx"' in cif
    assert "_iso_displacivemodenorm_value" in cif
    assert "_iso_displacivemodematrix_value" in cif
    assert "Undistorted superstructure:" in complete
    assert "Distorted superstructure:" in complete
    assert "  0.250000   0.500000   0.750000" in complete
    assert "  0.350000   0.500000   0.750000" in complete
    assert "0.25000   0.00000   0.00000" in isoviz
    assert "+  0.25000*a1" in topas

def test_exact_fixed_space_four_writers_use_real_irreps_without_fake_primary():
    data = _export_data(
        global_irreps=("GM1+", "GM2+", "GM3+"),
        context_kind="exact_fixed_space",
    )
    spec = _spec(data)
    spec.subgroup.basis_raw = "BOGUS-BASIS"
    spec.subgroup.origin_raw = "BOGUS-ORIGIN"

    cif = render_cif(spec.cif_structure, spec)
    complete = render_complete_modes(spec)
    isoviz = render_isoviz(spec)
    topas = render_topas(spec)

    assert "Exact fixed-space coupled embedding" in cif
    assert "basis=(1,0,0),(0,1,0),(0,0,1), origin=(0,0,0)" in cif
    assert "[0,0,0]" not in cif.split("# Order parameter values:", 1)[1]
    assert "exact fixed-space coupled embedding" in complete
    assert "IR=None" not in complete
    irrep_block = isoviz.split("!irreplist ", 1)[1].split("\n\n", 1)[0]
    assert "IR" not in irrep_block.split()
    assert all(irrep in irrep_block for irrep in ("GM1+", "GM2+", "GM3+"))
    assert "'{{{mode definitions" in topas
    assert all(
        "BOGUS-BASIS" not in payload and "BOGUS-ORIGIN" not in payload
        for payload in (cif, complete, isoviz, topas)
    )

    spec.subgroup.opd_symbol = "BOGUS"
    for writer in (
        lambda: render_cif(spec.cif_structure, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    ):
        with pytest.raises(
            DisplaciveModelError,
            match="incorrectly claims a primary route",
        ):
            writer()

    spec.subgroup.opd_symbol = ""
    spec.subgroup.k_point_label = "X"
    spec.subgroup.k_coordinates = ["1/2", "0", "0"]
    for writer in (
        lambda: render_cif(spec.cif_structure, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    ):
        with pytest.raises(
            DisplaciveModelError,
            match="incorrectly claims a primary route",
        ):
            writer()

def test_exact_fixed_space_export_rejects_mixed_primary_context_provenance():
    source = _export_data(context_kind="exact_fixed_space")
    first = source.model.canonical_modes[0]
    mixed_provenance = replace(
        first.microscopic_provenance,
        source_subgroup_context_kind="single_irrep",
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
    )
    mixed_modes = (
        replace(first, microscopic_provenance=mixed_provenance),
        *source.model.canonical_modes[1:],
    )
    mixed_model = build_displacive_cif_model(
        source.model.frame,
        source.model.atoms,
        tuple(orbit.representative_atom_id for orbit in source.model.orbits),
        mixed_modes,
        amplitudes=source.model.amplitudes,
    )

    with pytest.raises(
        DisplaciveModelError,
        match=r"query context differs|mix or mismatch",
    ):
        DisplaciveExportData(
            model=mixed_model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            final_structure=source.final_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

def test_isoviz_assigns_stable_numbers_to_two_explicit_displacive_irreps() -> None:
    spec = _spec(
        _export_data(global_irreps=("GM1+", "GM2+", "GM1+"))
    )
    text = render_isoviz(spec)

    irrep_block = text.split("!irreplist ", 1)[1].split("\n\n", 1)[0]
    irrep_rows = [
        (int(parts[0]), parts[1])
        for line in irrep_block.splitlines()
        if len(parts := line.split()) == 2 and parts[0].isdigit()
    ]
    assert irrep_rows == [(1, "GM1+"), (2, "GM2+")]

    mode_block = text.split("!displacivemodelist ", 1)[1].split("\n\n", 1)[0]
    mode_rows = [
        parts
        for line in mode_block.splitlines()
        if len(parts := line.split()) >= 6
        and parts[0].isdigit()
        and parts[1].isdigit()
    ]
    assert [int(parts[4]) for parts in mode_rows] == [1, 2, 1]

def test_isoviz_routes_same_species_modes_by_parent_orbit_contract() -> None:
    spec = _spec(_two_same_species_orbit_export_data())

    text = render_isoviz(spec)

    atom_types = text.split("!atomtypelist ", 1)[1].split("\n\n", 1)[0]
    assert "1 Fe1 Fe" in atom_types
    assert "2 Fe2 Fe" in atom_types
    mode_block = text.split("!displacivemodelist ", 1)[1].split("\n\n", 1)[0]
    mode_rows = [
        parts
        for line in mode_block.splitlines()
        if len(parts := line.split()) >= 6
        and parts[0].isdigit()
        and parts[1].isdigit()
    ]
    assert [int(parts[0]) for parts in mode_rows] == [1, 1, 1, 2, 2, 2]

def test_export_data_rejects_final_lattice_outside_validated_frame() -> None:
    data = _export_data()
    assert data.final_structure is not None
    wrong_final = Structure(
        Lattice.tetragonal(9.0, 11.0),
        [site.species for site in data.final_structure],
        data.final_structure.frac_coords,
        labels=[site.label for site in data.final_structure],
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=data.model,
            embedding=data.embedding,
            subgroup_identity=data.subgroup_identity,
            reference_structure=data.reference_structure,
            final_structure=wrong_final,
            representative_labels=data.representative_labels,
            parent_orbit_types=data.parent_orbit_types,
        )

    assert caught.value.code == "final_frame_mismatch"

def _pminus_one_export_data() -> DisplaciveExportData:
    identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    inversion = ((-1, 0, 0), (0, -1, 0), (0, 0, -1))
    operations = (
        ExactSeitzOperation.from_values(identity, (0, 0, 0)),
        ExactSeitzOperation.from_values(inversion, (0, 0, 0)),
    )
    frame = ExactChildFrame.from_values(
        "pminus1_exact",
        operations,
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atoms = (
        ChildAtom.from_values("first", (Fraction(1, 5),) * 3, "Fe"),
        ChildAtom.from_values("representative", (Fraction(4, 5),) * 3, "Fe"),
    )
    modes = tuple(
        _resolved_mode(
            axis=axis,
            frame_id=frame.frame_id,
            atom_ids=tuple(atom.atom_id for atom in atoms),
            displacements=(-np.eye(3)[axis], np.eye(3)[axis]),
            global_irrep="GM1+",
            parent_sg=2,
            wyckoff_letter="i",
            source_child_sg=2,
        )
        for axis in range(3)
    )
    model = build_displacive_cif_model(
        frame,
        atoms,
        ("representative",),
        modes,
        amplitudes=(0.1, 0.0, 0.0),
    )
    reference = Structure(
        Lattice.cubic(4.0),
        ["Fe", "Fe"],
        [[0.2, 0.2, 0.2], [0.8, 0.8, 0.8]],
        labels=["first", "representative"],
    )
    parent = Structure.from_sites(reference)
    embedding = ExactParentChildEmbedding(
        parent,
        2,
        _displacive_export_I3,
        _displacive_export_ZERO,
        (
            _site_mapping("first", 0, "orbit-i"),
            _site_mapping("representative", 1, "orbit-i"),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(
            embedding, target_space_group_number=2
        ),
        reference_structure=reference,
        representative_labels={"representative": "Fe1_1"},
        parent_orbit_types={
            "orbit-i": ParentOrbitType(
                1,
                "Fe1",
                ("first", "representative"),
            ),
        },
    )

def _pminus_one_nonzero_origin_export_data() -> DisplaciveExportData:
    """Return an exact child already expressed after a nonzero B,q embedding."""

    source = _pminus_one_export_data()
    source_origin = ("1/2", "0", "0")
    source_modes = tuple(
        replace(
            mode,
            microscopic_provenance=replace(
                mode.microscopic_provenance,
                source_subgroup_origin=source_origin,
            ),
        )
        for mode in source.model.canonical_modes
    )
    model = build_displacive_cif_model(
        source.model.frame,
        source.model.atoms,
        ("first",),
        source_modes,
        amplitudes=source.model.amplitudes,
    )
    parent = Structure(
        source.parent_structure.lattice,
        [site.species for site in source.parent_structure],
        [[0.7, 0.2, 0.2], [0.3, 0.8, 0.8]],
        labels=["first", "representative"],
    )
    embedding = ExactParentChildEmbedding(
        parent,
        2,
        _displacive_export_I3,
        (Fraction(1, 2), Fraction(0), Fraction(0)),
        (
            _site_mapping("first", 0, "orbit-i"),
            _site_mapping("representative", 1, "orbit-i", (1, 0, 0)),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(
            embedding, target_space_group_number=2
        ),
        reference_structure=source.reference_structure,
        representative_labels={"first": "Fe1_1"},
        parent_orbit_types=source.parent_orbit_types,
    )

def test_cif_uses_exact_contract_frame_without_reapplying_subgroup_origin() -> None:
    data = _pminus_one_nonzero_origin_export_data()
    spec = SubgroupExportSpec(
        subgroup=_subgroup_for_data(data),
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    setting = _setting(spec.subgroup, force_p1=False)
    _legacy_structure, legacy_shift = _apply_origin_choice(
        spec.cif_structure,
        setting,
        (-0.5, 0.0, 0.0),
    )

    assert legacy_shift == pytest.approx((0.5, 0.0, 0.0))
    _setting_value, shifted, _sites, writer_shift = _validated_subgroup_sites(
        spec,
        spec.cif_structure,
    )
    assert writer_shift == pytest.approx((0.0, 0.0, 0.0))
    assert shifted.frac_coords == pytest.approx(spec.cif_structure.frac_coords)

    cif = render_cif(spec.cif_structure, spec)
    assert "_iso_parent-to-child.transform_Pp_abc a,b,c;1/2,0,0" in cif

def test_cif_rejects_unrecorded_shift_after_exact_contract_frame() -> None:
    data = _pminus_one_nonzero_origin_export_data()
    spec = SubgroupExportSpec(
        subgroup=_subgroup_for_data(data),
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    setting, _shifted, sites, _writer_shift = _validated_subgroup_sites(
        spec,
        spec.cif_structure,
    )

    with pytest.raises(ValueError, match="unrecorded writer origin shift"):
        _validate_displacive_writer_frame(
            spec,
            setting,
            sites,
            np.asarray((0.5, 0.0, 0.0)),
        )

def test_cif_fails_closed_when_writer_selects_an_equivalent_nonrepresentative_atom() -> None:
    data = _pminus_one_export_data()
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=2,
        space_group_symbol="P-1",
        size=1,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        parent_sg=2,
    )
    subgroup._displacive_embedding_id = data.embedding.embedding_id
    spec = SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )

    with pytest.raises(ValueError, match="representative differs"):
        render_cif(spec.cif_structure, spec)

def test_cif_binds_contract_labels_when_coordinate_sort_reverses_orbit_order() -> None:
    source = _index_two_export_data()
    atoms = (
        ChildAtom.from_values(
            "a0", (Fraction(5, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
        ChildAtom.from_values(
            "a1", (Fraction(1, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
    )
    model = build_displacive_cif_model(
        source.model.frame,
        atoms,
        ("a0", "a1"),
        source.model.canonical_modes,
        amplitudes=source.model.amplitudes,
    )
    reference = Structure(
        source.reference_structure.lattice,
        ["Fe", "Fe"],
        [[0.625, 0.25, 0.25], [0.125, 0.25, 0.25]],
        labels=["Fe1", "Fe1"],
    )
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        _displacive_export_ZERO,
        (
            _site_mapping("a0", 0, "orbit-a", (1, 0, 0)),
            _site_mapping("a1", 0, "orbit-a", (0, 0, 0)),
        ),
    )
    data = DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels={"a0": "Fe1_1", "a1": "Fe1_2"},
        parent_orbit_types=source.parent_orbit_types,
    )
    spec = SubgroupExportSpec(
        subgroup=_subgroup_for_data(data),
        displacive_data=data,
        require_verified_displacive_data=True,
    )

    cif = render_cif(spec.cif_structure, spec)
    atom_rows = [
        line.split()
        for line in cif.splitlines()
        if line.startswith(("Fe1_1 ", "Fe1_2 ")) and len(line.split()) >= 8
    ]

    assert [(row[0], float(row[4])) for row in atom_rows] == [
        ("Fe1_2", 0.125),
        ("Fe1_1", 0.675),
    ]
    assert 'Fe1_1_x                      "5/8 + Fe1_1_dx"' in cif
    assert 'Fe1_2_x                      "1/8 + Fe1_2_dx"' in cif

def test_all_writers_revalidate_and_reject_changed_subgroup_embedding() -> None:
    data = _export_data()
    spec = _spec(data)
    spec.subgroup.origin = [0.5, 0, 0]

    writers = (
        lambda: render_cif(spec.cif_structure, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    )
    for writer in writers:
        with pytest.raises(DisplaciveModelError) as caught:
            writer()
        assert caught.value.code == "subgroup_embedding_mismatch"

def test_embedding_proves_reference_lattice_is_basis_times_parent_lattice() -> None:
    data = _export_data()
    wrong_embedding = ExactParentChildEmbedding(
        data.parent_structure,
        data.parent_space_group_number,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        _displacive_export_ZERO,
        data.embedding.atom_mappings,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=data.model,
            embedding=wrong_embedding,
            subgroup_identity=_subgroup_identity(wrong_embedding),
            reference_structure=data.reference_structure,
            representative_labels=data.representative_labels,
            parent_orbit_types=data.parent_orbit_types,
        )
    assert caught.value.code == "parent_child_lattice_mismatch"

def test_embedding_parent_space_group_must_match_canonical_mode_identity() -> None:
    data = _export_data()
    wrong_parent = ExactParentChildEmbedding(
        data.parent_structure,
        2,
        _displacive_export_I3,
        _displacive_export_ZERO,
        data.embedding.atom_mappings,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=data.model,
            embedding=wrong_parent,
            subgroup_identity=_subgroup_identity(wrong_parent),
            reference_structure=data.reference_structure,
            representative_labels=data.representative_labels,
            parent_orbit_types=data.parent_orbit_types,
        )
    assert caught.value.code == "mode_parent_space_group_mismatch"

def test_contract_snapshots_structures_mapping_and_rebinds_writer_inputs() -> None:
    source = _export_data()
    parent = Structure.from_sites(source.parent_structure)
    reference = Structure.from_sites(source.reference_structure)
    labels = dict(source.representative_labels)
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        _displacive_export_I3,
        _displacive_export_ZERO,
        source.embedding.atom_mappings,
    )
    data = DisplaciveExportData(
        model=source.model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels=labels,
        parent_orbit_types=source.parent_orbit_types,
    )
    parent.translate_sites([0], [0.2, 0, 0], frac_coords=True)
    reference.translate_sites([0], [0.2, 0, 0], frac_coords=True)
    labels["a0"] = "changed"

    assert isinstance(data.parent_structure, IStructure)
    assert isinstance(data.reference_structure, IStructure)
    assert isinstance(data.final_structure, IStructure)
    assert data.parent_structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert data.reference_structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert data.representative_label("a0") == "Fe1_1"
    with pytest.raises(TypeError):
        data.representative_labels["a0"] = "forbidden"  # type: ignore[index]
    with pytest.raises(TypeError):
        data.parent_orbit_types["orbit-a"] = ParentOrbitType(  # type: ignore[index]
            1,
            "changed",
            ("a0",),
        )

    spec = _spec(data)
    assert spec.mode_displacements_sc is not None
    first_key = data.mode_keys[0]
    spec.structure = Structure(Lattice.cubic(9), ["Fe"], [[0, 0, 0]])
    spec.mode_displacements_sc[first_key][0, 0] = 99.0
    spec.mode_labels[first_key] = "changed"
    spec.parent_symbol = "forged"
    render_complete_modes(spec)
    assert isinstance(spec.structure, IStructure)
    assert spec.structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert spec.mode_displacements_sc[first_key][0, 0] == pytest.approx(1.0)
    assert spec.mode_labels[first_key] == "mode-0"
    assert spec.parent_symbol == "P1"

@pytest.mark.parametrize("target", ["parent", "reference", "final"])
def test_contract_rejects_non_single_species_or_partial_occupancy(target: str) -> None:
    source = _export_data()
    disordered = Structure(
        source.reference_structure.lattice,
        [{"Fe": 0.5, "Co": 0.5}],
        source.reference_structure.frac_coords,
        labels=["Fe1"],
    )
    parent = disordered if target == "parent" else source.parent_structure
    if target == "parent":
        with pytest.raises(DisplaciveModelError) as caught:
            ExactParentChildEmbedding(
                parent,
                1,
                _displacive_export_I3,
                _displacive_export_ZERO,
                source.embedding.atom_mappings,
            )
    else:
        with pytest.raises(DisplaciveModelError) as caught:
            DisplaciveExportData(
                model=source.model,
                embedding=source.embedding,
                subgroup_identity=source.subgroup_identity,
                reference_structure=(
                    disordered if target == "reference" else source.reference_structure
                ),
                final_structure=(disordered if target == "final" else None),
                representative_labels=source.representative_labels,
                parent_orbit_types=source.parent_orbit_types,
            )
    assert caught.value.code == "unsupported_site_occupancy"

def test_spec_rejects_parallel_inputs_when_displacive_contract_is_present() -> None:
    data = _export_data()
    with pytest.raises(ValueError, match="only source for structure"):
        SubgroupExportSpec(
            subgroup=_spec(data).subgroup,
            structure=data.reference_structure,
            displacive_data=data,
        )

    with pytest.raises(ValueError, match="only source for parent_symbol"):
        SubgroupExportSpec(
            subgroup=_spec(data).subgroup,
            parent_symbol="forged",
            displacive_data=data,
        )

def test_direct_isoviz_entry_rejects_removed_contract_even_with_empty_views() -> None:
    spec = _spec(_export_data())
    spec.displacive_data = None
    spec.mode_displacements_sc = {}
    spec.mode_labels = {}
    spec.amplitudes = {}
    spec.require_verified_displacive_data = False

    with pytest.raises(ValueError, match="contract was removed"):
        render_isodistort_isoviz(spec)

def test_coordinate_formula_does_not_reduce_fixed_space_coefficients_modulo_one() -> None:
    identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    sheared_reflection = ((1, 0, 0), (4, -1, 0), (0, 0, 1))
    operations = (
        ExactSeitzOperation.from_values(identity, (0, 0, 0)),
        ExactSeitzOperation.from_values(sheared_reflection, (0, 0, 0)),
    )
    shear = np.array([[1.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    inverse = np.linalg.inv(shear)
    invariant_metric = inverse.T @ inverse
    lattice = np.linalg.cholesky(invariant_metric)
    frame = ExactChildFrame.from_values("sheared_order2", operations, lattice)
    site = site_fixed_coordinate_model("a", (0, 0, 0), frame.operations)

    assert site.symmform == "Dx,2*Dx,Dz"
    parameters = tuple(f"p{index}" for index in range(site.dimension))
    assert _coordinate_formula("0", site.fixed_basis[1], parameters) == "0 + 2*p0"

def _database_operations(hall_number: int) -> tuple[ExactSeitzOperation, ...]:
    database = spglib.get_symmetry_from_database(hall_number)
    assert database is not None
    return tuple(
        ExactSeitzOperation.from_values(rotation, translation)
        for rotation, translation in zip(
            database["rotations"],
            database["translations"],
            strict=True,
        )
    )

def _database_centering_index(operations: tuple[ExactSeitzOperation, ...]) -> int:
    identity = tuple(tuple(Fraction(value) for value in row) for row in _displacive_export_I3)
    return sum(operation.rotation == identity for operation in operations)

@pytest.mark.parametrize(
    "hall_number",
    (1, 424, 523, 185, 186, 188, 458, 459),
)
def test_emitted_space_group_accepts_exact_declared_hall_setting(
    hall_number: int,
) -> None:
    operations = _database_operations(hall_number)
    target = spglib.get_spacegroup_type(hall_number)
    assert target is not None

    _validate_emitted_space_group(
        operations,
        model_centering_index=_database_centering_index(operations),
        target_space_group_number=target.number,
        target_hall_number=hall_number,
    )

@pytest.mark.parametrize(
    ("target_hall", "emitted_hall", "error_code"),
    (
        (1, 2, "target_space_group_type_mismatch"),
        (424, 1, "target_centering_mismatch"),
        (523, 1, "target_centering_mismatch"),
        (185, 186, "target_space_group_setting_mismatch"),
        (186, 188, "target_space_group_setting_mismatch"),
        (188, 185, "target_space_group_setting_mismatch"),
        (458, 459, "target_centering_mismatch"),
        (459, 458, "target_centering_mismatch"),
    ),
)
def test_emitted_space_group_rejects_wrong_centering_type_or_setting(
    target_hall: int,
    emitted_hall: int,
    error_code: str,
) -> None:
    operations = _database_operations(emitted_hall)
    target = spglib.get_spacegroup_type(target_hall)
    assert target is not None

    with pytest.raises(DisplaciveModelError) as caught:
        _validate_emitted_space_group(
            operations,
            model_centering_index=_database_centering_index(operations),
            target_space_group_number=target.number,
            target_hall_number=target_hall,
        )

    assert caught.value.code == error_code

def test_export_contract_rejects_model_with_wrong_target_space_group_type() -> None:
    source = _export_data()
    identity = _subgroup_identity(
        source.embedding,
        target_space_group_number=2,
        target_hall_number=2,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "target_space_group_type_mismatch"

def test_subgroup_identity_rejects_hall_setting_from_another_space_group() -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        _subgroup_identity(
            source.embedding,
            target_space_group_number=1,
            target_hall_number=2,
        )

    assert caught.value.code == "invalid_target_hall_setting"

def test_export_contract_rejects_model_with_wrong_target_centering() -> None:
    source = _export_data()
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        _displacive_export_ZERO,
        source.embedding.atom_mappings,
    )
    identity = _subgroup_identity(
        embedding,
        target_space_group_number=5,
        target_hall_number=9,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=embedding,
            subgroup_identity=identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "target_centering_mismatch"

def test_periodic_cartesian_residual_uses_skew_lattice_nearest_image() -> None:
    lattice = np.asarray(
        ((1.0, 0.0, 0.0), (0.9, 0.1, 0.0), (0.0, 0.0, 1.0)),
        dtype=float,
    )
    first = np.asarray(((0.49, 0.49, 0.0),), dtype=float)
    second = np.zeros((1, 3), dtype=float)

    residual = _periodic_cartesian_residual(first, second, lattice)

    assert residual == pytest.approx(np.hypot(0.031, 0.051))
    componentwise_residual = np.linalg.norm(
        (first[0] - np.round(first[0])) @ lattice
    )
    assert componentwise_residual > 15 * residual

# Original coverage: test_export_fail_closed.py.

def _two_type_layout() -> _AtomLayout:
    return _AtomLayout(
        types=[("A1", "A"), ("B1", "B")],
        child_type_idx=[1, 2],
        image_source_idx=[0, 1],
        orbit_type_idx={"orbit-a": 1, "orbit-b": 2},
    )

def _authoritative_mode(
    rows: np.ndarray,
    *,
    key: str = "mode-1",
    orbit_id: str = "orbit-b",
    type_index: int = 2,
    type_label: str = "B1",
) -> SimpleNamespace:
    return SimpleNamespace(
        key=key,
        parent_orbit_id=orbit_id,
        parent_type_index=type_index,
        parent_type_label=type_label,
        normalized_fractional_per_angstrom=rows,
    )

def test_authoritative_isoviz_mode_rejects_cross_type_column() -> None:
    rows = np.array([[0.25, 0.0, 0.0], [0.0, 0.5, 0.0]])

    with pytest.raises(ValueError, match="spans multiple IsoVIZ parent-atom types"):
        _authoritative_mode_parentatom(
            _authoritative_mode(rows),
            _two_type_layout(),
        )

def test_authoritative_isoviz_mode_rejects_empty_column() -> None:
    with pytest.raises(ValueError, match="has no nonzero child-atom row"):
        _authoritative_mode_parentatom(
            _authoritative_mode(np.zeros((2, 3))),
            _two_type_layout(),
        )

def test_authoritative_isoviz_mode_uses_nonzero_rows_not_display_label() -> None:
    rows = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]])

    assert (
        _authoritative_mode_parentatom(
            _authoritative_mode(rows, key="misleading-[A1:a:dsp]"),
            _two_type_layout(),
        )
        == 2
    )

def test_authoritative_isoviz_mode_rejects_rows_on_wrong_orbit_type() -> None:
    rows = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="expected type 2 for orbit 'orbit-b'"):
        _authoritative_mode_parentatom(
            _authoritative_mode(rows),
            _two_type_layout(),
        )

def test_authoritative_isoviz_mode_rejects_unmapped_orbit_identity() -> None:
    rows = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="absent from the IsoVIZ type layout"):
        _authoritative_mode_parentatom(
            _authoritative_mode(
                rows,
                orbit_id="orbit-missing",
                type_index=1,
                type_label="A1",
            ),
            _two_type_layout(),
        )

# Original coverage: test_formats.py.

def test_parse_export_formats_aliases():
    assert parse_export_formats(None) == [FORMAT_CIF]
    assert parse_export_formats("cif,isoviz,modes,topas") == [
        FORMAT_CIF, FORMAT_ISOVIZ, FORMAT_MODES, FORMAT_TOPAS,
    ]
    assert parse_export_formats("structurefile,completemodesdetails,TOPAS.STR") == [
        FORMAT_CIF, FORMAT_MODES, FORMAT_TOPAS,
    ]
    assert parse_export_formats(["CIF file", "Save interactive distortion"]) == [
        FORMAT_CIF, FORMAT_ISOVIZ,
    ]

def test_parse_export_formats_rejects_unknown():
    with pytest.raises(ValueError, match="未知导出格式"):
        parse_export_formats("poscar")

def test_parse_export_method_single_only():
    assert parse_export_method(None) == 2
    assert parse_export_method("1") == 1
    assert parse_export_method("method2") == 2
    assert parse_export_method("M3") == 3
    with pytest.raises(ValueError, match="不能多选"):
        parse_export_method("1,2")
    with pytest.raises(ValueError, match="未知 Method"):
        parse_export_method("4")

def test_subgroup_label_and_filenames():
    sg = _sg("LD1", "C1")
    assert subgroup_label(sg) == "LD1 C1"
    assert format_filename("LD1 C1", FORMAT_CIF) == "subgroup.cif"
    assert format_filename("LD1 C1", FORMAT_ISOVIZ) == "data.isoviz"
    assert format_filename("LD1 C1", FORMAT_MODES) == "Complete modes details.txt"
    assert format_filename("LD1 C1", FORMAT_TOPAS) == "topas.str"

def test_safe_name_cleans_windows_characters_without_losing_irrep_signs():
    from backend.export.distortion_formats import safe_name

    assert safe_name("N1+/<bad>: 4D1*") == "N1+_bad_4D1"
    assert safe_name("N1- 4D1") == "N1-_4D1"
    assert safe_name("CON") == "_CON"

def test_opd_line_body_omits_irrep_for_cif_comment():
    from backend.utils.opd_format import format_opd_line_body

    body = format_opd_line_body(
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=71,
        space_group_symbol="Immm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=2,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
        pad_opd=True,
    )
    assert body.startswith("P1   (a)  71 Immm")
    assert not body.startswith("GM")

    m2 = format_opd_line_body(
        opd_symbol="C1",
        opd_dir_raw="(a,b)",
        space_group_number=99,
        space_group_symbol="P4mm",
        basis_raw="(1,0,0),(0,1,0),(0,0,6)",
        origin_raw="(0,0,0)",
        size=12,
        subgroup_index=24,
        k_coordinates=["0", "0", "1/6"],
        parent_sg=139,
        pad_opd=False,
    )
    assert m2.startswith("C1 (a,b)  99 P4mm")

def test_short_folder_names_follow_each_method_rule():
    subgroup = _sg("N1+", "4D1", number=2)

    assert folder_label_for_subgroup(subgroup, export_method=1) == "N1+_4D1_SG2"
    assert folder_label_for_subgroup(subgroup, export_method=2) == "N1+_4D1"
    assert (
        folder_label_for_subgroup(subgroup, export_method=3, sequence=7)
        == "C07_SG2"
    )

def test_unique_folder_name_disambiguates_with_stable_identity_digest():
    used: set[str] = set()
    first = _sg("LD1", "C1", index=0, symbol="P4/mmm")
    second = _sg("LD1", "C1", index=1, symbol="P4mm")
    second.origin = [0.5, 0, 0]

    a = unique_folder_name(first, used, export_method=1)
    b = unique_folder_name(second, used, export_method=1)

    assert a == "LD1_C1_SG123"
    assert re.fullmatch(r"LD1_C1_SG123_[0-9a-f]{10}", b)
    assert a != b

    case_insensitive_used = {"ld1_c1"}
    method2 = unique_folder_name(first, case_insensitive_used, export_method=2)
    assert re.fullmatch(r"LD1_C1_[0-9a-f]{10}", method2)

def test_method3_generated_case_id_is_generic_and_order_independent():
    first = _sg("GM1+", "P1", index=0, number=139)
    second = _sg("N1+", "4D1", index=1, number=2)

    generated = method3_case_folder([first, second])
    assert re.fullmatch(r"M3-[0-9a-f]{10}", generated)
    assert generated == method3_case_folder([second, first])
    assert method3_case_folder([first], "case:01") == "case_01"

def test_zip_contains_only_method2_subgroups_named_by_ir_opd():
    """压缩包只含给定子群，按子群文件夹分类；不混入其它文件名。"""
    parent = _cubic()
    disp = {"LD1": np.zeros((2, 3))}
    disp["LD1"][0, 2] = 1.0
    specs = [
        SubgroupExportSpec(
            subgroup=_sg("LD1", "C1", index=0),
            structure=parent,
            parent_structure=parent,
            parent_sg=139,
            parent_symbol="I4/mmm",
            mode_displacements_sc=disp,
            folder_name="LD1_C1",
        ),
        SubgroupExportSpec(
            subgroup=_sg("LD5", "P6", index=1, number=11, symbol="P2_1/m"),
            structure=parent,
            parent_structure=parent,
            parent_sg=139,
            parent_symbol="I4/mmm",
            folder_name="LD5_P6",
        ),
    ]
    raw = build_export_zip(
        specs,
        [FORMAT_CIF, FORMAT_ISOVIZ, FORMAT_MODES, FORMAT_TOPAS],
    )
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()

    # 各子群文件夹在 ZIP 根下（官网同款）
    assert any(n.startswith("LD1_C1/") for n in names)
    assert any(n.startswith("LD5_P6/") for n in names)
    assert "LD1_C1/subgroup.cif" in names
    assert "LD1_C1/data.isoviz" in names
    assert "LD1_C1/Complete modes details.txt" in names
    assert "LD1_C1/topas.str" in names
    assert "LD5_P6/subgroup.cif" in names
    # 不应出现旧版 output_dir 风格的无关文件名
    joined = "\n".join(names)
    assert "mixed_" not in joined
    assert "distorted_" not in joined
    assert ".json" not in joined
    assert len(names) == 8  # 2 子群 × 4 格式

def test_cif_only_zip_has_no_other_formats():
    spec = SubgroupExportSpec(
        subgroup=_sg("LD1", "C1"),
        structure=_cubic(),
        folder_name="LD1_C1",
    )
    raw = build_export_zip([spec], [FORMAT_CIF])
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()
    assert names == ["LD1_C1/subgroup.cif"]

def test_zip_never_emits_duplicate_paths_for_colliding_short_names():
    parent = _cubic()
    first = _sg("LD1", "C1", index=0)
    second = _sg("LD1", "C1", index=1)
    second.origin = [0.5, 0, 0]
    specs = [
        SubgroupExportSpec(first, parent, folder_name="LD1_C1"),
        SubgroupExportSpec(second, parent, folder_name="LD1_C1"),
    ]

    raw = build_export_zip(specs, [FORMAT_CIF])
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()

    assert len(names) == 2
    assert len({name.casefold() for name in names}) == 2
    assert "LD1_C1/subgroup.cif" in names
    assert any(
        re.fullmatch(r"LD1_C1_[0-9a-f]{10}/subgroup\.cif", name)
        for name in names
    )

def test_directory_export_never_overwrites_an_existing_short_folder(tmp_path):
    parent = _cubic()
    subgroup = _sg("N1+", "4D1", index=0, number=2)
    iso = _minimal_export_api(parent)

    first = iso.export_subgroups(
        tmp_path,
        formats=[FORMAT_CIF],
        subgroups=[subgroup],
        export_method=1,
    )
    first_content = first[0].read_bytes()
    second = iso.export_subgroups(
        tmp_path,
        formats=[FORMAT_CIF],
        subgroups=[subgroup],
        export_method=1,
    )

    assert first == [
        tmp_path / "N1+_4D1_SG2" / "subgroup.cif",
        tmp_path / "export_candidate_status.json",
        tmp_path / "export_candidate_status.txt",
    ]
    assert second[0].parent.name == (
        f"N1+_4D1_SG2_{subgroup_identity_digest(subgroup)}"
    )
    assert first[0].read_bytes() == second[0].read_bytes() == first_content
    assert json.loads(first[1].read_text(encoding="utf-8"))["status"] == "success"

def test_directory_writer_render_failure_leaves_no_destination(
    tmp_path,
    monkeypatch,
):
    spec = SubgroupExportSpec(
        subgroup=_sg("LD1", "C1"),
        structure=_cubic(),
        folder_name="LD1_C1",
    )
    dest_dir = tmp_path / "LD1_C1"

    def fail_on_second_format(fmt, _spec):
        if fmt == FORMAT_ISOVIZ:
            raise RuntimeError("later render failed")
        return b"first format"

    monkeypatch.setattr(
        "backend.export.distortion_formats.render_format",
        fail_on_second_format,
    )

    with pytest.raises(RuntimeError, match="later render failed"):
        write_subgroup_files(dest_dir, spec, [FORMAT_CIF, FORMAT_ISOVIZ])

    assert not dest_dir.exists()
    assert list(tmp_path.iterdir()) == []

def test_directory_writer_rejects_existing_destination_without_changes(tmp_path):
    spec = SubgroupExportSpec(
        subgroup=_sg("LD1", "C1"),
        structure=_cubic(),
        folder_name="LD1_C1",
    )
    dest_dir = tmp_path / "LD1_C1"
    dest_dir.mkdir()
    sentinel = dest_dir / "keep.txt"
    sentinel.write_bytes(b"unchanged")

    with pytest.raises(FileExistsError, match="目标目录已存在"):
        write_subgroup_files(dest_dir, spec, [FORMAT_CIF])

    assert sentinel.read_bytes() == b"unchanged"
    assert list(dest_dir.iterdir()) == [sentinel]

def test_directory_writer_rejects_duplicate_filenames_before_publication(tmp_path):
    spec = SubgroupExportSpec(
        subgroup=_sg("LD1", "C1"),
        structure=_cubic(),
        folder_name="LD1_C1",
    )
    dest_dir = tmp_path / "LD1_C1"

    with pytest.raises(ValueError, match="重复文件名"):
        write_subgroup_files(dest_dir, spec, [FORMAT_CIF, FORMAT_CIF])

    assert not dest_dir.exists()
    assert list(tmp_path.iterdir()) == []

def test_format_writers_contain_official_markers():
    parent = _cubic()
    disp = {"LD1": np.array([[0.0, 0.0, 0.2], [0.0, 0.0, 0.0]])}
    spec = SubgroupExportSpec(
        subgroup=_sg("LD1", "C1"),
        structure=parent,
        parent_structure=parent,
        parent_sg=139,
        parent_symbol="I4/mmm",
        mode_displacements_sc=disp,
        mode_labels={"LD1": "[Eu:a:dsp]A2u(a)"},
    )
    cif = render_cif(parent)
    assert "data_isodistort-output" in cif
    assert "_cell_length_a" in cif
    cif_spec = render_cif(parent, spec)
    assert "_iso_parentcell_length_a" in cif_spec
    assert "_space_group_symop_operation_xyz" in cif_spec

    isoviz = render_isoviz(spec)
    assert "!isoversion" in isoviz
    assert "!displacivemodelist" in isoviz
    assert "LD1" in isoviz or "[Eu:a:dsp]A2u(a)" in isoviz
    # IsoVIZ invariant: mode vector count == atomcoordlist rows for parentatom
    lines = isoviz.splitlines()
    coords: list[int] = []
    in_c = False
    for line in lines:
        if line.startswith("!atomcoordlist"):
            in_c = True
            continue
        if line.startswith("!") and in_c:
            in_c = False
        if in_c and line.strip():
            coords.append(int(line.split()[0]))
    i = next(k for k, line in enumerate(lines) if line.startswith("!displacivemodelist"))
    parentatom = None
    nvecs = 0
    for line in lines[i + 1 :]:
        if not line.strip():
            break
        if line.startswith("!"):
            break
        parts = line.split()
        try:
            pa = int(parts[0])
            float(parts[2])
            float(parts[3])
            int(parts[1])
            int(parts[4])
            parentatom = pa
            nvecs = 0
            continue
        except (ValueError, IndexError):
            pass
        if parentatom is not None:
            nvecs += 1
    assert parentatom is not None
    assert nvecs == coords.count(parentatom)

    details = render_complete_modes(spec)
    assert "Complete modes details" in details
    assert "normfactor" in details
    assert "As=" in details and "Ap=" in details

    topas = render_topas(spec)
    assert topas.startswith("'Topas .str file generated by ISODISTORT")
    assert "space_group" in topas
    assert "prm  !a1" in topas or "prm !a1" in topas
    assert "mode-amplitude to delta" in topas

def test_mode_maxamp_is_inverse_dmax_at_unit_as():
    """ISODISTORT slider/TOPAS bound makes dmax reach exactly 1 Angstrom.

    The eight equal displacements emulate a primitive-cell-normalized mode:
    after As normalization each atom moves 1/sqrt(8) Angstrom, so the allowed
    amplitude is sqrt(8), independently of the arbitrary fractional scale of
    the input eigenvector.
    """
    lattice = np.diag([4.0, 5.0, 10.0])
    mode = np.zeros((8, 3), dtype=float)
    mode[:, 2] = 0.2

    scaled, maxamp = cart_normalized_mode_matrix(mode, lattice)
    dmax_at_as_one = np.max(np.linalg.norm(scaled @ lattice, axis=1))

    assert np.sum((scaled @ lattice) ** 2) == pytest.approx(1.0)
    assert maxamp == pytest.approx(np.sqrt(8.0))
    assert maxamp * dmax_at_as_one == pytest.approx(1.0)

def test_centered_mode_normalization_uses_primitive_cell_and_physical_dmax():
    lattice = np.diag([4.0, 5.0, 10.0])
    mode = np.zeros((8, 3), dtype=float)
    mode[:, 2] = 1.0

    scaled, maxamp = cart_normalized_mode_matrix(
        mode, lattice, centering_mult=2
    )
    cartesian = scaled @ lattice

    assert np.sum(cartesian**2) / 2 == pytest.approx(1.0)
    assert maxamp * np.max(np.linalg.norm(cartesian, axis=1)) == pytest.approx(1.0)

@pytest.mark.parametrize("magnitude", [1e-20, 1e-300, 1e300])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_cart_normalization_preserves_direction_at_extreme_column_scales(magnitude, sign):
    lattice = np.array([[4.0, 0.4, 0.2], [0.1, 5.0, 0.7], [0.3, -0.2, 6.0]])
    mode = np.array([[1.0, -2.0, 0.5], [0.25, 0.0, 1.5], [-0.75, 0.5, 0.0]])
    # The reference is evaluated at ordinary magnitude from the declared
    # row-coordinate metric, independently of the implementation's scaling.
    ordinary_cart = mode @ lattice
    expected = sign * mode / np.sqrt(np.sum(ordinary_cart**2) / 2)
    expected_cart = expected @ lattice
    expected_bound = 1.0 / np.max(np.linalg.norm(expected_cart, axis=1))

    normalized, bound = cart_normalized_mode_matrix(
        sign * magnitude * mode, lattice, centering_mult=np.int64(2)
    )

    np.testing.assert_allclose(normalized, expected, rtol=2e-14, atol=0)
    cart = normalized @ lattice
    assert np.sum(cart**2) / 2 == pytest.approx(1.0, rel=2e-14)
    assert bound == pytest.approx(expected_bound, rel=2e-14)
    assert bound * np.max(np.linalg.norm(cart, axis=1)) == pytest.approx(1.0, rel=2e-14)

@pytest.mark.parametrize("lattice_scale", [1e-300, 1e300])
def test_cart_normalization_avoids_lattice_norm_overflow_and_underflow(lattice_scale):
    mode = np.array([[1.0, 2.0, 0.0], [-1.0, 0.5, 1.0]])
    lattice = lattice_scale * np.array([[4.0, 0.5, 0.0], [0.0, 5.0, 0.25], [0.5, 0.0, 6.0]])
    normalized, bound = cart_normalized_mode_matrix(mode, lattice)
    cart = normalized @ lattice
    assert np.isfinite(normalized).all()
    assert np.sum(cart**2) == pytest.approx(1.0, rel=2e-14)
    assert bound * np.max(np.linalg.norm(cart, axis=1)) == pytest.approx(1.0, rel=2e-14)

@pytest.mark.parametrize("mode,lattice", [
    # A/L loses the first Cartesian contribution while the third is rounded
    # into the subnormal range. Previously the returned norm squared was
    # 0.5325667378703975, yet the returned amplitude bound was exactly 1.
    ([[1.0, 0.0, 3e-24]], np.diag([2e-24, 1e300, 1.0])),
    # All scaled lattice/input components are normal, but their nonzero
    # product underflows; checking only the final dot product misses this.
    ([[1e-24, 0.0, 1.0]], np.diag([1e-200, 1e100, 1.0])),
    # The original column's component range cannot survive division by max|U|.
    ([[1e300, 1e-300, 1e300]], np.diag([1.0, 1e300, 1.0])),
    # Direction survives all intermediate products but not the returned V.
    ([[1.0, 1e-300, 0.0]], 1e300 * np.eye(3)),
])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_cart_normalization_rejects_component_underflow(mode, lattice, sign):
    with pytest.raises(ValueError, match="underflow"):
        cart_normalized_mode_matrix(sign * np.asarray(mode), lattice)

def test_cart_normalization_preserves_small_normal_cartesian_components():
    mode = np.array([[1.0, 0.0, 1e-300]])
    lattice = np.diag([2e-300, 1.0, 1.0])
    normalized, bound = cart_normalized_mode_matrix(mode, lattice)
    expected_cartesian = np.array([[2.0, 0.0, 1.0]]) / np.sqrt(5.0)
    np.testing.assert_allclose(normalized @ lattice, expected_cartesian, rtol=2e-14, atol=0)
    assert bound == pytest.approx(1.0, rel=2e-14)

def test_cart_normalization_ignores_lattice_components_absent_from_the_mode():
    # A subnormal entry in a row multiplied only by structural zeros cannot
    # affect this column. Do not classify it as a lost contribution.
    normalized, bound = cart_normalized_mode_matrix(
        np.array([[1.0, 0.0, 0.0]]), np.diag([1.0, 1e-310, 1.0]),
    )
    np.testing.assert_array_equal(normalized, [[1.0, 0.0, 0.0]])
    assert bound == 1.0

@pytest.mark.parametrize("multiplicity", [0, -1, True, np.bool_(True), 1.5, 2.0, "2", np.nan])
def test_cart_normalization_rejects_invalid_multiplicity(multiplicity):
    with pytest.raises(ValueError, match="positive integer"):
        cart_normalized_mode_matrix(np.array([[1.0, 0.0, 0.0]]), np.eye(3), centering_mult=multiplicity)

@pytest.mark.parametrize("mode", [
    np.zeros((2, 3)), np.zeros((0, 3)), np.ones(3), np.ones((2, 2)),
    np.array([[np.nan, 0.0, 0.0]]), np.array([[np.inf, 0.0, 0.0]]),
    np.array([[1.0 + 1.0j, 0.0, 0.0]]),
])
def test_cart_normalization_rejects_zero_nonfinite_or_invalid_mode(mode):
    with pytest.raises(ValueError):
        cart_normalized_mode_matrix(mode, np.eye(3))

@pytest.mark.parametrize("lattice", [
    np.ones((2, 3)), np.zeros((3, 3)),
    np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 0.0, 1.0]]),
    np.diag([1.0, np.nan, 1.0]), np.diag([1.0, np.inf, 1.0]),
])
def test_cart_normalization_rejects_invalid_lattice(lattice):
    with pytest.raises(ValueError):
        cart_normalized_mode_matrix(np.array([[1.0, 0.0, 0.0]]), lattice)

def test_cart_normalization_does_not_treat_small_nonzero_lattice_volume_as_singular():
    lattice = np.diag([1e-20, 1.0, 1.0])
    normalized, bound = cart_normalized_mode_matrix(np.array([[1.0, 0.0, 0.0]]), lattice)
    np.testing.assert_allclose(normalized @ lattice, [[1.0, 0.0, 0.0]], rtol=2e-14, atol=0)
    assert bound == pytest.approx(1.0, rel=2e-14)

def test_complete_modes_uses_primitive_norm_ap_and_dmax_definitions():
    structure = Structure(
        np.diag([4.0, 5.0, 10.0]),
        ["Eu"] * 4,
        [[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
    )
    subgroup = _sg("GM1+", "P1", number=139, symbol="I4/mmm")
    subgroup.size = 2
    mode = np.zeros((4, 3), dtype=float)
    mode[:, 2] = 1.0
    spec = SubgroupExportSpec(
        subgroup=subgroup,
        structure=structure,
        parent_structure=structure,
        mode_displacements_sc={"mode": mode},
        amplitudes={"mode": 2.0},
    )

    details = render_complete_modes(spec)

    # Σ_conventional |B'u|^2=400 A^2; I centring gives 200 A^2 per
    # primitive cell.  With As=2 A and s=2: norm=1/sqrt(200),
    # Ap=sqrt(2) A and dmax=As*norm*10=sqrt(2) A.
    assert "normfactor = 0.0707107 Angstrom^-1" in details
    assert "As=  2.000000" in details
    assert "Ap=  1.414214" in details
    assert "dmax=  1.414214 Angstrom" in details
    assert "over one primitive unit cell equal to 1.0" in details
    assert "centering-equivalent atoms" in details

def test_i4mmm_cif_keeps_origin_zero_wyckoff():
    """I4/mmm identity subgroup must keep Eu at 2a (0,0,0), not a quarter-shifted 8g."""
    from pymatgen.core import Lattice

    parent = Structure(
        Lattice.tetragonal(4.402, 11.163),
        ["Eu", "Eu", "Al", "Al", "Al", "Al", "Al", "Al", "Al", "Al"],
        [
            [0, 0, 0], [0.5, 0.5, 0.5],
            [0, 0.5, 0.25], [0.5, 0, 0.25], [0.5, 0, 0.75], [0, 0.5, 0.75],
            [0, 0, 0.38], [0, 0, 0.62], [0.5, 0.5, 0.12], [0.5, 0.5, 0.88],
        ],
    )
    sg = _formats_SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        subgroup_index=1,
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin=[0.0, 0.0, 0.0],
        origin_raw="(0,0,0)",
        k_point_label="GM",
        irrep_label="GM1+",
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    spec = SubgroupExportSpec(
        subgroup=sg,
        structure=parent,
        parent_structure=parent,
        parent_sg=139,
        parent_symbol="I4/mmm",
    )
    text = render_cif(parent, spec)
    assert 'H-M "I 4/m 2/m 2/m"' in text
    assert re.search(r"Eu\d*_1\s+Eu\s+2\s+a\s+0\.00000\s+0\.00000\s+0\.00000", text)
    assert re.search(r"Al\d*_1\s+Al\s+4\s+d\s+0\.00000\s+0\.50000\s+0\.25000", text)
    assert re.search(r"Al\d*_\d+\s+Al\s+4\s+e\s+0\.00000\s+0\.00000\s+0\.38000", text)
    assert "_iso_parent-to-child.transform_Pp_abc a,b,c;0,0,0" in text
    assert "8 g" not in text

def test_k_active_substitutes_ld_parameter():
    from backend.tables.kpoints_official import official_special_k_coords

    coords = official_special_k_coords(139, "LD", [], ["1/6"])
    assert coords == ["0", "0", "1/6"]

def test_subgroup_cif_origin_preserves_atom_count(tmp_path):
    """The subgroup's origin must make its symmetry operations valid."""
    from pymatgen.io.cif import CifParser

    from backend.input_cif import read_cif

    parent = read_cif(_formats_REPO / "experiment_data" / "EuAl4 Parent.cif")
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=137, space_group_symbol="P4_2/nmc",
        subgroup_index=2, size=2, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0.25, 0.75, 0.25], origin_raw="(1/4,3/4,1/4)",
        k_point_label="M", irrep_label="M4-", parent_sg=139,
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=parent, parent_structure=parent,
        parent_sg=139, parent_symbol="I4/mmm",
    )
    path = tmp_path / "subgroup.cif"
    path.write_text(render_cif(parent, spec), encoding="utf-8")
    parsed = CifParser(str(path), occupancy_tolerance=100).parse_structures(
        primitive=False
    )[0]
    assert len(parsed) == len(parent) == 10
    assert parsed.composition == parent.composition

def test_explicit_valid_origin_outweighs_asymmetric_site_count():
    """NdNiO2 A1+ P1 uses origin zero even though z=1/4 has fewer ASU rows."""
    from backend.utils.parent_header import parent_wyckoff_display
    from backend.export.isodistort_cif import _apply_origin_choice, _setting
    from backend.input_cif import SymmetryValidator, read_cif
    from backend.input_cif.coordinate_transform import build_supercell

    path = _formats_REPO / "experiment_data" / "NdNiO2 own.cif"
    parent = read_cif(path)
    basis = [[1, 1, 0], [-1, 1, 0], [0, 0, 2]]
    child = build_supercell(parent, basis)
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=139, space_group_symbol="I4/mmm",
        subgroup_index=1, size=2, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=basis, origin=[0, 0, 0], k_point_label="A",
        irrep_label="A1+", parent_sg=123,
    )
    _shifted, origin = _apply_origin_choice(
        child, _setting(subgroup, force_p1=False), [0, 0, 0]
    )
    assert np.allclose(origin, [0, 0, 0])
    wyckoff = SymmetryValidator().validate(parent)["wyckoff_sites"]
    labels = parent_wyckoff_display(parent, wyckoff, path)
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=child, parent_structure=parent,
        parent_sg=123, parent_symbol="P4/mmm",
        parent_wyckoff_sites=wyckoff, parent_wyckoff_lines=labels,
    )
    displayed = render_isoviz(spec).split("!atomcoordlist \n", 1)[1].split("\n\n", 1)[0]
    assert len(displayed.splitlines()) == 35
    assert "0.25000   0.25000   0.00000" in displayed

def test_topas_uses_target_subgroup_orbits():
    """A P1 TOPAS file must keep sites equivalent only in the parent group."""
    parent = Structure(
        [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=1, space_group_symbol="P1",
        subgroup_index=1, size=1, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM1+",
        parent_sg=229,
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=parent, parent_structure=parent,
        parent_sg=229, parent_symbol="Im-3m",
        parent_wyckoff_sites=[{
            "species": "Fe", "multiplicity": 2, "wyckoff_letter": "a",
            "representative_index": 0, "equivalent_indices": [0, 1],
            "display_label": "Fe", "display_order": 0,
        }],
    )
    text = render_topas(spec)
    sites = re.findall(r"^\s*site (\S+)", text, re.M)
    assert len(sites) == len(set(sites)) == 2

def test_isoviz_expands_boundary_images_and_mode_vectors():
    """Boundary atoms and their mode vectors need the same display images."""
    parent = Structure(
        [[3, 0, 0], [0, 3, 0], [0, 0, 3]], ["Fe"], [[0, 0, 0]],
    )
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=221, space_group_symbol="Pm-3m",
        subgroup_index=1, size=1, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM1+",
        parent_sg=221,
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=parent, parent_structure=parent,
        parent_sg=221, parent_symbol="Pm-3m",
        mode_displacements_sc={"test": np.array([[1.0, 0.0, 0.0]])},
    )
    text = render_isoviz(spec)
    atom_rows = text.split("!atomcoordlist \n", 1)[1].split("\n\n", 1)[0].splitlines()
    flags = text.split("!atomsinunitcell \n", 1)[1].split("\n\n", 1)[0].splitlines()
    mode_rows = text.split("!displacivemodelist \n", 1)[1].split("\n\n", 1)[0].splitlines()
    assert len(atom_rows) == len(flags) == 8
    assert [int(flag) for flag in flags] == [1, 0, 0, 0, 0, 0, 0, 0]
    assert len(mode_rows) == 9  # header plus one vector per displayed image

def test_isoviz_subtypes_follow_target_subgroup_orbits():
    """Atoms merged by parent symmetry remain separate P1 subatoms."""
    parent = Structure(
        [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=1, space_group_symbol="P1",
        subgroup_index=1, size=1, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM1+",
        parent_sg=229,
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=parent, parent_structure=parent,
        parent_sg=229, parent_symbol="Im-3m",
        parent_wyckoff_sites=[{
            "species": "Fe", "multiplicity": 2, "wyckoff_letter": "a",
            "representative_index": 0, "equivalent_indices": [0, 1],
            "display_label": "Fe", "display_order": 0,
        }],
    )
    text = render_isoviz(spec)
    subtypes = text.split("!atomsubtypelist \n", 1)[1].split("\n\n", 1)[0].splitlines()
    coord_rows = text.split("!atomcoordlist \n", 1)[1].split("\n\n", 1)[0].splitlines()
    assert len(subtypes) == 2
    assert {int(row.split()[1]) for row in coord_rows} == {1, 2}

def test_exports_preserve_parent_cif_atom_labels_and_order():
    """The NdNiO2 source spells its sites O / ND / NI in that order."""
    from backend.utils.parent_header import parent_wyckoff_display
    from backend.input_cif import SymmetryValidator, read_cif

    path = _formats_REPO / "experiment_data" / "NdNiO2 own.cif"
    parent = read_cif(path)
    sites = SymmetryValidator().validate(parent)["wyckoff_sites"]
    lines = parent_wyckoff_display(parent, sites, path)
    assert [line.split()[0] for line in lines] == ["O", "ND", "NI"]
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=123, space_group_symbol="P4/mmm",
        subgroup_index=1, size=1, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM1+",
        parent_sg=123,
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=parent, parent_structure=parent,
        parent_sg=123, parent_symbol="P4/mmm",
        parent_wyckoff_sites=sites, parent_wyckoff_lines=lines,
    )
    type_rows = render_isoviz(spec).split("!atomtypelist \n", 1)[1].split("\n\n", 1)[0]
    assert [row.split()[1] for row in type_rows.splitlines()] == ["O", "ND", "NI"]
    site_rows = render_cif(parent, spec).split("_atom_site_label", 1)[1]
    assert "O_1" in site_rows and "ND_1" in site_rows and "NI_1" in site_rows

def test_pmmm_wyckoff_letters_use_target_orbits_and_fixed_abc_setting():
    """Accidental parent symmetry must not relabel Pmmm sites after an axis swap."""
    from backend.input_cif import read_cif
    from backend.input_cif.coordinate_transform import build_supercell

    parent = read_cif(_formats_REPO / "experiment_data" / "NdNiO2 own.cif")
    basis = [[-3, 0, 0], [0, 0, 1], [0, 2, 0]]
    child = build_supercell(parent, basis)
    subgroup = _formats_SubgroupInfo(
        index=0, space_group_number=47, space_group_symbol="Pmmm",
        subgroup_index=6, size=6, opd_symbol="P1", opd_dir_raw="(a)",
        basis_vectors=basis, basis_raw="(-3,0,0),(0,0,1),(0,2,0)",
        origin=[0, 0, 0], origin_raw="(0,0,0)",
        k_point_label="Y", irrep_label="Y1",
        k_coordinates=["1/3", "1/2", "0"], parent_sg=123,
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup, structure=child, parent_structure=parent,
        parent_sg=123, parent_symbol="P4/mmm",
    )

    text = render_cif(child, spec)

    # These are standard-setting Pmmm coordinates and Wyckoff symbols.  A
    # metric-driven spglib axis permutation used to emit q/r/u/v instead.
    assert re.search(r"\s1\s+a\s+0\.00000\s+0\.00000\s+0\.00000", text)
    assert re.search(r"\s2\s+i\s+0\.16667\s+0\.00000\s+0\.00000", text)
    assert re.search(r"\s2\s+q\s+0\.00000\s+0\.00000\s+0\.25000", text)
    assert re.search(r"\s4\s+w\s+0\.66667\s+0\.00000\s+0\.25000", text)
    assert re.search(r"\s4\s+x\s+0\.16667\s+0\.50000\s+0\.25000", text)

def test_method1_x4_p3_cif_matches_official_setting(tmp_path):
    """X4- P3 CIF: subgroup Cmma #67, origin-2 ASU, same atoms as official."""
    from pymatgen.core import Lattice
    from pymatgen.io.cif import CifParser

    from tests.validation.backend.paths import ITEM_DIR, TRUE_DIR

    official = TRUE_DIR / "X4- P3 CIF.cif"
    local_sc = ITEM_DIR / "X4- P3 CIF.cif"
    if not official.is_file() or not local_sc.is_file():
        pytest.skip("official/local X4- P3 CIF not present")
    sc = CifParser(str(local_sc), occupancy_tolerance=100).parse_structures(primitive=False)[0]
    ref = CifParser(str(official), occupancy_tolerance=100).parse_structures(primitive=False)[0]
    parent = Structure(
        Lattice.tetragonal(4.402, 11.163),
        ["Eu", "Eu", "Al", "Al", "Al", "Al", "Al", "Al", "Al", "Al"],
        [
            [0, 0, 0], [0.5, 0.5, 0.5],
            [0, 0.5, 0.25], [0.5, 0, 0.25], [0.5, 0, 0.75], [0, 0.5, 0.75],
            [0, 0, 0.38], [0, 0, 0.62], [0.5, 0.5, 0.12], [0.5, 0.5, 0.88],
        ],
    )
    parent_wyckoff = [
        {"species": "Eu", "multiplicity": 2, "wyckoff_letter": "a",
         "representative_index": 0, "equivalent_indices": [0, 1]},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "d",
         "representative_index": 2, "equivalent_indices": [2, 3, 4, 5]},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "e",
         "representative_index": 6, "equivalent_indices": [6, 7, 8, 9]},
    ]
    sg = _formats_SubgroupInfo(
        index=0,
        space_group_number=67,
        space_group_symbol="Cmma",
        subgroup_index=4,
        size=2,
        opd_symbol="P3",
        opd_dir_raw="(a;0)",
        basis_vectors=[[0, 0, 1], [1, -1, 0], [1, 1, 0]],
        basis_raw="(0,0,1),(1,-1,0),(1,1,0)",
        origin=[0.0, 0.5, 0.5],
        origin_raw="(0,1/2,1/2)",
        k_point_label="X",
        irrep_label="X4-",
        k_coordinates=["1/2", "1/2", "0"],
        parent_sg=139,
    )
    spec = SubgroupExportSpec(
        subgroup=sg,
        structure=sc,
        parent_structure=parent,
        parent_sg=139,
        parent_symbol="I4/mmm",
        folder_name="X4- P3",
        parent_wyckoff_sites=parent_wyckoff,
    )
    text = render_cif(sc, spec)
    assert "data_isodistort-output" in text
    assert "_symmetry_Int_Tables_number 67" in text
    assert "C 2/m 2/m 2/a" in text
    assert "067:-C 2a 2" in text
    assert "_iso_parent-to-child.transform_Pp_abc" in text
    assert "c,a-b,a+b" in text
    assert re.search(r"Eu1_1\s+Eu\s+4\s+g\s+0\.00000\s+0\.25000\s+0\.25000", text)
    assert re.search(r"Al1_1\s+Al\s+4\s+a\s+0\.25000\s+0\.00000\s+0\.00000", text)
    assert re.search(r"Al1_2\s+Al\s+4\s+b\s+0\.25000\s+0\.00000\s+0\.50000", text)
    assert re.search(r"Al2_1\s+Al\s+8\s+n\s+0\.38000\s+0\.25000\s+0\.25000", text)
    assert "0,0,Dz" in text
    assert "Dx,0,Dz" in text
    assert "# Eu1 2a (0,0,0)" in text
    assert "# Al1 4d (0,1/2,1/4)" in text
    assert "# Al2 4e (0,0,z), z= 0.38000" in text
    tmp_cif = tmp_path / "x4_p3.cif"
    tmp_cif.write_text(text, encoding="utf-8")
    got = CifParser(str(tmp_cif), occupancy_tolerance=100).parse_structures(primitive=False)[0]
    assert len(got) == len(ref) == 20
    assert abs(got.lattice.a - ref.lattice.a) < 1e-4
    assert abs(got.lattice.b - ref.lattice.b) < 1e-4

    def _match(first, second, tol=1e-4) -> bool:
        used: set[int] = set()
        for site in first:
            found = False
            for j, other in enumerate(second):
                if j in used or site.species_string != other.species_string:
                    continue
                delta = np.abs(site.frac_coords - other.frac_coords)
                delta = np.minimum(delta, 1.0 - delta)
                if np.all(delta < tol):
                    used.add(j)
                    found = True
                    break
            if not found:
                return False
        return True

    assert _match(got, ref)

# Original coverage: test_microscopic_column_provenance.py.

def test_domain_extension_evidence_survives_mapping_and_raw_contract():
    parent, _reference, _canonical, selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    mapper = _microscopic_column_provenance_DistortionMapper()
    mapped_arrays = {
        mode.amplitude_key: np.asarray(
            [[float(index + 1), 0.0, 0.0], [0.0, float(index + 1), 0.0]]
        )
        for index, mode in enumerate(selected)
    }
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: mapped_arrays
    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        [],
        selected,
        np.eye(3).tolist(),
        subgroup_context=_microscopic_column_provenance_subgroup(),
        frame_id="child-frame",
        atom_ids=("atom-0", "atom-1"),
    )

    assert [column.domain_extension_evidence for column in mapped] == [
        mode.microscopic_domain_extension for mode in selected
    ]
    raw = RawModeColumn.from_values(
        mapped[0].mode_id,
        0,
        mapped[0].frame_id,
        mapped[0].atom_ids,
        mapped[0].displacements,
        ModeScaleProvenance(
            source="iso_microscopic_with_verified_bush_domain_extension",
            convention="mapped_unmixed_verified_full_domain_extension",
            direction_resolved=True,
        ),
        mode_identity=mapped[0].mode_identity,
        microscopic_provenance=mapped[0].provenance,
        microscopic_domain_extension=mapped[0].domain_extension_evidence,
    )
    assert raw.source_binding_issue() is None
    assert raw.microscopic_domain_extension == mapped[0].domain_extension_evidence

def test_public_raw_contract_revalidates_complete_extension_groups():
    parent, _reference, _canonical, selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    mapper = _microscopic_column_provenance_DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.asarray(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float,
        )
        for mode in selected
    }
    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        [],
        selected,
        np.eye(3).tolist(),
        subgroup_context=_microscopic_column_provenance_subgroup(),
        frame_id="child-frame",
        atom_ids=("atom-0", "atom-1"),
    )

    def _raw(column, index, *, extension=None, source=None):
        return RawModeColumn.from_values(
            column.mode_id,
            index,
            column.frame_id,
            column.atom_ids,
            column.displacements,
            ModeScaleProvenance(
                source=source or (
                    "iso_microscopic_with_verified_bush_domain_extension"
                ),
                convention="mapped_unmixed_verified_full_domain_extension",
                direction_resolved=True,
            ),
            mode_identity=column.mode_identity,
            microscopic_provenance=column.provenance,
            microscopic_domain_extension=(
                column.domain_extension_evidence
                if extension is None
                else extension
            ),
            label=f"mode-{index}",
        )

    raw = [_raw(column, index) for index, column in enumerate(mapped)]
    operation = ExactSeitzOperation.from_values(np.eye(3), (0, 0, 0))
    atoms = (
        ChildAtom.from_values("atom-0", (0, 0, 0), "Fe"),
        ChildAtom.from_values("atom-1", (Fraction(1, 2), 0, 0), "Fe"),
    )
    with pytest.raises(
        DisplaciveModelError,
        match=r"invalid_domain_extension_evidence.*source basis is incomplete",
    ):
        project_modes_to_free_rows(
            (),
            atoms,
            (operation,),
            [replace(raw[0], column_index=0)],
            frame_id="child-frame",
            reference_lattice=np.eye(3),
        )

    first_extension = raw[0].microscopic_domain_extension
    assert first_extension is not None
    tampered = replace(
        raw[0],
        microscopic_domain_extension=replace(
            first_extension,
            coefficient_column=(0.25, 0.5),
        ),
    )
    with pytest.raises(
        DisplaciveModelError,
        match=r"invalid_domain_extension_evidence.*change-of-basis digest",
    ):
        project_modes_to_free_rows(
            (),
            atoms,
            (operation,),
            [tampered, raw[1]],
            frame_id="child-frame",
            reference_lattice=np.eye(3),
        )

    assert replace(
        raw[0],
        microscopic_domain_extension=None,
    ).source_binding_issue() == "microscopic scale source claims missing domain extension"
    assert replace(
        raw[0],
        scale_provenance=ModeScaleProvenance(
            source="iso_microscopic_display_distortion",
            convention="mapped_unmixed_fractional_source_column",
            direction_resolved=True,
        ),
    ).source_binding_issue() == "microscopic scale source omits verified domain extension"

# Original coverage: test_production_displacive_core.py.

def _matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )

def _add(left, right):
    return tuple(left[index] + right[index] for index in range(3))

def _subtract(left, right):
    return tuple(left[index] - right[index] for index in range(3))

def test_exact_child_operations_conjugate_nondiagonal_basis_and_nonzero_origin():
    hall_number, _symbol = hall_for_number(3)
    database = spglib.get_symmetry_from_database(hall_number)
    basis = rational_matrix(((0, 1, 0), (1, 1, 0), (0, 0, -1)))
    origin = (Fraction(1, 3), Fraction(1, 5), Fraction(1, 7))
    parent_from_child = transpose(basis)
    child_from_parent = inverse(parent_from_child)
    parent_operations = []
    expected = []
    for rotation, translation in zip(
        database["rotations"],
        database["translations"],
        strict=True,
    ):
        child_rotation = rational_matrix(rotation.tolist())
        child_translation = tuple(as_fraction(float(value)) for value in translation)
        parent_rotation = multiply(
            multiply(parent_from_child, child_rotation),
            child_from_parent,
        )
        parent_translation = _add(
            _matrix_vector(parent_from_child, child_translation),
            _subtract(origin, _matrix_vector(parent_rotation, origin)),
        )
        parent_operations.append(AffineOperation(parent_rotation, parent_translation))
        expected.append(
            ExactSeitzOperation.from_values(child_rotation, child_translation)
        )
    subgroup = SubgroupInfo(
        index=0,
        parent_sg=3,
        space_group_number=3,
        space_group_symbol="P2",
        subgroup_index=1,
        size=1,
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM1+",
        opd_symbol="P1",
        basis_vectors=[list(row) for row in basis],
        origin=list(origin),
    )
    affine_child = SimpleNamespace(
        hall_number=hall_number,
        operations=tuple(parent_operations),
    )

    actual = exact_child_operations(subgroup, affine_child)

    assert {operation.key for operation in actual} == {
        operation.key for operation in expected
    }

# Original coverage: test_strain_export.py.

_strain_export_REPO = TESTS_ROOT.parent.parent

_F02 = _strain_export_REPO / "output_compare" / "EuAl4 Parent.cif" / "官网" / "Method4" / "F02"

def _p1_export_spec(amplitudes: tuple[float, ...] | None = None) -> SubgroupExportSpec:
    structure = Structure(Lattice.from_parameters(4.0, 5.0, 6.0, 78.0, 83.0, 71.0), ["H"], [[0.13, 0.27, 0.39]])
    fixed = compute_homogeneous_strain_modes(
        [np.eye(3)], structure.lattice.matrix
    )
    definitions = tuple(
        CanonicalStrainModeDefinition(
            label=f"GM_TESTstrain_{index + 1}(a)",
            q_raw=np.eye(6)[index],
            irrep_label="GM_TEST",
            irrep_direction="(a)",
        )
        for index in range(6)
    )
    result = apply_canonical_strain_basis(fixed, definitions)
    values = amplitudes if amplitudes is not None else (0.0,) * len(result.modes)
    return SubgroupExportSpec(
        subgroup=_strain_export_subgroup(1, "P1"),
        structure=structure,
        parent_structure=structure,
        parent_sg=1,
        parent_symbol="P1",
        distortion_types=["strain"],
        strain_data=StrainExportData(result=result, amplitudes=values),
    )

def _read_cif_strain_matrix(text: str, mode_count: int) -> np.ndarray:
    block = text.split("_iso_strainmodematrix_value", 1)[1].split(
        "_iso_parentcell_length_a", 1
    )[0]
    matrix = np.zeros((6, mode_count), dtype=float)
    for row, column, value in re.findall(
        r"^\s+(\d+)\s+(\d+)\s+([-+0-9.]+)\s*$", block, re.MULTILINE
    ):
        matrix[int(row) - 1, int(column) - 1] = float(value)
    return matrix

def _read_isoviz_strain_vectors(text: str) -> np.ndarray:
    block = text.split("!strainmodelist", 1)[1].split(
        "#parentatom/dispmodenum", 1
    )[0]
    rows = [
        [float(value) for value in line.split()]
        for line in block.splitlines()
        if len(line.split()) == 6
    ]
    return np.asarray(rows, dtype=float)

def _read_complete_strain_vectors(text: str, key: str) -> np.ndarray:
    rows = re.findall(rf"^  {re.escape(key)}\s*= \(([^)]*)\)$", text, re.MULTILINE)
    return np.asarray(
        [[float(value) for value in row.split()] for row in rows], dtype=float
    )

def _read_topas_fixed_cell(text: str) -> np.ndarray:
    values = {
        keyword: float(value)
        for keyword, value in re.findall(
            r"^\s*(a|b|c|al|be|ga)\s+([-+0-9.]+)\s*$", text, re.MULTILINE
        )
    }
    assert set(values) == {"a", "b", "c", "al", "be", "ga"}
    return np.asarray([values[key] for key in ("a", "b", "c", "al", "be", "ga")])

def _read_cif_lattice(text: str) -> Lattice:
    values = {
        key: float(value)
        for key, value in re.findall(
            r"^_cell_(length_[abc]|angle_(?:alpha|beta|gamma))\s+"
            r"([-+0-9.]+)\s*$",
            text,
            re.MULTILINE,
        )
    }
    assert set(values) == {
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    }
    return Lattice.from_parameters(
        values["length_a"],
        values["length_b"],
        values["length_c"],
        values["angle_alpha"],
        values["angle_beta"],
        values["angle_gamma"],
    )

def _read_complete_lattice(text: str, title: str) -> Lattice:
    block = text.split(title, 1)[1]
    match = re.search(
        r"a=([-+0-9.]+),?\s+b=([-+0-9.]+),?\s+c=([-+0-9.]+),?\s+"
        r"alpha=([-+0-9.]+),?\s+beta=([-+0-9.]+),?\s+gamma=([-+0-9.]+)",
        block,
    )
    assert match is not None
    return Lattice.from_parameters(*(float(value) for value in match.groups()))

def _read_isoviz_parentcell(text: str) -> np.ndarray:
    block = text.split("!parentcell", 1)[1]
    values = next(
        line.split()
        for line in block.splitlines()
        if len(line.split()) == 6
    )
    return np.asarray([float(value) for value in values])

def _read_isoviz_parentbasis(text: str) -> np.ndarray:
    block = text.split("!parentbasis", 1)[1].split(
        "!conv2primparentbasis", 1
    )[0]
    rows = [
        [float(value) for value in line.split()]
        for line in block.splitlines()
        if len(line.split()) == 3
    ]
    return np.asarray(rows, dtype=float)

def _read_isoviz_atom_coords(text: str) -> np.ndarray:
    block = text.split("!atomcoordlist", 1)[1].split("!atomsinunitcell", 1)[0]
    rows = [
        [float(value) for value in line.split()[-3:]]
        for line in block.splitlines()
        if len(line.split()) == 6
    ]
    return np.asarray(rows, dtype=float)

def _read_cif_strain_rows(
    text: str,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    mode_block = text.split("_iso_strainmode_value", 1)[1].split("loop_", 1)[0]
    mode_rows = re.findall(
        r"^\s*\d+\s+(\S*strain\S*)\s+([-+0-9.]+)\s*$",
        mode_block,
        re.MULTILINE,
    )
    norm_block = text.split("_iso_strainmodenorm_value", 1)[1].split(
        "loop_", 1
    )[0]
    norms = np.asarray(
        [
            float(value)
            for value in re.findall(
                r"^\s*\d+\s+([-+0-9.]+)\s*$", norm_block, re.MULTILINE
            )
        ]
    )
    value_block = text.split("_iso_strain_value", 1)[1].split(
        "# matrix conversion", 1
    )[0]
    raw_sum = np.asarray(
        [
            float(value)
            for value in re.findall(
                r"^\s*\d+\s+E_\d+\s+([-+0-9.]+)\s*$",
                value_block,
                re.MULTILINE,
            )
        ]
    )
    return (
        [label for label, _value in mode_rows],
        np.asarray([float(value) for _label, value in mode_rows]),
        norms,
        raw_sum,
    )

def _read_isoviz_strain_rows(text: str) -> tuple[list[str], np.ndarray, np.ndarray]:
    block = text.split("!strainmodelist", 1)[1].split(
        "!displacivemodelist", 1
    )[0]
    rows = re.findall(
        r"^\s*\d+\s+([-+0-9.]+)\s+[-+0-9.]+\s+\d+\s+(\S+)\s*$"
        r"\n\s*((?:[-+0-9.]+\s+){5}[-+0-9.]+)\s*$",
        block,
        re.MULTILINE,
    )
    return (
        [label for _amplitude, label, _vector in rows],
        np.asarray([float(amplitude) for amplitude, _label, _vector in rows]),
        np.asarray(
            [
                [float(value) for value in vector.split()]
                for _amplitude, _label, vector in rows
            ]
        ),
    )

def _read_official_complete_strain_rows(
    text: str,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    definition_block = text.split(
        "<b>Parent-cell strain mode definitions</b>", 1
    )[1].split("<b>Parent-cell strain mode amplitudes</b>", 1)[0]
    definitions = re.findall(
        r"^(\S*strain\S*)\s+normfactor\s*=\s*([-+0-9.]+)\s*$"
        r"\n\s*((?:[-+0-9.]+\s+){5}[-+0-9.]+)\s*$",
        definition_block,
        re.MULTILINE,
    )
    amplitude_block = text.split(
        "<b>Parent-cell strain mode amplitudes</b>", 1
    )[1]
    amplitudes = np.asarray(
        [
            float(value)
            for value in re.findall(
                r"^\S*strain\S*\s+([-+0-9.]+)\s*$",
                amplitude_block,
                re.MULTILINE,
            )
        ]
    )
    return (
        [label for label, _norm, _vector in definitions],
        np.asarray([float(norm) for _label, norm, _vector in definitions]),
        np.asarray(
            [
                [float(value) for value in vector.split()]
                for _label, _norm, vector in definitions
            ]
        ),
        amplitudes,
    )

def test_4310_n1_plus_4d1_has_six_embedded_strain_modes_in_strain_exports() -> None:
    parent = read_cif(_strain_export_REPO / "experiment_data" / "4310_tetra.cif")
    subgroup = SubgroupInfo(
        index=0,
        parent_sg=139,
        space_group_number=2,
        space_group_symbol="P-1",
        subgroup_index=64,
        size=8,
        opd_symbol="4D1",
        opd_dir_raw="(a;b;c;d)",
        basis_vectors=[[2, 0, 0], [0, 2, 0], [-1, -1, 1]],
        origin=[0.0, 0.0, 0.0],
        k_point_label="N",
        k_coordinates=["1/2", "0", "1/2"],
        irrep_label="N1+",
    )
    api = _bare_api(parent, ["strain"])
    api._iso = _TetragonalMacroIso()

    data = api._strain_export_data_for_subgroup(subgroup)

    assert data is not None
    assert len(data.result.modes) == 6
    assert data.result.diagnostics.point_group_order == 2

    spec = SubgroupExportSpec(
        subgroup=subgroup,
        structure=api._supercell_for_subgroup(subgroup),
        parent_structure=parent,
        parent_sg=139,
        parent_symbol="I4/mmm",
        distortion_types=["strain"],
        strain_data=data,
    )
    cif = render_cif(spec.structure, spec)
    isoviz = render_isoviz(spec)
    complete = render_complete_modes(spec)
    topas = render_topas(spec)
    expected_raw = np.column_stack(
        [np.asarray(mode.q_raw) for mode in data.result.modes]
    )
    expected_unit = np.vstack(
        [np.asarray(mode.q_unit) for mode in data.result.modes]
    )

    assert "_iso_strainmode_number    6" in cif
    assert all(
        label in cif
        for label in (
            "I4/mmm[0,0,0]GM1+(a)strain_1(a)",
            "I4/mmm[0,0,0]GM1+(a)strain_2(a)",
            "I4/mmm[0,0,0]GM2+(a)strain(a)",
            "I4/mmm[0,0,0]GM4+(a)strain(a)",
            "I4/mmm[0,0,0]GM5+(a,b)strain(a)",
            "I4/mmm[0,0,0]GM5+(a,b)strain(b)",
        )
    )
    assert _read_cif_strain_matrix(cif, 6) == pytest.approx(
        expected_raw, abs=5.1e-6
    )
    assert len(
        re.findall(
            r"^\s+\d+\s+[-+0-9.]+\s+0\.10000\s+\d+\s+\S*strain",
            isoviz,
            re.MULTILINE,
        )
    ) == 6
    assert _read_isoviz_strain_vectors(isoviz) == pytest.approx(
        expected_unit, abs=5.1e-6
    )
    assert complete.count("label_status=canonical_iso_macro") == 6
    assert _read_complete_strain_vectors(complete, "q_raw ") == pytest.approx(
        expected_raw.T, abs=5.1e-9
    )
    assert _read_complete_strain_vectors(complete, "q_unit") == pytest.approx(
        expected_unit, abs=5.1e-9
    )
    assert not re.search(r"prm\s+!s\d+\s", topas)
    assert "strain_q" not in topas
    assert "homogeneous strain modes" not in topas
    assert _read_topas_fixed_cell(topas) == pytest.approx(
        spec.structure.lattice.parameters, abs=5.1e-5
    )

def test_eual4_f02_cif_isoviz_complete_modes_close_same_strain_contract() -> None:
    official_cif = (_F02 / "subgroup.cif").read_text(encoding="utf-8")
    official_isoviz = (_F02 / "data.isoviz").read_text(encoding="utf-8")
    official_complete = (
        _F02 / "ISODISTORT_ complete modes details.html"
    ).read_text(encoding="utf-8")
    cif_labels, cif_amplitudes, cif_norms, cif_raw_sum = _read_cif_strain_rows(
        official_cif
    )
    isoviz_labels, isoviz_amplitudes, isoviz_q_unit = (
        _read_isoviz_strain_rows(official_isoviz)
    )
    complete_labels, complete_norms, complete_q_raw, complete_amplitudes = (
        _read_official_complete_strain_rows(official_complete)
    )
    official_q_raw = _read_cif_strain_matrix(official_cif, 6).T

    assert cif_labels == [
        "I4/mmm[0,0,0]GM1+(a)strain_1(a)",
        "I4/mmm[0,0,0]GM1+(a)strain_2(a)",
        "I4/mmm[0,0,0]GM2+(a)strain(a)",
        "I4/mmm[0,0,0]GM4+(a)strain(a)",
        "I4/mmm[0,0,0]GM5+(a,b)strain(a)",
        "I4/mmm[0,0,0]GM5+(a,b)strain(b)",
    ]
    assert complete_labels == cif_labels
    assert isoviz_labels == [
        "GM1+strain_1(a)",
        "GM1+strain_2(a)",
        "GM2+strain(a)",
        "GM4+strain(a)",
        "GM5+strain(a)",
        "GM5+strain(b)",
    ]
    assert isoviz_amplitudes == pytest.approx(cif_amplitudes, abs=5.1e-6)
    assert complete_amplitudes == pytest.approx(cif_amplitudes, abs=5.1e-6)
    assert complete_norms == pytest.approx(cif_norms, abs=5.1e-6)
    assert complete_q_raw == pytest.approx(official_q_raw, abs=5.1e-6)
    assert isoviz_q_unit == pytest.approx(
        official_q_raw * cif_norms[:, None], abs=5.1e-5
    )

    parent = read_cif(_strain_export_REPO / "experiment_data" / "EuAl4 Parent.cif")
    subgroup = SubgroupInfo(
        index=0,
        parent_sg=139,
        space_group_number=1,
        space_group_symbol="P1",
        subgroup_index=32,
        size=2,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_vectors=[[0, 1, 1], [1, 0, 0], [0, 0, -1]],
        origin=[0.0, 0.0, 0.0],
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM1+",
    )
    api = _bare_api(parent, ["strain"])
    api._iso = _TetragonalMacroIso()
    zero_data = api._strain_export_data_for_subgroup(subgroup)
    assert zero_data is not None
    strain_data = StrainExportData(
        result=zero_data.result,
        amplitudes=tuple(float(value) for value in cif_amplitudes),
    )
    spec = SubgroupExportSpec(
        subgroup=subgroup,
        structure=api._supercell_for_subgroup(subgroup),
        parent_structure=parent,
        parent_sg=139,
        parent_symbol="I4/mmm",
        distortion_types=["strain"],
        strain_data=strain_data,
    )
    generated_cif = render_cif(spec.structure, spec)
    generated_isoviz = render_isoviz(spec)
    generated_complete = render_complete_modes(spec)
    generated_topas = render_topas(spec)
    generated_labels, generated_amplitudes, generated_norms, generated_raw_sum = (
        _read_cif_strain_rows(generated_cif)
    )
    generated_isoviz_labels, generated_isoviz_amplitudes, generated_q_unit = (
        _read_isoviz_strain_rows(generated_isoviz)
    )

    assert generated_labels == cif_labels
    assert generated_amplitudes == pytest.approx(cif_amplitudes, abs=5.1e-6)
    assert generated_norms == pytest.approx(cif_norms, abs=5.1e-5)
    assert generated_raw_sum == pytest.approx(cif_raw_sum, abs=5.1e-6)
    assert _read_cif_strain_matrix(generated_cif, 6).T == pytest.approx(
        official_q_raw, abs=5.1e-6
    )
    assert generated_isoviz_labels == isoviz_labels
    assert generated_isoviz_amplitudes == pytest.approx(
        isoviz_amplitudes, abs=5.1e-6
    )
    assert generated_q_unit == pytest.approx(isoviz_q_unit, abs=5.1e-5)
    # Official IsoVIZ represents strain parametrically: the physical parent
    # cell, inverse reference basis and reference atom coordinates stay fixed;
    # amplitudes times q_unit reconstruct M interactively.
    assert _read_isoviz_parentcell(generated_isoviz) == pytest.approx(
        _read_isoviz_parentcell(official_isoviz), abs=5.1e-5
    )
    assert _read_isoviz_parentbasis(generated_isoviz) == pytest.approx(
        _read_isoviz_parentbasis(official_isoviz), abs=5.1e-5
    )
    generated_atom_coords = _read_isoviz_atom_coords(generated_isoviz)
    official_atom_coords = _read_isoviz_atom_coords(official_isoviz)
    generated_order = np.lexsort(
        (
            generated_atom_coords[:, 2],
            generated_atom_coords[:, 1],
            generated_atom_coords[:, 0],
        )
    )
    official_order = np.lexsort(
        (
            official_atom_coords[:, 2],
            official_atom_coords[:, 1],
            official_atom_coords[:, 0],
        )
    )
    assert generated_atom_coords[generated_order] == pytest.approx(
        official_atom_coords[official_order], abs=5.1e-5
    )
    generated_cif_lattice = _read_cif_lattice(generated_cif)
    official_cif_lattice = _read_cif_lattice(official_cif)
    assert generated_cif_lattice.parameters == pytest.approx(
        official_cif_lattice.parameters, abs=5.1e-5
    )
    official_metric = official_cif_lattice.matrix @ official_cif_lattice.matrix.T
    assert (
        np.linalg.norm(
            generated_cif_lattice.matrix @ generated_cif_lattice.matrix.T
            - official_metric
        )
        / np.linalg.norm(official_metric)
    ) < 1.0e-5
    generated_reference_lattice = _read_complete_lattice(
        generated_complete, "Undistorted superstructure"
    )
    official_reference_lattice = _read_complete_lattice(
        official_complete, "Undistorted superstructure"
    )
    generated_final_lattice = _read_complete_lattice(
        generated_complete, "Distorted superstructure"
    )
    official_final_lattice = _read_complete_lattice(
        official_complete, "Distorted superstructure"
    )
    assert generated_reference_lattice.parameters == pytest.approx(
        official_reference_lattice.parameters, abs=5.1e-5
    )
    assert generated_final_lattice.parameters == pytest.approx(
        official_final_lattice.parameters, abs=5.1e-5
    )
    assert _read_complete_strain_vectors(
        generated_complete, "q_raw "
    ) == pytest.approx(complete_q_raw, abs=5.1e-6)
    assert _read_complete_strain_vectors(
        generated_complete, "q_unit"
    ) == pytest.approx(isoviz_q_unit, abs=5.1e-5)
    assert generated_complete.count("  amplitude =") == 6
    official_topas = (_F02 / "topas.str").read_text(encoding="utf-8")
    assert "strain_q" not in official_topas
    assert "strain_q" not in generated_topas
    assert not re.search(r"prm\s+!s\d+\s", generated_topas)
    assert _read_topas_fixed_cell(generated_topas) == pytest.approx(
        _read_topas_fixed_cell(official_topas), abs=5.1e-5
    )

@pytest.mark.parametrize("distortion_types", ([], ["displacive"]))
def test_strain_is_omitted_when_type_is_not_selected(
    distortion_types: list[str],
) -> None:
    structure = Structure(Lattice.cubic(4.0), ["H"], [[0.0, 0.0, 0.0]])
    api = _bare_api(structure, distortion_types)

    assert api._strain_export_data_for_subgroup(_strain_export_subgroup(221, "Pm-3m")) is None

def test_contract_rejects_strain_data_without_strain_type() -> None:
    valid = _p1_export_spec()
    with pytest.raises(ValueError, match="distortion_types includes 'strain'"):
        SubgroupExportSpec(
            subgroup=valid.subgroup,
            structure=valid.structure,
            distortion_types=["displacive"],
            strain_data=valid.strain_data,
        )

def test_three_strain_exports_share_contract_and_topas_uses_fixed_cell() -> None:
    amplitudes = (0.011, -0.012, 0.013, -0.014, 0.015, -0.016)
    spec = _p1_export_spec(amplitudes)
    assert spec.strain_data is not None
    modes = spec.strain_data.result.modes

    cif = render_cif(spec.structure, spec)
    assert "_iso_strainmode_number    6" in cif
    expected_raw = np.column_stack([np.asarray(mode.q_raw) for mode in modes])
    assert _read_cif_strain_matrix(cif, 6) == pytest.approx(
        expected_raw, abs=5.1e-6
    )
    for mode in modes:
        assert f"   {mode.index}  {mode.normfactor:.5f}" in cif

    isoviz = render_isoviz(spec)
    strain_block = isoviz.split("!strainmodelist", 1)[1].split(
        "#parentatom/dispmodenum", 1
    )[0]
    expected_unit = np.vstack([np.asarray(mode.q_unit) for mode in modes])
    assert _read_isoviz_strain_vectors(isoviz) == pytest.approx(
        expected_unit, abs=5.1e-6
    )
    assert strain_block.count("0.10000") == 6
    assert "GM_TEST" in isoviz

    complete = render_complete_modes(spec)
    assert "Parent-lattice-basis engineering order q" in complete
    assert "q_unit = normfactor * q_raw" in complete
    assert complete.count("label_status=canonical_iso_macro") == 6
    assert complete.count("  normfactor =") == 6
    assert complete.count("  q_raw  = (") == 6
    assert complete.count("  q_unit = (") == 6
    assert complete.count("  amplitude =") == 6
    assert "(no displacive modes available for this subgroup)" in complete
    assert _read_complete_strain_vectors(complete, "q_raw ") == pytest.approx(
        expected_raw.T, abs=5.1e-9
    )
    assert _read_complete_strain_vectors(complete, "q_unit") == pytest.approx(
        expected_unit, abs=5.1e-9
    )

    topas = render_topas(spec)
    assert not re.search(r"prm\s+!s\d+\s", topas)
    assert "strain_q" not in topas
    assert "homogeneous strain modes" not in topas
    expected_lattice = Lattice(
        (np.eye(3) + engineering_voigt_to_tensor(spec.strain_data.applied_engineering_q()))
        @ spec.parent_structure.lattice.matrix
    )
    assert _read_cif_lattice(cif).parameters == pytest.approx(
        expected_lattice.parameters, abs=5.1e-5
    )
    assert _read_topas_fixed_cell(topas) == pytest.approx(
        expected_lattice.parameters, abs=5.1e-5
    )

def test_zero_strain_cif_preserves_reference_child_lattice() -> None:
    spec = _p1_export_spec()
    cif = render_cif(spec.structure, spec)

    assert _read_cif_lattice(cif).parameters == pytest.approx(
        spec.structure.lattice.parameters, abs=5.1e-5
    )

def test_nonzero_strain_cif_preserves_reference_fractional_coordinates() -> None:
    spec = _p1_export_spec((0.01, -0.02, 0.03, -0.04, 0.05, -0.06))

    parsed = CifParser.from_str(render_cif(spec.structure, spec)).parse_structures(
        primitive=False
    )[0]

    assert parsed.frac_coords == pytest.approx(spec.structure.frac_coords, abs=5.1e-6)

def test_consistent_explicit_final_cif_uses_its_fractional_coordinates() -> None:
    spec = _p1_export_spec((0.01, 0.0, 0.0, 0.0, 0.0, 0.0))
    assert spec.strain_data is not None
    final_lattice = Lattice(
        (
            np.eye(3)
            + engineering_voigt_to_tensor(spec.strain_data.applied_engineering_q())
        )
        @ spec.parent_structure.lattice.matrix
    )
    spec.cif_structure = Structure(final_lattice, ["H"], [[0.21, 0.32, 0.43]])

    parsed = CifParser.from_str(
        render_cif(spec.cif_structure, spec)
    ).parse_structures(primitive=False)[0]

    assert parsed.frac_coords == pytest.approx(
        spec.cif_structure.frac_coords, abs=5.1e-6
    )

@pytest.mark.parametrize(
    "renderer", (render_cif, render_isoviz, render_complete_modes, render_topas)
)
def test_strain_export_rejects_reference_child_lattice_mismatch(renderer) -> None:
    spec = _p1_export_spec((0.01, 0.0, 0.0, 0.0, 0.0, 0.0))
    wrong_lattice = Lattice(1.02 * spec.structure.lattice.matrix)
    spec.structure = Structure(
        wrong_lattice,
        [site.species for site in spec.structure],
        spec.structure.frac_coords,
    )

    with pytest.raises(ValueError, match=r"reference child lattice.*B @ P"):
        if renderer is render_cif:
            renderer(spec.structure, spec)
        else:
            renderer(spec)

@pytest.mark.parametrize(
    "renderer", (render_cif, render_isoviz, render_complete_modes, render_topas)
)
def test_strain_export_rejects_inconsistent_explicit_final_structure(renderer) -> None:
    spec = _p1_export_spec((0.01, 0.0, 0.0, 0.0, 0.0, 0.0))
    spec.cif_structure = spec.structure.copy()

    with pytest.raises(ValueError, match=r"explicit final CIF structure.*B @ M @ P"):
        if renderer is render_cif:
            renderer(spec.cif_structure, spec)
        else:
            renderer(spec)

def test_strain_contract_preserves_core_raw_normalized_relation() -> None:
    spec = _p1_export_spec()
    assert spec.strain_data is not None
    for mode in spec.strain_data.result.modes:
        assert mode.q_unit == pytest.approx(mode.normfactor * mode.q_raw, abs=1.0e-14)
        assert np.max(np.abs(mode.q_raw)) == pytest.approx(1.0, abs=1.0e-14)
        expected = 1.0 / math.sqrt(
            float(
                mode.q_raw
                @ np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])
                @ mode.q_raw
            )
        )
        assert mode.normfactor == pytest.approx(expected, abs=1.0e-14)

def _combined_displacive_strain_spec() -> SubgroupExportSpec:
    """Build a real verified displacive contract plus a nonzero P1 strain."""


    data = _export_data(final_x=0.35)
    base = _spec(data)
    parent = data.snapshot_parent()
    fixed = compute_homogeneous_strain_modes(
        [np.eye(3)],
        parent.lattice.matrix,
    )
    definitions = tuple(
        CanonicalStrainModeDefinition(
            label=f"GM1+strain_{index + 1}(a)",
            q_raw=np.eye(6)[index],
            irrep_label="GM1+",
            irrep_direction="(a)",
        )
        for index in range(6)
    )
    strain = StrainExportData(
        result=apply_canonical_strain_basis(fixed, definitions),
        amplitudes=(0.02, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
    return SubgroupExportSpec(
        subgroup=base.subgroup,
        distortion_types=["strain", "displacive"],
        strain_data=strain,
        displacive_data=data,
        require_verified_displacive_data=True,
    )

def test_nonzero_strain_combines_with_verified_displacive_final() -> None:
    spec = _combined_displacive_strain_spec()
    assert spec.strain_data is not None
    combined = resolve_strained_export_structure(spec)
    expected_lattice = (
        np.eye(3)
        + engineering_voigt_to_tensor(spec.strain_data.applied_engineering_q())
    ) @ spec.parent_structure.lattice.matrix

    assert combined.lattice.matrix == pytest.approx(expected_lattice, abs=1.0e-12)
    assert combined.frac_coords == pytest.approx(
        spec.displacive_data.snapshot_final().frac_coords,
        abs=1.0e-14,
    )

    parsed = CifParser.from_str(render_cif(spec.cif_structure, spec)).parse_structures(
        primitive=False
    )[0]
    assert parsed.frac_coords == pytest.approx(combined.frac_coords, abs=5.1e-6)
    assert _read_cif_lattice(render_cif(spec.cif_structure, spec)).parameters == (
        pytest.approx(combined.lattice.parameters, abs=5.1e-5)
    )
    assert "Distorted superstructure:" in render_complete_modes(spec)
    assert _read_topas_fixed_cell(render_topas(spec)) == pytest.approx(
        combined.lattice.parameters,
        abs=5.1e-5,
    )
    assert "!strainmode" in render_isoviz(spec)

@pytest.mark.parametrize(
    "renderer", (render_cif, render_isoviz, render_complete_modes, render_topas)
)
def test_combined_export_rejects_tampered_displacive_final_lattice(renderer) -> None:
    spec = _combined_displacive_strain_spec()
    data = spec.displacive_data
    assert data is not None
    verified_final = data.snapshot_final()
    tampered = Structure(
        Lattice(1.03 * verified_final.lattice.matrix),
        [site.species for site in verified_final],
        verified_final.frac_coords,
        labels=[site.label for site in verified_final],
    )
    object.__setattr__(data, "final_structure", tampered)

    with pytest.raises(
        ValueError,
        match=r"verified displacive final lattice.*B @ P reference frame",
    ):
        if renderer is render_cif:
            renderer(spec.cif_structure, spec)
        else:
            renderer(spec)

# Original coverage: test_structure.py.

def test_export_poscar_reread_same_sg(tmp_path):
    """导出 POSCAR -> 重读 -> 空间群一致（导出/导入闭环）。"""
    lattice = [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]]
    base = Structure(lattice, ["Fe"], [[0, 0, 0]])
    exporter = StructureExporter(tmp_path)
    out = exporter.to_poscar(base, "fe_export")
    s = read_structure(str(out))
    assert len(s) == 1
    sg = SpacegroupAnalyzer(s, symprec=1e-3).get_space_group_number()
    assert sg == 221  # 简单立方 Pm-3m

def test_export_filename_cannot_escape_output_directory(tmp_path):
    """分数型 k 标签含路径分隔符时仍必须导出到指定目录。"""
    base = Structure([[4, 0, 0], [0, 4, 0], [0, 0, 4]], ["Fe"], [[0, 0, 0]])
    exporter = StructureExporter(tmp_path)
    out = exporter.to_cif(base, r"distorted_R4+[1/2,1/2,1/2]\\mode")

    assert out.parent == tmp_path
    assert out.exists()
    assert "/" not in out.name and "\\" not in out.name
    assert out.name == "distorted_R4+[12,12,12]mode.cif"

# Original coverage: test_theory_contracts.py.

def test_unknown_wyckoff_letter_is_not_guessed_from_multiplicity():
    # P-1 has eight distinct multiplicity-one positions: multiplicity cannot
    # determine their letters. Unknown metadata must stay unknown in CIF.
    assert _wyckoff_letter(2, np.array([0.5, 0, 0]), 1,
                           SpaceGroup.from_int_number(2).symmetry_ops) == "?"

def test_sg139_special_c_position_is_not_mislabeled_as_free_d():
    ops = SpaceGroup.from_int_number(139).symmetry_ops
    assert _wyckoff_letter(139, np.array([0, 0.5, 0]), 4, ops) == "c"

@pytest.mark.parametrize("sg,letter,multiplicity,seed", [
    (sg, letter, mult, seed)
    for sg in (67, 69, 71, 99, 107, 123, 129, 139)
    for letter, mult, seed in cif_writer._wyckoff_table(sg)
])
def test_partial_wyckoff_table_agrees_with_database_in_actual_hall(sg, letter, multiplicity, seed):
    import spglib

    hall, _ = hall_for_number(sg)
    operations = spglib.get_symmetry_from_database(hall)
    coordinates, types = [], []
    for atom_type, point in ((1, seed), (2, [.123, .237, .361])):
        orbit = sorted({tuple(np.round((rotation @ point + translation) % 1, 10))
                        for rotation, translation in zip(operations['rotations'], operations['translations'])})
        coordinates.extend(orbit)
        types.extend([atom_type] * len(orbit))
        if atom_type == 1:
            assert len(orbit) == multiplicity
    # The independent generic orbit removes accidental supergroup symmetry.
    # Verify the setting rather than assuming that this construction suffices.
    ds = spglib.get_symmetry_dataset(
        (cif_writer._canonical_metric_lattice(sg), coordinates, types), hall_number=hall)
    assert ds is not None and ds.number == sg and ds.hall_number == hall
    assert np.allclose(ds.transformation_matrix, np.eye(3), atol=1e-8, rtol=0)
    assert np.allclose(ds.origin_shift, 0, atol=1e-8, rtol=0)
    assert ds.wyckoffs[0] == letter

@pytest.mark.parametrize("transform,origin,expected", [
    (np.eye(3), [0, 0, 0], ["a"]),
    (np.eye(3), [0.5, 0, 0], None),
    (np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), [0, 0, 0], None),
])
def test_spglib_wyckoff_metadata_requires_unchanged_setting(monkeypatch, transform, origin, expected):
    import spglib

    hall, _ = hall_for_number(2)
    dataset = SimpleNamespace(number=2, hall_number=hall, wyckoffs=["a"],
                              transformation_matrix=transform, origin_shift=np.array(origin))
    def query(cell, **kwargs):
        assert kwargs.get("hall_number") == hall
        return dataset
    monkeypatch.setattr(spglib, "get_symmetry_dataset", query)
    structure = Structure(Lattice.cubic(4), ["Si"], [[0, 0, 0]])
    assert cif_writer._spglib_wyckoff_letters(structure, 2) == expected
