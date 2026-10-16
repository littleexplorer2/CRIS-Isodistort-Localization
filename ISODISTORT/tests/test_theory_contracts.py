"""Independent algebraic counterexamples from the theoretical design audit.

These synthetic cases test contracts, not agreement with saved website answers.
"""
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from application.session import IsoDistort
from backend.method1 import OccupationalModeGenerator
from backend.method4 import DistortionEngine
from backend.method1.search_methods import _basis_is_same_lattice
from backend.utils.lattice import is_sublattice, rational_matrix, same_lattice
from backend.export.isodistort_cif import _wyckoff_letter
from backend.export import isodistort_cif as cif_writer
from pymatgen.symmetry.groups import SpaceGroup
from backend.wrappers import SubgroupInfo
from backend.wrappers import KPointInfo
from backend.tables.kpoints_official import official_special_k_coords, official_kparams_to_iso


def _session():
    iso = IsoDistort.__new__(IsoDistort)
    rotation = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
    iso._parent_rotations = lambda: [np.eye(3, dtype=int), rotation, rotation @ rotation]
    return iso, rotation


def test_row_lattice_orbit_uses_transposed_fractional_rotation():
    iso, rotation = _session()
    basis = np.diag([3, 1, 1])
    assert iso._same_lattice_orbit(basis @ rotation.T, basis)
    assert not iso._same_lattice_orbit(basis @ rotation, basis)
    rotations = iso._parent_rotations()
    assert _basis_is_same_lattice(basis @ rotation.T, basis, rotations)
    assert not _basis_is_same_lattice(basis @ rotation, basis, rotations)
    oriented = rational_matrix(basis @ rotation.T)
    assert is_sublattice(oriented, oriented)
    assert not is_sublattice(rational_matrix(basis @ rotation), oriented)


def test_gl3z_equivalence_has_no_relative_float_tolerance():
    basis = [[1, Fraction(1000000001, 10), 0], [0, 1, 0], [0, 0, 1]]
    assert not same_lattice(rational_matrix(basis), rational_matrix(np.eye(3, dtype=int)))
    basis[0][1] = 100000000
    assert same_lattice(rational_matrix(basis), rational_matrix(np.eye(3, dtype=int)))


def test_lattice_label_does_not_change_high_denominator_basis():
    label = IsoDistort._format_lattice(((Fraction(1, 17), 0, 0), (0, 1, 0), (0, 0, 1)))
    assert label == "(1/17,0,0),(0,1,0),(0,0,1)"


def test_method1_options_do_not_read_reference_archives(monkeypatch):
    iso, _ = _session()
    iso.structure = Structure(Lattice.cubic(4), ["Si"], [[0.13, 0.27, 0.39]])
    iso.symmetry_info = {"space_group_number": 1}
    subgroup = SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1",
                            basis_vectors=np.eye(3, dtype=int).tolist(),
                            k_point_label="GM", irrep_label="GM1")
    iso._method1_filtered_subgroups = lambda: [subgroup]
    original = Path.iterdir

    def guarded_iterdir(path):
        assert path.name not in {"webpage_info", "experiment_data", "output_compare"}
        return original(path)

    monkeypatch.setattr(Path, "iterdir", guarded_iterdir)
    assert iso.method1_options()["conventional_lattices"]


def test_unsupported_occupational_modes_do_not_partition_coordinates():
    parent = Structure(Lattice.cubic(4), ["Si"], [[0, 0, 0]])
    subgroup = SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1",
                            basis_vectors=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                            k_point_label="X", irrep_label="X1")
    sites = [{"wyckoff_letter": "a", "species": "Si", "equivalent_indices": [0]}]
    with pytest.raises(NotImplementedError, match="occupational"):
        OccupationalModeGenerator().generate(parent, sites, subgroup, {"Si"})


@pytest.mark.parametrize("amplitude", [0.0, 0.2, 2.0])
def test_legacy_binary_occupancy_contract_cannot_rewrite_parent(amplitude):
    parent = Structure(Lattice.cubic(4), [{"Si": 0.6}], [[0, 0, 0]])
    # The old path changed 0.6 to 1 even at zero amplitude and clipped at 0.01.
    with pytest.raises(NotImplementedError, match="occupational"):
        DistortionEngine().generate_modes(parent, occupancy_patterns=[(np.ones(1), amplitude)])


@pytest.mark.parametrize("kind", ["rotational", "occupational", "magnetic"])
def test_unproven_physical_types_are_not_inferred_from_displacement(kind):
    iso, _ = _session()
    iso.distortion_types = [kind]
    with pytest.raises(NotImplementedError, match=kind):
        iso._filter_method1_by_types([], [kind])


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

    hall, _ = cif_writer._hall_for_number(sg)
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

    hall, _ = cif_writer._hall_for_number(2)
    dataset = SimpleNamespace(number=2, hall_number=hall, wyckoffs=["a"],
                              transformation_matrix=transform, origin_shift=np.array(origin))
    def query(cell, **kwargs):
        assert kwargs.get("hall_number") == hall
        return dataset
    monkeypatch.setattr(spglib, "get_symmetry_dataset", query)
    structure = Structure(Lattice.cubic(4), ["Si"], [[0, 0, 0]])
    assert cif_writer._spglib_wyckoff_letters(structure, 2) == expected


def test_affine_wavevector_substitution_keeps_offsets_and_signs():
    assert official_special_k_coords(139, "Y", [], ["1/3"]) == ["1/3", "2/3", "0"]


def test_kvalue_conversion_solves_the_coordinate_equation():
    # Both affine coordinates, rather than their textual coefficients alone,
    # must denote the same reciprocal vector.
    kp = KPointInfo(label="LD", coordinates=["0", "0", "1-2a"],
                   parameters=["2a"], is_special=False)
    assert official_kparams_to_iso(139, "LD", ["1/6"], kp) == ["5/12"]


def test_one_free_variable_is_not_counted_as_two_coordinate_expressions():
    kp = KPointInfo(label="Y", coordinates=["a", "1-a", "0"],
                   parameters=[], is_special=True)
    assert kp.parameters == ["a"]
    assert not kp.is_special


def test_kvalue_conversion_rejects_inconsistent_coordinates():
    kp = KPointInfo(label="DT", coordinates=["a", "2a", "0"],
                   parameters=["a"], is_special=False)
    with pytest.raises(ValueError, match="Inconsistent"):
        official_kparams_to_iso(139, "DT", ["1/4"], kp)


def test_coupled_k_parameters_are_solved_as_one_linear_system():
    kp = KPointInfo(label="A", coordinates=["a+b", "a+b", "2a"],
                   parameters=[], is_special=False)
    values = official_kparams_to_iso(139, "A", ["1/7", "1/9"], kp)
    a, b = map(Fraction, values)
    assert a + b == Fraction(1, 7)
    assert 2 * a == Fraction(1, 9)
