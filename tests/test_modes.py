"""Algebraic counterexamples independent of official case results."""
from fractions import Fraction as F
from types import SimpleNamespace as NS
import unittest

import numpy as np

from cris_rea.backend.modes import construct_modes, distort, ModeError
from cris_rea.backend.modes.core import inverse, mm, rational_matrix, IDENTITY, close_group
from cris_rea.backend.modes.export import cif_text, checked_cif
from cris_rea.backend.structure import (
    Lattice, AtomicSite, CrystalStructure, SymmetryOperation, parse_cif, expand_structure,
)

I = rational_matrix(IDENTITY)
Z = (F(0),) * 3
ONE = SymmetryOperation(IDENTITY, Z)


def fixture(ops=(ONE,), cell=I, sites=None, parent_ops=None, hprimitive=I,
            parentprimitive=I, index=1, size=1, lattice=None):
    parent_ops = parent_ops or ops
    structure = CrystalStructure("algebra", lattice or Lattice(4, 5, 6, 80, 85, 75),
                                 sites or (AtomicSite("A", "Si", (F(1, 7), F(2, 7), F(3, 7))),),
                                 parent_ops, 1)
    parent = NS(operations=parent_ops, primitive_basis=parentprimitive)
    subgroup = NS(operations=ops, basis=cell, origin=Z, index=index, size=size,
                  subgroup_number=1, subgroup_symbol="P1")
    hparent = NS(primitive_basis=hprimitive)
    return structure, parent, subgroup, hparent


class ModeTests(unittest.TestCase):
    def test_fraction_inverse_and_shear(self):
        a = rational_matrix(((1, F(1, 2), 2), (0, 2, 1), (0, 0, 3)))
        self.assertEqual(mm(a, inverse(a)), I)

    def test_triclinic_p1_all_nine_dof_and_six_strains(self):
        sites = tuple(AtomicSite(f"A{i}", "Si", (F(i, 7), F(i, 11), F(i, 13)))
                      for i in range(3))
        model = construct_modes(*fixture(sites=sites))
        self.assertEqual(model["basis"].shape, (9, 9))
        self.assertEqual(model["internal_modes"], 6)
        self.assertEqual(model["rigid_translations"], 3)
        self.assertEqual(model["strain_basis"].shape, (6, 6))
        np.testing.assert_allclose(model["basis"]@model["basis"].T, np.eye(9), atol=1e-13)

    def test_inversion_center_has_no_polar_dof(self):
        inversion = SymmetryOperation(tuple(tuple(-v for v in row) for row in IDENTITY), Z)
        model = construct_modes(*fixture(ops=(ONE, inversion), sites=(AtomicSite("A", "Si", Z),)))
        self.assertEqual(model["basis"].shape, (3, 0))
        self.assertEqual(model["rigid_translations"], 0)
        self.assertEqual(model["strain_basis"].shape, (6, 6))
        distort(model)

    def test_inversion_related_pair_is_opposite(self):
        inversion = SymmetryOperation(tuple(tuple(-v for v in row) for row in IDENTITY), Z)
        model = construct_modes(*fixture(ops=(ONE, inversion)))
        self.assertEqual(model["basis"].shape, (6, 3))
        np.testing.assert_allclose(model["basis"][:3], -model["basis"][3:], atol=1e-13)
        self.assertEqual(model["rigid_translations"], 0)

    def test_nonzero_translation_phase_is_not_lost(self):
        # Same parent point group, doubled x cell with H lacking parent translation x.
        c = rational_matrix(((2, 0, 0), (0, 1, 0), (0, 0, 1)))
        model = construct_modes(*fixture(cell=c, index=2, size=2))
        self.assertEqual(len(model["atoms"]), 2)
        self.assertEqual(model["basis"].shape, (6, 6))
        self.assertEqual(model["parent_primitive_cells"], 2)
        np.testing.assert_allclose(model["basis"][:3, :3], -model["basis"][3:, :3], atol=1e-13)

    def test_non_diagonal_supercell_count(self):
        c = rational_matrix(((2, 1, 0), (0, 2, 1), (0, 0, 1)))
        model = construct_modes(*fixture(cell=c, index=4, size=4))
        self.assertEqual(len(model["atoms"]), 4)
        self.assertEqual(model["basis"].shape, (12, 12))

    def test_centering_volume_and_mode_count(self):
        p = rational_matrix(((0, F(1, 2), F(1, 2)),
                             (F(1, 2), 0, F(1, 2)), (F(1, 2), F(1, 2), 0)))
        ops = close_group([SymmetryOperation(IDENTITY, tuple(p[i][j] for i in range(3)))
                           for j in range(3)])
        model = construct_modes(*fixture(ops=ops, hprimitive=p, parentprimitive=p))
        self.assertEqual(len(model["atoms"]), 4)
        self.assertEqual(model["basis"].shape, (12, 3))
        self.assertEqual(model["parent_primitive_cells"], 4)
        np.testing.assert_allclose(model["basis"][0:3], model["basis"][9:12], atol=1e-13)

    def test_c4_axis_and_kelvin_space(self):
        rotation = SymmetryOperation(((0, -1, 0), (1, 0, 0), (0, 0, 1)), Z)
        ops = close_group([rotation])
        model = construct_modes(*fixture(ops=ops, sites=(AtomicSite("A", "Si", Z),),
                                         lattice=Lattice(4, 4, 6, 90, 90, 90)))
        self.assertEqual(model["basis"].shape, (3, 1))
        self.assertEqual(model["strain_basis"].shape, (6, 2))
        np.testing.assert_allclose(model["basis"][:, 0], (0, 0, 1), atol=1e-13)
        state = distort(model, [0.2], [0.01, -0.02])
        self.assertAlmostEqual(state["strain_tensor"][0, 0], state["strain_tensor"][1, 1])

    def test_amplitude_inverse_and_rss(self):
        model = construct_modes(*fixture(sites=(AtomicSite("A", "Si", Z),
                                                AtomicSite("B", "O", (F(1, 3),)*3))))
        coefficients = np.linspace(-0.06, 0.06, 6)
        state = distort(model, coefficients, [0.01, 0.02, -0.01, 0.003, 0.004, -0.005])
        u = state["displacements_angstrom"].ravel()
        np.testing.assert_allclose(model["basis"].T@u, coefficients, atol=1e-14)
        self.assertAlmostEqual(float(u@u), float(coefficients@coefficients), places=14)
        self.assertAlmostEqual(float(np.sum(state["strain_tensor"]**2)),
                               float(np.sum(state["strain_amplitudes"]**2)), places=14)

    def test_cif_asymmetric_unit_no_duplicates(self):
        inv = SymmetryOperation(tuple(tuple(-v for v in row) for row in IDENTITY), Z)
        args = fixture(ops=(ONE, inv))
        model = construct_modes(*args)
        state = distort(model, [0.01, -0.02, 0.03])
        exported = parse_cif(cif_text(model, state, args[2], "roundtrip"))
        self.assertEqual(len(exported.sites), 1)
        self.assertEqual(len(expand_structure(exported)), 2)
        self.assertAlmostEqual(exported.lattice.volume, np.linalg.det(state["lattice"]), places=10)

    def test_threefold_special_position_uses_verified_full_cell_cif(self):
        rotation = SymmetryOperation(((0, -1, 0), (1, -1, 0), (0, 0, 1)), Z)
        ops = close_group([rotation])
        args = fixture(ops=ops, sites=(AtomicSite("A", "Si", (F(1, 3), F(2, 3), F(1, 7))),),
                       lattice=Lattice(4, 4, 6, 90, 90, 120))
        args[2].subgroup_number = 143
        args[2].subgroup_symbol = "P3"
        model = construct_modes(*args)
        text, info = checked_cif(model, distort(model), args[2], "special_position")
        self.assertEqual(info["representation"], "P1_full_cell")
        self.assertEqual(info["preserved_subgroup_number"], 143)
        self.assertEqual(len(expand_structure(parse_cif(text))), 1)

    def test_incompatible_cell_rejected(self):
        c = rational_matrix(((F(1, 2), 0, 0), (0, 1, 0), (0, 0, 1)))
        with self.assertRaisesRegex(ModeError, "not a parent translation"):
            construct_modes(*fixture(cell=c))

    def test_wrong_size_rejected(self):
        c = rational_matrix(((2, 0, 0), (0, 1, 0), (0, 0, 1)))
        with self.assertRaisesRegex(ModeError, "primitive volume"):
            construct_modes(*fixture(cell=c, size=1, index=2))

    def test_wrong_index_rejected(self):
        with self.assertRaisesRegex(ModeError, "ISO index"):
            construct_modes(*fixture(index=2))

    def test_partial_occupancy_rejected(self):
        with self.assertRaisesRegex(ModeError, "full occupancy"):
            construct_modes(*fixture(sites=(AtomicSite("A", "Si", Z, F(1, 2)),)))

    def test_nonfinite_amplitude_and_bad_strain_rejected(self):
        model = construct_modes(*fixture())
        with self.assertRaises(ModeError):
            distort(model, [float("nan"), 0, 0])
        with self.assertRaisesRegex(ModeError, "positive definite"):
            distort(model, strain_amplitudes=[-2, 0, 0, 0, 0, 0])

    def test_mismatched_parent_setting_rejected(self):
        structure, parent, subgroup, hp = fixture()
        parent.operations = (ONE, SymmetryOperation(IDENTITY, (F(1, 2), F(1, 2), 0)))
        parent.primitive_basis = rational_matrix(((F(1, 2), F(-1, 2), 0),
                                                 (F(1, 2), F(1, 2), 0), (0, 0, 1)))
        with self.assertRaisesRegex(ModeError, "setting"):
            construct_modes(structure, parent, subgroup, hp)


if __name__ == "__main__":
    unittest.main()
