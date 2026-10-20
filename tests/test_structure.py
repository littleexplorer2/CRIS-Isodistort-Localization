"""Independent synthetic counterexamples; no website result or scientific case answer."""
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path
import unittest

from cris_rea.backend.structure import (
    AtomicSite, CrystalStructure, Lattice, StructureError, expand_structure,
    parse_cif, parse_cif_number, parse_symmetry_operation as op,
    validate_symmetry_operations,
)


def cif(operations=("x,y,z",), rows=("Si1 Si .137 .231 .319 1",),
        cell=("4.123", "5.287", "6.431", "78.2", "83.7", "74.6"),
        extra="", atom_extra_tags=""):
    tags = ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")
    lattice = "\n".join(f"_cell_{tag} {value}" for tag, value in zip(tags, cell))
    symmetry = ("loop_\n_space_group_symop_id\n_space_group_symop_operation_xyz\n"
                + "\n".join(f"{i+1} '{value}'" for i, value in enumerate(operations))) if operations else ""
    return (f"data_independent\n{lattice}\n{extra}\n{symmetry}\nloop_\n"
            "_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
            "_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
            f"{atom_extra_tags}\n" + "\n".join(rows) + "\n")


class StructureTests(unittest.TestCase):
    def error(self, code, call):
        with self.assertRaises(StructureError) as raised:
            call()
        self.assertEqual(raised.exception.code, code)
        return raised.exception

    def test_exact_numbers_and_uncertainty_exponent(self):
        number = parse_cif_number("-1.234e-2(5)")
        self.assertEqual(number.value, F(-617, 50000))
        self.assertEqual(number.standard_uncertainty, F(1, 20000))
        self.assertEqual(parse_cif_number("12(3)").standard_uncertainty, 3)
        self.assertEqual(parse_cif_number(".25").value, F(1, 4))
        for bad in ("?", ".", "NaN", "inf", "1/3", "1.2(3)e2"):
            self.error("invalid_number", lambda bad=bad: parse_cif_number(bad))
        self.error("number_limit", lambda: parse_cif_number("1e99999"))

    def test_parse_p1_preserves_input_and_metadata(self):
        source = cif(rows=("Fe1 Fe3+ -0.125(2) 1.25 .3 0.75(1) a 1",),
                     extra="_space_group_IT_number 1\n_space_group_name_H-M_alt 'P 1'",
                     atom_extra_tags="_atom_site_Wyckoff_symbol\n_atom_site_symmetry_multiplicity")
        before = source.encode("ascii")
        structure = parse_cif(before)
        self.assertEqual(structure.source_sha256, hashlib.sha256(before).hexdigest())
        self.assertEqual(structure.space_group_number, 1)
        site = structure.sites[0]
        self.assertEqual(site.symbol, "Fe3+")
        self.assertEqual(site.fractional, (F(-1, 8), F(5, 4), F(3, 10)))
        self.assertEqual(site.coordinate_uncertainties[0], F(1, 500))
        self.assertEqual(site.occupancy, F(3, 4))
        self.assertEqual(site.occupancy_uncertainty, F(1, 100))
        self.assertEqual(site.wyckoff_symbol, "a")
        expanded = expand_structure(structure)[0]
        self.assertEqual(expanded.fractional, (F(7, 8), F(1, 4), F(3, 10)))
        self.assertEqual(expanded.asymmetric_index, 0)
        self.assertEqual(before, source.encode("ascii"))
        self.assertIn("partial_single_component_occupancy:Fe1", structure.notices)

    def test_known_independent_converter_input(self):
        path = Path(__file__).resolve().parents[1] / "evidence/native/cif-converter-positive/synthetic_p1.cif"
        source = path.read_bytes()
        structure = parse_cif(source)
        self.assertEqual([site.symbol for site in structure.sites], ["Si", "O"])
        self.assertEqual(len(expand_structure(structure)), 2)
        self.assertEqual(path.read_bytes(), source)

    def test_cif_case_quotes_comments_and_multiline_metadata(self):
        source = cif(extra="_audit_note\n;a dog's unrelated text\nloop_ here is text\n;\n")
        source = source.replace("_cell_length_a", "_CELL_LENGTH_A").replace("data_independent", "DATA_Test")
        source = source.replace("_cell_length_b 5.287", "_cell_length_b 5.287 # comment")
        self.assertEqual(parse_cif(source.replace("\n", "\r\n")).data_name, "Test")

    def test_legacy_operation_alias(self):
        source = cif().replace("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
        self.assertEqual(len(parse_cif(source).symmetry_operations), 1)

    def test_missing_ops_even_with_group_one(self):
        self.error("missing_symmetry_operations",
                   lambda: parse_cif(cif(operations=(), extra="_space_group_IT_number 1")))

    def test_declared_group_is_metadata_not_lookup(self):
        structure = parse_cif(cif(extra="_space_group_IT_number 221"))
        self.assertEqual(structure.space_group_number, 221)
        self.assertEqual(len(structure.symmetry_operations), 1)
        # Caller must independently check that declared number matches operations.

    def test_invalid_and_ambiguous_cif(self):
        source = cif()
        self.error("multiple_data_blocks", lambda: parse_cif(source + "\ndata_another\n_a 1\n"))
        self.error("duplicate_cif_tag", lambda: parse_cif(source + "_cell_length_a 4\n"))
        self.error("cif_loop_cardinality", lambda: parse_cif(source.replace("Si1 Si .137 .231 .319 1", "Si1 Si .137")))
        self.error("quoted_numeric_value", lambda: parse_cif(source.replace("4.123", "'4.123'")))
        self.error("ambiguous_cif_alias", lambda: parse_cif(source + "_space_group_IT_number 1\n_symmetry_Int_Tables_number 1"))
        self.error("cif_syntax", lambda: parse_cif(source + "_unknown 'unterminated\n"))
        self.error("cif_syntax", lambda: parse_cif(source + "_unknown\n;unterminated\n"))
        self.error("unsupported_cif_version", lambda: parse_cif("#\\#CIF_2.0\n" + source))
        self.error("unsupported_cif_encoding", lambda: parse_cif(source + "# 中文"))

    def test_magnetic_modulated_disordered_and_missing_coordinates_rejected(self):
        for tag in ("_atom_site_moment_crystalaxis_x", "_cell_modulation_dimension",
                    "_atom_site_displace_Fourier_param_cos", "_space_group_ssg_name"):
            self.error("magnetic_or_modulated_unsupported",
                       lambda tag=tag: parse_cif(cif(extra=tag + " 0")))
        self.error("disorder_unsupported", lambda: parse_cif(cif(extra="_atom_site_disorder_group 1")))
        self.error("missing_cif_field", lambda: parse_cif(cif().replace("_atom_site_fract_x", "_atom_site_Cartn_x")))
        self.error("invalid_number", lambda: parse_cif(cif(rows=("Si1 Si ? 0 0 1",))))
        self.error("invalid_occupancy", lambda: parse_cif(cif(rows=("Si1 Si 0 0 0 1.2",))))

    def test_lattice_cartesian_and_gram(self):
        lattice = Lattice(4, 5, 6, 70, 80, 75)
        matrix = lattice.matrix
        for i, expected in enumerate((4, 5, 6)):
            self.assertAlmostEqual(math.sqrt(sum(matrix[j][i] ** 2 for j in range(3))), expected, places=12)
        metric = lattice.metric
        self.assertAlmostEqual(metric[0][1], 4 * 5 * math.cos(math.radians(75)), places=12)
        self.assertAlmostEqual(metric[0][2], 4 * 6 * math.cos(math.radians(80)), places=12)
        self.assertAlmostEqual(metric[1][2], 5 * 6 * math.cos(math.radians(70)), places=12)
        x = (F(1, 7), F(-3, 11), F(4, 5))
        restored = lattice.cartesian_to_fractional(lattice.fractional_to_cartesian(x))
        for actual, expected in zip(restored, x):
            self.assertAlmostEqual(actual, float(expected), places=13)
        ca, cb, cg = (math.cos(math.radians(value)) for value in (70, 80, 75))
        expected_volume = 4 * 5 * 6 * math.sqrt(1 + 2*ca*cb*cg - ca*ca - cb*cb - cg*cg)
        self.assertAlmostEqual(lattice.volume, expected_volume, places=11)

    def test_invalid_cells(self):
        for args in ((0, 1, 1, 90, 90, 90), (1, 1, 1, 0, 90, 90)):
            self.error("invalid_lattice", lambda args=args: Lattice(*args))
        self.error("singular_lattice", lambda: Lattice(1, 1, 1, 10, 10, 170))
        self.error("invalid_lattice", lambda: Lattice("1e300", "1e300", "1e300", 90, 90, 90))
        self.error("invalid_lattice", lambda: Lattice("1e-300", "1e-300", "1e-300", 90, 90, 90))

    def test_symmetry_composition_and_exact_involution(self):
        glide = op("x, 1/2-y, z+1/2")
        self.assertEqual(glide.compose(glide), op("x,y,z"))
        position = (F(1, 7), F(2, 9), F(3, 11))
        self.assertEqual(glide.compose(glide).apply(position), position)
        g = op("-x,y+1/2,-z")
        h = op("x+1/2,-y,-z")
        self.assertEqual(g.compose(h).apply(position), g.apply(h.apply(position)))

    def test_no_eval_and_nonlinear_or_nonlattice_operation(self):
        for value in ("x*y,y,z", "x/(y+1),y,z", "x**2,y,z"):
            self.error("nonaffine_symmetry_operation", lambda value=value: op(value))
        self.error("invalid_symmetry_operation", lambda: op("__import__('os'),y,z"))
        self.error("unsupported_symmetry_dimension", lambda: op("x,y,z,+1"))
        self.error("nonintegral_rotation", lambda: op("x/2,y,z"))
        self.error("nonunimodular_rotation", lambda: op("2*x,y,z"))
        self.error("invalid_symmetry_operation", lambda: op("x+1e99999999,y,z"))

    def test_exact_group_and_decimal_translation_diagnostic(self):
        exact = tuple(map(op, ("x,y,z", "x+1/3,y,z", "x+2/3,y,z")))
        validate_symmetry_operations(exact, Lattice(1, 2, 3, 90, 90, 90))
        rounded = tuple(map(op, ("x,y,z", "x+.333333,y,z", "x+.666667,y,z")))
        error = self.error("nonclosed_symmetry_operations", lambda: validate_symmetry_operations(rounded))
        self.assertIn("not silently changed", str(error))
        self.error("duplicate_symmetry_operation", lambda: validate_symmetry_operations((op("x,y,z"), op("x+1,y,z"))))
        self.error("missing_identity", lambda: validate_symmetry_operations((op("-x,-y,-z"),)))

    def test_hexagonal_metric_and_wrong_orthogonal_cell_counterexample(self):
        operations = tuple(map(op, ("x,y,z", "-y,x-y,z", "-x+y,-x,z")))
        validate_symmetry_operations(operations, Lattice(3, 3, 5, 90, 90, 120))
        self.error("symmetry_metric_mismatch",
                   lambda: validate_symmetry_operations(operations, Lattice(3, 3, 5, 90, 90, 90)))

    def test_orbit_stabilizer_and_declared_multiplicity(self):
        operations = ("x,y,z", "-x,-y,-z")
        source = cif(operations, rows=("Si1 Si 0 0 0 1 a 1", "O1 O .1 .2 .3 1 i 2"),
                     atom_extra_tags="_atom_site_Wyckoff_symbol\n_atom_site_symmetry_multiplicity")
        structure = parse_cif(source)
        sites = expand_structure(structure)
        self.assertEqual(len(sites), 3)
        self.assertEqual(sites[0].operation_indices, (0, 1))
        self.assertEqual(sites[1].operation_indices, (0,))
        self.assertEqual(sites[2].fractional, (F(9, 10), F(4, 5), F(7, 10)))
        self.error("multiplicity_mismatch", lambda: parse_cif(source.replace("1 a 1", "1 a 2")))

    def test_centering_orbit_and_covariance(self):
        source = cif(("x,y,z", "x+1/2,y+1/2,z"))
        structure = parse_cif(source)
        orbit = expand_structure(structure)
        self.assertEqual(len(orbit), 2)
        positions = {site.fractional for site in orbit}
        for operation in structure.symmetry_operations:
            self.assertEqual({operation.apply(site.fractional) for site in orbit}, positions)

    def test_mixed_or_overlapping_representatives_rejected(self):
        self.error("mixed_occupancy_unsupported",
                   lambda: parse_cif(cif(rows=("Fe1 Fe 0 0 0 .5", "Co1 Co 1 0 0 .5"))))
        self.error("overlapping_independent_sites",
                   lambda: parse_cif(cif(("x,y,z", "-x,-y,-z"),
                                        rows=("Fe1 Fe .1 .2 .3 1", "Fe2 Fe -.1 -.2 -.3 1"))))
        self.error("duplicate_atom_label",
                   lambda: parse_cif(cif(rows=("A Si 0 0 0 1", "A Si .1 .2 .3 1"))))

    def test_public_objects_accept_exact_inputs_and_do_not_mutate(self):
        sites = [AtomicSite("A", "Si", (0, 0, 0))]
        operations = [op("x,y,z")]
        structure = CrystalStructure("manual", Lattice(1, 1, 1, 90, 90, 90), sites, operations)
        sites.clear()
        operations.clear()
        self.assertEqual(len(expand_structure(structure)), 1)

    def test_modern_and_legacy_multiplicity_are_traceable(self):
        modern = cif(rows=("Si1 Si 0 0 0 1 1",),
                     atom_extra_tags="_atom_site_site_symmetry_multiplicity")
        site = parse_cif(modern).sites[0]
        self.assertEqual(site.multiplicity_source_tag, "_atom_site_site_symmetry_multiplicity")
        legacy = modern.replace("_atom_site_site_symmetry_multiplicity", "_atom_site_symmetry_multiplicity")
        self.assertIn("legacy_multiplicity_tag:Si1", parse_cif(legacy).notices)
        bad = legacy.replace("Si1 Si 0 0 0 1 1", "Si1 Si 0 0 0 1 2")
        error = self.error("multiplicity_mismatch", lambda: parse_cif(bad))
        self.assertIn("historically ambiguous", str(error))

    def test_numeric_uncertainty_is_retained_not_used_for_snapping(self):
        source = cif(("x,y,z", "-x,-y,-z"), rows=("A Si .00001(10) 0 0 1",))
        structure = parse_cif(source)
        self.assertEqual(structure.sites[0].coordinate_uncertainties[0], F(1, 10000))
        self.assertEqual(len(expand_structure(structure)), 2)
        # An interval containing zero does not authorize changing the central value.

    def test_metrics_use_explicit_bounded_tolerance(self):
        operations = tuple(map(op, ("x,y,z", "y,x,z")))
        cell = Lattice(1, "1.00000001", 2, 90, 90, 90)
        self.error("symmetry_metric_mismatch", lambda: validate_symmetry_operations(operations, cell))
        validate_symmetry_operations(operations, cell, metric_relative_tolerance=1e-7)
        for tolerance in (-1, float("nan"), 1):
            self.error("invalid_tolerance", lambda tolerance=tolerance:
                       validate_symmetry_operations((op("x,y,z"),), cell,
                                                    metric_relative_tolerance=tolerance))
        self.error("symmetry_metric_mismatch", lambda: validate_symmetry_operations(
            operations, Lattice(1, 2, 10000, 90, 90, 90)))


if __name__ == "__main__":
    unittest.main()
