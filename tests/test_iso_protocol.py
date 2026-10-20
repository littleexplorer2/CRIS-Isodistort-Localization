"""Protocol tests against immutable official binary output, no website answers."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from cris_rea.backend.iso_protocol import (
    IsoProtocolError, build_direction_query, build_irrep_query, build_isotropy_query,
    build_parent_query, parse_irrep_output, parse_isotropy_output, parse_parent_output,
    parse_seitz_operation,
)

BASE = Path(__file__).resolve().parents[1] / "evidence/phase2/iso-contract"


def output(case):
    return (BASE / case / "stdout.txt").read_text(encoding="utf-8")


class IsoProtocolTests(unittest.TestCase):
    def test_actual_irrep_parameter_dimensions_zero_through_three(self):
        cases = (
            ("irrep-metadata-221-gm4-", 221, "GM4-", (), 3, (), (0, 0, 0)),
            ("irrep-metadata-221-x1+", 221, "X1+", (), 3, (), (0, F(1, 2), 0)),
            ("irrep-metadata-221-dt1", 221, "DT1", ("a",), 6, ("1/4",), (0, F(1, 4), 0)),
            ("irrep-metadata-225-ld1", 225, "LD1", ("a",), 8, ("1/5",), (F(1, 5),) * 3),
            ("irrep-metadata-221-c1", 221, "C1", ("a", "b"), 24,
             ("1/3", "1/4"), (F(1, 3), F(1, 3), F(1, 4))),
            ("irrep-metadata-1-gp1gq1", 1, "GP1GQ1", ("a", "b", "g"), 2,
             ("1/3", "1/4", 0), (F(1, 3), F(1, 4), 0)),
        )
        for case, parent, irrep, symbols, dimension, parameters, k in cases:
            with self.subTest(case=case):
                result = parse_irrep_output(output(case), parent, irrep)
                self.assertEqual(result.parameter_symbols, symbols)
                self.assertEqual(result.dimension, dimension)
                self.assertEqual(result.k_vector(parameters), k)
                self.assertTrue(all(isinstance(x, F) for x in result.k_vector(parameters)))

    def test_parameter_symbols_depend_on_expression_not_irrep_label(self):
        raw = output("irrep-metadata-1-gp1gq1").replace("(a,b,g)", "(b,a,1/2-g/2)")
        result = parse_irrep_output(raw)
        self.assertEqual(result.parameter_symbols, ("a", "b", "g"))
        self.assertEqual(result.wavevector_expression, "(b,a,1/2-g/2)")
        self.assertEqual(result.k_vector(["1/3", "1/4", "1/5"]),
                         (F(1, 4), F(1, 3), F(2, 5)))
        constant = parse_irrep_output(output("irrep-metadata-1-gp1gq1").replace(
            "(a,b,g)", "(a-a,0,1/2)"))
        self.assertEqual(constant.parameter_symbols, ())
        self.assertEqual(constant.k_vector(), (0, 0, F(1, 2)))

    def test_missing_extra_and_inexact_k_parameters_rejected(self):
        free = parse_irrep_output(output("irrep-metadata-1-gp1gq1"))
        for parameters in ((), [1], [1, 2], [1, 2, 3, 4], [1, 2, 0.5], [1, 2, "0.5"]):
            with self.subTest(parameters=parameters), self.assertRaises(IsoProtocolError):
                free.k_vector(parameters)
        special = parse_irrep_output(output("irrep-metadata-221-gm4-"))
        with self.assertRaises(IsoProtocolError):
            special.k_vector([0])

    def test_irrep_metadata_rejects_wrong_selection_or_non_linear_expression(self):
        raw = output("irrep-metadata-1-gp1gq1")
        for expression in ("(a*b,b,g)", "(x,b,g)", "(0.5a,b,g)",
                           "(a,b)", "(a+,b,g)", "(a,b,g/0)"):
            with self.subTest(expression=expression), self.assertRaises(IsoProtocolError):
                parse_irrep_output(raw.replace("(a,b,g)", expression))
        with self.assertRaises(IsoProtocolError):
            parse_irrep_output(raw, expected_parent=2)
        with self.assertRaises(IsoProtocolError):
            parse_irrep_output(raw, expected_irrep="GM1")
        with self.assertRaises(IsoProtocolError):
            parse_irrep_output(output("discover-221-dt"))

    def test_irrep_query_only_requests_symbolic_metadata(self):
        commands = build_irrep_query(221, "DT1").decode()
        self.assertIn("SHOW KPOINT\nSHOW DIMENSION\nDISPLAY IRREP", commands)
        self.assertNotIn("VALUE KVALUE", commands)
        with self.assertRaises(IsoProtocolError):
            build_irrep_query(221, "DT1\nQUIT")

    def test_parent_face_centered_basis_and_wrapped_operations(self):
        result = parse_parent_output(output("parent-225-primitive-basis-elements"), 225)
        self.assertEqual(len(result.operations), 48)
        self.assertEqual(result.primitive_basis,
                         ((0, F(1, 2), F(1, 2)), (F(1, 2), 0, F(1, 2)),
                          (F(1, 2), F(1, 2), 0)))
        self.assertTrue(all(isinstance(x, F) for row in result.primitive_basis for x in row))

    def test_rhombohedral_basis_has_printed_vectors_as_columns(self):
        result = parse_parent_output(output("parent-160-primitive-basis-elements"), 160)
        self.assertEqual(result.primitive_basis,
                         ((F(2, 3), F(-1, 3), F(-1, 3)),
                          (F(1, 3), F(1, 3), F(-2, 3)),
                          (F(1, 3), F(1, 3), F(1, 3))))
        self.assertEqual(len(result.operations), 6)

    def test_primary_basis_orientation_and_metadata(self):
        result = parse_isotropy_output(output("elements-221-gm4minus-p1"), 221, "GM4-", "P1")
        self.assertEqual(result.basis, ((0, 0, 1), (1, 0, 0), (0, 1, 0)))
        self.assertEqual((result.subgroup_number, result.index, result.size), (99, 6, 1))
        self.assertEqual(result.opd_vector, "(a,0,0)")
        self.assertEqual(len(result.operations), 8)

    def test_centered_size_is_not_conventional_determinant(self):
        result = parse_isotropy_output(output("subgroup-221-gm4minus-p3-centered"))
        self.assertEqual(result.size, 1)
        self.assertEqual(result.basis, ((-1, 0, -1), (0, 1, -1), (1, -1, -1)))
        self.assertEqual(result.subgroup_number, 160)

    def test_generated_rational_k_table_banner_allowed(self):
        result = parse_isotropy_output(output("gp-1-rational-generate"), 1, "GP1GQ1", "C1")
        self.assertEqual((result.size, result.index), (12, 12))
        self.assertEqual(result.basis, ((0, 0, -3), (0, 4, 0), (1, 0, 0)))

    def test_live_builder_outputs_and_explicit_setting_failure(self):
        for case in ("builder-primary-gamma221-default-setting",
                     "screen-1000-primary-default-setting"):
            self.assertEqual(parse_isotropy_output(output(case)).subgroup_number, 99)
        self.assertEqual(parse_isotropy_output(
            output("builder-primary-gp-k-default-setting")).size, 12)
        self.assertEqual(parse_parent_output(output("builder-parent-225-default-setting")).number, 225)
        for case in ("builder-primary-gamma221", "builder-decline-gp-k-default-setting"):
            with self.subTest(case=case), self.assertRaises(IsoProtocolError):
                parse_isotropy_output(output(case))
        with self.assertRaises(IsoProtocolError):
            parse_parent_output(output("builder-parent-225-default-setting").replace(
                "International (new ed.)", "International (old ed.)"))

    def test_live_wide_output_and_nonzero_subgroup_origin(self):
        self.assertEqual(parse_parent_output(output("builder-parent-225-wide")).operations,
                         parse_parent_output(output("parent-225-primitive-basis-elements")).operations)
        result = parse_isotropy_output(output("builder-dt-quarter-p2-wide"), 221, "DT1", "P2")
        self.assertEqual(result.origin, (0, F(-1, 2), 0))
        self.assertEqual(result.basis, ((1, 0, 0), (0, 0, -4), (0, 1, 0)))
        self.assertEqual((result.index, result.size, len(result.operations)), (12, 4, 16))
        self.assertIn((0, -1, 0), [operation.translation for operation in result.operations])

    def test_integer_translations_are_preserved_for_supercell_phases(self):
        op = parse_seitz_operation("(x,-y-1,-z+3/2)")
        self.assertEqual(op.rotation, ((1, 0, 0), (0, -1, 0), (0, 0, -1)))
        self.assertEqual(op.translation, (0, -1, F(3, 2)))
        hexagonal = parse_seitz_operation("(-y,x-y,z)")
        self.assertEqual(hexagonal.rotation, ((0, -1, 0), (1, -1, 0), (0, 0, 1)))

    def test_reject_non_exact_or_invalid_operations(self):
        for value in ("(x,y,0.5+z)", "(x,y,z+1/0)", "(x,y,0)", "(x,y,z);print(1)",
                      "(x,y,z--1)", "(x,y,z**2)", "(x,y,z,0)"):
            with self.subTest(value=value), self.assertRaises(IsoProtocolError):
                parse_seitz_operation(value)

    def test_multiple_rows_and_declined_generation_are_not_a_result(self):
        for case in ("dt-221-quarter-generate", "dt-221-quarter-subgroups",
                     "gamma-221-forbidden-scalar-irrep"):
            with self.subTest(case=case), self.assertRaises(IsoProtocolError):
                parse_isotropy_output(output(case))

    def test_missing_fields_wrong_selection_and_truncation(self):
        raw = output("elements-221-gm4minus-p1")
        invalid = (raw.replace("(x,-z,-y)", "(x,-z,-y"),
                   raw.replace("(0,0,0)", "(0.0,0,0)"),
                   raw.replace("(x,-z,-y)", "*"),
                   raw + "\n***Syntax error: this command\n")
        for value in invalid:
            with self.subTest(value=value[-80:]), self.assertRaises(IsoProtocolError):
                parse_isotropy_output(value)
        for kwargs in ({"expected_parent": 225}, {"expected_irrep": "GM1+"},
                       {"expected_opd": "P3"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(IsoProtocolError):
                parse_isotropy_output(raw, **kwargs)
        with self.assertRaises(IsoProtocolError):
            parse_parent_output(output("parent-225-primitive-basis-elements"), 221)

    def test_builder_exact_k_and_explicit_generation_confirmation(self):
        commands = build_isotropy_query(1, "gp1gq1", "c1", [F(1, 3), "1/4", 0]).decode()
        self.assertIn("VALUE KVALUE 3,1/3,1/4,0\n", commands)
        self.assertIn("LABEL ELEMENT INTERNATIONAL\n", commands)
        self.assertIn("DISPLAY ISOTROPY\n\nQUIT\n", commands)
        self.assertIn("VALUE IRREP GP1GQ1", commands)
        decline = build_isotropy_query(221, "GM4-", "P1", allow_generate=False).decode()
        self.assertTrue(decline.endswith("DISPLAY ISOTROPY\nQUIT\nQUIT\n"))

    def test_direction_builder_basis_transpose(self):
        commands = build_direction_query(
            221, 99, ((0, 0, 1), (1, 0, 0), (0, 1, 0)), (0, F(1, 2), 0)
        ).decode()
        self.assertIn("VALUE BASIS 0,1,0 0,0,1 1,0,0\n", commands)
        self.assertIn("VALUE ORIGIN 0,1/2,0\n", commands)
        self.assertIn("DISPLAY DIRECTION", commands)

    def test_command_injection_inexact_inputs_and_invalid_shapes_rejected(self):
        actions = [
            lambda: build_parent_query(True), lambda: build_parent_query(231),
            lambda: build_isotropy_query(221, "GM4-\nQUIT", "P1"),
            lambda: build_isotropy_query(221, "GM4-", "P1;QUIT"),
            lambda: build_isotropy_query(1, "GP1GQ1", "C1", [0.25]),
            lambda: build_isotropy_query(1, "GP1GQ1", "C1", ["0.25"]),
            lambda: build_isotropy_query(1, "GP1GQ1", "C1", ["1/0"]),
            lambda: build_isotropy_query(1, "GP1GQ1", "C1", [1, 2, 3, 4]),
            lambda: build_direction_query(1, 1, ((0, 0, 0),) * 3, (0, 0, 0)),
        ]
        for action in actions:
            with self.subTest(action=action), self.assertRaises(IsoProtocolError):
                action()


if __name__ == "__main__":
    unittest.main()
