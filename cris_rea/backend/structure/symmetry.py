"""Parse explicit affine xyz expressions using a restricted AST, never eval."""
from __future__ import annotations

import ast
from fractions import Fraction
import re

from .model import StructureError, SymmetryOperation


def parse_symmetry_operation(expression: str) -> SymmetryOperation:
    if not isinstance(expression, str) or len(expression) > 1024:
        raise StructureError("invalid_symmetry_operation", "A bounded xyz operation string is required.")
    components = expression.lower().split(",")
    if len(components) != 3:
        raise StructureError("unsupported_symmetry_dimension", "Only three-dimensional xyz operations are supported.")
    rows = []
    translation = []
    for component in components:
        component = component.strip()
        if not component or not re.fullmatch(r"[0-9xyz+\-*/.()\seE]+", component):
            raise StructureError("invalid_symmetry_operation", "Only affine xyz arithmetic is permitted.")
        try:
            tree = ast.parse(component, mode="eval")
        except (SyntaxError, ValueError, RecursionError) as exc:
            raise StructureError("invalid_symmetry_operation", "Malformed affine expression.") from exc
        if sum(1 for _ in ast.walk(tree)) > 64:
            raise StructureError("invalid_symmetry_operation", "Affine expression is too complex.")

        def affine(node: ast.AST) -> tuple[Fraction, ...]:
            if isinstance(node, ast.Name) and node.id in "xyz" and len(node.id) == 1:
                return tuple(Fraction(int(node.id == name)) for name in "xyz") + (Fraction(0),)
            if isinstance(node, ast.Constant) and type(node.value) in (int, float):
                raw = ast.get_source_segment(component, node)
                if raw is None or len(raw) > 64:
                    raise StructureError("invalid_symmetry_operation", "Invalid affine numeric constant.")
                exponent = re.search(r"[eE]([+-]?\d+)$", raw)
                if exponent and (len(exponent[1]) > 4 or abs(int(exponent[1])) > 300):
                    raise StructureError("invalid_symmetry_operation", "Affine exponent is outside supported bounds.")
                try:
                    value = Fraction(raw)
                except (ValueError, ZeroDivisionError, OverflowError) as exc:
                    raise StructureError("invalid_symmetry_operation", "Invalid affine numeric constant.") from exc
                if value.numerator.bit_length() > 1024 or value.denominator.bit_length() > 1024:
                    raise StructureError("invalid_symmetry_operation", "Affine numeric constant is too large.")
                return (Fraction(0),) * 3 + (value,)
            if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
                sign = -1 if isinstance(node.op, ast.USub) else 1
                return tuple(sign * v for v in affine(node.operand))
            if isinstance(node, ast.BinOp):
                left, right = affine(node.left), affine(node.right)
                if isinstance(node.op, (ast.Add, ast.Sub)):
                    sign = -1 if isinstance(node.op, ast.Sub) else 1
                    return tuple(a + sign * b for a, b in zip(left, right))
                if isinstance(node.op, ast.Mult):
                    if not any(right[:3]):
                        return tuple(v * right[3] for v in left)
                    if not any(left[:3]):
                        return tuple(v * left[3] for v in right)
                if isinstance(node.op, ast.Div) and not any(right[:3]) and right[3]:
                    return tuple(v / right[3] for v in left)
            raise StructureError("nonaffine_symmetry_operation", "Only linear coordinates plus rational translations are permitted.")

        coefficients = affine(tree.body)
        rows.append(coefficients[:3])
        translation.append(coefficients[3])
    return SymmetryOperation(tuple(rows), tuple(translation), expression)
