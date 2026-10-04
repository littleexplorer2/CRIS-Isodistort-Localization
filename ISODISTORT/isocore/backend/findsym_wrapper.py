"""
FINDSYM 封装 - 空间群识别与 Wyckoff 位置分析

对应阶段一，步骤2：母相空间群与 Wyckoff 位置识别

调用方式（重要）：findsym v6/v7 的关键字输入**必须**以输入文件参数方式调用
（``findsym inputfilename``），不能通过 stdin 传入。
"""
import ast
import re
import shlex
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np
from pymatgen.core import Structure

from ..utils import (
    OutputParseError,
    parse_space_group_number,
    parse_wyckoff_sites,
)
from .base_wrapper import BaseWrapper


@dataclass
class FindsymResult:
    """findsym 计算结果"""

    space_group_number: int          # 空间群号
    space_group_symbol: str = ""     # 空间群 Hermann-Mauguin 符号
    wyckoff_sites: list[dict] = field(default_factory=list)  # Wyckoff 位置信息
    lattice_params: tuple[float, ...] = ()  # 晶格参数 a,b,c,alpha,beta,gamma
    origin: tuple[float, float, float] = ()
    basis_vectors: tuple[tuple[float, float, float], ...] = ()
    standardized_sites: list[dict] = field(default_factory=list)
    raw_output: str = ""             # 原始输出（调试用）


class FindsymSettingMismatchError(OutputParseError):
    """FINDSYM found the same group in a different basis and/or origin."""

    def __init__(self, result: FindsymResult, message: str):
        self.result = result
        super().__init__("findsym", message)


class FindsymWrapper(BaseWrapper):
    """
    FINDSYM 程序封装

    功能：给定晶格参数与原子坐标，识别空间群、分配 Wyckoff 位置
    """

    def __init__(self) -> None:
        super().__init__()
        self.binary = self.cfg.findsym_bin

    def identify(self, lattice_params: list[float],
                 atom_types: list[str],
                 atom_positions: list[list[float]],
                 centering: str = "P",
                 title: str = "findsym input") -> FindsymResult:
        """
        调用 findsym 识别空间群

        Args:
            lattice_params: [a, b, c, alpha, beta, gamma] 晶格参数
            atom_types: 原子种类列表，如 ["Na", "Cl"]
            atom_positions: 原子分数坐标列表，与 atom_types 一一对应
            centering: 点阵中心类型 P/I/F/A/B/C/R
            title: 输入文件标题

        Returns:
            FindsymResult: 识别结果
        """
        if len(atom_types) != len(atom_positions):
            raise ValueError("atom_types 与 atom_positions 长度必须一致")
        if len(lattice_params) != 6:
            raise ValueError("lattice_params 必须为 [a,b,c,alpha,beta,gamma] 六元组")

        input_text = self._build_input(
            lattice_params, atom_types, atom_positions, centering, title
        )

        # findsym 关键字输入必须以输入文件参数方式调用
        stdout = self.run_input_file(self.binary, input_text)
        return self._parse_output(stdout)

    def standardize_wyckoff_orbits(
        self,
        structure: Structure,
        wyckoff_sites: list[dict],
        *,
        expected_space_group: int,
    ) -> list[dict]:
        """Return FINDSYM standard representatives for physical site orbits.

        Every inequivalent orbit receives its own temporary atom-type token.
        This keeps independent occurrences of one Wyckoff letter distinct while
        preserving all symmetry operations within each orbit.  The supplied
        ``Structure`` already contains every conventional-cell translation, so
        FINDSYM must not add centering translations a second time.

        Standard coordinates can only be combined with the current structure
        when FINDSYM keeps its basis and origin.  Report a typed setting
        mismatch so the session API can transform the complete structure and
        retry instead of silently mixing coordinate systems.
        """
        atom_types: list[str | None] = [None] * len(structure)
        for orbit_index, site in enumerate(wyckoff_sites, start=1):
            for atom_index in site.get("equivalent_indices") or []:
                index = int(atom_index)
                if index < 0 or index >= len(atom_types):
                    raise OutputParseError(
                        "findsym",
                        f"Wyckoff orbit contains invalid atom index {index}",
                    )
                if atom_types[index] is not None:
                    raise OutputParseError(
                        "findsym",
                        f"atom {index} belongs to more than one Wyckoff orbit",
                    )
                atom_types[index] = str(orbit_index)
        if any(value is None for value in atom_types):
            missing = [i for i, value in enumerate(atom_types) if value is None]
            raise OutputParseError(
                "findsym",
                f"Wyckoff orbits do not cover structure atoms {missing}",
            )

        lattice = structure.lattice
        result = self.identify(
            [
                float(lattice.a),
                float(lattice.b),
                float(lattice.c),
                float(lattice.alpha),
                float(lattice.beta),
                float(lattice.gamma),
            ],
            [str(value) for value in atom_types],
            [list(map(float, site.frac_coords)) for site in structure],
            centering="P",
            title="ISODISTORT parent Wyckoff orbit standardization",
        )
        if result.space_group_number != int(expected_space_group):
            raise OutputParseError(
                "findsym",
                "space-group disagreement while standardizing parent orbits: "
                f"spglib={expected_space_group}, FINDSYM={result.space_group_number}",
            )
        if (
            len(result.basis_vectors) != 3
            or any(len(vector) != 3 for vector in result.basis_vectors)
        ):
            raise OutputParseError(
                "findsym",
                "missing or invalid standard-setting basis transformation",
            )
        if len(result.origin) != 3:
            raise OutputParseError(
                "findsym",
                "missing or invalid standard-setting origin transformation",
            )
        if not np.allclose(
            np.asarray(result.basis_vectors, dtype=float),
            np.eye(3),
            atol=1e-8,
            rtol=0.0,
        ):
            raise FindsymSettingMismatchError(
                result,
                "parent CIF is not in the FINDSYM standard basis; automatic "
                "whole-structure setting conversion is required",
            )
        if not np.allclose(
            np.mod(np.asarray(result.origin, dtype=float), 1.0),
            np.zeros(3),
            atol=1e-8,
            rtol=0.0,
        ):
            raise FindsymSettingMismatchError(
                result,
                "parent CIF is not in the FINDSYM standard origin; automatic "
                "whole-structure setting conversion is required",
            )

        by_orbit: dict[int, dict] = {}
        for standard in result.standardized_sites:
            token = str(standard.get("atom_type") or "").strip()
            if not token.isdigit():
                raise OutputParseError(
                    "findsym", f"unexpected standardized orbit token {token!r}"
                )
            orbit_index = int(token) - 1
            if orbit_index in by_orbit:
                raise OutputParseError(
                    "findsym", f"duplicate standardized orbit token {token!r}"
                )
            by_orbit[orbit_index] = standard

        if set(by_orbit) != set(range(len(wyckoff_sites))):
            raise OutputParseError(
                "findsym",
                "standardized Wyckoff orbit set does not match the input orbit set",
            )

        metadata: list[dict] = []
        for orbit_index, site in enumerate(wyckoff_sites):
            standard = by_orbit[orbit_index]
            missing_fields = {
                "representative_frac_coords",
                "symmform",
                "parameters",
            } - standard.keys()
            if missing_fields:
                raise OutputParseError(
                    "findsym",
                    "incomplete standardized Wyckoff metadata for orbit "
                    f"{site.get('orbit_id', orbit_index)}: "
                    f"missing {sorted(missing_fields)}",
                )
            expected_letter = str(site.get("wyckoff_letter") or "")
            expected_multiplicity = int(site.get("multiplicity") or 0)
            if (
                standard["wyckoff_letter"] != expected_letter
                or int(standard["multiplicity"]) != expected_multiplicity
            ):
                raise OutputParseError(
                    "findsym",
                    "Wyckoff identity disagreement for orbit "
                    f"{site.get('orbit_id', orbit_index)}: "
                    f"spglib={expected_multiplicity}{expected_letter}, "
                    f"FINDSYM={standard['multiplicity']}"
                    f"{standard['wyckoff_letter']}",
                )
            metadata.append({
                "standard_representative_frac_coords": list(
                    standard["representative_frac_coords"]
                ),
                "standard_representative_symmform": standard["symmform"],
                "standard_representative_parameters": dict(
                    standard.get("parameters") or {}
                ),
                "standardization_source": "findsym",
            })
        return metadata

    # ---- 输入生成 ----

    @staticmethod
    def _build_input(lattice_params, atom_types, atom_positions,
                     centering, title) -> str:
        """生成 findsym 输入文件内容（关键字格式，见 isobyu/findsym.txt）"""
        a, b, c, alpha, beta, gamma = lattice_params
        n_atoms = len(atom_types)

        lines = [
            "!useKeyWords",
            "!title",
            title,
            "!latticeParameters",
            f"{a} {b} {c} {alpha} {beta} {gamma}",
            "!unitCellCentering",
            centering,
            "!atomCount",
            str(n_atoms),
            "!atomType",
            " ".join(atom_types),
            "!atomPosition",
        ]
        for pos in atom_positions:
            lines.append(f"{pos[0]} {pos[1]} {pos[2]}")

        return "\n".join(lines) + "\n"

    # ---- 输出解析 ----

    @staticmethod
    def _parse_output(text: str) -> FindsymResult:
        """解析 findsym 输出（兼容 v6 与 v7 格式）"""
        sg_num = parse_space_group_number(text)
        if sg_num is None:
            raise OutputParseError("findsym", "未找到空间群号")

        # 空间群符号：v6 "Space Group 225  Oh-5      Fm-3m"
        #            v7 "Space Group: 225  Oh-5      Fm-3m"
        sym_match = re.search(r"Space Group:?\s+\d+\s+(\S+)\s+(\S+)", text)
        sg_symbol = sym_match.group(2) if sym_match else ""

        # Wyckoff 位置（v7 中位点行带原子标签后缀，如 "Wyckoff position a (Na1)"）
        sites = parse_wyckoff_sites(text)

        # 晶格参数
        lat_match = re.search(
            r"Lattice parameters.*?:\s*\n\s*([\d.\s]+)", text
        )
        lattice_params: tuple[float, ...] = ()
        if lat_match:
            vals = [float(x) for x in lat_match.group(1).split()]
            if len(vals) == 6:
                lattice_params = tuple(vals)

        origin: tuple[float, float, float] = ()
        origin_match = re.search(
            r"^Origin at\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)",
            text,
            re.MULTILINE,
        )
        if origin_match:
            origin = tuple(float(origin_match.group(i)) for i in range(1, 4))

        basis_vectors: tuple[tuple[float, float, float], ...] = ()
        basis_match = re.search(
            r"^Vectors a,b,c:\s*\n"
            r"\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\n"
            r"\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\n"
            r"\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)",
            text,
            re.MULTILINE,
        )
        if basis_match:
            values = [float(basis_match.group(i)) for i in range(1, 10)]
            basis_vectors = tuple(
                tuple(values[start : start + 3]) for start in range(0, 9, 3)
            )

        standardized_sites = FindsymWrapper._parse_standardized_sites(text)
        # FINDSYM normally emits both sections in the same order, but that is
        # presentation detail rather than an identity contract.  Match the CIF
        # row label to the label printed in the Wyckoff section first; use a
        # same-letter, still-unmatched row only for older outputs without the
        # label column.  This is essential when one letter occurs repeatedly.
        unmatched = list(range(len(sites)))
        for standard in standardized_sites:
            output_label = str(standard.pop("_output_label", "") or "")
            match_index = next(
                (
                    index
                    for index in unmatched
                    if output_label
                    and str(sites[index].get("output_label") or "") == output_label
                ),
                None,
            )
            if match_index is None:
                match_index = next(
                    (
                        index
                        for index in unmatched
                        if sites[index].get("wyckoff_letter")
                        == standard["wyckoff_letter"]
                    ),
                    None,
                )
            if match_index is None:
                continue
            unmatched.remove(match_index)
            parsed = sites[match_index]
            standard["parameters"] = dict(parsed.get("parameters") or {})
            standard["input_atom_indices"] = [
                int(atom["index"]) - 1 for atom in parsed.get("atoms") or []
            ]
            standard["symmform"] = FindsymWrapper._absolute_symmform(
                standard["representative_frac_coords"],
                standard.pop("symmform_delta"),
                standard["parameters"],
            )

        return FindsymResult(
            space_group_number=sg_num,
            space_group_symbol=sg_symbol,
            wyckoff_sites=sites,
            lattice_params=lattice_params,
            origin=origin,
            basis_vectors=basis_vectors,
            standardized_sites=standardized_sites,
            raw_output=text,
        )

    @staticmethod
    def _parse_standardized_sites(text: str) -> list[dict]:
        """Parse the FINDSYM-generated CIF atom-site loop from stdout."""
        marker = "# CIF file created by FINDSYM"
        if marker not in text:
            return []
        lines = text.split(marker, 1)[1].splitlines()
        line_index = 0
        required = {
            "_atom_site_type_symbol",
            "_atom_site_symmetry_multiplicity",
            "_atom_site_wyckoff_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_fract_symmform",
        }
        while line_index < len(lines):
            if lines[line_index].strip().lower() != "loop_":
                line_index += 1
                continue
            line_index += 1
            tags: list[str] = []
            while (
                line_index < len(lines)
                and lines[line_index].strip().startswith("_")
            ):
                tags.append(lines[line_index].strip().split(maxsplit=1)[0].lower())
                line_index += 1
            if not required.issubset(tags):
                continue
            tokens: list[str] = []
            while line_index < len(lines):
                stripped = lines[line_index].strip()
                lowered = stripped.lower()
                if (
                    lowered == "loop_"
                    or lowered.startswith("data_")
                    or stripped.startswith("_")
                    or stripped.startswith("#")
                ):
                    break
                line_index += 1
                if not stripped:
                    continue
                tokens.extend(shlex.split(stripped, comments=True, posix=True))
            width = len(tags)
            index = {tag: i for i, tag in enumerate(tags)}
            rows: list[dict] = []
            for start in range(0, len(tokens) - width + 1, width):
                row = tokens[start : start + width]
                try:
                    symmform = row[index["_atom_site_fract_symmform"]]
                    rows.append({
                        "atom_type": row[index["_atom_site_type_symbol"]],
                        "_output_label": (
                            row[index["_atom_site_label"]]
                            if "_atom_site_label" in index
                            else ""
                        ),
                        "multiplicity": int(
                            row[index["_atom_site_symmetry_multiplicity"]]
                        ),
                        "wyckoff_letter": row[
                            index["_atom_site_wyckoff_symbol"]
                        ],
                        "representative_frac_coords": [
                            float(row[index["_atom_site_fract_x"]]),
                            float(row[index["_atom_site_fract_y"]]),
                            float(row[index["_atom_site_fract_z"]]),
                        ],
                        "symmform_delta": symmform,
                        "parameters": {},
                    })
                except (KeyError, ValueError, IndexError):
                    continue
            return rows
        return []

    @staticmethod
    def _absolute_symmform(
        representative: list[float],
        delta_form: str,
        parameters: dict[str, float],
    ) -> str:
        """Combine FINDSYM's displacement form with its standard base point.

        ``_atom_site_fract_symmform`` stores allowed coordinate changes, not
        absolute Wyckoff coordinates.  For example, the #139 ``8g`` site is
        reported as numeric ``(0,1/2,z)`` with delta form ``0,0,Dz``.
        """
        pieces = [piece.strip() for piece in delta_form.split(",")]
        if len(pieces) != 3:
            raise OutputParseError(
                "findsym", f"invalid _atom_site_fract_symmform {delta_form!r}"
            )
        absolute: list[str] = []
        for coordinate, delta in zip(representative, pieces, strict=True):
            if not re.search(r"D[xyz]", delta):
                absolute.append(FindsymWrapper._format_fraction(coordinate))
                continue
            expression = re.sub(r"D([xyz])", r"\1", delta)
            delta_value = FindsymWrapper._eval_symmform(expression, parameters)
            offset = (float(coordinate) - delta_value) % 1.0
            if np.isclose(offset, 1.0, atol=1e-8) or np.isclose(
                offset, 0.0, atol=1e-8
            ):
                offset = 0.0
            if offset:
                expression = f"{expression}+{FindsymWrapper._format_fraction(offset)}"
            absolute.append(expression)
        return ",".join(absolute)

    @staticmethod
    def _eval_symmform(expression: str, parameters: dict[str, float]) -> float:
        """Evaluate FINDSYM's small linear coordinate grammar safely."""
        # FINDSYM writes conventional crystallographic products such as
        # ``2Dx`` and ``1/2Dx`` without an explicit multiplication sign.
        python_expression = re.sub(r"(?<=\d)(?=[xyz])", "*", expression)
        try:
            tree = ast.parse(python_expression, mode="eval")
        except SyntaxError as exc:
            raise OutputParseError(
                "findsym", f"invalid symmetry form expression {expression!r}"
            ) from exc

        def evaluate(node: ast.AST) -> float:
            if isinstance(node, ast.Expression):
                return evaluate(node.body)
            if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
                return float(node.value)
            if isinstance(node, ast.Name) and node.id in {"x", "y", "z"}:
                if node.id not in parameters:
                    raise OutputParseError(
                        "findsym",
                        f"missing {node.id} parameter for symmetry form {expression!r}",
                    )
                return float(parameters[node.id])
            if isinstance(node, ast.UnaryOp) and isinstance(
                node.op, (ast.UAdd, ast.USub)
            ):
                value = evaluate(node.operand)
                return value if isinstance(node.op, ast.UAdd) else -value
            if isinstance(node, ast.BinOp) and isinstance(
                node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)
            ):
                left = evaluate(node.left)
                right = evaluate(node.right)
                if isinstance(node.op, ast.Add):
                    return left + right
                if isinstance(node.op, ast.Sub):
                    return left - right
                if isinstance(node.op, ast.Mult):
                    return left * right
                if right == 0.0:
                    raise OutputParseError(
                        "findsym",
                        f"division by zero in symmetry form {expression!r}",
                    )
                return left / right
            raise OutputParseError(
                "findsym", f"unsupported symmetry form expression {expression!r}"
            )

        return evaluate(tree)

    @staticmethod
    def _format_fraction(value: float) -> str:
        wrapped = float(value) % 1.0
        if np.isclose(wrapped, 1.0, atol=1e-8):
            wrapped = 0.0
        fraction = Fraction(wrapped).limit_denominator(48)
        if not np.isclose(float(fraction), wrapped, atol=1e-8):
            return f"{wrapped:.10g}"
        if fraction.denominator == 1:
            return str(fraction.numerator)
        return f"{fraction.numerator}/{fraction.denominator}"
