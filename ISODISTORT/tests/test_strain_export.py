from __future__ import annotations

import math
import re
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from backend.api.core_api import IsoDistort
from backend.models.iso_mode_models import (
    InvariantDirection,
    MacroscopicTensorBasis,
    MacroscopicTensorBlock,
    SymbolicTensorComponent,
)
from backend.utils import get_config
from backend.wrappers import SubgroupInfo
from features.export.distortion_formats import (
    StrainExportData,
    SubgroupExportSpec,
    render_cif,
    render_complete_modes,
    render_isoviz,
    render_topas,
    resolve_strained_export_structure,
)
from features.input_cif import read_cif
from features.method4.strain_modes import (
    CanonicalStrainModeDefinition,
    apply_canonical_strain_basis,
    compute_homogeneous_strain_modes,
    engineering_voigt_to_tensor,
)

_REPO = Path(__file__).resolve().parents[2]
_F02 = _REPO / "output_compare" / "EuAl4 Parent.cif" / "官网" / "Method4" / "F02"


def _subgroup(number: int, symbol: str, basis: list[list[float]] | None = None) -> SubgroupInfo:
    return SubgroupInfo(
        index=0,
        parent_sg=number,
        space_group_number=number,
        space_group_symbol=symbol,
        subgroup_index=1,
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_vectors=basis or np.eye(3).tolist(),
        origin=[0.0, 0.0, 0.0],
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM_TEST",
    )


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
        subgroup=_subgroup(1, "P1"),
        structure=structure,
        parent_structure=structure,
        parent_sg=1,
        parent_symbol="P1",
        distortion_types=["strain"],
        strain_data=StrainExportData(result=result, amplitudes=values),
    )


def _bare_api(structure: Structure, distortion_types: list[str]) -> IsoDistort:
    """Construct only the pure export context; no WSL backend is needed."""

    api = object.__new__(IsoDistort)
    api.cfg = get_config()
    api.structure = structure
    api.distortion_types = distortion_types
    api._strain_representation_cache = None
    return api


def _macro_component(expression: str, values: list[int]) -> SymbolicTensorComponent:
    return SymbolicTensorComponent(
        expression_raw=expression,
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        coefficients=tuple(Fraction(value) for value in values),
        quality="exact_rational",
    )


def _tetragonal_macro_basis() -> MacroscopicTensorBasis:
    return MacroscopicTensorBasis(
        rank_signature="[12]",
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        blocks=(
            MacroscopicTensorBlock(
                "GM1+", (_macro_component("xx+yy", [1, 1, 0, 0, 0, 0]),), 1, 1
            ),
            MacroscopicTensorBlock(
                "GM1+", (_macro_component("zz", [0, 0, 1, 0, 0, 0]),), 2, 2
            ),
            MacroscopicTensorBlock(
                "GM2+", (_macro_component("xx-yy", [1, -1, 0, 0, 0, 0]),), 3
            ),
            MacroscopicTensorBlock(
                "GM4+", (_macro_component("xy", [0, 0, 0, 0, 0, 1]),), 4
            ),
            MacroscopicTensorBlock(
                "GM5+",
                (
                    _macro_component("yz", [0, 0, 0, 1, 0, 0]),
                    _macro_component("-xz", [0, 0, 0, 0, -1, 0]),
                ),
                5,
            ),
        ),
        status="verified",
    )


class _TetragonalMacroIso:
    def get_macroscopic_tensor_basis(self, _parent_sg: int) -> MacroscopicTensorBasis:
        return _tetragonal_macro_basis()

    def list_invariant_directions(
        self, _parent_sg: int, subgroup: SubgroupInfo
    ) -> list[InvariantDirection]:
        if int(subgroup.space_group_number) not in {1, 2}:
            raise AssertionError("test ISO stub is only valid for P1/P-1 embeddings")
        return [
            InvariantDirection(label, direction, 2, "P-1", subgroup.size)
            for label, direction in (
                ("GM1+", "(a)"),
                ("GM2+", "(a)"),
                ("GM4+", "(a)"),
                ("GM5+", "(a,b)"),
            )
        ]


class _FlakyMacroIso(_TetragonalMacroIso):
    def __init__(self) -> None:
        self.macroscopic_calls = 0

    def get_macroscopic_tensor_basis(self, _parent_sg: int) -> MacroscopicTensorBasis:
        self.macroscopic_calls += 1
        if self.macroscopic_calls == 1:
            return MacroscopicTensorBasis.unresolved(
                rank_signature="[12]",
                coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
                reason="transient_test_backend_failure",
            )
        return _tetragonal_macro_basis()


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


@pytest.mark.parametrize(
    ("symbol", "lattice", "expected"),
    (
        ("P-1", Lattice.from_parameters(4.1, 5.2, 6.3, 76.0, 83.0, 69.0), 6),
        ("P2/m", Lattice.monoclinic(4.1, 5.2, 6.3, 107.0), 4),
        ("Pmmm", Lattice.orthorhombic(4.1, 5.2, 6.3), 3),
        ("P4/mmm", Lattice.tetragonal(4.1, 6.3), 2),
        ("P6/mmm", Lattice.hexagonal(4.1, 6.3), 2),
        ("Pm-3m", Lattice.cubic(4.1), 1),
    ),
)
def test_metric_fixed_space_gives_general_crystal_system_counts(
    symbol: str, lattice: Lattice, expected: int
) -> None:
    structure = Structure.from_spacegroup(symbol, lattice, ["H"], [[0.123, 0.234, 0.345]])
    rotations = [
        operation.rotation_matrix
        for operation in SpacegroupAnalyzer(
            structure, symprec=1.0e-3
        ).get_symmetry_operations(cartesian=False)
    ]
    result = compute_homogeneous_strain_modes(
        rotations, structure.lattice.matrix
    )

    assert len(result.modes) == expected
    assert result.diagnostics.point_group_order >= 2


def test_4310_n1_plus_4d1_has_six_embedded_strain_modes_in_strain_exports() -> None:
    parent = read_cif(_REPO / "experiment_data" / "4310_tetra.cif")
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

    parent = read_cif(_REPO / "experiment_data" / "EuAl4 Parent.cif")
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

    assert api._strain_export_data_for_subgroup(_subgroup(221, "Pm-3m")) is None


def test_contract_rejects_strain_data_without_strain_type() -> None:
    valid = _p1_export_spec()
    with pytest.raises(ValueError, match="distortion_types includes 'strain'"):
        SubgroupExportSpec(
            subgroup=valid.subgroup,
            structure=valid.structure,
            distortion_types=["displacive"],
            strain_data=valid.strain_data,
        )


def test_unresolved_macro_basis_is_not_cached_and_next_attempt_recovers() -> None:
    structure = Structure(
        Lattice.from_parameters(4.0, 5.0, 6.0, 78.0, 83.0, 71.0),
        ["H"],
        [[0.13, 0.27, 0.39]],
    )
    api = _bare_api(structure, ["strain"])
    backend = _FlakyMacroIso()
    api._iso = backend
    subgroup = _subgroup(1, "P1")

    with pytest.raises(ValueError, match="macroscopic strain basis is unresolved"):
        api._strain_export_data_for_subgroup(subgroup)

    recovered = api._strain_export_data_for_subgroup(subgroup)

    assert recovered is not None
    assert recovered.result.export_ready
    assert backend.macroscopic_calls == 2


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

    from tests.test_displacive_export import (
        _export_data,
        _spec,
    )

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
