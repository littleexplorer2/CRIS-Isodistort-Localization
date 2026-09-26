"""OPD line format and Distortion Page batch export formats."""
from __future__ import annotations

import io
import re
import zipfile
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from isocore.api.core_api import IsoDistort
from isocore.backend import SubgroupInfo
from isocore.io.distortion_formats import (
    FORMAT_CIF,
    FORMAT_ISOVIZ,
    FORMAT_MODES,
    FORMAT_TOPAS,
    SubgroupExportSpec,
    build_export_zip,
    cart_normalized_mode_matrix,
    format_filename,
    parse_export_formats,
    parse_export_method,
    render_cif,
    render_complete_modes,
    render_isoviz,
    render_topas,
    subgroup_label,
    unique_folder_name,
)
from isocore.utils import format_opd_line, parse_subgroup_table
from isocore.utils.opd_format import format_k_active, format_number, k_star_tuples

_REPO = Path(__file__).resolve().parents[2]
_OFFICIAL_HTML = (
    _REPO / "webpage_info" / "a. ISODISTORT_ order parameter direction.html"
)


# --- from test_opd_format.py ---


def test_format_number_preserves_exact_fractions_without_rounding_to_zero():
    assert format_number(Fraction(1, 97)) == "1/97"
    assert format_number(Fraction(2_000_001, 2_000_000)) == "2000001/2000000"
    assert format_number("1/2000001") == "1/2000001"
    assert format_number(0.25) == "1/4"
    approximate = format_number(1 / 97)
    assert approximate == repr(1 / 97)
    assert approximate != "0"
    assert format_number(0.5000000001) == repr(0.5000000001)

def test_parse_subgroup_table_keeps_letter_opd_and_raw_basis():
    text = (
        "Subgroup  Max Index Size Dir         Basis Vectors                Origin\n"
        "139 I4/mmm yes 1     1    P1  (a) (1,0,0),(0,1,0),(0,0,1) (0,0,0)\n"
        "123 P4/mmm yes 4     4    P1  (a;a) (1,-1,0),(1,1,0),(0,0,1) (0,0,0)\n"
        " 12 C2/m  no  4     1    P1  (a,0) (0,1,1),(-1,0,0),(0,-1,0) (1/2,-1/2,0)\n"
    )
    rows = parse_subgroup_table(text)
    assert len(rows) == 3
    assert rows[0]["opd_dir_raw"] == "(a)"
    assert rows[0]["basis_raw"] == "(1,0,0),(0,1,0),(0,0,1)"
    assert rows[0]["origin_raw"] == "(0,0,0)"
    assert rows[1]["opd_dir_raw"] == "(a;a)"
    assert rows[2]["origin_raw"] == "(1/2,-1/2,0)"
    assert rows[2]["origin"] == [0.5, -0.5, 0.0]


def test_format_opd_line_gm1_matches_official_visible_text():
    line = format_opd_line(
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=1,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert line == (
        "GM1+     P1   (a) 139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
        "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
    )


def test_format_opd_line_pads_sg_number_and_x_point():
    line = format_opd_line(
        irrep_label="GM2+",
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
    )
    assert line.startswith("GM2+     P1   (a)  71 Immm")


def test_x_star_and_k_active_two_arms():
    star = k_star_tuples(["1/2", "1/2", "0"], 139)
    assert star[0] == "(1/2,1/2,0)"
    assert len(star) == 2
    assert "(1/2,1/2,1)" in star
    both = format_k_active("(a;a)", ["1/2", "1/2", "0"], 139)
    one = format_k_active("(a;0)", ["1/2", "1/2", "0"], 139)
    assert both == " (1/2,1/2,0),(1/2,1/2,1)"
    assert one == " (1/2,1/2,0)"


def test_n_star_four_arms_partial_active():
    star = k_star_tuples(["1/2", "0", "1/2"], 139)
    assert star == [
        "(1/2,0,1/2)",
        "(1/2,1,1/2)",
        "(0,1/2,1/2)",
        "(1,1/2,1/2)",
    ]
    assert format_k_active("(a;a;0;0)", ["1/2", "0", "1/2"], 139) == (
        " (1/2,0,1/2),(1/2,1,1/2)"
    )
    assert format_k_active("(a;0;a;0)", ["1/2", "0", "1/2"], 139) == (
        " (1/2,0,1/2),(0,1/2,1/2)"
    )


def test_m_point_k_active_keeps_cdml_coords():
    star = k_star_tuples(["1", "1", "1"], 139)
    assert star == ["(1,1,1)"]
    assert format_k_active("(a)", ["1", "1", "1"], 139) == " (1,1,1)"


def test_y_star_keeps_exact_thirds_after_reciprocal_rotation():
    star = k_star_tuples(["1/3", "1/2", "0"], 123)

    assert star[:2] == ["(1/3,1/2,0)", "(1/2,1/3,0)"]
    assert len(star) == 4
    assert all("0.333333" not in arm for arm in star)
    assert format_k_active("(a;a)", ["1/3", "1/2", "0"], 123) == (
        " (1/3,1/2,0),(1/2,1/3,0)"
    )


def test_strain_only_gamma_keeps_gm4_not_gm3():
    api = IsoDistort()
    api.structure = Structure.from_spacegroup(
        "I4/mmm", Lattice.tetragonal(4.0, 10.0), ["Eu"], [[0, 0, 0]],
    )
    gm4 = SimpleNamespace(subgroup=SubgroupInfo(
        index=0,
        parent_sg=139,
        k_point_label="GM",
        irrep_label="GM4+",
        size=1,
        subgroup_index=2,
        space_group_number=69,
        space_group_symbol="Fmmm",
        basis_vectors=[[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
    ))
    gm3 = SimpleNamespace(subgroup=SubgroupInfo(
        index=1,
        parent_sg=139,
        k_point_label="GM",
        irrep_label="GM3+",
        size=1,
        subgroup_index=2,
        space_group_number=87,
        space_group_symbol="I4/m",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
    ))
    mirror_only = SimpleNamespace(subgroup=SubgroupInfo(
        index=2,
        parent_sg=139,
        k_point_label="GM",
        # Synthetic even label isolates the representation test: Cs fixes
        # more tensors than D4h, but inversion acts trivially on every strain
        # tensor, so its exact pointwise stabilizer is strictly larger than Cs.
        irrep_label="GM_FAKE+",
        size=1,
        subgroup_index=8,
        space_group_number=8,
        space_group_symbol="Cm",
        basis_vectors=[[-1, 1, 0], [0, 0, 1], [1, 0, 0]],
        origin=[0, 0, 0],
    ))
    assert api._keep_strain_only_irrep(gm4, 139)
    assert not api._keep_strain_only_irrep(gm3, 139)
    assert not api._keep_strain_only_irrep(mirror_only, 139)


def test_subgroup_from_row_uses_iso_basis_raw():
    """Display tokens come from iso output, not a memorized official OPD table."""
    from isocore.backend.iso_wrapper import IsoWrapper

    row = {
        "space_group_number": 12,
        "space_group_symbol": "C2/m",
        "is_maximal": True,
        "subgroup_index": 4,
        "size": 1,
        "opd_symbol": "P1",
        "opd_vector": [1.0, 0.0],
        "opd_dir_raw": "(a,0)",
        "basis_vectors": [[0, 1, 1], [1, 0, 0], [0, 0, -1]],
        "basis_raw": "(0,1,1),(1,0,0),(0,0,-1)",
        "origin": [0.0, 0.0, 0.0],
        "origin_raw": "(0,0,0)",
    }
    sg = IsoWrapper._subgroup_from_row(
        row, index=0, k_point="GM", irrep_label="GM5+",
        k_coordinates=["0", "0", "0"], parent_sg=139,
    )
    assert sg.basis_raw == "(0,1,1),(1,0,0),(0,0,-1)"
    assert "basis={(0,1,1),(1,0,0),(0,0,-1)}" in sg.opd_line()


def test_subgroupinfo_opd_line_no_maximal_asterisk():
    sg = SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        subgroup_index=1,
        size=1,
        is_maximal=True,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        irrep_label="GM1+",
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    line = sg.opd_line()
    assert "*" not in line
    assert "k-active= (0,0,0)" in line


def test_official_html_radio_lines_parseable():
    if not _OFFICIAL_HTML.is_file():
        return
    html = _OFFICIAL_HTML.read_text(encoding="utf-8", errors="replace")
    lines = re.findall(r'name="orderparam"[^>]*>([^<]+)<br>', html)
    assert len(lines) >= 50
    assert lines[0].startswith("GM1+     P1   (a) 139 I4/mmm")
    assert "basis={(1,0,0),(0,1,0),(0,0,1)}" in lines[0]
    assert "k-active=" in lines[0]
    formatted = format_opd_line(
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=1,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert formatted == lines[0].strip()
    x_p1 = next(ln for ln in lines if ln.startswith("X1+") and "(a;a)" in ln)
    assert "k-active= (1/2,1/2,0),(1/2,1/2,1)" in x_p1
    gm4 = next(ln for ln in lines if ln.startswith("GM4+"))
    assert "69 Fmmm" in gm4
    assert "basis={(1,1,0),(-1,1,0),(0,0,1)}" in gm4
    gm5_bases = [
        ln for ln in lines if ln.startswith("GM5+")
    ]
    assert any("{(0,1,1),(-1,0,0),(0,-1,0)}" in ln for ln in gm5_bases)
    assert any("{(1,-1,0),(1,1,0),(-1/2,1/2,1/2)}" in ln for ln in gm5_bases)
    assert any("{(1,0,0),(0,1,0),(-1/2,-1/2,1/2)}" in ln for ln in gm5_bases)
    m1 = next(ln for ln in lines if ln.startswith("M1+"))
    assert "k-active= (1,1,1)" in m1


def test_official_method1_fields_match_radio_tokens():
    from isocore.utils.opd_format import official_method1_fields

    fields = official_method1_fields(
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=1,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert fields["irrep"] == "GM1+"
    assert fields["opd"] == "P1"
    assert fields["dir"] == "(a)"
    assert fields["sg"] == "139 I4/mmm"
    assert fields["basis"] == "{(1,0,0),(0,1,0),(0,0,1)}"
    assert fields["origin"] == "(0,0,0)"
    assert fields["s"] == 1
    assert fields["i"] == 1
    assert fields["k_active"] == "(0,0,0)"
    sg = SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        subgroup_index=1,
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        irrep_label="GM1+",
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert sg.official_fields()["k_active"] == "(0,0,0)"
    assert "k-active= (0,0,0)" in sg.opd_line()


# --- from test_distortion_export.py ---

def _sg(irrep: str, opd: str, index: int = 0, number: int = 123,
        symbol: str = "P4/mmm") -> SubgroupInfo:
    return SubgroupInfo(
        index=index,
        space_group_number=number,
        space_group_symbol=symbol,
        subgroup_index=12,
        size=12,
        opd_symbol=opd,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
        origin=[0, 0, 0],
        k_point_label="LD",
        irrep_label=irrep,
        k_parameters=["1/6"],
    )


def _cubic() -> Structure:
    return Structure(
        [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
        ["Eu", "Al"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )


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


def test_safe_name_deletes_slash_like_official_windows_download():
    """官网 Windows 下载删除 ``/``：``I4/mmm``→``I4mmm``，``1/2``→``12``。"""
    from isocore.io.distortion_formats import safe_name

    assert safe_name("I4/mmm") == "I4mmm"
    assert safe_name("C2/m") == "C2m"
    assert (
        safe_name(
            "GM5+ C1 (a,b) 2 P-1, basis={(1,0,0),(0,1,0),(-1/2,-1/2,1/2)}, "
            "origin=(0,0,0), s=1, i=8, k-active= (0,0,0)"
        )
        == (
            "GM5+ C1 (a,b) 2 P-1, basis={(1,0,0),(0,1,0),(-12,-12,12)}, "
            "origin=(0,0,0), s=1, i=8, k-active= (0,0,0)"
        )
    )


def test_opd_line_body_omits_irrep_for_cif_comment():
    from isocore.utils.opd_format import format_opd_line_body

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


def test_unique_folder_name_disambiguates():
    used: set[str] = set()
    a = unique_folder_name(_sg("LD1", "C1", index=0, symbol="P4/mmm"), used)
    b = unique_folder_name(_sg("LD1", "C1", index=1, symbol="P4mm"), used)
    assert a == "LD1 C1"
    assert b == "LD1 C1 P4mm"
    assert a != b


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
            folder_name="LD1 C1",
        ),
        SubgroupExportSpec(
            subgroup=_sg("LD5", "P6", index=1, number=11, symbol="P2_1/m"),
            structure=parent,
            parent_structure=parent,
            parent_sg=139,
            parent_symbol="I4/mmm",
            folder_name="LD5 P6",
        ),
    ]
    raw = build_export_zip(
        specs,
        [FORMAT_CIF, FORMAT_ISOVIZ, FORMAT_MODES, FORMAT_TOPAS],
    )
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()

    # 各子群文件夹在 ZIP 根下（官网同款）
    assert any(n.startswith("LD1 C1/") for n in names)
    assert any(n.startswith("LD5 P6/") for n in names)
    assert "LD1 C1/subgroup.cif" in names
    assert "LD1 C1/data.isoviz" in names
    assert "LD1 C1/Complete modes details.txt" in names
    assert "LD1 C1/topas.str" in names
    assert "LD5 P6/subgroup.cif" in names
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
        folder_name="LD1 C1",
    )
    raw = build_export_zip([spec], [FORMAT_CIF])
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = zf.namelist()
    assert names == ["LD1 C1/subgroup.cif"]


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
    sg = SubgroupInfo(
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
    from isocore.data.kpoints_official import official_special_k_coords

    coords = official_special_k_coords(139, "LD", [], ["1/6"])
    assert coords == ["0", "0", "1/6"]


def test_subgroup_cif_origin_preserves_atom_count(tmp_path):
    """The subgroup's origin must make its symmetry operations valid."""
    from pymatgen.io.cif import CifParser

    from isocore.structure import read_cif

    parent = read_cif(_REPO / "experiment_data" / "EuAl4 Parent.cif")
    subgroup = SubgroupInfo(
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
    from isocore.structure import SymmetryValidator, read_cif
    from isocore.structure.coordinate_transform import build_supercell
    from isocore.io.isodistort_cif import _apply_origin_choice, _setting
    from isocore.utils.parent_header import parent_wyckoff_display

    path = _REPO / "experiment_data" / "NdNiO2 own.cif"
    parent = read_cif(path)
    basis = [[1, 1, 0], [-1, 1, 0], [0, 0, 2]]
    child = build_supercell(parent, basis)
    subgroup = SubgroupInfo(
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
    subgroup = SubgroupInfo(
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
    subgroup = SubgroupInfo(
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
    subgroup = SubgroupInfo(
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
    from isocore.structure import SymmetryValidator, read_cif
    from isocore.utils.parent_header import parent_wyckoff_display

    path = _REPO / "experiment_data" / "NdNiO2 own.cif"
    parent = read_cif(path)
    sites = SymmetryValidator().validate(parent)["wyckoff_sites"]
    lines = parent_wyckoff_display(parent, sites, path)
    assert [line.split()[0] for line in lines] == ["O", "ND", "NI"]
    subgroup = SubgroupInfo(
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
    from isocore.structure import read_cif
    from isocore.structure.coordinate_transform import build_supercell

    parent = read_cif(_REPO / "experiment_data" / "NdNiO2 own.cif")
    basis = [[-3, 0, 0], [0, 0, 1], [0, 2, 0]]
    child = build_supercell(parent, basis)
    subgroup = SubgroupInfo(
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

    official = _REPO / "ISODISTORT_VALIDATE" / "compare" / "true" / "X4- P3 CIF.cif"
    local_sc = _REPO / "ISODISTORT_VALIDATE" / "compare" / "item" / "X4- P3 CIF.cif"
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
    sg = SubgroupInfo(
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
