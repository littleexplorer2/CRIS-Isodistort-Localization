"""OPD line format and Distortion Page batch export formats."""
from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Structure

from isocore.backend import SubgroupInfo
from isocore.io.distortion_formats import (
    FORMAT_CIF,
    FORMAT_ISOVIZ,
    FORMAT_MODES,
    FORMAT_TOPAS,
    SubgroupExportSpec,
    build_export_zip,
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
from isocore.utils.opd_format import format_k_active, k_star_tuples

_REPO = Path(__file__).resolve().parents[2]
_OFFICIAL_HTML = (
    _REPO / "webpage_info" / "a. ISODISTORT_ order parameter direction.html"
)


# --- from test_opd_format.py ---

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
    assert len(star) == 4
    kact = format_k_active("(a;a;0;0)", ["1/2", "0", "1/2"], 139)
    assert kact.startswith(" ")
    assert kact.count("(") == 2


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
    assert format_filename("LD1 C1", FORMAT_CIF) == "LD1 C1 CIF.cif"
    assert format_filename("LD1 C1", FORMAT_ISOVIZ) == (
        "LD1 C1 Save interactive distortion.isoviz"
    )
    assert format_filename("LD1 C1", FORMAT_MODES) == (
        "LD1 C1 Complete modes details.txt"
    )
    assert format_filename("LD1 C1", FORMAT_TOPAS) == "LD1 C1 TOPAS.STR"


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

    assert all(n.startswith("isodistort_outputs/") for n in names)
    # 每个子群一个文件夹
    assert any(n.startswith("isodistort_outputs/LD1 C1/") for n in names)
    assert any(n.startswith("isodistort_outputs/LD5 P6/") for n in names)
    assert "isodistort_outputs/LD1 C1/LD1 C1 CIF.cif" in names
    assert "isodistort_outputs/LD1 C1/LD1 C1 Save interactive distortion.isoviz" in names
    assert "isodistort_outputs/LD1 C1/LD1 C1 Complete modes details.txt" in names
    assert "isodistort_outputs/LD1 C1/LD1 C1 TOPAS.STR" in names
    assert "isodistort_outputs/LD5 P6/LD5 P6 CIF.cif" in names
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
    assert names == ["isodistort_outputs/LD1 C1/LD1 C1 CIF.cif"]


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
    assert "data_" in cif or "_cell_length_a" in cif

    isoviz = render_isoviz(spec)
    assert "origintype=isovizdistortion" in isoviz
    assert "atomicradius=" in isoviz
    assert "BEGIN MODES" in isoviz
    assert "LD1" in isoviz

    details = render_complete_modes(spec)
    assert "Complete modes details" in details
    assert "normfactor" in details
    assert "As=" in details and "Ap=" in details

    topas = render_topas(spec)
    assert topas.startswith("'Topas .str file generated by ISODISTORT")
    assert "space_group 123" in topas
    assert "prm !a1" in topas
    assert "mode-amplitude to delta" in topas

