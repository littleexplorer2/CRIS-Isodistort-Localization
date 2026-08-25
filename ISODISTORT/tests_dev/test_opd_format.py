"""Official Method 1 order-parameter line format and subgroup-table raw tokens."""
from __future__ import annotations

import re
from pathlib import Path

from isocore.backend import SubgroupInfo
from isocore.utils import format_opd_line, parse_subgroup_table
from isocore.utils.opd_format import format_k_active, k_star_tuples

_REPO = Path(__file__).resolve().parents[2]
_OFFICIAL_HTML = (
    _REPO / "webpage_info" / "a. ISODISTORT_ order parameter direction.html"
)


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
