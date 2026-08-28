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


def test_strain_only_gamma_keeps_gm4_not_gm3():
    from types import SimpleNamespace

    from isocore.api.core_api import IsoDistort

    gm4 = SimpleNamespace(subgroup=SimpleNamespace(
        k_point_label="GM", irrep_label="GM4+", size=1, space_group_number=69,
    ))
    gm3 = SimpleNamespace(subgroup=SimpleNamespace(
        k_point_label="GM", irrep_label="GM3+", size=1, space_group_number=87,
    ))
    assert IsoDistort._keep_strain_only_irrep(gm4, 139)
    assert not IsoDistort._keep_strain_only_irrep(gm3, 139)


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


