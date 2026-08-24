"""Distortion Page 批量导出：命名、ZIP 结构、四种官网格式、不混入无关文件。"""
from __future__ import annotations

import io
import zipfile

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
