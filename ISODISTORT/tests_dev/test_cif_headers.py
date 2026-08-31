"""Regression tests for Distortion CIF header metadata (HM / k / IR)."""
from __future__ import annotations

from isocore.backend import SubgroupInfo
from isocore.data.irreps_cdml import lookup_irrep_kovalev
from isocore.data.kpoints_official import format_k_point_display
from isocore.io.isodistort_cif import (
    _HM_FULL_1992,
    _irrep_kovalev_tag,
    _k_point_comment,
    _space_hm_symbol,
)


def test_space_hm_symbol_splits_mm2_and_1m1():
    assert _space_hm_symbol("Pmm2") == "P m m 2"
    assert _space_hm_symbol("P1m1") == "P 1 m 1"
    assert _space_hm_symbol("Imm2") == "I m m 2"
    assert _space_hm_symbol("I2/m2/m2/m") == "I 2/m 2/m 2/m"


def test_hm_full_table_covers_export_subgroups():
    assert _HM_FULL_1992[25] == "P m m 2"
    assert _HM_FULL_1992[6] == "P 1 m 1"
    assert _HM_FULL_1992[38] == "A m m 2"


def test_kpoint_display_p4mmm_and_i4mmm():
    assert format_k_point_display(123, "A") == "A, k20 (1/2,1/2,1/2)"
    assert format_k_point_display(123, "Y", ["1/3"]) == "Y, k6 (a,1/2,0), a=1/3"
    assert format_k_point_display(139, "LD", ["1/6"]) == "LD, k10 (0,0,g), g=1/6"


def test_irrep_tags_scoped_by_parent_sg():
    assert lookup_irrep_kovalev(
        "GM2+", k_kovalev="k14", parent_sg=139, k_point_label="GM"
    ) == "k14t5"
    assert lookup_irrep_kovalev(
        "M3-", k_kovalev="k15", parent_sg=139, k_point_label="M"
    ) == "k15t4"
    assert lookup_irrep_kovalev(
        "A1+", k_kovalev="k20", parent_sg=123, k_point_label="A"
    ) == "k20t1"
    assert lookup_irrep_kovalev(
        "Y2", k_kovalev="k6", parent_sg=123, k_point_label="Y"
    ) == "k6t2"
    assert lookup_irrep_kovalev(
        "LD2", k_kovalev="k10", parent_sg=139, k_point_label="LD"
    ) == "k10t3"


def test_cif_header_helpers_wire_parent_and_k():
    sg = SubgroupInfo(
        index=0,
        space_group_number=25,
        space_group_symbol="Pmm2",
        irrep_label="A1+",
        k_point_label="A",
        parent_sg=123,
        k_coordinates=["1/2", "1/2", "1/2"],
    )
    assert _k_point_comment(sg) == "k point: A, k20 (1/2,1/2,1/2)"
    assert _irrep_kovalev_tag(sg) == "k20t1"
