"""Tests for backend.tables; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.export.isodistort_cif import _HM_FULL_1992
from backend.tables.irreps_cdml import lookup_irrep_kovalev
from backend.tables.kpoints_official import format_k_point_display
from backend.tables.kpoints_official import official_kparams_to_iso
from backend.tables.kpoints_official import official_special_k_coords
from backend.wrappers import KPointInfo
from fractions import Fraction
import pytest


# Original coverage: test_basics.py.

# --- from test_kparams_official.py ---

def test_ld_g_to_iso():
    iso_kp = KPointInfo(
        label="LD", coordinates=["0", "0", "2a"], parameters=["2a"], is_special=False
    )
    out = official_kparams_to_iso(139, "LD", ["1/6"], iso_kp)
    assert out == ["1/12"]

def test_sm_a_to_iso():
    iso_kp = KPointInfo(
        label="SM", coordinates=["2a", "0", "0"], parameters=["2a"], is_special=False
    )
    out = official_kparams_to_iso(139, "SM", ["1/4"], iso_kp)
    assert out == ["1/8"]

def test_dt_unchanged():
    iso_kp = KPointInfo(
        label="DT", coordinates=["a", "a", "0"], parameters=["a"], is_special=False
    )
    out = official_kparams_to_iso(139, "DT", ["1/4"], iso_kp)
    assert out == ["1/4"]

def test_no_override_passthrough():
    iso_kp = KPointInfo(
        label="Z", coordinates=["0", "0", "1/2"], parameters=[], is_special=True
    )
    out = official_kparams_to_iso(221, "Z", ["1/6"], iso_kp)
    assert out == ["1/6"]

# Original coverage: test_cif_headers.py.

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

# Original coverage: test_theory_contracts.py.

def test_affine_wavevector_substitution_keeps_offsets_and_signs():
    assert official_special_k_coords(139, "Y", [], ["1/3"]) == ["1/3", "2/3", "0"]

def test_kvalue_conversion_solves_the_coordinate_equation():
    # Both affine coordinates, rather than their textual coefficients alone,
    # must denote the same reciprocal vector.
    kp = KPointInfo(label="LD", coordinates=["0", "0", "1-2a"],
                   parameters=["2a"], is_special=False)
    assert official_kparams_to_iso(139, "LD", ["1/6"], kp) == ["5/12"]

def test_one_free_variable_is_not_counted_as_two_coordinate_expressions():
    kp = KPointInfo(label="Y", coordinates=["a", "1-a", "0"],
                   parameters=[], is_special=True)
    assert kp.parameters == ["a"]
    assert not kp.is_special

def test_kvalue_conversion_rejects_inconsistent_coordinates():
    kp = KPointInfo(label="DT", coordinates=["a", "2a", "0"],
                   parameters=["a"], is_special=False)
    with pytest.raises(ValueError, match="Inconsistent"):
        official_kparams_to_iso(139, "DT", ["1/4"], kp)

def test_coupled_k_parameters_are_solved_as_one_linear_system():
    kp = KPointInfo(label="A", coordinates=["a+b", "a+b", "2a"],
                   parameters=[], is_special=False)
    values = official_kparams_to_iso(139, "A", ["1/7", "1/9"], kp)
    a, b = map(Fraction, values)
    assert a + b == Fraction(1, 7)
    assert 2 * a == Fraction(1, 9)
