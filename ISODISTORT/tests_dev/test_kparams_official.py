"""官网 k 点参数 ↔ iso KVALUE 换算单元测试。"""
from __future__ import annotations

from fractions import Fraction

from isocore.backend import KPointInfo
from isocore.data.kpoints_official import official_kparams_to_iso


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
