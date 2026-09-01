"""Basics: config, i18n, parent header, official k-params."""
from __future__ import annotations

import os
from fractions import Fraction

from pymatgen.core import Lattice, Structure

from isocore.backend import KPointInfo
from isocore.data.kpoints_official import official_kparams_to_iso
from isocore.i18n import MESSAGES, t
from isocore.utils import get_config
from isocore.utils.parent_header import (
    format_fixed_coord,
    format_wyckoff_site,
    format_wyckoff_sites,
)


# --- from test_config.py ---

def test_config_load():
    cfg = get_config()
    # isobyu 二进制路径可解析（存在性由部署保证；WSL 下为 Windows 侧路径）
    assert cfg.iso_bin.name == "iso"
    assert cfg.findsym_bin.name == "findsym"
    assert cfg.temp_dir.exists(), "临时目录已创建"
    assert cfg.output_dir.exists(), "输出目录已创建"
    assert "ISODATA" in os.environ, "ISODATA 环境变量已设置"
    print("配置加载测试通过")
    print(f"   ISODATA = {os.environ['ISODATA']}")




# --- from test_i18n.py ---

def test_messages_are_flat_english():
    assert isinstance(MESSAGES, dict)
    assert "zh" not in MESSAGES
    assert "en" not in MESSAGES
    assert t("load.done", sg=139, sym="I4/mmm", n=10).startswith("[Loaded]")
    assert "Space group" in t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "139" in t("load.done", sg=139, sym="I4/mmm", n=10)


def test_unknown_key_returns_key():
    assert t("no.such.key") == "no.such.key"


def test_no_language_switch_keys():
    assert "ui.menu.language" not in MESSAGES
    assert "ui.lang.current" not in MESSAGES
    assert "dist.gen" not in MESSAGES
    assert "dist.domainsBtn" not in MESSAGES


def test_distortion_and_method2_help_keys():
    assert "m2.genDbHelp" in MESSAGES
    assert "isotropy-subgroup" in MESSAGES["m2.genDbHelp"].lower() or (
        "isotropy" in MESSAGES["m2.genDbHelp"].lower()
        and "subgroup" in MESSAGES["m2.genDbHelp"].lower()
    )
    assert "cache" in MESSAGES["m2.genDbHelp"].lower()
    assert "\n" not in MESSAGES["m2.genDbHelp"]
    assert "m2.genDbManage" in MESSAGES
    assert "m2.genDbDelete" in MESSAGES
    assert "dist.method4" in MESSAGES
    assert "dist.tableLabel" in MESSAGES
    assert "not limited by the filter" not in MESSAGES["dist.tableNote"]
    assert "dist.zipNotM4" in MESSAGES
    assert "dist.zipWait" in MESSAGES
    assert "dist.zipDone" in MESSAGES
    assert "st.elapsed" in MESSAGES
    assert "st.busyHint" in MESSAGES
    assert "m2.nmodRemoved" in MESSAGES
    assert "not available locally" in MESSAGES["m2.nmodRemoved"].lower()
    assert "independent incommensurate" in MESSAGES["l.nmod"].lower()
    assert "ss.title" not in MESSAGES
    assert "ui.menu.superspace" not in MESSAGES
    assert "ok.nmod" not in MESSAGES


def test_numbered_web_placeholders():
    text = t("ok.nsup", 3)
    assert "3" in text


# --- from test_parent_header.py ---

def test_format_fixed_coord_fractions():
    assert format_fixed_coord(0.0) == "0"
    assert format_fixed_coord(0.5) == "1/2"
    assert format_fixed_coord(0.25) == "1/4"


def test_format_wyckoff_site_free_z():
    line = format_wyckoff_site("Al", 2, 4, "e", [0.0, 0.0, 0.38])
    assert line == "Al2 4e (0,0,z), z= 0.38000"


def test_format_wyckoff_sites_eual4_like():
    lattice = Lattice.tetragonal(4.402, 11.163)
    structure = Structure(
        lattice,
        ["Eu", "Eu", "Al", "Al", "Al", "Al", "Al", "Al", "Al", "Al"],
        [
            [0, 0, 0], [0.5, 0.5, 0.5],
            [0, 0.5, 0.25], [0.5, 0, 0.25], [0.5, 0, 0.75], [0, 0.5, 0.75],
            [0, 0, 0.38], [0, 0, 0.62], [0.5, 0.5, 0.12], [0.5, 0.5, 0.88],
        ],
    )
    sites = [
        {"species": "Eu", "multiplicity": 2, "wyckoff_letter": "a",
         "representative_index": 0},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "d",
         "representative_index": 2},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "e",
         "representative_index": 6},
    ]
    lines = format_wyckoff_sites(structure, sites)
    assert lines == [
        "Eu1 2a (0,0,0)",
        "Al1 4d (0,1/2,1/4)",
        "Al2 4e (0,0,z), z= 0.38000",
    ]


def test_parent_wyckoff_from_cif_ndnio2():
    """CIF atom_site labels/order drive the web header (not memorized EuAl4)."""
    from pathlib import Path

    from isocore.api import IsoDistort

    cris = Path(__file__).resolve().parents[3]
    cif = cris / "experiment_data" / "NdNiO2 own.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("NdNiO2 own.cif not in experiment_data")
    iso = IsoDistort()
    iso.load_structure(cif)
    lines = iso.parent_wyckoff_display()
    assert lines[0].startswith("O 2f")
    assert any(ln.startswith("ND 1d") for ln in lines)
    assert any(ln.startswith("NI 1a") for ln in lines)


def test_export_cif_uses_cif_wyckoff_lines_ndnio2():
    """Distortion CIF comments must reuse the same CIF-ordered parent header."""
    from pathlib import Path

    from isocore.api import IsoDistort
    from isocore.backend import SubgroupInfo
    from isocore.io.distortion_formats import render_cif

    cris = Path(__file__).resolve().parents[3]
    cif = cris / "experiment_data" / "NdNiO2 own.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("NdNiO2 own.cif not in experiment_data")
    iso = IsoDistort()
    iso.load_structure(cif)
    sg = SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        size=2,
        subgroup_index=2,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label="A1+",
        k_point_label="A",
        parent_sg=123,
        k_coordinates=["1/2", "1/2", "1/2"],
        basis_vectors=[[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
        origin=[0.0, 0.0, 0.0],
    )
    spec = iso._spec_for_subgroup(
        sg, use_current_modes=False, use_generated_structure=False
    )
    assert spec.parent_wyckoff_lines
    assert spec.parent_wyckoff_lines[0].startswith("O 2f")
    text = render_cif(spec.structure, spec)
    assert "# O 2f (0,1/2,0)" in text
    assert "# ND 1d (1/2,1/2,1/2)" in text
    assert "# NI 1a (0,0,0)" in text
    assert "# Nd1" not in text



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

