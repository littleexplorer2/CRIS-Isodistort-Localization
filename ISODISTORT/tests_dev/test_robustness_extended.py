"""
异常输入鲁棒性测试（第三层 §2）+ 可重复性测试（第三层 §3）。

原则：非法输入要么正确计算、要么给出清晰错误——绝不静默输出错误结果。
可重复性：同一输入多次运行结果完全一致（逐字节 / 机器精度内）。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort
from isocore.distortion import DistortionEngine
from isocore.io import StructureExporter
from isocore.structure import read_cif

from data_dir import experiment_data_dir

DATA_DIR = experiment_data_dir()
CIFS_DIR = Path(__file__).resolve().parent / "cifs_30"


# ------------------------------------------------------------------
# 第三层 §2：异常输入鲁棒性（不依赖 WSL 的部分）
# ------------------------------------------------------------------

def test_read_cif_missing_file_raises(tmp_path):
    with pytest.raises(OSError):
        read_cif(str(tmp_path / "nope.cif"))


def test_read_cif_empty_file_raises(tmp_path):
    bad = tmp_path / "empty.cif"
    bad.write_text("", encoding="utf-8")
    with pytest.raises((ValueError, KeyError, IndexError)):
        read_cif(str(bad))


def test_read_cif_garbage_raises(tmp_path):
    bad = tmp_path / "garbage.cif"
    bad.write_text("not a cif at all\n" * 10, encoding="utf-8")
    with pytest.raises((ValueError, KeyError, IndexError)):
        read_cif(str(bad))


def test_negative_amplitude_at_api_level(tmp_path):
    """API 层负振幅：数学上合法（位移方向翻转），不崩溃。"""
    cif = CIFS_DIR / "sg001.cif"
    if not cif.exists():
        pytest.skip("测试 CIF 不存在")
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    subs = iso.list_subgroups_at("GM", "GM1")
    found = False
    for sg in subs:
        try:
            iso.search_method_2(subgroup_idx=sg.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
            continue
        if iso.mode_displacements:
            label = next(iter(iso.mode_displacements))
            d = iso.generate_distortion(irrep_label=label, amplitude=-0.1)
            assert len(d) == len(iso.structure)
            found = True
            break
    if not found:
        pytest.skip("该母相无可计算位移模式")


def test_extreme_amplitude_no_crash():
    """极大振幅（1e6）：数值上允许，但必须不崩溃、不产生 NaN。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe", "Fe"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    d = DistortionEngine().generate_single_mode(struct, disp, amplitude=1e6)
    assert np.all(np.isfinite(d.frac_coords))
    assert len(d) == 2


def test_tiny_amplitude_returns_parent_within_tolerance():
    """极小振幅（1e-12）：位移在容差内，结构 ≈ 母相。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe", "Fe"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    d = DistortionEngine().generate_single_mode(struct, disp, amplitude=1e-12)
    diff = d.frac_coords - struct.frac_coords
    diff -= np.round(diff)
    assert np.max(np.abs(diff)) < 1e-6


def test_huge_supercell_atom_count():
    """超大超胞 (10,10,10)：原子数精确 = 1000×。"""
    struct = Structure(Lattice.cubic(3.0), ["Fe"], [[0, 0, 0]])
    d = DistortionEngine().generate_single_mode(
        struct, np.zeros((1, 3)), amplitude=0.0, supercell=[10, 10, 10]
    )
    assert len(d) == 1000


def test_many_atoms_structure_no_crash():
    """原子数极多的结构（200 原子）：引擎不崩溃。"""
    coords = []
    for i in range(200):
        coords.append([(i % 10) / 10.0, ((i // 10) % 10) / 10.0,
                       (i // 100) / 2.0])
    struct = Structure(Lattice.cubic(20.0), ["Fe"] * 200, coords)
    disp = np.tile([1.0, 0.0, 0.0], (200, 1))
    d = DistortionEngine().generate_single_mode(struct, disp, amplitude=0.1)
    assert len(d) == 200
    assert np.all(np.isfinite(d.frac_coords))


def test_nan_displacement_rejected():
    """含 NaN 的位移向量必须报错（不能静默产生非法结构）。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe"], [[0, 0, 0]])
    disp = np.array([[np.nan, 0.0, 0.0]])
    with pytest.raises((ValueError, FloatingPointError)):
        DistortionEngine().generate_single_mode(struct, disp, amplitude=0.1)


def test_zero_length_displacement_ok():
    """全零位移 + 非零振幅 = 母相（合法，不报错）。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe"], [[0, 0, 0]])
    d = DistortionEngine().generate_single_mode(
        struct, np.zeros((1, 3)), amplitude=0.5
    )
    diff = d.frac_coords - struct.frac_coords
    diff -= np.round(diff)
    assert np.max(np.abs(diff)) < 1e-9


# ------------------------------------------------------------------
# 第三层 §3：可重复性（确定性群论计算，结果必须逐次一致）
# ------------------------------------------------------------------

def test_engine_deterministic_repeat():
    """同一输入两次生成：逐字节一致（确定性）。"""
    struct = Structure(Lattice.cubic(4.0), ["Fe", "O"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    engine = DistortionEngine()
    d1 = engine.generate_single_mode(struct, disp, amplitude=0.1)
    d2 = engine.generate_single_mode(struct, disp, amplitude=0.1)
    assert d1.frac_coords.shape == d2.frac_coords.shape
    assert np.allclose(d1.frac_coords, d2.frac_coords, atol=1e-12)
    # 导出 CIF 逐字节一致
    out1 = StructureExporter().to_cif(d1, "rep_test_1")
    out2 = StructureExporter().to_cif(d2, "rep_test_2")
    assert out1.read_bytes() == out2.read_bytes(), \
        "同一输入两次导出的 CIF 必须逐字节一致"


def test_read_cif_roundtrip_stable(tmp_path):
    """CIF 读回两次：识别空间群一致（可重复性）。"""
    cif = CIFS_DIR / "sg123.cif"
    if not cif.exists():
        pytest.skip("测试 CIF 不存在")
    s1 = read_cif(str(cif))
    s2 = read_cif(str(cif))
    g1 = SpacegroupAnalyzer(s1, symprec=1e-3).get_space_group_number()
    g2 = SpacegroupAnalyzer(s2, symprec=1e-3).get_space_group_number()
    assert g1 == g2 == 123
    assert np.allclose(s1.frac_coords, s2.frac_coords, atol=1e-12)


def test_method1_options_repeatable(tmp_path):
    """API 层 method1_options 两次调用：选项列表一致（可重复性）。"""
    cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(cif)
    o1 = iso.method1_options()
    o2 = iso.method1_options()
    assert [g["number"] for g in o1["space_groups"]] == \
        [g["number"] for g in o2["space_groups"]]
    assert [g["label"] for g in o1["conventional_lattices"]] == \
        [g["label"] for g in o2["conventional_lattices"]]
