"""
测试：鲁棒性（异常输入）与可重复性。

异常输入原则：要么正确计算，要么给出清晰错误提示——绝不静默输出错误结果。
"""
import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from isocore.backend import SubgroupInfo
from isocore.distortion import DistortionEngine, OccupationalModeGenerator
from isocore.structure import read_cif


def test_negative_amplitude_is_allowed():
    """负振幅：位移方向翻转，不崩溃（合法数学输入）。"""
    struct = Structure(Lattice.cubic(5.0), ["A", "B"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    engine = DistortionEngine()
    d = engine.generate_single_mode(struct, disp, amplitude=-0.1)
    assert len(d) == 2
    # 原子 0 沿 -x 位移（≠ 原坐标）
    assert abs(d[0].frac_coords[0]) > 1e-6


def test_wrong_sized_displacements_raise():
    """位移向量与原子数不一致必须报错，不能静默。"""
    struct = Structure(Lattice.cubic(5.0), ["A", "B"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    engine = DistortionEngine()
    with pytest.raises(ValueError):
        engine.generate_single_mode(struct, np.zeros((3, 3)), amplitude=0.1)


def test_invalid_basis_matrix_raises():
    """非法基矢（非 3x3）必须报错。"""
    struct = Structure(Lattice.cubic(5.0), ["A"], [[0, 0, 0]])
    engine = DistortionEngine()
    with pytest.raises(ValueError):
        engine.generate_single_mode(struct, np.zeros((1, 3)),
                                    amplitude=0.1, supercell=[[1, 0], [0, 1]])


def test_corrupt_cif_raises_clear_error(tmp_path):
    """损坏的 CIF 必须给出清晰错误而非崩溃。"""
    bad = tmp_path / "bad.cif"
    bad.write_text("this is not a cif file at all\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_cif(str(bad))


def test_occupancy_amplitude_clamped():
    """occupational 幅度 >1 会被夹取，不产生非法占据率。"""
    parent = Structure(Lattice.tetragonal(4.0, 8.0),
                       ["Al", "Al"], [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]])
    wyckoff = [
        {"wyckoff_letter": "b", "multiplicity": 1, "species": "Al",
         "representative_index": 0, "equivalent_indices": [0, 1]},
    ]
    subgroup = SubgroupInfo(index=0, space_group_number=2,
                            space_group_symbol="P-1",
                            basis_vectors=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                            k_point_label="X", irrep_label="X1+")
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, wyckoff, subgroup, {"Al"})
    if modes:
        sc = DistortionEngine().generate_modes(
            parent, modes[0].basis_vectors,
            occupancy_patterns=[(modes[0].pattern, 5.0)],
        )
        for site in sc:
            for amt in site.species.values():
                assert 0.0 < amt <= 1.0 + 1e-6
