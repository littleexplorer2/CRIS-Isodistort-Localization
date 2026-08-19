"""
测试：畸变模式映射器（DistortionMapper）—— 带自由参数代表点的正确分配。

回归保护：BUSH 输出中带参数的代表点（如 ``(0,0,z)`` 与 ``(0,0,-z)``）必须
按结构真实坐标解析并按“模母相格点等价”分配，否则相反位移会落到错误的原子上
（实测曾导致全对称 GM1+ 模式把 I4/mmm 畸变成 P4mm 等错误对称性）。
"""
import numpy as np
from pymatgen.core import Lattice, Structure

from isocore.backend import BushMode, DistortionMode
from isocore.distortion import DistortionMapper


def _parent() -> Structure:
    """四方母相：Eu 2a + Al 2e（(0,0,0.38)/(0,0,0.62)）。"""
    return Structure(
        Lattice.tetragonal(4.4, 11.2),
        ["Eu", "Eu", "Al", "Al"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0, 0, 0.38], [0, 0, 0.62]],
    )


def _wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 2, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0, 1]},
        {"wyckoff_letter": "e", "multiplicity": 2, "species": "Al",
         "representative_index": 2, "equivalent_indices": [2, 3]},
    ]


def _mode() -> DistortionMode:
    mode = DistortionMode(irrep_label="GM1+", opd_symbol="P1",
                          k_point_label="GM")
    mode.bush_modes = [
        BushMode("GM1+", "P1", "a", [0, 0, 0], ["0", "0", "0"], [[0, 0, 0]]),
        # e 位点：两个相反位移的代表（(0,0,z) 与 (0,0,-z)）
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "-z"], [[0, 0, -1]]),
    ]
    return mode


def test_mapper_parameterized_z_points():
    """(0,0,z) 与 (0,0,-z) 必须解析为不同位置：z=0.38 -> +z，z=0.62 -> -z。"""
    mapper = DistortionMapper()
    result = mapper.map_modes_to_atoms(_parent(), _wyckoff(), [_mode()])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[0], 0) and np.allclose(disp[1], 0)  # Eu 不动
    assert np.allclose(disp[2], [0, 0, 1])   # (0,0,0.38) -> +z
    assert np.allclose(disp[3], [0, 0, -1])  # (0,0,0.62) -> -z


def test_mapper_uniform_single_rep():
    """单代表行（均匀模式）：同一位移作用于该位点全部原子。"""
    mode = _mode()
    mode.bush_modes = [
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
    ]
    mapper = DistortionMapper()
    result = mapper.map_modes_to_atoms(_parent(), _wyckoff(), [mode])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[2], [0, 0, 1])
    assert np.allclose(disp[3], [0, 0, 1])
