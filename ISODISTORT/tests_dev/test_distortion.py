"""
测试：畸变引擎（不依赖二进制）
"""
import numpy as np
from pymatgen.core import Lattice, Structure

from isocore.distortion import DistortionEngine
from isocore.distortion.phase_path import (
    DEFAULT_DISTORTION_TYPES,
    normalize_distortion_types,
)


def test_default_distortion_modes():
    """默认畸变模式应为 displacive + strain 两种（对齐官网 Displacive/Strain）。"""
    assert DEFAULT_DISTORTION_TYPES == ["displacive", "strain"]
    assert normalize_distortion_types(None) == ["displacive", "strain"]
    assert normalize_distortion_types(["strain", "displacement", "strain"]) == [
        "strain",
        "displacive",  # 旧名 displacement 自动映射为 displacive
    ]
    assert normalize_distortion_types("order") == ["occupational"]


def test_distortion_engine_single_mode():
    # 简单立方
    """Relative path: tests_dev/test_distortion.py"""

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A", "B"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    # 人工构造位移：A 原子沿 x 方向位移
    displacements = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])

    engine = DistortionEngine()
    distorted = engine.generate_single_mode(
        struct, displacements, amplitude=0.05
    )

    # 验证位移
    dx = distorted[0].frac_coords[0] - struct[0].frac_coords[0]
    assert abs(dx - 0.05) < 1e-6, f"期望位移 0.05，实际 {dx}"

    # B 原子不动
    dy = distorted[1].frac_coords[1] - struct[1].frac_coords[1]
    assert abs(dy) < 1e-6

    print("✅ 单模式畸变测试通过")


def test_distortion_engine_supercell():
    """Relative path: tests_dev/test_distortion.py"""

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A"], [[0, 0, 0]])

    displacements = np.array([[0.01, 0.0, 0.0]])
    engine = DistortionEngine()
    super_distorted = engine.generate_single_mode(
        struct, displacements, amplitude=1.0, supercell=[2, 2, 2]
    )

    assert len(super_distorted) == 8
    print("✅ 超胞畸变测试通过")


def test_mixed_mode():
    """Relative path: tests_dev/test_distortion.py"""

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A"], [[0, 0, 0]])

    disp_x = np.array([[1.0, 0.0, 0.0]])
    disp_y = np.array([[0.0, 1.0, 0.0]])

    all_disp = {"mode_x": disp_x, "mode_y": disp_y}
    contributions = {"mode_x": 0.1, "mode_y": 0.2}

    engine = DistortionEngine()
    mixed = engine.generate_mixed_mode(struct, contributions, all_disp)

    dx = mixed[0].frac_coords[0]
    dy = mixed[0].frac_coords[1]
    assert abs(dx - 0.1) < 1e-6
    assert abs(dy - 0.2) < 1e-6

    print("✅ 多模式混合畸变测试通过")


if __name__ == "__main__":
    test_distortion_engine_single_mode()
    test_distortion_engine_supercell()
    test_mixed_mode()
