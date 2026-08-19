"""
测试：物理自洽性自检模块（不依赖 WSL 二进制）。
"""
import numpy as np
from pymatgen.core import Lattice, Structure

from isocore.utils.self_check import (
    check_linearity,
    check_mode_orthogonality,
    check_subgroup_rule,
    check_symmetry_conservation,
    check_zero_amplitude,
)


class _StubIso:
    """最小桩：仅提供自检所需的 structure / generate_distortion。"""

    def __init__(self, struct: Structure, disp: np.ndarray):
        self.structure = struct
        self._disp = np.asarray(disp, dtype=float)
        self.mode_displacements = {}

    def generate_distortion(self, irrep_label: str, amplitude: float):
        coords = self.structure.frac_coords + amplitude * self._disp
        return Structure(self.structure.lattice, self.structure.species,
                         coords % 1.0, coords_are_cartesian=False)


def _parent() -> Structure:
    return Structure(Lattice.cubic(4.0), ["Fe", "Fe"],
                     [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])


def test_symmetry_conservation_same_structure():
    p = _parent()
    ok, detail = check_symmetry_conservation(p, p)
    assert ok, detail


def test_symmetry_conservation_extra_ops_breaks():
    """畸变结构出现母相没有的操作（如单斜母相被误建成更高对称）-> 不守恒。"""
    tetragonal = Structure(Lattice.tetragonal(4.0, 8.0), ["Fe", "Fe"],
                           [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    cubic = Structure(Lattice.cubic(4.0), ["Fe", "Fe"],
                      [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    ok, _ = check_symmetry_conservation(tetragonal, cubic)
    assert not ok  # 立方对称操作 ⊄ 四方对称操作


def test_subgroup_rule_same_sg():
    p = _parent()
    ok, detail = check_subgroup_rule(229, 229, p, p)
    assert ok, detail


def test_zero_amplitude_returns_parent():
    stub = _StubIso(_parent(), np.array([[0.2, 0, 0], [0, 0, 0]]))
    ok, detail = check_zero_amplitude(stub, "X")
    assert ok, detail


def test_linearity_doubling():
    stub = _StubIso(_parent(), np.array([[0.2, 0, 0], [0, 0, 0]]))
    ok, detail = check_linearity(stub, "X")
    assert ok, detail


def test_mode_orthogonality():
    ok, _ = check_mode_orthogonality({
        "A": {"displacements": np.array([[1, 0, 0], [0, 0, 0]])},
        "B": {"displacements": np.array([[0, 1, 0], [0, 0, 0]])},
    })
    assert ok
    bad, _ = check_mode_orthogonality({
        "A": {"displacements": np.array([[1, 0, 0], [0, 0, 0]])},
        "B": {"displacements": np.array([[2, 0, 0], [0, 0, 0]])},
    })
    assert not bad
