"""
测试：结构层工具（不依赖二进制）
"""
import numpy as np
from isocore.structure import (
    SymmetryValidator,
    build_supercell,
    coordinates_are_equal,
    wrap_to_unit_cell,
)
from pymatgen.core import Lattice, Structure


def test_coordinate_utils():
    # wrap_to_unit_cell
    """Relative path: tests/test_structure.py"""

    coords = np.array([[1.2, -0.3, 0.5]])
    wrapped = wrap_to_unit_cell(coords)
    assert abs(wrapped[0, 0] - 0.2) < 1e-6
    assert abs(wrapped[0, 1] - 0.7) < 1e-6

    # coordinates_are_equal
    assert coordinates_are_equal(np.array([0.1, 0.2, 0.3]),
                                np.array([0.1, 0.2, 0.3]))
    assert coordinates_are_equal(np.array([0.0, 0.0, 0.0]),
                                np.array([1.0, 1.0, 1.0]))
    assert not coordinates_are_equal(np.array([0.0, 0.0, 0.0]),
                                    np.array([0.5, 0.0, 0.0]))

    print("✅ 坐标工具测试通过")


def test_supercell():
    # 简单立方原胞
    """Relative path: tests/test_structure.py"""

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["Na"], [[0, 0, 0]])
    supercell = build_supercell(struct, [2, 2, 2])
    assert len(supercell) == 8
    print("✅ 超胞构建测试通过")


def test_symmetry_validator():
    # NaCl 结构（P 原胞：Na 与 Cl 各一个，晶系为 Pm-3m #221）
    """Relative path: tests/test_structure.py"""

    lattice = Lattice.cubic(5.63)
    struct = Structure(
        lattice,
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    validator = SymmetryValidator()
    result = validator.validate(struct)

    assert result["space_group_number"] == 221
    assert len(result["wyckoff_sites"]) == 2
    assert not result["has_disorder"]

    print("✅ 对称性校验测试通过")
    print(f"   空间群: #{result['space_group_number']} ({result['space_group_symbol']})")
    print(f"   Wyckoff 位点: {[s['wyckoff_letter'] for s in result['wyckoff_sites']]}")


if __name__ == "__main__":
    test_coordinate_utils()
    test_supercell()
    test_symmetry_validator()
