"""
测试：结构层工具（不依赖二进制）
"""
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.structure import (
    SymmetryValidator,
    build_supercell,
    coordinates_are_equal,
    wrap_to_unit_cell,
)


def test_coordinate_utils():
    # wrap_to_unit_cell
    """Relative path: tests_dev/test_structure.py"""

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
    """Relative path: tests_dev/test_structure.py"""

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["Na"], [[0, 0, 0]])
    supercell = build_supercell(struct, [2, 2, 2])
    assert len(supercell) == 8
    print("✅ 超胞构建测试通过")


def test_fractional_supercell_subcell():
    """回归：子群基矢行列式 < 1（带心母相的亚胞）不能走 pymatgen 整数路径。

    对应真实 bug：Fm-3m（4 原子惯用胞）的 I4/mmm 子群基矢
    det=0.5，pymatgen ``Structure * matrix`` 只接受整数矩阵，
    曾抛 LinAlgError: Singular matrix。
    """
    lattice = Lattice.cubic(3.6)
    struct = Structure(
        lattice,
        ["Fe", "Fe", "Fe", "Fe"],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )
    basis = [[-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.0, 0.0, 1.0]]
    sub = build_supercell(struct, basis)
    # 亚胞：4 原子 -> 2 原子（重合位点按周期合并）
    assert len(sub) == 2
    sg = SpacegroupAnalyzer(sub, symprec=1e-3).get_space_group_number()
    # 零畸变的 I4/mmm 亚胞仍为立方 Fm-3m（c/a=√2 的 bct 就是 fcc 晶格）
    assert sg == 225

    # 整数矩阵仍走原路径（回归保护）
    assert len(build_supercell(struct, np.eye(3, dtype=int))) == 4
    print("✅ 分数基矢亚胞构建测试通过")


def test_fractional_supercell_subcell_multi_species_merge():
    """回归：多物种亚胞的周期合并不得因索引错位漏并（原子重叠）。

    对应真实 bug：NaCl（Fm-3m，4 Na + 4 Cl）的 Imm2 子群基矢 det=0.5，
    4 个 Cl 映射到亚胞后两两重合。旧合并逻辑用 ``structure[j]`` 取
    「已收集坐标」对应的物种——当 Na 先被合并跳过（coords 索引与
    structure 索引错位）后，Cl 的物种检查会拿 Cl 与 Na 比较，重合的
    Cl 无法合并，畸变结构出现原子重叠，spglib 抛
    SymmetryUndeterminedError（NaCl/MgO 外部真实结构复现）。
    """
    lattice = Lattice.cubic(5.64)
    struct = Structure(
        lattice,
        ["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
         [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0]],
    )
    # Imm2 子群基矢（det=0.5）：8 原子 -> 4 原子（2 Na + 2 Cl）
    basis = [[0.0, 0.0, -1.0], [-0.5, 0.5, 0.0], [0.5, 0.5, 0.0]]
    sub = build_supercell(struct, basis)
    assert len(sub) == 4, f"亚胞应合并为 4 原子，实际 {len(sub)}"
    # 无原子重叠：合并后恰为 2 Na + 2 Cl（str 为 "Na1"/"Cl1" 形式）
    species_list = [str(s.species) for s in sub]
    assert species_list.count("Cl1") == 2 and species_list.count("Na1") == 2
    print("✅ 多物种亚胞合并回归测试通过")


def test_symmetry_validator():
    # NaCl 结构（P 原胞：Na 与 Cl 各一个，晶系为 Pm-3m #221）
    """Relative path: tests_dev/test_structure.py"""

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
    print(f"   Wyckoff 位置: {[s['wyckoff_letter'] for s in result['wyckoff_sites']]}")


if __name__ == "__main__":
    test_coordinate_utils()
    test_supercell()
    test_symmetry_validator()
