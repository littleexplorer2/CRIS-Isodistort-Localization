"""
测试：Occupational 占据率畸变模式生成器（不依赖 WSL 二进制）
"""
import numpy as np
from pymatgen.core import Lattice, Structure

from isocore.backend import SubgroupInfo
from isocore.distortion import DistortionEngine, OccupationalModeGenerator


def _parent() -> Structure:
    """简单四方晶胞：Eu(0,0,0) + Al(0.25,0.25,0.25)。"""
    return Structure(Lattice.tetragonal(4.0, 8.0),
                     ["Eu", "Al"], [[0, 0, 0], [0.25, 0.25, 0.25]])


def _wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 1, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0]},
        {"wyckoff_letter": "b", "multiplicity": 1, "species": "Al",
         "representative_index": 1, "equivalent_indices": [1]},
    ]


def _subgroup(basis) -> SubgroupInfo:
    return SubgroupInfo(
        index=0, space_group_number=2, space_group_symbol="P-1",
        basis_vectors=basis, k_point_label="X", irrep_label="X1+",
    )


def test_occupational_mode_doubled_cell():
    """2x1x1 超胞：Al 位点分裂为两类 -> 产生 +1/-1 占据率模式。"""
    parent = _parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _wyckoff(), _subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
                         {"Al"})
    assert modes, "倍胞下 Al 位点应产生占据率模式"
    m = modes[0]
    assert m.species == "Al"
    assert m.label == "occ-Al-b"
    # 模式值只含 +1/-1/0，且 +1 与 -1 数量相等（两类各半）
    assert set(np.unique(m.pattern)) <= {1.0, -1.0, 0.0}
    assert np.count_nonzero(m.pattern == 1) == np.count_nonzero(m.pattern == -1)


def test_occupational_mode_no_split_identity_cell():
    """单位基矢（t 子群，无超胞）：位点不分裂 -> 无占据率模式。"""
    parent = _parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _wyckoff(),
                         _subgroup([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Al"})
    assert modes == []


def test_occupational_scope_filters_species():
    """作用域只含 Eu 时，Al 位点不产生模式。"""
    parent = _parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _wyckoff(),
                         _subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Eu"})
    assert all(m.species == "Eu" for m in modes)


def test_occupancy_exportable_structure():
    """生成的占据率结构可被 pymatgen 构建（部分占据位点）。"""
    parent = _parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _wyckoff(),
                         _subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Al"})
    engine = DistortionEngine()
    m = modes[0]
    sc = engine.generate_modes(
        parent, m.basis_vectors,
        parent_displacements=None,
        occupancy_patterns=[(m.pattern, 0.2)],
    )
    assert len(sc) == 4  # 2 原子 x 2 倍胞
    occs = {s.species_string for s in sc if "Al" in s.species_string}
    assert any(o != "Al" for o in occs)  # 至少一个 Al 位点变为部分占据
