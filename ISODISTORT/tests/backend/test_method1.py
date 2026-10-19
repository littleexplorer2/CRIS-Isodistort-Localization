"""Tests for backend.method1; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.method1 import Method1Query
from backend.method1 import Method1Solver
from backend.method1 import OccupationalModeGenerator
from backend.method1.search import _basis_is_same_lattice
from backend.utils.lattice import is_sublattice
from backend.utils.lattice import rational_matrix
from backend.wrappers import SubgroupInfo
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.support_backend_method3 import _DummyIsoWrapper
import numpy as np
import pytest


# Original coverage: test_distortion.py.

# --- from test_occupational_modes.py ---

def _occ_parent() -> Structure:
    """简单四方晶胞：Eu(0,0,0) + Al(0.25,0.25,0.25)。"""
    return Structure(Lattice.tetragonal(4.0, 8.0),
                     ["Eu", "Al"], [[0, 0, 0], [0.25, 0.25, 0.25]])

def _occ_wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 1, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0]},
        {"wyckoff_letter": "b", "multiplicity": 1, "species": "Al",
         "representative_index": 1, "equivalent_indices": [1]},
    ]

def _occ_subgroup(basis) -> SubgroupInfo:
    return SubgroupInfo(
        index=0, space_group_number=2, space_group_symbol="P-1",
        basis_vectors=basis, k_point_label="X", irrep_label="X1+",
    )

@pytest.mark.parametrize("basis,scope", [
    ([[2, 0, 0], [0, 1, 0], [0, 0, 1]], {"Al"}),
    ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], {"Al"}),
    ([[2, 0, 0], [0, 1, 0], [0, 0, 1]], {"Eu"}),
])
def test_occupational_generation_requires_scalar_representation(basis, scope):
    """Neither supercell size nor species selection proves a scalar mode."""
    parent = _occ_parent()
    with pytest.raises(NotImplementedError, match="occupational"):
        OccupationalModeGenerator(tolerance=1e-4).generate(
            parent, _occ_wyckoff(), _occ_subgroup(basis), scope,
        )

# Original coverage: test_search_methods.py.

def test_method1_selector_matches_one_lattice_class():
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    doubled = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    assert _basis_is_same_lattice(identity, identity)
    assert not _basis_is_same_lattice(doubled, identity)

def test_sublattice_check_parent_point_group_orbit():
    """官网 lattice 选项含母相点群旋转：B 相对 S@R 亦算命中。"""
    s = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    # 右乘点群旋转：把 (2,1,1) 型子格转到 (1,2,1) 型（一般非 GL 等价）
    r = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    b = (np.asarray(s, dtype=float) @ r.T).tolist()
    assert is_sublattice(rational_matrix(b), rational_matrix(s)) is False
    assert is_sublattice(rational_matrix(b), rational_matrix(np.asarray(s) @ r.T)) is True

def test_method_1_filters():
    engine = Method1Solver(_DummyIsoWrapper())

    all_items = engine.search(225, Method1Query())
    assert len(all_items) == 3

    cubic_only = engine.search(225, Method1Query(crystal_system="cubic"))
    assert len(cubic_only) == 1
    assert cubic_only[0].subgroup.space_group_number == 225

    sg_filtered = engine.search(
        225, Method1Query(subgroup_space_group=62)
    )
    assert len(sg_filtered) == 1
    assert sg_filtered[0].subgroup.space_group_number == 62

    maximal_only = engine.search(
        225, Method1Query(maximal_subgroup_only=True)
    )
    assert all(item.is_maximal for item in maximal_only)
    assert len(maximal_only) == 2

def test_method_1_lattice_filter():
    engine = Method1Solver(_DummyIsoWrapper())
    # 请求 (2,2,2) 子格（对角阵）：仅基矢可被 2 整除的子群保留
    items = engine.search(
        225, Method1Query(lattice=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    )
    # 桩数据中 index=2 (Pnma, 单位基矢) 被过滤
    assert len(items) == 0 or all(
        item.subgroup.index != 2 for item in items
    )

# Original coverage: test_theory_contracts.py.

def test_unsupported_occupational_modes_do_not_partition_coordinates():
    parent = Structure(Lattice.cubic(4), ["Si"], [[0, 0, 0]])
    subgroup = SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1",
                            basis_vectors=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                            k_point_label="X", irrep_label="X1")
    sites = [{"wyckoff_letter": "a", "species": "Si", "equivalent_indices": [0]}]
    with pytest.raises(NotImplementedError, match="occupational"):
        OccupationalModeGenerator().generate(parent, sites, subgroup, {"Si"})
