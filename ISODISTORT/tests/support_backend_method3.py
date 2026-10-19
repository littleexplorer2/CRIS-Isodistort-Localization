"""Tests for support_backend_method3; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.wrappers import BushMode as _search_methods_BushMode
from backend.wrappers import DistortionMode as _search_methods_DistortionMode
from backend.wrappers import SubgroupInfo as _search_methods_SubgroupInfo
from dataclasses import dataclass
from pymatgen.core import Lattice
from pymatgen.core import Structure


# Original coverage: test_search_methods.py.

def _i4mmm_parent() -> Structure:
    return Structure.from_spacegroup(
        "I4/mmm",
        Lattice.tetragonal(4.0, 10.0),
        ["Eu"],
        [[0, 0, 0]],
    )

@dataclass
class _DummyIsoWrapper:
    """最小桩对象：提供 search engine 所需接口。"""

    def list_k_points(self, parent_sg):
        return []

    def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
        return [
            _search_methods_SubgroupInfo(index=0, space_group_number=225, space_group_symbol="Fm-3m",
                         subgroup_index=1, size=1, is_maximal=True,
                         opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
                         basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         k_point_label="GM", irrep_label="GM1+"),
            _search_methods_SubgroupInfo(index=1, space_group_number=139, space_group_symbol="I4/mmm",
                         subgroup_index=4, size=1, is_maximal=True,
                         opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
                         basis_vectors=[[0, 1, 0], [0, 0, 1], [1, 0, 0]],
                         k_point_label="X", irrep_label="X1"),
            _search_methods_SubgroupInfo(index=2, space_group_number=62, space_group_symbol="Pnma",
                         subgroup_index=8, size=1, is_maximal=False,
                         opd_symbol="C1", opd_vector=[1.0, 1.0, 0.0],
                         basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         k_point_label="L", irrep_label="L2"),
        ]

    def calc_distortion_modes(self, parent_sg, subgroup, wyckoff_letters):
        return [
            _search_methods_DistortionMode(
                irrep_label=subgroup.irrep_label,
                dimension=1,
                opd_symbol=subgroup.opd_symbol,
                bush_modes=[
                    _search_methods_BushMode(irrep_label=subgroup.irrep_label,
                             opd_symbol=subgroup.opd_symbol,
                             wyckoff_letter="a",
                             point=[0.0, 0.0, 0.0],
                             displacements=[[1.0, 0.0, 0.0]]),
                ],
            ),
        ]

def _method3_route(
    irrep: str,
    *,
    index: int = 0,
    origin: list[float] | None = None,
    k_parameters: list[str] | None = None,
) -> _search_methods_SubgroupInfo:
    return _search_methods_SubgroupInfo(
        index=index,
        space_group_number=99,
        space_group_symbol="P4mm",
        subgroup_index=24,
        size=12,
        opd_symbol="C1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
        origin=origin or [0, 0, 0],
        k_point_label="LD",
        irrep_label=irrep,
        k_parameters=k_parameters or ["1/6"],
        parent_sg=139,
    )
