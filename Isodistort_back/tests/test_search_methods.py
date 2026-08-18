"""
测试：ISODISTORT Method 1-4 搜索流程（使用桩对象，不依赖真实二进制）。
"""

from dataclasses import dataclass

import numpy as np
from pymatgen.core import Lattice, Structure

from isocore.backend import BushMode, DistortionMode, SubgroupInfo
from isocore.distortion import (
    IsoSearchEngine,
    Method1Query,
    Method2Query,
    Method3Query,
    Method4Query,
)


@dataclass
class _DummyIsoWrapper:
    """最小桩对象：提供 search engine 所需接口。"""

    def enumerate_all_special_subgroups(self, parent_sg, distortion_types=None):
        return [
            SubgroupInfo(index=0, space_group_number=225, space_group_symbol="Fm-3m",
                         subgroup_index=1, size=1, is_maximal=True,
                         opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
                         basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         k_point_label="GM", irrep_label="GM1+"),
            SubgroupInfo(index=1, space_group_number=139, space_group_symbol="I4/mmm",
                         subgroup_index=4, size=1, is_maximal=True,
                         opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
                         basis_vectors=[[0, 1, 0], [0, 0, 1], [1, 0, 0]],
                         k_point_label="X", irrep_label="X1"),
            SubgroupInfo(index=2, space_group_number=62, space_group_symbol="Pnma",
                         subgroup_index=8, size=1, is_maximal=False,
                         opd_symbol="C1", opd_vector=[1.0, 1.0, 0.0],
                         basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         k_point_label="L", irrep_label="L2"),
        ]

    def calc_distortion_modes(self, parent_sg, subgroup, wyckoff_letters):
        return [
            DistortionMode(
                irrep_label=subgroup.irrep_label,
                dimension=1,
                opd_symbol=subgroup.opd_symbol,
                bush_modes=[
                    BushMode(irrep_label=subgroup.irrep_label,
                             opd_symbol=subgroup.opd_symbol,
                             wyckoff_letter="a",
                             point=[0.0, 0.0, 0.0],
                             displacements=[[1.0, 0.0, 0.0]]),
                ],
            ),
        ]


def _make_engine() -> IsoSearchEngine:
    return IsoSearchEngine(_DummyIsoWrapper())


def test_method_1_filters():
    engine = _make_engine()

    all_items = engine.method_1_search(225, Method1Query())
    assert len(all_items) == 3

    cubic_only = engine.method_1_search(225, Method1Query(crystal_system="cubic"))
    assert len(cubic_only) == 1
    assert cubic_only[0].subgroup.space_group_number == 225

    sg_filtered = engine.method_1_search(
        225, Method1Query(subgroup_space_group=62)
    )
    assert len(sg_filtered) == 1
    assert sg_filtered[0].subgroup.space_group_number == 62

    maximal_only = engine.method_1_search(
        225, Method1Query(maximal_subgroup_only=True)
    )
    assert all(item.is_maximal for item in maximal_only)
    assert len(maximal_only) == 2


def test_method_1_direct_sublattice_filter():
    engine = _make_engine()
    # 请求 (2,2,2) 直接子格：仅基矢可被 2 整除的子群保留
    items = engine.method_1_search(
        225, Method1Query(direct_sublattice=[2, 2, 2])
    )
    # 桩数据中 index=2 (Pnma, 单位基矢) 被过滤
    assert len(items) == 0 or all(
        item.subgroup.index != 2 for item in items
    )


def test_method_2_metadata_and_modes():
    engine = _make_engine()
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)

    result = engine.method_2_search(
        225,
        subs,
        Method2Query(
            subgroup_idx=1,
            distortion_type="displacement",
            k_point_label="X",
            k_point_coordinates=["1/2", "0", "0"],
            number_of_superposed_irs=2,
        ),
        wyckoff_letters=["a", "b"],
    )

    assert result.subgroup.index == 1
    assert result.subgroup.space_group_number == 139
    assert len(result.modes) == 1
    assert result.modes[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]
    assert result.metadata["k_point_coordinates"] == [0.5, 0.0, 0.0]


def test_method_2_unknown_index_raises():
    engine = _make_engine()
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    try:
        engine.method_2_search(
            225, subs, Method2Query(subgroup_idx=99),
            wyckoff_letters=["a"],
        )
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_method_3_space_group_precedence():
    engine = _make_engine()

    items = engine.method_3_search(
        225,
        Method3Query(
            point_group="m-3m",
            space_group_type=62,
            supercell_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        ),
    )

    assert len(items) == 1
    assert items[0].subgroup.space_group_number == 62


def test_method_4_decomposition_recovery():
    engine = _make_engine()

    parent = Structure(Lattice.cubic(5.0), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    distorted = parent.copy()
    distorted.translate_sites([0], [0.03, 0.0, 0.0], frac_coords=True)

    mode_displacements = {
        "GM1+": np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "R5-": np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
    }

    result = engine.method_4_decompose(
        parent,
        distorted,
        mode_displacements,
        Method4Query(atom_matching_method="nearest-site"),
    )

    assert abs(result.amplitudes["GM1+"] - 0.03) < 1e-8
    assert abs(result.amplitudes["R5-"]) < 1e-8
    assert result.rms_residual < 1e-10
