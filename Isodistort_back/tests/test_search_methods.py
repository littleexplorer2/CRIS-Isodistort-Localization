"""
测试：ISODISTORT Method 1-4 搜索流程（不依赖真实二进制）。
"""

from dataclasses import dataclass

import numpy as np
from pymatgen.core import Lattice, Structure

from isocore.backend import DistortionMode, SubgroupInfo
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

    def list_subgroups(self, parent_sg, distortion_types=None):
        return [
            SubgroupInfo(index=1, space_group_number=225, space_group_symbol="Fm-3m"),
            SubgroupInfo(index=2, space_group_number=139, space_group_symbol="I4/mmm"),
            SubgroupInfo(index=3, space_group_number=62, space_group_symbol="Pnma"),
        ]

    def calc_distortion_modes(self, parent_sg, subgroup_idx, distortion_type="displacement"):
        return [
            DistortionMode(irrep_label="GM1+", dimension=1, basis_vectors=[[1.0, 0.0, 0.0]]),
            DistortionMode(irrep_label="R5-", dimension=3, basis_vectors=[[0.0, 1.0, 0.0]]),
        ]


def test_method_1_filters():
    engine = IsoSearchEngine(_DummyIsoWrapper())

    all_items = engine.method_1_search(225, Method1Query())
    assert len(all_items) == 3

    cubic_only = engine.method_1_search(225, Method1Query(crystal_system="cubic"))
    assert len(cubic_only) == 1
    assert cubic_only[0].subgroup.space_group_number == 225



def test_method_2_metadata_and_modes():
    engine = IsoSearchEngine(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().list_subgroups(225)

    result = engine.method_2_search(
        225,
        subs,
        Method2Query(
            subgroup_idx=2,
            distortion_type="displacement",
            k_point_label="X",
            k_point_coordinates=["1/2", "0", "0"],
            number_of_superposed_irs=2,
        ),
    )

    assert result.subgroup.index == 2
    assert len(result.modes) == 2
    assert result.metadata["k_point_coordinates"] == [0.5, 0.0, 0.0]



def test_method_3_space_group_precedence():
    engine = IsoSearchEngine(_DummyIsoWrapper())

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
    engine = IsoSearchEngine(_DummyIsoWrapper())

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
