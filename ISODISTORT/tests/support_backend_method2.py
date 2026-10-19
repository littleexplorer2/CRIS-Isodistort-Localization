"""Tests for support_backend_method2; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.method2.superspace import _select_verified_microscopic_modes
from backend.models.iso_mode_models import InvariantDirection
from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import MicroscopicVectorBlock
from backend.models.iso_mode_models import MicroscopicVectorRow
from backend.models.iso_mode_models import ModeIdentity
from backend.wrappers.iso_wrapper import BushMode
from backend.wrappers.iso_wrapper import BushMode as _partial_microscopic_mapping_BushMode
from backend.wrappers.iso_wrapper import DistortionMode
from backend.wrappers.iso_wrapper import DistortionMode as _partial_microscopic_mapping_DistortionMode
from backend.wrappers.iso_wrapper import IsoWrapper as _microscopic_column_provenance_IsoWrapper
from backend.wrappers.iso_wrapper import SubgroupInfo
from backend.wrappers.iso_wrapper import SubgroupInfo as _microscopic_column_provenance_SubgroupInfo
from fractions import Fraction
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.support_paths import TESTS_ROOT
from types import SimpleNamespace
import hashlib
import numpy as np


# Original coverage: test_microscopic_column_provenance.py.

_microscopic_column_provenance_QUERY_DIGEST = hashlib.sha256(b"synthetic microscopic query").hexdigest()

def _block(
    *,
    irrep: str = "X1+",
    source_order: int = 0,
) -> MicroscopicVectorBlock:
    return MicroscopicVectorBlock(
        global_irrep=irrep,
        wyckoff_letter="a",
        site_irrep="A",
        rows=(
            MicroscopicVectorRow(
                point_raw=("0", "0", "0"),
                displacements=(
                    (Fraction(1), Fraction(0), Fraction(0)),
                    (Fraction(0), Fraction(1), Fraction(0)),
                ),
            ),
            MicroscopicVectorRow(
                point_raw=("1", "0", "0"),
                displacements=(
                    (Fraction(-1), Fraction(0), Fraction(0)),
                    (Fraction(0), Fraction(-1), Fraction(0)),
                ),
            ),
        ),
        source_order=source_order,
        direction_symbol="VECTOR,A,B",
    )

def _provenance(column: int, *, block: MicroscopicVectorBlock | None = None):
    return MicroscopicColumnProvenance.from_exact_column(
        query_digest=_microscopic_column_provenance_QUERY_DIGEST,
        query_order=0,
        query_irrep_label=(block or _block()).global_irrep,
        source_parent_sg=62,
        source_k_coordinates=("0", "1/2", "0"),
        source_subgroup_space_group_number=1,
        source_subgroup_basis=np.eye(3, dtype=int).tolist(),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block or _block(),
        column_index=column,
    )

def _microscopic_column_provenance_subgroup(*, irrep: str = "GM1+") -> _microscopic_column_provenance_SubgroupInfo:
    return _microscopic_column_provenance_SubgroupInfo(
        index=0,
        space_group_number=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label=irrep,
        parent_sg=62,
        k_coordinates=["0", "0", "0"],
        basis_vectors=np.eye(3).tolist(),
        origin=[0.0, 0.0, 0.0],
    )

def _exact_fixed_space_wrapper() -> _microscopic_column_provenance_IsoWrapper:
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._direction_transport_for_invariant = lambda *_args: SimpleNamespace(
        symbol="P1",
        domain_number=1,
    )
    wrapper.list_invariant_directions = lambda *_args, **_kwargs: [
        InvariantDirection(
            "GM1+", "(a,b,c)", 1, "P1", 1,
            k_coordinates=("0", "0", "0"),
        )
    ]
    wrapper._run_session = lambda _commands: (
        "********Irrep (ML) Wyckoff Irrep Point Projected Vectors\n"
        "GM1+ a A (0,0,0) (1,0,0) (0,1,0) (0,0,1)\n"
    )
    return wrapper

def _exact_fixed_space_subgroup() -> _microscopic_column_provenance_SubgroupInfo:
    subgroup = _microscopic_column_provenance_SubgroupInfo(
        index=0,
        space_group_number=1,
        space_group_symbol="P1",
        subgroup_index=1,
        size=1,
        parent_sg=1,
        basis_vectors=np.eye(3, dtype=int).tolist(),
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin=[0, 0, 0],
        origin_raw="(0,0,0)",
    )
    subgroup.method3_route_resolution = "exact_fixed_space"
    return subgroup

def _mode_from_block(block: MicroscopicVectorBlock, column: int) -> DistortionMode:
    component = chr(ord("a") + column)
    provenance = _provenance(column, block=block)
    bushes = [
        BushMode(
            irrep_label=block.global_irrep,
            opd_symbol="P1",
            wyckoff_letter=block.wyckoff_letter,
            point=[float(Fraction(value)) for value in row.point_raw],
            point_raw=list(row.point_raw),
            displacements=[
                [float(value) for value in row.displacements[column]]
            ],
        )
        for row in block.rows
    ]
    return DistortionMode(
        irrep_label=block.global_irrep,
        dimension=block.column_count,
        wyckoff_site=block.wyckoff_letter,
        wyckoff_orbit_id="orbit-a",
        opd_symbol="P1",
        opd_dir_raw=(
            "(" + ",".join(
                chr(ord("a") + index) for index in range(block.column_count)
            ) + ")"
        ),
        bush_modes=bushes,
        amplitude_key=f"{block.global_irrep}__{component}",
        site_irrep=block.site_irrep,
        k_coords_label="0,1/2,0",
        opd_component=component,
        mode_identity=ModeIdentity(
            parent_sg=62,
            global_irrep=block.global_irrep,
            k_coordinates=("0", "1/2", "0"),
            wyckoff_letter=block.wyckoff_letter,
            orbit_id="orbit-a",
            site_irrep=block.site_irrep,
            component_index=column,
            component_label=component,
            source="iso_microscopic",
            status="partial",
            reason="awaiting_bush_or_smodes_subspace_validation",
        ),
        microscopic_provenance=provenance,
    )

def _current_mode(vector: tuple[float, float, float], index: int) -> DistortionMode:
    return DistortionMode(
        irrep_label="X1+",
        dimension=2,
        wyckoff_site="a",
        bush_modes=[
            BushMode(
                irrep_label="X1+",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[0.0, 0.0, 0.0],
                point_raw=["0", "0", "0"],
                displacements=[list(vector)],
            ),
            BushMode(
                irrep_label="X1+",
                opd_symbol="P1",
                wyckoff_letter="a",
                point=[1.0, 0.0, 0.0],
                point_raw=["1", "0", "0"],
                displacements=[[-value for value in vector]],
            ),
        ],
        amplitude_key=f"current-{index}",
    )

def _sparse_domain_extension_fixture():
    block = MicroscopicVectorBlock(
        global_irrep="X1+",
        wyckoff_letter="a",
        site_irrep="A",
        rows=(MicroscopicVectorRow(
            point_raw=("0", "0", "0"),
            displacements=(
                (Fraction(1), Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(1), Fraction(0)),
            ),
        ),),
        source_order=0,
        direction_symbol="VECTOR,A,B",
    )
    canonical = [_mode_from_block(block, 0), _mode_from_block(block, 1)]
    reference = [
        _current_mode((1.0, 1.0, 0.0), 0),
        _current_mode((1.0, -1.0, 0.0), 1),
    ]
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, canonical, parent_sg=62,
    )
    return parent, reference, canonical, selected, diagnostics

# Original coverage: test_partial_microscopic_mapping.py.

_partial_microscopic_mapping_QUERY_DIGEST = hashlib.sha256(
    b"synthetic exact secondary-irrep mapping query"
).hexdigest()

def _active_arm_mapping_fixture(
    *,
    irrep: str,
    k_coordinates: tuple[str, str, str],
    direction: str,
    active_k_coordinates: tuple[tuple[str, str, str], ...],
    basis: list[list[float]],
    points: tuple[tuple[str, str, str], ...],
    arrows: tuple[tuple[int, int, int], ...],
) -> tuple[Structure, list[dict], SubgroupInfo, _partial_microscopic_mapping_DistortionMode, list[list[float]]]:
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    orbit_id = "hydrogen-a"
    sites = [{
        "wyckoff_letter": "a",
        "orbit_id": orbit_id,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=1,
        subgroup_index=round(abs(np.linalg.det(np.asarray(basis)))),
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label=irrep,
        k_point_label=irrep.rstrip("0123456789+-"),
        k_coordinates=list(k_coordinates),
        parent_sg=123,
        basis_vectors=basis,
        basis_raw=",".join(
            "(" + ",".join(str(int(value)) for value in row) + ")"
            for row in basis
        ),
        origin=[0.0, 0.0, 0.0],
        origin_raw="(0,0,0)",
    )
    rows = tuple(
        MicroscopicVectorRow(
            point_raw=point,
            displacements=(tuple(Fraction(value) for value in arrow),),
        )
        for point, arrow in zip(points, arrows, strict=True)
    )
    block = MicroscopicVectorBlock(
        global_irrep=irrep,
        direction_symbol=direction,
        wyckoff_letter="a",
        site_irrep="A1",
        source_order=0,
        rows=rows,
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=_partial_microscopic_mapping_QUERY_DIGEST,
        query_order=0,
        query_irrep_label=irrep,
        source_parent_sg=123,
        source_k_coordinates=k_coordinates,
        source_active_k_coordinates=active_k_coordinates,
        source_subgroup_space_group_number=1,
        source_subgroup_basis=basis,
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label=irrep,
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    bushes = [
        _partial_microscopic_mapping_BushMode(
            irrep_label=irrep,
            opd_symbol="P1",
            wyckoff_letter="a",
            point=[float(Fraction(value)) for value in point],
            point_raw=list(point),
            displacements=[[float(value) for value in arrow]],
        )
        for point, arrow in zip(points, arrows, strict=True)
    ]
    identity = ModeIdentity(
        parent_sg=123,
        global_irrep=irrep,
        k_coordinates=k_coordinates,
        wyckoff_letter="a",
        orbit_id=orbit_id,
        site_irrep="A1",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="partial",
        reason="awaiting_bush_or_smodes_subspace_validation",
    )
    mode = _partial_microscopic_mapping_DistortionMode(
        irrep_label=irrep,
        dimension=1,
        mode_type="displacive",
        basis_vectors=[list(arrows[0])],
        wyckoff_site="a",
        wyckoff_orbit_id=orbit_id,
        k_point_label=subgroup.k_point_label,
        opd_symbol="P1",
        opd_dir_raw="(" + direction.split(",", 1)[1] + ")",
        bush_modes=bushes,
        amplitude_key=f"{irrep}__a__A1(a)",
        site_irrep="A1",
        k_coords_label=",".join(k_coordinates),
        opd_component="a",
        mode_identity=identity,
        microscopic_provenance=provenance,
    )
    return parent, sites, subgroup, mode, basis

# Original coverage: test_site_irreps.py.

ISODISTORT_ROOT = TESTS_ROOT.parent
