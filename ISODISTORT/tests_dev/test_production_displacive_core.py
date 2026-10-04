"""Production assembly of the shared nonempty displacive export contract."""

from __future__ import annotations

import hashlib
from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import spglib
from pymatgen.core import Lattice, Structure

from isocore.api import IsoDistort
from isocore.backend import BushMode, DistortionMode, SubgroupInfo
from isocore.backend.iso_mode_models import (
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
)
from isocore.distortion import DistortionMapper, embedding_from_identity, parent_affine_group
from isocore.distortion.affine_embeddings import AffineOperation
from isocore.io.cif_displacive_model import ExactSeitzOperation
from isocore.io.isodistort_cif import _hall_for_number
from isocore.utils import get_config
from isocore.utils.lattice import as_fraction, inverse, multiply, rational_matrix, transpose


def _matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    )


def _add(left, right):
    return tuple(left[index] + right[index] for index in range(3))


def _subtract(left, right):
    return tuple(left[index] - right[index] for index in range(3))


def test_exact_child_operations_conjugate_nondiagonal_basis_and_nonzero_origin():
    hall_number, _symbol = _hall_for_number(3)
    database = spglib.get_symmetry_from_database(hall_number)
    basis = rational_matrix(((0, 1, 0), (1, 1, 0), (0, 0, -1)))
    origin = (Fraction(1, 3), Fraction(1, 5), Fraction(1, 7))
    parent_from_child = transpose(basis)
    child_from_parent = inverse(parent_from_child)
    parent_operations = []
    expected = []
    for rotation, translation in zip(
        database["rotations"],
        database["translations"],
        strict=True,
    ):
        child_rotation = rational_matrix(rotation.tolist())
        child_translation = tuple(as_fraction(float(value)) for value in translation)
        parent_rotation = multiply(
            multiply(parent_from_child, child_rotation),
            child_from_parent,
        )
        parent_translation = _add(
            _matrix_vector(parent_from_child, child_translation),
            _subtract(origin, _matrix_vector(parent_rotation, origin)),
        )
        parent_operations.append(AffineOperation(parent_rotation, parent_translation))
        expected.append(
            ExactSeitzOperation.from_values(child_rotation, child_translation)
        )
    subgroup = SubgroupInfo(
        index=0,
        parent_sg=3,
        space_group_number=3,
        space_group_symbol="P2",
        subgroup_index=1,
        size=1,
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM1+",
        opd_symbol="P1",
        basis_vectors=[list(row) for row in basis],
        origin=list(origin),
    )
    affine_child = SimpleNamespace(
        hall_number=hall_number,
        operations=tuple(parent_operations),
    )

    actual = IsoDistort._exact_child_operations(subgroup, affine_child)

    assert {operation.key for operation in actual} == {
        operation.key for operation in expected
    }


def _p1_api_with_canonical_modes() -> tuple[IsoDistort, SubgroupInfo]:
    parent = Structure(
        Lattice.from_parameters(4.1, 5.2, 6.3, 77.0, 83.0, 71.0),
        ["H", "He"],
        [[0.1, 0.2, 0.3], [0.37, 0.41, 0.53]],
    )
    sites = [
        {
            "orbit_id": "orbit-h",
            "wyckoff_letter": "a",
            "species": "H",
            "representative_index": 0,
            "equivalent_indices": [0],
            "multiplicity": 1,
            "display_label": "H1",
            "display_order": 0,
        },
        {
            "orbit_id": "orbit-he",
            "wyckoff_letter": "a",
            "species": "He",
            "representative_index": 1,
            "equivalent_indices": [1],
            "multiplicity": 1,
            "display_label": "He1",
            "display_order": 1,
        },
    ]
    subgroup = SubgroupInfo(
        index=0,
        parent_sg=1,
        space_group_number=1,
        space_group_symbol="P1",
        subgroup_index=1,
        size=1,
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a,b,c)",
        opd_vector=[1],
        basis_vectors=[[1, 1, 0], [0, 1, 0], [0, 0, 1]],
        origin=[Fraction(1, 3), Fraction(1, 5), Fraction(1, 7)],
    )
    exact_vectors = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    block = MicroscopicVectorBlock(
        "GM1+",
        "a",
        "A",
        (
            MicroscopicVectorRow(
                point_raw=("x", "y", "z"),
                displacements=exact_vectors,
            ),
        ),
        source_order=0,
        direction_symbol="VECTOR,A,B,C",
    )
    digest = hashlib.sha256(b"synthetic-p1-production-contract").hexdigest()
    modes = []
    for orbit_id in ("orbit-h", "orbit-he"):
        for component_index, (component, vector) in enumerate(
            zip(("a", "b", "c"), exact_vectors, strict=True)
        ):
            provenance = MicroscopicColumnProvenance.from_exact_column(
                query_digest=digest,
                query_order=0,
                query_irrep_label="GM1+",
                source_parent_sg=1,
                source_k_coordinates=("0", "0", "0"),
                source_subgroup_space_group_number=1,
                source_subgroup_basis=subgroup.basis_vectors,
                source_subgroup_origin=subgroup.origin,
                source_subgroup_irrep_label="GM1+",
                source_subgroup_opd_symbol="P1",
                source_subgroup_primary_direction_selector="VECTOR,A,B,C",
                block=block,
                column_index=component_index,
            )
            identity = ModeIdentity(
                parent_sg=1,
                global_irrep="GM1+",
                k_coordinates=("0", "0", "0"),
                wyckoff_letter="a",
                orbit_id=orbit_id,
                site_irrep="A",
                component_index=component_index,
                component_label=component,
                source="iso_microscopic",
                status="verified",
            )
            key = f"GM1+[0,0,0]__a@{orbit_id}__A({component})"
            modes.append(
                DistortionMode(
                    irrep_label="GM1+",
                    dimension=3,
                    mode_type="displacive",
                    wyckoff_site="a",
                    k_point_label="GM",
                    opd_symbol="P1",
                    opd_dir_raw="(a,b,c)",
                    bush_modes=[
                        BushMode(
                            irrep_label="GM1+",
                            opd_symbol="P1",
                            wyckoff_letter="a",
                            point=[0.0, 0.0, 0.0],
                            point_raw=["x", "y", "z"],
                            displacements=[
                                [float(value) for value in vector]
                            ],
                        )
                    ],
                    amplitude_key=key,
                    site_irrep="A",
                    k_coords_label="0,0,0",
                    opd_component=component,
                    wyckoff_orbit_id=orbit_id,
                    mode_identity=identity,
                    microscopic_provenance=provenance,
                )
            )

    api = object.__new__(IsoDistort)
    api.structure = parent
    api.structure_path = None
    api.symmetry_info = {
        "space_group_number": 1,
        "space_group_symbol": "P1",
        "wyckoff_sites": sites,
    }
    api.cfg = get_config()
    api._dist_mapper = DistortionMapper()
    api.distortion_modes = modes
    api.mode_displacements = {}
    api.mode_occupancies = {}
    api.mode_displacements_sc = {}
    api._mode_label_overrides = {}
    api.phase_path = None
    api._selected_subgroup = subgroup
    api.distortion_types = ["displacive"]
    return api, subgroup


def test_production_contract_preserves_nonzero_amplitude_across_origin_frame():
    api, subgroup = _p1_api_with_canonical_modes()
    embedding = embedding_from_identity(subgroup, parent_affine_group(api.structure))
    base = api._supercell_for_subgroup(subgroup)
    atom_ids = tuple(f"test-{index}" for index in range(len(base)))
    mapped = api._dist_mapper.map_microscopic_columns_to_supercell(
        api.structure,
        api.symmetry_info["wyckoff_sites"],
        api.distortion_modes,
        subgroup.basis_vectors,
        subgroup_context=subgroup,
        frame_id="test-frame",
        atom_ids=atom_ids,
        subgroup_operations=embedding.operations,
        subgroup_translation_lattice=[
            [float(value) for value in row] for row in embedding.lattice
        ],
    )
    lifted = {column.amplitude_key: column.displacements for column in mapped}
    selected_key = api.distortion_modes[0].amplitude_key
    final = Structure(
        base.lattice,
        base.species,
        (
            np.asarray(base.frac_coords, dtype=float)
            + 0.125 * lifted[selected_key]
        )
        % 1.0,
        coords_are_cartesian=False,
    )

    data = api._displacive_export_data_for_subgroup(
        subgroup,
        lifted,
        final_structure=final,
    )

    assert data is not None
    assert data.subgroup_identity.target_hall_number == embedding.hall_number
    assert len(data.model.canonical_modes) == 6
    assert np.count_nonzero(np.abs(np.asarray(data.model.amplitudes)) > 1.0e-10) == 1
    assert data.final_structure is not None
    assert np.allclose(
        np.asarray(data.final_structure.frac_coords, dtype=float),
        (
            np.asarray(final.frac_coords, dtype=float)
            - np.asarray([float(value) for value in subgroup.origin])
            @ np.linalg.inv(np.asarray(subgroup.basis_vectors, dtype=float))
        )
        % 1.0,
        atol=1.0e-12,
    )
