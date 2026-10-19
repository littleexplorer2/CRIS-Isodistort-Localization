"""Tests for backend.method2; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.export.isodistort_isoviz import _isoviz_mode_label
from backend.input_cif import build_supercell
from backend.input_cif import read_cif
from backend.input_cif.symmetry_validator import SymmetryValidator
from backend.method2 import DistortionMapper as _distortion_DistortionMapper
from backend.method2 import Method2Query
from backend.method2 import Method2Solver
from backend.method2 import superspace as superspace_module
from backend.method2.distortion_mapper import DistortionMapper as _microscopic_column_provenance_DistortionMapper
from backend.method2.distortion_mapper import DistortionMapper as _partial_microscopic_mapping_DistortionMapper
from backend.method2.distortion_mapper import MappedMicroscopicColumn
from backend.method2.distortion_mapper import _has_exact_unit_phase
from backend.method2.distortion_mapper import _make_exact_translation_lattice
from backend.method2.distortion_mapper import _snap_to_exact_translation
from backend.method2.distortion_mapper import verify_mapped_microscopic_columns
from backend.method2.superspace import ParametricModeResult
from backend.method2.superspace import _KCandidate
from backend.method2.superspace import _child_origin_candidates
from backend.method2.superspace import _collect_smodes_blocks
from backend.method2.superspace import _common_point_rows
from backend.method2.superspace import _displacement_actions
from backend.method2.superspace import _display_site_label
from backend.method2.superspace import _extend_microscopic_modes_to_reference_domain
from backend.method2.superspace import _harmonic_order
from backend.method2.superspace import _instantiate_bush_modes_by_orbit
from backend.method2.superspace import _is_self_conjugate_parent_k
from backend.method2.superspace import _linearly_independent_mode_items
from backend.method2.superspace import _map_smodes_block_to_supercell
from backend.method2.superspace import _match_candidate_points_to_reference
from backend.method2.superspace import _ordered_complete_microscopic_space
from backend.method2.superspace import _parent_little_group_orbit
from backend.method2.superspace import _project_invariant
from backend.method2.superspace import _proven_real_transport_phase
from backend.method2.superspace import _quadrature_supercell_shift
from backend.method2.superspace import _replace_rootless_supplement_basis
from backend.method2.superspace import _select_verified_microscopic_modes
from backend.method2.superspace import _symbolic_subspace_proof
from backend.method2.superspace import _unresolved_modes
from backend.method2.superspace import _with_unresolved_identity_reason
from backend.method2.superspace import compute_parametric_modes
from backend.method2.superspace import folding_k_candidates
from backend.method2.superspace import format_k_token
from backend.method2.superspace import is_harmonic_of
from backend.method2.superspace import opd_free_components
from backend.method2.superspace import smodes_k_spec
from backend.method2.superspace import validate_symbolic_mode_subspaces
from backend.models.iso_mode_models import InvariantDirection
from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import MicroscopicVectorBlock
from backend.models.iso_mode_models import MicroscopicVectorRow
from backend.models.iso_mode_models import ModeIdentity
from backend.models.iso_mode_models import microscopic_atom_order_id
from backend.utils.k_lattice import k_compatible_with_supercell
from backend.utils.lattice import centering_primitive_matrix
from backend.utils.opd_format import _g_allowed
from backend.utils.opd_format import _k_equivalent_exact
from backend.utils.opd_format import k_star_fraction_vectors
from backend.utils.parent_header import parent_wyckoff_display
from backend.wrappers import BushMode as _distortion_BushMode
from backend.wrappers import BushMode as _search_methods_BushMode
from backend.wrappers import DistortionMode as _distortion_DistortionMode
from backend.wrappers import DistortionMode as _search_methods_DistortionMode
from backend.wrappers.iso_wrapper import BushMode as _partial_microscopic_mapping_BushMode
from backend.wrappers.iso_wrapper import BushMode as _site_irreps_BushMode
from backend.wrappers.iso_wrapper import DistortionMode as _partial_microscopic_mapping_DistortionMode
from backend.wrappers.iso_wrapper import DistortionMode as _site_irreps_DistortionMode
from backend.wrappers.iso_wrapper import DistortionMode as _superspace_DistortionMode
from backend.wrappers.iso_wrapper import IsoWrapper as _microscopic_column_provenance_IsoWrapper
from backend.wrappers.iso_wrapper import IsoWrapper as _site_irreps_IsoWrapper
from backend.wrappers.iso_wrapper import KPointInfo
from backend.wrappers.iso_wrapper import SubgroupInfo
from backend.wrappers.iso_wrapper import SubgroupInfo as _site_irreps_SubgroupInfo
from backend.wrappers.iso_wrapper import microscopic_extension_digest
from backend.wrappers.smodes_wrapper import SmodesAtom
from backend.wrappers.smodes_wrapper import SmodesModeBlock
from backend.wrappers.smodes_wrapper import parse_smodes_mode_blocks
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from pymatgen.core import Lattice
from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp
from tests.support_application_session import _mode_core_token
from tests.support_backend_method2 import ISODISTORT_ROOT
from tests.support_backend_method2 import _active_arm_mapping_fixture
from tests.support_backend_method2 import _block
from tests.support_backend_method2 import _current_mode
from tests.support_backend_method2 import _exact_fixed_space_subgroup
from tests.support_backend_method2 import _exact_fixed_space_wrapper
from tests.support_backend_method2 import _microscopic_column_provenance_QUERY_DIGEST
from tests.support_backend_method2 import _microscopic_column_provenance_subgroup
from tests.support_backend_method2 import _mode_from_block
from tests.support_backend_method2 import _partial_microscopic_mapping_QUERY_DIGEST
from tests.support_backend_method2 import _provenance
from tests.support_backend_method2 import _sparse_domain_extension_fixture
from tests.support_backend_method3 import _DummyIsoWrapper
from tests.support_backend_method3 import _method3_route
from tests.support_paths import TESTS_ROOT
from types import SimpleNamespace
import hashlib
import numpy as np
import pytest
import re
import shutil
import subprocess


# Original coverage: test_distortion.py.

# --- from test_distortion_mapper.py ---

def _mapper_parent() -> Structure:
    """四方母相：Eu 2a + Al 2e（(0,0,0.38)/(0,0,0.62)）。"""
    return Structure(
        Lattice.tetragonal(4.4, 11.2),
        ["Eu", "Eu", "Al", "Al"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0, 0, 0.38], [0, 0, 0.62]],
    )

def _mapper_wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 2, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0, 1]},
        {"wyckoff_letter": "e", "multiplicity": 2, "species": "Al",
         "representative_index": 2, "equivalent_indices": [2, 3]},
    ]

def _mapper_mode() -> _distortion_DistortionMode:
    mode = _distortion_DistortionMode(irrep_label="GM1+", opd_symbol="P1",
                          k_point_label="GM")
    mode.bush_modes = [
        _distortion_BushMode("GM1+", "P1", "a", [0, 0, 0], ["0", "0", "0"], [[0, 0, 0]]),
        _distortion_BushMode("GM1+", "P1", "a", [0.5, 0.5, 0.5], ["1/2", "1/2", "1/2"], [[0, 0, 0]]),
        # e 位点：两个相反位移的代表（(0,0,z) 与 (0,0,-z)）
        _distortion_BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
        _distortion_BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "-z"], [[0, 0, -1]]),
    ]
    return mode

def test_mapper_parameterized_z_points():
    """(0,0,z) 与 (0,0,-z) 必须解析为不同位置：z=0.38 -> +z，z=0.62 -> -z。"""
    mapper = _distortion_DistortionMapper()
    result = mapper.map_modes_to_atoms(_mapper_parent(), _mapper_wyckoff(), [_mapper_mode()])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[0], 0) and np.allclose(disp[1], 0)  # Eu 不动
    assert np.allclose(disp[2], [0, 0, 1])   # (0,0,0.38) -> +z
    assert np.allclose(disp[3], [0, 0, -1])  # (0,0,0.62) -> -z

def test_mapper_coordinate_token_parser_consumes_the_complete_expression():
    evaluate = _distortion_DistortionMapper._eval_token
    assert evaluate("-2a + 1", {"a": 0.25}) == pytest.approx(0.5)
    assert evaluate("x-y+1/2", {"x": 0.2, "y": 0.3}) == pytest.approx(0.4)

    for invalid in ("x*y", "x/2", "sin(x)", "x?1"):
        with pytest.raises(ValueError, match="BUSH coordinate"):
            evaluate(invalid, {"x": 0.2, "y": 0.3})
    with pytest.raises(ValueError, match="unresolved BUSH coordinate parameter"):
        evaluate("q", {})

def test_mapper_incomplete_single_rep_is_not_assumed_uniform():
    """单代表行不足以证明另一个位置具有相同位移。"""
    mode = _mapper_mode()
    mode.bush_modes = [
        _distortion_BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
    ]
    mapper = _distortion_DistortionMapper()
    with pytest.raises(ValueError, match="does not cover child atom"):
        mapper.map_modes_to_atoms(_mapper_parent(), _mapper_wyckoff(), [mode])

def test_mapper_periodic_distance_uses_true_skew_cell_minimum_image():
    """Component-wise fractional wrapping can miss the closest lattice image."""
    structure = Structure(
        Lattice.from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 10.0),
        ["Fe"],
        [[0.0, 0.0, 0.0]],
    )
    from backend.input_cif import map_supercell_sites

    first = np.asarray([0.0, 0.0, 0.0])
    second = np.asarray([0.49, 0.49, 0.0])
    component_wrapped = second - first
    component_wrapped -= np.round(component_wrapped)
    component_wrapped_distance = float(
        np.linalg.norm(component_wrapped @ structure.lattice.matrix)
    )

    child = Structure(structure.lattice, ["Fe"], [second])
    # This point has two equally close periodic images. Mapping it with a
    # radius covering both cannot determine a unique translation phase.
    with pytest.raises(ValueError, match="compatible parent sites"):
        map_supercell_sites(structure, child, np.eye(3), tolerance=1.0)
    expected = structure.lattice.get_distance_and_image(first, second)[0]
    assert expected < 1.0
    assert component_wrapped_distance > 9.0

def test_bush_supercell_mapping_preserves_explicit_phase_pattern():
    """BUSH representatives, not a guessed primary-k cosine, define phases."""
    parent = Structure(Lattice.cubic(4.0), ["H"], [[0, 0, 0]])
    sites = [{
        "wyckoff_letter": "a",
        "multiplicity": 1,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    mode = _distortion_DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="X1+__a__P1(1)__a",
        bush_modes=[
            _distortion_BushMode("X1+", "P1(1)", "a", [0, 0, 0],
                     ["0", "0", "0"], [[1, 0, 0]]),
            _distortion_BushMode("X1+", "P1(1)", "a", [1, 0, 0],
                     ["1", "0", "0"], [[-1, 0, 0]]),
        ],
    )

    mapped = _distortion_DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
    )[mode.amplitude_key]

    # Parent-fractional arrows +/-1 become child-fractional +/-1/2.
    assert mapped.shape == (2, 3)
    assert np.allclose(sorted(mapped[:, 0]), [-0.5, 0.5])
    assert np.allclose(mapped[:, 1:], 0.0)

def test_bush_supercell_mapping_is_invariant_to_unreduced_translation_basis():
    """A GL(3,Z) basis change must not change nearest-image coverage."""
    parent = Structure(Lattice.cubic(1.0), ["H"], [[0, 0, 0]])
    sites = [{
        "wyckoff_letter": "a",
        "multiplicity": 1,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    mode = _distortion_DistortionMode(
        irrep_label="GM1+",
        wyckoff_site="a",
        amplitude_key="GM1+__a__P1__a",
        bush_modes=[_distortion_BushMode(
            "GM1+", "P1", "a", [-0.98, -0.049, 0],
            ["-0.98", "-0.049", "0"], [[1, 0, 0]],
        )],
    )
    # This unimodular, strongly sheared basis generates exactly Z^3.  In its
    # fractional coordinates the representative residual is
    # (0.49, 0.049, 0): component-wise rounding gives 0.981 Å, while the true
    # nearest lattice image is only sqrt(0.02^2 + 0.049^2) Å away.
    translation_basis = [[1, 0, 0], [10, 1, 0], [0, 0, 1]]

    mapped = _distortion_DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        [mode],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        cartesian_tolerance=0.1,
        subgroup_translation_lattice=translation_basis,
    )[mode.amplitude_key]

    assert np.allclose(mapped, [[1, 0, 0]])

def test_bush_supercell_mapping_rejects_incomplete_representatives():
    parent = Structure(Lattice.cubic(4.0), ["H"], [[0, 0, 0]])
    sites = [{
        "wyckoff_letter": "a",
        "multiplicity": 1,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    mode = _distortion_DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="X1+__a__P1(1)__a",
        bush_modes=[_distortion_BushMode(
            "X1+", "P1(1)", "a", [0, 0, 0],
            ["0", "0", "0"], [[1, 0, 0]],
        )],
    )

    with pytest.raises(ValueError, match="does not cover child atom"):
        _distortion_DistortionMapper().map_bush_modes_to_supercell(
            parent, sites, [mode], [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
        )

def test_zone_center_bush_rows_repeat_over_centered_parent_translations():
    """A Gamma column printed for one primitive parent cell covers a larger child."""
    parent = Structure(
        Lattice.cubic(4.0),
        ["H"] * 4,
        [[0, 0, 0.2], [0, 0, 0.8], [0.5, 0.5, 0.7], [0.5, 0.5, 0.3]],
    )
    orbit_id = "centered-orbit"
    sites = [{
        "wyckoff_letter": "e",
        "orbit_id": orbit_id,
        "equivalent_indices": [0, 1, 2, 3],
        "standard_representative_parameters": {"z": 0.2},
    }]
    mode = _distortion_DistortionMode(
        irrep_label="GM",
        wyckoff_site="e",
        wyckoff_orbit_id=orbit_id,
        amplitude_key="gamma-column",
        mode_identity=ModeIdentity(
            # The fixture explicitly contains the I-centering translation
            # (1/2,1/2,1/2); its source identity must declare that lattice.
            parent_sg=229,
            global_irrep="GM",
            k_coordinates=("0", "0", "0"),
            wyckoff_letter="e",
            orbit_id=orbit_id,
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[
            _distortion_BushMode("GM", "P1", "e", [0, 0, 0],
                     ["0", "0", "z"], [[0, 0, 1]]),
            _distortion_BushMode("GM", "P1", "e", [0, 0, 0],
                     ["0", "0", "-z"], [[0, 0, -1]]),
        ],
    )
    basis = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mapped = _distortion_DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], basis,
    )[mode.amplitude_key]
    child = build_supercell(parent, basis)
    assert mapped.shape == (len(child), 3)
    for index, site in enumerate(child):
        parent_z = (site.frac_coords @ np.asarray(basis))[2] % 1
        expected = 0.5 if np.isclose(parent_z, (0.2, 0.7)).any() else -0.5
        assert np.allclose(mapped[index], [0, 0, expected])

def test_bush_rows_repeat_only_along_translations_invariant_under_k_star():
    """An in-plane boundary k permits c repetition but changes a/b phase."""
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    sites = [{"wyckoff_letter": "a", "equivalent_indices": [0]}]
    mode = _distortion_DistortionMode(
        irrep_label="X",
        wyckoff_site="a",
        amplitude_key="boundary-column",
        mode_identity=ModeIdentity(
            parent_sg=123,
            global_irrep="X",
            k_coordinates=("1/2", "1/2", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[
            _distortion_BushMode("X", "P1", "a", [0, 0, 0],
                     ["0", "0", "0"], [[1, 0, 0]]),
            _distortion_BushMode("X", "P1", "a", [1, 0, 0],
                     ["1", "0", "0"], [[-1, 0, 0]]),
            _distortion_BushMode("X", "P1", "a", [0, 1, 0],
                     ["0", "1", "0"], [[-1, 0, 0]]),
            _distortion_BushMode("X", "P1", "a", [1, 1, 0],
                     ["1", "1", "0"], [[1, 0, 0]]),
        ],
    )
    basis = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mapped = _distortion_DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], basis,
    )[mode.amplitude_key]
    child = build_supercell(parent, basis)
    for index, site in enumerate(child):
        x, y, _z = site.frac_coords @ np.asarray(basis)
        expected = 0.5 * (-1) ** (round(x) + round(y))
        assert np.allclose(mapped[index], [expected, 0, 0])

def test_explicit_multiarm_bush_rows_can_select_one_star_arm():
    """Complete BUSH rows can prove that one arm of a two-arm star is active."""
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    sites = [{"wyckoff_letter": "a", "equivalent_indices": [0]}]
    mode = _distortion_DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="single-arm-column",
        opd_symbol="C1",
        opd_component="a",
        mode_identity=ModeIdentity(
            parent_sg=123,
            global_irrep="X1+",
            k_coordinates=("0", "1/2", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[
            _distortion_BushMode("X1+", "C1", "a", [0, 0, 0],
                     ["0", "0", "0"], [[1, 0, 0]]),
            _distortion_BushMode("X1+", "C1", "a", [1, 0, 0],
                     ["1", "0", "0"], [[-1, 0, 0]]),
            _distortion_BushMode("X1+", "C1", "a", [0, 1, 0],
                     ["0", "1", "0"], [[1, 0, 0]]),
            _distortion_BushMode("X1+", "C1", "a", [1, 1, 0],
                     ["1", "1", "0"], [[-1, 0, 0]]),
        ],
    )
    basis = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
    mapped = _distortion_DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], basis,
    )[mode.amplitude_key]
    child = build_supercell(parent, basis)
    for index, site in enumerate(child):
        x, _y, _z = site.frac_coords @ np.asarray(basis)
        assert np.allclose(mapped[index], [0.5 * (-1) ** round(x), 0, 0])

def test_sparse_multiarm_bush_without_column_arm_weights_fails_closed():
    """A star identity and local column label do not identify its active arm.

    At the one printed translation, both legitimate X-star completions
    ``u=(-1)^x e_x`` and ``u=(-1)^y e_x`` have the same arrow.  The C1 OPD
    name and local ``a`` component do not bind this microscopic column to
    either arm, so repeating it would invent missing source information.
    """
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    mode = _distortion_DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="ambiguous-single-arm-column",
        opd_symbol="C1",
        opd_component="a",
        mode_identity=ModeIdentity(
            parent_sg=123,
            global_irrep="X1+",
            k_coordinates=("0", "1/2", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[_distortion_BushMode(
            "X1+", "C1", "a", [0, 0, 0],
            ["0", "0", "0"], [[1, 0, 0]],
        )],
    )
    with pytest.raises(ValueError, match="does not cover child atom"):
        _distortion_DistortionMapper().map_bush_modes_to_supercell(
            parent,
            [{"wyckoff_letter": "a", "equivalent_indices": [0]}],
            [mode],
            [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
        )

def test_integer_conventional_k_does_not_hide_centering_phase():
    """An integer conventional k can still be nonzero on a centered lattice."""
    parent = Structure(
        Lattice.cubic(4.0), ["H", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    mode = _distortion_DistortionMode(
        irrep_label="boundary",
        wyckoff_site="a",
        amplitude_key="centered-boundary",
        mode_identity=ModeIdentity(
            parent_sg=229,
            global_irrep="boundary",
            k_coordinates=("1", "0", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[_distortion_BushMode(
            "boundary", "P1", "a", [0, 0, 0],
            ["0", "0", "0"], [[1, 0, 0]],
        )],
    )
    with pytest.raises(ValueError, match="does not cover child atom"):
        _distortion_DistortionMapper().map_bush_modes_to_supercell(
            parent,
            [{"wyckoff_letter": "a", "equivalent_indices": [0, 1]}],
            [mode],
            [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        )

def test_bush_mapping_uses_standard_representative_for_one_physical_orbit():
    """A centered equivalent input point must use the standard free parameter."""
    parent = Structure(
        Lattice.cubic(4.0),
        ["H", "He"],
        [[0.5, 0.5, 0.25], [0.0, 0.0, 0.1]],
    )
    sites = [
        {
            "orbit_id": "hydrogen-e",
            "wyckoff_letter": "e",
            "multiplicity": 1,
            "species": "H",
            "representative_index": 0,
            "equivalent_indices": [0],
            "standard_representative_frac_coords": [0.0, 0.0, -0.25],
            "standard_representative_symmform": "0,0,z",
            "standard_representative_parameters": {"z": -0.25},
        },
        {
            "orbit_id": "helium-e",
            "wyckoff_letter": "e",
            "multiplicity": 1,
            "species": "He",
            "representative_index": 1,
            "equivalent_indices": [1],
        },
    ]
    mode = _distortion_DistortionMode(
        irrep_label="N1+",
        wyckoff_site="e",
        wyckoff_orbit_id="hydrogen-e",
        amplitude_key="N1+__e@hydrogen-e__P1__a",
        bush_modes=[_distortion_BushMode(
            "N1+", "P1", "e", [0, 0, 0],
            ["0", "0", "z"], [[1, 0, 0]],
        )],
    )

    mapped = _distortion_DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        [mode],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        subgroup_translation_lattice=[
            [1, 0, 0], [0, 1, 0], [0.5, 0.5, 0.5],
        ],
    )[mode.amplitude_key]

    assert np.allclose(mapped[0], [1, 0, 0])
    assert np.allclose(mapped[1], [0, 0, 0])

# Original coverage: test_microscopic_column_provenance.py.

def test_parametric_primary_direction_requires_independent_subspace_match():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.list_invariant_directions = lambda *_args, **_kwargs: [
        InvariantDirection(
            "LD5", "(0,a)", 1, "P1", 1,
            k_parameters=("1/12",), k_coordinates=("0", "0", "1/6"),
        )
    ]
    subgroup = _microscopic_column_provenance_subgroup(irrep="LD5")
    subgroup.k_coordinates = ["0", "0", "1/6"]
    subgroup.opd_dir_raw = "(a,0)"

    with pytest.raises(ValueError, match="disagrees with the selected subgroup"):
        _microscopic_column_provenance_IsoWrapper.calc_parametric_microscopic_distortion_modes(
            wrapper,
            139,
            subgroup,
            ["a"],
            primary_k_parameters=("1/12",),
        )

def test_exact_fixed_space_compute_installs_verified_canonical_columns():
    parent = Structure(Lattice.cubic(4.0), ["H"], [[0, 0, 0]])
    smodes_blocks = [
        SmodesModeBlock(
            irrep="GM1+",
            degeneracy=1,
            k_label="GM",
            k_coords=("0", "0", "0"),
            superlattice_cart=np.eye(3) * 4.0,
            atoms=[SmodesAtom(1, "H", np.zeros(3))],
            atom_indices=[1],
            displacements_cart=np.asarray([vector], dtype=float),
        )
        for vector in np.eye(3)
    ]

    class _Smodes:
        @staticmethod
        def compute_mode_blocks(*_args, **_kwargs):
            return smodes_blocks

    result = compute_parametric_modes(
        parent,
        {
            "space_group_number": 1,
            "space_group_symbol": "P1",
            "wyckoff_sites": [{
                "orbit_id": "orbit-a",
                "wyckoff_letter": "a",
                "species": "H",
                "representative_index": 0,
                "equivalent_indices": [0],
            }],
        },
        _exact_fixed_space_subgroup(),
        ["a"],
        _Smodes(),
        iso=_exact_fixed_space_wrapper(),
        kpoints=[KPointInfo("GM", ["0", "0", "0"], [], True)],
        nmod=0,
        wyckoff_orbit_ids=["orbit-a"],
    )

    assert len(result.modes) == 3
    assert len(result.supercell_displacements) == 3
    assert all(mode.mode_identity.status == "verified" for mode in result.modes)
    assert all(
        mode.microscopic_provenance.source_subgroup_context_kind
        == "exact_fixed_space"
        for mode in result.modes
    )
    assert "complete exact-fixed-space ISO basis installed" in result.note

def test_exact_fixed_space_compute_canonicalizes_non_gamma_folded_space():
    parent = Structure(Lattice.cubic(4.0), ["H"], [[0, 0, 0]])
    subgroup = _exact_fixed_space_subgroup()
    subgroup.size = 2
    subgroup.basis_vectors = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    subgroup.basis_raw = "(2,0,0),(0,1,0),(0,0,1)"

    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.list_invariant_directions = lambda *_args, **_kwargs: [
        InvariantDirection(
            "GM1+", "(a,b,c)", 1, "P1", 1,
            k_coordinates=("0", "0", "0"),
        ),
        InvariantDirection(
            "X1", "(a,b,c)", 1, "P1", 2,
            k_coordinates=("1/2", "0", "0"),
        ),
    ]
    wrapper._direction_transport_for_invariant = lambda *_args: SimpleNamespace(
        symbol="P1",
        domain_number=1,
    )
    wrapper._run_session = lambda _commands: (
        "********Irrep (ML) Wyckoff Irrep Point Projected Vectors\n"
        "GM1+ a A (0,0,0) (1,0,0) (0,1,0) (0,0,1)\n"
        "********Irrep (ML) Wyckoff Irrep Point Projected Vectors\n"
        "X1 a A (0,0,0) (1,0,0) (0,1,0) (0,0,1)\n"
        "          (1,0,0) (-1,0,0) (0,-1,0) (0,0,-1)\n"
    )

    gamma_atom = SmodesAtom(1, "H", np.zeros(3))
    x_atoms = [
        SmodesAtom(1, "H", np.asarray([0.0, 0.0, 0.0])),
        SmodesAtom(2, "H", np.asarray([4.0, 0.0, 0.0])),
    ]
    smodes_blocks = []
    for irrep, k_label, k_coords in (
        ("GM1+", "GM", ("0", "0", "0")),
        ("X1", "X", ("1/2", "0", "0")),
    ):
        for vector in np.eye(3):
            if k_label == "GM":
                atoms = [gamma_atom]
                atom_indices = [1]
                displacements = np.asarray([vector], dtype=float)
                superlattice = np.eye(3) * 4.0
            else:
                atoms = x_atoms
                atom_indices = [1, 2]
                displacements = np.asarray([vector, -vector], dtype=float)
                superlattice = np.diag([8.0, 4.0, 4.0])
            smodes_blocks.append(SmodesModeBlock(
                irrep=irrep,
                degeneracy=1,
                k_label=k_label,
                k_coords=k_coords,
                superlattice_cart=superlattice,
                atoms=atoms,
                atom_indices=atom_indices,
                displacements_cart=displacements,
            ))

    class _Smodes:
        @staticmethod
        def compute_mode_blocks(*_args, **_kwargs):
            return smodes_blocks

    result = compute_parametric_modes(
        parent,
        {
            "space_group_number": 1,
            "space_group_symbol": "P1",
            "wyckoff_sites": [{
                "orbit_id": "orbit-a",
                "wyckoff_letter": "a",
                "species": "H",
                "representative_index": 0,
                "equivalent_indices": [0],
            }],
        },
        subgroup,
        ["a"],
        _Smodes(),
        iso=wrapper,
        kpoints=[
            KPointInfo("GM", ["0", "0", "0"], [], True),
            KPointInfo("X", ["1/2", "0", "0"], [], True),
        ],
        nmod=0,
        wyckoff_orbit_ids=["orbit-a"],
    )

    assert len(result.modes) == 6
    assert len(result.supercell_displacements) == 6
    assert {
        mode.mode_identity.k_coordinates for mode in result.modes
    } == {("0", "0", "0"), ("1/2", "0", "0")}
    assert all(mode.mode_identity.status == "verified" for mode in result.modes), (
        result.note,
        [
            (mode.irrep_label, mode.mode_identity.status, mode.mode_identity.reason)
            for mode in result.modes
        ],
    )
    assert "complete exact-fixed-space ISO basis installed" in result.note

def _mapped(
    component_label: str,
    provenance: MicroscopicColumnProvenance,
    displacements: np.ndarray,
    *,
    canonical: bool,
    frame_id: str = "child-frame",
    atom_ids: tuple[str, ...] = ("a", "b"),
    amplitude_key: str | None = None,
    mode_identity: ModeIdentity | None = None,
) -> MappedMicroscopicColumn:
    bound_provenance = provenance.bind_mapping(
        frame_id=frame_id,
        atom_ids=atom_ids,
    )
    identity = mode_identity or ModeIdentity(
        parent_sg=62,
        global_irrep=provenance.source_global_irrep,
        k_coordinates=("0", "1/2", "0"),
        wyckoff_letter=provenance.source_wyckoff_letter,
        orbit_id="orbit-a",
        site_irrep=provenance.source_site_irrep,
        component_index=provenance.source_column_index,
        component_label=component_label,
        source="iso_microscopic",
        status="verified",
    )
    return MappedMicroscopicColumn(
        mode_id=identity.stable_token,
        amplitude_key=amplitude_key or f"amp-{component_label}",
        frame_id=frame_id,
        atom_order_id=microscopic_atom_order_id(atom_ids),
        atom_ids=atom_ids,
        displacements=displacements,
        mode_identity=identity,
        provenance=bound_provenance,
        relation_status="canonical" if canonical else "unverified",
        current_over_canonical_signed_scale=1.0 if canonical else None,
        cartesian_residual_angstrom=0.0 if canonical else None,
    )

def test_superspace_installs_canonical_iso_basis_after_nondiagonal_change():
    block = _block()
    canonical = [_mode_from_block(block, 0), _mode_from_block(block, 1)]
    mixed = [
        _current_mode((1.0, 1.0, 0.0), 0),
        _current_mode((1.0, -1.0, 0.0), 1),
    ]
    parent = Structure(Lattice.cubic(4.0), ["Fe", "Fe"], [[0, 0, 0], [0, 0, 0]])
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, mixed, list(reversed(canonical)), parent_sg=62
    )
    assert diagnostics[("X1+", "a")].matched
    assert [mode.amplitude_key for mode in selected] == [
        canonical[0].amplitude_key,
        canonical[1].amplitude_key,
    ]
    assert [mode.microscopic_provenance for mode in selected] == [
        canonical[0].microscopic_provenance,
        canonical[1].microscopic_provenance,
    ]
    assert all(mode.mode_identity.status == "verified" for mode in selected)

    selected, diagnostics = _select_verified_microscopic_modes(
        parent, [_current_mode((1.0, 0.0, 0.0), 0)], [canonical[0]], parent_sg=62
    )
    assert diagnostics[("X1+", "a")].reason == "microscopic_source_block_column_missing"
    assert selected[0].mode_identity.status == "unresolved"

def test_sparse_canonical_rows_are_uniquely_extended_with_separate_evidence():
    _parent, _reference, canonical, selected, diagnostics = (
        _sparse_domain_extension_fixture()
    )

    assert diagnostics[("X1+", "a")].matched
    assert len(selected) == 2
    assert all(len(mode.bush_modes) == 2 for mode in selected)
    assert all(mode.mode_identity.status == "verified" for mode in selected)
    assert [mode.microscopic_provenance for mode in selected] == [
        mode.microscopic_provenance for mode in canonical
    ]
    extensions = [mode.microscopic_domain_extension for mode in selected]
    assert all(extension is not None for extension in extensions)
    assert [extension.coefficient_column for extension in extensions] == [
        pytest.approx((0.5, 0.5)),
        pytest.approx((0.5, -0.5)),
    ]
    assert all(
        extension.relation_status == "verified_bush_domain_extension"
        and extension.comparison_domain_kind == "exact_common_points"
        and extension.reference_rank == 2
        and extension.canonical_condition_number <= extension.maximum_condition_number
        and extension.matrix_convention
        == "R_common@T=C_common;C_full=R_full@T"
        for extension in extensions
    )
    assert np.allclose(
        [mode.bush_modes[1].displacements[0] for mode in selected],
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
    )
    assert all(len(mode.bush_modes) == 1 for mode in canonical)
    unresolved = _unresolved_modes(selected, parent_sg=62, reason="forced_failure")
    assert all(mode.microscopic_domain_extension is None for mode in unresolved)

def test_sparse_extension_is_invariant_to_bush_basis_column_order():
    parent, reference, canonical, baseline, _diagnostics = (
        _sparse_domain_extension_fixture()
    )

    selected, diagnostics = _select_verified_microscopic_modes(
        parent,
        list(reversed(reference)),
        list(reversed(canonical)),
        parent_sg=62,
    )

    assert diagnostics[("X1+", "a")].matched
    assert [mode.amplitude_key for mode in selected] == [
        mode.amplitude_key for mode in canonical
    ]
    assert [mode.microscopic_provenance for mode in selected] == [
        mode.microscopic_provenance for mode in canonical
    ]
    for actual, expected in zip(selected, baseline, strict=True):
        assert np.allclose(
            [bush.displacements[0] for bush in actual.bush_modes],
            [bush.displacements[0] for bush in expected.bush_modes],
        )
        assert actual.mode_identity == expected.mode_identity
        assert actual.microscopic_domain_extension is not None
        assert (
            actual.microscopic_domain_extension.canonical_source_order_key
            == actual.microscopic_provenance.source_order_key
        )

def test_domain_extension_proof_rejects_rebound_columns_and_basis():
    parent, reference, canonical, _selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    lattice = np.asarray(parent.lattice.matrix, dtype=float)
    proof = _symbolic_subspace_proof(reference, canonical, lattice)
    assert proof.validation.matched

    with pytest.raises(ValueError, match="proof_candidate_binding_mismatch"):
        _extend_microscopic_modes_to_reference_domain(
            reference, list(reversed(canonical)), lattice, proof,
        )

    rebound = replace(canonical[0], bush_modes=list(canonical[1].bush_modes))
    with pytest.raises(
        ValueError, match="proof_candidate_source_vector_digest_mismatch",
    ):
        _extend_microscopic_modes_to_reference_domain(
            reference, [rebound, canonical[1]], lattice, proof,
        )

    with pytest.raises(ValueError, match="proof_reference_basis_mismatch"):
        _extend_microscopic_modes_to_reference_domain(
            list(reversed(reference)), canonical, lattice, proof,
        )

def test_active_arm_domain_extension_maps_with_verified_derived_rows():
    """A proven full-domain extension is not a forged direct ISO row table."""

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
    source_basis = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    canonical = []
    for column in range(2):
        mode = _mode_from_block(block, column)
        provenance = MicroscopicColumnProvenance.from_exact_column(
            query_digest=_microscopic_column_provenance_QUERY_DIGEST,
            query_order=0,
            query_irrep_label=block.global_irrep,
            source_parent_sg=62,
            source_k_coordinates=("1/2", "0", "0"),
            source_active_k_coordinates=(("1/2", "0", "0"),),
            source_subgroup_space_group_number=1,
            source_subgroup_basis=source_basis,
            source_subgroup_origin=(0, 0, 0),
            source_subgroup_irrep_label="GM1+",
            source_subgroup_opd_symbol="P1",
            source_subgroup_primary_direction_selector="VECTOR,A",
            block=block,
            column_index=column,
        )
        canonical.append(replace(
            mode,
            k_coords_label="1/2,0,0",
            mode_identity=replace(
                mode.mode_identity,
                k_coordinates=("1/2", "0", "0"),
            ),
            microscopic_provenance=provenance,
        ))
    reference = [
        _current_mode((1.0, 1.0, 0.0), 0),
        _current_mode((1.0, -1.0, 0.0), 1),
    ]
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, canonical, parent_sg=62,
    )
    assert diagnostics[("X1+", "a")].matched
    assert all(mode.microscopic_domain_extension is not None for mode in selected)

    subgroup = _microscopic_column_provenance_subgroup()
    subgroup.basis_vectors = source_basis
    subgroup.basis_raw = "(2,0,0),(0,1,0),(0,0,1)"
    sites = [{
        "wyckoff_letter": "a",
        "orbit_id": "orbit-a",
        "equivalent_indices": [0],
    }]
    mapped = _microscopic_column_provenance_DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        selected,
        source_basis,
        subgroup_translation_lattice=source_basis,
        subgroup_context=subgroup,
    )

    assert set(mapped) == {mode.amplitude_key for mode in selected}
    assert all(array.shape == (2, 3) for array in mapped.values())

    damaged_rows = list(selected[0].bush_modes)
    damaged_rows[1] = replace(
        damaged_rows[1], displacements=[[9.0, 0.0, 0.0]],
    )
    with pytest.raises(ValueError, match="verified full-domain column"):
        _microscopic_column_provenance_DistortionMapper().map_bush_modes_to_supercell(
            parent,
            sites,
            [replace(selected[0], bush_modes=damaged_rows), selected[1]],
            source_basis,
            subgroup_translation_lattice=source_basis,
            subgroup_context=subgroup,
        )

def test_mapper_rejects_tampered_or_incomplete_domain_extension_evidence():
    parent, _reference, _canonical, selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    mapper = _microscopic_column_provenance_DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.ones((2, 3), dtype=float) for mode in selected
    }

    def _map(modes, subgroup=None):
        return mapper.map_microscopic_columns_to_supercell(
            parent,
            [],
            modes,
            np.eye(3).tolist(),
            subgroup_context=subgroup or _microscopic_column_provenance_subgroup(),
            frame_id="child-frame",
            atom_ids=("atom-0", "atom-1"),
        )

    first_extension = selected[0].microscopic_domain_extension
    second_extension = selected[1].microscopic_domain_extension
    assert first_extension is not None and second_extension is not None
    wrong_coefficient = replace(
        selected[0],
        microscopic_domain_extension=replace(
            first_extension,
            coefficient_column=(0.25, 0.5),
        ),
    )
    with pytest.raises(ValueError, match="change-of-basis digest"):
        _map([wrong_coefficient, selected[1]])

    wrong_output_order = replace(
        selected[1],
        microscopic_domain_extension=replace(
            second_extension,
            canonical_output_column_index=0,
        ),
    )
    with pytest.raises(ValueError, match="canonical column order"):
        _map([selected[0], wrong_output_order])

    swapped_output_order = [
        replace(
            selected[column],
            microscopic_domain_extension=replace(
                selected[column].microscopic_domain_extension,
                canonical_output_column_index=1 - column,
            ),
        )
        for column in range(2)
    ]
    with pytest.raises(ValueError, match="differs from exact source order"):
        _map(swapped_output_order)

    wrong_source_index = replace(
        selected[0],
        microscopic_domain_extension=replace(
            first_extension,
            canonical_source_order_key=(0, 9, 0),
        ),
    )
    with pytest.raises(ValueError, match="another source index"):
        _map([wrong_source_index, selected[1]])

    damaged_bush = list(selected[0].bush_modes)
    damaged_bush[0] = replace(
        damaged_bush[0], displacements=[[9.0, 0.0, 0.0]],
    )
    with pytest.raises(ValueError, match="current full column"):
        _map([replace(selected[0], bush_modes=damaged_bush), selected[1]])

    with pytest.raises(ValueError, match="source basis is incomplete"):
        _map([selected[0]])

    wrong_transport = replace(
        selected[0],
        microscopic_domain_extension=replace(
            first_extension,
            transport_digest="0" * 64,
        ),
    )
    with pytest.raises(ValueError, match="group evidence is inconsistent"):
        _map([wrong_transport, selected[1]])

    with pytest.raises(ValueError, match="numerical proof is invalid"):
        replace(first_extension, maximum_normalized_residual=1.01)

    singular_change = ((0.5, 0.5), (0.5, 0.5))
    singular_digest = microscopic_extension_digest(singular_change)
    singular_modes = [
        replace(
            mode,
            microscopic_domain_extension=replace(
                mode.microscopic_domain_extension,
                coefficient_column=tuple(
                    singular_change[row][column] for row in range(2)
                ),
                change_of_basis_digest=singular_digest,
            ),
        )
        for column, mode in enumerate(selected)
    ]
    with pytest.raises(ValueError, match="change-of-basis rank loss"):
        _map(singular_modes)

    wrong_condition = [
        replace(
            mode,
            microscopic_domain_extension=replace(
                mode.microscopic_domain_extension,
                change_of_basis_condition_number=(
                    mode.microscopic_domain_extension
                    .change_of_basis_condition_number
                    * 1.01
                ),
            ),
        )
        for mode in selected
    ]
    with pytest.raises(ValueError, match="change-of-basis condition mismatch"):
        _map(wrong_condition)

    shifted_context = replace(_microscopic_column_provenance_subgroup(), origin=[1.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="query context differs"):
        _map(selected, shifted_context)

def test_domain_extension_groups_are_independent_per_physical_orbit():
    parent, _reference, _canonical, first_orbit, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    second_orbit = [
        replace(
            mode,
            amplitude_key=f"{mode.amplitude_key}__orbit-b",
            wyckoff_orbit_id="orbit-b",
            mode_identity=replace(mode.mode_identity, orbit_id="orbit-b"),
        )
        for mode in first_orbit
    ]
    modes = [*first_orbit, *second_orbit]
    mapper = _microscopic_column_provenance_DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.ones((2, 3), dtype=float) for mode in modes
    }

    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        [],
        modes,
        np.eye(3).tolist(),
        subgroup_context=_microscopic_column_provenance_subgroup(),
        frame_id="child-frame",
        atom_ids=("atom-0", "atom-1"),
    )
    assert len(mapped) == 4
    assert len({column.mode_id for column in mapped}) == 4
    assert len({
        column.provenance.exact_source_token for column in mapped
    }) == 2

    with pytest.raises(ValueError, match="source basis is incomplete"):
        mapper.map_microscopic_columns_to_supercell(
            parent,
            [],
            [*first_orbit, second_orbit[0]],
            np.eye(3).tolist(),
            subgroup_context=_microscopic_column_provenance_subgroup(),
            frame_id="child-frame",
            atom_ids=("atom-0", "atom-1"),
        )

def test_domain_extension_rejects_rank_loss_residual_and_ill_conditioned_change():
    parent, reference, canonical, _selected, _diagnostics = (
        _sparse_domain_extension_fixture()
    )
    lattice = np.asarray(parent.lattice.matrix, dtype=float)
    proof = _symbolic_subspace_proof(reference, canonical, lattice)
    assert proof.validation.matched

    rank_lost = replace(
        proof,
        reference_comparison=(
            np.asarray([[1.0, 0.0, 0.0]]),
            np.asarray([[0.0, 0.0, 0.0]]),
        ),
    )
    with pytest.raises(ValueError, match="reference_common_rank_loss"):
        _extend_microscopic_modes_to_reference_domain(
            reference, canonical, lattice, rank_lost,
        )

    residual = replace(
        proof,
        candidate_comparison=(
            np.asarray([[1.0, 0.0, 1.0]]),
            np.asarray([[0.0, 1.0, 0.0]]),
        ),
    )
    with pytest.raises(ValueError, match="residual_exceeds_tolerance"):
        _extend_microscopic_modes_to_reference_domain(
            reference, canonical, lattice, residual,
        )

    ill_conditioned = replace(
        proof,
        reference_comparison=(
            np.asarray([[1.0, 0.0, 0.0]]),
            np.asarray([[0.0, 1.0e-5, 0.0]]),
        ),
        candidate_comparison=(
            np.asarray([[1.0e-5, 0.0, 0.0]]),
            np.asarray([[0.0, 1.0, 0.0]]),
        ),
    )
    with pytest.raises(ValueError, match="change_of_basis_ill_conditioned"):
        _extend_microscopic_modes_to_reference_domain(
            reference, canonical, lattice, ill_conditioned,
        )

def test_mapper_preserves_source_provenance_and_signed_scalars():
    block = _block()
    modes = [_mode_from_block(block, 0), _mode_from_block(block, 1)]
    modes = [
        replace(mode, mode_identity=replace(mode.mode_identity, status="verified", reason=None))
        for mode in modes
    ]
    first = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    second = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    mapper = _microscopic_column_provenance_DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        modes[0].amplitude_key: first,
        modes[1].amplitude_key: second,
    }
    source = mapper.map_microscopic_columns_to_supercell(
        Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]]),
        [],
        modes,
        np.eye(3).tolist(),
        subgroup_context=_microscopic_column_provenance_subgroup(),
        frame_id="child-frame",
        atom_ids=("a", "b"),
    )
    assert [column.mode_identity for column in source] == [
        mode.mode_identity for mode in modes
    ]
    assert [column.provenance.exact_vector_digest for column in source] == [
        mode.microscopic_provenance.exact_vector_digest for mode in modes
    ]
    assert [column.provenance.source_float_column_digest for column in source] == [
        mode.microscopic_provenance.source_float_column_digest for mode in modes
    ]
    assert all(column.provenance.unresolved_fields == () for column in source)
    assert all(
        column.provenance.source_frame_id == "child-frame"
        and column.provenance.source_atom_order_id
        == microscopic_atom_order_id(("a", "b"))
        for column in source
    )

    current = [
        _mapped(
            "b",
            source[1].provenance,
            -3.0 * second,
            canonical=False,
            amplitude_key=source[1].amplitude_key,
            mode_identity=source[1].mode_identity,
        ),
        _mapped(
            "a",
            source[0].provenance,
            2.0 * first,
            canonical=False,
            amplitude_key=source[0].amplitude_key,
            mode_identity=source[0].mode_identity,
        ),
    ]
    verified = verify_mapped_microscopic_columns(
        source,
        current,
        [[2.0, 0.3, 0.0], [0.0, 3.0, 0.2], [0.0, 0.0, 4.0]],
    )
    assert [column.mode_id for column in verified] == [
        source[0].mode_id,
        source[1].mode_id,
    ]
    assert [column.current_over_canonical_signed_scale for column in verified] == pytest.approx([2.0, -3.0])
    assert all(column.relation_status == "scalar_verified" for column in verified)

def test_mapper_rejects_tampered_direct_displacement_and_outer_routing_fields():
    block = _block()
    modes = [
        replace(
            _mode_from_block(block, column),
            mode_identity=replace(
                _mode_from_block(block, column).mode_identity,
                status="verified",
                reason=None,
            ),
        )
        for column in range(2)
    ]
    mapper = _microscopic_column_provenance_DistortionMapper()
    mapper.map_bush_modes_to_supercell = lambda *_args, **_kwargs: {
        mode.amplitude_key: np.ones((2, 3), dtype=float) for mode in modes
    }

    def _map(candidate_modes, subgroup=None):
        return mapper.map_microscopic_columns_to_supercell(
            Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]]),
            [],
            candidate_modes,
            np.eye(3).tolist(),
            subgroup_context=subgroup or _microscopic_column_provenance_subgroup(),
            frame_id="child-frame",
            atom_ids=("atom-0", "atom-1"),
        )

    damaged_rows = list(modes[0].bush_modes)
    damaged_rows[0] = replace(
        damaged_rows[0],
        displacements=[[9.0, 0.0, 0.0]],
    )
    with pytest.raises(ValueError, match="direct microscopic displacements"):
        _map([replace(modes[0], bush_modes=damaged_rows), modes[1]])

    with pytest.raises(ValueError, match="direction differs"):
        _map([replace(modes[0], opd_dir_raw="(a,a)"), modes[1]])

    with pytest.raises(ValueError, match="mode OPD differs"):
        _map([replace(modes[0], opd_symbol="P2"), modes[1]])

    wrong_bush_opd = list(modes[0].bush_modes)
    wrong_bush_opd[0] = replace(wrong_bush_opd[0], opd_symbol="P2")
    with pytest.raises(ValueError, match="row routing fields differ"):
        _map([replace(modes[0], bush_modes=wrong_bush_opd), modes[1]])

    routing_changes = (
        {"irrep_label": "WRONG"},
        {"wyckoff_site": "b"},
        {"wyckoff_orbit_id": "orbit-b"},
        {"site_irrep": "B"},
        {"opd_component": "z"},
        {"dimension": 9},
        {"k_coords_label": "0,0,0"},
    )
    for changes in routing_changes:
        with pytest.raises(ValueError, match="routing fields differ"):
            _map([replace(modes[0], **changes), modes[1]])

    forged_parent_identity = replace(modes[0].mode_identity, parent_sg=63)
    with pytest.raises(ValueError, match="source provenance disagree"):
        _map([
            replace(modes[0], mode_identity=forged_parent_identity),
            modes[1],
        ])

    forged_k_identity = replace(
        modes[0].mode_identity,
        k_coordinates=("0", "0", "0"),
    )
    with pytest.raises(ValueError, match="source provenance disagree"):
        _map([
            replace(
                modes[0],
                k_coords_label="0,0,0",
                mode_identity=forged_k_identity,
            ),
            modes[1],
        ])

    damaged_row_identity = list(modes[0].bush_modes)
    damaged_row_identity[0] = replace(damaged_row_identity[0], wyckoff_letter="b")
    with pytest.raises(ValueError, match="row routing fields differ"):
        _map([replace(modes[0], bush_modes=damaged_row_identity), modes[1]])

    context_changes = (
        {"space_group_number": 2},
        {"basis_vectors": [[2, 0, 0], [0, 1, 0], [0, 0, 1]]},
        {"origin": [1, 0, 0]},
        {"irrep_label": "GM2+"},
        {"opd_symbol": "P2"},
        {"opd_dir_raw": "(0,a)"},
    )
    for changes in context_changes:
        with pytest.raises(ValueError, match=r"query context|primary direction"):
            _map(modes, replace(_microscopic_column_provenance_subgroup(), **changes))

def test_constant_point_raw_is_authoritative_and_numeric_tamper_fails_closed():
    block = MicroscopicVectorBlock(
        global_irrep="GM1+",
        wyckoff_letter="a",
        site_irrep="A",
        rows=(
            MicroscopicVectorRow(
                point_raw=("0", "0", "0"),
                displacements=((Fraction(1), Fraction(0), Fraction(0)),),
            ),
        ),
        source_order=0,
        direction_symbol="VECTOR,A",
    )
    mode = _mode_from_block(block, 0)
    gamma_provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=_microscopic_column_provenance_QUERY_DIGEST,
        query_order=0,
        query_irrep_label="GM1+",
        source_parent_sg=62,
        source_k_coordinates=("0", "0", "0"),
        source_subgroup_space_group_number=1,
        source_subgroup_basis=np.eye(3, dtype=int).tolist(),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    mode = replace(
        mode,
        k_coords_label="0,0,0",
        mode_identity=replace(
            mode.mode_identity,
            k_coordinates=("0", "0", "0"),
            status="verified",
            reason=None,
        ),
        microscopic_provenance=gamma_provenance,
        bush_modes=[replace(mode.bush_modes[0], point=[0.25, 0.0, 0.0])],
    )
    mapper = _microscopic_column_provenance_DistortionMapper()
    with pytest.raises(ValueError, match="differs from authoritative point_raw"):
        mapper.map_microscopic_columns_to_supercell(
            Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]]),
            [{
                "equivalent_indices": [0],
                "wyckoff_letter": "a",
                "orbit_id": "orbit-a",
            }],
            [mode],
            np.eye(3).tolist(),
            subgroup_context=_microscopic_column_provenance_subgroup(irrep="GM1+"),
            frame_id="child-frame",
            atom_ids=("atom-0",),
        )

def test_mapper_rejects_mixing_missing_frame_and_atom_order_mismatches():
    p0, p1 = _provenance(0), _provenance(1)
    first = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    second = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    canonical = [
        _mapped("m0", p0, first, canonical=True),
        _mapped("m1", p1, second, canonical=True),
    ]
    lattice = np.eye(3)

    mixed = [
        _mapped("m0", p0, first + second, canonical=False),
        _mapped("m1", p1, first - second, canonical=False),
    ]
    with pytest.raises(ValueError, match="non_diagonal_column_mixing"):
        verify_mapped_microscopic_columns(canonical, mixed, lattice)
    with pytest.raises(ValueError, match="missing or unexpected"):
        verify_mapped_microscopic_columns(canonical, mixed[:1], lattice)

    wrong_frame = [
        _mapped("m0", p0, first, canonical=False, frame_id="other-frame"),
        _mapped("m1", p1, second, canonical=False),
    ]
    with pytest.raises(ValueError, match="frame mismatch"):
        verify_mapped_microscopic_columns(canonical, wrong_frame, lattice)

    wrong_order = [
        _mapped("m0", p0, first[::-1], canonical=False, atom_ids=("b", "a")),
        _mapped("m1", p1, second, canonical=False),
    ]
    with pytest.raises(ValueError, match="atom-order mismatch"):
        verify_mapped_microscopic_columns(canonical, wrong_order, lattice)

def test_mapped_column_rejects_incomplete_or_mismatched_mode_identity():
    provenance = _provenance(0)
    vector = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    verified = ModeIdentity(
        parent_sg=62,
        global_irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        wyckoff_letter="a",
        orbit_id="orbit-a",
        site_irrep="A",
        component_index=0,
        component_label="a",
        source="iso_microscopic",
        status="verified",
    )

    with pytest.raises(ValueError, match="not verified"):
        _mapped(
            "a",
            provenance,
            vector,
            canonical=True,
            mode_identity=replace(verified, status="unresolved", reason="missing"),
        )
    with pytest.raises(ValueError, match="source provenance disagree"):
        _mapped(
            "a",
            provenance,
            vector,
            canonical=True,
            mode_identity=replace(verified, global_irrep="WRONG"),
        )

def test_mapped_column_keeps_legacy_positional_fields_and_checks_actual_atom_order():
    original = _mapped(
        "a",
        _provenance(0),
        np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        canonical=True,
    )
    positional = MappedMicroscopicColumn(
        original.mode_id,
        original.amplitude_key,
        original.frame_id,
        original.atom_order_id,
        original.displacements,
        original.provenance,
        original.relation_status,
        original.current_over_canonical_signed_scale,
        original.cartesian_residual_angstrom,
        atom_ids=original.atom_ids,
        mode_identity=original.mode_identity,
    )

    assert positional.mode_identity == original.mode_identity
    assert positional.atom_ids == original.atom_ids
    with pytest.raises(ValueError, match="atom-order ID does not match"):
        replace(original, atom_order_id="forged-order")

def test_verifier_rejects_display_identity_change_even_when_stable_token_matches():
    provenance = _provenance(0)
    vector = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    canonical = _mapped("a", provenance, vector, canonical=True)
    changed_identity = replace(canonical.mode_identity, component_label="renamed")
    assert changed_identity.stable_token == canonical.mode_id
    current = _mapped(
        "renamed",
        provenance,
        vector,
        canonical=False,
        mode_identity=changed_identity,
        amplitude_key=canonical.amplitude_key,
    )

    with pytest.raises(ValueError, match="mode identity mismatch"):
        verify_mapped_microscopic_columns((canonical,), (current,), np.eye(3))

# Original coverage: test_partial_microscopic_mapping.py.

def test_parent_translation_snap_uses_metric_threshold_before_exact_phase():
    lattice = _make_exact_translation_lattice(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        np.diag([4.0, 5.0, 6.0]),
    )

    assert _snap_to_exact_translation(
        np.asarray([1.0 + 2e-5, 0.0, 0.0]), lattice, 1e-4,
    ) == (Fraction(1), Fraction(0), Fraction(0))
    assert _snap_to_exact_translation(
        np.asarray([1.0 + 3e-5, 0.0, 0.0]), lattice, 1e-4,
    ) is None
    assert _snap_to_exact_translation(
        np.asarray([0.25, 0.0, 0.0]), lattice, 1e-4,
    ) is None

def test_parent_translation_snap_is_gl_invariant_and_exact_for_i_centering():
    metric = np.diag([4.0, 4.0, 6.0])
    conventional = _make_exact_translation_lattice(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)), metric,
    )
    sheared = _make_exact_translation_lattice(
        ((1, 0, 0), (10, 1, 0), (0, 0, 1)), metric,
    )
    delta = np.asarray([1.0 + 1e-6, 1.0 - 1e-6, 0.0])
    expected = (Fraction(1), Fraction(1), Fraction(0))
    assert _snap_to_exact_translation(delta, conventional, 1e-4) == expected
    assert _snap_to_exact_translation(delta, sheared, 1e-4) == expected

    body_centered = _make_exact_translation_lattice(
        centering_primitive_matrix("I"), metric,
    )
    centered = _snap_to_exact_translation(
        np.asarray([0.5 + 1e-6, 0.5 - 1e-6, 0.5]),
        body_centered,
        1e-4,
    )
    assert centered == (Fraction(1, 2),) * 3
    assert _has_exact_unit_phase(((Fraction(1), Fraction(1), Fraction(0)),), centered)
    assert not _has_exact_unit_phase(
        ((Fraction(1), Fraction(0), Fraction(0)),), centered,
    )

def _partial_secondary_mapping_fixture(
    *,
    k_coordinates: tuple[str, str, str] = ("0", "0", "0"),
) -> tuple[Structure, list[dict], SubgroupInfo, _partial_microscopic_mapping_DistortionMode, list[list[float]]]:
    """Return one sparse secondary column bound to a complete query context."""

    parent = Structure(
        Lattice(
            [
                [3.0, 0.2, 0.1],
                [0.1, 4.0, 0.3],
                [0.2, 0.4, 5.0],
            ]
        ),
        ["H", "He"],
        [[0.0, 0.0, 0.0], [0.23, 0.31, 0.41]],
    )
    orbit_id = "hydrogen-a"
    sites = [{
        "wyckoff_letter": "a",
        "orbit_id": orbit_id,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    basis = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=1,
        subgroup_index=2,
        size=2,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label="X1-",
        k_point_label="X",
        k_coordinates=["1/2", "0", "0"],
        parent_sg=1,
        basis_vectors=basis,
        basis_raw="(2,0,0),(0,1,0),(0,0,1)",
        origin=[0.0, 0.0, 0.0],
        origin_raw="(0,0,0)",
    )
    block = MicroscopicVectorBlock(
        global_irrep="GM1+",
        direction_symbol="VECTOR,A",
        wyckoff_letter="a",
        site_irrep="A1",
        source_order=0,
        rows=(MicroscopicVectorRow(
            point_raw=("0", "0", "0"),
            displacements=((Fraction(1), Fraction(0), Fraction(0)),),
        ),),
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=_partial_microscopic_mapping_QUERY_DIGEST,
        query_order=0,
        query_irrep_label="GM1+",
        source_parent_sg=1,
        source_k_coordinates=k_coordinates,
        source_subgroup_space_group_number=1,
        source_subgroup_basis=basis,
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="X1-",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    bush = _partial_microscopic_mapping_BushMode(
        irrep_label="GM1+",
        opd_symbol="P1",
        wyckoff_letter="a",
        point=[0.0, 0.0, 0.0],
        point_raw=["0", "0", "0"],
        displacements=[[1.0, 0.0, 0.0]],
    )
    identity = ModeIdentity(
        parent_sg=1,
        global_irrep="GM1+",
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
        irrep_label="GM1+",
        dimension=1,
        mode_type="displacive",
        basis_vectors=[[1.0, 0.0, 0.0]],
        wyckoff_site="a",
        wyckoff_orbit_id=orbit_id,
        k_point_label="GM",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        bush_modes=[bush],
        amplitude_key="GM1+[0,0,0]__a@hydrogen-a__A1(a)",
        site_irrep="A1",
        k_coords_label=",".join(k_coordinates),
        opd_component="a",
        mode_identity=identity,
        microscopic_provenance=provenance,
    )
    return parent, sites, subgroup, mode, basis

def _map_partial_fixture(
    parent: Structure,
    sites: list[dict],
    subgroup: SubgroupInfo,
    mode: _partial_microscopic_mapping_DistortionMode,
    basis: list[list[float]],
) -> np.ndarray:
    return _partial_microscopic_mapping_DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        [mode],
        basis,
        subgroup_translation_lattice=basis,
        subgroup_context=subgroup,
    )[mode.amplitude_key]

def test_exact_partial_secondary_column_can_cover_unit_phase_parent_translations():
    parent, sites, subgroup, mode, basis = _partial_secondary_mapping_fixture()

    mapped = _map_partial_fixture(parent, sites, subgroup, mode, basis)

    child_species = [
        site.species_string for site in build_supercell(parent, basis)
    ]
    hydrogen_rows = mapped[
        np.asarray([species == "H" for species in child_species], dtype=bool)
    ]
    assert hydrogen_rows.shape == (2, 3)
    assert np.allclose(hydrogen_rows, [[0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
    assert mode.mode_identity.status == "partial"

@pytest.mark.parametrize(
    "damage",
    ["missing-provenance", "identity-k", "context-q", "raw-row"],
)
def test_partial_mapping_rejects_missing_or_tampered_source_evidence(damage):
    parent, sites, subgroup, mode, basis = _partial_secondary_mapping_fixture()
    if damage == "missing-provenance":
        mode = replace(mode, microscopic_provenance=None)
    elif damage == "identity-k":
        mode = replace(
            mode,
            mode_identity=replace(
                mode.mode_identity,
                k_coordinates=("1/2", "0", "0"),
            ),
            k_coords_label="1/2,0,0",
        )
    elif damage == "context-q":
        subgroup = replace(
            subgroup,
            origin=[0.5, 0.0, 0.0],
            origin_raw="(1/2,0,0)",
        )
    else:
        mode = replace(
            mode,
            bush_modes=[
                replace(mode.bush_modes[0], point_raw=["1", "0", "0"])
            ],
        )

    with pytest.raises(ValueError, match="partial microscopic mode"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)

def test_partial_mapping_does_not_repeat_a_nonunit_k_phase():
    parent, sites, subgroup, mode, basis = _partial_secondary_mapping_fixture(
        k_coordinates=("1/2", "0", "0"),
    )

    with pytest.raises(ValueError, match="does not cover child atom"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)

def test_x_active_arm_provenance_allows_only_inactive_axis_repetition():
    fixture = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )

    mapped = _map_partial_fixture(*fixture)

    assert np.allclose(
        sorted(mapped[:, 0]),
        [-1 / 3, -1 / 3, -1 / 3, 1 / 3, 1 / 3, 1 / 3],
    )

def test_dt_active_conjugate_pair_allows_transverse_repetition():
    fixture = _active_arm_mapping_fixture(
        irrep="DT1",
        k_coordinates=("0", "1/3", "0"),
        direction="VECTOR,0,0,A,0",
        active_k_coordinates=(("1/3", "0", "0"), ("2/3", "0", "0")),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("1", "0", "0"), ("2", "0", "0")),
        arrows=((0, 0, 0), (1, 0, 0), (-1, 0, 0)),
    )

    mapped = _map_partial_fixture(*fixture)

    assert np.allclose(
        sorted(mapped[:, 0]),
        [-1 / 3, -1 / 3, 0.0, 0.0, 1 / 3, 1 / 3],
    )

def test_sparse_x_column_without_active_arm_provenance_stays_fail_closed():
    parent, sites, subgroup, mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_coordinates=(),
            source_active_k_resolution_kind=None,
        ),
    )

    with pytest.raises(ValueError, match="does not cover child atom"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)

def test_tampered_active_arm_provenance_is_rejected_before_repetition():
    parent, sites, subgroup, mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_coordinates=(("1/2", "0", "0"),),
        ),
    )

    with pytest.raises(ValueError, match="disagrees with its exact VECTOR"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)

def test_tampered_active_arm_resolution_kind_is_rejected_before_repetition():
    parent, sites, subgroup, mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_resolution_kind="compatible_envelope",
        ),
    )

    with pytest.raises(ValueError, match="disagrees with its exact VECTOR"):
        _map_partial_fixture(parent, sites, subgroup, mode, basis)

def test_microscopic_column_mapping_forwards_and_revalidates_active_arm_context():
    """The high-level mapping adapter must not drop or trust active-arm evidence."""

    parent, sites, subgroup, partial_mode, basis = _active_arm_mapping_fixture(
        irrep="X1+",
        k_coordinates=("0", "1/2", "0"),
        direction="VECTOR,A,0",
        active_k_coordinates=(("0", "1/2", "0"),),
        basis=[[3, 0, 0], [0, 2, 0], [0, 0, 1]],
        points=(("0", "0", "0"), ("0", "1", "0")),
        arrows=((1, 0, 0), (-1, 0, 0)),
    )
    mode = replace(
        partial_mode,
        mode_identity=replace(
            partial_mode.mode_identity,
            status="verified",
            reason=None,
        ),
    )
    child = build_supercell(parent, basis)
    kwargs = {
        "subgroup_context": subgroup,
        "frame_id": "active-arm-child-frame",
        "atom_ids": tuple(f"atom-{index}" for index in range(len(child))),
        "subgroup_translation_lattice": basis,
    }
    mapper = _partial_microscopic_mapping_DistortionMapper()

    mapped = mapper.map_microscopic_columns_to_supercell(
        parent,
        sites,
        [mode],
        basis,
        **kwargs,
    )

    assert len(mapped) == 1
    assert mapped[0].mode_id == mode.mode_identity.stable_token
    assert mapped[0].provenance.source_active_k_coordinates == (
        ("0", "1/2", "0"),
    )

    tampered = replace(
        mode,
        microscopic_provenance=replace(
            mode.microscopic_provenance,
            source_active_k_coordinates=(("1/2", "0", "0"),),
        ),
    )
    with pytest.raises(ValueError, match="disagrees with its exact VECTOR"):
        mapper.map_microscopic_columns_to_supercell(
            parent,
            sites,
            [tampered],
            basis,
            **kwargs,
        )

# Original coverage: test_search_methods.py.

def test_method_2_metadata_and_modes():
    engine = Method2Solver(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)

    result = engine.search(
        225,
        subs,
        Method2Query(subgroup_idx=1, distortion_type="displacement"),
        wyckoff_letters=["a", "b"],
    )

    assert result.subgroup.index == 1
    assert result.subgroup.space_group_number == 139
    assert len(result.modes) == 1
    assert result.modes[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]
    # k / IR / OPD 元数据来自所选子群（k 点坐标等不再重复传参）
    assert result.metadata["k_point_label"] == "X"
    assert result.metadata["irrep_label"] == "X1"
    assert result.metadata["number_of_independent_modulations"] == 0

def test_method_2_nmod_accepted():
    engine = Method2Solver(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    result = engine.search(
        225,
        subs,
        Method2Query(subgroup_idx=1, number_of_independent_modulations=1),
        wyckoff_letters=["a"],
    )
    assert result.metadata["number_of_independent_modulations"] == 1
    assert result.modes

def test_method_2_negative_nmod_rejected():
    engine = Method2Solver(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    try:
        engine.search(
            225,
            subs,
            Method2Query(subgroup_idx=1, number_of_independent_modulations=-1),
            wyckoff_letters=["a"],
        )
        raise AssertionError("expected ValueError")
    except ValueError as exc:
        assert "number_of_independent_modulations" in str(exc)

def test_method_2_superposed_ir_request_fails_closed_before_backend_use():
    engine = Method2Solver(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)

    with pytest.raises(NotImplementedError, match="cannot be combined as a union"):
        engine.search(
            225,
            subs,
            Method2Query(subgroup_idx=1, number_of_superposed_irreps=2),
            wyckoff_letters=["a"],
        )

@pytest.mark.parametrize("invalid_count", [0, -1, 1.5, True, "1"])
def test_method_2_superposed_ir_count_must_be_a_positive_integer(invalid_count):
    engine = Method2Solver(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)

    with pytest.raises(ValueError, match="positive integer"):
        engine.search(
            225,
            subs,
            Method2Query(
                subgroup_idx=1,
                number_of_superposed_irreps=invalid_count,
            ),
            wyckoff_letters=["a"],
        )

def test_method_2_missing_subgroup_rejected():
    engine = Method2Solver(_DummyIsoWrapper())
    subs = _DummyIsoWrapper().enumerate_all_special_subgroups(225)
    try:
        engine.search(
            225, subs, Method2Query(subgroup_idx=99),
            wyckoff_letters=["a"],
        )
        raise AssertionError("expected ValueError")
    except ValueError:
        pass

def test_special_modes_complete_only_wholly_rootless_orbits(monkeypatch) -> None:
    subgroup = _method3_route("P2")
    subgroup.k_parameters = []
    bush = _search_methods_DistortionMode(
        irrep_label="P2",
        wyckoff_site="d",
        amplitude_key="P2__d__C1__a",
        bush_modes=[
            _search_methods_BushMode(
                irrep_label="P2",
                opd_symbol="C1",
                wyckoff_letter="d",
                point=[0, 0, 0],
                displacements=[[1, 0, 0]],
            )
        ],
    )
    secondary = _search_methods_DistortionMode(
        irrep_label="GM1+",
        wyckoff_site="e",
        amplitude_key="GM1+[0,0,0]__e__A1(a)",
    )

    class _Iso:
        def calc_distortion_modes(self, parent_sg, selected, wyckoff_letters):
            assert wyckoff_letters == ["d", "e"]
            return [bush]

    def _fixed_space(parent, symmetry_info, selected, letters, smodes, **kwargs):
        assert letters == ["e"]
        return superspace_module.ParametricModeResult(
            modes=[secondary],
            supercell_displacements={secondary.amplitude_key: np.zeros((2, 3))},
            labels={secondary.amplitude_key: "secondary"},
            nmod=0,
            note="fixed-space diagnostic retained",
        )

    monkeypatch.setattr(superspace_module, "compute_parametric_modes", _fixed_space)
    result = superspace_module.compute_special_modes_with_rootless_supplement(
        object(),
        {"space_group_number": 139},
        subgroup,
        ["d", "e"],
        _Iso(),
        object(),
    )

    assert [mode.irrep_label for mode in result.modes] == ["P2", "GM1+"]
    # Without physical-orbit metadata the numerical BUSH column remains
    # available, but its scientific identity must fail closed.
    assert result.modes[0].site_irrep == ""
    assert result.modes[0].mode_identity.status == "unresolved"
    assert (
        result.modes[0].mode_identity.reason
        == "physical_orbit_metadata_unavailable"
    )
    assert result.modes[1] is secondary
    assert set(result.supercell_displacements) == {secondary.amplitude_key}
    assert "rootless orbit(s) e" in result.note
    assert "fixed-space diagnostic retained" in result.note
    assert "ISO microscopic query failed: AttributeError" in result.note

def test_special_modes_split_repeated_letters_into_physical_orbits(
    monkeypatch,
) -> None:
    subgroup = _method3_route("P2")
    subgroup.k_parameters = []
    bush = _search_methods_DistortionMode(
        irrep_label="P2",
        wyckoff_site="e",
        amplitude_key="P2__e__C1__a",
        bush_modes=[_search_methods_BushMode(
            irrep_label="P2",
            opd_symbol="C1",
            wyckoff_letter="e",
            point=[0, 0, 0],
            displacements=[[1, 0, 0]],
        )],
    )
    secondary = _search_methods_DistortionMode(
        irrep_label="GM1+",
        wyckoff_site="g",
        wyckoff_orbit_id="orbit-g",
        amplitude_key="GM1+[0,0,0]__g__A1(a)",
    )
    sites = [
        {"orbit_id": "orbit-e-1", "wyckoff_letter": "e"},
        {"orbit_id": "orbit-e-2", "wyckoff_letter": "e"},
        {"orbit_id": "orbit-g", "wyckoff_letter": "g"},
    ]

    class _Iso:
        def calc_distortion_modes(self, parent_sg, selected, wyckoff_letters):
            assert wyckoff_letters == ["e", "g"]
            return [bush]

    def _fixed_space(parent, symmetry_info, selected, letters, smodes, **kwargs):
        assert letters == ["g"]
        assert kwargs["wyckoff_orbit_ids"] == ["orbit-g"]
        return superspace_module.ParametricModeResult(
            modes=[secondary],
            supercell_displacements={secondary.amplitude_key: np.zeros((2, 3))},
            labels={secondary.amplitude_key: "secondary"},
            nmod=0,
        )

    monkeypatch.setattr(superspace_module, "compute_parametric_modes", _fixed_space)
    result = superspace_module.compute_special_modes_with_rootless_supplement(
        object(),
        {"space_group_number": 139, "wyckoff_sites": sites},
        subgroup,
        ["e", "g"],
        _Iso(),
        object(),
    )

    rooted = [mode for mode in result.modes if mode.wyckoff_site == "e"]
    assert [mode.wyckoff_orbit_id for mode in rooted] == [
        "orbit-e-1", "orbit-e-2",
    ]
    assert len({mode.amplitude_key for mode in rooted}) == 2
    assert result.modes[-1] is secondary

def test_special_modes_complete_all_rootless_physical_orbits(monkeypatch) -> None:
    subgroup = _method3_route("P2")
    subgroup.k_parameters = []
    sites = [
        {"orbit_id": "orbit-e-1", "wyckoff_letter": "e"},
        {"orbit_id": "orbit-e-2", "wyckoff_letter": "e"},
    ]
    secondary = [
        _search_methods_DistortionMode(
            irrep_label="GM1+",
            wyckoff_site="e",
            wyckoff_orbit_id=orbit_id,
            amplitude_key=f"GM1+__{orbit_id}",
        )
        for orbit_id in ("orbit-e-1", "orbit-e-2")
    ]

    class _Iso:
        def calc_distortion_modes(self, parent_sg, selected, wyckoff_letters):
            return []

        def calc_microscopic_distortion_modes(
            self, parent_sg, selected, wyckoff_letters
        ):
            return secondary

    def _fixed_space(parent, symmetry_info, selected, letters, smodes, **kwargs):
        assert letters == ["e"]
        assert kwargs["wyckoff_orbit_ids"] == ["orbit-e-1", "orbit-e-2"]
        return superspace_module.ParametricModeResult(
            modes=secondary,
            supercell_displacements={
                mode.amplitude_key: np.zeros((2, 3)) for mode in secondary
            },
            labels={mode.amplitude_key: mode.amplitude_key for mode in secondary},
            nmod=0,
            note="rootless fixed-space note",
        )

    monkeypatch.setattr(superspace_module, "compute_parametric_modes", _fixed_space)
    def _reject_replacement(*_args, **_kwargs):
        raise ValueError("whole-space mismatch")

    monkeypatch.setattr(
        superspace_module,
        "_canonicalize_rootless_supplement",
        _reject_replacement,
    )
    result = superspace_module.compute_special_modes_with_rootless_supplement(
        object(),
        {"space_group_number": 139, "wyckoff_sites": sites},
        subgroup,
        ["e"],
        _Iso(),
        object(),
    )

    assert result.modes == secondary
    assert set(result.supercell_displacements) == {
        "GM1+__orbit-e-1", "GM1+__orbit-e-2",
    }
    assert "rootless fixed-space note" in result.note
    assert "canonical ISO microscopic replacement rejected: whole-space mismatch" in result.note

def test_exact_fixed_space_method3_row_uses_complete_folded_mode_space(
    monkeypatch,
) -> None:
    target = _method3_route("GM2+")
    target.k_point_label = ""
    target.irrep_label = ""
    target.k_parameters = []
    target.method3_route_resolution = "exact_fixed_space"
    mode = _search_methods_DistortionMode(
        irrep_label="GM2+",
        wyckoff_site="a",
        amplitude_key="GM2+[0,0,0]__a__A(a)",
    )
    seen = {}

    def _fixed_space(parent, symmetry_info, selected, letters, smodes, **kwargs):
        seen.update(kwargs)
        assert selected is target
        assert letters == ["a"]
        return superspace_module.ParametricModeResult(
            modes=[mode],
            supercell_displacements={mode.amplitude_key: np.zeros((1, 3))},
            labels={mode.amplitude_key: "complete coupled mode"},
            nmod=0,
            note="complete fixed space",
        )

    monkeypatch.setattr(superspace_module, "compute_parametric_modes", _fixed_space)
    result = Method2Solver(_DummyIsoWrapper()).search(
        3,
        [target],
        Method2Query(subgroup_idx=target.index, number_of_independent_modulations=3),
        wyckoff_letters=["a"],
        structure=Structure(Lattice.cubic(4), ["Na"], [[0, 0, 0]]),
        wyckoff_sites=[],
        smodes=object(),
        kpoints=[],
        symmetry_info={"space_group_number": 3},
    )

    assert result.modes == [mode]
    assert seen["nmod"] == 0
    assert result.metadata["coupled_mode_note"] == "complete fixed space"
    assert set(result.metadata["supercell_displacements"]) == {mode.amplitude_key}

# Original coverage: test_site_irreps.py.

WORKSPACE = ISODISTORT_ROOT.parent

def _mode(
    points: list[tuple[str, str, str]],
    vectors: list[tuple[float, float, float]],
    *,
    irrep: str = "N1+",
    letter: str = "e",
    site_irrep: str = "A1",
    k_coordinates: tuple[str, str, str] = ("1/2", "0", "1/2"),
    source_block_order: int = 0,
) -> _site_irreps_DistortionMode:
    bushes = [
        _site_irreps_BushMode(
            irrep_label=irrep,
            opd_symbol="4D1",
            wyckoff_letter=letter,
            point=[0.0, 0.0, 0.0],
            point_raw=list(point),
            displacements=[list(vector)],
        )
        for point, vector in zip(points, vectors, strict=True)
    ]
    exact_block = MicroscopicVectorBlock(
        irrep,
        letter,
        site_irrep,
        tuple(
            MicroscopicVectorRow(
                point_raw=point,
                displacements=(tuple(Fraction(str(value)) for value in vector),),
            )
            for point, vector in zip(points, vectors, strict=True)
        ),
        source_order=source_block_order,
        direction_symbol="VECTOR,A",
    )
    try:
        exact_source_k = tuple(str(Fraction(value)) for value in k_coordinates)
    except (ValueError, ZeroDivisionError):
        exact_source_k = ()
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=hashlib.sha256(f"synthetic:{irrep}".encode()).hexdigest(),
        query_order=0,
        query_irrep_label=irrep,
        source_parent_sg=139,
        source_k_coordinates=exact_source_k,
        source_subgroup_space_group_number=2,
        source_subgroup_basis=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="N1+",
        source_subgroup_opd_symbol="4D1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=exact_block,
        column_index=0,
    )
    return _site_irreps_DistortionMode(
        irrep_label=irrep,
        dimension=1,
        wyckoff_site=letter,
        bush_modes=bushes,
        amplitude_key=f"{irrep}__{letter}__test",
        site_irrep=site_irrep,
        k_coords_label=",".join(k_coordinates),
        mode_identity=ModeIdentity(
            parent_sg=139,
            global_irrep=irrep,
            k_coordinates=k_coordinates,
            wyckoff_letter=letter,
            site_irrep=site_irrep,
            component_index=0,
            component_label="a",
            source="iso_microscopic",
            status="verified",
        ),
        microscopic_provenance=provenance,
    )

def _rootless_basis_fixture():
    orbit_id = "physical-orbit-a"
    x = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    y = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    reference_modes = [
        _site_irreps_DistortionMode(
            irrep_label="GM1+",
            wyckoff_site="e",
            wyckoff_orbit_id=orbit_id,
            amplitude_key=key,
            mode_identity=ModeIdentity.unresolved(
                parent_sg=139,
                global_irrep="GM1+",
                k_coordinates=("0", "0", "0"),
                wyckoff_letter="e",
                orbit_id=orbit_id,
                reason="smodes_basis_has_no_iso_microscopic_identity",
            ),
        )
        for key in ("smodes-1", "smodes-2")
    ]
    supplement = ParametricModeResult(
        modes=reference_modes,
        supercell_displacements={"smodes-1": x + y, "smodes-2": x - y},
        labels={"smodes-1": "unresolved 1", "smodes-2": "unresolved 2"},
        nmod=0,
    )
    microscopic = []
    for key, site_irrep, source_order in (
        ("canonical-x", "A1", 0),
        ("canonical-y", "B1", 1),
    ):
        mode = _mode(
            [("x", "0", "z")],
            [(1.0, 0.0, 0.0)],
            irrep="GM1+",
            site_irrep=site_irrep,
            k_coordinates=("0", "0", "0"),
            source_block_order=source_order,
        )
        microscopic.append(replace(
            mode,
            amplitude_key=key,
            wyckoff_orbit_id=orbit_id,
            mode_identity=replace(mode.mode_identity, orbit_id=orbit_id),
        ))
    return supplement, microscopic, {"canonical-x": x, "canonical-y": y}

def test_rootless_whole_space_replacement_accepts_nondiagonal_basis_change():
    supplement, microscopic, arrays = _rootless_basis_fixture()
    lattice = np.asarray([[2.0, 0.4, 0.0], [0.0, 3.0, 0.3], [0.0, 0.0, 4.0]])

    replaced = _replace_rootless_supplement_basis(
        supplement,
        microscopic,
        arrays,
        lattice,
    )

    assert [mode.amplitude_key for mode in replaced.modes] == [
        "canonical-x",
        "canonical-y",
    ]
    assert set(replaced.supercell_displacements) == {"canonical-x", "canonical-y"}
    assert np.array_equal(replaced.supercell_displacements["canonical-x"], arrays["canonical-x"])
    assert all(mode.mode_identity.status == "verified" for mode in replaced.modes)

def test_rootless_whole_space_replacement_rejects_cross_orbit_and_missing_column():
    supplement, microscopic, arrays = _rootless_basis_fixture()
    crossed = [
        microscopic[0],
        replace(
            microscopic[1],
            wyckoff_orbit_id="physical-orbit-b",
            mode_identity=replace(
                microscopic[1].mode_identity,
                orbit_id="physical-orbit-b",
            ),
        ),
    ]
    with pytest.raises(ValueError, match="cross or omit physical orbits"):
        _replace_rootless_supplement_basis(
            supplement,
            crossed,
            arrays,
            np.eye(3),
        )
    with pytest.raises(ValueError, match="complete child-fixed space"):
        _replace_rootless_supplement_basis(
            supplement,
            microscopic[:1],
            {"canonical-x": arrays["canonical-x"]},
            np.eye(3),
        )

def test_rootless_microscopic_queries_split_same_digest_by_query_order():
    _supplement, microscopic, _arrays = _rootless_basis_fixture()
    first = microscopic[0]
    second = microscopic[1]
    shared_digest = first.microscopic_provenance.query_digest
    second_provenance = replace(
        second.microscopic_provenance,
        query_digest=shared_digest,
        query_order=1,
        query_irrep_label="GM2+",
        source_global_irrep="GM2+",
    )
    second_identity = replace(second.mode_identity, global_irrep="GM2+")
    second = replace(
        second,
        irrep_label="GM2+",
        microscopic_provenance=second_provenance,
        mode_identity=second_identity,
    )

    ordered = _ordered_complete_microscopic_space([second, first])

    assert [mode.irrep_label for mode in ordered] == ["GM1+", "GM2+"]

    duplicate_order = replace(
        second,
        microscopic_provenance=replace(second_provenance, query_order=0),
    )
    with pytest.raises(ValueError, match="multiple irreps"):
        _ordered_complete_microscopic_space([first, duplicate_order])

    missing_source_column = replace(
        first,
        microscopic_provenance=replace(
            first.microscopic_provenance,
            source_block_column_count=2,
        ),
    )
    with pytest.raises(ValueError, match="source_block_column_missing"):
        _ordered_complete_microscopic_space([missing_source_column])

def test_symbolic_point_alignment_handles_permutation_and_gamma_integer_shifts():
    reference = [
        _mode(
            [("x", "0", "-z"), ("y", "0", "z")],
            [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        )
    ]
    permuted = [
        _mode(
            [("y", "0", "z"), ("x", "0", "-z")],
            [(0.0, 1.0, 0.0), (1.0, 0.0, 0.0)],
        )
    ]
    assert validate_symbolic_mode_subspaces(
        reference, permuted, np.eye(3)
    ).matched

    shifted_reference = [
        _mode(
            [("x", "0", "-z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    equivalent = [
        _mode(
            [("x+1", "0", "1-z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    assert validate_symbolic_mode_subspaces(
        shifted_reference, equivalent, np.eye(3)
    ).matched

def test_non_gamma_exact_common_rows_prove_nondiagonal_whole_space_match():
    reference_points = [
        ("x", "0", "z"),
        ("y", "0", "z"),
        ("y+1", "0", "z"),
    ]
    candidate_points = [
        ("x", "0", "z"),
        ("y", "0", "z"),
        ("y+2", "0", "z"),
    ]
    reference = [
        _mode(
            reference_points,
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
        ),
        _mode(
            reference_points,
            [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            site_irrep="B1",
        ),
    ]
    candidate = [
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        ),
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)],
            site_irrep="B1",
        ),
    ]

    validation = validate_symbolic_mode_subspaces(
        reference, candidate, np.eye(3)
    )

    assert validation.matched
    assert validation.reference_rank == validation.candidate_rank == 2

def test_non_gamma_rank_losing_common_rows_still_require_phase_proof():
    reference_points = [("x", "0", "z"), ("y", "0", "z")]
    candidate_points = [("x", "0", "z"), ("y+1", "0", "z")]
    reference = [
        _mode(reference_points, [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)]),
        _mode(
            reference_points,
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            site_irrep="B1",
        ),
    ]
    candidate = [
        _mode(candidate_points, [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)]),
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            site_irrep="B1",
        ),
    ]

    with pytest.raises(ValueError, match="star-arm phase evidence"):
        validate_symbolic_mode_subspaces(reference, candidate, np.eye(3))

def test_x_star_common_minus_phase_extends_sparse_microscopic_column():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z")],
            [(1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x+1", "0", "z")],
            [(-1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]

    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139,
    )

    assert diagnostics[("X1+", "e")].matched
    assert selected[0].mode_identity.status == "verified"
    assert selected[0].microscopic_domain_extension is not None
    assert (
        selected[0].microscopic_domain_extension.comparison_domain_kind
        == "phase_proven_transported_points"
    )
    assert selected[0].bush_modes[0].point_raw == ["x", "0", "z"]
    assert selected[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]

def test_star_phase_matching_assigns_distinct_equivalent_bush_images():
    reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x-1", "0", "z"), ("x", "0", "z")],
            [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            irrep="X1+",
            k_coordinates=("1/2", "1/2", "0"),
        )
    ]

    validation = validate_symbolic_mode_subspaces(
        reference, candidate, np.eye(3)
    )
    assert validation.matched
    assert validation.reference_rank == validation.candidate_rank == 1
    reference_points, reference_rows = _common_point_rows(reference)
    candidate_points, _candidate_rows = _common_point_rows(candidate)
    transports = _match_candidate_points_to_reference(
        reference_points,
        reference_rows,
        candidate_points,
        parent_sg=139,
        k_coordinates=(Fraction(1, 2), Fraction(1, 2), Fraction(0)),
    )
    assert [point for point, _translation in transports] == [
        reference_points[1],
        reference_points[0],
    ]
    assert [translation for _point, translation in transports] == [
        (-2, 0, 0),
        (0, 0, 0),
    ]

def test_symbolic_point_matching_fails_when_no_complete_injection_exists():
    reference_points = [("x", "0", "z"), ("x+1", "0", "z")]
    candidate_points = [
        ("x-1", "0", "z"),
        ("x+2", "0", "z"),
        ("x+3", "0", "z"),
    ]
    reference = [
        _mode(
            reference_points,
            [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        ),
        _mode(
            reference_points,
            [(0.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            site_irrep="B1",
            source_block_order=1,
            k_coordinates=("0", "0", "0"),
        ),
    ]
    candidate = [
        _mode(
            candidate_points,
            [(1.0, 0.0, 0.0)] * 3,
            k_coordinates=("0", "0", "0"),
        ),
        _mode(
            candidate_points,
            [(0.0, 1.0, 0.0)] * 3,
            site_irrep="B1",
            source_block_order=1,
            k_coordinates=("0", "0", "0"),
        ),
    ]

    with pytest.raises(ValueError, match="non-injectively"):
        validate_symbolic_mode_subspaces(reference, candidate, np.eye(3))

@pytest.mark.parametrize(
    ("irrep", "k_coordinates"),
    [
        ("N1+", ("1/2", "0", "1/2")),
        ("DT1", ("1/4", "0", "0")),
    ],
)
def test_non_gamma_integer_translation_without_star_arm_identity_fails_closed(
    irrep, k_coordinates,
):
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z"), ("y", "0", "z")],
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            irrep=irrep,
            k_coordinates=k_coordinates,
        )
    ]
    candidate = [
        _mode(
            [("x+1", "0", "z"), ("y", "0", "z")],
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
            irrep=irrep,
            k_coordinates=k_coordinates,
        )
    ]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert diagnostics[(irrep, "e")].reason == "subspace_validation_error"
    assert selected[0].mode_identity.status == "unresolved"

def test_r_hexagonal_reciprocal_centering_keeps_exact_star_arms_distinct():
    zero = Fraction(0)
    half = Fraction(1, 2)

    assert _g_allowed(1, 1, 0, "R")
    assert not _g_allowed(1, 0, 0, "R")
    assert _k_equivalent_exact((zero, zero, zero), (Fraction(1), Fraction(1), zero), "R")
    assert not _k_equivalent_exact(
        (zero, zero, zero), (Fraction(1), zero, zero), "R"
    )
    assert k_star_fraction_vectors((zero, half, zero), 146) == [
        (zero, half, zero),
        (half, zero, Fraction(1)),
        (Fraction(3, 2), half, zero),
    ]

def test_r_hexagonal_star_phase_proof_rejects_mixed_real_arm_phases():
    with pytest.raises(ValueError, match="star-arm phase evidence"):
        _proven_real_transport_phase(
            (-1, 0, 0),
            parent_sg=146,
            k_coordinates=(Fraction(0), Fraction(1, 2), Fraction(0)),
        )

def test_unknown_reciprocal_centering_fails_closed():
    with pytest.raises(ValueError, match="unsupported reciprocal-lattice centering"):
        _g_allowed(0, 0, 0, "Q")

def test_unparseable_symbolic_transport_k_fails_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("a", "0", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x+1", "0", "z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("a", "0", "0"),
        )
    ]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert diagnostics[("N1+", "e")].reason == "microscopic_k_identity_not_exact"
    assert selected[0].mode_identity.status == "unresolved"

def test_gamma_multiple_integer_transporters_require_column_consistency():
    consistent_reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    candidate = [
        _mode(
            [("x+2", "0", "z")],
            [(1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    assert validate_symbolic_mode_subspaces(
        consistent_reference, candidate, np.eye(3)
    ).matched

    conflicting_reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
            k_coordinates=("0", "0", "0"),
        )
    ]
    with pytest.raises(ValueError, match="unresolved phase"):
        validate_symbolic_mode_subspaces(
            conflicting_reference, candidate, np.eye(3)
        )

def test_unproven_symbolic_point_equivalence_fails_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [_mode([("x", "0", "z")], [(1.0, 0.0, 0.0)])]
    candidate = [_mode([("y", "0", "z")], [(1.0, 0.0, 0.0)])]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert diagnostics[("N1+", "e")].reason == "subspace_validation_error"
    assert len(selected) == len(reference)
    assert all(mode.site_irrep == "" for mode in selected)
    assert all(mode.mode_identity.status == "unresolved" for mode in selected)

def test_one_proven_star_phase_edge_is_sufficient_for_sparse_extension():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    reference = [
        _mode(
            [("x", "0", "z"), ("x+1", "0", "z")],
            [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)],
        )
    ]
    candidate = [_mode([("x+2", "0", "z")], [(1.0, 0.0, 0.0)])]

    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )

    assert diagnostics[("N1+", "e")].matched
    assert selected[0].mode_identity.status == "verified"
    assert selected[0].microscopic_domain_extension is not None
    assert [bush.displacements for bush in selected[0].bush_modes] == [
        [[1.0, 0.0, 0.0]],
        [[-1.0, 0.0, 0.0]],
    ]

def test_subspace_rank_mismatch_and_missing_copy_fail_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    point = [("x", "0", "z")]
    reference = [
        _mode(point, [(1.0, 0.0, 0.0)]),
        _mode(point, [(0.0, 1.0, 0.0)], site_irrep="B1"),
    ]
    candidate = [_mode(point, [(1.0, 0.0, 0.0)])]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, candidate, parent_sg=139
    )
    assert not diagnostics[("N1+", "e")].matched
    assert diagnostics[("N1+", "e")].reason == "rank_mismatch"
    assert len(selected) == 2
    assert all(mode.site_irrep == "" for mode in selected)
    assert all(mode.mode_identity.reason == "rank_mismatch" for mode in selected)

def test_duplicate_microscopic_column_fails_closed():
    parent = Structure(Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])
    point = [("x", "0", "z")]
    reference = [_mode(point, [(1.0, 0.0, 0.0)])]
    duplicate = [
        _mode(point, [(1.0, 0.0, 0.0)]),
        _mode(point, [(1.0, 0.0, 0.0)], source_block_order=1),
    ]
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, reference, duplicate, parent_sg=139
    )
    assert not diagnostics[("N1+", "e")].matched
    assert (
        diagnostics[("N1+", "e")].reason
        == "candidate_columns_not_independent"
    )
    assert len(selected) == 1
    assert selected[0].mode_identity.status == "unresolved"
    assert selected[0].mode_identity.reason == "candidate_columns_not_independent"

def test_five_physical_e_orbits_receive_distinct_stable_identities():
    mode = _mode([("0", "0", "z")], [(0.0, 0.0, 1.0)])
    sites = [
        {
            "wyckoff_letter": "e",
            "species": species,
            "orbit_id": f"physical-e-{index}",
            "display_label": label,
        }
        for index, (species, label) in enumerate(
            [("La", "La1"), ("La", "La2"), ("Ni", "Ni2"), ("O", "O4"), ("O", "O2")]
        )
    ]
    clones = _instantiate_bush_modes_by_orbit([mode], sites, ["e"])
    assert len(clones) == 5
    assert {item.wyckoff_orbit_id for item in clones} == {
        f"physical-e-{index}" for index in range(5)
    }
    assert {item.mode_identity.orbit_id for item in clones} == {
        f"physical-e-{index}" for index in range(5)
    }
    assert len({item.mode_identity.stable_token for item in clones}) == 5
    assert len({item.amplitude_key for item in clones}) == 5

def test_missing_or_duplicate_physical_orbit_metadata_fails_closed():
    mode = _mode([("0", "0", "z")], [(0.0, 0.0, 1.0)])
    duplicate = [
        {"wyckoff_letter": "e", "orbit_id": "same"},
        {"wyckoff_letter": "e", "orbit_id": "same"},
    ]
    with pytest.raises(ValueError, match="not unique"):
        _instantiate_bush_modes_by_orbit([mode], duplicate, ["e"])

    sites = [{"wyckoff_letter": "e", "orbit_id": "present"}]
    with pytest.raises(ValueError, match="metadata is missing"):
        _instantiate_bush_modes_by_orbit(
            [mode], sites, ["e"], requested_orbit_ids=["absent"]
        )
    with pytest.raises(ValueError, match="contain duplicates"):
        _instantiate_bush_modes_by_orbit(
            [mode], sites, ["e"], requested_orbit_ids=["present", "present"]
        )

    unresolved = _instantiate_bush_modes_by_orbit([mode], [], ["e"])
    assert len(unresolved) == 1
    assert unresolved[0].site_irrep == ""
    assert unresolved[0].mode_identity.status == "unresolved"
    assert unresolved[0].mode_identity.reason == "physical_orbit_metadata_unavailable"

def _wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(
            ["wsl.exe", "--status"],  # noqa: S607 - platform executable
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return result.returncode == 0

@pytest.mark.skipif(not _wsl_available(), reason="WSL/ISO binary unavailable")
def test_live_4310_n1plus_4d1_all_192_site_irrep_labels_match_official():
    parent_path = WORKSPACE / "experiment_data" / "4310_tetra.cif"
    official_path = (
        WORKSPACE / "output_compare" / "4310_tetra.cif" / "官网" / "Method1"
        / "N1+_4D1_SG2" / "subgroup.cif"
    )
    if not parent_path.is_file() or not official_path.is_file():
        pytest.skip("validated 4310 parent/official output unavailable")

    subgroup = _site_irreps_SubgroupInfo(
        index=0,
        space_group_number=2,
        space_group_symbol="P-1",
        subgroup_index=64,
        size=8,
        opd_symbol="4D1",
        opd_vector=[1.0, 1.0, 1.0, 1.0],
        basis_vectors=[[2, 0, 0], [0, 2, 0], [-1, -1, 1]],
        origin=[0.0, 0.0, 0.0],
        k_point_label="N",
        irrep_label="N1+",
        k_coordinates=["1/2", "0", "1/2"],
        parent_sg=139,
        opd_dir_raw="(a;b;c;d)",
        basis_raw="(2,0,0),(0,2,0),(-1,-1,1)",
        origin_raw="(0,0,0)",
    )
    parent = read_cif(parent_path)
    symmetry_info = SymmetryValidator().validate(parent)
    sites = symmetry_info["wyckoff_sites"]
    parent_wyckoff_display(parent, sites, parent_path)
    assert sum(site["wyckoff_letter"] == "e" for site in sites) == 5

    wrapper = _site_irreps_IsoWrapper()
    bush = wrapper.calc_distortion_modes(139, subgroup, ["e", "g", "c"])
    microscopic = wrapper.calc_microscopic_distortion_modes(
        139, subgroup, ["e", "g", "c"]
    )
    selected, diagnostics = _select_verified_microscopic_modes(
        parent, bush, microscopic, parent_sg=139
    )
    assert diagnostics
    assert all(item.matched for item in diagnostics.values())
    assert len(selected) == 96
    assert all(mode.mode_identity.status == "verified" for mode in selected)
    modes = _instantiate_bush_modes_by_orbit(selected, sites, ["e", "g", "c"])
    assert len(modes) == 192
    assert len({mode.wyckoff_orbit_id for mode in modes if mode.wyckoff_site == "e"}) == 5

    official_pattern = re.compile(
        r"^\s*\d+\s+\S+\[[^]]+\](?P<irrep>\S+?)\([^)]*\)"
        r"\[(?P<site>[^:]+):(?P<letter>[a-z]):dsp\]"
        r"(?P<site_irrep>[^\s(]+)\((?P<component>[^)]+)\)"
    )
    expected: Counter[tuple[str, str, str, str, str]] = Counter()
    for line in official_path.read_text(encoding="utf-8-sig").splitlines():
        match = official_pattern.match(line)
        if match is not None:
            expected[(
                match.group("irrep"),
                match.group("site"),
                match.group("letter"),
                match.group("site_irrep"),
                match.group("component"),
            )] += 1
    assert sum(expected.values()) == 192

    actual: Counter[tuple[str, str, str, str, str]] = Counter()
    for mode in modes:
        site_label, _species = _display_site_label(
            sites, mode.wyckoff_site, mode.wyckoff_orbit_id
        )
        actual[(
            mode.irrep_label,
            site_label,
            mode.wyckoff_site,
            mode.site_irrep,
            mode.opd_component,
        )] += 1
    assert actual == expected

# Original coverage: test_superspace.py.

_SAMPLE = (
    TESTS_ROOT / "fixtures" / "smodes_sample.out"
)

def _kp(label, coords, *, special=True, parameters=None):
    return SimpleNamespace(
        label=label,
        coordinates=list(coords),
        is_special=special,
        parameters=list(parameters or []),
    )

def test_parse_smodes_sample_blocks():
    text = _SAMPLE.read_text(encoding="utf-8", errors="replace")
    blocks = parse_smodes_mode_blocks(text)
    assert blocks
    irreps = {b.irrep for b in blocks}
    assert "GM3-" in irreps
    assert "GM4-" in irreps
    gm3 = next(b for b in blocks if b.irrep == "GM3-")
    assert gm3.k_label == "GM"
    dense = gm3.dense_displacements()
    assert dense.shape[1] == 3
    assert float(np.max(np.abs(dense))) > 0
    z_blocks = [b for b in blocks if b.k_label == "Z"]
    assert z_blocks
    assert z_blocks[0].superlattice_cart.shape == (3, 3)

def test_is_harmonic_of_line_and_gamma():
    q = np.array([0.0, 0.0, 1.0 / 6.0])
    assert is_harmonic_of(np.array([0.0, 0.0, 1.0 / 6.0]), [q])
    assert is_harmonic_of(np.array([0.0, 0.0, 1.0 / 3.0]), [q])
    assert is_harmonic_of(np.array([0.0, 0.0, 0.5]), [q])
    assert is_harmonic_of(np.zeros(3), [q])  # 3D average structure
    # X of I4/mmm is not a harmonic of LD along c*
    assert not is_harmonic_of(np.array([0.5, 0.5, 0.0]), [q])
    assert _harmonic_order(np.array([0.0, 0.0, 5.0 / 6.0]), q) == 5
    assert _harmonic_order(np.array([0.0, 0.0, 1.0 / 6.0]), q) == 1

def test_harmonic_order_uses_the_exact_finite_period_not_a_fixed_cap():
    q = np.array([0.0, 0.0, 1.0 / 30.0])
    thirteenth = np.array([0.0, 0.0, 13.0 / 30.0])

    assert is_harmonic_of(thirteenth, [q])
    assert _harmonic_order(thirteenth, q) == 13
    # Reciprocal equivalence uses the primitive translations of the centered
    # parent.  In an I lattice (5/4,0,0) is the fifth, not the first, harmonic
    # of (1/4,0,0), because (1,0,0) is not an allowed reciprocal vector.
    assert _harmonic_order(
        np.array([5.0 / 4.0, 0.0, 0.0]),
        np.array([1.0 / 4.0, 0.0, 0.0]),
        parent_sg=229,
    ) == 5

def test_format_k_token_preserves_long_exact_commensurate_periods():
    assert format_k_token(1.0 / 49.0) == "1/49"
    assert format_k_token(Fraction(1, 2_000_001)) == "1/2000001"
    assert format_k_token(1e-13) == "1/10000000000000"
    assert format_k_token(-1e-13) == "-1/10000000000000"

def test_opd_free_components_c1():
    assert opd_free_components("(a,b)", "C1") == ["a", "b"]
    assert opd_free_components("(a)", "P1") == ["a"]
    assert opd_free_components("", "C1") == ["a", "b"]

def test_parametric_rejection_reason_is_bound_to_every_unresolved_identity():
    mode = _superspace_DistortionMode(
        irrep_label="Y1",
        amplitude_key="Y1__f__unresolved(a)",
        wyckoff_site="f",
        wyckoff_orbit_id="sg123:2f:O:test",
        mode_identity=ModeIdentity.unresolved(
            parent_sg=123,
            global_irrep="Y1",
            k_coordinates=("1/3", "1/2", "0"),
            wyckoff_letter="f",
            orbit_id="sg123:2f:O:test",
            reason="smodes_basis_has_no_iso_microscopic_identity",
        ),
    )
    result = ParametricModeResult(
        modes=[mode],
        supercell_displacements={mode.amplitude_key: np.ones((1, 3))},
        labels={},
        nmod=0,
    )
    reason = (
        "parameter_k_microscopic_canonicalization_rejected; "
        "criterion=ValueError: candidate_columns_not_independent; "
        "parent_SG=123; k=(1/3,1/2,0); KVALUE=(2/3); "
        "IR=Y1; direction_selector=(a,0;0,0)"
    )

    updated = _with_unresolved_identity_reason(result, reason)

    assert updated.note == reason
    assert updated.modes[0].mode_identity.reason == reason
    assert updated.modes[0].mode_identity.status == "unresolved"

def test_cached_displacement_actions_equal_reynolds_projection():
    structure = Structure(
        Lattice.cubic(4.0),
        ["Fe", "Fe"],
        [[0.25, 0, 0], [0.75, 0, 0]],
    )
    identity = SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))
    inversion = SymmOp.from_rotation_and_translation(-np.eye(3), np.zeros(3))
    ops = [identity, inversion]
    disp = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

    actions = _displacement_actions(structure, ops)
    assert actions is not None
    projected = _project_invariant(disp, structure, ops, actions=actions)
    expected_inverted = np.asarray([[-4.0, -5.0, -6.0], [-1.0, -2.0, -3.0]])
    assert np.allclose(projected, 0.5 * (disp + expected_inverted))

def test_displacement_action_site_tolerance_uses_cartesian_angstrom_metric():
    """A fractional residual is judged in the actual direct-cell metric."""
    structure = Structure(
        Lattice.orthorhombic(10.0, 2.0, 3.0), ["Fe"], [[0.0, 0.0, 0.0]],
    )
    within_one_millangstrom = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.00005, 0.0, 0.0],
    )
    outside_one_millangstrom = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.00015, 0.0, 0.0],
    )

    assert _displacement_actions(
        structure, [within_one_millangstrom], atol=0.001,
    ) is not None
    assert _displacement_actions(
        structure, [outside_one_millangstrom], atol=0.001,
    ) is None

def test_displacement_action_uses_true_minimum_image_in_skew_cell():
    """Fractional component wrapping is not a nearest-image algorithm."""
    structure = Structure(
        Lattice.from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 10.0),
        ["Fe"],
        [[0.0, 0.0, 0.0]],
    )
    operation = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.49, 0.49, 0.0],
    )

    # The component-wrapped vector is 9.76 Å, while the true periodic image
    # is 0.894 Å because the a and b axes are nearly parallel.
    assert _displacement_actions(structure, [operation], atol=1.0) is not None

def test_displacement_action_solves_global_species_bijection():
    """A perfect matching must not fail because of greedy target ownership."""
    structure = Structure(
        Lattice.cubic(10.0),
        ["Fe", "Fe"],
        [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]],
    )
    operation = SymmOp.from_rotation_and_translation(
        np.eye(3), [0.11, 0.0, 0.0],
    )

    actions = _displacement_actions(structure, [operation], atol=1.2)
    assert actions is not None
    targets, _rotation = actions[0]
    assert targets.tolist() == [0, 1]

def test_child_origin_candidates_preserve_supercell_integer_lifts():
    basis = np.diag([1.0, 1.0, 6.0])
    translations = np.asarray([[0.0, 0.0, z] for z in range(6)])
    candidates = _child_origin_candidates(
        [0.75, 0.75, 0.25], basis, translations,
    )
    # The website's equivalent parent origin (-1/4,-1/4,9/4) maps to
    # child z=3/8.  Reducing 9/4 modulo the parent cell must not lose it.
    assert any(np.allclose(item, [0.75, 0.75, 0.375]) for item in candidates)
    assert len(candidates) == 6

def test_mode_rank_keeps_near_parallel_independent_and_drops_dependency():
    mode = SimpleNamespace()
    first = np.asarray([[1.0, 0.0, 0.0]])
    near_parallel = np.asarray([[1.0, 1e-3, 0.0]])
    dependent = first + near_parallel
    items = [
        ("a", mode, first, "a"),
        ("b", mode, near_parallel, "b"),
        ("c", mode, dependent, "c"),
    ]
    kept = _linearly_independent_mode_items(items, np.eye(3))
    assert [item[0] for item in kept] == ["a", "b"]

def test_parent_little_group_orbit_completes_vector_irrep():
    structure = Structure(Lattice.tetragonal(4.0, 5.0), ["Fe"], [[0, 0, 0]])
    seed = np.asarray([[1.0, 0.0, 0.0]])
    orbit = _parent_little_group_orbit(
        seed,
        structure,
        np.eye(3),
        123,
        np.zeros(3),
        2,
    )
    assert len(orbit) == 2
    cart = [item @ structure.lattice.matrix for item in orbit]
    assert abs(float(np.vdot(cart[0].ravel(), cart[1].ravel()))) < 1e-8

def test_self_conjugate_k_respects_parent_centering():
    assert _is_self_conjugate_parent_k(np.array([0.5, 0.5, 0.0]), 123)
    assert _is_self_conjugate_parent_k(np.array([0.5, 0.0, 0.5]), 139)
    assert not _is_self_conjugate_parent_k(np.array([0.0, 0.0, 1 / 6]), 139)
    # (0,0,1) is not a reciprocal vector of an I-centred conventional cell.
    assert not _is_self_conjugate_parent_k(np.array([0.0, 0.0, 0.5]), 139)

def test_folding_k_candidates_i4mmm_1x1x6():
    basis = np.diag([1.0, 1.0, 6.0])
    kpoints = [
        _kp("GM", ["0", "0", "0"], special=True),
        _kp("M", ["1", "1", "1"], special=True),
        _kp("X", ["1/2", "1/2", "0"], special=True),
        _kp("LD", ["0", "0", "g"], special=False, parameters=["g"]),
    ]
    primary = SimpleNamespace(
        label="LD",
        coords=np.array([0.0, 0.0, 1.0 / 6.0]),
        params=["1/6"],
        coord_tokens=("0", "0", "1/6"),
    )
    prim = _KCandidate("LD", primary.coords, ["1/6"], ("0", "0", "1/6"))
    lockin = folding_k_candidates(139, kpoints, basis, prim, nmod=0)
    labels = {(c.label, tuple(format_k_token(x) for x in c.coords)) for c in lockin}
    assert ("GM", ("0", "0", "0")) in labels
    assert ("M", ("1", "1", "1")) in labels
    assert ("LD", ("0", "0", "1/6")) in labels
    assert ("LD", ("0", "0", "1/3")) in labels
    assert ("LD", ("0", "0", "5/6")) in labels
    assert ("LD", ("0", "0", "1/2")) in labels
    assert ("X", ("1/2", "1/2", "0")) not in labels

    harm = folding_k_candidates(139, kpoints, basis, prim, nmod=1)
    harm_labels = {(c.label, tuple(format_k_token(x) for x in c.coords)) for c in harm}
    assert ("LD", ("0", "0", "1/6")) in harm_labels
    assert ("LD", ("0", "0", "1/3")) in harm_labels
    assert ("GM", ("0", "0", "0")) in harm_labels
    assert ("X", ("1/2", "1/2", "0")) not in harm_labels

def test_folding_k_candidates_uses_parent_primitive_lattice_for_fractional_basis():
    """An I-parent primitive basis is nonsingular although det(B_conv)=1/2."""
    half = 0.5
    basis = np.array([
        [-half, half, half],
        [half, -half, half],
        [half, half, -half],
    ])
    kpoints = [
        _kp("GM", ["0", "0", "0"], special=True),
        _kp("M", ["1", "1", "1"], special=True),
    ]
    folded = folding_k_candidates(139, kpoints, basis, None, nmod=0)
    assert [(item.label, item.coord_tokens) for item in folded] == [
        ("GM", ("0", "0", "0")),
    ]

def test_folding_candidates_canonicalize_i_boundary_star_before_early_stop():
    """Equivalent P-star boundary arms must not displace folded GM/M classes."""
    basis = np.array([
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 2.0],
    ])
    kpoints = [
        _kp("GM", ["0", "0", "0"], special=True),
        _kp("M", ["1", "1", "1"], special=True),
        _kp("P", ["1/2", "1/2", "1/2"], special=True),
    ]
    primary = _KCandidate(
        "P", np.array([0.5, 0.5, 0.5]), None, ("1/2", "1/2", "1/2"),
    )

    folded = folding_k_candidates(
        139, kpoints, basis, primary, nmod=0, child_centering="I",
    )

    assert [item.label for item in folded] == ["P", "GM", "M"]
    assert len(folded) == 3

def test_k_compatibility_uses_row_basis_for_axis_permutation():
    """A row-basis cell must test each child translation against k.

    NdNiO2 Method 2 uses this representative for Y=(1/3,1/2,0).  The
    old transpose convention rejected Y even though every child translation
    has an integral Bloch phase.
    """
    basis = np.array([
        [-3.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 2.0, 0.0],
    ])
    y = np.array([1.0 / 3.0, 1.0 / 2.0, 0.0])
    assert np.allclose(basis @ y, [-1.0, 0.0, 1.0])
    assert k_compatible_with_supercell(y, basis)

    incompatible = np.array([1.0 / 4.0, 1.0 / 2.0, 0.0])
    assert not k_compatible_with_supercell(incompatible, basis)

def test_folding_candidates_keep_compatible_rotated_star_arm():
    """A displayed line representative may differ from the active star arm."""
    basis = np.array([
        [-3.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 2.0, 0.0],
    ])
    kpoints = [
        _kp("DT", ["0", "b", "0"], special=False, parameters=["b"]),
    ]
    primary = _KCandidate(
        "Y", np.array([1.0 / 3.0, 1.0 / 2.0, 0.0]),
        ["1/3"], ("1/3", "1/2", "0"),
    )
    candidates = folding_k_candidates(123, kpoints, basis, primary, nmod=0)
    dt_thirds = [
        item for item in candidates
        if item.label == "DT" and item.coord_tokens == ("0", "1/3", "0")
    ]
    assert dt_thirds
    assert any(
        k_compatible_with_supercell(item.active_coords, basis)
        for item in dt_thirds
    )

def test_rotated_star_arm_keeps_physical_k_for_phase_generation():
    """The website representative must not replace a rotated arm's phase k."""
    parent = Structure(Lattice.cubic(4.0), ["O"], [[0, 0, 0]])
    raw = SmodesModeBlock(
        irrep="DT1",
        degeneracy=1,
        k_label="DT",
        k_coords=("0", "1/3", "0"),
        superlattice_cart=np.eye(3) * 4.0,
        atoms=[SmodesAtom(1, "O", np.zeros(3))],
        atom_indices=[1],
        displacements_cart=np.array([[1.0, 0.0, 0.0]]),
    )

    class _FakeSmodes:
        @staticmethod
        def compute_mode_blocks(*_args, **_kwargs):
            return [raw]

    candidate = _KCandidate(
        "DT",
        np.array([0.0, 1.0 / 3.0, 0.0]),
        ["1/3"],
        ("0", "1/3", "0"),
        active_coords=np.array([1.0 / 3.0, 0.0, 0.0]),
        arm_rotation=np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
    )
    blocks = _collect_smodes_blocks(
        _FakeSmodes(), parent, 123, [], [(candidate, ("DT", ["1/3"]))],
    )
    assert len(blocks) == 1
    assert blocks[0].k_coords == ("0", "1/3", "0")
    assert blocks[0].active_k_coords == ("1/3", "0", "0")

def test_smodes_mapping_rejects_an_uncovered_child_atom():
    """An unmatched row must not be emitted as a genuine zero displacement."""
    child = Structure(
        Lattice.cubic(10.0),
        ["Fe", "Fe"],
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )
    incomplete = SmodesModeBlock(
        irrep="X1",
        degeneracy=1,
        k_label="X",
        k_coords=("1/2", "0", "0"),
        superlattice_cart=np.eye(3) * 10.0,
        atoms=[SmodesAtom(1, "Fe", np.zeros(3))],
        atom_indices=[1],
        displacements_cart=np.array([[1.0, 0.0, 0.0]]),
    )

    assert _map_smodes_block_to_supercell(incomplete, child) is None

def test_smodes_mapping_preserves_a_covered_zero_displacement_atom():
    """Coverage is positional; a symmetry-required zero vector remains valid."""
    child = Structure(
        Lattice.cubic(10.0),
        ["Fe", "Fe"],
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
    )
    complete = SmodesModeBlock(
        irrep="X1",
        degeneracy=1,
        k_label="X",
        k_coords=("1/2", "0", "0"),
        superlattice_cart=np.eye(3) * 10.0,
        atoms=[
            SmodesAtom(1, "Fe", np.zeros(3)),
            SmodesAtom(2, "Fe", np.array([5.0, 0.0, 0.0])),
        ],
        atom_indices=[1],
        displacements_cart=np.array([[1.0, 0.0, 0.0]]),
    )

    mapped = _map_smodes_block_to_supercell(complete, child)
    assert mapped is not None
    assert np.allclose(mapped[0], [0.1, 0.0, 0.0])
    assert np.allclose(mapped[1], 0.0)

def test_smodes_k_spec_converts_official_ld_g():
    iso_kp = _kp("LD", ["0", "0", "2a"], special=False, parameters=["2a"])
    cand = _KCandidate(
        "LD", np.array([0.0, 0.0, 1.0 / 6.0]), ["1/6"], ("0", "0", "1/6"),
    )
    label, params = smodes_k_spec(139, cand, [iso_kp])
    assert label == "LD"
    assert params == ["1/12"]

def test_isoviz_mode_label_keeps_parametric_k_prefix():
    pretty = "I4/mmm[0,0,1/6]LD1(a,b)[Eu0:a:dsp]A2u(a)"
    assert _isoviz_mode_label(pretty, "LD1") == "[0,0,1/6]LD1[Eu0:a:dsp]A2u(a)"
    gamma = "I4/mmm[0,0,0]GM1+(a)[Al2:e:dsp]A1(a)"
    assert _isoviz_mode_label(gamma, "GM1+") == "GM1+[Al2:e:dsp]A1(a)"
    mpoint = "I4/mmm[1,1,1]M3-(a)[Eu0:a:dsp]A2u(a)"
    assert _isoviz_mode_label(mpoint, "M3-") == "[1,1,1]M3-[Eu0:a:dsp]A2u(a)"

def test_quadrature_shift_covers_long_primitive_supercell():
    """A fixed 0..5 parent-cell window misses k·T≈1/4 on a long P lattice."""
    basis = np.diag([1.0, 1.0, 24.0])
    k = np.array([0.0, 0.0, 1.0 / 24.0])
    shift = _quadrature_supercell_shift(k, basis, parent_sg=1)
    assert shift is not None
    phase = float(np.dot(k, np.asarray(shift, dtype=float) @ basis))
    phase -= round(phase)
    assert abs(abs(phase) - 0.25) < 0.02

def test_quadrature_shift_solves_long_period_exactly_without_a_span_cap():
    basis = np.diag([1.0, 1.0, 400.0])
    k = np.array([0.0, 0.0, 1.0 / 400.0])

    shift = _quadrature_supercell_shift(k, basis, parent_sg=1)
    assert shift is not None
    parent_translation = np.asarray(shift) @ basis
    phase = float(np.dot(k, parent_translation)) % 1.0
    assert phase == 0.25

def test_quadrature_shift_fails_closed_when_the_character_has_only_sign_phase():
    basis = np.diag([1.0, 1.0, 2.0])
    k = np.array([0.0, 0.0, 1.0 / 2.0])

    assert _quadrature_supercell_shift(k, basis, parent_sg=1) is None
    # I centering adds quarter-phase primitive translations for the same
    # conventional k point, so the centered-parent case remains resolvable.
    assert _quadrature_supercell_shift(k, basis, parent_sg=139) is not None

def test_gold_compact_m_labels_are_not_gamma():
    """Official IsoVIZ omits ``[1,1,1]`` on I4/mmm M compact labels."""
    compact = "M1+[Al2:d:dsp]B2(a)"
    full = "I4/mmm[1,1,1]M1+(a)[Al2:d:dsp]B2(a)"
    assert _mode_core_token(compact) == "1,1,1|M1+|d|B2(a)"
    assert _mode_core_token(full) == _mode_core_token(compact)
    assert _mode_core_token("M3-[Eu0:a:dsp]A2u(a)") == "1,1,1|M3-|a|A2u(a)"
    assert _mode_core_token("GM1+[Al2:e:dsp]A1(a)") == "0,0,0|GM1+|e|A1(a)"
