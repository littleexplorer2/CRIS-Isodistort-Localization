"""An unavailable parent-setting verification must not enter ISO-frame output.

The structure is a real nonstandard A-centered setting of SG 65. Only the
FINDSYM environment failure and downstream calls are mocked; no scientific
mode answer is fabricated by these boundary tests.
"""

from __future__ import annotations

import copy
import io
import zipfile
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifParser

from backend.api import IsoDistort, core_api
from backend.utils import IsodistortError, get_config
from backend.wrappers import SubgroupInfo
from features.input_cif import SymmetryValidator
from features.method2 import superspace


def _diagnostic_api(parent=None):
    if parent is None:
        parent = Structure(
            Lattice.orthorhombic(3, 4, 5), ["Fe", "Fe"],
            [[0, 0, 0], [0, 0.5, 0.5]],
        )
    api = object.__new__(IsoDistort)
    api.cfg = get_config()
    api.structure = None
    api.structure_path = None
    api.symmetry_info = None
    api._sym_val = SymmetryValidator()
    api._findsym = SimpleNamespace(standardize_wyckoff_orbits=Mock(
        side_effect=OSError("intentional FINDSYM environment failure"),
    ))
    api.set_structure(parent)
    api.distortion_scope = {}
    api.distortion_types = ["displacive"]
    api.number_of_independent_modulations = 7
    return api


def test_failed_setting_verification_keeps_the_upload_available_for_diagnostics():
    api = _diagnostic_api()
    assert api.symmetry_info["space_group_number"] == 65
    assert api.symmetry_info["parent_orbit_standardization"]["status"] == "unavailable"
    assert np.array_equal(api.structure.frac_coords, [[0, 0, 0], [0, 0.5, 0.5]])


def _invoke_mode_entry(api, entry):
    target = SimpleNamespace(index=0, _method3_route_resolution="exact_fixed_space")
    if entry == "calc":
        return api._calc_displacive_modes(65, target, ["a"])
    if entry == "scoped_displacive":
        return api._compute_scoped_modes(65, target, ["displacive"])
    if entry == "scoped_occupational":
        return api._compute_scoped_modes(65, target, ["occupational"])
    if entry == "scoped_supplied_raw":
        return api._compute_scoped_modes(65, target, ["displacive"], raw_modes=[object()])
    if entry == "select":
        return api.select_path(0)
    return api.search_method_2(0, number_of_independent_modulations=1)


@pytest.mark.parametrize("entry", [
    "calc", "scoped_displacive", "scoped_occupational", "scoped_supplied_raw",
    "select", "search",
])
def test_unavailable_setting_blocks_all_mode_entries_before_calls_or_mutation(monkeypatch, entry):
    api = _diagnostic_api()
    api._iso = Mock(name="ISO")
    api._smodes = Mock(name="SMODES")
    api._search = Mock(name="Search")
    api._dist_mapper = Mock(name="Mapper")
    api.phase_path = object()
    api._selected_subgroup = object()
    api.distortion_modes = [object()]
    api.mode_displacements = {"existing": object()}
    api.mode_displacements_sc = {"existing": object()}
    api.mode_occupancies = {"existing": object()}
    api._mode_label_overrides = {"existing": "existing label"}
    api._mode_cache_key = ("existing",)
    api.distorted_structure = api.structure.copy()
    api._generated_structure_cache_key = ("existing",)
    parametric = Mock(name="Parametric", side_effect=AssertionError("must not calculate modes"))
    special = Mock(name="Special", side_effect=AssertionError("must not calculate modes"))
    occupational = Mock(name="Occupational", side_effect=AssertionError("must not create generator"))
    monkeypatch.setattr(superspace, "compute_parametric_modes", parametric)
    monkeypatch.setattr(superspace, "compute_special_modes_with_rootless_supplement", special)
    monkeypatch.setattr(core_api, "OccupationalModeGenerator", occupational)
    before = dict(vars(api))
    before_metadata = copy.deepcopy(api.symmetry_info)

    with pytest.raises(IsodistortError, match=r"Parent setting has not been verified.*reload"):
        _invoke_mode_entry(api, entry)

    assert set(vars(api)) == set(before)
    for key, value in before.items():
        assert vars(api)[key] is value
    assert api.symmetry_info == before_metadata
    for backend in (api._iso, api._smodes, api._search, api._dist_mapper,
                    parametric, special, occupational):
        assert backend.mock_calls == []


@pytest.mark.parametrize("status", ["available", "canonicalized"])
@pytest.mark.parametrize("entry", ["calc", "scoped_occupational"])
def test_verified_status_is_not_blocked_by_private_entry_guard(monkeypatch, status, entry):
    # Empty downstream stubs prove only that the gate allows these statuses.
    # They do not purport to certify the deliberately nonstandard test parent.
    api = _diagnostic_api()
    api.symmetry_info["parent_orbit_standardization"] = {"status": status}
    api._iso = SimpleNamespace(list_k_points=Mock(return_value=[]))
    api._smodes = object()
    parametric = Mock(return_value=superspace.ParametricModeResult(
        modes=[], supercell_displacements={}, labels={}, nmod=0,
    ))
    generator = SimpleNamespace(generate=Mock(return_value=[]))
    factory = Mock(return_value=generator)
    monkeypatch.setattr(superspace, "compute_parametric_modes", parametric)
    monkeypatch.setattr(core_api, "OccupationalModeGenerator", factory)
    assert _invoke_mode_entry(api, entry) == []
    if entry == "calc":
        assert parametric.call_count == 1
        assert api._iso.list_k_points.call_count == 1
    else:
        assert factory.call_count == 1
        assert generator.generate.call_count == 1


@pytest.mark.parametrize("status", ["available", "canonicalized"])
@pytest.mark.parametrize("entry", ["select", "search"])
def test_verified_status_reaches_the_normal_public_entry_after_guard(status, entry):
    api = _diagnostic_api()
    api.symmetry_info["parent_orbit_standardization"] = {"status": status}
    downstream = RuntimeError("normal candidate listing reached")
    api.list_subgroups = Mock(side_effect=downstream)
    with pytest.raises(RuntimeError, match="normal candidate listing reached"):
        _invoke_mode_entry(api, entry)
    assert api.list_subgroups.call_count == 1


def test_legacy_missing_setting_metadata_is_not_reclassified_as_verified():
    api = _diagnostic_api()
    del api.symmetry_info["parent_orbit_standardization"]
    before = copy.deepcopy(api.symmetry_info)
    api._require_verified_parent_setting()
    assert api.symmetry_info == before
    assert "parent_orbit_standardization" not in api.symmetry_info


def _setting_boundary_candidate(parent_sg=65, size=2, group_index=8):
    return SubgroupInfo(
        index=0, space_group_number=2, space_group_symbol="P-1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0], parent_sg=parent_sg, size=size,
        subgroup_index=group_index,
    )


_SETTING_BOUNDARIES = ["search1", "search3", "zip", "disk", "spec", "collect", "supercell"]


def _invoke_setting_boundary(api, entry, dest, progress, cancel):
    candidate = _setting_boundary_candidate()
    if entry == "search1":
        return api.search_method_1(distortion_types=["occupational"])
    if entry == "search3":
        return api.search_method_3(
            distortion_types=["occupational"], space_group_type=2,
            supercell_basis=candidate.basis_vectors,
        )
    if entry == "zip":
        return api.export_subgroups_zip(
            formats=["cif"], subgroups=[candidate], compute_missing_modes=False,
            export_method=3, progress_callback=progress, cancel_check=cancel,
        )
    if entry == "disk":
        return api.export_subgroups(
            dest, formats=["cif"], subgroups=[candidate], compute_missing_modes=False,
            export_method=3, progress_callback=progress, cancel_check=cancel,
        )
    if entry == "spec":
        return api._spec_for_subgroup(
            candidate, use_current_modes=False, use_generated_structure=False,
        )
    if entry == "collect":
        return api._collect_export_specs(
            [candidate], ["cif"], False, export_method=3,
            progress_callback=progress, cancel_check=cancel,
        )
    return api._supercell_for_subgroup(candidate)


def _intercept_first_downstream(monkeypatch, api, entry, downstream):
    if entry in {"search1", "search3"}:
        monkeypatch.setattr(api, "_parent_rotations", downstream)
    elif entry == "zip":
        monkeypatch.setattr(api, "_collect_export_specs", downstream)
    elif entry == "disk":
        # The public disk entry must stop even before inspecting an export
        # destination, not only before writing or before preparing candidates.
        monkeypatch.setattr(core_api, "_published_export_folder_names", downstream)
    elif entry == "spec":
        monkeypatch.setattr(api, "_supercell_for_subgroup", downstream)
    elif entry == "collect":
        monkeypatch.setattr(api, "_snapshot_distortion_state", downstream)
    else:
        monkeypatch.setattr(core_api, "build_supercell", downstream)


@pytest.mark.parametrize("entry", _SETTING_BOUNDARIES)
def test_unavailable_setting_blocks_search_and_cif_exports_before_any_effect(
    monkeypatch, tmp_path, entry,
):
    api = _diagnostic_api()
    api._iso = Mock(name="ISO")
    api._smodes = Mock(name="SMODES")
    api._search = Mock(name="Search")
    api._dist_mapper = Mock(name="Mapper")
    downstream = Mock(side_effect=AssertionError("unverified setting reached downstream"))
    _intercept_first_downstream(monkeypatch, api, entry, downstream)
    progress = Mock(side_effect=AssertionError("must not call progress"))
    cancel = Mock(side_effect=AssertionError("must not call cancellation observer"))
    dest = tmp_path / "never-created"
    before = dict(vars(api))
    before_metadata = copy.deepcopy(api.symmetry_info)
    before_structure = api.structure.copy()

    with pytest.raises(IsodistortError, match=r"Parent setting has not been verified.*reload"):
        _invoke_setting_boundary(api, entry, dest, progress, cancel)

    assert not dest.exists()
    assert list(tmp_path.iterdir()) == []
    assert set(vars(api)) == set(before)
    assert all(vars(api)[key] is value for key, value in before.items())
    assert api.symmetry_info == before_metadata
    assert api.structure == before_structure
    for observer in (downstream, progress, cancel, api._iso, api._smodes, api._search, api._dist_mapper):
        assert observer.mock_calls == []


@pytest.mark.parametrize("status", ["available", "canonicalized"])
@pytest.mark.parametrize("entry", _SETTING_BOUNDARIES)
def test_verified_setting_reaches_the_normal_search_or_export_boundary(
    monkeypatch, tmp_path, status, entry,
):
    # Deliberate downstream stops prove gate reachability only; assigning a
    # status here is not a scientific certification of the diagnostic parent.
    api = _diagnostic_api()
    api.symmetry_info["parent_orbit_standardization"] = {"status": status}
    downstream = Mock(side_effect=RuntimeError("normal downstream reached"))
    _intercept_first_downstream(monkeypatch, api, entry, downstream)
    progress = Mock()
    cancel = Mock()

    with pytest.raises(RuntimeError, match="normal downstream reached"):
        _invoke_setting_boundary(api, entry, tmp_path / "never-created", progress, cancel)

    assert downstream.call_count == 1
    assert progress.mock_calls == []
    assert cancel.mock_calls == []
    assert list(tmp_path.iterdir()) == []


def test_unverified_shifted_parent_cannot_publish_a_cif_that_adds_inversion_images():
    """A valid origin-shifted parent is not in the standard ISO subgroup frame.

    Pmmm has one Ni at this shifted inversion center. If the unchecked writer
    applies standard P-1 at O=0, its CIF creates a second Ni at -x. Only the
    FINDSYM runtime failure is simulated; structure detection and any reached
    ZIP/CIF rendering and parsing are real, with no filesystem publication.
    """
    parent = Structure(Lattice.orthorhombic(3, 4, 5), ["Ni"], [[0.1, 0.2, 0.3]])
    api = _diagnostic_api(parent)
    assert api.symmetry_info["space_group_number"] == 47
    assert api.symmetry_info["parent_orbit_standardization"]["status"] == "unavailable"
    candidate = _setting_boundary_candidate(parent_sg=47, size=1, group_index=4)
    try:
        archive = api.export_subgroups_zip(
            formats=["cif"], subgroups=[candidate], compute_missing_modes=False,
            export_method=3,
        )
    except IsodistortError as exc:
        assert "Parent setting has not been verified" in str(exc)
        assert "reload" in str(exc)
        assert api.structure is parent
        return

    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        cif_names = [name for name in zipped.namelist() if name.endswith("subgroup.cif")]
        assert len(cif_names) == 1
        text = zipped.read(cif_names[0]).decode("utf-8")
    daughter = CifParser.from_str(text).parse_structures(primitive=False)[0]
    pytest.fail(
        "Unverified parent frame published a CIF: "
        f"{parent.composition} ({len(parent)} atoms) -> "
        f"{daughter.composition} ({len(daughter)} atoms)"
    )
