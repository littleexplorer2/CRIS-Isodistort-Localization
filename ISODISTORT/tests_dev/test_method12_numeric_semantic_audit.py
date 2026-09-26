"""Regression tests for the signed Method 1/2 numeric semantic audit."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from tests_dev.manual import audit_method12_numeric_semantics as numeric


def _dataset(
    atoms: list[numeric.Atom],
    columns: dict[str, np.ndarray],
    *,
    normfactor: float | None = None,
    bound: float | None = None,
) -> numeric.ModeDataset:
    lattice = np.diag([4.0, 5.0, 6.0])
    return numeric.ModeDataset(
        atoms,
        lattice,
        {
            label: numeric.Mode(
                label,
                np.asarray(vector, dtype=float),
                normfactor=normfactor,
                amplitude_bound=bound,
            )
            for label, vector in columns.items()
        },
    )


def test_global_sign_and_degenerate_rotation_are_semantically_equivalent() -> None:
    atoms = [
        numeric.Atom("Eu", np.array([0.0, 0.0, 0.0])),
        numeric.Atom("Al", np.array([0.5, 0.5, 0.5])),
    ]
    eu_x = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])
    eu_y = np.array([[0.0, 0.2, 0.0], [0.0, 0.0, 0.0]])
    al_z = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0 / 6.0]])
    official = _dataset(
        atoms,
        {
            "x1[eu:dsp]eu(a)": eu_x,
            "x1[eu:dsp]eu(b)": eu_y,
            "gm1+[al:dsp]a1(a)": al_z,
        },
    )
    # Reorder atoms, rotate the two-dimensional Eu basis, and reverse the
    # independent Al mode.  All operations leave the physical subspaces fixed.
    reordered = [atoms[1], atoms[0]]
    local = _dataset(
        reordered,
        {
            "x1[eu:dsp]eu(a)": np.array([[0, 0, 0], [0.25 / np.sqrt(2), 0.2 / np.sqrt(2), 0]]),
            "x1[eu:dsp]eu(b)": np.array([[0, 0, 0], [0.25 / np.sqrt(2), -0.2 / np.sqrt(2), 0]]),
            "gm1+[al:dsp]a1(a)": np.array([[0, 0, -1.0 / 6.0], [0, 0, 0]]),
        },
    )

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "pass"
    assert result["max_principal_angle_degrees"] < 1.0e-5


def test_missing_mode_is_not_hidden_by_common_label_comparison() -> None:
    atoms = [numeric.Atom("Eu", np.array([0.0, 0.0, 0.0]))]
    x = np.array([[0.25, 0.0, 0.0]])
    y = np.array([[0.0, 0.2, 0.0]])
    official = _dataset(atoms, {"x1[eu:dsp]eu(a)": x, "x1[eu:dsp]eu(b)": y})
    local = _dataset(atoms, {"x1[eu:dsp]eu(a)": x})

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "fail"
    assert result["overall_cartesian_subspace"]["official_rank"] == 2
    assert result["overall_cartesian_subspace"]["local_rank"] == 1


def test_component_names_are_noncanonical_within_fixed_site_irrep_family() -> None:
    atoms = [numeric.Atom("O", np.array([0.0, 0.0, 0.0]))]
    x = np.array([[0.25, 0.0, 0.0]])
    y = np.array([[0.0, 0.2, 0.0]])
    official = _dataset(
        atoms,
        {
            "dt1[o:f:dsp]b2u(a)": x,
            "dt1[o:f:dsp]b3u(b)": y,
        },
    )
    local = _dataset(
        atoms,
        {
            "dt1[o:f:dsp]b2u(d)": -x,
            "dt1[o:f:dsp]b3u(c)": y,
        },
    )

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "pass"
    assert not result["missing_local_labels"]
    assert not result["extra_local_labels"]
    assert len(result["component_copy_label_reassignments"]) == 2


def test_site_irrep_name_is_not_relaxed_as_component_permutation() -> None:
    atoms = [numeric.Atom("O", np.array([0.0, 0.0, 0.0]))]
    vector = np.array([[0.25, 0.0, 0.0]])
    official = _dataset(atoms, {"dt1[o:f:dsp]b2u(a)": vector})
    local = _dataset(atoms, {"dt1[o:f:dsp]b3u(a)": vector})

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "fail"
    assert result["missing_local_labels"] == ["dt1[o:f:dsp]b2u(a)"]
    assert result["extra_local_labels"] == ["dt1[o:f:dsp]b3u(a)"]


def test_optional_compact_k_prefix_is_removed_only_when_unambiguous() -> None:
    atoms = [numeric.Atom("Eu", np.array([0.0, 0.0, 0.0]))]
    vector = np.array([[0.25, 0.0, 0.0]])
    official = _dataset(atoms, {"x4-[eu:dsp]eu(a)": vector})
    local = _dataset(atoms, {"[1/2,1/2,0]x4-[eu:dsp]eu(a)": -vector})

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "pass"
    assert not result["missing_local_labels"]
    assert not result["extra_local_labels"]


def test_atom_matching_is_species_preserving_and_periodic() -> None:
    lattice = np.diag([4.0, 4.0, 4.0])
    left = [
        numeric.Atom("Eu", np.array([0.0, 0.0, 0.0])),
        numeric.Atom("Al", np.array([0.25, 0.5, 0.75])),
    ]
    right = [
        numeric.Atom("Al", np.array([1.25, -0.5, 0.75])),
        numeric.Atom("Eu", np.array([1.0, 0.0, 0.0])),
    ]

    mapping, maximum = numeric._atom_bijection(left, right, lattice, 1.0e-8)

    assert mapping == [1, 0]
    assert maximum < 1.0e-12


def test_atom_matching_can_prove_one_global_origin_shift() -> None:
    lattice = np.diag([4.0, 5.0, 6.0])
    left = [
        numeric.Atom("Eu", np.array([0.0, 0.0, 0.0])),
        numeric.Atom("Al", np.array([0.25, 0.5, 0.75])),
    ]
    shift = np.array([0.25, -0.25, 0.5])
    right = [
        numeric.Atom("Al", np.mod(left[1].frac + shift, 1.0)),
        numeric.Atom("Eu", np.mod(left[0].frac + shift, 1.0)),
    ]

    mapping, maximum, inferred = numeric._atom_bijection_with_origin_shift(
        left,
        right,
        lattice,
        1.0e-8,
    )

    assert mapping == [1, 0]
    assert maximum < 1.0e-12
    delta = np.asarray(inferred) - shift
    assert np.allclose(delta - np.rint(delta), 0.0)


def test_checkpoint_rejects_changed_signature_and_atomic_write_leaves_no_tmp(
    tmp_path: Path,
) -> None:
    checkpoint = tmp_path / "checkpoint.json"
    payload = {
        "schema": numeric.CHECKPOINT_SCHEMA,
        "signature": "first",
        "created_at": "time",
        "updated_at": "time",
        "results": {},
    }
    numeric._atomic_json(checkpoint, payload)

    assert json.loads(checkpoint.read_text(encoding="utf-8"))["signature"] == "first"
    assert not checkpoint.with_name(f"{checkpoint.name}.tmp").exists()
    with pytest.raises(RuntimeError, match="signature mismatch"):
        numeric._load_checkpoint(checkpoint, "second", restart=False)


def test_local_modes_parser_ignores_prose_mode_headings(tmp_path: Path) -> None:
    path = tmp_path / "Complete modes details.txt"
    path.write_text(
        """Lattice parameters of the supercell:
  a=4 b=4 c=4 alpha=90 beta=90 gamma=90
Superstructure in the traditional atomic-xyz-coordinate basis
  atom el x y z
     1 Eu 0 0 0
Mode definitions (every atom in the unit cell)
Mode vectors are given in fractional coordinates
Mode X1__a__A1(a)  P1[0,0,0]X1(a)[Eu:a:dsp]A1(a)
  normfactor = 0.25 Angstrom^-1
  atom 1 Eu (1, 0, 0)
Mode amplitudes (As)
""",
        encoding="utf-8",
    )

    parsed = numeric.parse_local_modes_text(path)

    assert len(parsed.modes) == 1
    assert next(iter(parsed.modes.values())).normfactor == pytest.approx(0.25)


def test_local_amplitude_triplet_uses_scientific_ap_and_dmax_formula(tmp_path: Path) -> None:
    path = tmp_path / "Complete modes details.txt"
    path.write_text(
        """Complete modes details
Basis = {(1,0,0),(0,1,0),(0,0,1)} origin=(0,0,0) s=4, i=4
Lattice parameters of the supercell:
  a=4 b=5 c=6 alpha=90 beta=90 gamma=90
Superstructure in the traditional atomic-xyz-coordinate basis
  atom el x y z
     1 Eu 0 0 0
Mode definitions (every atom in the unit cell)
Mode X1__a__A1(a)  P1[0,0,0]X1(a)[Eu:a:dsp]A1(a)
  normfactor = 0.25 Angstrom^-1
  atom 1 Eu (1, 0, 0)
Mode amplitudes (As = supercell-normalized, Ap = parent-cell-normalized)
  X1__a__A1(a)  As=2.000000  Ap=1.000000  dmax=2.000000 Angstrom
""",
        encoding="utf-8",
    )

    parsed = numeric.parse_local_modes_text(path)
    result = numeric._amplitude_formula_audit(parsed)

    assert result["status"] == "pass"
    assert result["nonzero_amplitude_count"] == 1
    assert result["max_Ap_formula_abs_error"] == pytest.approx(0.0)
    assert result["max_dmax_formula_abs_error_angstrom"] == pytest.approx(0.0)


def test_duplicate_pretty_labels_preserve_every_exported_column(tmp_path: Path) -> None:
    path = tmp_path / "Complete modes details.txt"
    path.write_text(
        """Lattice parameters of the supercell:
  a=4 b=4 c=4 alpha=90 beta=90 gamma=90
Superstructure in the traditional atomic-xyz-coordinate basis
  atom el x y z
     1 Nd 0 0 0
Mode definitions (every atom in the unit cell)
Mode DT1__d__Eu(a)  P4/mmm[0,1/3,0]DT1(a,b,c,d)[Nd:d:dsp]Eu(a)
  normfactor = 0.25 Angstrom^-1
  atom 1 Nd (1, 0, 0)
Mode DT1__d__Eu(a)_2  P4/mmm[0,1/3,0]DT1(a,b,c,d)[Nd:d:dsp]Eu(a)
  normfactor = 0.25 Angstrom^-1
  atom 1 Nd (0, 1, 0)
""",
        encoding="utf-8",
    )

    parsed = numeric.parse_local_modes_text(path)

    assert len(parsed.modes) == 2
    assert any("#duplicate2" in label for label in parsed.modes)


def test_mode_label_diagnostics_exposes_collision_without_losing_columns() -> None:
    labels = {
        "DT1__d__Eu(a)": "P4/mmm[0,1/3,0]DT1(a,b,c,d)[Nd:d:dsp]Eu(a)",
        "DT1__d__Eu(a)_2": "P4/mmm[0,1/3,0]DT1(a,b,c,d)[Nd:d:dsp]Eu(a)",
        "DT1__d__Eu(b)": "P4/mmm[0,1/3,0]DT1(a,b,c,d)[Nd:d:dsp]Eu(b)",
    }

    result = numeric._mode_label_diagnostics(labels)

    assert result["internal_mode_count"] == 3
    assert result["unique_canonical_pretty_label_count"] == 2
    assert result["collision_group_count"] == 1
    assert result["colliding_column_count"] == 2


def test_amplitude_bound_swap_is_accepted_as_family_basis_permutation() -> None:
    atoms = [numeric.Atom("Al", np.array([0.0, 0.0, 0.0]))]
    x = np.array([[0.25, 0.0, 0.0]])
    y = np.array([[0.0, 0.2, 0.0]])
    official = _dataset(
        atoms,
        {
            "ld1[al:e:dsp]a1_1(a)": x,
            "ld1[al:e:dsp]a1_2(a)": y,
        },
    )
    local = _dataset(
        atoms,
        {
            "ld1[al:e:dsp]a1_1(a)": x,
            "ld1[al:e:dsp]a1_2(a)": y,
        },
    )
    official.modes["ld1[al:e:dsp]a1_1(a)"].amplitude_bound = 2.0
    official.modes["ld1[al:e:dsp]a1_2(a)"].amplitude_bound = 3.0
    local.modes["ld1[al:e:dsp]a1_1(a)"].amplitude_bound = 3.0
    local.modes["ld1[al:e:dsp]a1_2(a)"].amplitude_bound = 2.0

    result = numeric.compare_mode_datasets(
        official,
        local,
        compare_amplitude_bounds=True,
    )

    assert result["status"] == "pass"
    assert result["amplitude_bound_family_multiset"]["permutation_equivalent"] is True
    assert result["amplitude_bound_family_multiset"]["max_abs_error"] == pytest.approx(0.0)
    assert "allowed component/copy basis permutation" in result["diagnostics"][0]


def test_amplitude_bound_family_multiset_difference_still_fails() -> None:
    atoms = [numeric.Atom("Al", np.array([0.0, 0.0, 0.0]))]
    x = np.array([[0.25, 0.0, 0.0]])
    y = np.array([[0.0, 0.2, 0.0]])
    official = _dataset(
        atoms,
        {
            "ld1[al:e:dsp]a1_1(a)": x,
            "ld1[al:e:dsp]a1_2(a)": y,
        },
    )
    local = _dataset(
        atoms,
        {
            "ld1[al:e:dsp]a1_1(a)": x,
            "ld1[al:e:dsp]a1_2(a)": y,
        },
    )
    official.modes["ld1[al:e:dsp]a1_1(a)"].amplitude_bound = 2.0
    official.modes["ld1[al:e:dsp]a1_2(a)"].amplitude_bound = 3.0
    local.modes["ld1[al:e:dsp]a1_1(a)"].amplitude_bound = 2.0
    local.modes["ld1[al:e:dsp]a1_2(a)"].amplitude_bound = 4.0

    result = numeric.compare_mode_datasets(
        official,
        local,
        compare_amplitude_bounds=True,
    )

    assert result["status"] == "fail"
    assert result["amplitude_bound_family_multiset"]["permutation_equivalent"] is False
    assert result["amplitude_bound_family_multiset"]["max_abs_error"] == pytest.approx(1.0)


def test_live_candidate_selector_uses_full_identity_and_periodic_origin() -> None:
    expected = {
        "irrep": "X4-",
        "opd": "C1",
        "space_group_number": 49,
        "basis": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
        "origin": [0, 0.5, 0],
    }
    correct = SimpleNamespace(
        index=3,
        k_point_label="X",
        k_parameters=[],
        irrep_label="X4-",
        opd_symbol="C1",
        space_group_number=49,
        space_group_symbol="Pccm",
        basis_vectors=((1, 1, 0), (-1, 1, 0), (0, 0, 1)),
        origin=(0, 1.5, 0),
        supercell_size=4,
        subgroup_index=8,
    )
    wrong_setting = SimpleNamespace(**{**vars(correct), "index": 4, "origin": (0, 0.25, 0)})

    selected = numeric._select_live_candidate([wrong_setting, correct], expected)

    assert selected is correct


def test_live_generation_cache_is_content_addressed(tmp_path: Path) -> None:
    case_root = tmp_path / "CASE"
    export_dir = case_root / "runs" / "one" / "exports" / "candidate"
    export_dir.mkdir(parents=True)
    exported = export_dir / "subgroup.cif"
    exported.write_text("first", encoding="utf-8")
    manifest_path = case_root / "runs" / "one" / "manifest.json"
    manifest = {
        "status": "complete",
        "source_signature": "source",
        "export_directory": str(export_dir),
        "export_files": [
            {
                "path": str(exported),
                "sha256": numeric._sha256(exported),
                "bytes": exported.stat().st_size,
            }
        ],
    }
    numeric._atomic_json(manifest_path, manifest)
    numeric._atomic_json(case_root / "latest.json", {"manifest_path": str(manifest_path)})

    assert numeric._read_valid_live_generation(case_root, "source") is not None
    exported.write_text("tampered", encoding="utf-8")
    assert numeric._read_valid_live_generation(case_root, "source") is None


def test_live_directory_audit_preserves_explicit_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    official = tmp_path / "official"
    local = tmp_path / "live"
    official.mkdir()
    local.mkdir()
    for name in numeric.CORE_FILENAMES:
        (official / name).write_text("placeholder", encoding="utf-8")
        (local / name).write_text("placeholder", encoding="utf-8")
    (official / "Complete modes details.html").write_text("placeholder", encoding="utf-8")
    (local / "Complete modes details.txt").write_text("placeholder", encoding="utf-8")
    atoms = [numeric.Atom("Eu", np.array([0.0, 0.0, 0.0]))]
    dataset = _dataset(atoms, {"gm1+[eu:dsp]a(a)": np.array([[0.25, 0.0, 0.0]])})
    passed = {"status": "pass", "reasons": []}
    monkeypatch.setattr(numeric, "parse_official_modes_html", lambda _path: dataset)
    monkeypatch.setattr(numeric, "parse_local_modes_text", lambda _path: dataset)
    monkeypatch.setattr(numeric, "parse_isoviz", lambda _path, _lattice: dataset)
    monkeypatch.setattr(numeric, "parse_topas", lambda _path, _lattice: dataset)
    monkeypatch.setattr(numeric, "compare_mode_datasets", lambda *_args, **_kwargs: passed)
    monkeypatch.setattr(numeric, "compare_cif", lambda *_args, **_kwargs: passed)
    case = {
        "id": "LIVE",
        "coverage": ["live-current-source"],
        "parent": "parent.cif",
        "method": "Method2",
        "folder": "candidate",
        "nmod": 0,
    }
    provenance = {"kind": "live_current_source", "source_signature": "abc"}

    result = numeric._audit_directories(
        case,
        official,
        local,
        evidence_scope="generated now",
        provenance=provenance,
    )

    assert result["status"] == "pass"
    assert result["provenance"] == provenance
    assert result["evidence_scope"] == "generated now"
