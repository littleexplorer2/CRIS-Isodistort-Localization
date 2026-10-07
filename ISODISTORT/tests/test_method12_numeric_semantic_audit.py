"""Regression tests for the signed Method 1/2 numeric semantic audit."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from tests.manual import audit_method12_numeric_semantics as numeric


def _official_complete_modes_html(preformatted: str = "") -> str:
    title = "ISODISTORT: complete modes details"
    return (
        f"<html><head><title>{title}</title></head><body><h1>{title}</h1>"
        "<strong>Subgroup details</strong>"
        "<strong>Undistorted superstructure</strong>"
        "<strong>Distorted superstructure</strong>"
        f"<pre>{preformatted}</pre></body></html>"
    )


def _audit_case(case_id: str = "TEST", folder: str = "candidate") -> dict[str, object]:
    return {
        "id": case_id,
        "coverage": ["official-html-resolution"],
        "parent": "parent.cif",
        "method": "Method2",
        "folder": folder,
        "nmod": 0,
    }


def _pass_record(case: dict[str, object]) -> dict[str, object]:
    exports = {
        name: {"status": "pass", "reasons": []}
        for name in ("complete_modes", "isoviz", "topas", "cif")
    }
    return {
        "id": case["id"],
        "status": "pass",
        "coverage": case["coverage"],
        "exports": exports,
        "failed_exports": [],
        "inconclusive_exports": [],
        "resumed_from_checkpoint": False,
    }


def _saved_audit_tree(
    root: Path,
) -> tuple[dict[str, object], Path, Path, str]:
    case = _audit_case()
    official, local = numeric._case_paths(root, case)
    official.mkdir(parents=True)
    local.mkdir(parents=True)
    for name in numeric.CORE_FILENAMES:
        (official / name).write_text(f"official {name}", encoding="utf-8")
        (local / name).write_text(f"local {name}", encoding="utf-8")
    (local / "Complete modes details.txt").write_text("local modes", encoding="utf-8")
    html = _official_complete_modes_html("captured A")
    html_path = official / "renamed evidence.htm"
    html_path.write_text(html, encoding="utf-8")
    return case, official, html_path, html


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


def test_atomic_json_serializes_numpy_scalars_and_nonfinite_values(tmp_path: Path) -> None:
    path = tmp_path / "report.json"

    numeric._atomic_json(
        path,
        {
            "ok": np.bool_(True),
            "count": np.int64(3),
            "finite": np.float64(1.25),
            "nonfinite": np.float64(np.inf),
        },
    )

    assert json.loads(path.read_text(encoding="utf-8")) == {
        "count": 3,
        "finite": 1.25,
        "nonfinite": None,
        "ok": True,
    }


def test_topas_parser_composes_shared_displacement_aliases(tmp_path: Path) -> None:
    official = tmp_path / "official.str"
    local = tmp_path / "local.str"
    header = "prm !a1 0 min -2 max 2 'X4+[Al1:e:dsp]A1(a)\n"
    official.write_text(
        header
        + "prm Al1_dx = + 0.25*a1;: 0\n"
        + "prm Al1_x = 0 + Al1_dx;: 0\n"
        + "prm Al1_y = 0 - Al1_dx;: 0\n"
        + "prm !Al1_z = 0;: 0\n",
        encoding="utf-8",
    )
    local.write_text(
        header
        + "prm Al1_dx = + 0.25*a1;: 0\n"
        + "prm Al1_dy = - 0.25*a1;: 0\n"
        + "prm Al1_x = 0 + Al1_dx;: 0\n"
        + "prm Al1_y = 0 + Al1_dy;: 0\n"
        + "prm !Al1_z = 0;: 0\n",
        encoding="utf-8",
    )

    left = numeric.parse_topas(official, np.eye(3))
    right = numeric.parse_topas(local, np.eye(3))

    assert len(left.atoms) == len(right.atoms) == 1
    assert np.array_equal(
        next(iter(left.modes.values())).vectors,
        next(iter(right.modes.values())).vectors,
    )


def test_topas_different_asu_representatives_are_proven_through_complete_modes() -> None:
    label = "x1[al:e:dsp]a1(a)"
    complete_atoms = [
        numeric.Atom("Al", np.array([0.0, 0.0, 0.0])),
        numeric.Atom("Al", np.array([0.5, 0.0, 0.0])),
        numeric.Atom("Eu", np.array([0.25, 0.0, 0.0])),
    ]
    complete_vector = np.array(
        [[0.25, 0.0, 0.0], [0.25, 0.0, 0.0], [0.1, 0.0, 0.0]]
    )
    official_complete = _dataset(complete_atoms, {label: complete_vector})
    local_complete = _dataset(complete_atoms, {label: complete_vector})
    official_topas = _dataset(
        [complete_atoms[0], complete_atoms[2]],
        {label: complete_vector[[0, 2]]},
        bound=2.0,
    )
    local_topas = _dataset(
        [complete_atoms[1], complete_atoms[2]],
        {label: complete_vector[[1, 2]]},
        bound=2.0,
    )

    direct = numeric.compare_mode_datasets(
        official_topas,
        local_topas,
        compare_amplitude_bounds=True,
    )
    result = numeric.combine_topas_evidence(
        direct,
        numeric.compare_topas_to_complete_modes(official_topas, official_complete),
        numeric.compare_topas_to_complete_modes(local_topas, local_complete),
        numeric.compare_mode_datasets(official_complete, local_complete),
    )

    assert direct["status"] == "inconclusive"
    assert result["status"] == "pass"
    assert result["transitive_complete_modes_proof"] is True
    assert result["comparison_route"] == "same-side-TOPAS-to-Complete-modes-transitive"


def test_exact_setting_rotation_origin_and_vector_map_are_semantically_equivalent() -> None:
    official_lattice = np.diag([4.0, 5.0, 6.0])
    matrix = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=float)
    inverse = np.linalg.inv(matrix)
    shift = np.array([0.25, 0.5, 0.0])
    official_atoms = [
        numeric.Atom("La", np.array([0.1, 0.2, 0.3])),
        numeric.Atom("O", np.array([0.6, 0.4, 0.8])),
    ]
    official_vector = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])
    label = "I4/mmm[1/2,0,1/2]N1+[La1:e:dsp]A1(a)"
    official = numeric.ModeDataset(
        official_atoms,
        official_lattice,
        {label: numeric.Mode(label, official_vector)},
    )
    # x_official = x_local @ U + q and d_official = d_local @ U.
    local_atoms_in_official_order = [
        numeric.Atom(atom.species, np.mod((atom.frac - shift) @ inverse, 1.0))
        for atom in official_atoms
    ]
    local = numeric.ModeDataset(
        list(reversed(local_atoms_in_official_order)),
        matrix @ official_lattice,
        {
            label: numeric.Mode(
                label,
                np.asarray([[0.0, 0.0, 0.0], *(-official_vector @ inverse)[:1]]),
            )
        },
    )
    transform = numeric.DatasetSettingTransform(
        matrix=tuple(tuple(numeric.Fraction(int(value)) for value in row) for row in matrix),
        origin_shift=tuple(numeric.Fraction(str(value)) for value in shift),
    )

    result = numeric.compare_mode_datasets(
        official,
        local,
        setting_transform=transform,
    )

    assert result["status"] == "pass"
    assert result["max_species_periodic_mapping_error_angstrom"] < 1.0e-12
    assert result["max_principal_angle_degrees"] < 1.0e-5
    assert result["setting_transform"]["local_to_official_matrix"] == [
        ["0", "1", "0"],
        ["-1", "0", "0"],
        ["0", "0", "1"],
    ]


def test_raw_normfactor_is_not_compared_across_sheared_setting() -> None:
    atoms = [numeric.Atom("Al", np.array([0.0, 0.0, 0.0]))]
    official_lattice = np.eye(3)
    label = "X1[Al:a:dsp]A1(a)"
    official = numeric.ModeDataset(
        atoms,
        official_lattice,
        {
            label: numeric.Mode(
                label,
                np.array([[1.0, -1.0, 0.0]]),
                normfactor=1.0 / np.sqrt(2.0),
            )
        },
    )
    matrix = np.array([[1, 2, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
    local = numeric.ModeDataset(
        atoms,
        matrix @ official_lattice,
        {
            label: numeric.Mode(
                label,
                np.array([[1.0 / 3.0, -1.0, 0.0]]),
                normfactor=3.0 / np.sqrt(2.0),
            )
        },
    )
    transform = numeric.DatasetSettingTransform(
        matrix=tuple(
            tuple(numeric.Fraction(int(value)) for value in row) for row in matrix
        ),
        origin_shift=(numeric.Fraction(0),) * 3,
    )

    result = numeric.compare_mode_datasets(
        official,
        local,
        setting_transform=transform,
        compare_normfactors=True,
    )

    assert result["status"] == "pass"
    assert result["max_gram_eigenvalue_abs_error"] == pytest.approx(0.0)
    assert result["max_normfactor_abs_error_angstrom_inverse"] is None
    assert result["normfactor_comparison"]["status"] == "skipped"


def test_duplicate_column_is_not_hidden_by_unchanged_rank() -> None:
    atoms = [numeric.Atom("La", np.array([0.0, 0.0, 0.0]))]
    vector = np.array([[0.25, 0.0, 0.0]])
    label = "I4/mmm[1/2,0,1/2]N1+[La1:e:dsp]A1(a)"
    official = _dataset(atoms, {label: vector})
    local = _dataset(
        atoms,
        {
            label: vector,
            f"{label}#duplicate2": vector,
        },
    )

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "fail"
    assert result["overall_cartesian_subspace"]["official_rank"] == 1
    assert result["overall_cartesian_subspace"]["local_rank"] == 1
    assert result["group_multiplicity_mismatches"]


def test_swapped_physical_orbits_fail_even_when_whole_space_is_unchanged() -> None:
    atoms = [
        numeric.Atom("La", np.array([0.0, 0.0, 0.0])),
        numeric.Atom("La", np.array([0.0, 0.0, 0.25])),
    ]
    first = np.array([[0.25, 0.0, 0.0], [0.0, 0.0, 0.0]])
    second = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    label1 = "I4/mmm[1/2,0,1/2]N1+[La1:e:dsp]A1(a)"
    label2 = "I4/mmm[1/2,0,1/2]N1+[La2:e:dsp]A1(a)"
    official = _dataset(atoms, {label1: first, label2: second})
    local = _dataset(atoms, {label1: second, label2: first})

    result = numeric.compare_mode_datasets(official, local)

    assert result["status"] == "fail"
    assert result["overall_cartesian_subspace"]["max_principal_angle_degrees"] < 1.0e-5
    assert result["max_principal_angle_degrees"] == pytest.approx(90.0)


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


def test_live_api_uses_matching_displacive_context_and_validation_cache(
    tmp_path: Path,
) -> None:
    class FakeWrapper:
        _mode = "native"
        _stage_dir = "old"

    class FakeApi:
        def __init__(self) -> None:
            self._iso = FakeWrapper()
            self.distortion_types: list[str] = []

        def set_distortion_types(self, values: list[str]) -> None:
            self.distortion_types = list(values)

    api = FakeApi()
    cache = tmp_path / "validation" / "isotropy_cache"

    state = numeric._configure_live_api(api, cache)

    assert api.distortion_types == ["displacive"]
    assert Path(api._iso._stage_dir) == cache.resolve()
    assert state == {
        "distortion_types": ["displacive"],
        "cache_directory": str(cache.resolve()),
        "backend_stage_directory": str(cache.resolve()),
    }


@pytest.mark.parametrize("route", ["list_subgroups_at_kpoint", "list_subgroups_at"])
def test_live_candidate_inventory_explicitly_generates_missing_parameter_db(
    route: str,
) -> None:
    calls: list[tuple[str, tuple, dict]] = []

    class FakeApi:
        def list_subgroups_at_kpoint(self, *args, **kwargs):
            calls.append(("list_subgroups_at_kpoint", args, kwargs))
            return ["candidate"]

        def list_subgroups_at(self, *args, **kwargs):
            calls.append(("list_subgroups_at", args, kwargs))
            return ["candidate"]

    selector = {
        "route": route,
        "k_point": "Y",
        "k_parameters": ["1/3"],
        "irrep": "Y1",
        "opd": "C1",
    }

    assert numeric._list_live_candidates(FakeApi(), selector) == ["candidate"]
    assert calls[0][0] == route
    assert calls[0][2]["generate_if_missing"] is True


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
    (official / "manually renamed evidence.hTm").write_text(
        _official_complete_modes_html("placeholder"),
        encoding="utf-8",
    )
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
    monkeypatch.setattr(
        numeric,
        "setting_transform_from_cifs",
        lambda *_args: numeric.DatasetSettingTransform(
            matrix=(
                (numeric.Fraction(1), numeric.Fraction(0), numeric.Fraction(0)),
                (numeric.Fraction(0), numeric.Fraction(1), numeric.Fraction(0)),
                (numeric.Fraction(0), numeric.Fraction(0), numeric.Fraction(1)),
            ),
            origin_shift=(numeric.Fraction(0),) * 3,
        ),
    )
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
        official_html_inventory=numeric.scan_official_html(official, recursive=False),
    )

    assert result["status"] == "pass"
    assert result["provenance"] == provenance
    assert result["evidence_scope"] == "generated now"
    assert result["official_html"]["selected_page"]["relative_path"] == (
        "manually renamed evidence.hTm"
    )


def test_official_parser_uses_captured_content_after_arbitrary_htm_is_changed(
    tmp_path: Path,
) -> None:
    path = tmp_path / "任意手动文件名.hTm"
    captured_html = _official_complete_modes_html(
        """
a=4, b=4, c=4, alpha=90, beta=90, gamma=90
Displacive mode definitions
P1[0,0,0]GM1+[Eu:a:dsp]A1(a) normfactor = 0.25
Eu1 0 0 0 1 0 0
Parent-cell strain mode definitions
"""
    )
    path.write_bytes(captured_html.encode("utf-8"))
    inventory = numeric.scan_official_html(tmp_path, recursive=False)
    page = inventory.require_unique("complete_modes_details")

    path.write_text("<html><title>changed after capture</title></html>", encoding="utf-8")
    parsed = numeric.parse_official_modes_html(page)

    assert page.raw_html == captured_html
    assert len(parsed.modes) == 1
    assert next(iter(parsed.modes.values())).normfactor == pytest.approx(0.25)


def test_audit_reports_duplicate_complete_modes_pages_as_structured_inconclusive(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    official = tmp_path / "official"
    official.mkdir()
    (official / "first arbitrary.html").write_text(
        _official_complete_modes_html(), encoding="utf-8"
    )
    (official / "second arbitrary.hTm").write_text(
        _official_complete_modes_html(), encoding="utf-8"
    )
    monkeypatch.setattr(
        numeric,
        "parse_official_modes_html",
        lambda _page: pytest.fail("an ambiguous inventory must never select a page"),
    )

    result = numeric._audit_directories(
        _audit_case(),
        official,
        tmp_path / "local",
        evidence_scope="test",
        provenance={"kind": "test"},
        official_html_inventory=numeric.scan_official_html(official, recursive=False),
    )

    assert result["status"] == "inconclusive"
    assert result["reason_code"] == "official_html_ambiguous_match"
    resolution = result["official_html"]["resolution"]
    assert resolution["code"] == "ambiguous_match"
    assert len(resolution["matched_paths"]) == 2
    assert result["official_html"]["inventory"]["page_count"] == 2


def test_audit_rejects_misleading_complete_modes_link_as_missing_role(
    tmp_path: Path,
) -> None:
    official = tmp_path / "official"
    official.mkdir()
    (official / "Complete modes details.html").write_text(
        "<html><head><title>Unrelated</title></head><body><h1>Unrelated</h1>"
        '<a href="real.html">ISODISTORT: complete modes details</a></body></html>',
        encoding="utf-8",
    )

    result = numeric._audit_directories(
        _audit_case(),
        official,
        tmp_path / "local",
        evidence_scope="test",
        provenance={"kind": "test"},
        official_html_inventory=numeric.scan_official_html(official, recursive=False),
    )

    assert result["status"] == "inconclusive"
    assert result["reason_code"] == "official_html_missing_match"
    assert result["official_html"]["resolution"]["code"] == "missing_match"
    pages = result["official_html"]["inventory"]["pages"]
    assert [page["role"] for page in pages] == ["unknown"]


def test_numeric_signature_covers_full_html_inventory_but_ignores_basename(
    tmp_path: Path,
) -> None:
    cases = [_audit_case("A", "one"), _audit_case("B", "two")]
    first_dir = tmp_path / "parent.cif" / "官网" / "Method2" / "one"
    second_dir = tmp_path / "parent.cif" / "官网" / "Method2" / "two"
    first_dir.mkdir(parents=True)
    second_dir.mkdir(parents=True)
    complete = first_dir / "browser download.html"
    unknown = first_dir / "unrelated.html"
    complete.write_text(_official_complete_modes_html(), encoding="utf-8")
    unknown.write_text("<html><title>unknown one</title></html>", encoding="utf-8")

    original = numeric._signature_state(tmp_path, cases)
    resolver_path = str(
        (
            numeric.PROJECT_ROOT
            / "tests"
            / "manual"
            / "official_html_resolver.py"
        ).resolve()
    )
    assert any(record["path"] == resolver_path for record in original["files"])

    renamed = first_dir / "用户修改后的任意名.hTm"
    complete.rename(renamed)
    after_rename = numeric._signature_state(tmp_path, cases)
    assert after_rename["signature"] == original["signature"]

    unknown.write_text("<html><title>unknown two</title></html>", encoding="utf-8")
    after_unknown_change = numeric._signature_state(tmp_path, cases)
    assert after_unknown_change["signature"] != after_rename["signature"]

    duplicate = first_dir / "duplicate.html"
    duplicate.write_text(_official_complete_modes_html(), encoding="utf-8")
    after_duplicate = numeric._signature_state(tmp_path, cases)
    assert after_duplicate["signature"] != after_unknown_change["signature"]

    duplicate.rename(second_dir / "moved duplicate.html")
    after_cross_directory_move = numeric._signature_state(tmp_path, cases)
    assert after_cross_directory_move["signature"] != after_duplicate["signature"]


def test_saved_audit_drift_invalidates_pass_and_checkpoint(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "output"
    case, _official, html_path, _captured_html = _saved_audit_tree(output_root)
    monkeypatch.setattr(numeric, "CASE_MATRIX", (case,))

    def audit_case(
        _output_root: Path,
        selected: dict[str, object],
        *,
        official_html_inventory: numeric.OfficialHtmlInventory | None,
    ) -> dict[str, object]:
        assert official_html_inventory is not None
        html_path.write_text(
            _official_complete_modes_html("changed after start signature"),
            encoding="utf-8",
        )
        return _pass_record(selected)

    monkeypatch.setattr(numeric, "audit_case", audit_case)
    report_path = tmp_path / "report.json"
    checkpoint_path = tmp_path / "checkpoint.json"

    report = numeric.run_audit(
        output_root,
        report_path,
        checkpoint_path,
        restart=True,
    )

    assert report["schema"] == numeric.SCHEMA_VERSION == 6
    assert report["status"] == "invalidated_by_concurrent_change"
    assert report["validation_status"] == "inconclusive"
    assert report["provenance_stability"]["complete_invocation_signature_stable"] is False
    assert report["summary"]["pass_case_count"] == 0
    assert report["summary"]["inconclusive_case_count"] == 1
    assert all(
        export["status"] == "inconclusive"
        for export in report["cases"][0]["exports"].values()
    )
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["schema"] == numeric.CHECKPOINT_SCHEMA == 3
    assert checkpoint["results"] == {}
    assert checkpoint["reusable"] is False
    assert checkpoint["run_state"] == "invalidated"
    assert checkpoint["signature"].startswith("invalidated:")
    with pytest.raises(RuntimeError, match="signature mismatch"):
        numeric._load_checkpoint(
            checkpoint_path,
            report["signature_state_at_start"]["signature"],
            restart=False,
        )


def test_saved_audit_aba_parses_start_capture_not_transient_html(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "output"
    case, official, html_path, captured_html = _saved_audit_tree(output_root)
    monkeypatch.setattr(numeric, "CASE_MATRIX", (case,))
    consumed: list[str] = []

    def audit_case(
        _output_root: Path,
        selected: dict[str, object],
        *,
        official_html_inventory: numeric.OfficialHtmlInventory | None,
    ) -> dict[str, object]:
        assert official_html_inventory is not None
        page = official_html_inventory.require_unique_in_directory(
            "complete_modes_details",
            official,
        )
        html_path.write_text(
            _official_complete_modes_html("transient unsigned B"),
            encoding="utf-8",
        )
        consumed.append(page.raw_html)
        html_path.write_text(captured_html, encoding="utf-8")
        return _pass_record(selected)

    monkeypatch.setattr(numeric, "audit_case", audit_case)
    report = numeric.run_audit(
        output_root,
        tmp_path / "report.json",
        tmp_path / "checkpoint.json",
        restart=True,
    )

    assert consumed == [captured_html]
    assert "transient unsigned B" not in consumed[0]
    assert report["status"] == "complete"
    assert report["provenance_stability"]["complete_invocation_signature_stable"] is True
    assert report["cases"][0]["status"] == "pass"


def test_saved_audit_resign_error_is_inconclusive_and_clears_checkpoint(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "output"
    case, _official, _html_path, _captured_html = _saved_audit_tree(output_root)
    monkeypatch.setattr(numeric, "CASE_MATRIX", (case,))
    monkeypatch.setattr(
        numeric,
        "audit_case",
        lambda _root, selected, **_kwargs: _pass_record(selected),
    )
    real_signature_state = numeric._signature_state
    calls = 0

    def signature_state(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated final re-sign failure")
        return real_signature_state(*args, **kwargs)

    monkeypatch.setattr(numeric, "_signature_state", signature_state)
    checkpoint_path = tmp_path / "checkpoint.json"
    report = numeric.run_audit(
        output_root,
        tmp_path / "report.json",
        checkpoint_path,
        restart=True,
    )

    assert report["status"] == "inconclusive_due_to_resign_error"
    assert report["validation_status"] == "inconclusive"
    assert report["signature_state_at_end"]["error_type"] == "OSError"
    assert report["cases"][0]["status"] == "inconclusive"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["results"] == {}
    assert checkpoint["reusable"] is False


def test_live_audit_drift_invalidates_pass_and_uses_start_inventory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_root = tmp_path / "output"
    case, official, _html_path, _captured_html = _saved_audit_tree(output_root)
    case["selector"] = {"route": "test"}
    monkeypatch.setattr(numeric, "LIVE_CASE_MATRIX", (case,))
    captured_inventory = numeric.scan_official_html(official, recursive=False)
    signature_calls = 0

    def signature_state(
        _output_root: Path,
        _cases: list[dict[str, object]],
        *,
        captured_html_inventories: (
            dict[str, numeric.OfficialHtmlInventory | None] | None
        ) = None,
    ) -> dict[str, object]:
        nonlocal signature_calls
        signature_calls += 1
        if captured_html_inventories is not None:
            captured_html_inventories[str(case["id"])] = captured_inventory
        suffix = "start" if signature_calls == 1 else "end"
        return {
            "signature": suffix,
            "source_signature": "source",
            "evidence_signature": f"evidence-{suffix}",
            "missing_files": [],
        }

    live_export = tmp_path / "live export"
    live_export.mkdir()

    def generate(*_args, **_kwargs) -> dict[str, object]:
        return {
            "status": "complete",
            "source_signature": "source",
            "manifest_sha256": "manifest",
            "manifest_path": str(tmp_path / "manifest.json"),
            "export_directory": str(live_export),
            "reused_immutable_generation": False,
            "mode_count": 1,
            "mode_label_diagnostics": {},
            "selected_candidate": {},
        }

    def audit_directories(
        selected: dict[str, object],
        _official: Path,
        _local: Path,
        **kwargs,
    ) -> dict[str, object]:
        assert kwargs["official_html_inventory"] is captured_inventory
        return _pass_record(selected)

    monkeypatch.setattr(numeric, "_live_signature_state", signature_state)
    monkeypatch.setattr(numeric, "_generate_live_case", generate)
    monkeypatch.setattr(numeric, "_audit_directories", audit_directories)
    checkpoint_path = tmp_path / "live-checkpoint.json"

    report = numeric.run_live_audit(
        output_root,
        tmp_path / "live",
        tmp_path / "live-report.json",
        checkpoint_path,
        restart=True,
    )

    assert report["status"] == "invalidated_by_concurrent_change"
    assert report["validation_status"] == "inconclusive"
    assert report["summary"]["inconclusive_case_count"] == 1
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["results"] == {}
    assert checkpoint["reusable"] is False
    assert checkpoint["signature"] == "invalidated:start"


def test_checkpoint_without_certified_final_signature_is_not_reusable(
    tmp_path: Path,
) -> None:
    checkpoint_path = tmp_path / "pending.json"
    numeric._atomic_json(
        checkpoint_path,
        {
            "schema": numeric.CHECKPOINT_SCHEMA,
            "signature": "same",
            "created_at": "time",
            "updated_at": "time",
            "run_state": "pending_final_signature",
            "reusable": False,
            "results": {"UNSIGNED": {"status": "pass"}},
        },
    )

    with pytest.raises(RuntimeError, match="not reusable"):
        numeric._load_checkpoint(checkpoint_path, "same", restart=False)


@pytest.mark.parametrize(
    ("report", "expected"),
    [
        (
            {
                "status": "complete",
                "validation_status": "pass",
                "provenance_stability": {
                    "complete_invocation_signature_stable": True,
                    "final_resign_succeeded": True,
                },
                "summary": {
                    "case_count": 1,
                    "pass_case_count": 1,
                    "fail_case_count": 0,
                    "inconclusive_case_count": 0,
                    "numeric_semantics_complete": True,
                    "all_numeric_checks_pass": True,
                },
                "cases": [],
            },
            0,
        ),
        (
            {
                "status": "complete",
                "validation_status": "fail",
                "provenance_stability": {
                    "complete_invocation_signature_stable": True,
                    "final_resign_succeeded": True,
                },
                "summary": {
                    "case_count": 1,
                    "pass_case_count": 0,
                    "fail_case_count": 1,
                    "inconclusive_case_count": 0,
                    "numeric_semantics_complete": True,
                    "all_numeric_checks_pass": False,
                },
                "cases": [],
            },
            1,
        ),
        (
            {
                "status": "invalidated_by_concurrent_change",
                "validation_status": "inconclusive",
                "provenance_stability": {
                    "complete_invocation_signature_stable": False,
                    "final_resign_succeeded": True,
                },
                "summary": {
                    "case_count": 1,
                    "pass_case_count": 1,
                    "fail_case_count": 0,
                    "inconclusive_case_count": 0,
                    "numeric_semantics_complete": True,
                    "all_numeric_checks_pass": True,
                },
                "cases": [],
            },
            2,
        ),
        (
            {
                "status": "complete",
                "validation_status": "inconclusive",
                "provenance_stability": {
                    "complete_invocation_signature_stable": True,
                    "final_resign_succeeded": True,
                },
                "summary": {
                    "case_count": 1,
                    "pass_case_count": 0,
                    "fail_case_count": 0,
                    "inconclusive_case_count": 1,
                    "numeric_semantics_complete": False,
                    "all_numeric_checks_pass": False,
                },
                "cases": [],
            },
            2,
        ),
        (
            {
                "status": "complete",
                "validation_status": "pass",
                "provenance_stability": {
                    "complete_invocation_signature_stable": True,
                    "final_resign_succeeded": True,
                },
                "summary": {
                    "case_count": 1,
                    "pass_case_count": True,
                    "fail_case_count": -1,
                    "inconclusive_case_count": 0,
                    "numeric_semantics_complete": True,
                    "all_numeric_checks_pass": True,
                },
                "cases": [],
            },
            2,
        ),
    ],
)
def test_main_rejects_false_success_reports(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    report: dict[str, object],
    expected: int,
) -> None:
    monkeypatch.setattr(numeric, "run_audit", lambda *_args, **_kwargs: report)
    monkeypatch.setattr(numeric, "_print_summary", lambda *_args, **_kwargs: None)

    result = numeric.main(
        [
            "--report",
            str(tmp_path / "report.json"),
            "--checkpoint",
            str(tmp_path / "checkpoint.json"),
        ]
    )

    assert result == expected
