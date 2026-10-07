"""Tests for the read-only Method output audit helpers."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import backend.api
import tests.manual.audit_method1_4310_live as method1_4310_module
import tests.manual.validate_method_outputs as validation_module
from tests.manual.audit_method1_4310_live import (
    _atomic_json,
    _classify_cif_result,
    _frozen_pair_preflight,
    _frozen_preflight_summary,
    _frozen_report_status,
    _preparse_official_mode_references,
    _resign_source_and_invalidate,
    _source_record_set_integrity,
)
from tests.manual.official_html_resolver import scan_official_html
from tests.manual.validate_method_outputs import (
    _count_local_displacive_modes,
    _live12_checkpoint_policy,
    _live12_signature,
    _resolve_official_displacive_modes,
    audit_output_root,
    bases_generate_same_lattice,
    compare_cif_alternate_settings,
    parse_candidate_header,
    run_live_method12_audit,
)


def _complete_modes_html(
    mode_count: int = 1,
    *,
    subgroup_details: str | None = None,
) -> str:
    modes = "".join(
        f"<pre>:dsp] mode-{index} normfactor = 1</pre>"
        for index in range(mode_count)
    )
    return (
        "<html><head><title>ISODISTORT: complete modes details</title></head>"
        "<body><h1>ISODISTORT: complete modes details</h1>"
        "<strong>Subgroup details</strong>"
        f"{subgroup_details or ''}"
        "<strong>Undistorted superstructure</strong>"
        "<strong>Distorted superstructure</strong>"
        f"{modes}</body></html>"
    )


def _basis(rows: tuple[tuple[int, int, int], ...]):
    return tuple(tuple(Fraction(value) for value in row) for row in rows)


def test_4310_audit_json_converts_numpy_scalars_and_rejects_unknowns(
    tmp_path: Path,
) -> None:
    output = tmp_path / "report.json"
    _atomic_json(
        output,
        {
            "passed": np.bool_(True),
            "rank": np.int64(3),
            "nonfinite": [
                float("inf"),
                float("-inf"),
                float("nan"),
                np.float64("inf"),
            ],
        },
    )
    assert json.loads(output.read_text(encoding="utf-8")) == {
        "passed": True,
        "rank": 3,
        "nonfinite": [None, None, None, None],
    }

    with pytest.raises(TypeError, match="is not JSON serializable"):
        _atomic_json(output, {"unsupported": object()})


def test_4310_end_resign_invalidates_checkpoint_records_and_blocks_aba_reuse(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    checkpoint = {
        "signature": "source-a",
        "records": {"candidate": {"status": "passed"}},
    }
    monkeypatch.setattr(
        method1_4310_module,
        "_signature_state",
        lambda: {"signature": "source-b"},
    )

    assert not _resign_source_and_invalidate(checkpoint)
    assert checkpoint["records"] == {}
    assert checkpoint["invalidated_record_count"] == 1
    assert checkpoint["source_signature"] == {
        "start": "source-a",
        "end": "source-b",
        "stable": False,
        "completion_error": None,
    }

    # Even if the source later returns to A, no result from the invalidated
    # A -> B run remains available for an ABA checkpoint reuse.
    monkeypatch.setattr(
        method1_4310_module,
        "_signature_state",
        lambda: {"signature": "source-a"},
    )
    assert _resign_source_and_invalidate(checkpoint)
    assert checkpoint["records"] == {}


def test_equivalent_basis_accepts_integer_unimodular_change() -> None:
    reference = _basis(((1, 0, 0), (0, 1, 0), (0, 0, 2)))
    equivalent = _basis(((1, 1, 0), (0, 1, 0), (0, 0, 2)))
    assert bases_generate_same_lattice(reference, equivalent)


def test_equivalent_basis_rejects_different_sublattice() -> None:
    reference = _basis(((1, 0, 0), (0, 1, 0), (0, 0, 2)))
    different = _basis(((2, 0, 0), (0, 1, 0), (0, 0, 2)))
    assert not bases_generate_same_lattice(reference, different)


def _write_setting_cif(
    path: Path,
    *,
    basis: str = "(1,0,0),(0,1,0),(0,0,1)",
    origin: str = "(0,0,0)",
    cell: tuple[str, str, str, str, str, str] = (
        "1.00000",
        "1.00000",
        "1.00000",
        "90.00000",
        "90.00000",
        "90.00000",
    ),
    species: str = "C",
    coordinates: tuple[str, str, str] = ("0.10000", "0.20000", "0.30000"),
    occupancy: str = "1.00000",
    atom_rows: tuple[tuple[str, str, str, str, str, str], ...] | None = None,
    symmetry_operations: tuple[str, ...] | None = ("x,y,z",),
    magnetic_moments: tuple[tuple[str, str, str], ...] | None = None,
    space_group_symbol: str = "P 1",
    space_group_number: str | None = "1",
) -> None:
    a, b, c, alpha, beta, gamma = cell
    x, y, z = coordinates
    number_line = (
        f"_symmetry_Int_Tables_number {space_group_number}\n"
        if space_group_number is not None
        else ""
    )
    if atom_rows is None:
        atom_rows = ((f"{species}1", species, x, y, z, occupancy),)
    if magnetic_moments is not None and len(magnetic_moments) != len(atom_rows):
        raise ValueError("magnetic moment rows must align with atom rows")
    symmetry_loop = ""
    if symmetry_operations is not None:
        symmetry_loop = (
            "loop_\n_symmetry_equiv_pos_as_xyz\n"
            + "\n".join(f"'{operation}'" for operation in symmetry_operations)
            + "\n"
        )
    moment_columns = ""
    if magnetic_moments is not None:
        moment_columns = (
            "_atom_site_moment.crystalaxis_x\n"
            "_atom_site_moment.crystalaxis_y\n"
            "_atom_site_moment.crystalaxis_z\n"
        )
    atom_lines = []
    for index, row in enumerate(atom_rows):
        fields = list(row)
        if magnetic_moments is not None:
            fields.extend(magnetic_moments[index])
        atom_lines.append(" ".join(fields))
    path.write_text(
        f"""# IR: GM1+
# P1 (a) 1 P1, basis={{{basis}}}, origin={origin}, s=1, i=1, k-active= (0,0,0)
data_test
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha {alpha}
_cell_angle_beta {beta}
_cell_angle_gamma {gamma}
_symmetry_space_group_name_H-M '{space_group_symbol}'
{number_line}{symmetry_loop}loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
{moment_columns}{chr(10).join(atom_lines)}
""",
        encoding="utf-8",
    )


def test_precision_cif_comparison_applies_exact_basis_origin_and_rounding(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        basis="(1,1,0),(0,1,0),(0,0,1)",
        origin="(1/2,0,0)",
        cell=(
            "1.41421",
            "1.00000",
            "1.00000",
            "90.00000",
            "90.00000",
            "45.00000",
        ),
        coordinates=("0.10000", "0.20000", "0.30000"),
    )
    _write_setting_cif(
        official,
        coordinates=("0.59999", "0.30001", "0.30000"),
    )

    comparison = compare_cif_alternate_settings(local, official)

    assert comparison.equivalent
    assert comparison.conclusive
    assert comparison.details["transform"] == {
        "local_to_official_matrix": [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
        "origin_shift_in_official_cell": ["1/2", 0, 0],
        "determinant": 1,
        "integer": True,
        "unimodular": True,
        "unique": True,
        "coordinate_convention": "x_official = x_local @ U + q (mod 1)",
    }
    assert comparison.details["precision"][
        "coordinate_bound_by_official_component"
    ] == pytest.approx([1e-5, 1.5e-5, 1e-5])
    assert comparison.details["lattice"]["metric_equal_at_output_precision"]
    assert comparison.details["sites"]["coordinates_equal_at_output_precision"]


@pytest.mark.parametrize(
    ("change", "expected_issue"),
    [
        ({"coordinates": ("0.10002", "0.20000", "0.30000")}, "coordinate-mismatch-outside-output-precision"),
        ({"species": "O"}, "species-count-mismatch"),
        ({"occupancy": "0.99998"}, "occupancy-mismatch-outside-output-precision"),
    ],
)
def test_precision_cif_comparison_rejects_semantic_changes(
    tmp_path: Path,
    change: dict[str, object],
    expected_issue: str,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local, **change)
    _write_setting_cif(official)

    comparison = compare_cif_alternate_settings(local, official)

    assert not comparison.equivalent
    assert expected_issue in comparison.details["issues"]
    if expected_issue == "coordinate-mismatch-outside-output-precision":
        assert "occupancy-mismatch-outside-output-precision" not in comparison.details["issues"]


def test_precision_cif_comparison_accepts_occupancy_rounding_boundary(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local, occupancy="0.99999")
    _write_setting_cif(official, occupancy="1.00000")

    comparison = compare_cif_alternate_settings(local, official)

    assert comparison.equivalent
    assert not comparison.details["sites"]["composition_equal"]
    assert comparison.details["sites"]["occupancies_equal_at_output_precision"]


def test_precision_cif_comparison_rejects_gram_difference_beyond_precision(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        basis="(1,1,0),(0,1,0),(0,0,1)",
        cell=(
            "1.41430",
            "1.00000",
            "1.00000",
            "90.00000",
            "90.00000",
            "45.00000",
        ),
    )
    _write_setting_cif(
        official,
        coordinates=("0.10000", "0.30000", "0.30000"),
    )

    comparison = compare_cif_alternate_settings(local, official)

    assert not comparison.equivalent
    assert "lattice-metric-outside-output-precision" in comparison.details["issues"]
    assert comparison.details["lattice"]["failed_components"]


def test_precision_cif_comparison_requires_one_joint_metric_witness(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        basis="(1,1,0),(0,1,0),(0,0,1)",
        cell=(
            "1.58113883",
            "1.00000000",
            "1.00000000",
            "90.00000000",
            "90.00000000",
            "18.43494882",
        ),
        coordinates=("0.10000000", "0.20000000", "0.30000000"),
    )
    # Every transformed Gram-component interval overlaps independently because
    # the printed official b length is intentionally coarse.  The same b cannot
    # satisfy local G22=1 and local G12=1.5, so there is no joint cell witness.
    _write_setting_cif(
        official,
        cell=(
            "1.00000000",
            "1",
            "1.00000000",
            "90.00000000",
            "90.00000000",
            "90.00000000",
        ),
        coordinates=("0.10000000", "0.30000000", "0.30000000"),
    )

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert not comparison.details["lattice"]["failed_components"]
    assert not comparison.details["lattice"]["joint_feasibility_witness"]["found"]
    assert not comparison.equivalent
    assert not comparison.conclusive
    assert "lattice-output-precision-joint-feasibility-not-proven" in (
        comparison.details["inconclusive_reasons"]
    )


def test_precision_cif_comparison_rejects_non_unimodular_basis(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local, basis="(2,0,0),(0,1,0),(0,0,1)")
    _write_setting_cif(official)

    comparison = compare_cif_alternate_settings(local, official)

    assert not comparison.equivalent
    assert comparison.details["transform"]["determinant"] == 2
    assert comparison.details["issues"] == ["non-unimodular-basis-transform"]


def test_precision_cif_comparison_marks_multiple_setting_headers_inconclusive(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local)
    _write_setting_cif(official)
    text = local.read_text(encoding="utf-8")
    duplicate = next(
        line
        for line in text.splitlines()
        if line.startswith("#") and "basis={" in line
    )
    local.write_text(f"{duplicate}\n{text}", encoding="utf-8")

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert not comparison.equivalent
    assert not comparison.conclusive
    assert not comparison.details["transform"]["unique"]
    assert comparison.details["inconclusive_reasons"] == [
        "setting-transform-not-uniquely-declared"
    ]


def test_precision_cif_comparison_conjugates_full_seitz_set_exactly(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        basis="(1,1,0),(0,1,0),(0,0,1)",
        origin="(1/2,1/3,1/4)",
        cell=(
            "1.41421356",
            "1.00000000",
            "1.00000000",
            "90.00000000",
            "90.00000000",
            "45.00000000",
        ),
        coordinates=("0.10000000", "0.20000000", "0.30000000"),
        symmetry_operations=(
            "x,y,z",
            "x,y,z",  # duplicate declarations are removed before set comparison
            "y+1/4,x+3/4,z+1/2",
            "y+1/4,x+3/4,z+1/2",
        ),
    )
    _write_setting_cif(
        official,
        cell=(
            "1.00000000",
            "1.00000000",
            "1.00000000",
            "90.00000000",
            "90.00000000",
            "90.00000000",
        ),
        coordinates=("0.60000000", "0.63333333", "0.55000000"),
        symmetry_operations=(
            "x,y,z",
            "-x+y+11/12,y+1,z+1/2",
        ),
    )

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert comparison.equivalent
    operation = comparison.details["space_group"]["operation_comparison"]
    assert operation["status"] == "equivalent"
    assert operation["local_declared_operation_count"] == 4
    assert operation["local_unique_operation_count"] == 2
    assert operation["official_declared_operation_count"] == 2
    assert operation["official_unique_operation_count"] == 2
    witness = comparison.details["lattice"]["joint_feasibility_witness"]
    assert witness["all_six_gram_components_jointly_validated"]
    assert all(witness["local_components_inside_printed_box_intervals"])
    assert all(witness["official_components_inside_printed_box_intervals"])


def test_precision_cif_comparison_rejects_tampered_seitz_translation(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    common_local = {
        "basis": "(1,1,0),(0,1,0),(0,0,1)",
        "origin": "(1/2,1/3,1/4)",
        "cell": (
            "1.41421356",
            "1.00000000",
            "1.00000000",
            "90.00000000",
            "90.00000000",
            "45.00000000",
        ),
        "coordinates": ("0.10000000", "0.20000000", "0.30000000"),
        "symmetry_operations": ("x,y,z", "y+1/4,x+3/4,z+1/2"),
    }
    _write_setting_cif(local, **common_local)
    _write_setting_cif(
        official,
        cell=(
            "1.00000000",
            "1.00000000",
            "1.00000000",
            "90.00000000",
            "90.00000000",
            "90.00000000",
        ),
        coordinates=("0.60000000", "0.63333333", "0.55000000"),
        symmetry_operations=("x,y,z", "-x+y+1/12,y,z+1/2"),
    )

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert not comparison.equivalent
    assert comparison.conclusive
    assert "declared-space-group-operation-mismatch" in comparison.details["issues"]
    operation = comparison.details["space_group"]["operation_comparison"]
    assert operation["missing_operation_count"] == 1
    assert operation["extra_operation_count"] == 1


def test_precision_cif_comparison_keeps_inferred_symmetry_diagnostic(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local)
    _write_setting_cif(official)
    direct_details = {
        "declared_space_group_equal": True,
        "inferred_space_group_equal": False,
        "space_group": {
            "local_inferred": {"number": 123, "symbol": "P4/mmm"},
            "reference_inferred": {"number": 139, "symbol": "I4/mmm"},
        },
        "issues": [
            "spglib-inferred space group differs",
            "declared space group does not match spglib inference",
        ],
    }

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details=direct_details,
    )

    assert comparison.equivalent
    assert comparison.details["space_group"]["inference_is_diagnostic_only"]
    assert comparison.details["diagnostics"] == [
        "spglib-inferred-space-group-differs",
        "declared-space-group-differs-from-spglib-inference",
    ]


def test_precision_cif_comparison_prefers_declared_it_number_over_symbol(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local, space_group_symbol="P 1")
    _write_setting_cif(official, space_group_symbol="P1 alternate spelling")

    comparison = compare_cif_alternate_settings(local, official)

    assert comparison.equivalent
    assert comparison.details["space_group"]["declared_comparison_status"] == (
        "same-it-number"
    )
    assert "declared-symbol-differs-with-same-it-number" in comparison.details[
        "diagnostics"
    ]


def test_precision_cif_comparison_rejects_different_declared_it_numbers(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local, space_group_number="1")
    _write_setting_cif(official, space_group_number="2")

    comparison = compare_cif_alternate_settings(local, official)

    assert not comparison.equivalent
    assert "declared-space-group-it-number-mismatch" in comparison.details["issues"]


def test_precision_cif_comparison_marks_unresolved_symbols_inconclusive(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        space_group_symbol="P 1",
        space_group_number=None,
        symmetry_operations=None,
    )
    _write_setting_cif(
        official,
        space_group_symbol="P1 alternate setting",
        space_group_number=None,
        symmetry_operations=None,
    )

    comparison = compare_cif_alternate_settings(local, official, direct_details={})

    assert not comparison.equivalent
    assert not comparison.conclusive
    assert comparison.details["space_group"]["declared_comparison_status"] == (
        "inconclusive-without-comparable-it-numbers"
    )
    assert "declared-space-group-comparison-inconclusive" in comparison.details[
        "diagnostics"
    ]
    assert comparison.details["inconclusive_reasons"] == [
        "declared-space-group-operation-equivalence-not-proven"
    ]


def test_precision_cif_comparison_allows_disordered_site_total_occupancy_one(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    rows = (
        ("C1", "C", "0.10000", "0.20000", "0.30000", "0.40000"),
        ("O1", "O", "0.10000", "0.20000", "0.30000", "0.60000"),
    )
    _write_setting_cif(local, atom_rows=rows)
    _write_setting_cif(official, atom_rows=rows)

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert comparison.equivalent
    assert comparison.details["sites"]["local_atom_count"] == 2
    assert comparison.details["sites"]["local_composition"] == {"C": 0.4, "O": 0.6}


def test_precision_cif_comparison_rejects_raw_occupancy_above_one(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local, occupancy="1.10000")
    _write_setting_cif(official)

    with pytest.raises(ValueError, match="invalid occupancy"):
        compare_cif_alternate_settings(local, official, direct_details={})


def test_precision_cif_comparison_rejects_total_site_occupancy_above_one(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    rows = (
        ("C1", "C", "0.10000", "0.20000", "0.30000", "0.60000"),
        ("O1", "O", "0.10000", "0.20000", "0.30000", "0.50000"),
    )
    _write_setting_cif(local, atom_rows=rows)
    _write_setting_cif(official)

    with pytest.raises(ValueError, match="exceed total occupancy 1"):
        compare_cif_alternate_settings(local, official, direct_details={})


def test_precision_cif_comparison_requires_one_aligned_raw_atom_loop(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(local)
    _write_setting_cif(official)
    text = local.read_text(encoding="utf-8")
    normal_loop = """loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 0.10000 0.20000 0.30000 1.00000
"""
    split_loops = """loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_z
_atom_site_occupancy
C1 C 0.10000 0.30000 1.00000
loop_
_atom_site_fract_y
0.20000
"""
    local.write_text(text.replace(normal_loop, split_loops), encoding="utf-8")

    with pytest.raises(ValueError, match="exactly one atom-site loop"):
        compare_cif_alternate_settings(local, official, direct_details={})


def test_precision_cif_comparison_propagates_precision_per_atom_row(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        atom_rows=(
            ("C1", "C", "0.10002", "0.20000", "0.30000", "1.00000"),
            ("O1", "O", "0.5", "0.6", "0.7", "1.0"),
        ),
    )
    _write_setting_cif(
        official,
        atom_rows=(
            ("C1", "C", "0.10000", "0.20000", "0.30000", "1.00000"),
            ("O1", "O", "0.5", "0.6", "0.7", "1.0"),
        ),
    )

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert not comparison.equivalent
    assert comparison.conclusive
    assert "coordinate-mismatch-outside-output-precision" in comparison.details["issues"]
    assert comparison.details["precision"][
        "local_coordinate_quantum_ranges_by_component"
    ][0] == pytest.approx([1e-5, 0.1])


def test_precision_cif_comparison_marks_magnetic_cif_inconclusive(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    _write_setting_cif(
        local,
        magnetic_moments=(("1.0", "0.0", "0.0"),),
    )
    _write_setting_cif(
        official,
        magnetic_moments=(("1.0", "0.0", "0.0"),),
    )

    comparison = compare_cif_alternate_settings(
        local,
        official,
        direct_details={},
    )

    assert not comparison.equivalent
    assert not comparison.conclusive
    assert "magnetic-moment-setting-transform-not-implemented" in (
        comparison.details["inconclusive_reasons"]
    )


def test_direct_match_cannot_override_precision_aware_failure() -> None:
    setting = SimpleNamespace(equivalent=False, conclusive=True)

    assert _classify_cif_result(
        setting,
        direct_structure_equal=True,
        equivalent_setting_label="precision-equivalent-setting",
    ) == "failed"


def test_frozen_report_integrity_closes_count_uniqueness_and_drift_gates() -> None:
    source = {
        "summary": {"completed_records": 2},
        "identity": {
            "paired_candidates": 2,
            "official_candidates": 2,
            "live_candidates": 2,
            "duplicate_live": [],
            "duplicate_official": [],
            "extra_live": [],
            "missing_live": [],
        },
        "signature_unchanged_at_completion": True,
    }
    duplicate_records = [
        {"identity": "same", "position": 1, "errors": []},
        {"identity": "same", "position": 1, "errors": []},
    ]

    failures, expected_count, identities, positions = _source_record_set_integrity(
        source,
        duplicate_records,
    )
    codes = {failure["code"] for failure in failures}

    assert expected_count == 2
    assert identities == ["same", "same"]
    assert positions == [1, 1]
    assert "source-identities-not-unique-and-complete" in codes
    assert "source-positions-not-unique-and-complete" in codes
    accepted_records = [
        {"status": "direct"},
        {"status": "precision-equivalent-setting"},
    ]
    assert _frozen_report_status(accepted_records, failures) == "complete-failed"
    assert _frozen_report_status(
        accepted_records,
        [{"code": "frozen-cif-corpus-changed-during-reanalysis"}],
    ) == "complete-failed"
    assert _frozen_report_status(accepted_records, []) == "complete-passed"


def test_frozen_pair_preflight_records_raw_loops_operations_and_unique_transform(
    tmp_path: Path,
) -> None:
    local = tmp_path / "local.cif"
    official = tmp_path / "official.cif"
    rows = (
        ("C1", "C", "0.10000", "0.20000", "0.30000", "0.40000"),
        ("O1", "O", "0.10000", "0.20000", "0.30000", "0.60000"),
    )
    _write_setting_cif(local, atom_rows=rows)
    _write_setting_cif(official, atom_rows=rows)

    preflight = _frozen_pair_preflight(local, official)
    summary = _frozen_preflight_summary([{"preflight": preflight}])

    assert preflight["transform"] == {
        "unique": True,
        "matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "origin_shift": [0, 0, 0],
        "determinant": 1,
        "integer": True,
        "unimodular": True,
    }
    assert preflight["local"]["atom_site_raw_row_count"] == 2
    assert preflight["local"]["space_group_operation_tags"] == [
        "_symmetry_equiv_pos_as_xyz"
    ]
    assert preflight["local"]["space_group_declared_operation_count"] == 1
    assert preflight["local"]["space_group_unique_operation_count"] == 1
    assert summary["record_pairs_parsed"] == 1
    assert summary["parse_error_count"] == 0
    assert summary["unique_transform_count"] == 1
    assert summary["non_unimodular_transform_count"] == 0


def test_candidate_identity_reduces_origin_modulo_integer_translation(tmp_path: Path) -> None:
    template = """# IR: GM1+
# P1 (a) 1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=ORIGIN, s=1, i=1, k-active= (0,0,0)
data_test
"""
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    first.write_text(template.replace("ORIGIN", "(5/4,-3/2,2)"), encoding="utf-8")
    second.write_text(template.replace("ORIGIN", "(1/4,1/2,0)"), encoding="utf-8")

    assert parse_candidate_header(first).identity == parse_candidate_header(second).identity


def test_local_mode_count_excludes_prose_headings() -> None:
    text = """Mode definitions (every atom)
Mode vectors are given in unitless coordinates
Mode LD1[0,0,1/6]__a__A2u(a)  pretty label
Mode GM5-[0,0,0]__e__Eu(b)  another label
Mode amplitudes (As = supercell-normalized)
"""
    assert _count_local_displacive_modes(text) == 2


def test_local_mode_count_excludes_preceding_strain_modes() -> None:
    text = """Strain mode definitions
Mode GM1+strain_1(a)
Mode GM1+strain_2(a)
Displacive mode definitions (every atom in the unit cell)
Mode GM1+__e__A1(a)  displacive label
Mode amplitudes (As = supercell-normalized)
"""

    assert _count_local_displacive_modes(text) == 1


def test_official_mode_reference_uses_content_not_basename(tmp_path: Path) -> None:
    page = tmp_path / "user renamed this page.HTM"
    page.write_text(_complete_modes_html(2), encoding="utf-8")

    reference = _resolve_official_displacive_modes(tmp_path)

    assert reference["count"] == 2
    provenance = reference["page_provenance"]
    assert provenance["selected_page"]["path"] == str(page.resolve())
    assert provenance["selected_page"]["role"] == "complete_modes_details"
    assert provenance["selected_page"]["content_sha256"]


def test_official_mode_reference_parses_only_the_captured_recursive_page(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    page = candidate / "renamed.htm"
    page.write_text(_complete_modes_html(2), encoding="utf-8")
    inventory = scan_official_html(tmp_path, recursive=True)

    # A post-capture replacement must neither enter this run nor trigger a
    # second read hidden inside the candidate parser.
    page.write_text(_complete_modes_html(9), encoding="utf-8")
    reference = _resolve_official_displacive_modes(
        candidate,
        inventory=inventory,
    )

    assert reference["count"] == 2
    assert reference["page_provenance"]["inventory_signature"] == inventory.signature
    page.write_text(_complete_modes_html(2), encoding="utf-8")
    assert scan_official_html(tmp_path, recursive=True).signature == inventory.signature
    # The A -> B -> A interval is safe because the candidate used captured A,
    # never transient B; the end re-sign therefore agrees with what was parsed.


@pytest.mark.parametrize(
    "subgroup_details",
    [
        (
            "2 P-1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
        ),
        (
            "1 P1, basis={(2,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
        ),
        (
            "1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(1/2,0,0), s=1, i=1, k-active= (0,0,0)"
        ),
        (
            "1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=2, i=1, k-active= (0,0,0)"
        ),
        (
            "1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=2, k-active= (0,0,0)"
        ),
        (
            "1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1, k-active= (1/2,0,0)"
        ),
    ],
)
def test_4310_reference_rejects_same_count_page_with_wrong_subgroup_identity(
    tmp_path: Path,
    subgroup_details: str,
) -> None:
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    _write_setting_cif(candidate / "subgroup.cif")
    (candidate / "details.html").write_text(
        _complete_modes_html(
            1,
            subgroup_details=subgroup_details,
        ),
        encoding="utf-8",
    )
    inventory = scan_official_html(tmp_path, recursive=True)

    references, errors = _preparse_official_mode_references(
        {"candidate": candidate},
        inventory,
    )

    assert references == {}
    assert len(errors) == 1
    assert errors[0]["candidate"] == "candidate"
    assert "official-subgroup-identity-mismatch" in errors[0]["message"]


@pytest.mark.parametrize("shape", ["missing", "fake", "duplicate"])
def test_official_mode_reference_fails_closed(
    tmp_path: Path,
    shape: str,
) -> None:
    if shape == "fake":
        (tmp_path / "complete modes details.html").write_text(
            "<html><title>unrelated</title><body>not an official result</body></html>",
            encoding="utf-8",
        )
    elif shape == "duplicate":
        for name in ("one.html", "two.htm"):
            (tmp_path / name).write_text(_complete_modes_html(), encoding="utf-8")

    with pytest.raises(ValueError, match=r"missing_match|ambiguous_match"):
        _resolve_official_displacive_modes(tmp_path)


def test_method1_reference_preflight_captures_errors_and_page_provenance(
    tmp_path: Path,
) -> None:
    valid = tmp_path / "valid"
    invalid = tmp_path / "invalid"
    valid.mkdir()
    invalid.mkdir()
    _write_setting_cif(valid / "subgroup.cif")
    _write_setting_cif(invalid / "subgroup.cif")
    page = valid / "arbitrary.htm"
    page.write_text(
        _complete_modes_html(
            3,
            subgroup_details=(
                "1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
                "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
            ),
        ),
        encoding="utf-8",
    )
    (invalid / "misleading complete modes details.html").write_text(
        "<html><title>not official</title></html>",
        encoding="utf-8",
    )

    inventory = scan_official_html(tmp_path, recursive=True)
    references, errors = _preparse_official_mode_references(
        {"valid identity": valid, "invalid identity": invalid},
        inventory,
    )

    assert references["valid identity"]["count"] == 3
    assert references["valid identity"]["page_provenance"]["selected_page"][
        "path"
    ] == str(page.resolve())
    assert [error["candidate"] for error in errors] == ["invalid identity"]
    assert errors[0]["resolution_code"] == "missing_match"
    assert errors[0]["page_provenance"]["pages"][0]["role"] == "unknown"


def test_saved_output_audit_records_incomplete_official_mode_reference(
    tmp_path: Path,
) -> None:
    official = tmp_path / "official"
    local = tmp_path / "local"
    official.mkdir()
    local.mkdir()
    (official / "convincing filename complete modes details.html").write_text(
        "<html><title>unrelated</title></html>",
        encoding="utf-8",
    )
    (local / "Complete modes details.txt").write_text(
        "Displacive mode definitions\nMode GM1+__a(a)\nMode amplitudes\n",
        encoding="utf-8",
    )
    summary = validation_module.AuditSummary(output_root=str(tmp_path))

    validation_module._compare_mode_counts(
        summary,
        "case",
        "candidate",
        official,
        local,
    )

    assert summary.failures == 1
    assert summary.mode_files_compared == 0
    assert summary.issues[0].code == "official-mode-reference-incomplete"
    reference = summary.official_mode_references[0]
    assert reference["status"] == "incomplete_reference"
    assert reference["error"]["resolution_code"] == "missing_match"
    assert reference["error"]["page_provenance"]["pages"][0]["role"] == "unknown"


def test_saved_method12_audit_delegates_method3_affine_rows(tmp_path: Path) -> None:
    for side in ("官网", "现有网页版交互"):
        candidate = tmp_path / "parent.cif" / side / "Method3" / "affine candidate"
        candidate.mkdir(parents=True)
        # An affine/coupled Method 3 CIF need not carry a single-IR header.
        (candidate / "subgroup.cif").write_text(
            "# Subgroup:  1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1\n",
            encoding="utf-8",
        )

    summary = audit_output_root(tmp_path)

    assert summary.failures == 0
    assert summary.groups == {}


def test_live12_signature_covers_scientific_inputs_but_not_generated_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workspace = tmp_path / "workspace"
    package_root = workspace / "ISODISTORT"
    config_dir = package_root / "resources" / "config"
    runtime_dir = package_root / "runtime"
    backend_dir = package_root / "backend"
    input_dir = workspace / "experiment_data"
    output_root = workspace / "output_compare"
    for directory in (config_dir, runtime_dir, backend_dir, input_dir):
        directory.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "settings.yaml"
    config_path.write_text(
        """isobyu:
  bin_dir: ../../runtime
  data_dir: ../../runtime
  iso_bin: iso
  smodes_bin: smodes
""",
        encoding="utf-8",
    )
    (backend_dir / "engine.py").write_text("VALUE = 1\n", encoding="utf-8")
    for parent in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        (input_dir / parent).write_text(f"parent {parent}\n", encoding="utf-8")
    iso_path = runtime_dir / "iso"
    smodes_path = runtime_dir / "smodes"
    data_path = runtime_dir / "data_isotropy.txt"
    iso_path.write_bytes(b"iso-v1")
    smodes_path.write_bytes(b"smodes-v1")
    data_path.write_bytes(b"data-v1")
    reference_dir = output_root / "EuAl4 Parent.cif" / "官网" / "Method1" / "candidate"
    reference_dir.mkdir(parents=True)
    (reference_dir / "subgroup.cif").write_text("reference-v1\n", encoding="utf-8")
    details_path = reference_dir / "ISODISTORT complete modes details.html"
    details_path.write_text(_complete_modes_html(), encoding="utf-8")

    runtime = {"python": {"version": "runtime-v1"}}
    monkeypatch.setattr(validation_module, "_live12_runtime_versions", lambda: dict(runtime))

    def signature() -> str:
        return _live12_signature(
            output_root,
            workspace=workspace,
            package_root=package_root,
        )

    initial = signature()
    renamed_details = reference_dir / "manually renamed evidence.HTM"
    details_path.rename(renamed_details)
    assert signature() == initial
    details_path = renamed_details
    generated_cache = runtime_dir / "i123.iso"
    generated_cache.write_bytes(b"generated-v1")
    assert signature() == initial
    generated_cache.write_bytes(b"generated-v2")
    assert signature() == initial
    unrelated_method3 = output_root / "EuAl4 Parent.cif" / "官网" / "Method3" / "candidate"
    unrelated_method3.mkdir(parents=True)
    (unrelated_method3 / "subgroup.cif").write_text("method3-only\n", encoding="utf-8")
    assert signature() == initial

    signatures = [initial]
    config_path.write_text(config_path.read_text(encoding="utf-8") + "# changed\n", encoding="utf-8")
    signatures.append(signature())
    iso_path.write_bytes(b"iso-v2")
    signatures.append(signature())
    smodes_path.write_bytes(b"smodes-v2")
    signatures.append(signature())
    data_path.write_bytes(b"data-v2")
    signatures.append(signature())
    runtime["python"] = {"version": "runtime-v2"}
    signatures.append(signature())
    details_path.write_text(_complete_modes_html() + "<!-- changed -->", encoding="utf-8")
    signatures.append(signature())
    assert len(set(signatures)) == len(signatures)

    content_baseline = signature()
    unknown = reference_dir / "unrelated.html"
    unknown.write_text("<html><title>unknown</title></html>", encoding="utf-8")
    with_unknown = signature()
    assert with_unknown != content_baseline
    duplicate = reference_dir / "second complete page.htm"
    duplicate.write_text(_complete_modes_html(), encoding="utf-8")
    with_duplicate = signature()
    assert with_duplicate != with_unknown
    moved_dir = output_root / "EuAl4 Parent.cif" / "官网" / "Method1" / "other candidate"
    moved_dir.mkdir()
    moved = moved_dir / details_path.name
    details_path.rename(moved)
    assert signature() != with_duplicate
    assert not _live12_checkpoint_policy()["generated_wsl_i_star_iso_cache"][
        "included_in_global_signature"
    ]


class _Live12Subgroup:
    irrep_label = "GM1+"
    origin = (0, 0, 0)
    index = 0
    basis_scale = 1

    @classmethod
    def opd_line_body(cls) -> str:
        return (
            f"P1 (a) 1 P1, basis={{({cls.basis_scale},0,0),(0,1,0),(0,0,1)}}, "
            "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
        )


def _write_live12_reference(output_root: Path, case_name: str, method: str) -> None:
    folder = output_root / case_name / "官网" / method / "GM1+ P1"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "subgroup.cif").write_text(
        """# IR: GM1+
# P1 (a) 1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=(0,0,0), s=1, i=1, k-active= (0,0,0)
data_test
""",
        encoding="utf-8",
    )
    (folder / "ISODISTORT complete modes details.html").write_text(
        _complete_modes_html(),
        encoding="utf-8",
    )


def test_live12_checkpoint_marks_interruption_and_safely_reuses_identity_counts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workspace = tmp_path / "workspace"
    output_root = workspace / "output_compare"
    input_root = workspace / "experiment_data"
    input_root.mkdir(parents=True)
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        (input_root / case_name).write_text("parent\n", encoding="utf-8")
        for method in ("Method1", "Method2"):
            _write_live12_reference(output_root, case_name, method)

    state = {"mode_calls": 0, "fail_after": 2}
    monkeypatch.setattr(_Live12Subgroup, "basis_scale", 1)

    class FakeIsoDistort:
        def __init__(self) -> None:
            self.mode_displacements: list[object] = []

        def load_structure(self, _path: Path) -> None:
            return None

        def search_method_1(self, **_kwargs: object) -> list[SimpleNamespace]:
            return [SimpleNamespace(subgroup=_Live12Subgroup())]

        def list_subgroups_at_kpoint(self, *_args: object, **_kwargs: object) -> list[_Live12Subgroup]:
            return [_Live12Subgroup()]

        def search_method_2(self, **_kwargs: object) -> None:
            if state["fail_after"] is not None and state["mode_calls"] >= state["fail_after"]:
                raise RuntimeError("simulated interruption")
            state["mode_calls"] += 1
            self.mode_displacements = [object()]

    monkeypatch.setattr(validation_module, "WORKSPACE", workspace)
    monkeypatch.setattr(backend.api, "IsoDistort", FakeIsoDistort)
    checkpoint_path = workspace / "output" / "validation" / "checkpoint.json"
    report_path = workspace / "output" / "validation" / "report.json"

    with pytest.raises(RuntimeError, match="simulated interruption"):
        run_live_method12_audit(output_root, checkpoint_path, report_path)
    interrupted_checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    interrupted_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert interrupted_checkpoint["status"] == "incomplete_error"
    assert not interrupted_checkpoint["execution_complete"]
    assert interrupted_report["status"] == "incomplete_error"
    assert interrupted_report["progress"]["completed_groups"] == 2
    assert interrupted_report["last_invocation_elapsed_seconds"] >= 0

    state["fail_after"] = None
    completed = run_live_method12_audit(output_root, checkpoint_path, report_path)
    assert completed["status"] == "complete_passed"
    assert completed["execution_complete"]
    assert completed["validation_complete"]
    assert completed["passed"]
    assert completed["resumed_from_checkpoint"]
    assert sum(report["checkpoint_counts_reused"] for report in completed["reports"]) == 2
    assert sum(report["mode_counts_computed"] for report in completed["reports"]) == 2
    assert state["mode_calls"] == 4
    saved_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert saved_report["status"] == "complete_passed"
    assert saved_report["progress"]["completed_groups"] == 4

    # The basis-independent official identity stays the same, but a changed
    # live basis must invalidate each saved mode count rather than mixing runs.
    _Live12Subgroup.basis_scale = 2
    changed_candidate = run_live_method12_audit(output_root, checkpoint_path, report_path)
    assert sum(report["checkpoint_counts_reused"] for report in changed_candidate["reports"]) == 0
    assert sum(report["mode_counts_computed"] for report in changed_candidate["reports"]) == 4
    assert state["mode_calls"] == 8


def test_live12_end_source_drift_invalidates_all_checkpoint_results(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workspace = tmp_path / "workspace"
    output_root = workspace / "output_compare"
    input_root = workspace / "experiment_data"
    input_root.mkdir(parents=True)
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        (input_root / case_name).write_text("parent\n", encoding="utf-8")
        for method in ("Method1", "Method2"):
            _write_live12_reference(output_root, case_name, method)

    changed_page = (
        output_root
        / "EuAl4 Parent.cif"
        / "官网"
        / "Method1"
        / "GM1+ P1"
        / "ISODISTORT complete modes details.html"
    )
    original_page = changed_page.read_text(encoding="utf-8")
    state = {"mutate": True, "mode_calls": 0}

    class DriftIsoDistort:
        def __init__(self) -> None:
            self.mode_displacements: list[object] = []

        def load_structure(self, _path: Path) -> None:
            return None

        def search_method_1(self, **_kwargs: object) -> list[SimpleNamespace]:
            return [SimpleNamespace(subgroup=_Live12Subgroup())]

        def list_subgroups_at_kpoint(
            self,
            *_args: object,
            **_kwargs: object,
        ) -> list[_Live12Subgroup]:
            return [_Live12Subgroup()]

        def search_method_2(self, **_kwargs: object) -> None:
            state["mode_calls"] += 1
            self.mode_displacements = [object()]
            if state["mutate"]:
                changed_page.write_text(_complete_modes_html(2), encoding="utf-8")
                state["mutate"] = False

    monkeypatch.setattr(_Live12Subgroup, "basis_scale", 1)
    monkeypatch.setattr(validation_module, "WORKSPACE", workspace)
    monkeypatch.setattr(backend.api, "IsoDistort", DriftIsoDistort)
    checkpoint_path = workspace / "output" / "validation" / "checkpoint.json"
    report_path = workspace / "output" / "validation" / "report.json"

    invalidated = run_live_method12_audit(
        output_root,
        checkpoint_path,
        report_path,
    )

    assert invalidated["status"] == "invalidated_source_drift"
    assert invalidated["source_drift_detected"]
    assert not invalidated["execution_complete"]
    assert not invalidated["validation_complete"]
    assert not invalidated["passed"]
    assert invalidated["reports"] == []
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert checkpoint["counts"] == {}
    assert checkpoint["reports"] == {}
    assert checkpoint["invalidated_result_counts"] == {
        "candidate_counts": 4,
        "group_reports": 4,
    }

    changed_page.write_text(original_page, encoding="utf-8")
    rerun = run_live_method12_audit(output_root, checkpoint_path, report_path)
    assert rerun["status"] == "complete_passed"
    assert not rerun["resumed_from_checkpoint"]
    assert sum(report["checkpoint_counts_reused"] for report in rerun["reports"]) == 0
    assert sum(report["mode_counts_computed"] for report in rerun["reports"]) == 4
    assert state["mode_calls"] == 8
