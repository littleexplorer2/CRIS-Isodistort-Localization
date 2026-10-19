"""Tests for test_audit_method1_4310_live; migrated without changing scientific expectations."""
from __future__ import annotations

from pathlib import Path
from tests.audit_method1_4310_live import _atomic_json
from tests.audit_method1_4310_live import _frozen_pair_preflight
from tests.audit_method1_4310_live import _frozen_preflight_summary
from tests.audit_method1_4310_live import _frozen_report_status
from tests.audit_method1_4310_live import _preparse_official_mode_references
from tests.audit_method1_4310_live import _resign_source_and_invalidate
from tests.audit_method1_4310_live import _source_record_set_integrity
from tests.official_html_resolver import scan_official_html
from tests.support_validate_method_outputs import _complete_modes_html
from tests.support_validate_method_outputs import _write_setting_cif
import json
import numpy as np
import pytest
import tests.audit_method1_4310_live as method1_4310_module


# Original coverage: test_method_validation.py.

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
