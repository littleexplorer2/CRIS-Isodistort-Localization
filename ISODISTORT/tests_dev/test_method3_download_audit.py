"""Tests for the read-only Method 3 download auditor."""

from __future__ import annotations

import json
from pathlib import Path

from tests_dev.manual import audit_method3_downloads as audit_module
from tests_dev.manual.audit_method3_downloads import (
    audit_method3_root,
    parse_details_page_identity,
    parse_result_page,
    parse_result_page_context,
    parse_subgroup_cif,
)

IDENTITY = "12 C2/m, basis={(0,1,1),(-1,0,0),(0,-1,0)}, origin=(-1/2,0,1/4), s=1, i=4"
SECOND_IDENTITY = "8 Cm, basis={(-1,1,0),(0,0,1),(1,0,0)}, origin=(0,0,0), s=1, i=8"


def _write_candidate(
    directory: Path,
    *,
    identity: str = IDENTITY,
    omit: str | None = None,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    files = {
        "subgroup.cif": f"# Subgroup: {identity}\ndata_test\n",
        "data.isoviz": "isoviz\n",
        "topas.str": "topas\n",
        "ISODISTORT_ complete modes details.html": (
            "<html><h1>complete modes details</h1><pre><b>Subgroup details</b>\n"
            f"{identity}\n</pre></html>\n"
        ),
    }
    for name, content in files.items():
        if name != omit:
            (directory / name).write_text(content, encoding="utf-8")


def _write_context_result(
    path: Path,
    *,
    identity: str = IDENTITY,
    parent_sg: int = 139,
    target_sg: int = 12,
    target_symbol: str = "C2/m",
    basis: str = "(0,1,1),(-1,0,0),(0,-1,0)",
    centering: str = "C",
    distortion_types: tuple[str, ...] = ("strain", "displacive"),
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    effective_centering = "no" if centering.upper() == "P" else centering.upper()
    inputs = [f'<input type="hidden" name="spacegroup" value="{parent_sg} Parent">']
    if "strain" in distortion_types:
        inputs.append('<input type="hidden" name="includestrain" value="true">')
    input_names = {
        "displacive": "includedisplacive",
        "occupational": "includescalar",
        "magnetic": "includemagnetic",
        "rotational": "includerotational",
    }
    for distortion_type, input_name in input_names.items():
        if distortion_type in distortion_types:
            inputs.append(
                f'<input type="hidden" name="{input_name}001" value="true">'
            )
    path.write_text(
        "<html><body>"
        f"Space Group: {parent_sg} Parent<br>"
        f"Include {', '.join(distortion_types)} ALL distortions<br>"
        f"Subgroup: {target_sg} {target_symbol}, basis of lattice = {basis}, "
        f"{effective_centering} centering<br>"
        + "".join(inputs)
        + "Finish selecting the distortion mode:"
        + f'<input type="radio" name="orderparam" value="x *Subgroup: {identity}">'
        + f" {identity}<br></body></html>",
        encoding="utf-8",
    )


def _manifest(*cases: dict, parent_sg: int = 139) -> dict:
    return {
        "parents": [
            {
                "parent_cif": "Parent.cif",
                "parent_space_group_type": parent_sg,
                "cases": list(cases),
            }
        ]
    }


def _case(
    folder: str,
    *,
    target_sg: int = 12,
    basis: list[list[int]] | None = None,
    centering: str = "d",
    distortion_types: list[str] | None = None,
) -> dict:
    return {
        "id": folder,
        "readable_case_folder": folder,
        "space_group_type": target_sg,
        "supercell_basis": basis or [[0, 1, 1], [-1, 0, 0], [0, -1, 0]],
        "direct_sublattice_centering": centering,
        "distortion_types": distortion_types or ["strain", "displacive"],
    }


def test_result_page_and_cif_parse_to_same_fractional_identity(tmp_path: Path) -> None:
    page = tmp_path / "result.html"
    page.write_text(
        "<html>Finish selecting the distortion mode:"
        f'<input type="RADIO" name="orderparam" value="ignored *Subgroup: {IDENTITY}">'
        f" {IDENTITY}<br></html>",
        encoding="utf-8",
    )
    cif = tmp_path / "subgroup.cif"
    cif.write_text(f"# Subgroup:  {IDENTITY}\n", encoding="utf-8")

    rows = parse_result_page(page)
    assert rows == [parse_subgroup_cif(cif)]
    assert rows[0].origin[0].numerator == -1
    assert rows[0].origin[0].denominator == 2


def test_audit_matches_result_table_and_checks_four_core_file_classes(tmp_path: Path) -> None:
    case = tmp_path / "M3-TEST-01"
    case.mkdir()
    (case / "result.html").write_text(
        f'Finish selecting<input name="orderparam" type="radio" value="x *Subgroup: {IDENTITY}">{IDENTITY}<br>',
        encoding="utf-8",
    )
    _write_candidate(case / "candidate")

    audit = audit_method3_root(tmp_path)

    assert audit.error_count == 0
    assert audit.cases[0].status == "verified"
    assert len(audit.cases[0].candidates[0].core_files) == 4


def test_audit_pairs_swapped_web_packages_by_exact_identity(tmp_path: Path) -> None:
    """Content-complete downloads remain evidence despite a browser folder swap."""
    case = tmp_path / "M3-TEST-SWAPPED-PACKAGES"
    candidate = case / "candidate"
    _write_candidate(candidate, omit="ISODISTORT_ complete modes details.html")
    _write_context_result(candidate / "method3-result.html")
    details = case / "ISODISTORT_ complete modes details.html"
    details.write_text(
        "<html><h1>Complete modes details</h1><pre><b>Subgroup details</b>\n"
        f"{IDENTITY}\n</pre></html>",
        encoding="utf-8",
    )

    audit = audit_method3_root(tmp_path)
    result = audit.cases[0]
    candidate_result = result.candidates[0]
    warning_codes = {
        issue.code for issue in result.issues + candidate_result.issues
        if issue.severity == "warning"
    }

    assert audit.error_count == 0
    assert result.status == "verified"
    assert parse_details_page_identity(details) == parse_subgroup_cif(
        candidate / "subgroup.cif"
    )
    assert candidate_result.core_files["complete_modes_details"] == str(details)
    assert {
        "noncanonical_result_page_location",
        "noncanonical_details_page_location",
    } <= warning_codes


def test_explicit_zero_candidate_result_is_verified_without_downloads(tmp_path: Path) -> None:
    case = tmp_path / "M3-TEST-EMPTY"
    case.mkdir()
    page = case / "ISODISTORT_ distorted structure_ subgroup.html"
    page.write_text(
        "<html><body>Finish selecting the distortion mode by selecting one of the following choices:"
        "<p>There are <strong>no subgroups</strong> with the space-group symmetry you selected.</p>"
        "<p>Try again.</p></body></html>",
        encoding="utf-8",
    )

    audit = audit_method3_root(tmp_path)
    result = audit.cases[0]
    codes = {issue.code for issue in result.issues}

    assert audit.error_count == 0
    assert result.status == "verified"
    assert result.html_pages[page.name] == "empty_result"
    assert result.result_pages == [str(page)]
    assert result.official_candidates == []
    assert result.downloaded_candidates == []
    assert "missing_result_page" not in codes
    assert "search_page_saved_instead_of_result" not in codes
    assert "no_downloaded_candidates" not in codes
    assert parse_result_page(page) == []


def test_audit_flags_search_page_saved_instead_of_result_and_missing_core(tmp_path: Path) -> None:
    case = tmp_path / "M3-TEST-02"
    case.mkdir()
    (case / "ISODISTORT_ search.html").write_text(
        "<h2>Method 3: Search over arbitrary k points for specified space group and lattice</h2>",
        encoding="utf-8",
    )
    _write_candidate(case, omit="topas.str")

    audit = audit_method3_root(tmp_path)
    codes = {issue.code for issue in audit.cases[0].issues}
    candidate_codes = {issue.code for issue in audit.cases[0].candidates[0].issues}

    assert audit.cases[0].status == "unverified_result_page"
    assert "search_page_saved_instead_of_result" in codes
    assert "missing_core_file" in candidate_codes


def test_manifest_binds_result_context_to_intended_case(tmp_path: Path) -> None:
    case = tmp_path / "case-one"
    _write_context_result(case / "result.html")
    _write_candidate(case / "candidate")
    manifest = _manifest(_case("case-one"))

    audit = audit_method3_root(tmp_path, manifest=manifest, parent_cif="Parent.cif")

    assert audit.error_count == 0
    assert audit.cases[0].status == "verified"
    context = parse_result_page_context(case / "result.html")
    assert context.parent_space_group_number == 139
    assert context.target_space_group_number == 12
    assert context.centering == "C"
    assert context.distortion_types == frozenset({"strain", "displacive"})
    assert audit.cases[0].result_page_contexts[str(case / "result.html")]["centering"] == "C"


def test_result_context_recognizes_official_scalar_name_as_occupational(tmp_path: Path) -> None:
    page = tmp_path / "result.html"
    _write_context_result(
        page,
        distortion_types=("occupational",),
    )

    assert parse_result_page_context(page).distortion_types == frozenset({"occupational"})


def test_manifest_detects_complete_result_and_exports_swapped_between_cases(tmp_path: Path) -> None:
    first = tmp_path / "case-one"
    second = tmp_path / "case-two"
    # Each directory is internally self-consistent, so the generic audit alone
    # cannot know that both complete payloads were swapped.
    _write_context_result(
        first / "result.html",
        identity=SECOND_IDENTITY,
        target_sg=8,
        target_symbol="Cm",
        basis="(-1,1,0),(0,0,1),(1,0,0)",
    )
    _write_candidate(first / "candidate", identity=SECOND_IDENTITY)
    _write_context_result(second / "result.html")
    _write_candidate(second / "candidate")
    manifest = _manifest(
        _case("case-one"),
        _case("case-two", target_sg=8, basis=[[-1, 1, 0], [0, 0, 1], [1, 0, 0]]),
    )

    assert audit_method3_root(tmp_path).error_count == 0
    bound = audit_method3_root(tmp_path, manifest=manifest, parent_cif="Parent.cif")
    first_codes = {issue.code for issue in bound.cases[0].issues}
    second_codes = {issue.code for issue in bound.cases[1].issues}

    assert bound.error_count >= 4
    assert {"result_target_space_group_mismatch", "result_basis_mismatch"} <= first_codes
    assert {"result_target_space_group_mismatch", "result_basis_mismatch"} <= second_codes
    assert all(case.status == "failed" for case in bound.cases)


def test_manifest_checks_parent_centering_and_distortion_types(tmp_path: Path) -> None:
    case = tmp_path / "case-one"
    _write_context_result(
        case / "result.html",
        parent_sg=123,
        centering="P",
        distortion_types=("strain",),
    )
    _write_candidate(case / "candidate")

    audit = audit_method3_root(
        tmp_path,
        manifest=_manifest(_case("case-one")),
        parent_cif="Parent.cif",
    )
    codes = {issue.code for issue in audit.cases[0].issues}

    assert "result_parent_space_group_mismatch" in codes
    assert "result_centering_mismatch" in codes
    assert "result_distortion_types_mismatch" in codes
    assert audit.cases[0].status == "failed"


def test_explicit_manifest_centering_is_compared_to_effective_page_centering(tmp_path: Path) -> None:
    case = tmp_path / "case-one"
    _write_context_result(case / "result.html")
    _write_candidate(case / "candidate")

    matching = audit_method3_root(
        tmp_path,
        manifest=_manifest(_case("case-one", centering="C")),
        parent_cif="Parent.cif",
    )
    mismatching = audit_method3_root(
        tmp_path,
        manifest=_manifest(_case("case-one", centering="P")),
        parent_cif="Parent.cif",
    )

    assert matching.error_count == 0
    assert "result_centering_mismatch" in {issue.code for issue in mismatching.cases[0].issues}


def test_cli_uses_project_manifest_by_default_but_api_can_omit_it(
    tmp_path: Path,
    monkeypatch,
) -> None:
    case = tmp_path / "case-one"
    _write_context_result(case / "result.html", parent_sg=123)
    _write_candidate(case / "candidate")
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(_manifest(_case("case-one"))), encoding="utf-8")
    monkeypatch.setattr(audit_module, "DEFAULT_MANIFEST_PATH", manifest_path)

    assert audit_module.main([str(tmp_path), "--parent-cif", "Parent.cif"]) == 1
    assert audit_module.main([str(tmp_path), "--no-manifest"]) == 0
    assert audit_method3_root(tmp_path).error_count == 0
