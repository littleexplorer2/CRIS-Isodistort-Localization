from __future__ import annotations

from pathlib import Path

import pytest

from tests.manual.official_html_resolver import (
    OfficialHtmlResolutionError,
    classify_official_html,
    scan_official_html,
)


def _form(*inputs: tuple[str, str]) -> str:
    controls = "".join(
        f'<input type="hidden" name="{name}" value="{value}">' for name, value in inputs
    )
    return f'<form method="POST" action="isodistortform.php">{controls}</form>'


def _page(title: str, body: str) -> str:
    return f"<html><head><title>{title}</title></head><body><h1>{title}</h1>{body}</body></html>"


@pytest.mark.parametrize(
    ("title", "body", "expected_role"),
    [
        (
            " ISODISTORT:\n complete   modes details ",
            "<b>Subgroup details</b><b>Undistorted superstructure</b>"
            "<b>Distorted superstructure</b>",
            "complete_modes_details",
        ),
        (
            "ISODISTORT: distortion",
            _form(("input", "displaydistort")),
            "distortion_result",
        ),
        (
            "ISODISTORT: distorted structure (basis)",
            _form(("input", "distort"), ("origintype", "method4")),
            "distorted_structure_basis",
        ),
        (
            "ISODISTORT: search",
            _form(("input", "changesearch")) + _form(("input", "kvector")),
            "search",
        ),
    ],
)
def test_classifies_pages_from_normalized_content(
    tmp_path: Path,
    title: str,
    body: str,
    expected_role: str,
) -> None:
    path = tmp_path / "arbitrary renamed page.HTM"
    path.write_text(_page(title, body), encoding="utf-8")

    page = classify_official_html(path)

    assert page.role == expected_role
    assert page.path == path
    assert page.raw_html == page.text
    assert page.evidence


@pytest.mark.parametrize(
    "body",
    [
        '<a href="details.html">ISODISTORT: complete modes details</a>',
        '<a href="distortion.html">ISODISTORT: distortion</a>',
        '<a href="basis.html">ISODISTORT: distorted structure (basis)</a>',
        '<a href="search.html">ISODISTORT: search</a>',
    ],
)
def test_link_text_alone_never_classifies_page(tmp_path: Path, body: str) -> None:
    path = tmp_path / "misleading-name.html"
    path.write_text(_page("Unrelated page", body), encoding="utf-8")

    assert classify_official_html(path).role == "unknown"


def test_title_without_required_structure_is_unknown(tmp_path: Path) -> None:
    path = tmp_path / "complete modes details.html"
    path.write_text(
        _page("ISODISTORT: complete modes details", "<p>Subgroup details</p>"),
        encoding="utf-8",
    )

    assert classify_official_html(path).role == "unknown"


def test_distortion_rejection_page_without_result_form_is_classified(tmp_path: Path) -> None:
    path = tmp_path / "server-rejection.HTML"
    path.write_text(
        _page(
            "ISODISTORT: distortion",
            "<p>Types of atoms in subgroup do not match types of atoms in parent.</p>",
        ),
        encoding="utf-8",
    )

    page = classify_official_html(path)

    assert page.role == "distortion_result"
    assert "result_form:absent_expected_rejection_shape" in page.evidence


def test_scan_scope_and_suffixes_are_explicit(tmp_path: Path) -> None:
    direct = tmp_path / "renamed.HTML"
    nested_dir = tmp_path / "candidate"
    nested_dir.mkdir()
    nested = nested_dir / "also-renamed.hTm"
    ignored = tmp_path / "looks-like.html.txt"
    html = _page("ISODISTORT: distortion", _form(("input", "displaydistort")))
    direct.write_text(html, encoding="utf-8")
    nested.write_text(html, encoding="utf-8")
    ignored.write_text(html, encoding="utf-8")

    shallow = scan_official_html(tmp_path, recursive=False)
    recursive = scan_official_html(tmp_path, recursive=True)

    assert [page.path for page in shallow.pages] == [direct]
    assert [page.path for page in recursive.pages] == [nested, direct]
    assert shallow.require_unique("distortion_result").path == direct
    assert recursive.in_directory(nested_dir) == (classify_official_html(nested, root=tmp_path),)
    assert (
        recursive.require_unique_in_directory("distortion_result", nested_dir).path
        == nested
    )


def test_directory_unique_resolution_never_leaks_descendant_pages(
    tmp_path: Path,
) -> None:
    child = tmp_path / "auxiliary-run"
    child.mkdir()
    html = _page("ISODISTORT: distortion", _form(("input", "displaydistort")))
    primary = tmp_path / "primary.html"
    auxiliary = child / "auxiliary.html"
    primary.write_text(html, encoding="utf-8")
    auxiliary.write_text(html, encoding="utf-8")

    inventory = scan_official_html(tmp_path, recursive=True)

    assert inventory.require_unique_in_directory("distortion_result", tmp_path).path == primary
    assert inventory.require_unique_in_directory("distortion_result", child).path == auxiliary


def test_require_unique_missing_error_lists_inspected_paths(tmp_path: Path) -> None:
    unknown = tmp_path / "unknown.html"
    unknown.write_text(_page("Unrelated", ""), encoding="utf-8")
    inventory = scan_official_html(tmp_path, recursive=False)

    with pytest.raises(OfficialHtmlResolutionError) as exc_info:
        inventory.require_unique("complete_modes_details")

    error = exc_info.value
    assert error.code == "missing_match"
    assert error.matched_paths == ()
    assert error.inspected_paths == (unknown,)
    assert str(unknown) in str(error)


def test_require_unique_ambiguous_error_lists_every_match(tmp_path: Path) -> None:
    html = _page("ISODISTORT: distortion", _form(("input", "displaydistort")))
    paths = (tmp_path / "first.html", tmp_path / "second.htm")
    for path in paths:
        path.write_text(html, encoding="utf-8")
    inventory = scan_official_html(tmp_path, recursive=False)

    with pytest.raises(OfficialHtmlResolutionError) as exc_info:
        inventory.require_unique("distortion_result")

    error = exc_info.value
    assert error.code == "ambiguous_match"
    assert error.matched_paths == paths
    assert all(str(path) in str(error) for path in paths)


def test_signature_and_parsing_share_captured_content(tmp_path: Path) -> None:
    path = tmp_path / "download.html"
    original = _page("ISODISTORT: distortion", _form(("input", "displaydistort")))
    path.write_text(original, encoding="utf-8")
    inventory = scan_official_html(tmp_path, recursive=False)
    captured = inventory.require_unique("distortion_result")

    path.write_text(_page("Changed", ""), encoding="utf-8")

    assert captured.raw_html == original
    assert captured.role == "distortion_result"
    assert inventory.require_unique("distortion_result") is captured


def test_signature_ignores_basename_but_tracks_content_and_parent(tmp_path: Path) -> None:
    html = _page("ISODISTORT: distortion", _form(("input", "displaydistort")))
    path = tmp_path / "original.html"
    path.write_text(html, encoding="utf-8")
    original = scan_official_html(tmp_path, recursive=True)

    renamed = tmp_path / "manually renamed.HTML"
    path.rename(renamed)
    after_rename = scan_official_html(tmp_path, recursive=True)
    assert after_rename.signature == original.signature

    child = tmp_path / "different-candidate"
    child.mkdir()
    moved = child / renamed.name
    renamed.rename(moved)
    after_move = scan_official_html(tmp_path, recursive=True)
    assert after_move.signature != original.signature

    moved.write_text(html + "<!-- changed -->", encoding="utf-8")
    after_change = scan_official_html(tmp_path, recursive=True)
    assert after_change.signature != after_move.signature


def test_unknown_pages_and_duplicates_contribute_to_signature(tmp_path: Path) -> None:
    known = tmp_path / "known.html"
    unknown = tmp_path / "unknown.html"
    html = _page("ISODISTORT: distortion", _form(("input", "displaydistort")))
    known.write_text(html, encoding="utf-8")
    unknown.write_text(_page("Unrelated", "one"), encoding="utf-8")
    first = scan_official_html(tmp_path, recursive=False)

    unknown.write_text(_page("Unrelated", "two"), encoding="utf-8")
    changed_unknown = scan_official_html(tmp_path, recursive=False)
    assert changed_unknown.signature != first.signature

    duplicate = tmp_path / "duplicate.htm"
    duplicate.write_text(html, encoding="utf-8")
    with_duplicate = scan_official_html(tmp_path, recursive=False)
    assert with_duplicate.signature != changed_unknown.signature


def test_non_html_file_and_invalid_root_fail_closed(tmp_path: Path) -> None:
    text_file = tmp_path / "page.txt"
    text_file.write_text("<html></html>", encoding="utf-8")

    with pytest.raises(ValueError, match="not an HTML file"):
        classify_official_html(text_file)
    with pytest.raises(NotADirectoryError, match="not a directory"):
        scan_official_html(text_file, recursive=False)
