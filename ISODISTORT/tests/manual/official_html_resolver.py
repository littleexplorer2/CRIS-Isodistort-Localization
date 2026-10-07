"""Content-driven discovery of archived official ISODISTORT HTML pages.

The archive is user-supplied validation evidence. Browser-assigned filenames
are therefore not part of a page's identity: a page is classified from its
HTML structure, and callers must explicitly choose whether discovery is local
to one directory or recursive. An :class:`OfficialHtmlInventory` retains the
exact decoded content that was hashed and classified so downstream parsers do
not need to reopen a possibly changed file.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path
from typing import Literal, TypeAlias

HtmlPageRole: TypeAlias = Literal[
    "complete_modes_details",
    "distortion_result",
    "distorted_structure_basis",
    "search",
    "unknown",
]

_KNOWN_ROLES = frozenset(
    {
        "complete_modes_details",
        "distortion_result",
        "distorted_structure_basis",
        "search",
        "unknown",
    }
)
_HTML_SUFFIXES = frozenset({".html", ".htm"})


def _normalize_text(value: str) -> str:
    return " ".join(value.split()).casefold()


@dataclass(slots=True)
class _Form:
    method: str
    action: str
    inputs: list[tuple[str, str]] = field(default_factory=list)


class _PageStructureParser(HTMLParser):
    """Collect only the structural fields used to identify an official page."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self._title: list[str] | None = None
        self._h1: list[str] | None = None
        self._strong: list[str] | None = None
        self._current_form: _Form | None = None
        self.titles: list[str] = []
        self.h1s: list[str] = []
        self.strong_text: list[str] = []
        self.forms: list[_Form] = []

    def handle_starttag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        tag = tag.casefold()
        attributes = {name.casefold(): value or "" for name, value in attrs}
        if tag == "title":
            self._title = []
        elif tag == "h1":
            self._h1 = []
        elif tag in {"b", "strong"}:
            self._strong = []
        elif tag == "form":
            form = _Form(
                method=_normalize_text(attributes.get("method", "")),
                action=attributes.get("action", "").strip(),
            )
            self.forms.append(form)
            self._current_form = form
        elif tag == "input" and self._current_form is not None:
            name = _normalize_text(attributes.get("name", ""))
            value = _normalize_text(attributes.get("value", ""))
            if name:
                self._current_form.inputs.append((name, value))

    def handle_startendtag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        self.handle_starttag(tag, attrs)

    def handle_data(self, data: str) -> None:
        if self._title is not None:
            self._title.append(data)
        if self._h1 is not None:
            self._h1.append(data)
        if self._strong is not None:
            self._strong.append(data)

    def handle_endtag(self, tag: str) -> None:
        tag = tag.casefold()
        if tag == "title" and self._title is not None:
            self.titles.append(_normalize_text("".join(self._title)))
            self._title = None
        elif tag == "h1" and self._h1 is not None:
            self.h1s.append(_normalize_text("".join(self._h1)))
            self._h1 = None
        elif tag in {"b", "strong"} and self._strong is not None:
            self.strong_text.append(_normalize_text("".join(self._strong)))
            self._strong = None
        elif tag == "form":
            self._current_form = None


@dataclass(frozen=True, slots=True)
class OfficialHtmlPage:
    """One HTML page, captured and classified from one immutable byte read."""

    path: Path
    relative_path: Path
    role: HtmlPageRole
    title: str | None
    raw_html: str = field(repr=False)
    content_sha256: str
    evidence: tuple[str, ...]

    @property
    def text(self) -> str:
        """Return the captured raw decoded HTML (not visible-text extraction)."""

        return self.raw_html


class OfficialHtmlResolutionError(ValueError):
    """A structured fail-closed error for a missing or ambiguous page role."""

    def __init__(
        self,
        *,
        code: Literal["missing_match", "ambiguous_match"],
        role: HtmlPageRole,
        root: Path,
        matched_paths: tuple[Path, ...],
        inspected_paths: tuple[Path, ...],
    ) -> None:
        self.code = code
        self.role = role
        self.root = root
        self.matched_paths = matched_paths
        self.inspected_paths = inspected_paths
        matched = ", ".join(str(path) for path in matched_paths) or "<none>"
        inspected = ", ".join(str(path) for path in inspected_paths) or "<none>"
        super().__init__(
            f"{code}: expected exactly one {role!r} HTML page under {root}; "
            f"matched=[{matched}]; inspected=[{inspected}]"
        )


@dataclass(frozen=True, slots=True)
class OfficialHtmlInventory:
    """Stable inventory used for both provenance signatures and page parsing."""

    root: Path
    recursive: bool
    pages: tuple[OfficialHtmlPage, ...]
    signature: str

    def by_role(self, role: HtmlPageRole) -> tuple[OfficialHtmlPage, ...]:
        _validate_role(role)
        return tuple(page for page in self.pages if page.role == role)

    def in_directory(self, directory: Path) -> tuple[OfficialHtmlPage, ...]:
        """Return captured pages stored directly in ``directory``.

        Filtering already-captured objects prevents a second filesystem read
        between provenance signing and scientific parsing. Descendant
        directories are excluded so auxiliary runs cannot make a primary
        candidate's page role ambiguous.
        """

        resolved = Path(directory).resolve()
        return tuple(page for page in self.pages if page.path.parent.resolve() == resolved)

    def require_unique(self, role: HtmlPageRole) -> OfficialHtmlPage:
        """Return the sole matching page, otherwise raise a structured error."""

        _validate_role(role)
        matches = self.by_role(role)
        if len(matches) == 1:
            return matches[0]
        code: Literal["missing_match", "ambiguous_match"] = (
            "missing_match" if not matches else "ambiguous_match"
        )
        raise OfficialHtmlResolutionError(
            code=code,
            role=role,
            root=self.root,
            matched_paths=tuple(page.path for page in matches),
            inspected_paths=tuple(page.path for page in self.pages),
        )

    def require_unique_in_directory(
        self,
        role: HtmlPageRole,
        directory: Path,
    ) -> OfficialHtmlPage:
        """Return one role from one physical directory of this capture."""

        _validate_role(role)
        pages = self.in_directory(directory)
        matches = tuple(page for page in pages if page.role == role)
        if len(matches) == 1:
            return matches[0]
        code: Literal["missing_match", "ambiguous_match"] = (
            "missing_match" if not matches else "ambiguous_match"
        )
        raise OfficialHtmlResolutionError(
            code=code,
            role=role,
            root=Path(directory),
            matched_paths=tuple(page.path for page in matches),
            inspected_paths=tuple(page.path for page in pages),
        )


def _validate_role(role: str) -> None:
    if role not in _KNOWN_ROLES:
        raise ValueError(f"unsupported official HTML role: {role!r}")


def _form_has_input(parser: _PageStructureParser, name: str, value: str) -> bool:
    expected = (_normalize_text(name), _normalize_text(value))
    return any(
        form.method == "post" and expected in form.inputs
        for form in parser.forms
    )


def _classify_structure(parser: _PageStructureParser) -> tuple[HtmlPageRole, tuple[str, ...]]:
    title = parser.titles[0] if len(parser.titles) == 1 else None
    h1s = set(parser.h1s)
    strong = set(parser.strong_text)

    complete_title = "isodistort: complete modes details"
    if (
        title == complete_title
        and complete_title in h1s
        and {
            "subgroup details",
            "undistorted superstructure",
            "distorted superstructure",
        }.issubset(strong)
    ):
        return (
            "complete_modes_details",
            (
                "title:isodistort_complete_modes_details",
                "h1:isodistort_complete_modes_details",
                "sections:subgroup_and_superstructures",
            ),
        )

    distortion_title = "isodistort: distortion"
    if title == distortion_title and distortion_title in h1s:
        evidence = [
            "title:isodistort_distortion",
            "h1:isodistort_distortion",
        ]
        if _form_has_input(parser, "input", "displaydistort"):
            evidence.append("post_input:displaydistort")
        else:
            evidence.append("result_form:absent_expected_rejection_shape")
        return "distortion_result", tuple(evidence)

    basis_title = "isodistort: distorted structure (basis)"
    if (
        title == basis_title
        and basis_title in h1s
        and _form_has_input(parser, "input", "distort")
        and _form_has_input(parser, "origintype", "method4")
    ):
        return (
            "distorted_structure_basis",
            (
                "title:isodistort_distorted_structure_basis",
                "h1:isodistort_distorted_structure_basis",
                "post_input:distort",
                "post_origintype:method4",
            ),
        )

    search_title = "isodistort: search"
    method_inputs = {"isosubgroup", "kvector", "parentbasis"}
    search_inputs = {
        value
        for form in parser.forms
        if form.method == "post"
        for name, value in form.inputs
        if name == "input"
    }
    if (
        title == search_title
        and search_title in h1s
        and "changesearch" in search_inputs
        and method_inputs.intersection(search_inputs)
    ):
        return (
            "search",
            (
                "title:isodistort_search",
                "h1:isodistort_search",
                "post_input:changesearch",
                "post_input:search_method",
            ),
        )

    return "unknown", ()


def classify_official_html(path: Path, *, root: Path | None = None) -> OfficialHtmlPage:
    """Read and classify one HTML file without consulting its basename."""

    path = Path(path)
    if path.suffix.casefold() not in _HTML_SUFFIXES:
        raise ValueError(f"not an HTML file: {path}")
    payload = path.read_bytes()
    raw_html = payload.decode("utf-8-sig", errors="replace")
    parser = _PageStructureParser()
    parser.feed(raw_html)
    parser.close()
    role, evidence = _classify_structure(parser)
    relative_path = path.name if root is None else path.relative_to(root)
    title = parser.titles[0] if len(parser.titles) == 1 else None
    return OfficialHtmlPage(
        path=path,
        relative_path=Path(relative_path),
        role=role,
        title=title,
        raw_html=raw_html,
        content_sha256=hashlib.sha256(payload).hexdigest(),
        evidence=evidence,
    )


def _inventory_signature(pages: tuple[OfficialHtmlPage, ...]) -> str:
    # Basenames are excluded so a browser/manual rename cannot invalidate
    # otherwise identical evidence. The relative parent remains significant.
    entries = sorted(
        (
            page.relative_path.parent.as_posix().casefold(),
            page.role,
            page.content_sha256,
        )
        for page in pages
    )
    canonical = json.dumps(entries, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def scan_official_html(root: Path, *, recursive: bool) -> OfficialHtmlInventory:
    """Scan ``root`` with explicit current-directory or recursive scope."""

    root = Path(root)
    if not root.is_dir():
        raise NotADirectoryError(f"official HTML root is not a directory: {root}")
    iterator = root.rglob("*") if recursive else root.iterdir()
    paths = sorted(
        (
            path
            for path in iterator
            if path.is_file() and path.suffix.casefold() in _HTML_SUFFIXES
        ),
        key=lambda path: (
            path.relative_to(root).as_posix().casefold(),
            path.relative_to(root).as_posix(),
        ),
    )
    pages = tuple(classify_official_html(path, root=root) for path in paths)
    return OfficialHtmlInventory(
        root=root,
        recursive=recursive,
        pages=pages,
        signature=_inventory_signature(pages),
    )


__all__ = [
    "HtmlPageRole",
    "OfficialHtmlInventory",
    "OfficialHtmlPage",
    "OfficialHtmlResolutionError",
    "classify_official_html",
    "scan_official_html",
]
