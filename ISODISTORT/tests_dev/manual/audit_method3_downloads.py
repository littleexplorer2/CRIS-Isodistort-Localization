"""Audit saved official Method 3 result pages and candidate downloads.

The auditor is deliberately read-only and independent of any crystal-specific
answer list.  A case's authoritative candidate list is read from the saved
ISODISTORT result-table HTML (the page containing both ``Finish selecting`` and
``orderparam`` radio controls), then compared with identities in downloaded
``subgroup.cif`` files.

When a Method 3 manifest is supplied, the auditor also checks that each saved
result page belongs to the intended query.  This closes a gap where a result
page and all four matching exports could be moved together into the wrong case
directory and still pass the file-to-page comparison.

Usage::

    python tests_dev/manual/audit_method3_downloads.py PATH/TO/Method3
    python tests_dev/manual/audit_method3_downloads.py PATH/TO/Method3 \
        --json-output audit.json

The CLI uses the project's Method 3 manifest by default.  Library callers can
omit ``manifest`` to retain the generic, manifest-independent audit.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

CORE_FILENAMES = ("subgroup.cif", "data.isoviz", "topas.str")
RESULT_MARKER = "Finish selecting"
EMPTY_RESULT_MARKER = "There are no subgroups"
TRY_AGAIN_MARKER = "Try again"
PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method3_download_manifest.json"


def _read_text(path: Path) -> str:
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "cp1252"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _vector_text(values: Iterable[Fraction]) -> str:
    return "(" + ",".join(_fraction_text(value) for value in values) + ")"


@dataclass(frozen=True, order=True)
class SubgroupIdentity:
    """Identity printed by ISODISTORT for one embedded subgroup."""

    space_group_number: int
    symbol: str
    basis: tuple[tuple[Fraction, Fraction, Fraction], ...]
    origin: tuple[Fraction, Fraction, Fraction]
    sublattice_index: int
    group_index: int

    @property
    def text(self) -> str:
        basis = "{" + ",".join(_vector_text(row) for row in self.basis) + "}"
        return (
            f"{self.space_group_number} {self.symbol}, basis={basis}, "
            f"origin={_vector_text(self.origin)}, s={self.sublattice_index}, i={self.group_index}"
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "space_group_number": self.space_group_number,
            "symbol": self.symbol,
            "basis": [[_fraction_text(value) for value in row] for row in self.basis],
            "origin": [_fraction_text(value) for value in self.origin],
            "s": self.sublattice_index,
            "i": self.group_index,
            "text": self.text,
        }


@dataclass(frozen=True)
class ResultPageContext:
    """Method 3 query context echoed by one official result page."""

    parent_space_group_number: int
    target_space_group_number: int
    target_space_group_symbol: str
    basis: tuple[tuple[Fraction, Fraction, Fraction], ...]
    centering: str
    distortion_types: frozenset[str]

    def to_json(self) -> dict[str, Any]:
        return {
            "parent_space_group_number": self.parent_space_group_number,
            "target_space_group_number": self.target_space_group_number,
            "target_space_group_symbol": self.target_space_group_symbol,
            "basis": [[_fraction_text(value) for value in row] for row in self.basis],
            "centering": self.centering,
            "distortion_types": sorted(self.distortion_types),
        }


@dataclass(frozen=True)
class Method3CaseExpectation:
    """Manifest fields that can be proved from an official result page."""

    case: str
    parent_cif: str
    parent_space_group_number: int
    target_space_group_number: int
    basis: tuple[tuple[Fraction, Fraction, Fraction], ...]
    centering: str
    distortion_types: frozenset[str]


_IDENTITY_RE = re.compile(
    r"^\s*(?:Subgroup:\s*)?"
    r"(?P<number>\d+)\s+(?P<symbol>.+?),\s*"
    r"basis=\{(?P<basis>.*?)\},\s*"
    r"origin=\((?P<origin>[^()]*)\),\s*"
    r"s=(?P<s>\d+)\s*,\s*i=(?P<i>\d+)\s*$",
    re.IGNORECASE,
)
_VECTOR_RE = re.compile(r"\(([^()]*)\)")
_PARENT_SPACE_GROUP_RE = re.compile(r"\bSpace\s+Group:\s*(\d+)\b", re.IGNORECASE)
_RESULT_CONTEXT_RE = re.compile(
    r"\bSubgroup:\s*(?P<number>\d+)\s+(?P<symbol>.*?),\s*"
    r"basis\s+of\s+lattice\s*=\s*"
    r"(?P<basis>\([^()]+\)\s*,\s*\([^()]+\)\s*,\s*\([^()]+\))\s*,\s*"
    r"(?P<center>no|[PABCIFR])\s+centering\b",
    re.IGNORECASE,
)
_DISTORTION_INPUT_PREFIXES = {
    "includedisplacive": "displacive",
    "includescalar": "occupational",
    "includemagnetic": "magnetic",
    "includerotational": "rotational",
}


def _parse_vector(text: str, *, size: int = 3) -> tuple[Fraction, ...]:
    values = tuple(Fraction(part.strip()) for part in text.split(","))
    if len(values) != size:
        raise ValueError(f"expected {size} coordinates, got {len(values)} in {text!r}")
    return values


def parse_subgroup_identity(text: str) -> SubgroupIdentity:
    """Parse one visible Method 3 result row or a CIF ``# Subgroup`` line."""

    cleaned = re.sub(r"^\s*#\s*Subgroup:\s*", "", text.strip(), flags=re.IGNORECASE)
    match = _IDENTITY_RE.fullmatch(cleaned)
    if match is None:
        raise ValueError(f"not an ISODISTORT subgroup identity: {text!r}")

    rows = tuple(_parse_vector(row) for row in _VECTOR_RE.findall(match.group("basis")))
    if len(rows) != 3:
        raise ValueError(f"expected three basis vectors in {text!r}")
    basis = tuple(tuple(row) for row in rows)
    origin = tuple(_parse_vector(match.group("origin")))
    symbol = " ".join(match.group("symbol").split())
    return SubgroupIdentity(
        space_group_number=int(match.group("number")),
        symbol=symbol,
        basis=basis,  # type: ignore[arg-type]
        origin=origin,  # type: ignore[arg-type]
        sublattice_index=int(match.group("s")),
        group_index=int(match.group("i")),
    )


class _OrderParamHTMLParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.rows: list[str] = []
        self.visible_text: list[str] = []
        self.inputs: list[dict[str, str | None]] = []
        self._current: list[str] | None = None
        self._fallback: str | None = None

    def _finish_current(self) -> None:
        if self._current is None:
            return
        visible = " ".join("".join(self._current).split())
        row = visible
        if not row and self._fallback:
            row = self._fallback.split("*Subgroup:", maxsplit=1)[-1].strip()
        if row:
            self.rows.append(row)
        self._current = None
        self._fallback = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "input":
            values = {name.lower(): value for name, value in attrs}
            self.inputs.append(values)
            if (values.get("name") or "").lower() == "orderparam":
                self._finish_current()
                self._current = []
                self._fallback = values.get("value")
        elif tag.lower() == "br":
            self._finish_current()

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.handle_starttag(tag, attrs)

    def handle_data(self, data: str) -> None:
        self.visible_text.append(data)
        if self._current is not None:
            self._current.append(data)

    def close(self) -> None:
        super().close()
        self._finish_current()


def _parsed_page(path: Path) -> tuple[str, _OrderParamHTMLParser, str]:
    """Read one HTML page and return raw, parsed and normalized visible text."""

    text = _read_text(path)
    parser = _OrderParamHTMLParser()
    parser.feed(text)
    parser.close()
    visible = " ".join(" ".join(parser.visible_text).split())
    return text, parser, visible


def _is_empty_result_page(text: str, parser: _OrderParamHTMLParser, visible: str) -> bool:
    """Recognize the official, successful Method 3 response with zero rows."""

    folded_raw = text.casefold()
    folded_visible = visible.casefold()
    return (
        RESULT_MARKER.casefold() in folded_raw
        and EMPTY_RESULT_MARKER.casefold() in folded_visible
        and TRY_AGAIN_MARKER.casefold() in folded_visible
        and not parser.rows
    )


def parse_result_page(path: Path) -> list[SubgroupIdentity]:
    """Return candidates from a genuine Method 3 result page.

    The official site has two valid terminal responses: a table containing one
    or more ``orderparam`` rows, or an explicit ``There are no subgroups`` page.
    The latter is an authoritative empty result, not a missing/search page.
    """

    text, parser, visible = _parsed_page(path)
    if _is_empty_result_page(text, parser, visible):
        return []
    if RESULT_MARKER.casefold() not in text.casefold() or not parser.rows:
        raise ValueError("page is not a Method 3 result table")
    return [parse_subgroup_identity(row) for row in parser.rows]


def _active_distortion_types(parser: _OrderParamHTMLParser) -> frozenset[str]:
    result: set[str] = set()
    for values in parser.inputs:
        name = (values.get("name") or "").casefold()
        input_type = (values.get("type") or "").casefold()
        active = input_type == "hidden" or "checked" in values
        if not active or (values.get("value") or "true").casefold() not in {"1", "true", "yes", "on"}:
            continue
        if name == "includestrain":
            result.add("strain")
            continue
        for prefix, distortion_type in _DISTORTION_INPUT_PREFIXES.items():
            if name.startswith(prefix) and not name.endswith("none"):
                result.add(distortion_type)
                break
    return frozenset(result)


def _parse_context_basis(text: str) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    rows = tuple(_parse_vector(row) for row in _VECTOR_RE.findall(text))
    if len(rows) != 3:
        raise ValueError(f"expected three basis vectors in result context, got {len(rows)}")
    return tuple(tuple(row) for row in rows)  # type: ignore[return-value]


def parse_result_page_context(path: Path) -> ResultPageContext:
    """Parse the official page's echoed Method 3 inputs.

    ``centering`` is the effective centering printed by ISODISTORT: ``P`` for
    ``no centering`` and otherwise one of A/B/C/I/F/R.  For a manifest's
    ``Default`` (``d``) request this is later compared with the conventional
    centering letter of the echoed target Hermann--Mauguin symbol.
    """

    _text, parser, visible = _parsed_page(path)
    parent_numbers = {int(value) for value in _PARENT_SPACE_GROUP_RE.findall(visible)}
    for values in parser.inputs:
        if (values.get("name") or "").casefold() != "spacegroup":
            continue
        match = re.match(r"\s*(\d+)\b", values.get("value") or "")
        if match:
            parent_numbers.add(int(match.group(1)))
    if len(parent_numbers) != 1:
        raise ValueError(
            "result page must identify exactly one parent space-group number; "
            f"found {sorted(parent_numbers)}"
        )

    context_match = _RESULT_CONTEXT_RE.search(visible)
    if context_match is None:
        raise ValueError("result page is missing the Method 3 subgroup/basis/centering context")
    raw_centering = context_match.group("center").upper()
    centering = "P" if raw_centering == "NO" else raw_centering
    return ResultPageContext(
        parent_space_group_number=next(iter(parent_numbers)),
        target_space_group_number=int(context_match.group("number")),
        target_space_group_symbol=" ".join(context_match.group("symbol").split()),
        basis=_parse_context_basis(context_match.group("basis")),
        centering=centering,
        distortion_types=_active_distortion_types(parser),
    )


def parse_subgroup_cif(path: Path) -> SubgroupIdentity:
    """Read the first ``# Subgroup:`` identity from a downloaded CIF."""

    for line in _read_text(path).splitlines():
        if re.match(r"^\s*#\s*Subgroup\s*:", line, flags=re.IGNORECASE):
            return parse_subgroup_identity(line)
    raise ValueError(f"missing '# Subgroup:' header in {path}")


def parse_details_page_identity(path: Path) -> SubgroupIdentity:
    """Read the exact subgroup identity printed by a Complete modes page.

    Browser downloads are sometimes placed one directory above or below the
    corresponding CIF.  Pairing such a page by filename or download time is
    not scientific evidence; the page itself must print the same complete
    ``(SG, basis, origin, s, i)`` identity as the candidate CIF.
    """

    _text, parser, _visible = _parsed_page(path)
    identities: list[SubgroupIdentity] = []
    for chunk in parser.visible_text:
        candidate = " ".join(chunk.split())
        if not candidate:
            continue
        try:
            identity = parse_subgroup_identity(candidate)
        except ValueError:
            continue
        if identity not in identities:
            identities.append(identity)
    if len(identities) != 1:
        raise ValueError(
            "Complete modes page must print exactly one subgroup identity; "
            f"found {len(identities)}"
        )
    return identities[0]


def _classify_html(path: Path) -> str:
    text, parser, visible = _parsed_page(path)
    folded = text.casefold()
    if RESULT_MARKER.casefold() in folded and parser.rows:
        return "result_table"
    if _is_empty_result_page(text, parser, visible):
        return "empty_result"
    if "method 3: search over arbitrary" in folded or "isodistort: search" in folded:
        return "search_input"
    if "complete modes details" in path.name.casefold() or "complete mode details" in folded:
        return "complete_modes_details"
    return "other"


@dataclass
class AuditIssue:
    severity: str
    code: str
    message: str
    path: str | None = None


@dataclass
class CandidateAudit:
    directory: str
    identity: dict[str, Any] | None
    core_files: dict[str, str]
    issues: list[AuditIssue] = field(default_factory=list)


@dataclass
class CaseAudit:
    case: str
    status: str
    html_pages: dict[str, str]
    result_pages: list[str]
    result_page_contexts: dict[str, dict[str, Any]]
    official_candidates: list[dict[str, Any]]
    downloaded_candidates: list[dict[str, Any]]
    missing_downloads: list[dict[str, Any]]
    unexpected_downloads: list[dict[str, Any]]
    candidates: list[CandidateAudit]
    issues: list[AuditIssue] = field(default_factory=list)


@dataclass
class Method3Audit:
    root: str
    cases: list[CaseAudit]
    issues: list[AuditIssue] = field(default_factory=list)

    @property
    def error_count(self) -> int:
        return sum(issue.severity == "error" for issue in self.issues) + sum(
            issue.severity == "error"
            for case in self.cases
            for issue in case.issues + [item for candidate in case.candidates for item in candidate.issues]
        )

    @property
    def warning_count(self) -> int:
        return sum(issue.severity == "warning" for issue in self.issues) + sum(
            issue.severity == "warning"
            for case in self.cases
            for issue in case.issues + [item for candidate in case.candidates for item in candidate.issues]
        )

    def to_json(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["summary"] = {
            "case_count": len(self.cases),
            "verified_case_count": sum(case.status == "verified" for case in self.cases),
            "unverified_case_count": sum(case.status != "verified" for case in self.cases),
            "official_candidate_count": sum(len(case.official_candidates) for case in self.cases),
            "downloaded_candidate_count": sum(len(case.downloaded_candidates) for case in self.cases),
            "error_count": self.error_count,
            "warning_count": self.warning_count,
        }
        return payload


def _counter_items(counter: Counter[SubgroupIdentity]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for identity, count in sorted(counter.items(), key=lambda item: item[0].text):
        result.extend(identity.to_json() for _ in range(count))
    return result


def _manifest_basis(value: Any) -> tuple[tuple[Fraction, Fraction, Fraction], ...]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError("supercell_basis must contain exactly three rows")
    rows: list[tuple[Fraction, Fraction, Fraction]] = []
    for row in value:
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError("each supercell_basis row must contain exactly three values")
        rows.append(tuple(Fraction(str(item)) for item in row))  # type: ignore[arg-type]
    return tuple(rows)


def _identity_basis(basis: tuple[tuple[Fraction, Fraction, Fraction], ...]) -> bool:
    return basis == (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )


def manifest_parent_space_group_number(parent: Mapping[str, Any]) -> int:
    """Return the parent SG encoded by one manifest parent record."""

    explicit = parent.get("parent_space_group_type") or parent.get("parent_space_group_number")
    if explicit is not None:
        return int(explicit)

    # The frozen project manifest predates an explicit parent-SG field.  Its
    # totally symmetric Gamma query is nevertheless an exact, scientific
    # encoding of the parent group: GM1+ P1 on the identity lattice preserves
    # the full parent symmetry.  Prefer an explicit field in new manifests,
    # while retaining compatibility with that existing batch.
    candidates: list[int] = []
    for case in parent.get("cases", []):
        expected = case.get("expected_candidate", {})
        irrep = str(expected.get("irrep", "")).replace("Γ", "GM").upper()
        try:
            basis = _manifest_basis(case.get("supercell_basis"))
        except (TypeError, ValueError):
            continue
        if irrep == "GM1+" and str(expected.get("opd", "")).upper() == "P1" and _identity_basis(basis):
            candidates.append(int(case["space_group_type"]))
    unique = sorted(set(candidates))
    if len(unique) != 1:
        raise ValueError(
            "manifest parent must define parent_space_group_type or contain exactly one "
            f"identity-lattice GM1+ P1 query; found {unique}"
        )
    return unique[0]


def manifest_case_expectation(
    *,
    parent_cif: str,
    parent_space_group_number: int,
    case: Mapping[str, Any],
) -> Method3CaseExpectation:
    """Build a page-check expectation from one Method 3 manifest case."""

    case_name = str(case.get("readable_case_folder", ""))
    if not case_name:
        raise ValueError(f"{parent_cif}: manifest case has no readable_case_folder")
    return Method3CaseExpectation(
        case=case_name,
        parent_cif=parent_cif,
        parent_space_group_number=parent_space_group_number,
        target_space_group_number=int(case["space_group_type"]),
        basis=_manifest_basis(case.get("supercell_basis")),
        centering=str(case.get("direct_sublattice_centering", "")).upper(),
        distortion_types=frozenset(str(item).casefold() for item in case.get("distortion_types", [])),
    )


def _read_manifest(manifest: Path | Mapping[str, Any]) -> Mapping[str, Any]:
    if isinstance(manifest, Mapping):
        return manifest
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping):
        raise ValueError("Method 3 manifest root must be a JSON object")
    return payload


def _select_manifest_parent(
    payload: Mapping[str, Any],
    root: Path,
    case_names: set[str],
    parent_cif: str | None,
) -> Mapping[str, Any]:
    parents = payload.get("parents")
    if not isinstance(parents, list) or not parents:
        raise ValueError("Method 3 manifest has no parents")
    records = [item for item in parents if isinstance(item, Mapping)]
    if parent_cif is not None:
        matches = [item for item in records if str(item.get("parent_cif", "")).casefold() == parent_cif.casefold()]
    else:
        root_parts = {part.casefold() for part in root.parts}
        matches = [item for item in records if str(item.get("parent_cif", "")).casefold() in root_parts]
        if not matches:
            scores: list[tuple[int, Mapping[str, Any]]] = []
            for item in records:
                expected_names = {
                    str(case.get("readable_case_folder", ""))
                    for case in item.get("cases", [])
                    if isinstance(case, Mapping)
                }
                scores.append((len(case_names & expected_names), item))
            best = max((score for score, _item in scores), default=0)
            matches = [item for score, item in scores if score == best and score > 0]
        if not matches and len(records) == 1:
            matches = records
    if len(matches) != 1:
        names = [str(item.get("parent_cif", "<unnamed>")) for item in matches]
        raise ValueError(f"could not select exactly one manifest parent for {root}; matches={names}")
    return matches[0]


def _manifest_expectations(
    manifest: Path | Mapping[str, Any],
    root: Path,
    case_names: set[str],
    parent_cif: str | None,
) -> dict[str, Method3CaseExpectation]:
    parent = _select_manifest_parent(_read_manifest(manifest), root, case_names, parent_cif)
    parent_name = str(parent.get("parent_cif", ""))
    parent_sg = manifest_parent_space_group_number(parent)
    expectations: dict[str, Method3CaseExpectation] = {}
    for raw_case in parent.get("cases", []):
        if not isinstance(raw_case, Mapping):
            raise ValueError(f"{parent_name}: manifest case is not an object")
        case_name = str(raw_case.get("readable_case_folder", ""))
        if not case_name:
            raise ValueError(f"{parent_name}: manifest case has no readable_case_folder")
        if case_name in expectations:
            raise ValueError(f"{parent_name}: duplicate manifest case folder {case_name!r}")
        expectations[case_name] = manifest_case_expectation(
            parent_cif=parent_name,
            parent_space_group_number=parent_sg,
            case=raw_case,
        )
    return expectations


def _default_centering(symbol: str) -> str:
    match = re.match(r"\s*([PABCIFR])", symbol.upper())
    if match is None:
        raise ValueError(f"cannot infer Default centering from target symbol {symbol!r}")
    return match.group(1)


def _context_issues(
    path: Path,
    context: ResultPageContext,
    expected: Method3CaseExpectation,
) -> list[AuditIssue]:
    issues: list[AuditIssue] = []
    checks = [
        (
            context.parent_space_group_number == expected.parent_space_group_number,
            "result_parent_space_group_mismatch",
            f"parent SG is {context.parent_space_group_number}; expected {expected.parent_space_group_number} "
            f"for {expected.parent_cif}",
        ),
        (
            context.target_space_group_number == expected.target_space_group_number,
            "result_target_space_group_mismatch",
            f"target SG is {context.target_space_group_number}; expected {expected.target_space_group_number}",
        ),
        (
            context.basis == expected.basis,
            "result_basis_mismatch",
            "input basis is "
            + ",".join(_vector_text(row) for row in context.basis)
            + "; expected "
            + ",".join(_vector_text(row) for row in expected.basis),
        ),
        (
            context.distortion_types == expected.distortion_types,
            "result_distortion_types_mismatch",
            f"distortion types are {sorted(context.distortion_types)}; "
            f"expected {sorted(expected.distortion_types)}",
        ),
    ]
    try:
        expected_centering = (
            _default_centering(context.target_space_group_symbol)
            if expected.centering == "D"
            else expected.centering
        )
        checks.append(
            (
                context.centering == expected_centering,
                "result_centering_mismatch",
                f"effective centering is {context.centering}; expected {expected_centering} "
                f"from manifest value {expected.centering}",
            )
        )
    except ValueError as exc:
        issues.append(AuditIssue("error", "invalid_result_context", str(exc), str(path)))
    issues.extend(AuditIssue("error", code, message, str(path)) for passed, code, message in checks if not passed)
    return issues


def _candidate_core_files(
    directory: Path,
    *,
    matched_details_page: Path | None = None,
) -> tuple[dict[str, str], list[AuditIssue]]:
    files: dict[str, str] = {}
    issues: list[AuditIssue] = []
    for name in CORE_FILENAMES:
        path = directory / name
        if not path.is_file():
            issues.append(AuditIssue("error", "missing_core_file", f"missing {name}", str(path)))
        elif path.stat().st_size == 0:
            issues.append(AuditIssue("error", "empty_core_file", f"empty {name}", str(path)))
        else:
            files[name] = str(path)

    details = [path for path in directory.glob("*.htm*") if "complete modes details" in path.name.casefold()]
    if not details and matched_details_page is not None:
        if matched_details_page.stat().st_size > 0:
            files["complete_modes_details"] = str(matched_details_page)
            issues.append(
                AuditIssue(
                    "warning",
                    "noncanonical_details_page_location",
                    "Complete modes details was paired by its exact subgroup "
                    "identity but is stored outside the candidate directory",
                    str(matched_details_page),
                )
            )
        else:
            issues.append(
                AuditIssue(
                    "error",
                    "empty_core_file",
                    "Complete modes details HTML is empty",
                    str(matched_details_page),
                )
            )
    elif not details:
        issues.append(
            AuditIssue(
                "error",
                "missing_core_file",
                "missing Complete modes details HTML",
                str(directory),
            )
        )
    else:
        nonempty = [path for path in details if path.stat().st_size > 0]
        if not nonempty:
            issues.append(
                AuditIssue(
                    "error",
                    "empty_core_file",
                    "Complete modes details HTML is empty",
                    str(details[0]),
                )
            )
        else:
            files["complete_modes_details"] = str(nonempty[0])
        if len(details) > 1:
            issues.append(
                AuditIssue(
                    "warning",
                    "multiple_details_pages",
                    f"found {len(details)} Complete modes details HTML files in one candidate directory",
                    str(directory),
                )
            )
    return files, issues


def audit_case(
    case_dir: Path,
    expectation: Method3CaseExpectation | None = None,
) -> CaseAudit:
    """Audit one Method 3 case, optionally against its manifest query."""

    all_html = sorted(case_dir.rglob("*.htm*"))
    html_pages = {
        str(path.relative_to(case_dir)): _classify_html(path)
        for path in all_html
    }
    result_paths = [
        case_dir / relative
        for relative, kind in html_pages.items()
        if kind in {"result_table", "empty_result"}
    ]
    issues: list[AuditIssue] = []
    official: list[SubgroupIdentity] = []
    result_contexts: dict[str, dict[str, Any]] = {}

    if not result_paths:
        search_pages = [name for name, kind in html_pages.items() if kind == "search_input"]
        if search_pages:
            issues.append(
                AuditIssue(
                    "error",
                    "search_page_saved_instead_of_result",
                    "saved Method 3 search/input page instead of the candidate result table: "
                    + ", ".join(search_pages),
                    str(case_dir),
                )
            )
        else:
            issues.append(
                AuditIssue(
                    "error",
                    "missing_result_page",
                    "no root-level HTML is a candidate table or an explicit zero-candidate result",
                    str(case_dir),
                )
            )
    else:
        if len(result_paths) > 1:
            issues.append(
                AuditIssue(
                    "error",
                    "multiple_result_pages",
                    f"found {len(result_paths)} result-table pages; expected exactly one",
                    str(case_dir),
                )
            )
        for path in result_paths:
            if path.parent != case_dir:
                issues.append(
                    AuditIssue(
                        "warning",
                        "noncanonical_result_page_location",
                        "the unique result table is stored below the case root; "
                        "its query context and complete candidate inventory are "
                        "still audited from page content",
                        str(path),
                    )
                )
            try:
                official.extend(parse_result_page(path))
            except ValueError as exc:
                issues.append(AuditIssue("error", "invalid_result_row", str(exc), str(path)))
            if expectation is not None:
                try:
                    context = parse_result_page_context(path)
                except ValueError as exc:
                    issues.append(AuditIssue("error", "invalid_result_context", str(exc), str(path)))
                else:
                    result_contexts[str(path)] = context.to_json()
                    issues.extend(_context_issues(path, context, expectation))

    official_counter = Counter(official)
    for identity, count in official_counter.items():
        if count > 1:
            issues.append(
                AuditIssue(
                    "warning",
                    "duplicate_result_row",
                    f"result table repeats {identity.text} {count} times",
                    str(case_dir),
                )
            )

    details_by_identity: dict[SubgroupIdentity, list[Path]] = {}
    for path in all_html:
        if html_pages[str(path.relative_to(case_dir))] != "complete_modes_details":
            continue
        try:
            identity = parse_details_page_identity(path)
        except ValueError:
            continue
        details_by_identity.setdefault(identity, []).append(path)

    candidates: list[CandidateAudit] = []
    downloaded: list[SubgroupIdentity] = []
    for cif_path in sorted(case_dir.rglob("subgroup.cif")):
        candidate_issues: list[AuditIssue] = []
        identity: SubgroupIdentity | None = None
        try:
            identity = parse_subgroup_cif(cif_path)
            downloaded.append(identity)
        except ValueError as exc:
            candidate_issues.append(AuditIssue("error", "invalid_subgroup_cif", str(exc), str(cif_path)))
        matched_details_page: Path | None = None
        local_details = [
            path for path in cif_path.parent.glob("*.htm*")
            if "complete modes details" in path.name.casefold()
        ]
        if identity is not None:
            for details_path in local_details:
                try:
                    details_identity = parse_details_page_identity(details_path)
                except ValueError as exc:
                    candidate_issues.append(
                        AuditIssue(
                            "error",
                            "invalid_details_page_identity",
                            str(exc),
                            str(details_path),
                        )
                    )
                else:
                    if details_identity != identity:
                        candidate_issues.append(
                            AuditIssue(
                                "error",
                                "details_page_identity_mismatch",
                                "Complete modes page identity does not match its "
                                f"candidate CIF: {details_identity.text} != {identity.text}",
                                str(details_path),
                            )
                        )
        if identity is not None and not local_details:
            identity_pages = details_by_identity.get(identity, [])
            if len(identity_pages) == 1:
                matched_details_page = identity_pages[0]
            elif len(identity_pages) > 1:
                candidate_issues.append(
                    AuditIssue(
                        "error",
                        "ambiguous_details_page_identity",
                        f"found {len(identity_pages)} Complete modes pages for "
                        f"the same exact candidate identity: {identity.text}",
                        str(case_dir),
                    )
                )
        core_files, core_issues = _candidate_core_files(
            cif_path.parent,
            matched_details_page=matched_details_page,
        )
        candidate_issues.extend(core_issues)
        candidates.append(
            CandidateAudit(
                directory=str(cif_path.parent),
                identity=identity.to_json() if identity else None,
                core_files=core_files,
                issues=candidate_issues,
            )
        )

    if not candidates and official:
        issues.append(AuditIssue("error", "no_downloaded_candidates", "no subgroup.cif was found", str(case_dir)))

    downloaded_counter = Counter(downloaded)
    for identity, count in downloaded_counter.items():
        if count > 1:
            issues.append(
                AuditIssue(
                    "error",
                    "duplicate_downloaded_identity",
                    f"downloaded identity appears {count} times: {identity.text}",
                    str(case_dir),
                )
            )

    missing: Counter[SubgroupIdentity] = Counter()
    unexpected: Counter[SubgroupIdentity] = Counter()
    if result_paths:
        missing = official_counter - downloaded_counter
        unexpected = downloaded_counter - official_counter
        for identity, count in missing.items():
            issues.append(
                AuditIssue(
                    "error",
                    "missing_candidate_download",
                    f"missing {count} download(s) for {identity.text}",
                    str(case_dir),
                )
            )
        for identity, count in unexpected.items():
            issues.append(
                AuditIssue(
                    "error",
                    "unexpected_candidate_download",
                    f"downloaded {count} candidate(s) absent from the result table: {identity.text}",
                    str(case_dir),
                )
            )

    candidate_errors = any(issue.severity == "error" for candidate in candidates for issue in candidate.issues)
    verified = bool(result_paths) and not missing and not unexpected and not candidate_errors
    if verified and not any(issue.severity == "error" for issue in issues):
        status = "verified"
    elif not result_paths:
        status = "unverified_result_page"
    else:
        status = "failed"

    return CaseAudit(
        case=case_dir.name,
        status=status,
        html_pages=html_pages,
        result_pages=[str(path) for path in result_paths],
        result_page_contexts=result_contexts,
        official_candidates=[identity.to_json() for identity in official],
        downloaded_candidates=[identity.to_json() for identity in downloaded],
        missing_downloads=_counter_items(missing),
        unexpected_downloads=_counter_items(unexpected),
        candidates=candidates,
        issues=issues,
    )


def audit_method3_root(
    root: Path,
    *,
    manifest: Path | Mapping[str, Any] | None = None,
    parent_cif: str | None = None,
) -> Method3Audit:
    """Audit every immediate case directory beneath ``root`` without writing it.

    ``manifest=None`` preserves the generic file/page audit.  Supplying a
    manifest additionally binds every case directory to its intended parent
    SG, target SG, basis, effective centering and distortion-type selection.
    """

    root = root.resolve()
    issues: list[AuditIssue] = []
    if not root.is_dir():
        issues.append(AuditIssue("error", "invalid_root", "Method 3 root is not a directory", str(root)))
        return Method3Audit(str(root), [], issues)

    case_dirs = sorted(path for path in root.iterdir() if path.is_dir())
    if not case_dirs:
        issues.append(AuditIssue("error", "no_cases", "Method 3 root contains no case directories", str(root)))
    expectations: dict[str, Method3CaseExpectation] = {}
    if manifest is not None:
        try:
            expectations = _manifest_expectations(
                manifest,
                root,
                {path.name for path in case_dirs},
                parent_cif,
            )
        except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            issues.append(AuditIssue("error", "invalid_manifest", str(exc), str(manifest)))
        else:
            actual_names = {path.name for path in case_dirs}
            for missing in sorted(expectations.keys() - actual_names):
                issues.append(
                    AuditIssue(
                        "error",
                        "missing_manifest_case",
                        f"manifest case directory is missing: {missing}",
                        str(root / missing),
                    )
                )
            for unexpected in sorted(actual_names - expectations.keys()):
                issues.append(
                    AuditIssue(
                        "error",
                        "unexpected_case_directory",
                        f"case directory is absent from the selected manifest parent: {unexpected}",
                        str(root / unexpected),
                    )
                )
    cases = [audit_case(path, expectations.get(path.name)) for path in case_dirs]
    return Method3Audit(str(root), cases, issues)


def _print_summary(audit: Method3Audit) -> None:
    payload = audit.to_json()
    summary = payload["summary"]
    print(
        "Method 3 download audit: "
        f"{summary['verified_case_count']}/{summary['case_count']} cases verified, "
        f"{summary['official_candidate_count']} official rows, "
        f"{summary['downloaded_candidate_count']} downloaded candidates, "
        f"{summary['error_count']} errors, {summary['warning_count']} warnings"
    )
    for case in audit.cases:
        marker = "OK" if case.status == "verified" else "CHECK"
        print(
            f"[{marker}] {case.case}: {case.status}; "
            f"official={len(case.official_candidates)}, downloaded={len(case.downloaded_candidates)}"
        )
        for issue in case.issues:
            print(f"  {issue.severity.upper()} {issue.code}: {issue.message}")
        for candidate in case.candidates:
            for issue in candidate.issues:
                print(f"  {issue.severity.upper()} {issue.code}: {issue.message}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Read-only audit of official ISODISTORT Method 3 HTML and candidate downloads."
    )
    parser.add_argument("root", type=Path, help="directory whose immediate children are Method 3 cases")
    parser.add_argument("--json-output", type=Path, help="write the complete audit report as UTF-8 JSON")
    manifest_group = parser.add_mutually_exclusive_group()
    manifest_group.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST_PATH,
        help=f"query manifest used to bind pages to cases (default: {DEFAULT_MANIFEST_PATH})",
    )
    manifest_group.add_argument(
        "--no-manifest",
        action="store_true",
        help="run only the generic page/export audit",
    )
    parser.add_argument(
        "--parent-cif",
        help="manifest parent_cif to select explicitly (normally inferred from the root path/case names)",
    )
    args = parser.parse_args(argv)

    manifest = None if args.no_manifest else args.manifest
    audit = audit_method3_root(args.root, manifest=manifest, parent_cif=args.parent_cif)
    _print_summary(audit)
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(audit.to_json(), ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"JSON: {args.json_output.resolve()}")
    return 1 if audit.error_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
