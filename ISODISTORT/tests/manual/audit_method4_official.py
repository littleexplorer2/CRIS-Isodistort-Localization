"""Audit archived official ISODISTORT Method 4 runs without modifying them.

The archive under ``output_compare`` is protected and read-only.  This tool
hashes and parses it, compares the selected transformation and exported
structure with the frozen daughter input, and writes a machine report only to
``ISODISTORT/output/validation``.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import sys
from collections import Counter
from collections.abc import Mapping
from datetime import UTC, datetime
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

WORKSPACE = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = WORKSPACE / "ISODISTORT"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from features.input_cif.cif_io import read_cif  # noqa: E402
from tests.manual.method4_provenance import (  # noqa: E402
    PROVENANCE_KIND,
    PROVENANCE_SCHEMA,
    build_method4_provenance_state,
)
from tests.manual.official_html_resolver import (  # noqa: E402
    OfficialHtmlInventory,
    OfficialHtmlPage,
    scan_official_html,
)

MANIFEST_PATH = PACKAGE_ROOT / "docs" / "manifests" / "method4_download_manifest.json"
LOCAL_REPORT_PATH = PACKAGE_ROOT / "output" / "validation" / "method4_local_validation.json"
OFFICIAL_ROOT = WORKSPACE / "output_compare"
DEFAULT_REPORT = PACKAGE_ROOT / "output" / "validation" / "method4_official_audit.json"
SOURCE_PATHS = [
    Path(__file__).resolve(),
    PACKAGE_ROOT / "tests" / "manual" / "method4_provenance.py",
    PACKAGE_ROOT / "tests" / "manual" / "official_html_resolver.py",
    MANIFEST_PATH,
    LOCAL_REPORT_PATH,
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_signature(paths: list[Path]) -> dict[str, Any]:
    records = []
    combined = hashlib.sha256()
    for path in paths:
        digest = _sha256(path)
        try:
            label = path.resolve().relative_to(WORKSPACE.resolve()).as_posix()
        except ValueError:
            label = str(path.resolve())
        records.append({"path": label, "sha256": digest})
        combined.update(label.encode("utf-8"))
        combined.update(b"\0")
        combined.update(digest.encode("ascii"))
        combined.update(b"\n")
    return {"sha256": combined.hexdigest(), "files": records}


def _html_page_record(
    page: OfficialHtmlPage,
    *,
    root: Path,
) -> dict[str, Any]:
    """Serialize evidence from the resolver's captured read without reopening it."""

    try:
        relative_path = page.path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        relative_path = str(page.path.resolve())
    return {
        "path": str(page.path.resolve()),
        "relative_path": relative_path,
        "role": page.role,
        "title": page.title,
        "sha256": page.content_sha256,
        "classification_evidence": list(page.evidence),
        "captured_character_count": len(page.raw_html),
    }


def _direct_html_pages(
    inventory: OfficialHtmlInventory,
    directory: Path,
) -> tuple[OfficialHtmlPage, ...]:
    """Return only pages physically stored in ``directory``, never descendants."""

    resolved = directory.resolve()
    return tuple(
        page for page in inventory.pages
        if page.path.parent.resolve() == resolved
    )


def _resolve_directory_html(
    inventory: OfficialHtmlInventory,
    directory: Path,
    *,
    require_details: bool,
    require_basis: bool = False,
) -> tuple[dict[str, Any], dict[str, OfficialHtmlPage]]:
    """Resolve the scientific HTML roles in one directory, fail-closed.

    The shared resolver classifies content without consulting basenames.  This
    layer merely applies Method 4's path-specific cardinality contract: result
    is always required, details is required for accepted decompositions, and a
    basis page may be optional on paths with independent control evidence.  An
    optional role still rejects duplicates.
    """

    pages = _direct_html_pages(inventory, directory)
    selected: dict[str, OfficialHtmlPage] = {}
    issues: list[str] = []
    role_contract = (
        ("result_page", "distortion_result", True),
        ("details_page", "complete_modes_details", require_details),
        ("basis_page", "distorted_structure_basis", require_basis),
    )
    for key, role, required in role_contract:
        matches = tuple(page for page in pages if page.role == role)
        if len(matches) == 1:
            selected[key] = matches[0]
        elif len(matches) > 1:
            issues.append(
                f"expected at most one {role} HTML page in {directory}, "
                f"found {len(matches)}"
            )
        elif required:
            issues.append(
                f"missing required {role} HTML page in {directory}"
            )
    record = {
        "directory": str(directory.resolve()),
        "scope": "current_directory_only",
        "inspected_pages": [
            _html_page_record(page, root=directory) for page in pages
        ],
        "selected_pages": {
            key: _html_page_record(page, root=directory)
            for key, page in selected.items()
        },
        "unknown_pages": [
            _html_page_record(page, root=directory)
            for page in pages
            if page.role == "unknown"
        ],
        "issues": issues,
    }
    return record, selected


def _capture_selected_html_evidence(
    manifest: Mapping[str, Any],
    *,
    parent_filter: str | None,
    case_filter: str | None,
) -> tuple[
    dict[str, Any],
    dict[tuple[str, str], OfficialHtmlInventory | None],
    dict[tuple[str, str], str],
]:
    """Capture every HTML page recursively for the selected scientific cases."""

    records: list[dict[str, Any]] = []
    inventories: dict[tuple[str, str], OfficialHtmlInventory | None] = {}
    errors: dict[tuple[str, str], str] = {}
    for parent_name, parent in manifest["parents"].items():
        if parent_filter is not None and parent_name != parent_filter:
            continue
        for case_id in parent["cases"]:
            if case_filter is not None and case_id != case_filter:
                continue
            key = (parent_name, case_id)
            readable_folder = manifest["readable_case_folders"][case_id]
            directory = (
                OFFICIAL_ROOT / parent_name / "官网" / "Method4" / readable_folder
            )
            record: dict[str, Any] = {
                "parent": parent_name,
                "case_id": case_id,
                "directory": str(directory.resolve()),
                "exists": directory.is_dir(),
                "recursive": True,
                "signature": None,
                "pages": [],
            }
            inventory: OfficialHtmlInventory | None = None
            if directory.is_dir():
                try:
                    inventory = scan_official_html(directory, recursive=True)
                except (OSError, ValueError) as exc:
                    error = f"could not capture official HTML inventory: {exc}"
                    record["scan_error"] = error
                    errors[key] = error
                else:
                    record["signature"] = inventory.signature
                    record["pages"] = [
                        _html_page_record(page, root=directory)
                        for page in inventory.pages
                    ]
            inventories[key] = inventory
            records.append(record)

    canonical_records = [
        {
            "parent": record["parent"],
            "case_id": record["case_id"],
            "exists": record["exists"],
            "signature": record["signature"],
            "scan_error": record.get("scan_error"),
        }
        for record in records
    ]
    canonical = json.dumps(
        canonical_records,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    )
    evidence = {
        "schema": 1,
        "sha256": hashlib.sha256(canonical.encode("ascii")).hexdigest(),
        "selected_case_count": len(records),
        "scan_error_count": len(errors),
        "case_inventories": records,
    }
    return evidence, inventories, errors


def _selection_filter_gate(
    manifest: Mapping[str, Any],
    *,
    parent_filter: str | None,
    case_filter: str | None,
) -> dict[str, Any]:
    """Validate CLI filters and describe the exact selected case set.

    An explicitly supplied empty filter is not equivalent to omitting the
    filter. Unknown and empty values are kept as an inconclusive selection
    error so ``main`` can write a machine-readable report before returning 2.
    """

    parents = manifest.get("parents")
    if not isinstance(parents, Mapping):
        parents = {}
    issues: list[str] = []
    if parent_filter is not None:
        if not parent_filter.strip():
            issues.append("--parent must not be empty")
        elif parent_filter not in parents:
            issues.append(f"unknown --parent value: {parent_filter!r}")

    all_case_ids = {
        str(case_id)
        for parent in parents.values()
        if isinstance(parent, Mapping) and isinstance(parent.get("cases"), Mapping)
        for case_id in parent["cases"]
    }
    if case_filter is not None:
        if not case_filter.strip():
            issues.append("--case-id must not be empty")
        elif case_filter not in all_case_ids:
            issues.append(f"unknown --case-id value: {case_filter!r}")

    selected_cases = [
        {"parent": str(parent_name), "case_id": str(case_id)}
        for parent_name, parent in parents.items()
        if isinstance(parent, Mapping)
        and (parent_filter is None or parent_name == parent_filter)
        and isinstance(parent.get("cases"), Mapping)
        for case_id in parent["cases"]
        if case_filter is None or case_id == case_filter
    ]
    if not selected_cases:
        issues.append("the requested filters selected no Method 4 cases")
    return {
        "status": "pass" if not issues else "inconclusive",
        "issues": issues,
        "requested": {"parent": parent_filter, "case_id": case_filter},
        "selected_cases": selected_cases,
        "selected_case_count": len(selected_cases),
    }


def _local_report_provenance_gate(
    local: Mapping[str, Any],
    manifest: Mapping[str, Any],
    manifest_path: Path,
) -> dict[str, Any]:
    """Independently reproduce and verify the local Method 4 evidence state."""

    issues: list[str] = []
    schema = local.get("schema")
    if not isinstance(schema, int) or schema < 3:
        issues.append("local report schema predates layered Method 4 evidence (schema 3)")
    if local.get("kind") != "method4_local_inverse_consistency":
        issues.append("local report kind is not Method 4 inverse consistency")
    if local.get("provenance_failure_count") != 0:
        issues.append("local report records a provenance failure")

    raw_selection = local.get("selection")
    if isinstance(raw_selection, Mapping):
        selection = {
            "parent": raw_selection.get("parent"),
            "context": raw_selection.get("context"),
            "case_id": raw_selection.get("case_id"),
        }
        if dict(raw_selection) != selection:
            issues.append("local report selection is malformed")
    else:
        selection = {"parent": None, "context": None, "case_id": None}
        issues.append("local report selection is missing")

    recomputed = build_method4_provenance_state(manifest_path, selection)
    current_manifest_sha256 = _sha256(manifest_path)
    if local.get("manifest_sha256") != current_manifest_sha256:
        issues.append("local report manifest SHA-256 does not match the current manifest")
    if local.get("manifest_source_signature") != manifest.get("signature"):
        issues.append("local report manifest source signature does not match the manifest")
    if local.get("tolerances") != manifest.get("local_validation_tolerances"):
        issues.append("local report tolerances do not match the current manifest")

    reported_case_count = sum(
        len(parent.get("cases", {}))
        for parent in local.get("parents", {}).values()
        if isinstance(parent, Mapping) and isinstance(parent.get("cases"), Mapping)
    ) if isinstance(local.get("parents"), Mapping) else 0
    if local.get("case_count") != reported_case_count:
        issues.append("local report case_count does not match its case records")

    provenance = local.get("provenance")
    if not isinstance(provenance, Mapping):
        issues.append("local report has no reproducible provenance object")
        provenance = {}
    if provenance.get("schema") != PROVENANCE_SCHEMA:
        issues.append("local provenance schema is unsupported")
    if provenance.get("kind") != PROVENANCE_KIND:
        issues.append("local provenance kind is unsupported")
    if provenance.get("stable_during_run") is not True:
        issues.append("local provenance was not stable during validation")
    if provenance.get("changed_surfaces_during_run") != []:
        issues.append("local provenance records changed surfaces during validation")

    recorded_state = provenance.get("state")
    recorded_signature = provenance.get("signature")
    if recorded_signature != recomputed["signature"]:
        issues.append("local provenance signature does not match independent recomputation")
    if recorded_state != recomputed:
        issues.append("local provenance state does not match independent recomputation")
    if recorded_signature != (
        recorded_state.get("signature") if isinstance(recorded_state, Mapping) else None
    ):
        issues.append("local provenance signature and embedded state disagree")

    start_signature = provenance.get("run_start_signature")
    pre_update_signature = provenance.get("run_end_signature_before_manifest_update")
    if start_signature != pre_update_signature:
        issues.append("local provenance changed while Method 4 validation was running")
    manifest_self_update = provenance.get("manifest_self_update")
    if not isinstance(manifest_self_update, Mapping):
        issues.append("local provenance has no manifest self-update record")
        manifest_self_update = {}
    performed_update = manifest_self_update.get("performed") is True
    changed_surfaces = manifest_self_update.get("changed_surfaces")
    unexpected_surfaces = manifest_self_update.get("unexpected_changed_surfaces")
    if unexpected_surfaces != []:
        issues.append("local provenance records unexpected post-validation changes")
    if performed_update:
        if not isinstance(changed_surfaces, list) or any(
            item != "input-manifest" for item in changed_surfaces
        ):
            issues.append("local manifest self-update changed a non-manifest surface")
    elif pre_update_signature != recorded_signature:
        issues.append("local provenance changed after validation without a manifest self-update")

    return {
        "status": "pass" if not issues else "fail",
        "issues": issues,
        "recorded_signature": recorded_signature,
        "recomputed_signature": recomputed["signature"],
        "current_manifest_sha256": current_manifest_sha256,
        "selection": selection,
        "reported_case_count": reported_case_count,
        "recomputed_state": recomputed,
    }


def _html_text(raw: str) -> str:
    """Extract visible text from HTML already captured by the resolver."""

    raw = re.sub(r"(?i)<br\s*/?>", "\n", raw)
    raw = re.sub(r"(?i)</(?:p|pre|div|tr|h\d)>", "\n", raw)
    raw = re.sub(r"(?s)<[^>]+>", "", raw)
    return html.unescape(raw).replace("\xa0", " ")


def _expected_failure_matches(expectation: str, page_text: str) -> bool:
    """Recognize an expected scientific rejection, never a generic server error."""
    patterns = {
        "reject_species_or_stoichiometry": (
            r"types? of atoms?.*do not match",
            r"species.*(?:do not match|mismatch)",
            r"stoichiometr(?:y|ic).*(?:do not match|mismatch)",
        ),
        "reject_incompatible_lattice": (
            r"lattice.*(?:incompatible|not compatible|does not match)",
            r"(?:incompatible|not compatible).*lattice",
        ),
        "reject_unmatched_atom": (
            r"failed to find match",
            r"larger value for dmax",
            r"unable to match.*atom",
        ),
    }
    return any(
        re.search(pattern, page_text, flags=re.IGNORECASE | re.DOTALL)
        for pattern in patterns.get(expectation, ())
    )


def _number(token: str) -> float:
    return float(Fraction(token.strip()))


def _parse_vector(text: str) -> list[float]:
    return [_number(token) for token in text.split(",")]


def _parse_basis(text: str) -> list[list[float]]:
    rows = re.findall(r"\(([^()]*)\)", text)
    if len(rows) != 3:
        raise ValueError(f"Expected three basis vectors, found {rows!r}")
    return [_parse_vector(row) for row in rows]


def _parse_subgroup(text: str) -> dict[str, Any]:
    patterns = [
        r"Subgroup:\s*(\d+)\s+([^,\n]+),\s*basis=\{([^}]+)\},\s*"
        r"origin=\(([^)]+)\),\s*s=(\d+),\s*i=(\d+)",
        r"Subgroup details\s*(\d+)\s+([^,\n]+),\s*basis=\{([^}]+)\},\s*"
        r"origin=\(([^)]+)\),\s*s=(\d+),\s*i=(\d+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return {
                "space_group_number": int(match.group(1)),
                "space_group_symbol": match.group(2).strip(),
                "basis": _parse_basis(match.group(3)),
                "origin": _parse_vector(match.group(4)),
                "s": int(match.group(5)),
                "i": int(match.group(6)),
            }
    raise ValueError("Could not parse official Method 4 subgroup identity")


def _amplitude_rows(details_text: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    displacement_block = details_text.split("Displacive mode amplitudes", 1)[1]
    displacement_block, strain_tail = displacement_block.split(
        "Parent-cell strain mode definitions", 1,
    )
    displacement_rows = []
    for line in displacement_block.splitlines():
        match = re.match(
            r"^\s*(\[[^\n]+?\))\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$",
            line,
        )
        if match and "  all" not in line:
            displacement_rows.append({
                "label": match.group(1).strip(),
                "As": float(match.group(2)),
                "Ap": float(match.group(3)),
                "dmax_angstrom": float(match.group(4)),
            })
    strain_rows = []
    if "Parent-cell strain mode amplitudes" in strain_tail:
        strain_block = strain_tail.split("Parent-cell strain mode amplitudes", 1)[1]
        for line in strain_block.splitlines():
            match = re.match(r"^\s*(\[[^\n]+?\))\s+([-+0-9.Ee]+)\s*$", line)
            if match:
                strain_rows.append({
                    "label": match.group(1).strip(),
                    "amplitude": float(match.group(2)),
                })
    return displacement_rows, strain_rows


def _definition_normfactors(details_text: str) -> dict[str, float]:
    block = details_text.split("Displacive mode definitions", 1)[1]
    block = block.split("Displacive mode amplitudes", 1)[0]
    result = {}
    for line in block.splitlines():
        match = re.match(
            r"^\s*[^\s\[]+(\[[^\n]+?\))\s+normfactor\s*=\s*([-+0-9.Ee]+)\s*$",
            line,
        )
        if match:
            result[match.group(1).strip()] = float(match.group(2))
    return result


def _rounding_half_unit(token: str) -> float:
    mantissa, _, exponent_text = token.lower().partition("e")
    decimals = len(mantissa.rsplit(".", 1)[1]) if "." in mantissa else 0
    exponent = int(exponent_text) if exponent_text else 0
    return 0.5 * 10.0 ** (exponent - decimals)


def _isoviz_applied_strain_with_rounding(
    text: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the six Voigt components actually applied by IsoVIZ.

    The displayed/CIF strain amplitudes are symmetry-mode coordinates.  The
    viewer applies the following six-component mode vector before converting
    engineering shear components to a symmetric 3x3 matrix, so comparing the
    amplitude numbers alone is not a tensor comparison.
    """

    marker = "!strainmodelist"
    if marker not in text:
        raise ValueError("IsoVIZ export lacks a strain-mode list")
    block = text.split(marker, 1)[1]
    block = block.split("!displacivemodelist", 1)[0]
    lines = [line.strip() for line in block.splitlines() if line.strip()]
    applied = np.zeros(6, dtype=float)
    rounding_bound = np.zeros(6, dtype=float)
    index = 0
    parsed = 0
    while index < len(lines):
        header = re.match(
            r"^\d+\s+([-+0-9.Ee]+)\s+[-+0-9.Ee]+\s+\d+\s+\S+\s*$",
            lines[index],
        )
        if not header:
            index += 1
            continue
        if index + 1 >= len(lines):
            raise ValueError("IsoVIZ strain mode is missing its six-component vector")
        vector_tokens = lines[index + 1].split()
        if len(vector_tokens) != 6:
            raise ValueError("IsoVIZ strain mode vector does not have six components")
        amplitude_token = header.group(1)
        amplitude = float(amplitude_token)
        amplitude_error = _rounding_half_unit(amplitude_token)
        vector = np.asarray([float(token) for token in vector_tokens], dtype=float)
        vector_error = np.asarray(
            [_rounding_half_unit(token) for token in vector_tokens], dtype=float
        )
        applied += amplitude * vector
        rounding_bound += (
            np.abs(vector) * amplitude_error
            + abs(amplitude) * vector_error
            + amplitude_error * vector_error
        )
        parsed += 1
        index += 2
    if parsed == 0:
        raise ValueError("IsoVIZ strain-mode list contains no parseable modes")
    return applied, rounding_bound


def _isoviz_applied_strain(text: str) -> np.ndarray:
    return _isoviz_applied_strain_with_rounding(text)[0]


def _canonical_wavevector(text: str) -> str:
    values = [Fraction(token.strip()) for token in text.split(",")]
    if len(values) != 3:
        raise ValueError(f"Expected a three-component wave vector, got {text!r}")
    if all(value == 0 for value in values):
        return ""
    return ",".join(str(value) for value in values)


def _local_mode_identity(label: str) -> tuple[str, str, str] | None:
    parts = label.split("__")
    if len(parts) < 3:
        return None
    ir_match = re.match(r"^(.+)\[([^]]+)\]$", parts[0])
    if ir_match:
        ir = ir_match.group(1)
        wavevector = _canonical_wavevector(ir_match.group(2))
    else:
        ir = parts[0]
        wavevector = ""
    tail_match = re.match(r"^(.+)\(([^()]*)\)$", parts[-1])
    if tail_match:
        site_ir, component = tail_match.groups()
    else:
        site_ir, component = "", parts[-1]
    # At Gamma the local subgroup-restricted basis and the official P1 basis
    # can use different site-irrep names for the same Cartesian column, so the
    # component is the stable identity.  Away from Gamma, retain both k and the
    # site-irrep name: LD1 harmonics and repeated A1_1/A1_2 branches otherwise
    # collide even though they are distinct physical modes.
    ir_key = f"{ir}[{wavevector}]" if wavevector else ir
    component_key = f"{site_ir}({component})" if wavevector else component
    return ir_key, parts[1], component_key


def _official_mode_identity(label: str) -> tuple[str, str, str] | None:
    match = re.match(
        r"^\[([^]]+)\]([^(]+)\([^)]*\)\[[^:\]]+:([^:\]]+):dsp\]"
        r"([^()]+)\(([^)]+)\)$",
        label,
    )
    if not match:
        return None
    wavevector = _canonical_wavevector(match.group(1))
    ir_key = f"{match.group(2)}[{wavevector}]" if wavevector else match.group(2)
    component_key = (
        f"{match.group(4)}({match.group(5)})" if wavevector else match.group(5)
    )
    return ir_key, match.group(3), component_key


def _parameter_mode_subspace_key(
    identity: tuple[str, str, str] | None,
) -> tuple[str, str, str] | None:
    """Return the basis-independent branch identity for a parameter-k mode.

    Equivalent multidimensional irrep bases may rotate components such as
    ``A2u(a)`` and ``A2u(b)``. Their individual amplitudes are basis
    dependent, while the Euclidean norm of a complete, equally normalized
    subspace is invariant. Gamma-mode identities intentionally return None:
    their local identity deliberately omits the official site-irrep name and
    cannot establish a complete official subspace by itself.
    """
    if identity is None or "[" not in identity[0]:
        return None
    component_match = re.fullmatch(r"(.+)\(([^()]*)\)", identity[2])
    if component_match is None:
        return None
    return identity[0], identity[1], component_match.group(1)


def _gamma_site_subspace_key(
    identity: tuple[str, str, str] | None,
) -> tuple[str, str] | None:
    """Return the complete Gamma irrep/site block for basis-invariant checks.

    A repeated site representation need not have the same component names in
    the two programs.  For example, the four local columns ``a,b,c,d`` for the
    NdNiO2 oxygen ``f`` orbit are the same space that official ISODISTORT calls
    ``B3u(a,b) + B2u(a,b)``.  Individual columns are not identifiable from the
    labels alone, but a complete equally normalized block has invariant As/Ap
    norms.  Parameter-k modes deliberately keep their existing, stricter key.
    """
    if identity is None or "[" in identity[0]:
        return None
    return identity[0], identity[1]


def _complete_gamma_site_subspaces(
    local_rows: list[dict[str, Any]],
    official_rows: list[dict[str, Any]],
) -> list[tuple[tuple[str, str], list[dict[str, Any]], list[dict[str, Any]]]]:
    """Return complete ambiguous Gamma site blocks from both mode bases."""
    local_groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    official_groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in local_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _gamma_site_subspace_key(identity)
        if key is not None:
            local_groups.setdefault(key, []).append(row)
    for row in official_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _gamma_site_subspace_key(identity)
        if key is not None:
            official_groups.setdefault(key, []).append(row)

    complete = []
    for key, rows in local_groups.items():
        official = official_groups.get(key, [])
        # Use a block comparison only when one-to-one component labels are
        # genuinely ambiguous and both sides span the entire same-sized block.
        official_identities = [tuple(row["identity"]) for row in official]
        if (
            len(rows) > 1
            and len(rows) == len(official)
            and len(set(official_identities)) < len(official_identities)
        ):
            complete.append((key, rows, official))
    return complete


def _complete_required_amplitude_subspaces(
    rows: list[dict[str, Any]],
    required_labels: set[str] | None,
    official_identities: set[tuple[str, str, str]],
) -> list[list[dict[str, Any]]]:
    """Select only fully covered multidimensional subspaces for norm checks."""
    if not required_labels:
        return []

    local_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    required_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    official_groups: dict[
        tuple[str, str, str], set[tuple[str, str, str]]
    ] = {}
    for row in rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _parameter_mode_subspace_key(identity)
        if key is None:
            continue
        local_groups.setdefault(key, []).append(row)
        if row["local_label"] in required_labels:
            required_groups.setdefault(key, []).append(row)
    for identity in official_identities:
        key = _parameter_mode_subspace_key(identity)
        if key is not None:
            official_groups.setdefault(key, set()).add(identity)

    complete = []
    for key, required_rows in required_groups.items():
        required_identities = {tuple(row["identity"]) for row in required_rows}
        local_identities = {tuple(row["identity"]) for row in local_groups[key]}
        if (
            len(required_rows) > 1
            and required_identities == local_identities
            and required_identities == official_groups.get(key, set())
            and all(row["official_label"] is not None for row in required_rows)
        ):
            complete.append(required_rows)
    return complete


def _parameter_origin_phase_turns(
    subspace_key: tuple[str, str, str],
    origin_shift_parent_fractional: list[float] | tuple[float, ...],
) -> Fraction | None:
    """Return ``k dot origin_shift`` modulo one for a parameter-k branch."""
    match = re.search(r"\[([^]]+)\]$", subspace_key[0])
    if match is None:
        return None
    coordinates = [
        Fraction(value.strip()) for value in match.group(1).split(",")
    ]
    if len(coordinates) != len(origin_shift_parent_fractional):
        return None
    shift = [
        Fraction(str(float(value))).limit_denominator(1_000_000)
        for value in origin_shift_parent_fractional
    ]
    return sum(
        (coordinate * value for coordinate, value in zip(coordinates, shift, strict=True)),
        Fraction(0),
    ) % 1


def _complete_origin_phase_amplitude_subspaces(
    local_rows: list[dict[str, Any]],
    required_labels: set[str] | None,
    official_rows: list[dict[str, Any]],
    origin_shift_parent_fractional: list[float] | tuple[float, ...],
) -> list[
    tuple[
        tuple[str, str, str],
        list[dict[str, Any]],
        list[dict[str, Any]],
        Fraction,
    ]
]:
    """Return complete parameter branches rotated by an official origin shift.

    A parent-lattice origin shift changes a non-Gamma complex mode by the phase
    ``exp(-2*pi*i*k.dot(t))``. The real ``a,b`` components may therefore rotate
    or change sign even when both programs reconstruct the same daughter. The
    Euclidean norm is invariant, but it is safe to use only when every local
    component in the branch is an explicit case target and the full official
    branch is present in Complete modes details.
    """
    if not required_labels:
        return []

    local_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    official_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in local_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _parameter_mode_subspace_key(identity)
        if key is not None:
            local_groups.setdefault(key, []).append(row)
    for row in official_rows:
        identity = tuple(row["identity"]) if row.get("identity") is not None else None
        key = _parameter_mode_subspace_key(identity)
        if key is not None:
            official_groups.setdefault(key, []).append(row)

    complete = []
    for key, rows in local_groups.items():
        official = official_groups.get(key, [])
        local_identities = {tuple(row["identity"]) for row in rows}
        official_identities = {tuple(row["identity"]) for row in official}
        phase_turns = _parameter_origin_phase_turns(
            key, origin_shift_parent_fractional,
        )
        if (
            official
            and phase_turns not in {None, Fraction(0)}
            and all(row["local_label"] in required_labels for row in rows)
            and len(local_identities) == len(rows)
            and len(official_identities) == len(official)
            and local_identities <= official_identities
            and all(row["official_label"] is not None for row in rows)
        ):
            complete.append((key, rows, official, phase_turns))
    return complete


def _species_counts(structure) -> Counter[str]:
    return Counter(site.species_string for site in structure)


def _compare_structures(expected_path: Path, actual_path: Path) -> dict[str, Any]:
    expected = read_cif(expected_path)
    actual = read_cif(actual_path)
    expected_metric = np.asarray(expected.lattice.metric_tensor, dtype=float)
    actual_metric = np.asarray(actual.lattice.metric_tensor, dtype=float)
    metric_residual = float(
        np.linalg.norm(expected_metric - actual_metric)
        / max(float(np.linalg.norm(expected_metric)), 1.0)
    )
    result: dict[str, Any] = {
        "expected_atom_count": len(expected),
        "actual_atom_count": len(actual),
        "expected_species": dict(_species_counts(expected)),
        "actual_species": dict(_species_counts(actual)),
        "lattice_metric_relative_residual": metric_residual,
        "direct_max_species_matched_distance_angstrom": None,
        "max_species_matched_distance_angstrom": None,
        "origin_aligned_root_sum_squared_distance_angstrom": None,
        "origin_aligned_fractional_shift": None,
    }
    if len(expected) != len(actual) or _species_counts(expected) != _species_counts(actual):
        return result
    assignments = [-1] * len(expected)
    direct_max_distance = 0.0
    expected_species = [site.species_string for site in expected]
    actual_species = [site.species_string for site in actual]
    for species in sorted(set(expected_species)):
        left = [index for index, value in enumerate(expected_species) if value == species]
        right = [index for index, value in enumerate(actual_species) if value == species]
        distances = expected.lattice.get_all_distances(
            np.asarray(expected.frac_coords)[left],
            np.asarray(actual.frac_coords)[right],
        )
        rows, columns = linear_sum_assignment(distances)
        for row, column in zip(rows, columns, strict=True):
            assignments[left[int(row)]] = right[int(column)]
            direct_max_distance = max(
                direct_max_distance, float(distances[row, column]),
            )
    result["direct_max_species_matched_distance_angstrom"] = direct_max_distance

    # Official automatic-origin detection can translate the complete P1 atom
    # list while preserving the same physical structure.  Compare every
    # species-compatible anchor translation and retain the globally optimal
    # species-wise Hungarian assignment.  This is a single global origin
    # shift, not an independent per-atom adjustment.
    expected_frac = np.asarray(expected.frac_coords, dtype=float)
    actual_frac = np.asarray(actual.frac_coords, dtype=float)
    best: tuple[float, float, list[int], np.ndarray] | None = None
    for species in sorted(set(expected_species)):
        expected_anchors = [
            index for index, value in enumerate(expected_species) if value == species
        ]
        actual_anchors = [
            index for index, value in enumerate(actual_species) if value == species
        ]
        for expected_anchor in expected_anchors:
            for actual_anchor in actual_anchors:
                shift = expected_frac[expected_anchor] - actual_frac[actual_anchor]
                # Refine the anchor translation to the least-squares global
                # origin shift for the resulting periodic assignment.  An
                # anchor alone forces one rounded CIF coordinate to match
                # exactly and can leave a larger artificial maximum residual.
                for _ in range(5):
                    shifted_actual = actual_frac + shift
                    paired: list[tuple[int, int]] = []
                    for trial_species in sorted(set(expected_species)):
                        left = [
                            index for index, value in enumerate(expected_species)
                            if value == trial_species
                        ]
                        right = [
                            index for index, value in enumerate(actual_species)
                            if value == trial_species
                        ]
                        distances = expected.lattice.get_all_distances(
                            expected_frac[left], shifted_actual[right],
                        )
                        rows, columns = linear_sum_assignment(distances)
                        paired.extend(
                            (left[int(row)], right[int(column)])
                            for row, column in zip(rows, columns, strict=True)
                        )
                    fractional_deltas = np.asarray([
                        expected_frac[left] - shifted_actual[right]
                        for left, right in paired
                    ])
                    fractional_deltas -= np.round(fractional_deltas)
                    correction = np.mean(fractional_deltas, axis=0)
                    shift = shift + correction
                    if float(np.linalg.norm(correction)) < 1e-14:
                        break
                shifted_actual = actual_frac + shift
                trial_assignments = [-1] * len(expected)
                squared_sum = 0.0
                max_distance = 0.0
                for trial_species in sorted(set(expected_species)):
                    left = [
                        index for index, value in enumerate(expected_species)
                        if value == trial_species
                    ]
                    right = [
                        index for index, value in enumerate(actual_species)
                        if value == trial_species
                    ]
                    distances = expected.lattice.get_all_distances(
                        expected_frac[left], shifted_actual[right],
                    )
                    rows, columns = linear_sum_assignment(distances)
                    for row, column in zip(rows, columns, strict=True):
                        value = float(distances[row, column])
                        trial_assignments[left[int(row)]] = right[int(column)]
                        squared_sum += value * value
                        max_distance = max(max_distance, value)
                candidate = (squared_sum, max_distance, trial_assignments, shift)
                if best is None or candidate[:2] < best[:2]:
                    best = candidate
    if best is not None:
        aligned_squared_sum, aligned_max, aligned_assignments, aligned_shift = best
        aligned_shift = aligned_shift - np.round(aligned_shift)
        result["max_species_matched_distance_angstrom"] = aligned_max
        result["origin_aligned_root_sum_squared_distance_angstrom"] = float(
            np.sqrt(aligned_squared_sum)
        )
        result["origin_aligned_fractional_shift"] = aligned_shift.tolist()
        result["assignment"] = aligned_assignments
    else:
        result["max_species_matched_distance_angstrom"] = direct_max_distance
        result["assignment"] = assignments
    return result


def _compare_exported_cif(
    expected_path: Path,
    actual_path: Path,
) -> tuple[dict[str, Any], float]:
    comparison = _compare_structures(expected_path, actual_path)
    cif_text = actual_path.read_text(encoding="utf-8", errors="replace")
    coordinate_rows = re.findall(
        r"(?m)^\s*\S+\s+\S+\s+\d+\s+\S+\s+"
        r"([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+",
        cif_text,
    )
    coordinate_decimals = (
        min(
            len(token.rsplit(".", 1)[1])
            for row in coordinate_rows
            for token in row
        )
        if coordinate_rows else None
    )
    expected_structure = read_cif(expected_path)
    coordinate_rounding_bound = (
        0.5 * 10.0 ** (-coordinate_decimals)
        * float(np.sum(np.linalg.norm(expected_structure.lattice.matrix, axis=1)))
        if coordinate_decimals is not None else 0.0
    )
    distance_tolerance = max(1e-4, 1.05 * coordinate_rounding_bound)
    comparison["official_fractional_coordinate_decimals"] = coordinate_decimals
    comparison["coordinate_rounding_bound_angstrom"] = coordinate_rounding_bound
    comparison["distance_tolerance_angstrom"] = distance_tolerance
    return comparison, distance_tolerance


def _amplitude_comparison_tolerance(
    case: dict[str, Any],
    structure_comparison: dict[str, Any],
    local_case: dict[str, Any] | None = None,
) -> tuple[float | None, dict[str, Any]]:
    """Return a case-derived bound for local/official shared-mode amplitudes.

    Exact inputs are limited by the website's five-decimal amplitude display.
    For a deliberately noisy input, official P1 decomposition can select a
    slightly different global origin and distribute noise into additional
    modes.  Unit column norms alone do not bound amplitude sensitivity: the
    aligned Cartesian RSS difference must be divided by the smallest singular
    value of the actual unit-mode matrix.  Source rounding can perturb that
    matrix, so Weyl's inequality supplies the denominator sigma_min-epsilon.
    The separately reported source bound already contains the local residual
    and matrix perturbation terms; only the structure term is added here.
    """
    display_rounding = 5.1e-6
    structure_bound = 0.0
    noisy = case.get("noise_sigma_fractional") is not None
    if noisy:
        value = structure_comparison.get(
            "origin_aligned_root_sum_squared_distance_angstrom"
        )
        if value is not None:
            structure_bound = float(value)
        else:
            maximum = structure_comparison.get(
                "max_species_matched_distance_angstrom"
            )
            count = structure_comparison.get("expected_atom_count")
            if maximum is not None and count is not None:
                structure_bound = float(maximum) * float(np.sqrt(int(count)))
    source_evidence = (
        (local_case or {}).get("mode_source_quantization_bound", {})
        if isinstance(local_case, Mapping) else {}
    )
    if not isinstance(source_evidence, Mapping):
        source_evidence = {}
    source_available = source_evidence.get("status") == "available"
    source_as = (
        float(source_evidence.get("As_vector_l2_error_bound_angstrom", 0.0))
        if source_available else 0.0
    )
    source_ap = (
        float(source_evidence.get("Ap_vector_l2_error_bound_angstrom", 0.0))
        if source_available else 0.0
    )
    evidence = {
        "status": "available",
        "display_rounding_angstrom": display_rounding,
        "aligned_structure_rss_bound_angstrom": structure_bound,
        "source_mode_quantization_status": source_evidence.get("status"),
        "source_mode_quantization_reason": source_evidence.get("reason"),
        "source_As_vector_l2_bound_angstrom": source_as,
        "source_Ap_vector_l2_bound_angstrom": source_ap,
        "derivation": (
            "official display rounding + aligned-export Cartesian RSS / "
            "(sigma_min(unit_mode_matrix)-source_matrix_error) + "
            "DISPLAY DISTORTION lexical-quantization bound propagated through "
            "the same measured singular value; the source term already includes "
            "the local residual and is added once; Ap divides both propagated "
            "As terms by sqrt(supercell_size)"
        ),
    }
    structure_as = 0.0
    structure_ap = 0.0
    if noisy:
        reason = None
        try:
            if not source_available:
                raise ValueError("unit-mode matrix source evidence is unavailable")
            if value is None and (maximum is None or count is None):
                raise ValueError("aligned structure RSS evidence is unavailable")
            sigma_min = float(source_evidence["unit_mode_matrix_smallest_singular_value"])
            matrix_error = float(source_evidence["unit_mode_matrix_l2_error_bound"])
            source_as = float(source_evidence["As_vector_l2_error_bound_angstrom"])
            source_ap = float(source_evidence["Ap_vector_l2_error_bound_angstrom"])
            supercell_size = float(local_case["metadata"]["supercell_size"])
            if (
                not all(np.isfinite(item) for item in (
                    structure_bound, sigma_min, matrix_error,
                    supercell_size, source_as, source_ap,
                ))
                or structure_bound < 0.0 or matrix_error < 0.0
                or sigma_min <= 0.0 or supercell_size <= 0.0
                or source_as < 0.0 or source_ap < 0.0
            ):
                raise ValueError("amplitude sensitivity evidence must be finite and physically valid")
            denominator = sigma_min - matrix_error
            evidence.update(
                unit_mode_matrix_smallest_singular_value=sigma_min,
                unit_mode_matrix_l2_error_bound=matrix_error,
                amplitude_sensitivity_denominator=denominator,
                supercell_size=supercell_size,
            )
            if denominator <= 0.0:
                raise ValueError("source matrix error reaches the smallest singular value")
            structure_as = structure_bound / denominator
            structure_ap = structure_as / float(np.sqrt(supercell_size))
            source_ap = source_as / float(np.sqrt(supercell_size))
            if not all(np.isfinite(item) for item in (
                structure_as, structure_ap, source_as + structure_as,
                source_ap + structure_ap,
            )):
                raise ValueError("propagated amplitude sensitivity is non-finite")
            evidence["source_As_vector_l2_bound_angstrom"] = source_as
            evidence["source_Ap_vector_l2_bound_angstrom"] = source_ap
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            reason = str(exc)
        if reason is not None:
            evidence.update(
                status="inconclusive",
                reason=reason,
                aligned_structure_As_vector_l2_bound_angstrom=None,
                aligned_structure_Ap_vector_l2_bound_angstrom=None,
                As_total_angstrom=None,
                Ap_total_angstrom=None,
                total_angstrom=None,
            )
            return None, evidence
    as_total = display_rounding + structure_as + source_as
    ap_total = display_rounding + structure_ap + source_ap
    total = max(as_total, ap_total)
    evidence.update(
        aligned_structure_As_vector_l2_bound_angstrom=structure_as,
        aligned_structure_Ap_vector_l2_bound_angstrom=structure_ap,
        As_total_angstrom=as_total,
        Ap_total_angstrom=ap_total,
        total_angstrom=total,
    )
    return total, evidence


def _subspace_amplitude_tolerances(
    evidence: Mapping[str, Any],
    official_dimension: int,
) -> tuple[float, float]:
    """Return As/Ap norm bounds for a displayed vector of given dimension."""

    if evidence.get("status") != "available":
        raise ValueError("amplitude sensitivity evidence is unavailable")
    dimension = max(int(official_dimension), 1)
    display = float(evidence["display_rounding_angstrom"]) * float(
        np.sqrt(dimension)
    )
    return (
        display
        + float(evidence["aligned_structure_As_vector_l2_bound_angstrom"])
        + float(evidence.get("source_As_vector_l2_bound_angstrom", 0.0)),
        display
        + float(evidence["aligned_structure_Ap_vector_l2_bound_angstrom"])
        + float(evidence.get("source_Ap_vector_l2_bound_angstrom", 0.0)),
    )


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _checked(raw: str, name: str, value: str) -> bool:
    pattern = (
        rf"<input\b(?=[^>]*\bname=[\"']{re.escape(name)}[\"'])"
        rf"(?=[^>]*\bvalue=[\"']{re.escape(value)}[\"'])[^>]*\bchecked(?:=[\"'][^\"']*[\"'])?[^>]*>"
    )
    return re.search(pattern, raw, flags=re.IGNORECASE) is not None


def _audit_diagnostic_html_subdirectories(
    case_directory: Path,
    inventory: OfficialHtmlInventory,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Resolve direct HTML evidence in auxiliary diagnostic run directories."""

    diagnostic_directories: set[Path] = {
        path
        for path in case_directory.iterdir()
        if path.is_dir()
        and path.name.casefold() != "auto-origin"
        and not path.name.casefold().endswith("_files")
    }
    for page in inventory.pages:
        try:
            relative = page.path.relative_to(case_directory)
        except ValueError:
            continue
        if len(relative.parts) != 2:
            continue
        child_name = relative.parts[0]
        if child_name.casefold() == "auto-origin" or child_name.casefold().endswith(
            "_files"
        ):
            continue
        diagnostic_directories.add(case_directory / child_name)

    records: list[dict[str, Any]] = []
    issues: list[str] = []
    for directory in sorted(
        diagnostic_directories,
        key=lambda path: (path.name.casefold(), path.name),
    ):
        resolution, _selected = _resolve_directory_html(
            inventory,
            directory,
            require_details=False,
            require_basis=False,
        )
        resolution["status"] = "pass" if not resolution["issues"] else "fail"
        records.append(resolution)
        issues.extend(
            f"diagnostic HTML directory {directory.name!r}: {issue}"
            for issue in resolution["issues"]
        )
    return records, issues


def _audit_g05_auto_origin(
    directory: Path,
    expected_daughter: Path,
    primary_result: dict[str, Any],
    html_inventory: OfficialHtmlInventory | None = None,
    html_inventory_error: str | None = None,
) -> dict[str, Any]:
    """Audit the independent G05 automatic-origin run.

    The explicit-origin run remains the primary case record.  This secondary
    record proves that website origin discovery converges to the same physical
    decomposition; missing screenshots are evidence warnings, not numerical
    failures.
    """
    result: dict[str, Any] = {
        "directory": str(directory.resolve()),
        "status": "not_downloaded",
        "issues": [],
        "warnings": [],
    }
    if not directory.is_dir():
        result["issues"].append("auto-origin directory is missing")
        result["status"] = "fail"
        return result
    if html_inventory_error is not None:
        result["issues"].append(html_inventory_error)
        result["status"] = "fail"
        return result
    if html_inventory is None:
        try:
            html_inventory = scan_official_html(directory, recursive=True)
        except (OSError, ValueError) as exc:
            result["issues"].append(
                f"could not capture auto-origin HTML inventory: {exc}"
            )
            result["status"] = "fail"
            return result
    core = {
        "cif": directory / "subgroup.cif",
        "topas": directory / "topas.str",
        "isoviz": directory / "data.isoviz",
    }
    missing = [name for name, path in core.items() if not path.is_file()]
    html_resolution, selected_pages = _resolve_directory_html(
        html_inventory,
        directory,
        require_details=True,
        require_basis=False,
    )
    result["html_evidence"] = html_resolution
    result["selected_html_pages"] = html_resolution["selected_pages"]
    result["issues"].extend(html_resolution["issues"])
    if missing:
        result["issues"].append(f"missing auto-origin core file(s): {missing}")
    if result["issues"]:
        result["status"] = "fail"
        return result
    files = {name: _file_record(path) for name, path in core.items()}
    files.update(html_resolution["selected_pages"])
    if "basis_page" not in selected_pages:
        result["warnings"].append("auto-origin basis/matching page was not archived")
    else:
        result["input_page_controls"] = {
            "automatic_origin_checked_in_saved_html": _checked(
                selected_pages["basis_page"].raw_html,
                "chooseorigin",
                "true",
            ),
            "nearest_site_checked_in_saved_html": _checked(
                selected_pages["basis_page"].raw_html,
                "trynearest",
                "true",
            ),
            "robust_checked_in_saved_html": _checked(
                selected_pages["basis_page"].raw_html,
                "trynearest",
                "false",
            ),
        }
    result["files"] = files
    control_captures = [
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}
    ]
    result["dynamic_control_captures"] = [
        _file_record(path) for path in control_captures
    ]
    if not control_captures:
        result["warnings"].append(
            "auto-origin selection and completed controls lack a screenshot/PDF"
        )

    result_text = _html_text(selected_pages["result_page"].raw_html)
    details_text = _html_text(selected_pages["details_page"].raw_html)
    result_subgroup = _parse_subgroup(result_text)
    details_subgroup = _parse_subgroup(details_text)
    result["official_subgroup"] = result_subgroup
    result["details_subgroup"] = details_subgroup
    if result_subgroup != details_subgroup:
        result["issues"].append(
            "auto-origin result and details subgroup identities differ"
        )
    if result_subgroup != primary_result.get("official_subgroup"):
        result["issues"].append(
            "auto-origin and explicit-origin subgroup identities differ"
        )

    displacement_rows, strain_rows = _amplitude_rows(details_text)
    result["official_amplitudes"] = {
        "displacive_mode_count": len(displacement_rows),
        "strain_mode_count": len(strain_rows),
        "displacive": displacement_rows,
        "strain": strain_rows,
    }
    primary_amplitudes = primary_result.get("official_amplitudes", {})
    if displacement_rows != primary_amplitudes.get("displacive"):
        result["issues"].append(
            "auto-origin and explicit-origin displacive amplitudes differ"
        )
    if strain_rows != primary_amplitudes.get("strain"):
        result["issues"].append(
            "auto-origin and explicit-origin strain amplitudes differ"
        )

    structure_comparison, structure_tolerance = _compare_exported_cif(
        expected_daughter, core["cif"],
    )
    result["exported_cif_vs_frozen_daughter"] = structure_comparison
    if structure_comparison["expected_atom_count"] != structure_comparison["actual_atom_count"]:
        result["issues"].append("auto-origin CIF atom count differs from frozen daughter")
    if structure_comparison["expected_species"] != structure_comparison["actual_species"]:
        result["issues"].append("auto-origin CIF species differ from frozen daughter")
    if structure_comparison["lattice_metric_relative_residual"] > 1e-6:
        result["issues"].append("auto-origin CIF lattice differs from frozen daughter")
    distance = structure_comparison["max_species_matched_distance_angstrom"]
    if distance is None or distance > structure_tolerance:
        result["issues"].append(
            "auto-origin CIF does not reproduce the frozen daughter within its "
            f"coordinate-precision tolerance {structure_tolerance:.8g} angstrom: {distance}"
        )

    result["status"] = (
        "pass" if not result["issues"] and not result["warnings"]
        else "pass_with_warnings" if not result["issues"]
        else "fail"
    )
    return result


def _audit_case(
    parent_name: str,
    case_id: str,
    case: dict[str, Any],
    parent: dict[str, Any],
    readable_folder: str,
    local_case: dict[str, Any] | None,
    *,
    html_inventory: OfficialHtmlInventory | None = None,
    html_inventory_error: str | None = None,
) -> dict[str, Any]:
    directory = OFFICIAL_ROOT / parent_name / "官网" / "Method4" / readable_folder
    result: dict[str, Any] = {
        "parent": parent_name,
        "case_id": case_id,
        "directory": str(directory.resolve()),
        "status": "not_downloaded",
        "issues": [],
        "warnings": [],
        "inconclusive_reasons": [],
    }
    if not directory.is_dir():
        return result
    if html_inventory_error is not None:
        result.update(status="fail", issues=[html_inventory_error])
        return result
    if html_inventory is None:
        try:
            html_inventory = scan_official_html(directory, recursive=True)
        except (OSError, ValueError) as exc:
            result.update(
                status="fail",
                issues=[f"could not capture official HTML inventory: {exc}"],
            )
            return result

    required = {
        "cif": directory / "subgroup.cif",
        "topas": directory / "topas.str",
        "isoviz": directory / "data.isoviz",
    }
    expected_rejection = str(case.get("expect", "")).startswith("reject_")
    archive_started = bool(html_inventory.pages) or any(
        path.is_file() for path in required.values()
    ) or any(
        path.is_file()
        for path in directory.iterdir()
        if path.suffix.casefold() in {".png", ".jpg", ".jpeg", ".pdf", ".cif"}
    ) or any(
        path.is_dir() and not path.name.casefold().endswith("_files")
        for path in directory.iterdir()
    )
    if not archive_started:
        return result

    html_resolution, selected_pages = _resolve_directory_html(
        html_inventory,
        directory,
        require_details=not expected_rejection,
        require_basis=False,
    )
    result["html_evidence"] = html_resolution
    result["selected_html_pages"] = html_resolution["selected_pages"]
    diagnostic_records, diagnostic_issues = _audit_diagnostic_html_subdirectories(
        directory,
        html_inventory,
    )
    result["diagnostic_html_subdirectories"] = diagnostic_records
    result["issues"].extend(html_resolution["issues"])
    result["issues"].extend(diagnostic_issues)
    if html_resolution["issues"]:
        result["status"] = "fail"
        return result

    if expected_rejection:
        result["files"] = dict(html_resolution["selected_pages"])
        for name, path in required.items():
            if path.is_file():
                result["files"][name] = _file_record(path)

        uploaded_copies = [
            path for path in directory.glob("*.cif")
            if path.name.lower() != "subgroup.cif"
        ]
        result["uploaded_daughter_copies"] = [
            _file_record(path) for path in uploaded_copies
        ]
        expected_daughter = WORKSPACE / case["daughter_cif"]
        actual_hash = _sha256(expected_daughter)
        result["frozen_daughter"] = {
            "path": str(expected_daughter.resolve()),
            "sha256_expected": case["sha256"],
            "sha256_actual": actual_hash,
        }
        if actual_hash != case["sha256"]:
            result["issues"].append("frozen daughter CIF hash mismatch")
        if not uploaded_copies:
            result["warnings"].append(
                "The exact uploaded daughter CIF was not copied into the official case directory"
            )
        elif not any(_sha256(path) == case["sha256"] for path in uploaded_copies):
            result["issues"].append(
                "no archived uploaded daughter CIF matches the frozen SHA-256"
            )

        control_captures = [
            path for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}
        ]
        result["dynamic_control_captures"] = [
            _file_record(path) for path in control_captures
        ]
        basis_page = selected_pages.get("basis_page")
        if basis_page is not None:
            basis_raw = basis_page.raw_html
            result["input_page_controls"] = {
                "automatic_origin_checked_in_saved_html": _checked(
                    basis_raw, "chooseorigin", "false",
                ),
                "nearest_site_checked_in_saved_html": _checked(
                    basis_raw, "trynearest", "true",
                ),
                "robust_checked_in_saved_html": _checked(
                    basis_raw, "trynearest", "false",
                ),
            }
        if not control_captures:
            result["warnings"].append(
                "expected-failure controls/error message lack a screenshot/PDF"
            )

        page_text = _html_text(selected_pages["result_page"].raw_html)
        expectation = str(case["expect"])
        observed = _expected_failure_matches(expectation, page_text)
        success_artifacts = [
            name for name in ("details_page", "cif", "topas", "isoviz")
            if (
                name in selected_pages
                if name == "details_page"
                else required[name].is_file()
            )
        ]
        result["expected_rejection"] = {
            "expectation": expectation,
            "recognized": observed,
            "success_artifacts": success_artifacts,
        }
        if success_artifacts:
            result["issues"].append(
                "official Method 4 accepted a case frozen as an expected rejection; "
                f"success artifacts present: {success_artifacts}"
            )
        if not observed:
            result["issues"].append(
                f"official page does not contain the expected rejection semantics {expectation!r}"
            )
        result["local_case"] = local_case
        if local_case is None or local_case.get("status") != "pass":
            result["issues"].append("matching local expected-rejection case is absent or not passing")
        result["status"] = (
            "pass" if not result["issues"] and not result["warnings"]
            else "pass_with_warnings" if not result["issues"]
            else "fail"
        )
        return result

    missing = [name for name, path in required.items() if not path.is_file()]
    if missing:
        result["issues"].append(f"missing required file(s): {missing}")
        result["status"] = "fail"
        return result
    result["files"] = {
        name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": _sha256(path)}
        for name, path in required.items()
    }
    result["files"].update(html_resolution["selected_pages"])
    uploaded_copies = [
        path for path in directory.glob("*.cif")
        if path.name.lower() not in {"subgroup.cif"}
    ]
    result["uploaded_daughter_copies"] = [
        {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in uploaded_copies
    ]
    if not uploaded_copies:
        result["warnings"].append(
            "The exact uploaded daughter CIF was not copied into the official case directory"
        )

    expected_daughter = WORKSPACE / case["daughter_cif"]
    actual_hash = _sha256(expected_daughter)
    result["frozen_daughter"] = {
        "path": str(expected_daughter.resolve()),
        "sha256_expected": case["sha256"],
        "sha256_actual": actual_hash,
    }
    if actual_hash != case["sha256"]:
        result["issues"].append("frozen daughter CIF hash mismatch")
    if uploaded_copies and not any(_sha256(path) == case["sha256"] for path in uploaded_copies):
        result["issues"].append(
            "no archived uploaded daughter CIF matches the frozen SHA-256"
        )

    basis_page = selected_pages.get("basis_page")
    basis_raw = basis_page.raw_html if basis_page is not None else ""
    basis_text = _html_text(basis_raw) if basis_page is not None else ""
    result_text = _html_text(selected_pages["result_page"].raw_html)
    details_text = _html_text(selected_pages["details_page"].raw_html)
    result_subgroup = _parse_subgroup(result_text)
    details_subgroup = _parse_subgroup(details_text)
    result["official_subgroup"] = result_subgroup
    result["details_subgroup"] = details_subgroup
    if result_subgroup != details_subgroup:
        result["issues"].append("result and Complete modes details subgroup identities differ")

    expected_context = parent["resolved_contexts"][case["context"]]
    expected_basis = np.asarray(expected_context["basis_vectors"], dtype=float)
    official_basis = np.asarray(result_subgroup["basis"], dtype=float)
    if not np.allclose(expected_basis, official_basis, atol=1e-12, rtol=0):
        result["issues"].append(
            f"official basis {official_basis.tolist()} differs from frozen context "
            f"{expected_basis.tolist()}"
        )
    result["input_page_controls"] = {
        "automatic_origin_checked_in_saved_html": _checked(
            basis_raw, "chooseorigin", "false",
        ),
        "nearest_site_checked_in_saved_html": _checked(basis_raw, "trynearest", "true"),
        "robust_checked_in_saved_html": _checked(basis_raw, "trynearest", "false"),
    }
    control_captures = [
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}
    ]
    result["dynamic_control_captures"] = [
        {
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in control_captures
    ]
    if basis_page is None:
        if control_captures:
            result["warnings"].append(
                "basis/matching HTML was not archived; the completed-controls screenshot/PDF "
                "is the surviving input-form evidence"
            )
        else:
            result["issues"].append(
                "neither basis/matching HTML nor a completed-controls screenshot/PDF was archived"
            )
    if "Specify basis as:" in basis_text and not control_captures:
        result["warnings"].append(
            "Saved HTML contains server-default form values, not proof of dynamic manual basis input; "
            "no screenshot/PDF was archived; the accepted basis is instead proven by the "
            "result/details/export artifacts"
        )

    displacement_rows, strain_rows = _amplitude_rows(details_text)
    definition_normfactors = _definition_normfactors(details_text)
    result["official_amplitudes"] = {
        "displacive_mode_count": len(displacement_rows),
        "strain_mode_count": len(strain_rows),
        "displacive": displacement_rows,
        "strain": strain_rows,
    }
    isoviz_text = required["isoviz"].read_text(
        encoding="utf-8", errors="replace"
    )
    try:
        (
            official_applied_strain,
            official_strain_rounding_bound,
        ) = _isoviz_applied_strain_with_rounding(isoviz_text)
    except ValueError as exc:
        official_applied_strain = None
        official_strain_rounding_bound = None
        result["issues"].append(str(exc))
    result["official_applied_strain_voigt_engineering"] = (
        None
        if official_applied_strain is None
        else {
            label: float(value)
            for label, value in zip(
                ("xx", "yy", "zz", "2yz", "2xz", "2xy"),
                official_applied_strain,
                strict=True,
            )
        }
    )
    result["official_applied_strain_rounding_bound"] = (
        None
        if official_strain_rounding_bound is None
        else official_strain_rounding_bound.tolist()
    )
    if not displacement_rows:
        result["issues"].append("no displacive amplitude rows parsed")
    is_zero_case = case_id in {"G01-zero", "P01-zero-supercell"}
    if is_zero_case:
        nonzero = [
            row for row in displacement_rows
            if max(abs(row["As"]), abs(row["Ap"]), abs(row["dmax_angstrom"])) > 5e-6
        ]
        nonzero_strain = [row for row in strain_rows if abs(row["amplitude"]) > 5e-6]
        if nonzero or nonzero_strain:
            result["issues"].append(
                f"zero case has nonzero official amplitudes: displacive={nonzero}, strain={nonzero_strain}"
            )

    cif_text = required["cif"].read_text(encoding="utf-8", errors="replace")
    cif_displacive = re.search(r"_iso_displacivemode_number\s+(\d+)", cif_text)
    cif_strain = re.search(r"_iso_strainmode_number\s+(\d+)", cif_text)
    result["export_consistency"] = {
        "cif_displacive_mode_count": int(cif_displacive.group(1)) if cif_displacive else None,
        "cif_strain_mode_count": int(cif_strain.group(1)) if cif_strain else None,
        "topas_mode_parameter_count": len(re.findall(
            r"(?m)^\s*prm\s+!a\d+\s+[-+0-9.Ee]+", required["topas"].read_text(
                encoding="utf-8", errors="replace",
            ),
        )),
        "isoviz_has_displacive_modes": "!displacivemodelist" in isoviz_text,
    }
    if result["export_consistency"]["cif_displacive_mode_count"] != len(displacement_rows):
        result["issues"].append("CIF and Complete modes details displacive-mode counts differ")
    if result["export_consistency"]["cif_strain_mode_count"] != len(strain_rows):
        result["issues"].append("CIF and Complete modes details strain-mode counts differ")
    if result["export_consistency"]["topas_mode_parameter_count"] != len(displacement_rows):
        result["issues"].append("TOPAS and Complete modes details displacive-mode counts differ")
    if not result["export_consistency"]["isoviz_has_displacive_modes"]:
        result["issues"].append("IsoVIZ export lacks a displacive-mode list")

    structure_comparison, structure_distance_tolerance = _compare_exported_cif(
        expected_daughter, required["cif"],
    )
    result["exported_cif_vs_frozen_daughter"] = structure_comparison
    _amplitude_tolerance, amplitude_tolerance_evidence = (
        _amplitude_comparison_tolerance(case, structure_comparison, local_case)
    )
    result["shared_mode_amplitude_tolerance"] = amplitude_tolerance_evidence
    amplitude_tolerance_ready = amplitude_tolerance_evidence["status"] == "available"
    if not amplitude_tolerance_ready:
        result["inconclusive_reasons"].append(
            "noisy As/Ap sensitivity is unproven; shared-mode comparisons were skipped: "
            + amplitude_tolerance_evidence["reason"]
        )
    if structure_comparison["expected_atom_count"] != structure_comparison["actual_atom_count"]:
        result["issues"].append("exported CIF atom count differs from frozen daughter")
    if structure_comparison["expected_species"] != structure_comparison["actual_species"]:
        result["issues"].append("exported CIF species differ from frozen daughter")
    lattice_metric_tolerance = max(
        1e-6,
        4.0 * float(np.linalg.norm(official_strain_rounding_bound))
        if official_strain_rounding_bound is not None else 1e-6,
    )
    structure_comparison["lattice_metric_relative_tolerance"] = (
        lattice_metric_tolerance
    )
    if (
        structure_comparison["lattice_metric_relative_residual"]
        > lattice_metric_tolerance
    ):
        result["issues"].append(
            "exported CIF lattice differs from frozen daughter beyond the "
            "archived strain-vector rounding bound"
        )
    distance = structure_comparison["max_species_matched_distance_angstrom"]
    if distance is None or distance > structure_distance_tolerance:
        result["issues"].append(
            "exported CIF does not reproduce the frozen daughter within its "
            f"coordinate-precision tolerance {structure_distance_tolerance:.8g} angstrom: "
            f"{distance}"
        )

    result["local_case"] = local_case
    inverse_consistency = (local_case or {}).get("inverse_consistency", {})
    identity_readiness = (local_case or {}).get("identity_readiness", {})
    inverse_passed = (
        isinstance(inverse_consistency, Mapping)
        and inverse_consistency.get("status") == "pass"
    )
    identity_status = (
        identity_readiness.get("status")
        if isinstance(identity_readiness, Mapping) else None
    )
    result["local_evidence_layers"] = {
        "inverse_consistency": inverse_consistency,
        "identity_readiness": identity_readiness,
    }
    if not inverse_passed:
        result["issues"].append("matching local inverse_consistency is absent or not passing")
    if identity_status in {"inconclusive", "unresolved"}:
        result["inconclusive_reasons"].append(
            "local scientific mode identities are unresolved; "
            "identity-based normfactor and As/Ap comparisons were skipped"
        )
    elif identity_status != "pass":
        result["issues"].append("matching local identity_readiness is absent or failed")
    if local_case is not None and local_case.get("status") == "fail":
        result["issues"].append("matching local validation case failed")
    result["local_identity_amplitude_comparison"] = {
        "status": (
            "ready"
            if inverse_passed and identity_status == "pass" and amplitude_tolerance_ready
            else "skipped"
        ),
        "requires": [
            "inverse_consistency=pass", "identity_readiness=pass",
            "amplitude_tolerance=available",
        ],
    }
    if inverse_passed:
        local_strain = local_case.get(
            "strain_applied_engineering_q_parent_basis",
            local_case.get("strain_voigt_engineering"),
        )
        if official_applied_strain is None:
            pass
        elif not isinstance(local_strain, dict):
            result["issues"].append(
                "local validation case lacks applied strain components"
            )
        else:
            labels = (
                ("11", "22", "33", "2*23", "2*13", "2*12")
                if "11" in local_strain
                else ("xx", "yy", "zz", "2yz", "2xz", "2xy")
            )
            local_vector = np.asarray(
                [float(local_strain[label]) for label in labels]
            )
            difference = np.abs(local_vector - official_applied_strain)
            strain_tolerance = official_strain_rounding_bound + 5e-10
            result["local_vs_official_applied_strain"] = {
                "local": local_vector.tolist(),
                "official_isoviz": official_applied_strain.tolist(),
                "absolute_difference": difference.tolist(),
                "max_absolute_difference": float(np.max(difference)),
                "componentwise_tolerance": strain_tolerance.tolist(),
            }
            if bool(np.any(difference > strain_tolerance)):
                result["issues"].append(
                    "local applied strain differs from the official IsoVIZ tensor: "
                    f"component differences {difference.tolist()} exceed archived "
                    f"rounding bounds {strain_tolerance.tolist()}"
                )
    if (
        local_case is not None
        and inverse_passed
        and identity_status == "pass"
        and amplitude_tolerance_ready
        and definition_normfactors
    ):
        official_normfactor_rows: list[dict[str, Any]] = []
        official_by_identity: dict[
            tuple[str, str, str], list[dict[str, Any]]
        ] = {}
        for label, value in definition_normfactors.items():
            identity = _official_mode_identity(label)
            if identity is not None:
                official_row = {
                    "identity": list(identity),
                    "label": label,
                    "value": float(value),
                }
                official_normfactor_rows.append(official_row)
                official_by_identity.setdefault(identity, []).append(official_row)
        normalization_matches = []
        required_normfactor_labels = (
            set(case.get("contributions", {}))
            if case.get("context") == "parameter"
            else None
        )
        for local_label, local_value in local_case.get(
            "mode_normfactors_inverse_angstrom", {},
        ).items():
            identity = _local_mode_identity(local_label)
            candidates = official_by_identity.get(identity, []) if identity else []
            official = candidates[0] if len(candidates) == 1 else None
            row = {
                "local_label": local_label,
                "identity": list(identity) if identity is not None else None,
                "local_normfactor_inverse_angstrom": float(local_value),
                "official_label": official["label"] if official else None,
                "official_normfactor_inverse_angstrom": (
                    official["value"] if official else None
                ),
                "absolute_difference": (
                    abs(float(local_value) - official["value"])
                    if official else None
                ),
                "comparison_required": (
                    required_normfactor_labels is None
                    or local_label in required_normfactor_labels
                ),
                "comparison_method": "componentwise",
            }
            normalization_matches.append(row)

        gamma_normalization_comparisons = []
        gamma_normalization_labels: set[str] = set()
        for key, local_rows, official_rows in _complete_gamma_site_subspaces(
            normalization_matches, official_normfactor_rows,
        ):
            local_values = sorted(
                row["local_normfactor_inverse_angstrom"] for row in local_rows
            )
            official_values = sorted(row["value"] for row in official_rows)
            differences = np.abs(
                np.asarray(local_values, dtype=float)
                - np.asarray(official_values, dtype=float)
            )
            comparison = {
                "subspace_identity": list(key),
                "local_labels": [row["local_label"] for row in local_rows],
                "official_labels": [row["label"] for row in official_rows],
                "sorted_local_normfactors_inverse_angstrom": local_values,
                "sorted_official_normfactors_inverse_angstrom": official_values,
                "max_absolute_difference": float(np.max(differences)),
                "tolerance_inverse_angstrom": (
                    5.1e-6
                    + max(
                        float(
                            (local_case or {}).get(
                                "mode_source_quantization_bound", {}
                            ).get("mode_records", {})
                            .get(row["local_label"], {})
                            .get("normfactor_error_bound_inverse_angstrom", 0.0)
                        )
                        for row in local_rows
                    )
                ),
            }
            gamma_normalization_comparisons.append(comparison)
            gamma_normalization_labels.update(
                row["local_label"] for row in local_rows
            )
            for row in local_rows:
                row["comparison_method"] = "complete_gamma_site_subspace"
                row["official_subspace_labels"] = comparison["official_labels"]
            if (
                any(row["comparison_required"] for row in local_rows)
                and comparison["max_absolute_difference"]
                > comparison["tolerance_inverse_angstrom"]
            ):
                result["issues"].append(
                    "normfactor multiset mismatch for complete Gamma site subspace "
                    f"{comparison['subspace_identity']!r}"
                )

        for row in normalization_matches:
            source_normfactor_bound = float(
                (local_case or {}).get("mode_source_quantization_bound", {})
                .get("mode_records", {})
                .get(row["local_label"], {})
                .get("normfactor_error_bound_inverse_angstrom", 0.0)
            )
            row["source_quantization_bound_inverse_angstrom"] = (
                source_normfactor_bound
            )
            row["tolerance_inverse_angstrom"] = 5.1e-6 + source_normfactor_bound
            if (
                row["official_label"] is None
                and row["local_label"] not in gamma_normalization_labels
            ):
                result["issues"].append(
                    f"no unique official mode matches local mode {row['local_label']!r}"
                )
            elif (
                row["official_label"] is not None
                and row["comparison_required"]
                and row["absolute_difference"] > row["tolerance_inverse_angstrom"]
            ):
                result["issues"].append(
                    f"normfactor mismatch for {row['local_label']!r}: local "
                    f"{row['local_normfactor_inverse_angstrom']}, official "
                    f"{row['official_normfactor_inverse_angstrom']}"
                )
        result["overlapping_mode_normalization"] = normalization_matches
        result["complete_gamma_site_normalization_comparisons"] = (
            gamma_normalization_comparisons
        )
        official_amplitude_rows: list[dict[str, Any]] = []
        official_amplitudes_by_identity: dict[
            tuple[str, str, str], list[dict[str, Any]]
        ] = {}
        for official_row in displacement_rows:
            identity = _official_mode_identity(official_row["label"])
            if identity is None:
                continue
            indexed_row = {**official_row, "identity": list(identity)}
            official_amplitude_rows.append(indexed_row)
            official_amplitudes_by_identity.setdefault(identity, []).append(
                indexed_row
            )
        amplitude_matches = []
        required_amplitude_labels = (
            set(case.get("contributions", {}))
            if (
                case.get("context") == "parameter"
                or case.get("noise_sigma_fractional") is not None
            )
            else None
        )
        for local_label, local_as in local_case.get("actual_As_angstrom", {}).items():
            identity = _local_mode_identity(local_label)
            candidates = (
                official_amplitudes_by_identity.get(identity, [])
                if identity is not None else []
            )
            official_row = candidates[0] if len(candidates) == 1 else None
            local_ap = float(local_case.get("actual_Ap_angstrom", {})[local_label])
            row = {
                "local_label": local_label,
                "identity": list(identity) if identity is not None else None,
                "local_As_angstrom": float(local_as),
                "local_Ap_angstrom": local_ap,
                "official_label": official_row["label"] if official_row else None,
                "official_As_angstrom": official_row["As"] if official_row else None,
                "official_Ap_angstrom": official_row["Ap"] if official_row else None,
                "As_absolute_difference": (
                    abs(float(local_as) - official_row["As"])
                    if official_row else None
                ),
                "Ap_absolute_difference": (
                    abs(local_ap - official_row["Ap"])
                    if official_row else None
                ),
                "comparison_required": (
                    required_amplitude_labels is None
                    or local_label in required_amplitude_labels
                ),
                "comparison_method": "componentwise",
            }
            amplitude_matches.append(row)
        complete_subspaces = _complete_required_amplitude_subspaces(
            amplitude_matches,
            required_amplitude_labels,
            set(official_amplitudes_by_identity),
        )
        official_parameter_groups: dict[
            tuple[str, str, str], list[dict[str, Any]]
        ] = {}
        for official_row in official_amplitude_rows:
            identity = tuple(official_row["identity"])
            key = _parameter_mode_subspace_key(identity)
            if key is not None:
                official_parameter_groups.setdefault(key, []).append(official_row)

        comparison_groups: list[
            tuple[
                tuple[str, str, str],
                list[dict[str, Any]],
                list[dict[str, Any]],
                str,
                Fraction | None,
            ]
        ] = []
        for rows in complete_subspaces:
            key = _parameter_mode_subspace_key(tuple(rows[0]["identity"]))
            if key is not None:
                comparison_groups.append((
                    key,
                    rows,
                    official_parameter_groups[key],
                    "complete_subspace_norm",
                    None,
                ))

        expected_origin = np.asarray(
            expected_context.get("origin", [0.0, 0.0, 0.0]), dtype=float,
        )
        official_origin = np.asarray(result_subgroup["origin"], dtype=float)
        origin_shift = (official_origin - expected_origin).tolist()
        already_grouped_labels = {
            row["local_label"]
            for _, rows, _, _, _ in comparison_groups
            for row in rows
        }
        for key, rows, official_rows, phase_turns in (
            _complete_origin_phase_amplitude_subspaces(
                amplitude_matches,
                required_amplitude_labels,
                official_amplitude_rows,
                origin_shift,
            )
        ):
            if any(row["local_label"] in already_grouped_labels for row in rows):
                continue
            comparison_groups.append((
                key,
                rows,
                official_rows,
                "origin_phase_subspace_norm",
                phase_turns,
            ))
            already_grouped_labels.update(row["local_label"] for row in rows)

        subspace_labels = set(already_grouped_labels)
        subspace_comparisons = []
        for key, rows, official_rows, method, phase_turns in comparison_groups:
            for row in rows:
                row["comparison_method"] = method
                row["official_subspace_labels"] = [
                    official_row["label"] for official_row in official_rows
                ]
            local_as_norm = float(np.linalg.norm([row["local_As_angstrom"] for row in rows]))
            official_as_norm = float(
                np.linalg.norm([row["As"] for row in official_rows])
            )
            local_ap_norm = float(np.linalg.norm([row["local_Ap_angstrom"] for row in rows]))
            official_ap_norm = float(
                np.linalg.norm([row["Ap"] for row in official_rows])
            )
            as_tolerance, ap_tolerance = _subspace_amplitude_tolerances(
                amplitude_tolerance_evidence,
                len(official_rows),
            )
            comparison = {
                "subspace_identity": list(key),
                "local_labels": [row["local_label"] for row in rows],
                "official_labels": [row["label"] for row in official_rows],
                "local_dimension": len(rows),
                "official_dimension": len(official_rows),
                "local_As_norm_angstrom": local_as_norm,
                "official_As_norm_angstrom": official_as_norm,
                "As_norm_absolute_difference": abs(local_as_norm - official_as_norm),
                "local_Ap_norm_angstrom": local_ap_norm,
                "official_Ap_norm_angstrom": official_ap_norm,
                "Ap_norm_absolute_difference": abs(local_ap_norm - official_ap_norm),
                "As_tolerance_angstrom": as_tolerance,
                "Ap_tolerance_angstrom": ap_tolerance,
                "tolerance_angstrom": max(as_tolerance, ap_tolerance),
                "comparison_method": method,
                "origin_phase_turns": (
                    None if phase_turns is None else str(phase_turns)
                ),
            }
            subspace_comparisons.append(comparison)
            if (
                comparison["As_norm_absolute_difference"] > as_tolerance
                or comparison["Ap_norm_absolute_difference"] > ap_tolerance
            ):
                result["issues"].append(
                    "As/Ap complete-subspace norm mismatch for "
                    f"{comparison['subspace_identity']!r}: local "
                    f"({local_as_norm}, {local_ap_norm}), official "
                    f"({official_as_norm}, {official_ap_norm})"
                )

        for key, local_rows, official_rows in _complete_gamma_site_subspaces(
            amplitude_matches, official_amplitude_rows,
        ):
            for row in local_rows:
                row["comparison_method"] = "complete_gamma_site_subspace_norm"
                row["official_subspace_labels"] = [
                    official_row["label"] for official_row in official_rows
                ]
            local_as_norm = float(np.linalg.norm([
                row["local_As_angstrom"] for row in local_rows
            ]))
            official_as_norm = float(np.linalg.norm([
                row["As"] for row in official_rows
            ]))
            local_ap_norm = float(np.linalg.norm([
                row["local_Ap_angstrom"] for row in local_rows
            ]))
            official_ap_norm = float(np.linalg.norm([
                row["Ap"] for row in official_rows
            ]))
            as_tolerance, ap_tolerance = _subspace_amplitude_tolerances(
                amplitude_tolerance_evidence,
                len(official_rows),
            )
            comparison = {
                "subspace_identity": list(key),
                "local_labels": [row["local_label"] for row in local_rows],
                "official_labels": [row["label"] for row in official_rows],
                "local_As_norm_angstrom": local_as_norm,
                "official_As_norm_angstrom": official_as_norm,
                "As_norm_absolute_difference": abs(local_as_norm - official_as_norm),
                "local_Ap_norm_angstrom": local_ap_norm,
                "official_Ap_norm_angstrom": official_ap_norm,
                "Ap_norm_absolute_difference": abs(local_ap_norm - official_ap_norm),
                "As_tolerance_angstrom": as_tolerance,
                "Ap_tolerance_angstrom": ap_tolerance,
                "tolerance_angstrom": max(as_tolerance, ap_tolerance),
                "comparison_method": "complete_gamma_site_subspace_norm",
            }
            subspace_comparisons.append(comparison)
            subspace_labels.update(row["local_label"] for row in local_rows)
            if (
                required_amplitude_labels is None
                or any(
                    row["local_label"] in required_amplitude_labels
                    for row in local_rows
                )
            ) and (
                comparison["As_norm_absolute_difference"] > as_tolerance
                or comparison["Ap_norm_absolute_difference"] > ap_tolerance
            ):
                result["issues"].append(
                    "As/Ap complete Gamma-site subspace norm mismatch for "
                    f"{comparison['subspace_identity']!r}: local "
                    f"({local_as_norm}, {local_ap_norm}), official "
                    f"({official_as_norm}, {official_ap_norm})"
                )

        for row in amplitude_matches:
            component_as_tolerance, component_ap_tolerance = (
                _subspace_amplitude_tolerances(
                    amplitude_tolerance_evidence,
                    1,
                )
            )
            row["As_tolerance_angstrom"] = component_as_tolerance
            row["Ap_tolerance_angstrom"] = component_ap_tolerance
            if (
                row["official_label"] is None
                and row["local_label"] not in subspace_labels
            ):
                result["issues"].append(
                    "no unique official amplitude row matches local mode "
                    f"{row['local_label']!r}"
                )
            if (
                row["official_label"] is not None
                and row["comparison_required"]
                and row["local_label"] not in subspace_labels
                and (
                    row["As_absolute_difference"] > component_as_tolerance
                    or row["Ap_absolute_difference"] > component_ap_tolerance
                )
            ):
                result["issues"].append(
                    f"As/Ap mismatch for {row['local_label']!r}: local "
                    f"({row['local_As_angstrom']}, {row['local_Ap_angstrom']}), official "
                    f"({row['official_As_angstrom']}, {row['official_Ap_angstrom']}), "
                    f"tolerances As={component_as_tolerance}, Ap={component_ap_tolerance}"
                )
        result["overlapping_mode_amplitudes"] = amplitude_matches
        result["complete_subspace_amplitude_comparisons"] = subspace_comparisons
    if result["inconclusive_reasons"]:
        result["scientific_scope"] = (
            f"{case_id} audits archived official identities, amplitudes, normalization, "
            "and exported structure independently. Local inverse consistency and the "
            "applied strain tensor remain separately checked; unresolved scientific mode "
            "identities or amplitude-sensitivity evidence prevent asserting complete "
            "local/official normfactor and As/Ap equivalence."
        )
    elif is_zero_case:
        result["scientific_scope"] = (
            f"{case_id} proves the selected basis, atom mapping, zero-amplitude decomposition, "
            "normalization factors, and export consistency at the website's displayed precision. "
            "It does not test nonzero As/Ap values. The official P1 basis can span more modes "
            "than the local subgroup-restricted basis."
        )
    else:
        result["scientific_scope"] = (
            f"{case_id} compares signed nonzero As/Ap values when the official and local mode "
            "components share a fixed basis, and compares the invariant Euclidean norm when a "
            "complete multidimensional parameter-k or repeated Gamma site subspace is present. "
            "For a nontrivial parameter-k phase caused by the archived official origin, the "
            "complete targeted branch norm is compared because its real components may rotate "
            "or change sign. "
            "Normalization factors are checked for target modes; extra official P1 modes are "
            "retained but are not treated as local columns."
        )
    if case_id == "G05-origin-shift":
        auto_origin = _audit_g05_auto_origin(
            directory / "auto-origin",
            expected_daughter,
            result,
            html_inventory=html_inventory,
        )
        result["auto_origin_run"] = auto_origin
        result["issues"].extend(
            f"auto-origin: {issue}" for issue in auto_origin["issues"]
        )
        result["warnings"].extend(
            f"auto-origin: {warning}" for warning in auto_origin["warnings"]
        )
    result["status"] = (
        "fail" if result["issues"]
        else "inconclusive" if result["inconclusive_reasons"]
        else "pass_with_warnings" if result["warnings"]
        else "pass"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument("--local-report", type=Path, default=LOCAL_REPORT_PATH)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--parent")
    parser.add_argument("--case-id")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    local = json.loads(args.local_report.read_text(encoding="utf-8"))
    selection_gate = _selection_filter_gate(
        manifest,
        parent_filter=args.parent,
        case_filter=args.case_id,
    )
    local_report_sha256_at_start = _sha256(args.local_report)
    (
        html_evidence_at_start,
        html_inventories,
        html_inventory_errors,
    ) = _capture_selected_html_evidence(
        manifest,
        parent_filter=args.parent,
        case_filter=args.case_id,
    )
    provenance_at_start = _local_report_provenance_gate(
        local,
        manifest,
        args.manifest,
    )
    selected_case_count = selection_gate["selected_case_count"]
    cases = []
    if (
        selection_gate["status"] == "pass"
        and provenance_at_start["status"] == "pass"
    ):
        for parent_name, parent in manifest["parents"].items():
            if args.parent is not None and parent_name != args.parent:
                continue
            for case_id, case in parent["cases"].items():
                if args.case_id is not None and case_id != args.case_id:
                    continue
                local_case = (
                    local.get("parents", {})
                    .get(parent_name, {})
                    .get("cases", {})
                    .get(case_id)
                )
                cases.append(_audit_case(
                    parent_name,
                    case_id,
                    case,
                    parent,
                    manifest["readable_case_folders"][case_id],
                    local_case,
                    html_inventory=html_inventories.get((parent_name, case_id)),
                    html_inventory_error=html_inventory_errors.get(
                        (parent_name, case_id)
                    ),
                ))
    html_evidence_at_completion, _completion_inventories, _completion_errors = (
        _capture_selected_html_evidence(
            manifest,
            parent_filter=args.parent,
            case_filter=args.case_id,
        )
    )
    html_evidence_issues: list[str] = []
    if html_evidence_at_start["scan_error_count"]:
        html_evidence_issues.append(
            "one or more selected official HTML inventories could not be captured at start"
        )
    if html_evidence_at_completion["scan_error_count"]:
        html_evidence_issues.append(
            "one or more selected official HTML inventories could not be captured at completion"
        )
    if html_evidence_at_start["sha256"] != html_evidence_at_completion["sha256"]:
        html_evidence_issues.append(
            "selected official HTML evidence changed while the audit was running"
        )
    html_evidence_gate = {
        "status": "pass" if not html_evidence_issues else "fail",
        "issues": html_evidence_issues,
        "stable_during_run": (
            html_evidence_at_start["sha256"]
            == html_evidence_at_completion["sha256"]
        ),
        "at_start": html_evidence_at_start,
        "at_completion": html_evidence_at_completion,
    }
    provenance_at_completion = (
        _local_report_provenance_gate(local, manifest, args.manifest)
        if provenance_at_start["status"] == "pass"
        else provenance_at_start
    )
    local_report_sha256_at_completion = _sha256(args.local_report)
    provenance_gate_issues = list(provenance_at_start["issues"])
    if (
        provenance_at_start["status"] == "pass"
        and provenance_at_completion["status"] != "pass"
    ):
        provenance_gate_issues.extend(
            f"completion: {issue}"
            for issue in provenance_at_completion["issues"]
            if f"completion: {issue}" not in provenance_gate_issues
        )
    if (
        provenance_at_start["recomputed_signature"]
        != provenance_at_completion["recomputed_signature"]
    ):
        provenance_gate_issues.append(
            "signed Method 4 inputs changed while the official audit was running"
        )
    if local_report_sha256_at_start != local_report_sha256_at_completion:
        provenance_gate_issues.append(
            "local report changed while the official audit was running"
        )
    provenance_gate = {
        "status": "pass" if not provenance_gate_issues else "fail",
        "issues": provenance_gate_issues,
        "at_start": provenance_at_start,
        "at_completion": provenance_at_completion,
        "local_report_sha256_at_start": local_report_sha256_at_start,
        "local_report_sha256_at_completion": local_report_sha256_at_completion,
    }
    downloaded = [case for case in cases if case.get("status") != "not_downloaded"]
    failures = [case for case in downloaded if case.get("status") == "fail"]
    warnings = [
        case for case in downloaded if case.get("status") == "pass_with_warnings"
    ]
    inconclusive = [
        case for case in downloaded if case.get("status") == "inconclusive"
    ]
    not_downloaded = [
        case for case in cases if case.get("status") == "not_downloaded"
    ]
    unexpected_statuses = sorted({
        str(case.get("status"))
        for case in cases
        if case.get("status") not in {
            "pass", "pass_with_warnings", "inconclusive", "fail", "not_downloaded",
        }
    })
    summary_status = (
        "fail" if (
            failures
            or provenance_gate["status"] != "pass"
            or html_evidence_gate["status"] != "pass"
        )
        else "inconclusive" if (
            selection_gate["status"] != "pass"
            or selected_case_count == 0
            or len(cases) != selected_case_count
            or len(downloaded) != selected_case_count
            or not_downloaded
            or inconclusive
            or unexpected_statuses
        )
        else "pass_with_warnings" if warnings
        else "pass"
    )
    report = {
        "schema": 5,
        "kind": "method4_official_archive_audit",
        "generated_at": datetime.now(UTC).isoformat(),
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": _sha256(args.manifest),
        "local_report": str(args.local_report.resolve()),
        "local_report_sha256": local_report_sha256_at_completion,
        "source_signature": _source_signature([
            Path(__file__).resolve(),
            SOURCE_PATHS[1],
            SOURCE_PATHS[2],
            args.manifest,
            args.local_report,
        ]),
        "selection_gate": selection_gate,
        "local_provenance_gate": provenance_gate,
        "official_html_evidence_gate": html_evidence_gate,
        "scope_note": (
            "Downloaded official runs are audited read-only. Zero cases validate routing and "
            "zero preservation, not nonzero amplitude normalization."
        ),
        "summary": {
            "selected_case_count": selected_case_count,
            "audited_case_count": len(cases),
            "blocked_case_count": selected_case_count - len(cases),
            "downloaded_case_count": len(downloaded),
            "passed_case_count": sum(
                case.get("status") == "pass" for case in downloaded
            ),
            "warning_case_count": len(warnings),
            "inconclusive_case_count": len(inconclusive),
            "failed_case_count": len(failures),
            "provenance_gate_failure_count": int(provenance_gate["status"] != "pass"),
            "html_evidence_gate_failure_count": int(
                html_evidence_gate["status"] != "pass"
            ),
            "not_downloaded_case_count": sum(
                case.get("status") == "not_downloaded" for case in cases
            ),
            "unexpected_case_statuses": unexpected_statuses,
            "status": summary_status,
        },
        "cases": cases,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["summary"], ensure_ascii=False))
    print(args.json_output.resolve())
    return {"fail": 1, "inconclusive": 2}.get(report["summary"]["status"], 0)


if __name__ == "__main__":
    raise SystemExit(main())
