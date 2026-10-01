"""Static evidence for a CSV -> IsoVIZ amplitude patch.

This module deliberately does not claim that the Java GUI accepted a file.
That final check needs an explicit GUI launch and a human observation because
IsoVIZ does not expose a machine-readable load acknowledgement.
"""
from __future__ import annotations

from collections import Counter
from typing import Any

from .amplitudes import PatchReport, _header_match, list_mode_headers


def _mask_mode_amplitudes(isoviz_text: str) -> str:
    """Return text with only recognized mode-amplitude fields masked."""
    section = ""
    out: list[str] = []
    for line in isoviz_text.splitlines(keepends=True):
        raw = line.rstrip("\r\n")
        ending = line[len(raw) :]
        if raw.startswith("!"):
            section = raw[1:].split(maxsplit=1)[0].lower()
        header = _header_match(section, raw)
        if header is None:
            out.append(line)
            continue
        start, end = header.span("amp")
        out.append(f"{raw[:start]}<amp>{raw[end:]}{ending}")
    return "".join(out)


def section_scalar(isoviz_text: str, section_name: str) -> str | None:
    """Return the first non-empty line following ``!section_name``."""
    marker = f"!{section_name}".lower()
    lines = isoviz_text.splitlines()
    for index, raw in enumerate(lines):
        stripped = raw.strip()
        lowered = stripped.lower()
        if lowered == marker:
            inline = ""
        elif (
            lowered.startswith(marker)
            and len(stripped) > len(marker)
            and stripped[len(marker)].isspace()
        ):
            inline = stripped[len(marker) :].strip()
        else:
            continue
        if inline:
            return inline
        for following_line in lines[index + 1 :]:
            value = following_line.strip()
            if value:
                return value
        return None
    return None


def validate_patched_text(
    original_text: str,
    patched_text: str,
    patch_report: PatchReport,
    *,
    amplitude_tolerance: float = 5.1e-6,
) -> dict[str, Any]:
    """Check that a patch preserved mode identity and wrote expected amplitudes."""
    original = list_mode_headers(original_text)
    patched = list_mode_headers(patched_text)
    identity_preserved = [item[:2] for item in original] == [item[:2] for item in patched]
    non_amplitude_content_preserved = (
        _mask_mode_amplitudes(original_text) == _mask_mode_amplitudes(patched_text)
    )

    expected = Counter(
        (label, round(float(amplitude), 5))
        for label, amplitude in patch_report.matched
    )
    observed = Counter(
        (label, round(float(amplitude), 5))
        for _section, label, amplitude in patched
    )
    missing_expected: list[dict[str, Any]] = []
    for (label, amplitude), count in expected.items():
        available = sum(
            found_count
            for (found_label, found_amplitude), found_count in observed.items()
            if found_label == label
            and abs(found_amplitude - amplitude) <= amplitude_tolerance
        )
        if available < count:
            missing_expected.append(
                {
                    "label": label,
                    "amplitude": amplitude,
                    "expected_occurrences": count,
                    "observed_occurrences": available,
                }
            )

    errors: list[str] = []
    warnings: list[str] = []
    if not original:
        errors.append("The source .isoviz contains no recognized strain/displacive modes.")
    if not identity_preserved:
        errors.append("Mode labels, sections, count, or order changed during patching.")
    if not non_amplitude_content_preserved:
        errors.append("Content outside recognized mode-amplitude fields changed during patching.")
    if not patch_report.matched:
        errors.append("No CSV amplitude matched an IsoVIZ mode.")
    if missing_expected:
        errors.append("At least one matched CSV amplitude is absent from the patched mode table.")
    if patch_report.unmatched_csv:
        warnings.append("Some CSV rows were not used.")
    if patch_report.unmatched_isoviz:
        warnings.append("Some IsoVIZ modes retained their original amplitude.")

    metadata_keys = (
        "isoversion",
        "numberOfModulations",
        "parent_space_group_number",
        "parent_space_group_name_short",
        "child_space_group_number",
        "child_space_group_name_short",
    )
    metadata = {key: section_scalar(patched_text, key) for key in metadata_keys}
    status = "fail" if errors else ("pass_with_warnings" if warnings else "pass")
    return {
        "status": status,
        "mode_count": len(patched),
        "strain_mode_count": sum(1 for section, _label, _amp in patched if section == "strainmodelist"),
        "displacive_mode_count": sum(
            1 for section, _label, _amp in patched if section == "displacivemodelist"
        ),
        "matched_count": len(patch_report.matched),
        "identity_preserved": identity_preserved,
        "non_amplitude_content_preserved": non_amplitude_content_preserved,
        "missing_expected": missing_expected,
        "unmatched_csv": list(patch_report.unmatched_csv),
        "unmatched_isoviz": list(patch_report.unmatched_isoviz),
        "metadata": metadata,
        "expected_visible_amplitudes": [
            {"label": label, "amp": float(amplitude)}
            for label, amplitude in patch_report.matched
        ],
        "errors": errors,
        "warnings": warnings,
    }
