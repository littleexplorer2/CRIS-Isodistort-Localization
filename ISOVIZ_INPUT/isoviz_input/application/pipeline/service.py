"""Coordinate the complete amplitude-input to IsoVIZ-launch use case."""
from __future__ import annotations

import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ...backend.amplitudes import PatchReport, apply_amplitudes, parse_amplitude_csv
from ...backend.inputs import TextDocument
from ...backend.launcher import find_isoviz_launcher, open_isoviz
from ...backend.validation import validate_patched_text


@dataclass(frozen=True, slots=True)
class PipelineRequest:
    """The two documents required by the ISOVIZ_INPUT workflow."""

    amplitudes: TextDocument
    structure: TextDocument


@dataclass(frozen=True, slots=True)
class PreparedIsoVIZ:
    """A reproducible static result, before any GUI acceptance claim."""

    request: PipelineRequest
    patched_text: str
    patch_report: PatchReport
    validation: dict[str, Any]
    output_path: Path | None

    @property
    def launchable(self) -> bool:
        return self.output_path is not None and self.validation["status"] != "fail"


@dataclass(frozen=True, slots=True)
class LaunchReceipt:
    """Records process dispatch only; it is not an IsoVIZ load acknowledgement."""

    prepared_path: Path
    launcher: Path | None


def prepare_isoviz(
    request: PipelineRequest,
    *,
    output_dir: Path | None = None,
) -> PreparedIsoVIZ:
    """Parse, patch, validate, and write a temporary launch file.

    Invalid results remain inspectable but are not written as launchable files.
    """
    modes = parse_amplitude_csv(
        request.amplitudes.text,
        source=request.amplitudes.source,
    )
    patched, patch_report = apply_amplitudes(request.structure.text, modes)
    validation = validate_patched_text(
        request.structure.text,
        patched,
        patch_report,
    )
    output_path = None
    if validation["status"] != "fail":
        output_path = _write_launch_copy(
            patched,
            stem=request.structure.stem,
            output_dir=output_dir,
        )
    return PreparedIsoVIZ(
        request=request,
        patched_text=patched,
        patch_report=patch_report,
        validation=validation,
        output_path=output_path,
    )


def dispatch_to_isoviz(
    prepared: PreparedIsoVIZ,
    *,
    launcher: Path | None = None,
) -> LaunchReceipt:
    """Dispatch a statically valid prepared file to the external IsoVIZ app."""
    if not prepared.launchable or prepared.output_path is None:
        errors = "; ".join(prepared.validation.get("errors", []))
        raise ValueError(f"Prepared input is not launchable: {errors or 'static validation failed'}")
    selected = launcher if launcher is not None else find_isoviz_launcher()
    open_isoviz(prepared.output_path, launcher=selected)
    return LaunchReceipt(prepared_path=prepared.output_path, launcher=selected)


def _write_launch_copy(text: str, *, stem: str, output_dir: Path | None) -> Path:
    safe = re.sub(r"[^\w.-]+", "_", stem).strip("._") or "isoviz"
    if output_dir is None:
        fd, name = tempfile.mkstemp(prefix=f"{safe}_", suffix=".isoviz")
    else:
        output_dir.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(
            prefix=f"{safe}_",
            suffix=".isoviz",
            dir=output_dir,
        )
    path = Path(name)
    with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    return path
