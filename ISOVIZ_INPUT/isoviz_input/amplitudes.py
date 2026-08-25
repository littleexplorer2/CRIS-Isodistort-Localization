"""Apply gradient-descent CSV amplitudes onto an official IsoVIZ ``.isoviz`` file."""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

_SECTION = re.compile(r"^!(\S+)")
_DISP_HEADER = re.compile(
    r"^(?P<pre>\s*)(?P<atom>\d+)(?P<w1>\s+)(?P<n>\d+)(?P<w2>\s+)"
    r"(?P<amp>[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?P<rest>\s+[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\s+\d+\s+(?P<label>\S.*))$"
)
_STRAIN_HEADER = re.compile(
    r"^(?P<pre>\s*)(?P<n>\d+)(?P<w1>\s+)"
    r"(?P<amp>[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?P<rest>\s+[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\s+\d+\s+(?P<label>\S.*))$"
)


@dataclass
class ModeAmplitude:
    name: str
    amplitude: float
    max_amplitude: float | None = None
    alias: str = ""


@dataclass
class PatchReport:
    matched: list[tuple[str, float]]
    unmatched_csv: list[str]
    unmatched_isoviz: list[str]


def read_amplitude_csv(path: str | Path) -> list[ModeAmplitude]:
    """Read GD-style or generic amplitude tables.

    Recognized headers (case-insensitive):
    Mode Name / modelabel / name; Best Model Parameter / amplitude / amp / value;
    Maximum Mode Amplitude / maxamp; Mode (a1, a2, ... as order fallback).

    Mode labels such as ``[0,0,1/6]...`` may contain commas. Quoted CSV (pandas
    default) is preferred; unquoted rows are repaired by merging extra fields
    back into the name column.
    """
    src = Path(path)
    if not src.is_file():
        raise FileNotFoundError(f"Data file not found: {src}")
    with src.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        try:
            raw_header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"CSV has no header row: {src}") from exc
        header = [(name or "").strip() for name in raw_header]
        if not any(header):
            raise ValueError(f"CSV has no header row: {src}")
        name_idx = _header_index(header, ("mode name", "modelabel", "label", "name"))
        amp_idx = _header_index(
            header,
            ("best model parameter", "amplitude", "amp", "value"),
        )
        max_idx = _header_index(header, ("maximum mode amplitude", "maxamp", "max amplitude"))
        alias_idx = _header_index(header, ("mode",))
        if amp_idx is None:
            raise ValueError(
                "CSV must contain an amplitude column "
                "(Best Model Parameter / amplitude / amp)."
            )
        rows: list[ModeAmplitude] = []
        for rec in reader:
            rec = _repair_row(rec, len(header), name_idx)
            if amp_idx >= len(rec):
                continue
            amp_raw = rec[amp_idx].strip()
            if not amp_raw:
                continue
            name = rec[name_idx].strip() if name_idx is not None and name_idx < len(rec) else ""
            alias = rec[alias_idx].strip() if alias_idx is not None and alias_idx < len(rec) else ""
            max_amp = None
            if max_idx is not None and max_idx < len(rec) and rec[max_idx].strip():
                try:
                    max_amp = float(rec[max_idx])
                except ValueError:
                    max_amp = None
            rows.append(
                ModeAmplitude(
                    name=name or alias,
                    amplitude=float(amp_raw),
                    max_amplitude=max_amp,
                    alias=alias,
                )
            )
    if not rows:
        raise ValueError(f"No amplitude rows in {src}")
    return rows


def list_mode_headers(isoviz_text: str) -> list[tuple[str, str, float]]:
    """Return ``(section, label, amp)`` for each IsoVIZ mode header."""
    found: list[tuple[str, str, float]] = []
    section = ""
    for raw in isoviz_text.splitlines():
        sec = _SECTION.match(raw)
        if sec:
            section = sec.group(1).lower()
            continue
        header = _header_match(section, raw)
        if header is not None:
            found.append((section, header.group("label").strip(), float(header.group("amp"))))
    return found


def apply_amplitudes(isoviz_text: str, modes: list[ModeAmplitude]) -> tuple[str, PatchReport]:
    """Return patched ``.isoviz`` text and a match report."""
    remaining = list(modes)
    matched: list[tuple[str, float]] = []
    unmatched_isoviz: list[str] = []
    sequential = 0
    section = ""
    out: list[str] = []
    for line in isoviz_text.splitlines(keepends=True):
        raw = line.rstrip("\r\n")
        ending = line[len(raw) :]
        sec = _SECTION.match(raw)
        if sec:
            section = sec.group(1).lower()
            out.append(line)
            continue
        header = _header_match(section, raw)
        if header is None:
            out.append(line)
            continue
        sequential += 1
        label = header.group("label").strip()
        hit = _consume(remaining, label, sequential)
        if hit is None:
            unmatched_isoviz.append(label)
            out.append(line)
            continue
        new_amp = _format_amp(hit.amplitude, header.group("amp"))
        raw = _rebuild_header(header, new_amp)
        matched.append((label, hit.amplitude))
        out.append(raw + ending)
    report = PatchReport(
        matched=matched,
        unmatched_csv=[m.name or m.alias for m in remaining],
        unmatched_isoviz=unmatched_isoviz,
    )
    return "".join(out), report


def patch_isoviz_file(
    structure_path: str | Path,
    data_path: str | Path,
    output_path: str | Path,
) -> PatchReport:
    src = Path(structure_path)
    if not src.is_file():
        raise FileNotFoundError(f"IsoVIZ structure file not found: {src}")
    dest = Path(output_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    text = src.read_text(encoding="utf-8", errors="replace")
    modes = read_amplitude_csv(data_path)
    patched, report = apply_amplitudes(text, modes)
    dest.write_text(patched, encoding="utf-8", newline="\n")
    return report


def _header_match(section: str, raw: str):
    if section == "displacivemodelist":
        return _DISP_HEADER.match(raw)
    if section == "strainmodelist":
        return _STRAIN_HEADER.match(raw)
    return None


def _rebuild_header(match: re.Match[str], new_amp: str) -> str:
    if "atom" in match.groupdict() and match.group("atom") is not None:
        return (
            f"{match.group('pre')}{match.group('atom')}{match.group('w1')}"
            f"{match.group('n')}{match.group('w2')}{new_amp}{match.group('rest')}"
        )
    return (
        f"{match.group('pre')}{match.group('n')}{match.group('w1')}"
        f"{new_amp}{match.group('rest')}"
    )


def _format_amp(value: float, original: str) -> str:
    text = f"{value:.5f}"
    if len(text) < len(original):
        return text.rjust(len(original))
    return text


def _consume(remaining: list[ModeAmplitude], label: str, order: int) -> ModeAmplitude | None:
    key = _norm(label)
    for i, mode in enumerate(remaining):
        if mode.name and _norm(mode.name) == key:
            return remaining.pop(i)
    for i, mode in enumerate(remaining):
        name = _norm(mode.name)
        if name and name != _norm(mode.alias) and (name in key or key in name):
            return remaining.pop(i)
    alias = f"a{order}"
    for i, mode in enumerate(remaining):
        if _norm(mode.alias) != alias:
            continue
        if not mode.name or _norm(mode.name) == alias:
            return remaining.pop(i)
    return None


def _header_index(header: list[str], names: tuple[str, ...]) -> int | None:
    lowered = [name.strip().lower() for name in header]
    for name in names:
        if name in lowered:
            return lowered.index(name)
    return None


def _repair_row(row: list[str], n_header: int, name_idx: int | None) -> list[str]:
    if len(row) <= n_header or name_idx is None:
        return row
    extra = len(row) - n_header
    end = name_idx + 1 + extra
    merged = ",".join(row[name_idx:end])
    return row[:name_idx] + [merged] + row[end:]


def _norm(text: str) -> str:
    return re.sub(r"\s+", "", text.strip().lower())
