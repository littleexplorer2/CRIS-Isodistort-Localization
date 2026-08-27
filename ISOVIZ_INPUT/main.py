"""Read amplitude CSV + subgroup .isoviz and launch IsoVIZ.

Usage (from the CRIS root, using CRIS/.venv):

  python ISOVIZ_INPUT/main.py --data path/to.csv --structure path/to.isoviz

Omit the flags to pick files from input_content/ (or be prompted).
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from isoviz_input.amplitudes import apply_amplitudes, list_mode_headers, read_amplitude_csv
from isoviz_input.launcher import find_isoviz_launcher, open_isoviz
from isoviz_input.paths import DATA_DIR, STRUCTURE_DIR, ensure_input_content


def _prompt(label: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{label}{suffix}: ").strip()
    return value or default


def _pick_from_dir(folder: Path, pattern: str) -> str:
    if not folder.is_dir():
        return ""
    files = sorted(folder.glob(pattern))
    if not files:
        return ""
    print(f"Files in {folder}:")
    for i, path in enumerate(files, start=1):
        print(f"  {i}. {path.name}")
    raw = input("Choose a number or paste a path: ").strip()
    if raw.isdigit() and 1 <= int(raw) <= len(files):
        return str(files[int(raw) - 1])
    return raw


def _write_launch_copy(text: str, stem: str) -> Path:
    """Write a temporary .isoviz that IsoVIZ can open (not a project output folder)."""
    safe = re.sub(r"[^\w.-]+", "_", stem).strip("._") or "isoviz"
    fd, name = tempfile.mkstemp(prefix=f"{safe}_", suffix=".isoviz")
    path = Path(name)
    with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    return path


def _print_report(report, structure_path: Path) -> None:
    print(f"[matched]   {len(report.matched)} mode(s)")
    for label, amp in report.matched:
        print(f"  {amp:10.5f}  {label}")
    if report.unmatched_csv:
        print("[csv leftover]")
        for name in report.unmatched_csv:
            print(f"  {name}")
    if report.unmatched_isoviz:
        print("[isoviz modes without CSV values; left at original amplitude]")
        for label in report.unmatched_isoviz:
            print(f"  {label}")
    if not report.matched:
        headers = list_mode_headers(structure_path.read_text(encoding="utf-8", errors="replace"))
        print("[hint] IsoVIZ mode labels in this file:")
        for section, label, amp in headers[:20]:
            print(f"  {section}: {label} (amp={amp:.5f})")
        print("Match CSV 'Mode Name' to those labels, or use Mode=a1,a2,... in file order.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Read an amplitude CSV and a subgroup .isoviz, apply amplitudes, "
            "and launch IsoVIZ. This tool does not write a project output folder."
        ),
    )
    parser.add_argument("--data", help="Path to the amplitude CSV (GD best-model parameters).")
    parser.add_argument("--structure", help="Path to the subgroup .isoviz file.")
    args = parser.parse_args(argv)

    ensure_input_content()
    data = args.data or _pick_from_dir(DATA_DIR, "*.csv") or _prompt("Amplitude CSV path")
    structure = (
        args.structure
        or _pick_from_dir(STRUCTURE_DIR, "*.isoviz")
        or _prompt("Subgroup .isoviz path")
    )
    if not data or not structure:
        print("Both a data CSV and a subgroup .isoviz path are required.")
        return 2

    data_path = Path(data).expanduser()
    structure_path = Path(structure).expanduser()
    if not data_path.is_file():
        print(f"Data file not found: {data_path}")
        return 2
    if not structure_path.is_file():
        print(f"IsoVIZ structure file not found: {structure_path}")
        return 2

    print(f"[data]      {data_path}")
    print(f"[structure] {structure_path}")
    text = structure_path.read_text(encoding="utf-8", errors="replace")
    modes = read_amplitude_csv(data_path)
    patched, report = apply_amplitudes(text, modes)
    _print_report(report, structure_path)

    launch_path = _write_launch_copy(patched, structure_path.stem)
    launcher = find_isoviz_launcher()
    if launcher is not None:
        print(f"[isoviz]    {launcher}")
    try:
        open_isoviz(launch_path, launcher=launcher)
    except (OSError, RuntimeError) as exc:
        print(f"[error] {exc}")
        return 2
    print("[done] IsoVIZ launched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
