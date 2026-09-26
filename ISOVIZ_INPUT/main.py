"""Read amplitude CSV + subgroup .isoviz and launch IsoVIZ.

Usage (from the CRIS root, using CRIS/.venv):

  python ISOVIZ_INPUT/main.py

You will be asked for:
  1) a folder name under Desktop/Best_Model_Parameters
  2) a CSV file name inside that folder
  3) the absolute path of the crystal .isoviz file (quotes are stripped)

  python ISOVIZ_INPUT/main.py --data path/to.csv --structure path/to.isoviz
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
from isoviz_input.paths import (
    ensure_best_model_root,
    ensure_input_content,
    list_csv_files,
    list_named_subfolders,
    resolve_best_model_csv,
    resolve_best_model_folder,
    resolve_isoviz_path,
    strip_user_path,
)


def _prompt(label: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{label}{suffix}: ").strip()
    return value or default


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


def _prompt_best_model_csv() -> Path:
    root = ensure_best_model_root()
    folders = list_named_subfolders(root)
    print(f"Folders in {root}:")
    if folders:
        for i, folder in enumerate(folders, start=1):
            print(f"  {i}. {folder.name}")
    else:
        print("  (empty — run the GD notebook save cell first)")
    folder_name = strip_user_path(_prompt("Best_Model_Parameters folder name"))
    folder = resolve_best_model_folder(folder_name, root=root)

    files = list_csv_files(folder)
    if files:
        print(f"CSV files in {folder}:")
        for i, path in enumerate(files, start=1):
            print(f"  {i}. {path.name}")
    else:
        print(f"No CSV files in {folder}")
    filename = strip_user_path(_prompt("CSV file name"))
    return resolve_best_model_csv(folder.name, filename, root=root)


def _prompt_isoviz_path() -> Path:
    raw = _prompt("Absolute path of the crystal .isoviz file")
    return resolve_isoviz_path(raw)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Read a GD amplitude CSV from Desktop/Best_Model_Parameters and a "
            "subgroup .isoviz, apply amplitudes, and launch IsoVIZ."
        ),
    )
    parser.add_argument("--data", help="Path to the amplitude CSV (skips folder/file prompts).")
    parser.add_argument("--structure", help="Path to the subgroup .isoviz file.")
    parser.add_argument("--folder", help="Subfolder name under Desktop/Best_Model_Parameters.")
    parser.add_argument("--csv-name", help="CSV file name inside that subfolder.")
    args = parser.parse_args(argv)

    ensure_input_content()
    ensure_best_model_root()

    try:
        if args.data:
            data_path = Path(strip_user_path(args.data)).expanduser()
            if not data_path.is_file():
                print(f"Data file not found: {data_path}")
                return 2
        elif args.folder and args.csv_name:
            data_path = resolve_best_model_csv(args.folder, args.csv_name)
        else:
            data_path = _prompt_best_model_csv()

        if args.structure:
            structure_path = resolve_isoviz_path(args.structure)
        else:
            structure_path = _prompt_isoviz_path()
    except (OSError, ValueError) as exc:
        print(f"[error] {exc}")
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
