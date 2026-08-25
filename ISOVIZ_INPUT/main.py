"""Write GD amplitudes into an IsoVIZ ``.isoviz`` file and open the visualizer.

Usage (from the CRIS root, using CRIS/.venv):

  python ISOVIZ_INPUT/main.py --data path/to.csv --structure path/to.isoviz

Omit the flags to be prompted for paths.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from isoviz_input.amplitudes import list_mode_headers, patch_isoviz_file
from isoviz_input.launcher import find_isoviz_launcher, open_isoviz


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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Write amplitude CSV values into a subgroup .isoviz file and open IsoVIZ.",
    )
    parser.add_argument("--data", help="Path to the amplitude CSV (GD best-model parameters).")
    parser.add_argument("--structure", help="Path to the subgroup .isoviz file.")
    parser.add_argument(
        "--output",
        help="Patched .isoviz output path (default: ISOVIZ_INPUT/output/<name>_patched.isoviz).",
    )
    parser.add_argument("--no-open", action="store_true", help="Write the file but do not launch IsoVIZ.")
    args = parser.parse_args(argv)

    data = args.data or _pick_from_dir(ROOT / "data.csv", "*.csv") or _prompt("Amplitude CSV path")
    structure = (
        args.structure
        or _pick_from_dir(ROOT / "subgroup.isoviz", "*.isoviz")
        or _prompt("Subgroup .isoviz path")
    )
    if not data or not structure:
        print("Both a data CSV and a subgroup .isoviz path are required.")
        return 2

    data_path = Path(data).expanduser()
    structure_path = Path(structure).expanduser()
    if args.output:
        output_path = Path(args.output).expanduser()
    else:
        output_dir = ROOT / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{structure_path.stem}_patched.isoviz"

    print(f"[data]      {data_path}")
    print(f"[structure] {structure_path}")
    print(f"[output]    {output_path}")
    report = patch_isoviz_file(structure_path, data_path, output_path)
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

    if args.no_open:
        print("[done] File written. IsoVIZ was not launched (--no-open).")
        return 0
    launcher = find_isoviz_launcher()
    if launcher is not None:
        print(f"[isoviz]    {launcher}")
    open_isoviz(output_path, launcher=launcher)
    print("[done] IsoVIZ launch requested.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
