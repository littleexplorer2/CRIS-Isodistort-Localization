"""Command-line adapter for scripted and compatibility workflows."""
from __future__ import annotations

import argparse
from pathlib import Path

from ...application.pipeline import PipelineRequest, dispatch_to_isoviz, prepare_isoviz
from ...application.input_selection import load_text_document
from ...application.environment import find_isoviz_launcher
from ...application.input_selection import (
    best_model_root,
    list_csv_files,
    list_named_subfolders,
    resolve_best_model_csv,
    resolve_best_model_folder,
    resolve_isoviz_path,
    resolve_data_path,
)


def _prompt(label: str) -> str:
    return input(f"{label}: ").strip()


def _prompt_best_model_csv() -> Path:
    root = best_model_root()
    if not root.is_dir():
        raise FileNotFoundError(
            f"Best_Model_Parameters folder not found: {root}. "
            "Run the GD export first or pass --data."
        )
    folders = list_named_subfolders(root)
    print(f"Folders in {root}:")
    for index, folder in enumerate(folders, start=1):
        print(f"  {index}. {folder.name}")
    if not folders:
        raise FileNotFoundError(f"No model folders found in: {root}")
    folder = resolve_best_model_folder(_prompt("Best_Model_Parameters folder name"), root=root)
    files = list_csv_files(folder)
    print(f"CSV files in {folder}:")
    for index, path in enumerate(files, start=1):
        print(f"  {index}. {path.name}")
    if not files:
        raise FileNotFoundError(f"No CSV files found in: {folder}")
    return resolve_best_model_csv(folder.name, _prompt("CSV file name"), root=root)


def _resolve_data_path(args: argparse.Namespace) -> Path:
    if args.data:
        return resolve_data_path(args.data)
    if args.folder or args.csv_name:
        if not (args.folder and args.csv_name):
            raise ValueError("--folder and --csv-name must be provided together.")
        return resolve_best_model_csv(args.folder, args.csv_name)
    return _prompt_best_model_csv()


def _resolve_structure_path(args: argparse.Namespace) -> Path:
    raw = args.structure or _prompt("Absolute path of the crystal .isoviz file")
    return resolve_isoviz_path(raw)


def _print_report(prepared) -> None:
    report = prepared.patch_report
    print(f"[static]    {prepared.validation['status']}")
    print(f"[matched]   {len(report.matched)} mode(s)")
    for label, amplitude in report.matched:
        print(f"  {amplitude:10.5f}  {label}")
    if report.unmatched_csv:
        print("[csv leftover]")
        for name in report.unmatched_csv:
            print(f"  {name}")
    if report.unmatched_isoviz:
        print("[isoviz modes without CSV values; left at original amplitude]")
        for label in report.unmatched_isoviz:
            print(f"  {label}")
    for error in prepared.validation["errors"]:
        print(f"[error]     {error}")
    if prepared.output_path is not None:
        print(f"[prepared]  {prepared.output_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Prepare amplitude input and optionally launch IsoVIZ (CLI mode).",
    )
    parser.add_argument("--data", help="Path to the amplitude CSV.")
    parser.add_argument("--structure", help="Path to the subgroup .isoviz file.")
    parser.add_argument("--folder", help="Subfolder under Best_Model_Parameters.")
    parser.add_argument("--csv-name", help="CSV filename inside --folder.")
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Create and validate the temporary file without launching IsoVIZ.",
    )
    args = parser.parse_args(argv)
    try:
        data_path = _resolve_data_path(args)
        structure_path = _resolve_structure_path(args)
        request = PipelineRequest(
            amplitudes=load_text_document(data_path, expected_suffix=".csv"),
            structure=load_text_document(
                structure_path,
                expected_suffix=".isoviz",
                errors="replace",
            ),
        )
        prepared = prepare_isoviz(request)
        _print_report(prepared)
        if not prepared.launchable:
            return 2
        if args.prepare_only:
            return 0
        launcher = find_isoviz_launcher()
        receipt = dispatch_to_isoviz(prepared, launcher=launcher)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"[error] {exc}")
        return 2
    print(f"[dispatch]  {receipt.launcher or 'Windows file association'}")
    print("[done] IsoVIZ launch requested; window opening and data recognition require observation.")
    return 0
