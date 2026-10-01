r"""Prepare and optionally perform the manual IsoVIZ GUI acceptance check.

Run from the CRIS root, for example::

    .\.venv\Scripts\python.exe ISOVIZ_INPUT\tests_dev\manual\validate_isoviz_gui.py \
        --data "C:\\...\\best_model_parameters.csv" \
        --structure "C:\\...\\subgroup.isoviz"

Add ``--launch`` only when a desktop GUI may be opened. A successful launch
call proves dispatch, not that IsoVIZ rendered or understood the file, so the
script asks the operator for those two observations separately.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from isoviz_input.amplitudes import apply_amplitudes, read_amplitude_csv  # noqa: E402
from isoviz_input.launcher import find_isoviz_launcher, open_isoviz  # noqa: E402
from isoviz_input.paths import resolve_isoviz_path, strip_user_path  # noqa: E402
from isoviz_input.validation import validate_patched_text  # noqa: E402


def _answer(prompt: str) -> bool:
    return input(f"{prompt} [y/N]: ").strip().lower() in {"y", "yes"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Create static patch evidence and optionally launch the IsoVIZ GUI for manual acceptance."
    )
    parser.add_argument("--data", required=True, help="Amplitude CSV path.")
    parser.add_argument("--structure", required=True, help="Original subgroup .isoviz path.")
    parser.add_argument(
        "--launch",
        action="store_true",
        help="Dispatch the prepared temporary file to IsoVIZ and ask for GUI confirmation.",
    )
    args = parser.parse_args(argv)

    data_path = Path(strip_user_path(args.data)).expanduser().resolve()
    structure_path = resolve_isoviz_path(args.structure)
    if not data_path.is_file():
        parser.error(f"CSV not found: {data_path}")

    original = structure_path.read_text(encoding="utf-8", errors="replace")
    patched, patch_report = apply_amplitudes(original, read_amplitude_csv(data_path))
    static = validate_patched_text(original, patched, patch_report)

    evidence_dir = Path(tempfile.mkdtemp(prefix="isoviz_gui_validation_"))
    launch_path = evidence_dir / f"{structure_path.stem}_patched.isoviz"
    report_path = evidence_dir / "validation_report.json"
    launch_path.write_text(patched, encoding="utf-8", newline="\n")
    launcher = find_isoviz_launcher()
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_path": str(data_path),
        "source_structure_path": str(structure_path),
        "prepared_structure_path": str(launch_path),
        "configured_launcher": str(launcher) if launcher else None,
        "static_validation": static,
        "gui_validation": {
            "requested": bool(args.launch),
            "dispatch_succeeded": False,
            "window_opened": "not_run",
            "data_recognized": "not_run",
            "status": "not_run",
        },
    }

    print(f"[static] {static['status']}")
    print(f"[file]   {launch_path}")
    print(f"[modes]  {static['mode_count']} total; {static['matched_count']} amplitudes patched")
    print(f"[meta]   {json.dumps(static['metadata'], ensure_ascii=False)}")
    print("[expected amplitudes visible in IsoVIZ]")
    for item in static["expected_visible_amplitudes"]:
        print(f"  {item['amp']: .5f}  {item['label']}")
    if static["unmatched_csv"]:
        print(f"[warning] unused CSV rows: {static['unmatched_csv']}")
    if static["unmatched_isoviz"]:
        print(f"[warning] modes retaining original amp: {static['unmatched_isoviz']}")

    if args.launch and static["status"] != "fail":
        try:
            open_isoviz(launch_path, launcher=launcher)
            report["gui_validation"]["dispatch_succeeded"] = True
            opened = _answer("Did an IsoVIZ window open with the prepared file")
            report["gui_validation"]["window_opened"] = "user_confirmed" if opened else "user_rejected"
            recognized = False
            if opened:
                recognized = _answer(
                    "Do the parent/child structure information, mode labels, and amp values match the list above"
                )
                report["gui_validation"]["data_recognized"] = (
                    "user_confirmed" if recognized else "user_rejected"
                )
            else:
                report["gui_validation"]["data_recognized"] = "not_observable"
            report["gui_validation"]["status"] = "pass" if opened and recognized else "fail"
        except (OSError, RuntimeError) as exc:
            report["gui_validation"]["status"] = "fail"
            report["gui_validation"]["error"] = str(exc)
            print(f"[launch error] {exc}")

    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[report] {report_path}")
    if static["status"] == "fail":
        return 2
    if args.launch and report["gui_validation"]["status"] != "pass":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
