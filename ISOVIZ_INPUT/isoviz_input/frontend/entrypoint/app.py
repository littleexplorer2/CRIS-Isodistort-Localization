"""Choose the desktop or command-line frontend without processing inputs."""
from __future__ import annotations

import sys


def _print_help() -> None:
    print(
        "ISOVIZ_INPUT\n\n"
        "  main.py                 Start the desktop GUI (default).\n"
        "  main.py --gui           Start the desktop GUI explicitly.\n"
        "  main.py --cli [options] Use the command-line workflow.\n"
        "  main.py --data ...      Legacy shorthand for CLI mode.\n\n"
        "Run `main.py --cli --help` for command-line options."
    )


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args == ["--gui"]:
        from ..gui import run_gui

        try:
            return run_gui()
        except RuntimeError as exc:
            print(f"[error] {exc}")
            return 2
    if args in (["-h"], ["--help"]):
        _print_help()
        return 0
    if args[0] == "--cli":
        args = args[1:]
    from ..cli import main as cli_main

    return cli_main(args)
