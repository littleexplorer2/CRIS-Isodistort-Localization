"""Run the shared comparison commands without an independently installed project."""

from .frontend.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
