"""Local ISODISTORT web UI.

    python main_web.py

Direct execution starts the server, opens the browser, and supports the same
WSL-backed calculations as the repository launcher. ``setup_cris.py`` repairs
the affected Windows venv launcher when needed. Port: resources/config/settings.yaml
``runtime.web_port`` (the next free port is used if that one is taken).

The Distortion panel downloads filtered result tables (Methods 1–4) and
subgroup files for exactly one Method (1 / 2 / 3) as CIF / ISOVIZ /
Complete modes details / TOPAS.STR. It does not scan output/. See README.md.
"""

import sys
from pathlib import Path

# 允许从任意工作目录直接运行 ISODISTORT 根目录入口。
_PROJECT_ROOT = Path(__file__).resolve().parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from frontend.web.server import main as web_main  # noqa: E402, I001


if __name__ == "__main__":
    sys.exit(web_main())
