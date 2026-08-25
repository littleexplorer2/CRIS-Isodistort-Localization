"""Local ISODISTORT web UI.

    python main_web.py

Same as ``python web\\server.py``. Port: config/settings.yaml ``runtime.web_port``
(the next free port is used if that one is taken).

The Distortion panel downloads filtered result tables (Methods 1–4) and
subgroup files for exactly one Method (1 / 2 / 3) as CIF / ISOVIZ /
Complete modes details / TOPAS.STR. It does not scan output/. See README.md.
"""

import sys
from pathlib import Path

# 允许直接以脚本方式运行（python main_web.py）
# 先注入项目根目录，再导入 web 服务（有意放在函数外但 sys.path 之后）
_PROJECT_ROOT = Path(__file__).resolve().parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from web.server import main as web_main  # noqa: E402, I001


if __name__ == "__main__":
    sys.exit(web_main())
