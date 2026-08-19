"""
ISODISTORT 网页交互入口（方式 A）

直接运行即可在本地端口启动网页界面并自动打开浏览器，
无需在终端执行任何命令：

    python main_web.py

等价于 ``python web\\server.py``（端口可在 config/settings.yaml 的
runtime.web_port 修改；被占用时自动顺延）。
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
