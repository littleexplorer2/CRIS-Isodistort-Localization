"""示例入口脚本 - 快速演示如何使用 Python API

说明:
- 本项目仅通过 Python API 使用（不再提供独立命令行脚本）。
- 运行示例: `python examples\01_basic_workflow.py`
"""

from pathlib import Path
import sys


def _print_quick_start():
    here = Path(__file__).parent
    example = here / "examples" / "01_basic_workflow.py"
    print("isodistort - Quick start")
    print()
    if example.exists():
        print(f"Run the example workflow: python {example}")
    else:
        print("Run your Python script and call the IsoDistort API from isodistort.api")


if __name__ == "__main__":
    _print_quick_start()
    sys.exit(0)