"""pytest 共享配置：确保测试可导入 VALIDATE 根目录模块（compare_cif/batch_compare）。

VALIDATE 为无包结构的脚本式工具，测试脚本直接 ``from batch_compare import ...``，
需将项目根目录加入 sys.path（由本 conftest 在收集阶段统一注入）。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
