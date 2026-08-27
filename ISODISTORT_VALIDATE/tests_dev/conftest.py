"""pytest 共享配置：把 ISODISTORT_VALIDATE 根目录加入 sys.path，以便导入 isodistort_validate。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
