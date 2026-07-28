"""
结果序列化 - 将完整计算结果保存为 JSON

对应阶段六：结果持久化与复现
实现方式：❌ 自研
"""
import json
import numpy as np
from pathlib import Path
from dataclasses import asdict
from typing import Any

from ..utils import get_config


class ResultSerializer:
    """
    计算结果序列化器

    支持将完整计算流程的结果保存为 JSON，便于：
    1. 断点复算
    2. 结果复现
    3. 与在线版 ISODISTORT 对标调试
    """

    def __init__(self, output_dir: str | Path = None):
        cfg = get_config()
        self.output_dir = Path(output_dir) if output_dir else cfg.output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def save(self, data: dict, filename: str) -> Path:
        """保存结果为 JSON"""
        path = self.output_dir / f"{filename}.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2,
                      default=self._numpy_default)
        return path

    def load(self, filename: str) -> dict:
        """加载 JSON 结果"""
        path = self.output_dir / f"{filename}.json"
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    @staticmethod
    def _numpy_default(obj: Any):
        """numpy 类型序列化兼容"""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
