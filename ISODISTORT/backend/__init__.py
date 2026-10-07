"""
ISODISTORT backend package - 通用底层计算与共用设施。

本包只放"所有 Method 与网页都可能用到"的底层能力，任何功能改动若属于通用底层，
就改这里；某一 Method 专有的实现放在 :mod:`features` 对应包内。

子包：
    wrappers - iso / findsym / smodes 二进制封装、WSL 进程与暂存、ISO 缓存
    models   - 模式与微观来源领域模型
    tables   - CDML / Kovalev / 官网 k 点与空间群查表
    utils    - 精确有理数晶格与群论工具、文本解析、OPD 文本、配置、异常、自检
    api      - IsoDistort 会话 API 入口
"""

__version__ = "0.4.0"
