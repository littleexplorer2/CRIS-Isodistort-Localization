"""
自定义异常类

异常层次：
    IsodistortError（基异常）
    ├── WrapperError        二进制封装层
    │   ├── WrapperRunError     运行失败
    │   ├── WrapperTimeoutError 超时
    │   └── OutputParseError    输出解析失败
    └── DistortionError     畸变业务层
        └── PhasePathError      相变路径参数错误
"""


class IsodistortError(Exception):
    """项目基异常"""


# ---- 底层封装层异常 ----

class WrapperError(IsodistortError):
    """二进制封装基异常"""


class WrapperRunError(WrapperError):
    """二进制程序运行失败"""

    def __init__(self, binary: str, returncode: int, stderr: str = ""):
        self.binary = binary
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(f"{binary} 执行失败 (exit code {returncode}): {stderr}")


class WrapperTimeoutError(WrapperError):
    """二进制程序运行超时"""

    def __init__(self, binary: str, timeout: float):
        super().__init__(f"{binary} 执行超时 ({timeout}s)")


class OutputParseError(WrapperError):
    """输出解析失败"""

    def __init__(self, binary: str, message: str):
        super().__init__(f"解析 {binary} 输出失败: {message}")


# ---- 畸变层异常 ----

class DistortionError(IsodistortError):
    """畸变计算基异常"""


class PhasePathError(DistortionError):
    """相变路径参数错误"""
