"""
自定义异常类

异常层次：
    IsodistortError（基异常）
    ├── InputError              输入参数错误（类型 / 取值不合法）
    ├── DimensionMismatchError  维度不匹配（矩阵 / 向量长度）
    ├── WrapperError        二进制封装层
    │   ├── WrapperRunError     运行失败
    │   ├── WrapperTimeoutError 超时
    │   └── OutputParseError    输出解析失败
    └── DistortionError     畸变业务层
        ├── PhasePathError              相变路径参数错误
        ├── SymmetryIncompatibleError   对称不兼容
        └── NumericalSingularError      数值奇异
"""


class IsodistortError(Exception):
    """项目基异常"""


class InputError(IsodistortError):
    """输入参数错误（非法 nmod / d、类型错误、缺字段）"""


class DimensionMismatchError(IsodistortError):
    """矩阵或向量维度不匹配"""


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


class SymmetryIncompatibleError(DistortionError):
    """对称不兼容：波矢不满足超空间群约束，或对称操作集合不闭合"""


class NumericalSingularError(DistortionError):
    """数值奇异（度量矩阵不可逆、旋转矩阵奇异等）"""
