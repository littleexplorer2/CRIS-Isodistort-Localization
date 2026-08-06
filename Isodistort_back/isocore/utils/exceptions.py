"""
自定义异常类
"""


class IsodistortError(Exception):
    """项目基异常"""
    pass


LocalIsodistortError = IsodistortError


# ---- 底层封装层异常 ----

class WrapperError(IsodistortError):
    """二进制封装基异常"""
    pass


class WrapperRunError(WrapperError):
    """二进制程序运行失败"""

    def __init__(self, binary: str, returncode: int, stderr: str = ""):
        """Relative path: isocore/utils/exceptions.py"""

        self.binary = binary
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(f"{binary} 执行失败 (exit code {returncode}): {stderr}")


class WrapperTimeoutError(WrapperError):
    """二进制程序运行超时"""

    def __init__(self, binary: str, timeout: float):
        """Relative path: isocore/utils/exceptions.py"""

        super().__init__(f"{binary} 执行超时 ({timeout}s)")


class OutputParseError(WrapperError):
    """输出解析失败"""

    def __init__(self, binary: str, message: str):
        """Relative path: isocore/utils/exceptions.py"""
        
        super().__init__(f"解析 {binary} 输出失败: {message}")


# ---- 结构层异常 ----

class StructureError(IsodistortError):
    """结构处理基异常"""
    pass


class SiteMappingError(StructureError):
    """位点匹配失败"""
    pass


class CoordinateTransformError(StructureError):
    """坐标变换失败"""
    pass


# ---- 畸变层异常 ----

class DistortionError(IsodistortError):
    """畸变计算基异常"""
    pass


class DistortionModeError(DistortionError):
    """畸变模式错误"""
    pass


class PhasePathError(DistortionError):
    """相变路径参数错误"""
    pass
