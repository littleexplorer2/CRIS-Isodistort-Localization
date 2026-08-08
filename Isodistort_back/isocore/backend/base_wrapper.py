"""
二进制封装基类 - 统一管理子进程调用、临时文件与错误处理
"""
import os
import subprocess
import tempfile
from pathlib import Path
from abc import ABC

from ..utils import get_config, WrapperRunError, WrapperTimeoutError


class BaseWrapper(ABC):
    """
    ISOTROPY 套件二进制程序封装基类

    职责：
    1. 统一管理临时输入/输出文件
    2. 子进程调用与超时控制
    3. ISODATA 环境变量透传
    4. 统一的错误处理

    注意：isobyu 中的二进制为 Linux ELF 格式，
    Windows 环境下需要通过 WSL 调用，本基类自动检测并适配。
    """

    def __init__(self):
        """Relative path: isocore/backend/base_wrapper.py"""
        
        self.cfg = get_config()
        self._use_wsl = self._detect_wsl_need()

    def _detect_wsl_need(self) -> bool:
        """检测：在windows平台运行则自动通过 WSL 调用 Linux 二进制

        Relative path: isocore/backend/base_wrapper.py"""

        import platform
        if platform.system() != "Windows":
            return False
        # Windows 下默认使用 WSL 调用 isobyu 的 Linux 二进制
        return True

    def _run_binary(self, binary_path: Path, input_text: str,
                    extra_args: list = None) -> str:
        """
        运行二进制程序，传入输入文本，返回标准输出

        Args:
            binary_path: 可执行文件路径
            input_text: 写入标准输入的内容（或输入文件内容）
            extra_args: 额外命令行参数

        Returns:
            str: 程序标准输出
        
        Relative path: isocore/backend/base_wrapper.py"""

        extra_args = extra_args or []

        # 创建临时输入文件并写入input_text，以及获取临时文件路径
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".in", delete=False,
            dir=str(self.cfg.temp_dir)
        ) as f_in:
            f_in.write(input_text)
            in_path = f_in.name  # 输入文件在windows下的路径
        out_path = in_path.replace(".in", ".log")

        try:
            if self._use_wsl:
                # WSL 模式：将 Windows 路径转为 WSL 路径
                wsl_in = self._win_to_wsl_path(in_path)
                wsl_bin = self._win_to_wsl_path(str(binary_path))
                # 输入参数运行wsl binary：wsl binary < input.in
                cmd = ["wsl", wsl_bin] + extra_args
                with open(in_path, "r") as f_stdin:
                    result = subprocess.run(
                        cmd,
                        stdin=f_stdin,
                        capture_output=True,
                        text=True,
                        timeout=self.cfg.timeout,
                        env={**os.environ, "ISODATA": self._win_to_wsl_path(
                            str(self.cfg.resolve_path(self.cfg._cfg["isobyu"]["data_dir"])
                        ))}
                    )
            else:
                raise NotImplementedError("暂不支持在非 Windows 平台运行本程序。")
            
            if result.returncode != 0:
                raise WrapperRunError(
                    binary_path.name, result.returncode, result.stderr
                )

            return result.stdout

        except subprocess.TimeoutExpired:
            raise WrapperTimeoutError(binary_path.name, self.cfg.timeout)
        finally:
            # 清理临时文件
            for p in [in_path, out_path]:
                if os.path.exists(p):
                    try:
                        os.remove(p)
                    except OSError:
                        pass

    @staticmethod
    def _win_to_wsl_path(win_path: str) -> str:
        """Windows 路径转 WSL 路径，如 C:\\xxx -> /mnt/c/xxx

        Relative path: isocore/backend/base_wrapper.py"""

        if not win_path:
            return win_path

        # 已经是 WSL/Linux 路径时直接返回
        if win_path.startswith("/"):
            return win_path

        abs_path = os.path.abspath(win_path)

        # 优先使用 wsl 自带的 wslpath 方法，兼容中文、空格、UNC 等复杂路径场景
        try:
            result = subprocess.run(  # 调用子进程 wsl
                ["wsl", "wslpath", "-a", abs_path],
                capture_output=True,
                text=True,
                check=True,
            )
            wsl_path = result.stdout.strip()
            if wsl_path:
                return wsl_path
        except (subprocess.SubprocessError, OSError):
            pass

        # 当 wslpath 调用失败，程序手动拼接 WSL 挂载路径 C:\\x -> /mnt/c/x
        path = Path(abs_path)
        drive = path.drive.rstrip(":").lower()
        rest = str(path).replace(path.drive, "", 1).replace("\\", "/")
        return f"/mnt/{drive}{rest}"
