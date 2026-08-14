"""
二进制封装基类 - 统一管理 WSL 调用、短路径暂存与错误处理

背景
----
isobyu 中的 `iso` / `findsym` / `comsubs` / `smodes` 均为 Linux ELF 二进制，
Windows 下必须通过 WSL 运行。实测发现两个关键约束：

1. `findsym`（v7.1.6）内部用定长缓冲区存放“输入文件路径”，
   深路径（如 ``/mnt/c/Users/.../OneDrive/...``）会被截断导致崩溃，
   因此输入文件必须放在 WSL 侧的短路径下；
2. 二进制通过 ``ISODATA`` 环境变量定位数据库目录，且拼接文件名时
   “不带分隔符”，因此 ``ISODATA`` 必须以 ``/`` 结尾。

本基类统一解决：WSL 命令构造、短路径暂存、ISODATA 符号链接与
两种二进制调用方式（stdin 模式 / 输入文件参数模式）。
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from ..utils import WrapperRunError, WrapperTimeoutError, get_config


class BaseWrapper:
    """
    ISOTROPY 套件二进制程序封装基类。

    职责：
    1. Windows 下自动通过 WSL 调用 Linux 二进制（Linux 下原生调用）
    2. 在 WSL 侧建立短路径暂存目录与 ISODATA 符号链接
    3. 提供 stdin 模式与输入文件参数模式两种调用方式
    4. 统一的超时与错误处理
    """

    #: 在 WSL 用户主目录下创建的私有工作目录名（必须短：findsym 输入文件路径
    #: 缓冲区仅约 48 字符，超长会被截断导致崩溃）
    _WSL_DIR_NAME = ".id"
    #: 暂存目录名（短）
    _WSL_TMP_NAME = "tmp"
    #: ISODATA 符号链接名（短）
    _WSL_DATA_LINK = "data"
    #: 暂存文件保留时间（秒），超过则清理
    _STALE_SECONDS = 7 * 24 * 3600

    def __init__(self) -> None:
        self.cfg = get_config()
        self._mode = self._detect_mode()
        self._wsl_home: str = ""
        self._stage_dir: str = ""
        self._isodata_link: str = ""
        self._binary_dir_wsl: str = ""
        if self._mode == "wsl":
            self._init_wsl_environment()

    # ================================================================
    # 平台与环境初始化
    # ================================================================

    @staticmethod
    def _detect_mode() -> str:
        """返回运行模式：'wsl'（Windows + WSL）或 'native'（Linux）。"""
        if sys.platform.startswith("win"):
            return "wsl"
        if sys.platform.startswith("linux"):
            return "native"
        raise RuntimeError(
            f"暂不支持平台 {sys.platform}：isobyu 二进制为 Linux ELF 格式，"
            "请使用 Windows + WSL 或 Linux 环境。"
        )

    def _init_wsl_environment(self) -> None:
        """探测 WSL 用户主目录，建立短路径暂存目录与 ISODATA 符号链接。"""
        home = self._wsl_get_home()
        if not home:
            raise WrapperRunError(
                "wsl", 1, "无法获取 WSL 用户主目录，请确认已安装 WSL 并可用。"
            )
        self._wsl_home = home
        base = f"{home}/{self._WSL_DIR_NAME}"
        self._stage_dir = f"{base}/{self._WSL_TMP_NAME}"
        self._isodata_link = f"{base}/{self._WSL_DATA_LINK}"

        self._wsl_mkdir(self._stage_dir)
        self._wsl_cleanup_stale(self._stage_dir)
        self._wsl_link_isodata()

    def _wsl(self, command: str, timeout: float | None = None) -> subprocess.CompletedProcess:
        """执行一条 WSL 命令（Linux 模式下直接执行）。

        Args:
            command: 要在 WSL 中执行的 shell 命令
            timeout: 超时秒数，None 使用配置默认值

        Returns:
            subprocess.CompletedProcess
        """
        timeout = timeout if timeout is not None else float(self.cfg.timeout)
        try:
            if self._mode == "wsl":
                # 不使用 check=True：调用方需读取 returncode 决定如何处理。
                # S603/S607：命令由内部构造，路径均经 shlex.quote 转义；
                # 依赖 PATH 定位 wsl.exe 是刻意行为。
                result = subprocess.run(  # noqa: PLW1510, S603
                    ["wsl.exe", "-e", "sh", "-c", command],  # noqa: S607
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    # 关键：子进程不得继承父进程 stdin，否则会吞掉交互输入
                    # （脚本化/管道输入时会导致父进程 input() 直接 EOF）
                    stdin=subprocess.DEVNULL,
                    timeout=timeout,
                )
            else:
                result = subprocess.run(  # noqa: PLW1510, S603
                    ["sh", "-c", command],  # noqa: S607
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    stdin=subprocess.DEVNULL,
                    timeout=timeout,
                )
        except subprocess.TimeoutExpired as exc:
            name = command.split(maxsplit=1)[0] if command else "wsl"
            raise WrapperTimeoutError(name, timeout) from exc
        except FileNotFoundError as exc:
            raise WrapperRunError(
                "wsl", 1, "未找到 wsl.exe，请先安装并初始化 WSL。"
            ) from exc
        return result

    def _wsl_get_home(self) -> str:
        result = self._wsl("echo $HOME")
        if result.returncode != 0:
            return ""
        return result.stdout.strip()

    def _wsl_mkdir(self, path: str) -> None:
        result = self._wsl(f"mkdir -p {shlex.quote(path)}")
        if result.returncode != 0:
            raise WrapperRunError("wsl", result.returncode, result.stderr)

    def _wsl_cleanup_stale(self, stage_dir: str) -> None:
        """清理超过保留期的暂存文件，避免长期堆积。"""
        days = self._STALE_SECONDS // 86400
        result = self._wsl(
            f"find {shlex.quote(stage_dir)} -maxdepth 1 -type f "
            f"-mtime +{days} -delete 2>/dev/null || true"
        )
        _ = result  # 清理失败不影响主流程

    def _wsl_link_isodata(self) -> None:
        """在 WSL 侧创建到 isobyu 数据目录的短符号链接（幂等）。

        isobyu 数据目录在 Windows 侧的绝对路径通常很深，而二进制要求
        ISODATA 为短路径，因此统一通过符号链接暴露。
        """
        data_dir_win = str(self.cfg.resolve_path(self.cfg._cfg["isobyu"]["data_dir"]))
        data_dir_wsl = self._win_to_wsl(data_dir_win)
        cmd = (
            f"ln -sfn {shlex.quote(data_dir_wsl)} {shlex.quote(self._isodata_link)}"
        )
        result = self._wsl(cmd)
        if result.returncode != 0:
            raise WrapperRunError("wsl", result.returncode, result.stderr)

    # ================================================================
    # 路径转换与暂存
    # ================================================================

    def _win_to_wsl(self, win_path: str) -> str:
        """Windows 路径转 WSL 路径（优先 wslpath，失败时手动拼接）。"""
        if not win_path:
            return win_path
        if win_path.startswith("/"):
            return win_path
        if self._mode != "wsl":
            return win_path
        abs_path = os.path.abspath(win_path)
        try:
            # 依赖 PATH 定位 wsl.exe；参数为绝对路径，无注入面
            result = subprocess.run(  # noqa: S603
                ["wsl.exe", "wslpath", "-a", abs_path],  # noqa: S607
                capture_output=True,
                text=True,
                check=True,
                stdin=subprocess.DEVNULL,
                timeout=30,
            )
            wsl_path = result.stdout.strip()
            if wsl_path:
                return wsl_path
        except (subprocess.SubprocessError, OSError):
            pass
        # 手动拼接：C:\xxx -> /mnt/c/xxx
        path = Path(abs_path)
        drive = path.drive.rstrip(":").lower()
        rest = str(path).replace(path.drive, "", 1).replace("\\", "/")
        return f"/mnt/{drive}{rest}"

    def _wsl_bin_path(self, binary_path: Path) -> str:
        """将二进制文件路径转换为 WSL 可执行路径。

        Windows 侧路径直接转 WSL；若给出的是 WSL 路径则原样返回。
        """
        text = str(binary_path)
        if text.startswith("/"):
            return text
        return self._win_to_wsl(text)

    def _stage_text(self, prefix: str, text: str) -> str:
        """将文本写入 WSL 短路径暂存文件，返回 WSL 路径。

        Args:
            prefix: 文件名前缀（如 ``fs`` / ``iso``，务必保持简短，
                使最终路径（WSL 主目录 + 暂存目录 + 文件名）不超过约 45 字符，
                否则 findsym 的定长缓冲区会截断路径）
            text: 文件内容

        Returns:
            str: 暂存文件的 WSL 绝对路径（短路径）
        """
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".in", prefix=f"{prefix}_", delete=False,
            dir=str(self.cfg.temp_dir),
        ) as f_in:
            f_in.write(text)
            win_path = f_in.name
        try:
            wsl_tmp = self._stage_dir
            # 只保留随机后缀（去掉 Windows 侧前缀与扩展名再重新拼接）
            stem = os.path.basename(win_path)
            stem = stem.removesuffix(".in")
            wsl_name = f"{prefix}_{stem.split('_', 1)[-1]}.in"
            wsl_dst = f"{wsl_tmp}/{wsl_name}"
            cmd = (
                f"cp {shlex.quote(self._win_to_wsl(win_path))} {shlex.quote(wsl_dst)} && "
                f"chmod 600 {shlex.quote(wsl_dst)}"
            )
            result = self._wsl(cmd)
            if result.returncode != 0:
                raise WrapperRunError("wsl", result.returncode, result.stderr)
            return wsl_dst
        finally:
            try:
                os.remove(win_path)
            except OSError:
                pass

    # ================================================================
    # 二进制调用
    # ================================================================

    def _run_program(self, binary_path: Path, input_text: str,
                     input_via_stdin: bool,
                     extra_args: list[str] | None = None,
                     timeout: float | None = None) -> str:
        """运行二进制程序并返回标准输出。

        Args:
            binary_path: 可执行文件路径（Windows 绝对路径或 WSL 路径）
            input_text: 输入内容（stdin 模式直接作为 stdin；文件模式写入暂存文件）
            input_via_stdin: True 时通过 stdin 传入（iso 的交互命令流），
                False 时作为输入文件参数传入（findsym 关键字输入必须用文件参数）
            extra_args: 额外命令行参数
            timeout: 子进程超时秒数；None 使用配置默认值

        Returns:
            str: 程序标准输出
        """
        extra_args = extra_args or []
        wsl_bin = self._wsl_bin_path(binary_path)

        if input_via_stdin:
            wsl_in = self._stage_text("iso", input_text)
            quoted_args = " ".join(shlex.quote(a) for a in extra_args)
            shell_cmd = (
                f"{shlex.quote(wsl_bin)} {quoted_args} < {shlex.quote(wsl_in)}"
            )
        else:
            wsl_in = self._stage_text("fs", input_text)
            shell_cmd = (
                f"{shlex.quote(wsl_bin)} {' '.join(shlex.quote(a) for a in extra_args)} "
                f"{shlex.quote(wsl_in)}"
            )

        full_cmd = (
            f"export ISODATA={shlex.quote(self._isodata_path())}; "
            f"cd {shlex.quote(self._stage_dir)}; "
            f"{shell_cmd}"
        )
        result = self._wsl(full_cmd, timeout=timeout)

        if result.returncode != 0:
            raise WrapperRunError(binary_path.name, result.returncode, result.stderr)
        return result.stdout

    def _isodata_path(self) -> str:
        """返回供二进制使用的 ISODATA 路径（必须以 / 结尾）。"""
        if self._mode == "wsl":
            return self._isodata_link + "/"
        data_dir_win = str(self.cfg.resolve_path(self.cfg._cfg["isobyu"]["data_dir"]))
        return str(Path(data_dir_win)) + os.sep

    # ---- 面向子类的便捷入口 ----

    def run_stdin(self, binary_path: Path, input_text: str,
                  extra_args: list[str] | None = None,
                  timeout: float | None = None) -> str:
        """stdin 模式调用（适用于 iso 等从标准输入读取命令的程序）。

        Args:
            timeout: 子进程超时秒数；None 使用配置默认值
        """
        return self._run_program(binary_path, input_text, input_via_stdin=True,
                                 extra_args=extra_args, timeout=timeout)

    def run_input_file(self, binary_path: Path, input_text: str,
                       extra_args: list[str] | None = None) -> str:
        """输入文件参数模式调用（适用于 findsym 关键字输入）。"""
        return self._run_program(binary_path, input_text, input_via_stdin=False,
                                 extra_args=extra_args)

    # ---- 兼容旧接口 ----

    def _run_binary(self, binary_path: Path, input_text: str,
                    extra_args: list | None = None) -> str:
        """旧接口：默认走 stdin 模式（iso 交互命令流）。"""
        return self.run_stdin(binary_path, input_text, extra_args)
