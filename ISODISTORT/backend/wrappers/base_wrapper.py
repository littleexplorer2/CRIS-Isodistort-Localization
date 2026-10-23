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
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import weakref
from pathlib import Path

from backend.utils import WrapperRunError, WrapperTimeoutError, get_config


_WSL_WORKSPACE_RE = re.compile(r"/var/tmp/id-[A-Za-z0-9]{8}\Z")
_UTF16_ASCII_RECORDS = (
    (re.compile(rb"(?:[\t\r\x20-\x7e]\x00)*(?:\n\x00|\Z)"), "utf-16-le"),
    (re.compile(rb"(?:\x00[\t\r\x20-\x7e])*(?:\x00\n|\Z)"), "utf-16-be"),
)


def _run_captured_command(
    args: list[str], *, timeout: float, env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[bytes]:
    """Capture a process without waiting for descendants to close output pipes.

    Windows ``run(PIPE, timeout=...)`` calls an unbounded ``communicate`` after
    killing the parent. A descendant can retain those pipe handles forever.
    Regular files need no reader threads or EOF handshake. Reaping only the
    process created here gets a separate one-second bound; no process tree or
    unrelated WSL process is killed. Process creation itself is an OS call,
    outside the wait timeout. Inherited temporary files may be deleted later
    when their last handle closes; their disappearance is never awaited.
    """
    with tempfile.TemporaryFile(mode="w+b") as stdout, tempfile.TemporaryFile(mode="w+b") as stderr:
        process = subprocess.Popen(  # noqa: S603 - caller constructs the argv
            args, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
            env=env,
        )
        try:
            returncode = process.wait(timeout=timeout)
        except BaseException as exc:
            try:
                process.kill()
                process.wait(timeout=1.0)
            except (OSError, subprocess.TimeoutExpired) as cleanup_error:
                raise exc from cleanup_error
            raise
        # Snapshot a finite length: descendants may still append after the
        # primary process exits, so an unbounded read could chase their writes.
        captured = []
        for stream in (stdout, stderr):
            length = os.fstat(stream.fileno()).st_size
            stream.seek(0)
            captured.append(stream.read(length))
        return subprocess.CompletedProcess(args, returncode, *captured)


def _run_wsl_command(args: list[str], *, timeout: float) -> subprocess.CompletedProcess[bytes]:
    """Use WSL's declared encoding; leave the parent and Linux locale unchanged.

    WSL 0.64.0 introduced WSL_UTF8=1 for the Windows client's stdout/stderr:
    https://github.com/microsoft/WSL/discussions/8641
    This also prevents UTF-16 client warnings and UTF-8 Linux errors from
    sharing an unframed byte stream. Only wsl.exe launches use this adapter.
    """
    env = os.environ.copy()
    env["WSL_UTF8"] = "1"
    return _run_captured_command(args, timeout=timeout, env=env)


def _remove_wsl_workspace(path: str, timeout: float) -> None:
    """Remove only one mktemp-owned short directory, never its symlink targets."""
    if not _WSL_WORKSPACE_RE.fullmatch(path):
        raise ValueError(f"Refusing to remove an unowned WSL workspace: {path!r}")
    result = _run_wsl_command(
        ["wsl.exe", "-e", "sh", "-c", f"rm -rf -- {shlex.quote(path)}"],
        timeout=timeout,
    )
    if result.returncode != 0:
        detail = "\n".join(
            part for raw in (result.stdout, result.stderr)
            if (part := _decode_process_output(raw))
        )
        raise WrapperRunError("wsl", result.returncode, detail)


def _decode_process_output(raw: bytes, *, diagnostic: bool = True) -> str:
    """Decode declared UTF-8 output, with bounded legacy diagnostic support.

    Success stdout is scientific/program text: no encoding inference, NUL
    deletion or replacement characters are allowed. Its existing outer-space
    trimming contract is retained.

    Diagnostic fallback supports BOM-declared UTF-16 streams and independent
    ASCII UTF-16LE/BE records terminated by an encoded LF (or stream end).
    Such records occur in redirected legacy WSL service messages. Matching
    their exact byte grammar, rather than a whole-stream NUL ratio, keeps a
    following UTF-8 Linux error readable. Arbitrary mid-line/interleaved or
    localized BOM-less UTF-16 is ambiguous and is not guessed: remaining bytes
    use UTF-8 with visible escapes. This fallback never feeds success stdout.
    """
    if not raw:
        return ""
    if not diagnostic:
        decoded = raw.decode("utf-8")
        if "\x00" in decoded:
            raise UnicodeError("WSL stdout contains NUL bytes; UTF-8 text is required")
        return decoded.strip()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        decoded = raw.decode("utf-16", errors="backslashreplace")
    else:
        parts = []
        position = 0
        while position < len(raw):
            for pattern, encoding in _UTF16_ASCII_RECORDS:
                match = pattern.match(raw, position)
                if match is not None and match.end() > position:
                    parts.append(match.group().decode(encoding))
                    position = match.end()
                    break
            else:
                end = raw.find(b"\n", position)
                end = len(raw) if end < 0 else end + 1
                parts.append(raw[position:end].decode("utf-8", errors="backslashreplace"))
                position = end
        decoded = "".join(parts)
    return decoded.replace("\x00", r"\x00").replace("\r\r\n", "\n").replace("\r\n", "\n").strip()


def _wsl_result_detail(result: subprocess.CompletedProcess[str]) -> str:
    """Return useful text regardless of which stream wsl.exe used."""
    return "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part and part.strip()
    )


def _is_wsl_access_denied(text: str) -> bool:
    compact = text.upper().replace("_", "")
    return "WSL/" in compact and "EACCESSDENIED" in compact


def _wsl_access_denied_hint() -> str:
    return (
        "Windows 拒绝当前 Python 进程访问 WSL (Wsl/E_ACCESSDENIED)。"
        "本机从 OneDrive 内的 .venv\\Scripts\\python.exe 启动时会稳定触发该限制；"
        "这不是缺少管理员权限，也不应通过反复重启 WSL 处理。请在 CRIS 根目录使用 "
        ".\\run_cris.ps1 <原 Python 参数>，例如 "
        ".\\run_cris.ps1 ISODISTORT\\main_web.py。"
    )


class BaseWrapper:
    """
    ISOTROPY 套件二进制程序封装基类。

    职责：
    1. Windows 下自动通过 WSL 调用 Linux 二进制（Linux 下原生调用）
    2. 在 WSL 侧建立短路径暂存目录与 ISODATA 符号链接
    3. 提供 stdin 模式与输入文件参数模式两种调用方式
    4. 统一的超时与错误处理
    """

    def __init__(self) -> None:
        self._init_runtime_state()
        self.cfg = get_config()
        self._mode = self._detect_mode()
        try:
            if self._mode == "wsl":
                self._init_wsl_environment()
            else:
                self._init_native_environment()
        except Exception as exc:
            try:
                self.close()
            except Exception as cleanup_error:
                # Explicit chaining retains both errors on Python 3.10 too;
                # failed cleanup must not replace the initialization error.
                raise exc from cleanup_error
            raise

    def _init_runtime_state(self) -> None:
        self._runtime_lock = threading.RLock()
        self._closed = False
        self._runtime_used = False
        self._runtime_isodata = None
        self._native_stage = None
        self._wsl_cleanup = None
        self._stage_dir = ""
        self._isodata_link = ""

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("Binary wrapper is closed")

    def close(self) -> None:
        """Release this instance's generated data and staging, after active calls."""
        with self._runtime_lock:
            if self._closed:
                return
            self._closed = True
            try:
                if self._wsl_cleanup is not None:
                    self._wsl_cleanup()
                if self._native_stage is not None:
                    self._native_stage.cleanup()
            finally:
                if self._runtime_isodata is not None:
                    self._runtime_isodata.cleanup()

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    def clone_runtime_from(self, source: BaseWrapper) -> None:
        """Give a fresh wrapper an independent copy of generated runtime state.

        Transaction snapshots must retain locally generated ISO tables after
        the preceding session closes. Only a generated private ISODATA and
        numeric ``i*.iso`` working-directory caches are copied. Packaged data
        is still read through this wrapper's own source link. No directories
        or mutable links are shared with the source instance.
        """
        if source is self:
            raise ValueError("Runtime cloning requires two distinct wrappers")
        # A fixed lock ordering prevents A<-B and B<-A deadlock. Object IDs
        # order infrastructure locks only; they never identify a candidate.
        first, second = sorted((self, source), key=id)
        with first._runtime_lock, second._runtime_lock:
            try:
                self._require_open()
                source._require_open()
                if self._mode != source._mode:
                    raise ValueError("Runtime cloning requires matching execution modes")
                if self._runtime_used or self._runtime_isodata is not None:
                    raise ValueError("Runtime cloning requires a fresh destination wrapper")
                if source._runtime_isodata is not None:
                    self._copy_writable_isodata(Path(source._runtime_isodata.name))
                self._copy_runtime_cache_from(source)
                self._runtime_used = True
            except Exception as exc:
                try:
                    self.close()
                except Exception as cleanup_error:
                    raise exc from cleanup_error
                raise

    def _copy_runtime_cache_from(self, source: BaseWrapper) -> None:
        if self._mode == "native":
            for cache in Path(source._stage_dir).iterdir():
                if not re.fullmatch(r"i\d+\.iso", cache.name, re.IGNORECASE):
                    continue
                if cache.is_symlink() or not cache.is_file():
                    raise ValueError(f"Generated ISO cache is not a regular file: {cache}")
                destination = Path(self._stage_dir) / cache.name
                if destination.exists():
                    raise ValueError(f"Runtime cache destination already exists: {destination}")
                shutil.copy2(cache, destination)
            return
        command = (
            f"for f in {shlex.quote(source._stage_dir)}/i*.iso; do "
            '[ -e "$f" ] || continue; name=${f##*/}; digits=${name#i}; '
            'digits=${digits%.iso}; case "$digits" in ""|*[!0-9]*) continue;; esac; '
            '[ -f "$f" ] && [ ! -L "$f" ] || exit 2; '
            f'dst={shlex.quote(self._stage_dir)}/"$name"; '
            '[ ! -e "$dst" ] || exit 3; cp -- "$f" "$dst" || exit; done'
        )
        result = self._wsl(command)
        if result.returncode != 0:
            raise WrapperRunError("wsl", result.returncode, "Cannot clone generated ISO cache: " + _wsl_result_detail(result))

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
        """Allocate an owner-only short directory for this wrapper instance.

        The 20-byte directory plus a staged input name stays below FINDSYM's
        approximately 48-byte filename buffer, independently of the username.
        /tmp may be tmpfs: natural WSL distro shutdown can erase it while the
        Windows owner is still alive. /var/tmp retains this private workspace
        across that idle lifecycle; no instance can retarget another's link.
        """
        self._require_open()
        result = self._wsl("umask 077; mktemp -d /var/tmp/id-XXXXXXXX")
        if result.returncode != 0:
            detail = _wsl_result_detail(result)
            if _is_wsl_access_denied(detail):
                detail = f"{detail}\n{_wsl_access_denied_hint()}"
            raise WrapperRunError(
                "wsl", result.returncode, detail or "Cannot create private WSL workspace"
            )
        stage = result.stdout.strip()
        if not _WSL_WORKSPACE_RE.fullmatch(stage):
            raise WrapperRunError("wsl", 1, f"Invalid private WSL workspace: {stage!r}")
        self._stage_dir = stage
        self._isodata_link = f"{stage}/data"
        self._wsl_cleanup = weakref.finalize(
            self, _remove_wsl_workspace, stage, float(self.cfg.timeout)
        )
        self._wsl_link_isodata()

    def _init_native_environment(self) -> None:
        """Keep native inputs and executable copies in a private short directory.

        Linux wheels need not preserve executable mode bits.  Packaged files
        remain read-only; only private copies receive execute permission.  A
        random directory also avoids trusting a predictable shared /tmp path.
        """
        try:
            self._native_stage = tempfile.TemporaryDirectory(prefix="id-", dir="/tmp")
        except OSError as exc:
            raise WrapperRunError("native", 1, f"Cannot create private staging directory: {exc}") from exc
        self._stage_dir = self._native_stage.name
        self._native_binary_cache: dict[Path, tuple[tuple[int, ...], Path]] = {}

    def _native_binary_path(self, binary_path: Path) -> str:
        """Return an executable path without changing the configured binary.

        Direct execution is preferred when the source already has permission.
        Otherwise a copy is reused only while the source stat identity matches;
        replacing a user-configured binary cannot silently reuse an old copy.
        """
        try:
            source = binary_path.resolve(strict=True)
            if not source.is_file():
                raise OSError(f"Not a regular binary file: {source}")
            if os.access(source, os.X_OK):
                return str(source)
            metadata = source.stat()
            signature = (
                metadata.st_dev, metadata.st_ino, metadata.st_size,
                metadata.st_mtime_ns, metadata.st_ctime_ns,
            )
            cached = self._native_binary_cache.get(source)
            if cached is not None and cached[0] == signature and os.access(cached[1], os.X_OK):
                return str(cached[1])
            descriptor, filename = tempfile.mkstemp(prefix="b-", dir=self._stage_dir)
            os.close(descriptor)
            destination = Path(filename)
            try:
                shutil.copyfile(source, destination)
                destination.chmod(0o700)
                if not os.access(destination, os.X_OK):
                    raise OSError(f"Private binary copy is not executable: {destination}")
            except OSError:
                destination.unlink(missing_ok=True)
                raise
            self._native_binary_cache[source] = (signature, destination)
            return str(destination)
        except OSError as exc:
            raise WrapperRunError(binary_path.name, 1, f"Cannot prepare native executable: {exc}") from exc

    def _wsl(self, command: str, timeout: float | None = None) -> subprocess.CompletedProcess:
        """执行一条 WSL 命令（Linux 模式下直接执行）。

        Args:
            command: 要在 WSL 中执行的 shell 命令
            timeout: 超时秒数，None 使用配置默认值

        Returns:
            subprocess.CompletedProcess
        """
        self._require_open()
        timeout = timeout if timeout is not None else float(self.cfg.timeout)
        try:
            if self._mode == "wsl":
                # 不使用 check=True：调用方需读取 returncode 决定如何处理。
                # S603/S607：命令由内部构造，路径均经 shlex.quote 转义；
                # 依赖 PATH 定位 wsl.exe 是刻意行为。
                raw_result = _run_wsl_command(
                    ["wsl.exe", "-e", "sh", "-c", command],  # noqa: S607
                    timeout=timeout,
                )
                try:
                    stdout = _decode_process_output(
                        raw_result.stdout or b"", diagnostic=raw_result.returncode != 0,
                    )
                except UnicodeError as exc:
                    detail = "\n".join(
                        part for raw in (raw_result.stdout, raw_result.stderr)
                        if (part := _decode_process_output(raw or b""))
                    )
                    raise WrapperRunError(
                        "wsl", raw_result.returncode,
                        "Unsupported WSL stdout encoding: UTF-8 text is required "
                        "(WSL_UTF8=1 support is required).\n" + detail,
                    ) from exc
                result = subprocess.CompletedProcess(
                    args=raw_result.args,
                    returncode=raw_result.returncode,
                    stdout=stdout,
                    stderr=_decode_process_output(raw_result.stderr or b""),
                )
            else:
                raw_result = _run_captured_command(
                    ["sh", "-c", command],  # noqa: S607
                    timeout=timeout,
                )
                result = subprocess.CompletedProcess(
                    args=raw_result.args,
                    returncode=raw_result.returncode,
                    stdout=raw_result.stdout.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n"),
                    stderr=raw_result.stderr.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n"),
                )
        except subprocess.TimeoutExpired as exc:
            name = command.split(maxsplit=1)[0] if command else "wsl"
            timeout_error = WrapperTimeoutError(name, timeout)
            # Killing wsl.exe does not prove every Linux descendant stopped.
            # Retire its private runtime before any future query/cache access.
            try:
                self.close()
            except Exception as cleanup_error:
                raise timeout_error from cleanup_error
            raise timeout_error from exc
        except FileNotFoundError as exc:
            raise WrapperRunError(
                "wsl", 1, "未找到 wsl.exe，请先安装并初始化 WSL。"
            ) from exc
        return result

    def _wsl_link_isodata(self, data_dir_win: str | None = None) -> None:
        """在 WSL 侧创建到 isobyu 数据目录的短符号链接（幂等）。

        isobyu 数据目录在 Windows 侧的绝对路径通常很深，而二进制要求
        ISODATA 为短路径，因此统一通过符号链接暴露。
        """
        if data_dir_win is None:
            data_dir_win = str(
                self.cfg.resolve_path(self.cfg._cfg["isobyu"]["data_dir"])
            )
        data_dir_wsl = self._win_to_wsl(data_dir_win)
        cmd = (
            f"ln -sfn {shlex.quote(data_dir_wsl)} {shlex.quote(self._isodata_link)}"
        )
        result = self._wsl(cmd)
        if result.returncode != 0:
            raise WrapperRunError("wsl", result.returncode, _wsl_result_detail(result))

    def _ensure_writable_isodata(self) -> None:
        """Switch this wrapper to a private disposable ISODATA copy.

        ISO appends on-demand parametric-k subgroup tables to its data
        directory. Packaged/reference ``resources/isobyu`` is read-only, so
        generation must never point there. The private copy lives under the
        configured runtime temp root and is owned by ``TemporaryDirectory``;
        it is removed at normal interpreter shutdown and is never a committed
        scientific result or answer database.
        """

        with self._runtime_lock:
            self._require_open()
            self._copy_writable_isodata()

    def _copy_writable_isodata(self, source: Path | None = None) -> None:
        if self._runtime_isodata is not None:
            return
        if source is None:
            source = Path(
                self.cfg.resolve_path(self.cfg._cfg["isobyu"]["data_dir"])
            ).resolve()
        if not source.is_dir():
            raise WrapperRunError(
                "iso",
                1,
                f"Cannot prepare private ISODATA: source is not a directory: {source}",
            )
        runtime_root = Path(self.cfg.temp_dir)
        runtime_root.mkdir(parents=True, exist_ok=True)
        holder = None
        try:
            holder = tempfile.TemporaryDirectory(
                prefix="isodata-", dir=str(runtime_root)
            )
            destination = Path(holder.name)
            shutil.copytree(source, destination, dirs_exist_ok=True)
            if self._mode == "wsl":
                self._wsl_link_isodata(str(destination))
            self._runtime_isodata = holder
        except (OSError, WrapperRunError, WrapperTimeoutError) as exc:
            if holder is not None:
                holder.cleanup()
            raise WrapperRunError(
                "iso", 1, f"Cannot prepare private writable ISODATA: {exc}"
            ) from exc

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
            result = _run_wsl_command(
                ["wsl.exe", "wslpath", "-a", abs_path],  # noqa: S607
                timeout=30,
            )
            wsl_path = _decode_process_output(result.stdout, diagnostic=False)
            if result.returncode == 0 and wsl_path:
                return wsl_path
        except (subprocess.SubprocessError, OSError, UnicodeError):
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
        if self._mode == "native":
            return self._native_binary_path(binary_path)
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
            dir=str(self.cfg.temp_dir), encoding="utf-8", newline="\n",
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
                raise WrapperRunError("wsl", result.returncode, _wsl_result_detail(result))
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
        with self._runtime_lock:
            self._require_open()
            self._runtime_used = True
            return self._run_program_locked(
                binary_path, input_text, input_via_stdin, extra_args, timeout
            )

    def _run_program_locked(self, binary_path: Path, input_text: str,
                            input_via_stdin: bool,
                            extra_args: list[str] | None,
                            timeout: float | None) -> str:
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
            f"export ISODATA={shlex.quote(self._isodata_path())} && "
            f"cd {shlex.quote(self._stage_dir)} && "
            f"{shell_cmd}"
        )
        result = self._wsl(full_cmd, timeout=timeout)

        if result.returncode != 0:
            raise WrapperRunError(binary_path.name, result.returncode, _wsl_result_detail(result))
        return result.stdout

    def _isodata_path(self) -> str:
        """返回供二进制使用的 ISODATA 路径（必须以 / 结尾）。"""
        runtime_isodata = getattr(self, "_runtime_isodata", None)
        if runtime_isodata is not None:
            if self._mode == "wsl":
                return self._isodata_link + "/"
            return str(Path(runtime_isodata.name)) + os.sep
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
