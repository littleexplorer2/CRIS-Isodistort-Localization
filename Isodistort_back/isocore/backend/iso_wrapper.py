"""
iso（Isotropy 9.6.1）封装：k 点枚举、不可约表示、各向同性子群、畴、模式基矢

实测命令接口（见 isobyu/iso 与 ISOTROPY 用户手册 isotropy_doc）：
    PAGE NOBREAK            关闭分页（否则大表会停在 “Enter RETURN to continue”）
    SCREEN 200              加宽屏幕列宽（默认 80 列会导致子群表折行）
    VALUE PARENT <sg>       选择母相空间群（数字或符号）
    VALUE KPOINT <label>    选择 k 点（Miller-Love 记号）
    VALUE KVALUE <n>,<v1>,.. 设置 k 点参数（n 为参数个数，如 "1,1/4"）
    VALUE IRREP <label>     选择不可约表示（Miller-Love 记号）
    VALUE DIRECTION <sym>   选择序参量方向（如 P1）
    VALUE WYCKOFF <letter>  选择 Wyckoff 位置（可多次）
    SHOW <flag>             控制 DISPLAY 输出内容
    DISPLAY KPOINT/IRREP/ISOTROPY/BUSH/PARENT
    QUIT

对应 ISODISTORT 官网阶段：
- 阶段二步骤4：枚举各向同性子群（DISPLAY ISOTROPY）
- 阶段三步步骤6/7：计算畸变模式基矢（DISPLAY BUSH + SHOW MODES）
- 阶段五步骤10：畴变体（SHOW DOMAIN）
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field

from ..utils import (
    OutputParseError,
    WrapperRunError,
    detect_blocked_generation,
    detect_missing_subgroup_db,
    extract_section,
    get_config,
    parse_bush_table,
    parse_domain_table,
    parse_fraction,
    parse_irrep_table,
    parse_kpoint_table,
    parse_subgroup_table,
)
from .base_wrapper import BaseWrapper

# ================================================================
# 数据模型
# ================================================================

@dataclass
class KPointInfo:
    """k 点信息（DISPLAY KPOINT）"""

    label: str                       # Miller-Love 记号，如 GM / DT
    coordinates: list[str]           # 坐标分量字符串，如 ["0","2a","0"]
    parameters: list[str]            # 自由参数字母，如 ["a"]
    is_special: bool                 # 无自由参数的特殊 k 点

    def __post_init__(self) -> None:
        self.parameters = sorted(
            {c for c in self.coordinates if re.search(r"[a-zA-Z]", c)}
        )
        self.is_special = not self.parameters


@dataclass
class IrrepInfo:
    """不可约表示信息（DISPLAY IRREP）"""

    label: str                       # Miller-Love 记号，如 GM1+
    dimension: int                   # 维度
    active: bool = True              # Landau active 标记（仅用于参考）


@dataclass
class SubgroupInfo:
    """各向同性子群信息（DISPLAY ISOTROPY）"""

    index: int                       # 本地枚举序号（用户选择句柄，0 起）
    space_group_number: int          # 子群空间群号
    space_group_symbol: str = ""     # 子群空间群短符号
    subgroup_index: int = 0          # iso 输出 Index（相对母相的子群指数 = 畴数）
    size: int = 1                    # 子群原胞相对母相的大小 s
    is_maximal: bool = False         # 是否为 maximal 子群
    opd_symbol: str = ""             # 序参量方向符号（如 P1）
    opd_vector: list[float] = field(default_factory=list)   # 序参量方向向量
    basis_vectors: list[list[float]] = field(default_factory=list)  # 超胞基矢（母相格单位）
    origin: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])  # 超胞原点
    k_point_label: str = ""          # 产生该子群的 k 点
    irrep_label: str = ""            # 产生该子群的不可约表示
    k_parameters: list[str] = field(default_factory=list)  # k 点参数值（带参数 k 点）

    def describe(self) -> str:
        return (
            f"#{self.index}: SG {self.space_group_number} {self.space_group_symbol} "
            f"k={self.k_point_label} IR={self.irrep_label} OPD={self.opd_symbol}"
        )


@dataclass
class DomainInfo:
    """畴变体信息（SHOW DOMAIN）"""

    domain_number: int               # 畴编号（1 起）
    generator: str = ""              # 域生成元，如 (C2y|0,0,0)
    space_group_number: int = 0
    space_group_symbol: str = ""
    subgroup_index: int = 0          # 母相中的子群指数（= 畴总数）
    opd_symbol: str = ""
    opd_vector: list[float] = field(default_factory=list)
    basis_vectors: list[list[float]] = field(default_factory=list)
    origin: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])


@dataclass
class BushMode:
    """模式基矢行（DISPLAY BUSH + SHOW MODES）"""

    irrep_label: str
    opd_symbol: str
    wyckoff_letter: str
    point: list[float]                       # 代表位点坐标
    displacements: list[list[float]] = field(default_factory=list)  # 位移向量（可多个，对应模式分量）


@dataclass
class DistortionMode:
    """畸变模式（序参量 + 原子位移基矢）"""

    irrep_label: str                 # 不可约表示标号，如 GM4-
    dimension: int = 1               # 模式维度
    mode_type: str = "displacement"  # 类型: displacement/order/strain/magnetic
    basis_vectors: list[list[float]] = field(default_factory=list)  # 兼容旧接口
    wyckoff_site: str = ""           # 对应 Wyckoff 位置
    k_point_label: str = ""          # k 点
    opd_symbol: str = ""             # 序参量方向
    bush_modes: list[BushMode] = field(default_factory=list)  # 原子级位移模式


# ================================================================
# iso 封装
# ================================================================

class IsoWrapper(BaseWrapper):
    """
    ISOTROPY (iso) 命令行程序封装

    iso 是交互式程序，通过关键字命令序列驱动（见模块 docstring）。
    本封装将 ISODISTORT 所需能力映射为真实的 iso 命令序列：

    1. list_k_points    枚举母相的全部 k 点（Method 2 的 k 点下拉列表）
    2. list_irreps      枚举指定 k 点下的不可约表示（Method 2 的 IR 下拉列表）
    3. list_subgroups   枚举（k 点, IR）对应的各向同性子群（Method 1/2/3 的核心）
    4. calc_distortion_modes  计算指定路径的畸变模式基矢（DISPLAY BUSH）
    5. get_domains      获取畴变体列表（Domains 输出）
    6. get_wyckoff_letters  获取空间群的 Wyckoff 位置（含代表坐标）
    """

    def __init__(self) -> None:
        super().__init__()
        self.binary = self.cfg.iso_bin

    # ================================================================
    # 通用命令流构造
    # ================================================================

    @staticmethod
    def _session(text: str) -> str:
        """包装一个 iso 会话：关闭分页、加宽屏幕，最后退出。"""
        return f"PAGE NOBREAK\nSCREEN 200\n{text}QUIT\n"

    def _run_session(self, commands: str, timeout: float | None = None) -> str:
        """运行一次 iso 会话并返回标准输出。

        Args:
            commands: 会话命令流（不含 QUIT）
            timeout: 子进程超时秒数；None 使用配置默认值
        """
        stdout = self.run_stdin(self.binary, self._session(commands), timeout=timeout)
        # 非零返回码下 run_stdin 已抛 WrapperRunError
        return stdout

    # ================================================================
    # k 点 / IR 枚举
    # ================================================================

    def list_k_points(self, parent_sg: int) -> list[KPointInfo]:
        """
        枚举母相空间群的全部 k 点（Method 2 的下拉列表数据源）。

        Args:
            parent_sg: 母相空间群号 (1-230)

        Returns:
            List[KPointInfo]
        """
        stdout = self._run_session(
            f"VALUE PARENT {parent_sg}\nSHOW KPOINT\nDISPLAY KPOINT\n"
        )
        try:
            rows = parse_kpoint_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析 k 点列表失败: {exc}") from exc
        return [KPointInfo(**row) for row in rows]

    def list_irreps(self, parent_sg: int, k_point: str,
                    k_parameters: Sequence[str] | None = None) -> list[IrrepInfo]:
        """
        枚举指定 k 点下的不可约表示（Method 2 的 IR 下拉列表数据源）。

        Args:
            parent_sg: 母相空间群号
            k_point: k 点标签（Miller-Love 记号）
            k_parameters: k 点参数值序列（如 ["1/4"]）；带参数 k 点必须提供

        Returns:
            List[IrrepInfo]
        """
        commands = [f"VALUE PARENT {parent_sg}", f"VALUE KPOINT {k_point}"]
        if k_parameters:
            commands.append(self._kvalue_command(k_parameters))
        commands += ["SHOW IRREP", "SHOW DIMENSION", "SHOW ACTIVE", "DISPLAY IRREP"]
        stdout = self._run_session("\n".join(commands) + "\n")
        try:
            rows = parse_irrep_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析不可约表示列表失败: {exc}") from exc
        return [IrrepInfo(**row) for row in rows]

    @staticmethod
    def _kvalue_command(k_parameters: Sequence[str]) -> str:
        """
        构造 VALUE KVALUE 命令。

        iso 语法：``VALUE KVALUE <参数个数>,<v1>,<v2>...``
        例如单个参数 1/4 -> "VALUE KVALUE 1,1/4"。
        """
        values = ",".join(str(v).strip() for v in k_parameters)
        return f"VALUE KVALUE {len(list(k_parameters))},{values}"

    # ================================================================
    # 子群枚举（Method 1/2/3 的核心）
    # ================================================================

    def list_subgroups(self, parent_sg: int,
                       k_point: str,
                       irrep_label: str,
                       k_parameters: Sequence[str] | None = None,
                       opd_symbol: str | None = None,
                       start_index: int = 0,
                       generate_if_missing: bool = False) -> list[SubgroupInfo]:
        """
        枚举指定 (k 点, IR) 下的各向同性子群。

        对应官网 Method 1（遍历全部特殊 k 点与 IR）与 Method 2
        （指定 k 点/IR）中“子群+序参量方向”列表。

        参数 k 点（如 LD、DT 等带 a/b/g 的点）说明：
        - iso 要求 **先选择 IR、再设置 KVALUE**（顺序颠倒会报
          “parameters not selected for k vector”）；
        - 参数 k 点的子群数据库默认不存在，iso 会询问是否在线生成
          （对应官网 “Generate isotropy subgroups”，可能耗时数分钟到数小时）。
          当 generate_if_missing=True 时自动应答并等待生成完成（生成的数据库
          会保存到暂存目录，后续查询立即返回）；否则抛出 WrapperRunError。

        Args:
            parent_sg: 母相空间群号
            k_point: k 点标签
            irrep_label: 不可约表示标签
            k_parameters: k 点参数（带参数 k 点必须提供）
            opd_symbol: 若指定则只返回该序参量方向对应的子群
            start_index: 本地枚举序号的起始值（用于多 IR 合并时连续编号）
            generate_if_missing: 子群数据库缺失时是否自动在线生成
                （默认 False；生成可能耗时很长，请谨慎开启）

        Returns:
            List[SubgroupInfo]

        Raises:
            OutputParseError: 解析失败
            WrapperRunError: 需要在线生成子群数据库等无法自动完成的情况
        """
        commands = [
            f"VALUE PARENT {parent_sg}",
            f"VALUE KPOINT {k_point}",
            f"VALUE IRREP {irrep_label}",
        ]
        if k_parameters:
            # 注意顺序：iso 要求先选 IR 再设置 KVALUE
            commands.append(self._kvalue_command(k_parameters))
        if opd_symbol:
            commands.append(f"VALUE DIRECTION {opd_symbol}")
        commands += [
            "SHOW SUBGROUP", "SHOW INDEX", "SHOW SIZE", "SHOW DIRECTION VEC",
            "SHOW BASIS", "SHOW ORIGIN", "SHOW MAXIMAL",
            "DISPLAY ISOTROPY",
        ]

        extra_input = ""
        timeout = None
        if generate_if_missing:
            # 应答 “Enter RETURN to continue” 提示，进入在线生成流程
            extra_input = "\n"
            timeout = float(get_config().generation_timeout)
        else:
            # 应答 “Enter any character to stop”：停止生成、回到命令提示符，
            # 避免程序等待更多输入导致 EOF 崩溃
            extra_input = "q\n"

        stdout = self._run_session("\n".join(commands) + "\n" + extra_input,
                                   timeout=timeout)

        # 先解析子群表：若在线生成成功（或数据库已存在），直接返回结果
        try:
            rows = parse_subgroup_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析子群表失败: {exc}") from exc

        if not rows:
            # 无子群行：区分“需要生成”与“参数错误”
            if detect_missing_subgroup_db(stdout):
                raise WrapperRunError(
                    "iso", 1,
                    f"k 点 {k_point} 的子群数据库在本地不存在，需要在线生成"
                    f"（官网对应 “Generate isotropy subgroups”，可能耗时数分钟到数小时）。"
                    f"如确认需要，请以 generate_if_missing=True 重试。",
                )
            if detect_blocked_generation(stdout):
                raise WrapperRunError(
                    "iso", 1, f"k 点 {k_point} 的参数未正确指定，无法计算子群。"
                )

        subgroups: list[SubgroupInfo] = []
        for i, row in enumerate(rows):
            subgroups.append(SubgroupInfo(
                index=start_index + i,
                space_group_number=row["space_group_number"],
                space_group_symbol=row["space_group_symbol"],
                subgroup_index=row["subgroup_index"],
                size=row["size"],
                is_maximal=row["is_maximal"],
                opd_symbol=row["opd_symbol"],
                opd_vector=row["opd_vector"],
                basis_vectors=row["basis_vectors"],
                origin=row["origin"],
                k_point_label=k_point,
                irrep_label=irrep_label,
                k_parameters=list(k_parameters) if k_parameters else [],
            ))
        return subgroups

    def enumerate_all_special_subgroups(self, parent_sg: int,
                                        distortion_types=None) -> list[SubgroupInfo]:
        """
        枚举母相的全部特殊 k 点 × 全部 IR 的各向同性子群。

        对应官网 Method 1（“Search over all special k points”）的完整候选集。
        只包含无自由参数的特殊 k 点（有参数 k 点属于 Method 2 范畴）。

        实现：每个特殊 k 点一个 iso 会话，用 ``DISPLAY SETTING`` 的输出行
        （“Current setting is International ...”）作为各 IR 子群表之间的
        分隔标记，大幅减少 WSL 进程启动次数。

        Args:
            parent_sg: 母相空间群号
            distortion_types: 畸变类型（保留参数；类型过滤在模式计算阶段执行，
                见 README 已知差异说明）

        Returns:
            List[SubgroupInfo]
        """
        _ = distortion_types
        subgroups: list[SubgroupInfo] = []
        kpoints = self.list_k_points(parent_sg)
        for kp in kpoints:
            if not kp.is_special:
                continue
            try:
                irreps = self.list_irreps(parent_sg, kp.label)
            except WrapperRunError:
                continue
            if not irreps:
                continue

            commands = [
                f"VALUE PARENT {parent_sg}",
                f"VALUE KPOINT {kp.label}",
                "SHOW SUBGROUP", "SHOW INDEX", "SHOW SIZE",
                "SHOW DIRECTION VEC", "SHOW BASIS", "SHOW ORIGIN", "SHOW MAXIMAL",
            ]
            for ir in irreps:
                commands.append(f"VALUE IRREP {ir.label}")
                commands.append("DISPLAY ISOTROPY")
                commands.append("DISPLAY SETTING")
            stdout = self._run_session("\n".join(commands) + "\n")

            # 分隔符：每个 IR 的子群表后紧跟 DISPLAY SETTING 输出
            chunks = stdout.split("Current setting is International")
            # chunks[0] = 程序横幅（含起始的 setting 行之前部分）；
            # chunks[1:] 依次对应 irreps[0], irreps[1], ...
            for ir, chunk in zip(irreps, chunks[1:], strict=False):
                try:
                    rows = parse_subgroup_table(chunk)
                except (ValueError, IndexError) as exc:
                    raise OutputParseError(
                        "iso", f"解析 {kp.label}/{ir.label} 子群表失败: {exc}"
                    ) from exc
                for row in rows:
                    subgroups.append(SubgroupInfo(
                        index=len(subgroups),
                        space_group_number=row["space_group_number"],
                        space_group_symbol=row["space_group_symbol"],
                        subgroup_index=row["subgroup_index"],
                        size=row["size"],
                        is_maximal=row["is_maximal"],
                        opd_symbol=row["opd_symbol"],
                        opd_vector=row["opd_vector"],
                        basis_vectors=row["basis_vectors"],
                        origin=row["origin"],
                        k_point_label=kp.label,
                        irrep_label=ir.label,
                    ))
        return subgroups

    # ================================================================
    # 模式基矢（Method 2 的 Distortion Page 数据源）
    # ================================================================

    def calc_distortion_modes(self, parent_sg: int,
                              subgroup: SubgroupInfo,
                              wyckoff_letters: Sequence[str]) -> list[DistortionMode]:
        """
        计算指定子群路径下的畸变模式基矢（DISPLAY BUSH + SHOW MODES）。

        Args:
            parent_sg: 母相空间群号
            subgroup: 目标子群（须含 k_point_label / irrep_label / opd_symbol）
            wyckoff_letters: 母相结构中各 Wyckoff 位置字母（来自 findsym）

        Returns:
            List[DistortionMode]：每个模式含 BushMode 原子位移基矢
        """
        if not wyckoff_letters:
            raise WrapperRunError("iso", 1, "未提供任何 Wyckoff 位置，无法计算模式。")

        if subgroup.k_parameters:
            # iso 的 DISPLAY BUSH 仅支持对称 k 点（“Selected irrep must belong
            # to a k point of symmetry”）；参数 k 点（如 LD/DT）的模式计算
            # 依赖官网的 (3+d) 维超空间机制，本地二进制无法完成。
            raise WrapperRunError(
                "iso", 1,
                f"k 点 {subgroup.k_point_label}（参数 {subgroup.k_parameters}）为参数"
                f"（非对称）k 点：iso 二进制只能枚举其子群，无法计算原子位移模式；"
                f"官网对该场景使用 (3+d) 维超空间机制，本地暂不支持。",
            )

        commands = [
            f"VALUE PARENT {parent_sg}",
            f"VALUE KPOINT {subgroup.k_point_label}",
            f"VALUE IRREP {subgroup.irrep_label}",
        ]
        if subgroup.k_parameters:
            # 带参数 k 点必须设置 KVALUE（顺序：先 IR 后 KVALUE）
            commands.append(self._kvalue_command(subgroup.k_parameters))
        commands.append(f"VALUE DIRECTION {subgroup.opd_symbol}")
        for letter in wyckoff_letters:
            commands.append(f"VALUE WYCKOFF {letter}")
        commands += [
            "SHOW MODES", "SHOW MICROSCOPIC",
            "SHOW SUBGROUP", "SHOW DIRECTION VEC", "SHOW INDEX",
            "SHOW BASIS", "SHOW ORIGIN",
            "DISPLAY BUSH",
        ]
        stdout = self._run_session("\n".join(commands) + "\n")

        if detect_blocked_generation(stdout):
            raise WrapperRunError(
                "iso", 1,
                f"k 点 {subgroup.k_point_label} 的参数未正确指定，无法计算模式；"
                f"请确认子群包含有效的 k 参数（{subgroup.k_parameters}）。",
            )

        try:
            rows = parse_bush_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析模式基矢表失败: {exc}") from exc

        # 按 IR + 序参量方向分组为 DistortionMode
        modes: dict[tuple, DistortionMode] = {}
        for row in rows:
            key = (row["irrep_label"], row["opd_symbol"])
            mode = modes.get(key)
            if mode is None:
                mode = DistortionMode(
                    irrep_label=row["irrep_label"],
                    opd_symbol=row["opd_symbol"],
                    mode_type="displacement",
                    k_point_label=subgroup.k_point_label,
                )
                modes[key] = mode
            mode.bush_modes.append(BushMode(**row))
            # 兼容旧接口：basis_vectors 取各代表原子的首个位移向量
            mode.basis_vectors = [b.displacements[0] for b in mode.bush_modes
                                  if b.displacements]
        return list(modes.values())

    # ================================================================
    # 畴变体
    # ================================================================

    def get_domains(self, parent_sg: int, subgroup: SubgroupInfo,
                    k_parameters: Sequence[str] | None = None) -> list[DomainInfo]:
        """
        获取指定子群的畴变体列表（SHOW DOMAIN）。

        畴总数等于子群在母相中的指数（subgroup_index），
        与官网 “Domains” 输出一致：每个畴含生成元、空间群、基矢与原点。

        Args:
            parent_sg: 母相空间群号
            subgroup: 目标子群（须含 k/IR/OPD）
            k_parameters: k 点参数

        Returns:
            List[DomainInfo]
        """
        commands = [
            f"VALUE PARENT {parent_sg}",
            f"VALUE KPOINT {subgroup.k_point_label}",
            f"VALUE IRREP {subgroup.irrep_label}",
        ]
        if k_parameters is None:
            k_parameters = subgroup.k_parameters
        if k_parameters:
            # 顺序：先 IR 后 KVALUE
            commands.append(self._kvalue_command(k_parameters))
        commands += [
            f"VALUE DIRECTION {subgroup.opd_symbol}",
            "SHOW DOMAIN", "SHOW DOMAIN GENERATORS",
            "SHOW SUBGROUP", "SHOW INDEX", "SHOW BASIS", "SHOW ORIGIN",
            "SHOW DIRECTION VEC",
            "DISPLAY ISOTROPY",
        ]
        stdout = self._run_session("\n".join(commands) + "\n")
        try:
            rows = parse_domain_table(stdout)
        except (ValueError, IndexError) as exc:
            raise OutputParseError("iso", f"解析畴表失败: {exc}") from exc
        return [DomainInfo(**row) for row in rows]

    # ================================================================
    # Wyckoff 位置
    # ================================================================

    def get_wyckoff_letters(self, parent_sg: int) -> list[dict]:
        """
        获取空间群的全部 Wyckoff 位置（字母 + 代表坐标）。

        Returns:
            list of dict: 每项含 wyckoff_letter, coordinates（含自由参数的
                位点坐标为 []，因为坐标由 x,y,z 参数化）
        """
        stdout = self._run_session(
            f"VALUE PARENT {parent_sg}\nSHOW WYCKOFF VECTOR\nDISPLAY PARENT\n"
        )
        # 只解析 "Wyckoff Points" 之后的区段，避免程序横幅中的
        # "International (new ed.)" 等文本被误认为位点
        section = extract_section(stdout, "Wyckoff Points")
        sites: list[dict] = []
        for line in section.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("Wyckoff"):
                continue
            # 形如: a (0,0,0), b (1/2,1/2,1/2), c (1/4,1/4,1/4), e (x,0,0), ...
            for m in re.finditer(r"([a-z])\s+\(([^()]*)\)", stripped):
                raw = [c.strip() for c in m.group(2).split(",") if c.strip()]
                coords: list[float] = []
                if len(raw) == 3 and all(
                    re.fullmatch(r"[+-]?\d+(?:/\d+)?", c) for c in raw
                ):
                    coords = [parse_fraction(c) for c in raw]
                sites.append({
                    "wyckoff_letter": m.group(1),
                    "coordinates": coords,
                })
        return sites

    # ================================================================
    # 兼容旧接口（已弃用但保留签名，抛错提示新用法）
    # ================================================================

    def get_site_splitting(self, *args, **kwargs):
        """已废弃：位点分裂由 findsym + 结构对称分析完成。"""
        raise NotImplementedError(
            "get_site_splitting 已废弃：请使用 SymmetryValidator 获取 Wyckoff 位置。"
        )

    def get_domain_operations(self, *args, **kwargs):
        """已废弃：请使用 get_domains。"""
        raise NotImplementedError("get_domain_operations 已废弃：请使用 get_domains。")
