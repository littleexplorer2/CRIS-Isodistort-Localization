"""
文本解析工具 - 用于解析 iso / findsym 的非结构化文本输出

包含两类解析器：
1. findsym 输出解析：空间群号、Wyckoff 位置、晶格参数
2. iso 输出解析（实测自 iso 9.6.1）：
   - k 点列表（DISPLAY KPOINT）
   - 不可约表示列表（DISPLAY IRREP）
   - 各向同性子群表（DISPLAY ISOTROPY）
   - 畴表（DISPLAY ISOTROPY + SHOW DOMAIN）
   - 模式基矢表（DISPLAY BUSH + SHOW MODES）

注意：解析前应通过 ``SCREEN 200`` 避免 80 列折行，本模块假设输入为
单行完整表格。
"""
import re
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from typing import Any


def extract_section(text: str, start_marker: str, end_marker: str | None = None) -> str:
    """
    从文本中提取两个标记之间的内容段

    Args:
        text: 原始文本
        start_marker: 起始标记（包含该行）
        end_marker: 结束标记（不包含该行），为 None 则提取到末尾
    """
    lines = text.splitlines()
    start_idx = None
    for i, line in enumerate(lines):
        if start_marker in line:
            start_idx = i
            break
    if start_idx is None:
        return ""

    if end_marker is None:
        return "\n".join(lines[start_idx:])

    for j in range(start_idx + 1, len(lines)):
        if end_marker in lines[j]:
            return "\n".join(lines[start_idx:j])
    return "\n".join(lines[start_idx:])


def parse_floats(line: str) -> list[float]:
    """从一行文本中提取所有浮点数"""
    return [float(x) for x in re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", line)]


# ================================================================
# findsym 输出解析
# ================================================================

def parse_space_group_number(text: str) -> int | None:
    """从 findsym 输出中提取空间群号（兼容 v6 与 v7 格式）。"""
    m = re.search(r"Space Group:?\s+(\d+)", text)
    if m:
        return int(m.group(1))
    return None


def parse_wyckoff_sites(text: str) -> list[dict]:
    """
    从 findsym 输出中解析 Wyckoff 位置信息。

    v6 格式：``Wyckoff position a``
    v7 格式：``Wyckoff position e (Na1), z = 0.12345``

    Returns:
        list of dict: 每个位点包含 wyckoff_letter, output_label,
        parameters, atoms
    """
    sites = []
    lines = text.splitlines()
    in_section = False
    current: dict | None = None
    for line in lines:
        stripped = line.strip()
        m = re.match(
            r"^Wyckoff position\s+([a-z])"
            r"(?:\s+\(([^)]*)\))?(.*)$",
            stripped,
        )
        if m:
            in_section = True
            parameters = {
                name: float(value)
                for name, value in re.findall(
                    r"\b([xyz])\s*=\s*"
                    r"([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)",
                    m.group(3),
                )
            }
            current = {
                "wyckoff_letter": m.group(1),
                "output_label": (m.group(2) or "").strip(),
                "parameters": parameters,
                "atoms": [],
            }
            sites.append(current)
            continue
        if in_section and current is not None and re.match(r"^\s*\d+\s+\w", stripped):
            # 格式: 1   0.00000   0.00000   0.00000     1.00000
            parts = stripped.split()
            if len(parts) >= 5:
                atom_idx = int(parts[0])
                pos = [float(parts[1]), float(parts[2]), float(parts[3])]
                occ = float(parts[4])
                current["atoms"].append({
                    "index": atom_idx,
                    "position": pos,
                    "occupancy": occ,
                })
        if in_section and stripped.startswith("---"):
            break
    return sites


def parse_subgroup_list(text: str) -> list[dict]:
    """
    从 iso 输出中解析子群列表（兼容旧启发式，保留用于历史数据）。

    Returns:
        list of dict: 每个子群包含 index, space_group_number, space_group_symbol
    """
    subgroups = []
    lines = text.splitlines()
    in_table = False
    for line in lines:
        stripped = line.strip()
        if re.match(r"^\s*\d+\s+\d+\s+[A-Z]", stripped):
            in_table = True
            parts = stripped.split()
            if len(parts) >= 3:
                subgroups.append({
                    "index": int(parts[0]),
                    "space_group_number": int(parts[1]),
                    "space_group_symbol": parts[2],
                })
        elif in_table and stripped == "":
            break
    return subgroups


# ================================================================
# iso 输出解析（实测格式）
# ================================================================

def parse_fraction(token: str) -> float:
    """解析单个数值 token，支持分数（如 '1/2'）与负分数。"""
    token = token.strip()
    if not token:
        raise ValueError("empty numeric token")
    if "/" in token:
        num, den = token.split("/", 1)
        return float(num) / float(den)
    return float(token)


def parse_vector_token(token: str, letter_value: float = 1.0) -> list[float]:
    """
    解析形如 ``(a,0,0)``、``(1/2,-1/2,0)``、``(a,0.577a)`` 的向量 token。

    参数中出现的字母（a/b/c/g 等 k 参数占位）按 ``letter_value`` 处理，
    系数支持整数、分数与小数（如 0.577a → 0.577×letter_value）。
    """
    token = token.strip()
    inner = token.strip("()").strip()
    if not inner:
        return []
    out = []
    for raw_part in inner.split(","):
        part = raw_part.strip()
        if not part:
            continue
        # 支持: 1/2, -1/2, 0.577, a, -a, 2a, -1/2a, 0.577a
        m = re.fullmatch(r"([+-]?(?:\d+(?:\.\d+)?(?:/\d+)?)?)([a-zA-Z]*)", part)
        if m and (m.group(1) or m.group(2)):
            sign = -1.0 if m.group(1).startswith("-") else 1.0
            coeff_text = m.group(1).lstrip("+-")
            if coeff_text:
                out.append(sign * parse_fraction(coeff_text) * letter_value)
            else:
                out.append(sign * letter_value)
    return out


def parse_coords_token(token: str) -> list[float]:
    """解析形如 ``(0,1/2,-1/2)`` 的坐标 token（仅数字与分数）。

    BUSH 输出中带自由参数的位点（如 ``(0,0,z)``）按参数=0 处理，
    保证解析不中断（该行仅用于代表位点匹配，参数值不影响方向模式）。
    """
    token = token.strip()
    inner = token.strip("()").strip()
    if not inner:
        return []
    out = []
    for raw_part in inner.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if re.search(r"[a-zA-Z]", part):
            out.append(0.0)  # 自由参数（如 x/y/z）按 0 处理
        else:
            out.append(parse_fraction(part))
    return out


def parse_basis_token(token: str) -> list[list[float]]:
    """解析形如 ``(1,0,0),(0,1,0),(0,0,1)`` 的 3x3 基矢 token。"""
    rows: list[list[float]] = []
    for m in re.finditer(r"\(([^()]*)\)", token):
        row = parse_coords_token(m.group(0))
        if row:
            rows.append(row)
    return rows


# k 点标签：1-3 个字母（GM/DT/LD/SM/L/X/W/Q/V/C/A/GP/Z/M/R/N/Y/T/S/H/U 等）
# k 矢量分量可含分数、字母参数与线性组合，如 "1/2"、"2a"、"-2a+1"、"a"
_KOORD_COMPONENT_RE = re.compile(r"^[+-]?[0-9a-zA-Z/+-]+$")


def _is_kpoint_label(token: str) -> bool:
    """判断 token 是否像 k 点标签（1-3 个字母）。"""
    return bool(re.fullmatch(r"[A-Za-z]{1,3}", token))


def parse_kpoint_table(text: str) -> list[dict[str, Any]]:
    """
    解析 DISPLAY KPOINT 输出。

    格式::

        GM  (0,0,0)
        DT  (0,2a,0)
        LD  (a,a,a)

    Returns:
        list of dict: 每个 k 点包含 label, coordinates(字符串列表), parameters(字母列表)
    """
    kpoints: list[dict[str, Any]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or not _is_kpoint_label(stripped.split()[0]):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        label = parts[0]
        coords_raw = "".join(parts[1:]).strip("()")
        coords = [c for c in coords_raw.split(",") if c]
        # 严格校验：恰好 3 个坐标分量，且每个分量符合 k 矢量格式
        if len(coords) != 3 or not all(_KOORD_COMPONENT_RE.match(c) for c in coords):
            continue
        parameters = sorted({c for c in coords if re.search(r"[a-zA-Z]", c)})
        kpoints.append({
            "label": label,
            "coordinates": coords,
            "parameters": parameters,
            "is_special": not parameters,
        })
    return kpoints


# Miller-Love 不可约表示标签，如 GM1+ / DT1 / Z5- / mGM1+ / A2A3（配对标签）
_IRREP_LABEL_RE = re.compile(
    r"^m?[A-Z]{1,3}\d{1,2}[+-]?(?:[A-Z]{1,3}\d{1,2}[+-]?)?$"
)


def parse_irrep_table(text: str) -> list[dict[str, Any]]:
    """
    解析 DISPLAY IRREP 输出。

    两种格式（取决于是否显示 k 矢量列）::

        GM1+       1    no
        DT1        (0,2a,0) 6

    Returns:
        list of dict: 每个 IR 包含 label, dimension, active
    """
    irreps: list[dict[str, Any]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        label = parts[0]
        if not _IRREP_LABEL_RE.match(label):
            continue
        dim = None
        active = True
        # 从右侧解析维度与 active 标记（数值或 yes/no）
        for tok in reversed(parts[1:]):
            if re.fullmatch(r"\d+", tok):
                dim = int(tok)
                break
            if tok.lower() in {"yes", "no"}:
                active = tok.lower() == "yes"
        if dim is None or not (1 <= dim <= 12):
            continue
        irreps.append({"label": label, "dimension": dim, "active": active})
    return irreps


# 子群表（DISPLAY ISOTROPY，SCREEN 200 单行格式）
# 注意：数据库查询的子群表带 Max 列（yes/no）；
# 在线生成（非特殊 k 点）的子群表不带 Max 数据，故 Max 列可选
_SUBGROUP_ROW_RE = re.compile(
    r"^\s*(?P<sg>\d+)\s+(?P<symbol>\S+)\s+(?:(?P<max>yes|no)\s+)?"
    r"(?P<index>\d+)\s+(?P<size>\d+)\s+(?P<dir>\S+)\s+"
    r"(?P<dirvec>\([^()]*\))\s+(?P<basis>\([^()]*\)(?:,\([^()]*\)){2})\s+"
    r"(?P<origin>\([^()]*\))\s*$"
)


def parse_subgroup_table(text: str) -> list[dict[str, Any]]:
    """
    解析 DISPLAY ISOTROPY 的子群表（SHOW SUBGROUP/INDEX/SIZE/DIRECTION/BASIS/ORIGIN/MAXIMAL）。

    两种格式（取决于数据来源）::

        # 特殊 k 点（数据库查询，带 Max 列）
        Subgroup  Max Index Size Dir         Basis Vectors                Origin
        121 I-42m yes 6     1    P1  (a,0,0) (-1/2,1/2,0),(-1/2,-1/2,0),(0,0,1) (0,0,0)

        # 非特殊 k 点（在线生成，无 Max 数据）
        Subgroup   Max Index Size Dir            Basis Vectors           Origin
        123 P4/mmm     6     6    P1  (a,0)      (1,0,0),(0,1,0),(0,0,3) (0,0,0)

    Returns:
        list of dict，每项含:
        space_group_number, space_group_symbol, is_maximal,
        subgroup_index, size, opd_symbol, opd_vector, opd_dir_raw,
        basis_vectors, basis_raw, origin, origin_raw
    """
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        m = _SUBGROUP_ROW_RE.match(line)
        if not m:
            continue
        dirvec_raw = m.group("dirvec")
        basis_raw = m.group("basis")
        origin_raw = m.group("origin")
        rows.append({
            "space_group_number": int(m.group("sg")),
            "space_group_symbol": m.group("symbol"),
            "is_maximal": (m.group("max") or "no") == "yes",
            "subgroup_index": int(m.group("index")),
            "size": int(m.group("size")),
            "opd_symbol": m.group("dir"),
            "opd_vector": parse_vector_token(dirvec_raw.replace(";", ",")),
            "opd_dir_raw": dirvec_raw,
            "basis_vectors": parse_basis_token(basis_raw),
            "basis_raw": basis_raw,
            "origin": parse_coords_token(origin_raw),
            "origin_raw": origin_raw,
        })
    return rows


# 畴表（DISPLAY ISOTROPY + SHOW DOMAIN）
# 第一行携带 index 与 dir 列，后续行这两列为空
_DOMAIN_ROW_RE = re.compile(
    r"^\s*(?P<domain>\d+)\s+(?P<gen>\([^()]*\))\s+(?P<sg>\d+)\s+(?P<symbol>\S+)"
    r"(?:\s+(?P<index>\d+))?\s+(?:(?P<dir>\S+)\s+)?"
    r"(?P<dirvec>\([^()]*\))\s+"
    r"(?P<basis>\([^()]*\)(?:,\([^()]*\)){2})\s+(?P<origin>\([^()]*\))\s*$"
)


def parse_domain_table(text: str) -> list[dict[str, Any]]:
    """
    解析畴表（DISPLAY ISOTROPY + SHOW DOMAIN + SHOW DOMAIN GENERATORS）。

    格式::

        Domain Gen          Subgroup Index Dir          Basis Vectors       Origin
        1      (E|0,0,0)    107 I4mm 6     P1  (a,0,0) (0,1/2,-1/2),...     (0,0,0)
        2      (C2y|0,0,0)  107 I4mm           (-a,0,0) ...                 (0,0,0)

    第一行携带 index 与 dir，后续行两者为空。
    """
    rows: list[dict[str, Any]] = []
    fallback_index: int | None = None
    fallback_dir: str | None = None
    fallback_dirvec: list[float] = []
    for line in text.splitlines():
        m = _DOMAIN_ROW_RE.match(line)
        if not m:
            continue
        index = int(m.group("index")) if m.group("index") else fallback_index
        dir_symbol = m.group("dir") or fallback_dir or ""
        dirvec = parse_vector_token(m.group("dirvec")) if m.group("dirvec") else fallback_dirvec
        if m.group("index"):
            fallback_index = index
            fallback_dir = dir_symbol
            fallback_dirvec = dirvec
        rows.append({
            "domain_number": int(m.group("domain")),
            "generator": m.group("gen"),
            "space_group_number": int(m.group("sg")),
            "space_group_symbol": m.group("symbol"),
            "subgroup_index": index,
            "opd_symbol": dir_symbol,
            "opd_vector": dirvec,
            "basis_vectors": parse_basis_token(m.group("basis")),
            "origin": parse_coords_token(m.group("origin")),
        })
    return rows


# 模式基矢表（DISPLAY BUSH + SHOW MODES）
# 首行带 irrep/dir/wyckoff 列，后续行为续行（仅 Point + Displacement），
# Displacement 列可含多个向量（高维模式），如 "(0,1,-1), (0,1,1)"
_BUSH_HEADER_RE = re.compile(
    r"^\s*(?P<irrep>\S+)\s+(?P<dir>\S+)\s+(?P<wyckoff>[a-z])\s+"
    r"\((?P<point>[^()]*)\)\s+(?P<disp>.*)$"
)
_BUSH_CONT_RE = re.compile(
    r"^\s+\((?P<point>[^()]*)\)\s+(?P<disp>.*)$"
)


def _parse_disp_tokens(disp_text: str) -> list[list[float]]:
    """解析位移列中的一组向量，如 "(0,1,-1), (0,1,1)"。"""
    vectors: list[list[float]] = []
    for m in re.finditer(r"\(([^()]*)\)", disp_text):
        vec = parse_coords_token(m.group(0))
        if vec:
            vectors.append(vec)
    return vectors


def parse_bush_table(text: str) -> list[dict[str, Any]]:
    """
    解析 DISPLAY BUSH + SHOW MODES 的模式基矢表。

    格式::

        Irrep (ML) Dir(dom) Wyckoff Point            Displacement
        GM4-       P1(1)    c       (1/4,1/4,1/4)    (1,0,0)
                                   (-1/4,-1/4,-1/4) (1,0,0)
        GM5-       P1(1)    d       (0,1/4,1/4)      (0,1,-1), (0,1,1)

    带自由参数的位点（如 ``(0,0,z)``）会同时保存原始 token（point_raw），
    供模式映射器用结构真实坐标解析；数值版 point 中参数按 0 处理。

    Returns:
        list of dict，每项（一个 Wyckoff 位置的一个代表原子）含:
        irrep_label, opd_symbol, wyckoff_letter, point, point_raw,
        displacements（可能含多个向量）
    """
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("Irrep"):
            continue
        m = _BUSH_HEADER_RE.match(line)
        if m:
            current = {
                "irrep_label": m.group("irrep"),
                "opd_symbol": m.group("dir"),
                "wyckoff_letter": m.group("wyckoff"),
                "point": parse_coords_token(m.group("point")),
                "point_raw": m.group("point").strip("()").split(","),
                "displacements": _parse_disp_tokens(m.group("disp")),
            }
            rows.append(current)
            continue
        m = _BUSH_CONT_RE.match(line)
        if m and current is not None:
            rows.append({
                "irrep_label": current["irrep_label"],
                "opd_symbol": current["opd_symbol"],
                "wyckoff_letter": current["wyckoff_letter"],
                "point": parse_coords_token(m.group("point")),
                "point_raw": m.group("point").strip("()").split(","),
                "displacements": _parse_disp_tokens(m.group("disp")),
            })
    return rows


# Exact child-embedding direction table (DISPLAY DIRECTION).
_INVARIANT_DIRECTION_WITH_K_RE = re.compile(
    r"^\s*(?P<irrep>\S+)\s+"
    r"(?:(?P<k_parameters>(?!\()[^\s]+)\s+)?"
    r"(?P<k_coordinates>\([^()]*\))\s+"
    r"(?P<direction>\([^()]*\))\s+"
    r"(?P<sg>\d+)\s+(?P<symbol>\S+)\s+(?P<size>\d+)\s*$"
)
_INVARIANT_DIRECTION_RE = re.compile(
    r"^\s*(?P<irrep>\S+)\s+(?P<direction>\([^()]*\))\s+"
    r"(?P<sg>\d+)\s+(?P<symbol>\S+)\s+(?P<size>\d+)\s*$"
)


def parse_invariant_direction_table(text: str) -> list[dict[str, Any]]:
    """Parse the parent irreps fixed by an exact child embedding.

    The command sequence is ``VALUE SUBGROUP/BASIS/ORIGIN`` followed by
    ``SHOW IRREP``, ``SHOW DIRECTION VEC``, ``DISPLAY DIRECTION``.  Direction
    coordinates are retained as text because semicolons distinguish star-arm
    groupings in ISO notation.
    """

    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        match = _INVARIANT_DIRECTION_WITH_K_RE.match(line)
        if match is None:
            match = _INVARIANT_DIRECTION_RE.match(line)
        if match is None:
            continue
        groups = match.groupdict()
        k_parameters = tuple(
            value.strip()
            for value in str(groups.get("k_parameters") or "").split(",")
            if value.strip()
        )
        k_coordinates = tuple(
            value.strip()
            for value in str(groups.get("k_coordinates") or "").strip("()").split(",")
            if value.strip()
        )
        row = {
            "irrep_label": match.group("irrep"),
            "direction_raw": match.group("direction"),
            "subgroup_number": int(match.group("sg")),
            "subgroup_symbol": match.group("symbol"),
            "size": int(match.group("size")),
        }
        if k_parameters or k_coordinates:
            row["k_parameters"] = k_parameters
            row["k_coordinates"] = k_coordinates
        rows.append(row)
    return rows


# Macroscopic tensor basis (VALUE RANK ... + SHOW MACROSCOPIC).
_ENGINEERING_STRAIN_BASIS = ("xx", "yy", "zz", "yz", "xz", "xy")


def parse_symbolic_tensor_expression(
    expression: str,
    coefficient_basis: tuple[str, ...] = _ENGINEERING_STRAIN_BASIS,
) -> dict[str, Any]:
    """Resolve one ISO tensor expression without inventing exact radicals.

    Integer and rational coefficients are retained as :class:`Fraction`.
    A finite decimal such as ISO's printed ``1.732`` is retained as
    :class:`Decimal` and marked ``printed_decimal``; it is deliberately not
    treated as an exact representation of ``sqrt(3)``.  Unsupported symbolic
    functions are preserved verbatim and marked ``unresolved``.
    """

    raw = str(expression).strip()
    compact = re.sub(r"\s+", "", raw)
    if not compact:
        return {
            "expression_raw": raw,
            "coefficient_basis": tuple(coefficient_basis),
            "coefficients": None,
            "quality": "unresolved",
            "reason": "empty_tensor_expression",
        }
    symbols = sorted(set(coefficient_basis), key=len, reverse=True)
    symbol_pattern = "|".join(re.escape(symbol) for symbol in symbols)
    number = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
    term_re = re.compile(
        rf"(?P<sign>[+-]?)(?P<number>{number}(?:/{number})?)?"
        rf"(?P<symbol>{symbol_pattern})"
    )

    values: dict[str, Fraction | Decimal] = {
        symbol: Fraction(0) for symbol in coefficient_basis
    }
    saw_decimal = False
    offset = 0
    while offset < len(compact):
        match = term_re.match(compact, offset)
        if match is None:
            return {
                "expression_raw": raw,
                "coefficient_basis": tuple(coefficient_basis),
                "coefficients": None,
                "quality": "unresolved",
                "reason": "unsupported_tensor_expression",
            }
        token = match.group("number") or "1"
        sign = -1 if match.group("sign") == "-" else 1
        try:
            if "." in token or "e" in token.lower():
                if "/" in token:
                    numerator, denominator = token.split("/", 1)
                    coefficient: Fraction | Decimal = (
                        Decimal(numerator) / Decimal(denominator)
                    )
                else:
                    coefficient = Decimal(token)
                saw_decimal = True
            else:
                coefficient = Fraction(token)
        except (InvalidOperation, ValueError, ZeroDivisionError):
            return {
                "expression_raw": raw,
                "coefficient_basis": tuple(coefficient_basis),
                "coefficients": None,
                "quality": "unresolved",
                "reason": "invalid_tensor_coefficient",
            }
        coefficient *= sign
        symbol = match.group("symbol")
        previous = values[symbol]
        if isinstance(previous, Decimal) or isinstance(coefficient, Decimal):
            if not isinstance(previous, Decimal):
                previous = Decimal(previous.numerator) / Decimal(previous.denominator)
            if not isinstance(coefficient, Decimal):
                coefficient = (
                    Decimal(coefficient.numerator) / Decimal(coefficient.denominator)
                )
            values[symbol] = previous + coefficient
            saw_decimal = True
        else:
            values[symbol] = previous + coefficient
        offset = match.end()

    if saw_decimal:
        coefficients = tuple(
            value
            if isinstance(value, Decimal)
            else Decimal(value.numerator) / Decimal(value.denominator)
            for value in (values[symbol] for symbol in coefficient_basis)
        )
        quality = "printed_decimal"
    else:
        coefficients = tuple(values[symbol] for symbol in coefficient_basis)
        quality = "exact_rational"
    return {
        "expression_raw": raw,
        "coefficient_basis": tuple(coefficient_basis),
        "coefficients": coefficients,
        "quality": quality,
        "reason": None,
    }


_MACROSCOPIC_ROW_RE = re.compile(
    r"^\s*(?P<irrep>\S+)\s{2,}(?P<expressions>\S.*)\s*$"
)


def parse_macroscopic_tensor_blocks(
    text: str,
    coefficient_basis: tuple[str, ...] = _ENGINEERING_STRAIN_BASIS,
) -> list[dict[str, Any]]:
    """Parse ordered irrep-copy rows from ``SHOW MACROSCOPIC`` output.

    Comma-separated expressions on one printed row are components of one
    physically irreducible real representation.  An indented continuation is
    a new copy of the preceding irrep.  Combined physical labels such as
    ``GM2GM3`` are kept intact.
    """

    blocks: list[dict[str, Any]] = []
    active = False
    current_irrep: str | None = None
    for line in text.splitlines():
        if "Irrep (ML) Basis Functions" in line:
            active = True
            current_irrep = None
            continue
        if not active:
            continue
        stripped = line.strip()
        if stripped == "*":
            active = False
            current_irrep = None
            continue
        if not stripped:
            continue

        match = _MACROSCOPIC_ROW_RE.match(line)
        if match is not None:
            current_irrep = match.group("irrep")
            payload = match.group("expressions")
        elif line[:1].isspace() and current_irrep is not None:
            payload = stripped
        else:
            raise ValueError(f"malformed macroscopic tensor row {line!r}")
        expressions = [item.strip() for item in payload.split(",")]
        if not expressions or any(not item for item in expressions):
            raise ValueError("macroscopic tensor row contains an empty component")
        blocks.append({
            "global_irrep": current_irrep,
            "source_order": len(blocks),
            "components": [
                parse_symbolic_tensor_expression(item, coefficient_basis)
                for item in expressions
            ],
        })
    return blocks


_MICRO_HEADER_WITH_DIRECTION_RE = re.compile(
    r"^\s*(?P<irrep>\S+)\s+(?P<direction>\S+)\s+"
    r"(?P<wyckoff>[a-z])\s+(?P<site_irrep>\S+)\s+"
    r"\((?P<point>[^()]*)\)\s+(?P<disp>.*)$"
)
_MICRO_HEADER_RE = re.compile(
    r"^\s*(?P<irrep>\S+)\s+(?P<wyckoff>[a-z])\s+"
    r"(?P<site_irrep>\S+)\s+\((?P<point>[^()]*)\)\s+(?P<disp>.*)$"
)
_MICRO_CONT_RE = re.compile(
    r"^\s+\((?P<point>[^()]*)\)\s+(?P<disp>.*)$"
)


def _printed_half_step(token: str) -> Fraction:
    """Return the rounding half-width encoded by one numeric source token."""

    value = token.strip()
    if "/" in value or ("." not in value and "e" not in value.lower()):
        return Fraction(0)
    decimal = Decimal(value)
    return Fraction(Decimal("0.5").scaleb(decimal.as_tuple().exponent))


def _exact_vector_tokens(
    text: str,
) -> tuple[
    tuple[tuple[Fraction, Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction, Fraction], ...],
]:
    vectors: list[tuple[Fraction, Fraction, Fraction]] = []
    half_steps: list[tuple[Fraction, Fraction, Fraction]] = []
    for match in re.finditer(r"\(([^()]*)\)", text):
        tokens = tuple(part.strip() for part in match.group(1).split(","))
        if len(tokens) != 3:
            raise ValueError("microscopic displacement vector must have three components")
        try:
            vector = tuple(Fraction(token) for token in tokens)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(
                f"non-rational microscopic displacement vector {match.group(0)!r}"
            ) from exc
        vectors.append(vector)  # type: ignore[arg-type]
        half_steps.append(tuple(_printed_half_step(token) for token in tokens))
    return tuple(vectors), tuple(half_steps)  # type: ignore[return-value]


def parse_microscopic_vector_blocks(text: str) -> list[dict[str, Any]]:
    """Parse ordered local-irrep blocks from ``DISPLAY DISTORTION``.

    ISO omits the ``Dir`` column when the calculation is conditioned directly
    on ``VALUE SUBGROUP/BASIS/ORIGIN``.  Both that form and the older explicit
    ``VALUE DIRECTION`` form are accepted.  Every rational vector component is
    retained as :class:`Fraction`; point expressions such as ``-x+1/2`` stay
    as tokens for later substitution against a concrete physical orbit.
    """

    blocks: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("Irrep") or stripped == "*":
            continue
        match = _MICRO_HEADER_WITH_DIRECTION_RE.match(line)
        if match is None:
            match = _MICRO_HEADER_RE.match(line)
        if match is not None:
            groups = match.groupdict()
            point = tuple(part.strip() for part in groups["point"].split(","))
            if len(point) != 3:
                raise ValueError("microscopic representative point must have three components")
            vectors, half_steps = _exact_vector_tokens(groups["disp"])
            if not vectors:
                raise ValueError("microscopic-vector block contains no displacement columns")
            current = {
                "global_irrep": groups["irrep"],
                "direction_symbol": groups.get("direction") or "",
                "wyckoff_letter": groups["wyckoff"],
                "site_irrep": groups["site_irrep"],
                "source_order": len(blocks),
                "rows": [{
                    "point_raw": point,
                    "displacements": vectors,
                    "displacement_half_steps": half_steps,
                }],
            }
            blocks.append(current)
            continue
        match = _MICRO_CONT_RE.match(line)
        if match is not None and current is not None:
            point = tuple(part.strip() for part in match.group("point").split(","))
            if len(point) != 3:
                raise ValueError("microscopic continuation point must have three components")
            vectors, half_steps = _exact_vector_tokens(match.group("disp"))
            expected = len(current["rows"][0]["displacements"])
            if len(vectors) != expected:
                raise ValueError(
                    "microscopic-vector continuation row has inconsistent column count"
                )
            current["rows"].append({
                "point_raw": point,
                "displacements": vectors,
                "displacement_half_steps": half_steps,
            })
    return blocks


_SITE_CHARACTER_ROW_RE = re.compile(
    r"^(?P<wyckoff>[a-z]+)\s+(?P<point_group>\S+)\s+(?P<body>\(.*)$"
)
_SITE_OPERATION_RE = re.compile(
    r"\((?P<label>[^|()]+)\|(?P<translation>[^()]*)\)"
)
_SITE_CHARACTER_VALUE_RE = re.compile(
    r"(?P<label>[^,\s]+)\s+(?P<value>[+-]?\d+(?:/\d+)?)\s*(?:,|$)"
)


def parse_wyckoff_character_tables(text: str) -> list[dict[str, Any]]:
    """Parse orientation-resolved ``SHOW WYCKOFF CHARACTER`` output."""

    logical_rows: list[str] = []
    current: str | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if _SITE_CHARACTER_ROW_RE.match(stripped):
            if current is not None:
                logical_rows.append(current)
            current = stripped
            continue
        # Long C6/D6 tables wrap at SCREEN 200.  A continuation begins with
        # indentation and belongs to the current Wyckoff row.  Embedded stars
        # in labels such as E1* are data, not truncation sentinels.
        if current is not None and line[:1].isspace() and stripped and stripped != "*":
            current += " " + stripped
    if current is not None:
        logical_rows.append(current)

    tables: list[dict[str, Any]] = []
    for logical in logical_rows:
        header = _SITE_CHARACTER_ROW_RE.match(logical)
        if header is None:
            continue
        body = header.group("body")
        operations = list(_SITE_OPERATION_RE.finditer(body))
        if not operations:
            raise ValueError("Wyckoff character row contains no site operations")
        operation_rows: list[dict[str, Any]] = []
        irrep_order: list[str] = []
        characters: dict[str, list[Fraction]] = {}
        for index, operation in enumerate(operations):
            translation_tokens = tuple(
                part.strip() for part in operation.group("translation").split(",")
            )
            if len(translation_tokens) != 3:
                raise ValueError("site-operation translation must have three components")
            translation = tuple(Fraction(token) for token in translation_tokens)
            operation_rows.append({
                "point_operation": operation.group("label").strip(),
                "translation": translation,
            })
            end = operations[index + 1].start() if index + 1 < len(operations) else len(body)
            payload = body[operation.end():end].strip(" ,")
            values = list(_SITE_CHARACTER_VALUE_RE.finditer(payload))
            if not values:
                raise ValueError("site operation contains no character values")
            labels = [value.group("label") for value in values]
            if index == 0:
                irrep_order = labels
                characters = {label: [] for label in labels}
            elif labels != irrep_order:
                raise ValueError("site character columns change order between operations")
            for value in values:
                characters[value.group("label")].append(Fraction(value.group("value")))
        tables.append({
            "wyckoff_letter": header.group("wyckoff"),
            "point_group": header.group("point_group"),
            "operations": operation_rows,
            "irreps": [
                {"label": label, "characters": tuple(characters[label])}
                for label in irrep_order
            ],
        })
    return tables


def detect_missing_subgroup_db(text: str) -> bool:
    """判断输出中是否出现“子群数据库缺失、需要生成”的提示。

    对应 iso 在非特殊 k 点/超胞下询问是否在线生成子群数据库的行为，
    与官网“Generate isotropy subgroups”的等待场景一致。
    """
    return "Should the data base be added?" in text


def detect_blocked_generation(text: str) -> bool:
    """判断输出是否表明“参数未选择/不可计算”等阻塞状态。"""
    return "parameters not selected for k vector" in text
