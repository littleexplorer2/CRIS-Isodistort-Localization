"""
文本解析工具 - 用于解析 iso / findsym 的非结构化文本输出
"""
import re
from typing import List, Dict, Tuple, Optional


def extract_section(text: str, start_marker: str, end_marker: str = None) -> str:
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


def parse_floats(line: str) -> List[float]:
    """从一行文本中提取所有浮点数"""
    return [float(x) for x in re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", line)]


def parse_space_group_number(text: str) -> Optional[int]:
    """从 findsym 输出中提取空间群号"""
    m = re.search(r"Space Group\s+(\d+)", text)
    if m:
        return int(m.group(1))
    return None


def parse_wyckoff_sites(text: str) -> List[Dict]:
    """
    从 findsym 输出中解析 Wyckoff 位点信息

    Returns:
        list of dict: 每个位点包含 label, multiplicity, letter, position, occupancy
    """
    sites = []
    lines = text.splitlines()
    in_section = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Wyckoff position"):
            in_section = True
            # 格式: Wyckoff position a
            letter = stripped.split()[-1]
            current = {"wyckoff_letter": letter, "atoms": []}
            sites.append(current)
            continue
        if in_section and stripped and re.match(r"^\s*\d+\s+\w", stripped):
            # 格式: 1   0.00000   0.00000   0.00000     1.00000
            parts = stripped.split()
            if len(parts) >= 5:
                atom_idx = int(parts[0])
                pos = [float(parts[1]), float(parts[2]), float(parts[3])]
                occ = float(parts[4])
                current["atoms"].append({
                    "index": atom_idx,
                    "position": pos,
                    "occupancy": occ
                })
        if in_section and stripped.startswith("---"):
            break
    return sites


def parse_subgroup_list(text: str) -> List[Dict]:
    """
    从 iso 输出中解析子群列表

    Returns:
        list of dict: 每个子群包含 number, symbol, index, origin
    """
    subgroups = []
    lines = text.splitlines()
    in_table = False
    for line in lines:
        stripped = line.strip()
        # 简单启发式：包含数字 + 空间群符号样式的行
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
