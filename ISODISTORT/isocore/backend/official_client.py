"""从 ISODISTORT 官网获取子群（实验性，需联网）。

本地 isobyu/iso 二进制对带参数 k 点（如 LD/DT）默认不自动生成子群数据库；
当本地枚举为空时，可把相同的 (k 点, 参数) 喂给 ISODISTORT 官网
（landau3.byu.edu），复现官网「Method 2 → 不可约表示 → 序参量方向」流程，
从返回的 HTML 中解析出各向同性子群。

注意：
- 本模块是**实验性**功能：依赖官网可达、表单结构未变、以及 k 点参数编码
  与官网一致；任一不满足都会抛出明确异常（绝不静默输出错误结果）。
- 官网的 k 点/IR/OPD 约定与本地 iso 二进制可能存在版本差异；解析出的子群
  用于**浏览与选择**。即使成功，本地仍无法对带参数 k 点计算**位移模式**
  （官网用 (3+d) 维超空间机制，本地 iso 无法完成，见 README 已知差异第 5 条）。
- 若官网不可达，请改用本地「缺失时在线生成子群数据库」（
  iso 二进制自带的 Generate isotropy subgroups，离线可用）。
"""

from __future__ import annotations

import html as _html
import re
import urllib.parse
import urllib.request
from http.cookiejar import CookieJar

from ..backend import SubgroupInfo

# 官网入口（官网页面保存在 webpage_info/，URL 来自其 saved-from 元信息）
_BASE = "https://landau3.byu.edu"
_UPLOAD_URL = _BASE + "/isodistortuploadfile.php"
_FORM_URL = _BASE + "/isodistortform.php"

# 官网允许的最大响应体（防止异常超大响应）
_MAX_HTML = 10 * 1024 * 1024


def _opener() -> urllib.request.OpenerDirector:
    """带 cookie 会话的 opener（官网用会话承载母相结构状态）。"""
    return urllib.request.build_opener(urllib.request.HTTPCookieProcessor(CookieJar()))


def _post(opener, url: str, data: bytes, ctype: str) -> str:
    req = urllib.request.Request(url, data=data, headers={"Content-Type": ctype})  # noqa: S310 - 仅访问官网站点（https），无自定义 scheme
    with opener.open(req, timeout=60) as resp:
        return resp.read(_MAX_HTML).decode("utf-8", errors="replace")


def _multipart(fields: list[tuple[str, str]], files: list[tuple[str, str, str]]) -> tuple[bytes, str]:
    """构造 multipart/form-data 请求体（文件上传用）。"""
    boundary = "----ISODISTORTBOUNDARY"
    parts: list[bytes] = []
    for name, value in fields:
        parts.append(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"\r\n\r\n{value}\r\n".encode())
    for name, fname, content in files:
        parts.append(
            (f"--{boundary}\r\nContent-Disposition: form-data; name=\"{name}\"; "
             f"filename=\"{fname}\"\r\nContent-Type: application/octet-stream\r\n\r\n").encode()
            + content.encode("utf-8")
            + b"\r\n"
        )
    parts.append(f"--{boundary}--\r\n".encode())
    return b"".join(parts), f"multipart/form-data; boundary={boundary}"


def _extract_hidden(html: str) -> dict[str, str]:
    """提取表单里的隐藏字段（上传母相后官网回填的会话状态）。"""
    return {
        name: value
        for name, value in re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', html)
    }


def _html_unescape(text: str) -> str:
    return _html.unescape(text)


def _parse_opd_radio(text: str, k_point_label: str, irrep_label: str, index: int) -> SubgroupInfo | None:
    """解析官网 order parameter direction 页的 radio 文案为 SubgroupInfo。

    官网文案形如：
        P1 (a,0) 123 P4/mmm, basis={(1,0,0),(0,1,0),(0,0,6)}, origin=(0,0,0), s=12, i=12, k-active= (0,0,1/6)
    """
    m = re.match(
        r"\s*([A-Za-z0-9_\-+]+)\s*\([^)]*\)\s*(\d{1,3})\s+(\S+),\s*"
        r"basis=\{([^}]*)\},\s*origin=\(([^)]*)\),\s*s=(\d+),\s*i=(\d+),\s*k-active=\s*\(([^)]*)\)",
        text,
    )
    if not m:
        return None
    opd, sg_num, sg_sym, basis_raw, origin_raw, s_str, i_str, _kactive = m.groups()

    def _parse_vec3(raw: str) -> list[float]:
        return [_fraction(x.strip()) for x in raw.split(",") if x.strip()]

    basis = [_parse_vec3(row) for row in re.findall(r"\(([^)]*)\)", basis_raw)]
    basis = [r for r in basis if len(r) == 3]
    origin = _parse_vec3(origin_raw)
    if len(origin) != 3:
        origin = [0.0, 0.0, 0.0]

    return SubgroupInfo(
        index=index,
        space_group_number=int(sg_num),
        space_group_symbol=_html_unescape(sg_sym),
        subgroup_index=int(i_str),
        size=int(s_str),
        is_maximal=False,  # 官网文案未直接给出 maximal 标记，保守置 False
        opd_symbol=opd,
        basis_vectors=basis if len(basis) == 3 else [],
        origin=origin,
        k_point_label=k_point_label,
        irrep_label=irrep_label,
    )


def _fraction(text: str) -> float:
    """解析分数/小数坐标，如 '1/2' -> 0.5、'-3/2' -> -1.5。"""
    text = text.strip()
    try:
        if "/" in text:
            num, den = text.split("/", 1)
            return float(num) / float(den)
        return float(text)
    except (ValueError, ZeroDivisionError):
        return 0.0


def fetch_official_subgroups(parent_cif_content: str, k_point_label: str,
                             k_parameters: list[str] | None = None,
                             language: str = "en") -> list[SubgroupInfo]:
    """把 (k 点, 参数) 喂给 ISODISTORT 官网，读取其各向同性子群。

    Args:
        parent_cif_content: 母相 CIF 原文
        k_point_label: k 点标签（Miller-Love 记号，如 LD）
        k_parameters: k 点参数值（按 a/b/g 顺序，如 ["1/6"]）；带参数 k 点必须提供
        language: 界面语言（zh/en），仅用于控制台提示

    Returns:
        List[SubgroupInfo]（按官网 IR 枚举顺序）

    Raises:
        RuntimeError: 官网不可达、上传失败、或 HTML 结构变化无法解析
    """
    if not parent_cif_content.strip():
        raise RuntimeError("母相 CIF 内容为空")
    opener = _opener()

    # 1) 上传母相 CIF -> 搜索页（回填会话隐藏字段）
    body, ctype = _multipart([("input", "uploadparentcif")],
                             [("toProcess", "parent.cif", parent_cif_content)])
    html = _post(opener, _UPLOAD_URL, body, ctype)
    fields = _extract_hidden(html)
    if not fields.get("spacegroup"):
        raise RuntimeError("官网上传母相 CIF 失败：未返回搜索页（官网不可达或已改版）")

    # 2) Method 2：POST kvector -> 不可约表示页
    kp = dict(fields)
    kp.update({
        "input": "kvector",
        "kvec1": k_point_label,
        "kparam11": k_parameters[0] if k_parameters and len(k_parameters) > 0 else "",
        "kparam21": k_parameters[1] if k_parameters and len(k_parameters) > 1 else "",
        "kparam31": k_parameters[2] if k_parameters and len(k_parameters) > 2 else "",
        "nmodstar1": "0",
    })
    ir_html = _post(opener, _FORM_URL, urllib.parse.urlencode(kp).encode(),
                    "application/x-www-form-urlencoded")
    # IR 下拉选项：<option value=" 6654 *LD1, k10t1">LD1, k10t1</option>
    ir_options = re.findall(r'<option value="\s*\d+\s+\*([^"]+)">([^<]+)</option>', ir_html)
    ir_labels: list[str] = []
    for _val, label in ir_options:
        lab = _html_unescape(label).strip()
        if lab and lab not in ir_labels:
            ir_labels.append(lab)
    if not ir_labels:
        raise RuntimeError("官网未返回任何不可约表示（参数 k 点需正确提供 a/b/g 值）")

    # 3) 逐 IR POST irrep -> order parameter direction 页 -> 解析子群
    subgroups: list[SubgroupInfo] = []
    for irrep_label in ir_labels:
        ip = dict(fields)
        ip.update({
            "input": "irrep",
            "kvec1": k_point_label,
            "kvecnumber1": "  3",  # 官网回填的 k 点编号（此处用默认）
            "nmodstar1": "0",
            "irrep1": irrep_label,
            "irrepopd1": "",
        })
        opd_html = _post(opener, _FORM_URL, urllib.parse.urlencode(ip).encode(),
                         "application/x-www-form-urlencoded")
        # radio 文案：<input type="RADIO" name="orderparam" value="...">TEXT
        for m in re.finditer(r'name="orderparam"[^>]*>\s*([^<]+)</?', opd_html):
            sg = _parse_opd_radio(m.group(1), k_point_label, irrep_label, len(subgroups))
            if sg is not None:
                subgroups.append(sg)
    return subgroups
