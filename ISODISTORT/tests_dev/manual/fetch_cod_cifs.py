"""从 COD（Crystallography Open Database）批量获取测试结构 CIF。

用途：为「与官网 Isodistort 输出差异对比」提供外部真实结构源（非程序内嵌）。
只下载数据、不修改任何只读目录；产物保存到 tests_dev/cifs_external/。

下载策略：
- 按 (元素组合, 空间群) 搜索 COD（result.php），取前 N 个候选；
- 独立校验（不经过本程序）：spglib 空间群 == 预期、无部分占据
  （_atom_site_occupancy 全为 1 或缺失）、原子数合理（<= 96）；
- 首个通过者保存为 tests_dev/cifs_external/<name>.cif，来源 URL 记入 SOURCES.md。

用法：python tests_dev/manual/fetch_cod_cifs.py
"""
from __future__ import annotations

import re
import urllib.request
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

OUT_DIR = Path(__file__).resolve().parents[1] / "cifs_external"
COD_BASE = "https://www.crystallography.net/cod/"

# (name, 元素列表, 预期空间群(集合), 说明)
TARGETS = [
    ("SrTiO3_Pm3m", ["Sr", "Ti", "O"], {221}, "官网示例：R4+/M3+/GM4- 金标准"),
    ("BaTiO3_Pm3m", ["Ba", "Ti", "O"], {221}, "钙钛矿高对称母相"),
    ("WO3_Pm3m", ["W", "O"], {221}, "官网多模式耦合（M3-+R4+）候选"),
    ("LaMnO3_Pbnm", ["La", "Mn", "O"], {62}, "官网 Method 4 分解示例（正交）"),
    ("Cu_Fm3m", ["Cu"], {225}, "面心立方单质（1 原子原胞边界）"),
    ("Si_Fd3m", ["Si"], {227}, "金刚石结构（Fd-3m 227）"),
    ("NaCl_Fm3m", ["Na", "Cl"], {225}, "岩盐结构"),
    ("BaTiO3_P4mmm", ["Ba", "Ti", "O"], {123, 99}, "四方低对称母相（P4/mmm 或 P4mm）"),
    ("BiFeO3_R3c", ["Bi", "Fe", "O"], {161}, "菱方极性体系（R3c 161）"),
    ("ZnO_P63mc", ["Zn", "O"], {186}, "六方纤锌矿（P6_3mc 186）"),
    ("GaN_P63mc", ["Ga", "N"], {186}, "六方纤锌矿（P6_3mc 186）"),
    ("CsCl_Pm3m", ["Cs", "Cl"], {221}, "简单立方（CsCl 型）"),
    ("SiO2_quartz", ["Si", "O"], {152, 154}, "三方石英（P3_121/P3_221，左右手对映）"),
    ("MgO_Fm3m", ["Mg", "O"], {225}, "岩盐（二元简单）"),
    # ---- 第二批（不同 COD 批次：其它晶系/空间群/结构类型）----
    ("Fe2O3_hematite", ["Fe", "O"], {167}, "赤铁矿（R-3c 167，刚玉型）"),
    ("CaCO3_calcite", ["Ca", "C", "O"], {167}, "方解石（R-3c 167，碳酸盐）"),
    ("TiO2_rutile", ["Ti", "O"], {136}, "金红石（P4_2/mnm 136）"),
    ("TiO2_anatase", ["Ti", "O"], {141}, "锐钛矿（I4_1/amd 141）"),
    ("CaF2_fluorite", ["Ca", "F"], {225}, "萤石（Fm-3m 225）"),
    ("BaSO4_barite", ["Ba", "S", "O"], {62}, "重晶石（Pnma 62）"),
    ("MgAl2O4_spinel", ["Mg", "Al", "O"], {227}, "尖晶石（Fd-3m 227）"),
    ("KBr_Fm3m", ["K", "Br"], {225}, "岩盐型（KBr）"),
    ("Al2O3_corundum", ["Al", "O"], {167}, "刚玉（R-3c 167）"),
    ("CaTiO3_Pbnm", ["Ca", "Ti", "O"], {62}, "正交钙钛矿（Pbnm 62）"),
]


def _cod_search(els: list[str], sg: int) -> list[str]:
    """COD 元素+空间群搜索，返回 CIF 文件 ID 列表。

    COD 表单字段为 space_group_number（非 sg）；排除标记错误/重复条目。
    """
    params = [f"el{i + 1}={e}" for i, e in enumerate(els[:3])] + [
        f"space_group_number={sg}",
        "include_errors=0", "include_duplicates=0",
    ]
    url = COD_BASE + "result.php?" + "&".join(params)
    with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310 - COD 公共 HTTPS 接口
        html = resp.read().decode("utf-8", errors="replace")
    return re.findall(r'href="(\d+)\.cif', html)


def _fetch_cif(cod_id: str) -> str:
    with urllib.request.urlopen(COD_BASE + f"{cod_id}.cif", timeout=60) as resp:  # noqa: S310
        return resp.read().decode("utf-8", errors="replace")


def _validate(text: str, expected_sgs: set[int]) -> tuple[bool, str]:
    """独立校验：spglib SG ∈ 预期 且 无部分占据 且 原子数合理。"""
    tmp = Path(__file__).resolve().parents[1] / "_cod_tmp.cif"
    tmp.write_text(text, encoding="utf-8")
    try:
        struct = Structure.from_file(str(tmp))
    except Exception as exc:  # noqa: BLE001
        return False, f"解析失败: {type(exc).__name__}: {exc}"
    sg = SpacegroupAnalyzer(struct, symprec=1e-3).get_space_group_number()
    if sg not in expected_sgs:
        return False, f"spglib SG={sg} 不在预期 {sorted(expected_sgs)}"
    if len(struct) > 96:
        return False, f"原子数 {len(struct)} 过多"
    # 部分占据检查
    for site in struct:
        comp = site.species
        if len(comp) > 1 or any(abs(v - 1.0) > 1e-6 for v in comp.values()):
            return False, f"部分占据: {comp}"
    return True, f"OK SG={sg} atoms={len(struct)}"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sources: list[str] = ["# 外部测试结构来源（COD, https://www.crystallography.net/cod/）\n"]
    ok, fail = 0, 0
    for name, els, sgs, note in TARGETS:
        out = OUT_DIR / f"{name}.cif"
        if out.exists():
            print(f"skip  {name}（已存在）", flush=True)
            ok += 1
            sources.append(f"- {name}: 已缓存")
            continue
        ids: list[str] = []
        for sg in sgs:
            ids.extend(_cod_search(els, sg))
        # 去重保序
        ids = list(dict.fromkeys(ids))
        print(f"search {name} ({'+'.join(els)}, sg={sorted(sgs)}): "
              f"{len(ids)} 候选", flush=True)
        picked = None
        for cod_id in ids[:8]:
            try:
                text = _fetch_cif(cod_id)
            except Exception:  # noqa: BLE001,S112 - 单个候选下载失败继续尝试下一个
                continue
            good, msg = _validate(text, sgs)
            if good:
                out.write_text(text, encoding="utf-8")
                picked = cod_id
                print(f"  -> {cod_id}.cif  {msg}  [{note}]", flush=True)
                sources.append(
                    f"- {name}: {COD_BASE}{cod_id}.cif （{msg}；{note}）"
                )
                ok += 1
                break
            else:
                print(f"  x  {cod_id}.cif  {msg}", flush=True)
        if picked is None:
            print(f"  !! {name}: 未找到合格结构", flush=True)
            fail += 1
    (OUT_DIR / "SOURCES.md").write_text("\n".join(sources), encoding="utf-8")
    print(f"\n=== 完成: 成功 {ok} / 失败 {fail} ===", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
