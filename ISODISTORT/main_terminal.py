"""Terminal menu for the local ISODISTORT workflow.

Search Page: load a parent CIF, set distortion types, run Method 1-4
(parameters + compute; filter/sort the result table).
Distortion: download filtered result tables (Methods 1-4) and export
subgroup structure files from Method 1, 2, or 3 (same files as the web ZIP).

Requires WSL on Windows (isobyu binaries are Linux ELF). First run creates
a short WSL staging path and an ISODATA symlink.

Method 1 enumerates isotropy subgroups over all special k points and can take
several seconds on the first query.
"""

from __future__ import annotations

import csv
import sys
from collections.abc import Callable, Sequence
from pathlib import Path

from isocore.api import IsoDistort
from isocore.distortion.search_methods import CRYSTAL_SYSTEMS
from isocore.i18n import t
from isocore.utils import IsodistortError, get_config
from isocore.utils.parent_header import format_wyckoff_sites
from isocore.utils.schoenflies import (
    POINT_GROUP_SCHOENFLIES,
    POINT_GROUP_SYSTEM,
    hm_symbol,
    schoenflies_symbol,
)

DISTORTION_TYPE_MAP = {
    1: "displacive",
    2: "occupational",
    3: "strain",
    4: "magnetic",
    5: "rotational",
}

# 旧名兼容（displacement -> displacive，order -> occupational）
_TYPE_ALIASES = {
    "displacement": "displacive",
    "order": "occupational",
}


def _line() -> None:
    """打印分隔线"""
    print("-" * 78)


def _prompt(text: str, default: str | None = None) -> str:
    """通用的读取用户输入的方法：提示用户输入，支持默认值"""
    hint = f" [{default}]" if default is not None else ""
    while True:  # 循环直到用户输入有效值（非空）
        value = input(f"{text}{hint}: ").strip()
        if value:
            return value
        if default is not None:
            return default


def _prompt_int(text: str, default: int | None = None) -> int:
    while True:
        raw = _prompt(text, str(default) if default is not None else None)
        try:
            return int(raw)
        except ValueError:
            print(t("ui.prompt.int"))


def _prompt_float(text: str, default: float | None = None) -> float:
    while True:
        raw = _prompt(text, str(default) if default is not None else None)
        try:
            return float(raw)
        except ValueError:
            print(t("ui.prompt.float"))


def _prompt_yes_no(text: str, default_yes: bool = False) -> bool:
    default = "y" if default_yes else "n"
    while True:
        raw = _prompt(f"{text} (y/n)", default).lower()
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False
        print(t("ui.prompt.yes_no"))


def _find_cif_candidates(root: Path, limit: int = 30) -> list[Path]:
    """在指定目录下查找 CIF 文件，返回前 limit 个结果，按路径排序。

    - 只扫描 ``root``（终端调用时为项目根目录）及其子目录下的 ``.cif``。
    - 自动排除配置的输出目录与临时目录（``output_dir`` / ``temp_dir``）中的
      CIF——这些多为程序生成的成品/上传暂存文件，不应出现在“输入母相 CIF”
      候选列表里。如需加载某个输出文件，请用菜单的“手动输入路径”选项。
    """
    cfg = get_config()
    exclude = {cfg.output_dir.resolve(), cfg.temp_dir.resolve()}
    found: list[Path] = []
    for path in root.rglob("*.cif"):
        if not path.is_file():
            continue
        try:
            if any(path.is_relative_to(ex) for ex in exclude):
                continue
        except ValueError:
            # 相对/绝对路径不一致时无法判定包含，按“不排除”处理
            pass
        found.append(path)
        if len(found) >= limit:
            break
    return sorted(found)


def _choose_cif(project_root: Path, title: str) -> str:
    """交互式选择 CIF 文件（列出候选或手动输入路径）。"""
    _line()
    print(title)
    candidates = _find_cif_candidates(project_root)

    if candidates:
        print(t("ui.msg.cif_found"))
        for idx, path in enumerate(candidates, start=1):
            try:
                display = str(path.relative_to(project_root))
            except ValueError:
                display = str(path)
            print(f"  {idx:2d}. {display}")
        print("   0. " + t("ui.msg.cif_manual"))

        while True:
            choice = _prompt_int(t("ui.msg.cif_choice"), 1)
            if choice == 0:
                break
            if 1 <= choice <= len(candidates):
                return str(candidates[choice - 1])
            print(t("ui.msg.cif_range"))

    return _prompt(t("ui.prompt.cif_manual"))


def _prompt_distortion_types(default: Sequence[str]) -> list[str]:
    """选择 Distortion Types（对应官网勾选框，逻辑 OR）。"""
    _line()
    print("Distortion Types")
    print("  1. displacive（" + t("ui.type.displacive") + "）")
    print("  2. occupational（" + t("ui.type.occupational") + "）")
    print("  3. strain（" + t("ui.type.strain") + "）")
    print("  4. magnetic（" + t("ui.type.magnetic") + "）")
    print("  5. rotational（" + t("ui.type.rotational") + "）")
    default_raw = ",".join(default)
    raw = _prompt(t("ui.prompt.types"), default_raw)

    labels: list[str] = []
    for token in [x.strip().lower() for x in raw.split(",") if x.strip()]:
        if token.isdigit() and int(token) in DISTORTION_TYPE_MAP:
            labels.append(DISTORTION_TYPE_MAP[int(token)])
        elif token in DISTORTION_TYPE_MAP.values() or token in _TYPE_ALIASES:
            labels.append(_TYPE_ALIASES.get(token, token))

    dedup: list[str] = []
    seen = set()
    for label in labels:
        if label not in seen:
            dedup.append(label)
            seen.add(label)

    if not dedup:
        return ["displacive", "strain"]
    return dedup


def _prompt_supercell(default_cell: Sequence[int]) -> list[int]:
    default_raw = ",".join(str(x) for x in default_cell)
    while True:
        raw = _prompt(t("ui.prompt.supercell"), default_raw)
        parts = [x.strip() for x in raw.split(",") if x.strip()]
        if len(parts) != 3:
            print("Supercell must be 3 integers")
            continue
        try:
            values = [int(x) for x in parts]
        except ValueError:
            print("Supercell values must be integers")
            continue
        if any(v < 1 for v in values):
            print("Supercell values must be >= 1")
            continue
        return values


def _prompt_basis_matrix() -> list[list[str]]:
    print("Enter a 3x3 basis matrix (3 values per row; fractions like 1/2 allowed)")
    rows: list[list[str]] = []
    for i in range(3):
        while True:
            row = _prompt(f"Row {i + 1}", "1 0 0" if i == 0 else "0 1 0" if i == 1 else "0 0 1")
            parts = [x for x in row.split() if x]
            if len(parts) != 3:
                print("Each row must have 3 values")
                continue
            rows.append(parts)
            break
    return rows


_CS_ORDER = {
    "triclinic": 0,
    "monoclinic": 1,
    "orthorhombic": 2,
    "tetragonal": 3,
    "trigonal": 4,
    "hexagonal": 5,
    "cubic": 6,
}


def _sg_text(sg) -> str:
    return f"{sg.space_group_number} {sg.space_group_symbol or ''}".strip()


def _empty_tbl(
    cols: list[tuple[str, str, Callable, Callable, bool]],
) -> dict:
    return {
        "rows": [],
        "filters": {},
        "sort": (None, 1),
        "show_only": False,
        "cols": cols,
    }


def _method1_cols() -> list[tuple[str, str, Callable, Callable, bool]]:
    return [
        ("idx", "idx", lambda r: str(r["index"]), lambda r: int(r["index"]), False),
        ("sg", "SG", lambda r: r["sg"], lambda r: r["_sort_sg"], True),
        ("k", "k", lambda r: r["k"], lambda r: r["k"], True),
        ("irrep", "Irrep", lambda r: r["irrep"], lambda r: r["irrep"], True),
        ("opd", "OPD", lambda r: r["opd"], lambda r: r["opd"], True),
        ("cs", "crystal system", lambda r: r["cs"], lambda r: r["_sort_cs"], True),
        ("max", "maximal", lambda r: r["max"], lambda r: r["_sort_max"], True),
    ]


def _method2_cols() -> list[tuple[str, str, Callable, Callable, bool]]:
    return [
        ("idx", "idx", lambda r: str(r["index"]), lambda r: int(r["index"]), False),
        ("sg", "SG", lambda r: r["sg"], lambda r: r["_sort_sg"], True),
        ("k", "k", lambda r: r["k"], lambda r: r["k"], True),
        ("irrep", "Irrep", lambda r: r["irrep"], lambda r: r["irrep"], True),
        ("opd", "OPD", lambda r: r["opd"], lambda r: r["opd"], True),
        ("s", "s", lambda r: str(r["s"]), lambda r: r["_sort_s"], True),
        ("i", "i", lambda r: str(r["i"]), lambda r: r["_sort_i"], True),
    ]


def _method3_cols() -> list[tuple[str, str, Callable, Callable, bool]]:
    return [
        ("idx", "idx", lambda r: str(r["index"]), lambda r: int(r["index"]), False),
        ("sg", "SG", lambda r: r["sg"], lambda r: r["_sort_sg"], True),
        ("k", "k", lambda r: r["k"], lambda r: r["k"], True),
        ("irrep", "Irrep", lambda r: r["irrep"], lambda r: r["irrep"], True),
        ("pg", "point group", lambda r: r["pg"], lambda r: r["pg"], True),
    ]


def _method4_cols() -> list[tuple[str, str, Callable, Callable, bool]]:
    return [
        ("mode", "mode", lambda r: r["mode"], lambda r: r["mode"], True),
        ("amp", "amplitude", lambda r: r["amp"], lambda r: r["_sort_amp"], True),
    ]


def _row_method1(item) -> dict:
    sg = item.subgroup
    return {
        "index": sg.index,
        "sg": _sg_text(sg),
        "k": sg.k_point_label or "",
        "irrep": sg.irrep_label or "",
        "opd": sg.opd_symbol or "",
        "cs": item.crystal_system or "",
        "max": "yes" if item.is_maximal else "",
        "_sg": sg,
        "_sort_sg": int(sg.space_group_number or 0),
        "_sort_cs": _CS_ORDER.get(item.crystal_system, 99),
        "_sort_max": 1 if item.is_maximal else 0,
    }


def _row_method2(sg) -> dict:
    return {
        "index": sg.index,
        "sg": _sg_text(sg),
        "k": sg.k_point_label or "",
        "irrep": sg.irrep_label or "",
        "opd": sg.opd_symbol or "",
        "s": sg.size,
        "i": sg.subgroup_index,
        "_sg": sg,
        "_sort_sg": int(sg.space_group_number or 0),
        "_sort_s": int(sg.size or 0),
        "_sort_i": int(sg.subgroup_index or 0),
    }


def _row_method3(item) -> dict:
    sg = item.subgroup
    return {
        "index": sg.index,
        "sg": _sg_text(sg),
        "k": sg.k_point_label or "",
        "irrep": sg.irrep_label or "",
        "pg": item.point_group or "",
        "_sg": sg,
        "_sort_sg": int(sg.space_group_number or 0),
    }


def _row_method4(mode: str, amp: float) -> dict:
    return {
        "mode": mode,
        "amp": f"{amp:.6f}",
        "_sort_amp": float(amp),
    }


def _matches_filters(row: dict, st: dict) -> bool:
    colmap = {c[0]: c for c in st["cols"]}
    for key, needle in st["filters"].items():
        text = str(needle or "").strip().lower()
        if not text or key not in colmap:
            continue
        if text not in str(colmap[key][2](row)).lower():
            return False
    return True


def _displayed_rows(st: dict, matching_only: bool = False) -> tuple[list[dict], list[dict]]:
    rows = list(st["rows"])
    key, direction = st["sort"]
    colmap = {c[0]: c for c in st["cols"]}
    if key and key in colmap:
        rows.sort(key=colmap[key][3], reverse=(direction == -1))
    matching = [r for r in rows if _matches_filters(r, st)]
    shown = matching if (matching_only or st["show_only"]) else rows
    return shown, matching


class IsoDistortConsoleApp:
    """终端交互应用：Search Page + Distortion Page。"""

    def __init__(self) -> None:
        self.project_root = Path(__file__).resolve().parent
        self.iso = IsoDistort()
        # 默认对齐网页/官网：strain + displacive（displacive 作用于全部物种）
        self.distortion_types: list[str] = ["strain", "displacive"]
        self.distortion_scope: dict[str, list[str]] = {
            "displacive": ["*"],
            "occupational": [],
            "strain": [],
            "magnetic": [],
            "rotational": [],
        }
        self.iso.set_distortion_scope(self.distortion_scope)
        self.iso.set_distortion_types(self.distortion_types)
        self.last_method1: list = []
        self.last_method2 = None
        self.last_method2_subgroups: list = []
        self.last_method3: list = []
        self.last_method4: list[dict] = []
        self.last_method4_meta: dict = {"rms": None, "max_abs": None}
        self.tbl: dict[int, dict] = {
            1: _empty_tbl(_method1_cols()),
            2: _empty_tbl(_method2_cols()),
            3: _empty_tbl(_method3_cols()),
            4: _empty_tbl(_method4_cols()),
        }

    def run(self) -> int:
        self._banner()
        self._load_parent_structure()
        while True:
            choice = self._main_menu()
            if choice == 0:
                print(t("ui.exit.done"))
                return 0
            if choice == 1:
                self._load_parent_structure()
            elif choice == 2:
                self._set_distortion_types()
            elif choice == 3:
                self._run_method_1()
            elif choice == 4:
                self._run_method_2()
            elif choice == 5:
                self._run_method_3()
            elif choice == 6:
                self._run_method_4()
            elif choice == 7:
                self._distortion_page_menu()
            elif choice == 8:
                self._show_state()

    # ----------------------------------------------------------------
    # 基础流程
    # ----------------------------------------------------------------

    def _banner(self) -> None:
        """打印使用前说明"""
        _line()
        print(t("ui.banner.title"))
        print(t("ui.banner.desc"))
        _line()

    def _load_parent_structure(self) -> None:
        cif_path = _choose_cif(self.project_root, t("ui.prompt.cif"))
        try:
            self.iso.load_structure(cif_path)
        except Exception as exc:  # noqa: BLE001 - CLI 边界：任何加载错误都转为提示
            print(t("load.failed", err=exc))
            return
        self.iso.set_distortion_scope(self.distortion_scope)
        self.iso.set_distortion_types(self.distortion_types)
        self.last_method1 = []
        self.last_method2 = None
        self.last_method2_subgroups = []
        self.last_method3 = []
        self.last_method4 = []
        self.last_method4_meta = {"rms": None, "max_abs": None}
        self.tbl = {
            1: _empty_tbl(_method1_cols()),
            2: _empty_tbl(_method2_cols()),
            3: _empty_tbl(_method3_cols()),
            4: _empty_tbl(_method4_cols()),
        }
        self._print_parent_header()

    def _print_parent_header(self) -> None:
        """Print the official Done / Space Group / lattice / prefs / Wyckoff block."""
        if self.iso.structure is None or not self.iso.symmetry_info:
            return
        sg = self.iso.symmetry_info["space_group_number"]
        sym = self.iso.symmetry_info["space_group_symbol"]
        sch = schoenflies_symbol(sg)
        lat = self.iso.structure.lattice
        print(t("web.done"))
        print(f"Space Group: {sg} {sym} {sch}")
        print(
            "Lattice parameters: "
            f"a= {lat.a:.5f}, b= {lat.b:.5f}, c= {lat.c:.5f}, "
            f"alpha= {lat.alpha:.5f}, beta= {lat.beta:.5f}, gamma= {lat.gamma:.5f}"
        )
        print("Default space-group preferences: "
              f"{self.iso.space_group_preferences()}")
        lines = format_wyckoff_sites(
            self.iso.structure, self.iso.symmetry_info["wyckoff_sites"]
        )
        print(",\n".join(lines))

    def _set_distortion_types(self) -> None:
        """选择 Distortion Types + 各类型的作用域物种（官网 all/none/Eu/Al）。"""
        self.distortion_types = _prompt_distortion_types(self.distortion_types)

        species = self.iso.species()
        enabled: list[str] = []
        for tp in self.distortion_types:
            if tp == "strain":
                enabled.append("strain")
                continue
            print(
                f"  {tp} scope (species in the loaded structure: "
                f"{', '.join(species) or 'none loaded'}; "
                "all = every species, none = disable this type, "
                "or comma-separated names such as Eu,Al)"
            )
            raw = _prompt(f"  {tp} scope", "all").strip().lower()
            if raw == "none":
                self.distortion_scope[tp] = []
                continue
            enabled.append(tp)
            if raw in ("all", "*", ""):
                self.distortion_scope[tp] = ["*"]
            else:
                self.distortion_scope[tp] = [s for s in raw.split(",") if s]

        self.distortion_types = enabled or ["strain"]
        self.iso.set_distortion_scope(self.distortion_scope)
        self.iso.set_distortion_types(self.distortion_types)
        scope_desc = ", ".join(
            f"{tp}={','.join(v) if v else 'none'}"
            for tp, v in self.distortion_scope.items()
            if tp in self.distortion_types
        )
        print(f"Enabled types: {', '.join(self.distortion_types)} | scope: {scope_desc}")

    def _main_menu(self) -> int:
        """Print the Search Page menu and read the next action."""
        _line()
        print(t("ui.search_page"))
        print(t("ui.menu.reload"))
        print(t("ui.menu.types"))
        print(t("ui.menu.method1"))
        print(t("ui.menu.method2"))
        print(t("ui.menu.method3"))
        print(t("ui.menu.method4"))
        print(t("ui.menu.distortion"))
        print(t("ui.menu.state"))
        print(t("ui.menu.exit"))
        return _prompt_int(t("ui.prompt.action"), 3)

    # ----------------------------------------------------------------
    # Method 1-4
    # ----------------------------------------------------------------

    def _run_method_1(self) -> None:
        _line()
        print("Method 1: Search over all special k points")

        print("Crystal system(s) (comma-separated; blank = no filter):")
        print("  " + ", ".join(sorted(CRYSTAL_SYSTEMS)))
        raw_cs = _prompt(t("m1.cs"), "").strip().lower()
        crystal_system = [x.strip() for x in raw_cs.split(",") if x.strip()] or None

        try:
            opts = self.iso.method1_options()
        except IsodistortError as exc:
            print(f"Failed to load Method 1 options: {exc}")
            return
        reachable = opts.get("space_groups") or []
        print("Space-group symmetry (reachable for current Types; blank = no filter):")
        for g in reachable:
            print(f"  {g['number']:3d} {g['symbol']} {g.get('schoenflies') or ''}")
        sg_raw = _prompt(t("m1.sg"), "").strip()
        subgroup_space_group = int(sg_raw) if sg_raw else None
        if subgroup_space_group is not None:
            allowed = {int(g["number"]) for g in reachable}
            if subgroup_space_group not in allowed:
                print("That space group is not in the reachable list; "
                      "search would return no rows. Continuing anyway.")

        maximal_only = _prompt_yes_no(t("m1.maximal"), False)
        lattice = self._prompt_lattice_selection(opts)
        lattice_matrix = None
        if lattice:
            lattice_matrix = self.iso.lattice_in_conventional_frame(
                lattice["matrix"], lattice["frame"]
            )

        print(t("m1.wait"))
        result = self.iso.search_method_1(
            distortion_types=self.distortion_types,
            crystal_system=crystal_system,
            subgroup_space_group=subgroup_space_group,
            lattice=lattice_matrix,
            maximal_subgroup_only=maximal_only,
        )
        self.last_method1 = result
        self.tbl[1] = _empty_tbl(_method1_cols())
        self.tbl[1]["rows"] = [_row_method1(item) for item in result]
        print(t("method1.result", n=len(result)))
        idx = self._review_result_table(1, allow_idx=True)
        if idx is not None:
            self._compute_modes(idx, "method1")

    def _prompt_lattice_selection(self, opts: dict | None = None) -> dict | None:
        """官网 Conventional/Primitive lattice 下拉的终端版（0 = 不选）。"""
        print("\nMethod 1 lattice filter (official Conventional/Primitive lattice; 0 = none)")
        if opts is None:
            try:
                opts = self.iso.method1_options()
            except IsodistortError as exc:
                print(f"Failed to load lattice options: {exc}")
                return None
        conv = opts.get("conventional_lattices") or []
        prim = opts.get("primitive_lattices") or []
        print("  Conventional lattice:")
        for i, lat in enumerate(conv, start=1):
            print(f"    C{i:2d}. {lat['label']}")
        print("  Primitive lattice:")
        for i, lat in enumerate(prim, start=1):
            print(f"    P{i:2d}. {lat['label']}")
        raw = _prompt("Choose (e.g. C3 or P2; 0 = none)", "0").strip().lower()
        if raw in ("", "0"):
            return None
        if len(raw) >= 2 and raw[0] in ("c", "p"):
            try:
                idx = int(raw[1:]) - 1
            except ValueError:
                idx = -1
            pool = conv if raw[0] == "c" else prim
            if 0 <= idx < len(pool):
                return {"matrix": pool[idx]["basis"], "frame": "conventional"}
        print("Invalid index; skipping lattice filter.")
        return None

    def _run_method_2(self) -> None:
        _line()
        print("Method 2: General method - search over specific k points")
        groups = self._prompt_kpoint_groups()
        if not groups:
            return

        print(t("m2.genDbHelp"))
        generate = _prompt_yes_no(t("lGenDb"), False)
        subs = self._enumerate_subgroups_for_groups(groups, generate_if_missing=generate)
        if not subs:
            return

        idx = self._review_result_table(2, allow_idx=True)
        if idx is None:
            return
        self._compute_modes(idx, "subgroups")

    def _prompt_kpoint_groups(self) -> list[dict]:
        """终端版 Method 2 k 点组输入（与网页 nsup + k vector 行一致）。"""
        kpoints = self.iso.list_k_points()
        print("\n--- Specify k point(s) ---")
        for i, kp in enumerate(kpoints, start=1):
            params = f" (parameters: {','.join(kp.parameters)})" if kp.parameters else ""
            kov = f" {kp.kovalev}" if getattr(kp, "kovalev", None) else ""
            print(f"  {i:2d}. {kp.label:<4s}{kov} {kp.coordinates}{params}")
        nsup = _prompt_int("Number of superposed IRs", 1)
        if nsup < 1:
            print(t("err.badNsup"))
            return []

        groups: list[dict] = []
        for i in range(1, nsup + 1):
            choice = _prompt_int(f"k-vector group {i}: k-point number", 1)
            if not (1 <= choice <= len(kpoints)):
                print(t("m2.range"))
                return []
            kp = kpoints[choice - 1]
            params: list[str] | None = None
            if kp.parameters:
                print(f"k point {kp.label} needs parameters (in order: {','.join(kp.parameters)})")
                vals = [_prompt(t("m2.param_value", p=p), "").strip() for p in kp.parameters]
                if any(not v for v in vals):
                    print(t("m2.param_empty"))
                    return []
                params = vals
            groups.append({"k": kp.label, "params": params})
        return groups

    def _enumerate_subgroups_for_groups(
        self, groups: list[dict], generate_if_missing: bool = False
    ) -> list:
        """枚举多组 k 点的全部 IR 子群（与网页 /api/subgroups 一致）。"""
        all_subs: list = []
        for i, grp in enumerate(groups, start=1):
            print(f"{t('m2.enumKp', grp['k'], i, len(groups))}")
            try:
                subs = self.iso.list_subgroups_at_kpoint(
                    grp["k"],
                    k_parameters=grp.get("params"),
                    generate_if_missing=generate_if_missing,
                )
            except IsodistortError as exc:
                text = str(exc)
                if ("在线生成" in text or "generate" in text.lower()) and _prompt_yes_no(
                    "The local subgroup database for this parametric k point is missing. "
                    "Generate it now (may take a long time)?",
                    False,
                ):
                    print(t("m2.genDbHelp"))
                    subs = self.iso.list_subgroups_at_kpoint(
                        grp["k"],
                        k_parameters=grp.get("params"),
                        generate_if_missing=True,
                    )
                else:
                    print(f"Enumeration failed: {exc}")
                    return []
            all_subs.extend(subs)

        for j, sg in enumerate(all_subs):
            sg.index = j
        self.iso.subgroups = all_subs
        self.last_method2_subgroups = list(all_subs)
        self.tbl[2] = _empty_tbl(_method2_cols())
        self.tbl[2]["rows"] = [_row_method2(sg) for sg in all_subs]
        if not all_subs:
            print(t("m2.noSubsAtKp"))
        return all_subs

    def _review_result_table(self, method: int, allow_idx: bool) -> int | None:
        """Filter/sort a Method result table. Optionally pick idx to compute modes."""
        st = self.tbl[method]
        if not st["rows"]:
            print(t("st.noSubs") if method != 4 else t("ui.export.table_none"))
            return None
        col_keys = [c[0] for c in st["cols"] if c[4]]
        col_help = ", ".join(col_keys)
        print(t("ui.tbl.cols", cols=col_help))
        print(t("ui.tbl.cmd") if allow_idx else t("ui.tbl.cmd4"))
        while True:
            shown, matching = _displayed_rows(st)
            print(
                f"\n{t('ui.tbl.shown', shown=len(shown), match=len(matching), total=len(st['rows']))}"
            )
            preview = shown[:80]
            for row in preview:
                self._print_table_row(method, row)
            if len(shown) > 80:
                print(f"  ... {len(shown) - 80} more")
            raw = _prompt("Next", "q").strip()
            if not raw or raw.lower() in {"q", "done"}:
                print(t("ui.tbl.done"))
                return None
            parts = raw.split()
            cmd = parts[0].lower()
            if allow_idx and cmd.isdigit():
                idx = int(cmd)
                if any(r.get("index") == idx for r in st["rows"]):
                    return idx
                print(t("m2.idx_range"))
                continue
            if cmd in {"f", "filter"}:
                spec = parts[1] if len(parts) > 1 else ""
                if "=" not in spec:
                    print("Use: f <col>=text   (blank text clears that column)")
                    continue
                col, _, val = spec.partition("=")
                col = col.strip().lower()
                if col not in col_keys:
                    print(t("ui.tbl.bad_col", cols=col_help))
                    continue
                st["filters"][col] = val
                _, matching = _displayed_rows(st)
                if not matching:
                    print(t("ui.tbl.no_match"))
                continue
            if cmd in {"c", "clear"}:
                st["filters"] = {}
                continue
            if cmd in {"s", "sort"}:
                if len(parts) < 3:
                    print("Use: s <col> a|d")
                    continue
                col = parts[1].lower()
                direction = parts[2].lower()
                sortable = {c[0] for c in st["cols"]}
                if col not in sortable:
                    print(t("ui.tbl.bad_col", cols=", ".join(sorted(sortable))))
                    continue
                if direction not in {"a", "d", "asc", "desc"}:
                    print("Use a (ascending) or d (descending).")
                    continue
                st["sort"] = (col, 1 if direction in {"a", "asc"} else -1)
                continue
            if cmd == "only":
                st["show_only"] = not st["show_only"]
                print("Show filtered rows only:" + (" on" if st["show_only"] else " off"))
                continue
            print(t("ui.tbl.bad_cmd"))
        return None

    def _print_table_row(self, method: int, row: dict) -> None:
        if method == 1:
            print(
                f"  idx={row['index']:3d} | SG {row['sg']:<16s} | k={row['k']:<4s} "
                f"IR={row['irrep']:<6s} OPD={row['opd']:<4s} "
                f"| crystal_system={row['cs']:<12s} | maximal={row['max'] or 'no'}"
            )
        elif method == 2:
            print(
                f"  idx={row['index']:3d} | SG {row['sg']:<16s} | k={row['k']:<4s} "
                f"IR={row['irrep']:<6s} OPD={row['opd']:<4s} "
                f"| s={row['s']:<3} i={row['i']:<3}"
            )
        elif method == 3:
            print(
                f"  idx={row['index']:3d} | SG {row['sg']:<16s} | k={row['k']:<4s} "
                f"IR={row['irrep']:<6s} | point_group={row['pg']}"
            )
        else:
            print(f"  {row['mode']:<16s} {row['amp']}")

    def _compute_modes(self, idx: int, source: str) -> None:
        if source == "method1":
            pool = [item.subgroup for item in self.last_method1]
            self.iso.subgroups = list(pool)
        elif source == "method3":
            pool = [item.subgroup for item in self.last_method3]
            self.iso.subgroups = list(pool)
        else:
            pool = list(self.last_method2_subgroups)
        target = next((s for s in pool if s.index == idx), None)
        if target is None:
            print(t("m2.idx_range"))
            return
        if getattr(target, "k_parameters", None):
            self.last_method2 = None
            self.iso.mode_displacements = {}
            self.iso.mode_occupancies = {}
            print(t("m2.paramKNote"))
            return
        result = self.iso.search_method_2(
            subgroup_idx=idx,
            distortion_type=self.distortion_types,
        )
        self.last_method2 = result
        print(f"Method 2: {len(result.modes) + len(self.iso.mode_occupancies)} mode(s)")
        for mode in result.modes:
            sites = sorted({b.wyckoff_letter for b in mode.bush_modes})
            print(
                f"  {mode.irrep_label:<8s} OPD={mode.opd_symbol:<6s} "
                f"dim={mode.dimension:<2d} sites={sites}"
            )
        for label, entry in self.iso.mode_occupancies.items():
            om = entry["mode"]
            flag = "" if entry["validated"] else " (approx)"
            print(f"  {label:<8s} occupational occupancy mode sites={om.wyckoff_letter}{flag}")

    def _run_method_3(self) -> None:
        _line()
        print("Method 3: Search over arbitrary k for point group and supercell")

        print("Point group (crystal class); blank = none:")
        for hm, sch in POINT_GROUP_SCHOENFLIES.items():
            system = POINT_GROUP_SYSTEM.get(hm, "")
            print(f"  {system}: {hm} {sch}")
        point_group = _prompt(t("m3.pg"), "").strip() or None

        print("Space-group symmetry: enter number 1-230 (same HM symbols as the web, "
              "e.g. 64 Cmca), blank = none, or type list to print all 230.")
        sg_raw = _prompt(t("m3.sg"), "").strip()
        if sg_raw.lower() == "list":
            for n in range(1, 231):
                print(f"  {n:3d} {hm_symbol(n)} {schoenflies_symbol(n)}")
            sg_raw = _prompt(t("m3.sg"), "").strip()
        space_group_type = int(sg_raw) if sg_raw else None

        lattice_type = _prompt("lattice_type (direct/reciprocal)", "direct").strip().lower()
        if lattice_type not in ("direct", "reciprocal"):
            lattice_type = "direct"

        centering = _prompt(
            "direct sublattice centering (d/P/A/B/C/I/F/R, blank = d)",
            "d",
        ).strip() or "d"
        if centering.upper() == "D":
            centering = "d"

        # Match web: always send a basis (default identity) so lattice filter is defined
        if _prompt_yes_no(t("m3.basis_q"), True):
            basis = _prompt_basis_matrix()
        else:
            basis = [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]]
            print("Using identity basis (same default as the web form).")

        try:
            result = self.iso.search_method_3(
                distortion_types=self.distortion_types,
                point_group=point_group,
                space_group_type=space_group_type,
                supercell_basis=basis,
                direct_sublattice_centering=centering,
                lattice_type=lattice_type,
            )
        except (IsodistortError, ValueError) as exc:
            print(f"Method 3 error: {exc}")
            return
        self.last_method3 = result
        self.tbl[3] = _empty_tbl(_method3_cols())
        self.tbl[3]["rows"] = [_row_method3(item) for item in result]
        print(t("method3.result", n=len(result)))
        idx = self._review_result_table(3, allow_idx=True)
        if idx is not None:
            self._compute_modes(idx, "method3")

    def _run_method_4(self) -> None:
        _line()
        print("Method 4: Mode decomposition of a distorted structure")
        daughter_cif = _choose_cif(self.project_root, "Choose a daughter CIF")

        # 与网页一致：使用官网默认匹配参数（nearest-site / 阈值 0.25）
        result = self.iso.search_method_4(
            distorted_cif_path=daughter_cif,
            atom_matching_method="nearest-site",
            robust_distance_threshold=0.25,
            provided_origin_shift=None,
        )

        ranked = sorted(result.amplitudes.items(), key=lambda kv: abs(kv[1]), reverse=True)
        self.last_method4 = [_row_method4(label, float(amp)) for label, amp in ranked]
        self.last_method4_meta = {
            "rms": result.rms_residual,
            "max_abs": result.max_abs_residual,
        }
        self.tbl[4] = _empty_tbl(_method4_cols())
        self.tbl[4]["rows"] = list(self.last_method4)
        print(t("m4.result_title"))
        print(
            f"RMS residual: {result.rms_residual:.8e}  "
            f"Max residual: {result.max_abs_residual:.8e}"
        )
        self._review_result_table(4, allow_idx=False)

    # ----------------------------------------------------------------
    # Distortion: export subgroup files
    # ----------------------------------------------------------------

    def _distortion_page_menu(self) -> None:
        while True:
            _line()
            print(t("ui.distortion_page"))
            print(t("ui.dist.single"))
            print(t("ui.dist.mixed"))
            print(t("ui.dist.export"))
            print(t("ui.dist.table"))
            print(t("ui.dist.back"))
            choice = _prompt_int(t("ui.dist.choice"), 1)
            if choice == 0:
                return
            if choice in (1, 2, 3):
                self._export_subgroups_flow(choice)
            elif choice == 4:
                self._export_filtered_table()

    def _subgroups_for_method(self, method: int) -> list:
        if method == 1:
            return [item.subgroup for item in self.last_method1]
        if method == 3:
            return [item.subgroup for item in self.last_method3]
        return list(self.last_method2_subgroups)

    def _export_subgroups_flow(self, method: int) -> None:
        subs = self._subgroups_for_method(method)
        if not subs:
            print(t("ui.export.none"))
            return
        st = self.tbl.get(method)
        if st and any(str(v).strip() for v in st.get("filters", {}).values()):
            _rows, matching = _displayed_rows(st, matching_only=True)
            want = {row["index"] for row in matching}
            subs = [sg for sg in subs if sg.index in want]
        if not subs:
            print(t("ui.export.none"))
            return
        formats_raw = _prompt(t("ui.export.formats"), "cif,isoviz,modes,topas")
        formats = [x.strip().lower() for x in formats_raw.split(",") if x.strip()]
        as_zip = _prompt_yes_no("Write a ZIP file (same as the web Download all)?", True)
        default_dest = str(
            self.iso.cfg.output_dir
            / (f"isodistort_method{method}.zip" if as_zip else f"isodistort_method{method}")
        )
        dest = _prompt(t("ui.export.dest"), default_dest)
        need_modes = any(fmt != "cif" for fmt in formats)
        saved = list(self.iso.subgroups)
        try:
            self.iso.subgroups = list(subs)
            if as_zip:
                body = self.iso.export_subgroups_zip(
                    formats=formats,
                    subgroups=subs,
                    compute_missing_modes=need_modes,
                    wrapping=f"isodistort_method{method}",
                )
                out = Path(dest)
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(body)
                print(t("ui.export.done", n=len(subs), dest=out))
                return
            paths = self.iso.export_subgroups(
                dest,
                formats=formats,
                subgroups=subs,
                compute_missing_modes=need_modes,
            )
        finally:
            self.iso.subgroups = saved
        print(t("ui.export.done", n=len(paths), dest=dest))
        for path in paths[:20]:
            print(f"  {path}")
        if len(paths) > 20:
            print(f"  ... {len(paths) - 20} more")

    def _export_filtered_table(self) -> None:
        method = _prompt_int(t("ui.dist.table_method"), 2)
        if method not in (1, 2, 3, 4):
            print(t("m2.range"))
            return
        st = self.tbl.get(method)
        if not st or not st["rows"]:
            print(t("ui.export.table_none"))
            return
        fmt = _prompt(t("ui.dist.table_fmt"), "txt").strip().lower()
        if fmt not in {"txt", "csv"}:
            print("Use txt or csv.")
            return
        rows, _matching = _displayed_rows(st, matching_only=True)
        export_headers = {
            1: ["index", "space_group", "k_point", "irrep", "OPD", "crystal_system", "maximal"],
            2: ["index", "space_group", "k_point", "irrep", "OPD", "size", "subgroup_index"],
            3: ["index", "space_group", "k_point", "irrep", "point_group"],
            4: ["mode", "amplitude"],
        }
        cols = st["cols"]
        headers = export_headers[method]
        data = [[c[2](row) for c in cols] for row in rows]
        default_dest = str(
            self.iso.cfg.output_dir / f"isodistort_method{method}_filtered.{fmt}"
        )
        dest = Path(_prompt(t("ui.export.dest"), default_dest))
        dest.parent.mkdir(parents=True, exist_ok=True)
        if fmt == "csv":
            with dest.open("w", encoding="utf-8-sig", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(headers)
                writer.writerows(data)
        else:
            lines = ["\t".join(headers)]
            lines.extend("\t".join(map(str, row)) for row in data)
            dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(t("ui.export.table_done", path=dest, n=len(rows)))

    # ----------------------------------------------------------------
    # State
    # ----------------------------------------------------------------
    def _show_state(self) -> None:
        _line()
        print(t("ui.msg.state"))
        if self.iso.structure is None:
            print(f"  Parent structure: {t('state.not_loaded')}")
        else:
            self._print_parent_header()

        print(f"  Distortion types: {', '.join(self.distortion_types)}")
        scope_desc = ", ".join(
            f"{tp}={','.join(v) if v else 'none'}"
            for tp, v in self.distortion_scope.items()
            if tp in self.distortion_types
        )
        print(f"  Scope: {scope_desc}")
        print(f"  Last Method1 count: {len(self.last_method1)}")
        print(f"  Last Method2 subgroup count: {len(self.last_method2_subgroups)}")
        print(f"  Last Method3 count: {len(self.last_method3)}")
        print(f"  Last Method4 mode count: {len(self.last_method4)}")

        mode_count = len(self.iso.mode_displacements) + len(self.iso.mode_occupancies)
        print(f"  Available mapped modes: {mode_count}")


def main() -> int:
    # Windows 控制台默认代码页可能不是 UTF-8，统一重配置避免中文乱码
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        return IsoDistortConsoleApp().run()
    except KeyboardInterrupt:
        print(f"\n{t('ui.exit.interrupt')}")
        return 130
    except IsodistortError as exc:
        print(f"\nFailed: {exc}")
        return 1
    except Exception as exc:  # noqa: BLE001 - CLI 入口兜底，防止未捕获异常直接崩溃
        print(f"\nFailed: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
