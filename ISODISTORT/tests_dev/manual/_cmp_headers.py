"""Diff official vs local Method1/2 exports after NdNiO2 refresh."""
from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3] / "output_compare"


def _grab(lines: list[str], pref: str) -> str | None:
    for ln in lines:
        if ln.startswith(pref):
            return ln
    return None


def _pair_dirs(off_root: Path, loc_root: Path) -> list[tuple[Path, Path]]:
    locs = {p.parent.name: p.parent for p in loc_root.rglob("subgroup.cif")}
    # also match modes-only folders
    for p in loc_root.iterdir() if loc_root.is_dir() else []:
        if p.is_dir():
            locs.setdefault(p.name, p)
    pairs = []
    for op in sorted(off_root.rglob("subgroup.cif")):
        name = op.parent.name
        lp = locs.get(name)
        if lp is None:
            for kn, kp in locs.items():
                if name[:28] == kn[:28] or kn[:28] == name[:28]:
                    lp = kp
                    break
        if lp is not None:
            pairs.append((op.parent, lp))
    return pairs


def _count_modes(text: str) -> tuple[int, int]:
    # rough: ndisp / nstrain style comments or mode amplitude lines
    nd = len(re.findall(r"\bndisp\b|\bdisplacive\b", text, re.I))
    ns = len(re.findall(r"\bnstrain\b|\bstrain\b", text, re.I))
    return nd, ns


def analyze(cif: str, method: str) -> None:
    off_root = ROOT / cif / "官网" / method
    loc_root = ROOT / cif / "现有网页版交互" / method
    print(f"\n==== {cif} {method} ====")
    if not off_root.is_dir() or not loc_root.is_dir():
        print(" missing dirs", off_root.exists(), loc_root.exists())
        return
    off_n = len(list(off_root.iterdir()))
    loc_n = len(list(loc_root.iterdir()))
    pairs = _pair_dirs(off_root, loc_root)
    print(f" folders off={off_n} loc={loc_n} paired={len(pairs)}")

    stats = Counter()
    samples: dict[str, list[str]] = defaultdict(list)

    for od, ld in pairs:
        for fname, key in [
            ("subgroup.cif", "cif"),
            ("data.isoviz", "isoviz"),
            ("topas.str", "topas"),
        ]:
            of, lf = od / fname, ld / fname
            if of.exists() and not lf.exists():
                stats[f"missing_local_{key}"] += 1
            elif lf.exists() and not of.exists():
                stats[f"missing_off_{key}"] += 1
            elif of.exists() and lf.exists():
                stats[f"both_{key}"] += 1

        # modes: official HTML or Complete modes; local Complete modes.txt
        off_modes = list(od.glob("*modes*")) + list(od.glob("*.html"))
        loc_modes = list(ld.glob("*modes*")) + list(ld.glob("*.txt"))
        if off_modes and loc_modes:
            stats["both_modes"] += 1

        oc = od / "subgroup.cif"
        lc = ld / "subgroup.cif"
        if not (oc.exists() and lc.exists()):
            continue
        ot = oc.read_text(encoding="utf-8", errors="replace").splitlines()
        lt = lc.read_text(encoding="utf-8", errors="replace").splitlines()

        checks = [
            ("k", "# k point:"),
            ("ir", "# IR:"),
            ("hm", '_symmetry_space_group_name_H-M'),
            ("it", "_symmetry_Int_Tables_number"),
            ("a", "_cell_length_a"),
            ("b", "_cell_length_b"),
            ("c", "_cell_length_c"),
            ("vol", "_cell_volume"),
            ("opd", "# Order parameter"),  # next lines vary
        ]
        for tag, pref in checks:
            oh, lh = _grab(ot, pref), _grab(lt, pref)
            if oh is None and lh is None:
                continue
            if oh != lh:
                stats[f"diff_{tag}"] += 1
                if len(samples[tag]) < 3:
                    samples[tag].append(
                        f"{od.name[:55]}\n  off={oh}\n  loc={lh}"
                    )

        # origin / transform
        oh = _grab(ot, "_space_group.transform_Pp_abc")
        lh = _grab(lt, "_space_group.transform_Pp_abc")
        if oh != lh:
            stats["diff_p2c"] += 1
            if len(samples["p2c"]) < 3:
                samples["p2c"].append(
                    f"{od.name[:55]}\n  off={oh}\n  loc={lh}"
                )

        # atom site count
        def n_atoms(lines):
            in_loop = False
            n = 0
            for ln in lines:
                if ln.startswith("_atom_site_"):
                    in_loop = True
                    continue
                if in_loop:
                    if not ln.strip() or ln.startswith("_") or ln.startswith("loop_"):
                        if n:
                            break
                        continue
                    if ln.startswith("#"):
                        continue
                    parts = ln.split()
                    if len(parts) >= 4:
                        n += 1
            return n

        no, nl = n_atoms(ot), n_atoms(lt)
        if no != nl:
            stats["diff_natom"] += 1
            if len(samples["natom"]) < 3:
                samples["natom"].append(
                    f"{od.name[:55]} off={no} loc={nl}"
                )

        # header parent wyckoff lines
        ow = [ln for ln in ot if re.match(r"^# [A-Za-z].*\(", ln) and "Space" not in ln][:6]
        lw = [ln for ln in lt if re.match(r"^# [A-Za-z].*\(", ln) and "Space" not in ln][:6]
        # better: lines after Lattice parameters until Include
        def wyck(lines):
            out = []
            grab = False
            for ln in lines:
                if ln.startswith("# Lattice parameters"):
                    grab = True
                    continue
                if grab:
                    if ln.startswith("# Include") or ln.startswith("# k point"):
                        break
                    if ln.startswith("# ") and not ln.startswith("# Default"):
                        out.append(ln)
            return out

        ow, lw = wyck(ot), wyck(lt)
        if ow != lw:
            stats["diff_wyck"] += 1
            if len(samples["wyck"]) < 2:
                samples["wyck"].append(
                    f"{od.name[:40]}\n  off={ow}\n  loc={lw}"
                )

    print("stats:", dict(stats))
    for k, items in samples.items():
        print(f"-- samples {k} --")
        for s in items:
            print(s)


if __name__ == "__main__":
    for cif in ("NdNiO2 own.cif", "EuAl4 Parent.cif"):
        for method in ("Method1", "Method2"):
            analyze(cif, method)
