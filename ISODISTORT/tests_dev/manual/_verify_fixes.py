"""Verify CIF header tags vs official after fixes."""
from __future__ import annotations

import re
from pathlib import Path

from isocore.backend import SubgroupInfo
from isocore.data.irreps_cdml import lookup_irrep_kovalev
from isocore.data.kpoints_official import format_k_point_display
from isocore.io.isodistort_cif import (
    _HM_FULL_1992,
    _irrep_kovalev_tag,
    _k_point_comment,
    _space_hm_symbol,
)

ROOT = Path(__file__).resolve().parents[3] / "output_compare"


def check_lookup_vs_official(cif_name: str, method: str) -> None:
    off = ROOT / cif_name / "官网" / method
    pairs: dict[tuple[str, str, int], str] = {}
    for p in off.rglob("subgroup.cif"):
        t = p.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"# IR:\s*(\S+),\s*(\S+)", t)
        kp = re.search(r"# k point:\s*([A-Za-z]+),\s*(k\d+)", t)
        sg = re.search(r"# Space Group:\s*(\d+)", t)
        if m and kp and sg:
            pairs[(m.group(1), kp.group(1), int(sg.group(1)))] = m.group(2)
    bad = 0
    for (ir, klab, parent), want in sorted(pairs.items()):
        # resolve kov from table via format
        from isocore.data.kpoints_official import KPOINT_OFFICIAL

        entry = KPOINT_OFFICIAL.get(parent, {}).get(klab)
        kk = entry[0] if entry else None
        got = lookup_irrep_kovalev(
            ir, k_kovalev=kk, parent_sg=parent, k_point_label=klab
        )
        if got != want:
            bad += 1
            if bad <= 12:
                print(f"  BAD {ir} {klab} parent={parent} got={got} want={want}")
    print(f"{cif_name} {method}: unique={len(pairs)} mismatches={bad}")


def check_hm() -> None:
    cases = {
        "Pmm2": "P m m 2",
        "P1m1": "P 1 m 1",
        "Amm2": "A m m 2",
        "Imm2": "I m m 2",
        "I4/mmm": "I 4/m m m",
        "I2/m2/m2/m": "I 2/m 2/m 2/m",
        "P4_2nmc": "P 42 n m c",
    }
    for raw, want in cases.items():
        got = _space_hm_symbol(raw)
        print(f"HM {raw!r} -> {got!r} {'OK' if got == want else 'BAD want '+want}")
    for num in (6, 25, 38, 115):
        print(f"table {num}", _HM_FULL_1992.get(num))


def sim_headers() -> None:
    cases = [
        SubgroupInfo(
            index=0, space_group_number=129, space_group_symbol="P4/nmm",
            irrep_label="M3-", k_point_label="M", parent_sg=139,
            k_coordinates=["1", "1", "1"],
        ),
        SubgroupInfo(
            index=0, space_group_number=123, space_group_symbol="P4/mmm",
            irrep_label="LD1", k_point_label="LD", parent_sg=139,
            k_coordinates=["0", "0", "g"], k_parameters=["1/6"],
        ),
        SubgroupInfo(
            index=0, space_group_number=25, space_group_symbol="Pmm2",
            irrep_label="A1+", k_point_label="A", parent_sg=123,
            k_coordinates=["1/2", "1/2", "1/2"],
        ),
        SubgroupInfo(
            index=0, space_group_number=6, space_group_symbol="Pm",
            irrep_label="Y1", k_point_label="Y", parent_sg=123,
            k_coordinates=["a", "1/2", "0"], k_parameters=["1/3"],
        ),
        SubgroupInfo(
            index=0, space_group_number=139, space_group_symbol="I4/mmm",
            irrep_label="GM2+", k_point_label="GM", parent_sg=139,
            k_coordinates=["0", "0", "0"],
        ),
    ]
    for sg in cases:
        print(
            sg.irrep_label,
            "|",
            _k_point_comment(sg),
            "|",
            _irrep_kovalev_tag(sg),
        )


if __name__ == "__main__":
    check_hm()
    print("--- sim ---")
    sim_headers()
    print("--- vs official ---")
    for cif in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        for method in ("Method1", "Method2"):
            check_lookup_vs_official(cif, method)
