from isocore.data.irreps_cdml import lookup_irrep_kovalev
from isocore.data.kpoints_official import format_k_point_display, KPOINT_OFFICIAL
from isocore.io.isodistort_cif import _space_hm_symbol, _HM_FULL_1992

print("LD1 k10 ->", lookup_irrep_kovalev("LD1", k_kovalev="k10"))
print("A1+ k20 ->", lookup_irrep_kovalev("A1+", k_kovalev="k20"))
print("Y1 k6 ->", lookup_irrep_kovalev("Y1", k_kovalev="k6"))
print("M3- k15 ->", lookup_irrep_kovalev("M3-", k_kovalev="k15"))
print("M3- k15 companions GM...", lookup_irrep_kovalev(
    "M3-", k_kovalev="k15", companion_irreps=["M1+", "M2+", "M3+", "M4+", "M5+", "M1-", "M2-", "M3-", "M4-", "M5-"]
))
print("123 in table", 123 in KPOINT_OFFICIAL)
print("format A", format_k_point_display(123, "A", [], ["1/2", "1/2", "1/2"]))
print("format Y", format_k_point_display(123, "Y", ["1/3"], ["a", "1/2", "0"]))

for raw in ["Pmm2", "P1m1", "Amm2", "Imm2", "Pm", "C2/m", "P4_2nmc", "I4/mmm"]:
    print(raw, "->", repr(_space_hm_symbol(raw)))

# missing from HM table among common
missing = [n for n in range(1, 143) if n not in _HM_FULL_1992]
print("missing HM count", len(missing), missing[:40])
