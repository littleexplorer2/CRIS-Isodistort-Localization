from isocore.data.irreps_cdml import _ml_kov_tables, lookup_irrep_kovalev
from isocore.io.isodistort_cif import _irrep_names_for_kovalev
import re

ml, kov = _ml_kov_tables()
k_kov = "k15"
ir = "M3-"
companions = [
    name for name, _tag in _irrep_names_for_kovalev(k_kov)
    if name.startswith("M")
][:24]
print("companions used", companions)

k_prefix = "k15"
candidates = []
for i, name in enumerate(ml):
    if name != ir:
        continue
    tag = kov[i]
    if k_prefix and not tag.lower().startswith(k_prefix + "t"):
        continue
    candidates.append(i)
print("candidates", [(i, kov[i]) for i in candidates])

companions_set = set(companions)
companions_set.discard(ir)
for i in candidates:
    lo, hi = max(0, i - 24), min(len(ml), i + 25)
    nearby = set(ml[lo:hi])
    score = len(companions_set & nearby)
    tag = kov[i]
    if tag.lower().startswith(k_prefix + "t"):
        score += 5
    print(i, tag, "score", score, "overlap", sorted(companions_set & nearby)[:12])

print("lookup", lookup_irrep_kovalev(ir, k_kovalev=k_kov, companion_irreps=companions))
print("lookup no companions", lookup_irrep_kovalev(ir, k_kovalev=k_kov))
