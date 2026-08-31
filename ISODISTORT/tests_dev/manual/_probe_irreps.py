from isocore.data.irreps_cdml import _ml_kov_tables, lookup_irrep_kovalev

ml, kov = _ml_kov_tables()
print("n", len(ml))
for label in ["LD1", "LD1+", "LD1-", "M3-", "M4-", "M5+", "GM1+", "A1+", "Y1"]:
    hits = [(ml[i], kov[i]) for i in range(len(ml)) if ml[i] == label]
    print(label, "hits", len(hits), hits[:8])

# neighbors around first M3-
idxs = [i for i, n in enumerate(ml) if n == "M3-"]
for i in idxs[:3]:
    print("--- idx", i, "tag", kov[i])
    for j in range(max(0, i - 8), min(len(ml), i + 9)):
        print(f"  {j}: {ml[j]:8s} {kov[j]}")
