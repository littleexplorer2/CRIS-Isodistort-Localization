from isocore.data.irreps_cdml import _ml_kov_tables

ml, kov = _ml_kov_tables()

# Find LD-like labels
ld = sorted({n for n in ml if n.upper().startswith("LD") or "LD" in n.upper()})
print("LD-like", ld[:50], "count", len(ld))

# Find labels near k10t for tetragonal / I4mmm style
for i, tag in enumerate(kov):
    if tag.strip().lower() == "k10t1":
        print("--- k10t1 at", i, "ml=", ml[i])
        for j in range(max(0, i - 6), min(len(ml), i + 10)):
            print(f"  {j}: {ml[j]:10s} {kov[j]}")

# M3- k15 neighborhoods
for i, (n, t) in enumerate(zip(ml, kov)):
    if n == "M3-" and t.strip() in {"k15t3", "k15t4"}:
        print("---", n, t, "at", i)
        for j in range(max(0, i - 10), min(len(ml), i + 12)):
            print(f"  {j}: {ml[j]:10s} {kov[j]}")
