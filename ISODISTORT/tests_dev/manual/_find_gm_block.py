from isocore.data.irreps_cdml import _ml_kov_tables

ml, kov = _ml_kov_tables()
for i, (n, t) in enumerate(zip(ml, kov)):
    if n == "GM2+" and t.strip().startswith("k14"):
        print("---", i, n, t)
        for j in range(max(0, i - 8), min(len(ml), i + 12)):
            print(f"  {j}: {ml[j]:10s} {kov[j]}")

print("\n=== want block for I4/mmm GM: GM1+..GM5+/- ===")
# Find block containing GM1+ k14t1 and GM2+ k14t5
for i, (n, t) in enumerate(zip(ml, kov)):
    if n == "GM1+" and t.strip() == "k14t1":
        window = list(zip(ml[i:i+12], kov[i:i+12]))
        labels = [x[0] for x in window]
        if "GM2+" in labels and "GM5-" in labels:
            print("candidate block at", i)
            for a, b in window:
                print(f"  {a:10s} {b}")
