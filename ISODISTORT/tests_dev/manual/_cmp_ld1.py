from pathlib import Path

root = Path(r"C:\Users\devou\OneDrive\Desktop\CRIS\output_compare")
offs = [p for p in root.rglob("subgroup.cif") if p.parent.name == "LD1 C1"]
locs = [p for p in root.rglob("*CIF.cif") if p.parent.name == "LD1 C1"]
print("off", offs[0])
print("loc", locs[0])
off = offs[0].read_text(encoding="utf-8", errors="replace")
loc = locs[0].read_text(encoding="utf-8", errors="replace")
print("=== OFFICIAL CIF head ===")
print("\n".join(off.splitlines()[:50]))
print("=== LOCAL CIF head ===")
print("\n".join(loc.splitlines()[:50]))
print("=== OFFICIAL first atom lines ===")
for line in off.splitlines():
    if line[:4] in ("Eu1_", "Al1_", "Al2_") or line.startswith("Eu1 ") or "Eu1_1" in line[:10]:
        print(line)
        if line.startswith("Al2_1") or line.startswith("Al2_1 "):
            break
print("=== LOCAL first atom lines ===")
n = 0
for line in loc.splitlines():
    if any(line.startswith(p) for p in ("Eu", "Al")) and " " in line and not line.startswith("#"):
        print(line)
        n += 1
        if n >= 6:
            break
print("=== OFFICIAL isoviz head ===")
iso = next(offs[0].parent.glob("data.isoviz")).read_text(encoding="utf-8", errors="replace")
print("\n".join(iso.splitlines()[:40]))
print("=== LOCAL isoviz head ===")
liso = next(locs[0].parent.glob("*.isoviz")).read_text(encoding="utf-8", errors="replace")
print("\n".join(liso.splitlines()[:40]))
print("=== OFFICIAL topas head ===")
print("\n".join(next(offs[0].parent.glob("topas.str")).read_text(encoding="utf-8", errors="replace").splitlines()[:35]))
print("=== LOCAL topas head ===")
print("\n".join(next(locs[0].parent.glob("*TOPAS*")).read_text(encoding="utf-8", errors="replace").splitlines()[:35]))
