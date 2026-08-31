from isocore.data.irreps_cdml import _ml_kov_tables, lookup_irrep_kovalev

ml, kov, sgs = _ml_kov_tables()
for label in ["Y1", "Y2", "Y3", "Y4"]:
    hits = [(i, kov[i], sgs[i]) for i, n in enumerate(ml) if n == label]
    k6 = [h for h in hits if h[1].lower().startswith("k6t")]
    print(label, "total", len(hits), "k6", k6[:5], "first", hits[:3])
    print("  lookup", lookup_irrep_kovalev(label, k_kovalev="k6", parent_sg=123, k_point_label="Y"))
