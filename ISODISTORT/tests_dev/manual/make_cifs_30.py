"""生成 30 个覆盖全晶系的测试母相 CIF（pymatgen from_spacegroup + 对称性自检）。

供 30-CIF 科学验证使用：每个 CIF 生成后用 SpacegroupAnalyzer 自检，
只保留识别空间群与目标一致的候选。
"""
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

OUT_DIR = Path(__file__).resolve().parents[1] / "cifs_30"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# (sg, lattice, species, coords) —— 常见原型结构（低对称空间群用一般位置坐标）
CANDIDATES = [
    (1, Lattice.from_parameters(4.0, 5.0, 6.0, 80, 70, 85), ["Fe", "O"], [[0.11, 0.22, 0.33], [0.37, 0.61, 0.82]]),
    (2, Lattice.from_parameters(4.5, 5.0, 6.0, 82, 75, 88), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
    (12, Lattice.monoclinic(5.0, 6.0, 7.0, 99.0), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.0]]),
    (14, Lattice.monoclinic(5.0, 6.0, 7.0, 101.0), ["Fe"], [[0.24, 0.13, 0.19]]),
    (33, Lattice.orthorhombic(5.0, 6.0, 7.0), ["Fe", "O"], [[0.22, 0.11, 0.31], [0.31, 0.42, 0.17]]),
    (62, Lattice.orthorhombic(5.0, 6.0, 7.0), ["Fe", "O"], [[0.23, 0.11, 0.19], [0.73, 0.61, 0.69]]),
    (63, Lattice.orthorhombic(5.0, 6.0, 7.0), ["Fe", "O"], [[0.0, 0.31, 0.25], [0.0, 0.81, 0.75]]),
    (65, Lattice.orthorhombic(5.0, 6.0, 7.0), ["Fe"], [[0.0, 0.0, 0.0]]),
    (74, Lattice.orthorhombic(5.0, 6.0, 7.0), ["Fe", "O"], [[0.13, 0.22, 0.28], [0.17, 0.61, 0.74]]),
    (99, Lattice.tetragonal(4.0, 8.0), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.25]]),
    (123, Lattice.tetragonal(4.0, 8.0), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
    (129, Lattice.tetragonal(4.0, 8.0), ["Fe", "O"], [[0.21, 0.32, 0.35], [0.71, 0.82, 0.15]]),
    (139, Lattice.tetragonal(4.4, 11.2), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.38]]),
    (140, Lattice.tetragonal(4.4, 11.2), ["Fe", "O"], [[0.21, 0.32, 0.12], [0.71, 0.18, 0.62]]),
    (141, Lattice.tetragonal(4.4, 11.2), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.25]]),
    (147, Lattice.hexagonal(5.0, 8.0), ["Fe"], [[0.2, 0.3, 0.15]]),
    (148, Lattice.hexagonal(5.0, 8.0), ["Fe", "O"], [[0.15, 0.25, 0.2], [0.65, 0.75, 0.7]]),
    (160, Lattice.hexagonal(5.0, 8.0), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.25]]),
    (166, Lattice.hexagonal(5.0, 8.0), ["Fe", "O"], [[0.0, 0.0, 0.2], [0.0, 0.0, 0.35]]),
    (167, Lattice.hexagonal(5.0, 8.0), ["Fe", "O"], [[0.15, 0.2, 0.22], [0.31, 0.47, 0.58]]),
    (186, Lattice.hexagonal(3.2, 5.2), ["Fe", "O"], [[1 / 3, 2 / 3, 0.0], [1 / 3, 2 / 3, 0.375]]),
    (191, Lattice.hexagonal(3.2, 5.2), ["Fe", "O"], [[0.0, 0.0, 0.0], [1 / 3, 2 / 3, 0.25]]),
    (194, Lattice.hexagonal(3.2, 5.2), ["Fe", "O"], [[1 / 3, 2 / 3, 0.25], [1 / 3, 2 / 3, 0.75]]),
    (200, Lattice.cubic(4.0), ["Fe", "O"], [[0.15, 0.2, 0.1], [0.65, 0.7, 0.6]]),
    (205, Lattice.cubic(4.0), ["Fe", "O"], [[0.15, 0.2, 0.3], [0.65, 0.7, 0.8]]),
    (216, Lattice.cubic(5.6), ["Fe", "O"], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
    (221, Lattice.cubic(3.2), ["Fe"], [[0.0, 0.0, 0.0]]),
    (225, Lattice.cubic(3.6), ["Fe"], [[0.0, 0.0, 0.0]]),
    (227, Lattice.cubic(3.57), ["Fe"], [[0.125, 0.125, 0.125]]),
    (229, Lattice.cubic(3.2), ["Fe"], [[0.0, 0.0, 0.0]]),
]

written = []
for sg, lattice, species, coords in CANDIDATES:
    try:
        s = Structure.from_spacegroup(sg, lattice, species, coords)
    except Exception as exc:  # noqa: BLE001
        print(f"SKIP sg={sg}: from_spacegroup failed: {exc}")
        continue
    # 重叠原子自检（坐标被空间群操作映射到同一位置时 CIF 无法回读）
    fc = np.round(s.frac_coords, 6)
    overlap = False
    for i in range(len(fc)):
        d = fc - fc[i]
        d -= np.round(d)
        if np.any((np.linalg.norm(d, axis=1) < 1e-4)
                  & (np.arange(len(fc)) != i)):
            overlap = True
            break
    if overlap:
        print(f"SKIP sg={sg}: 原子重叠（坐标被群操作映射到同一位置）")
        continue
    try:
        got = SpacegroupAnalyzer(s).get_space_group_number()
    except Exception as exc:  # noqa: BLE001
        print(f"SKIP sg={sg}: analysis failed: {exc}")
        continue
    if got != sg:
        print(f"SKIP sg={sg}: identified as {got} (supergroup/placement issue)")
        continue
    name = f"sg{sg:03d}"
    path = OUT_DIR / f"{name}.cif"
    CifWriter(s).write_file(str(path))
    # 回读校验：确保生成的 CIF 能被 CifParser 重新解析（科研验证数据可用性）
    try:
        from pymatgen.io.cif import CifParser
        reread = CifParser(str(path)).parse_structures()[0]
        assert len(reread) == len(s)
    except Exception as exc:  # noqa: BLE001
        path.unlink(missing_ok=True)
        print(f"SKIP sg={sg}: CIF 回读失败: {exc}")
        continue
    written.append((sg, len(s), str(path)))
    print(f"OK  sg={sg:3d} atoms={len(s):2d} -> {path.name}")

print(f"\nTotal valid: {len(written)}/30")
for sg, n, p in written:
    print(f"  {p} (sg {sg}, {n} atoms)")
