"""Read-only scientific audits of frozen phase-2 probes, not production logic."""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "evidence/phase2/iso-contract"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(case):
    return (BASE / case / "stdout.txt").read_text()


def tuples(text):
    return re.findall(r"\(([^()]+)\)", text)


def vector(text):
    return tuple(F(x.strip()) for x in text.split(","))


def det(m):
    return sum(m[0][j] * (
        m[1][(j + 1) % 3] * m[2][(j + 2) % 3]
        - m[1][(j + 2) % 3] * m[2][(j + 1) % 3]
    ) for j in range(3))


def mat(columns):
    return tuple(tuple(columns[j][i] for j in range(3)) for i in range(3))


def product(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def parent_group(case, number):
    raw = read(case)
    row = next(line for line in raw.splitlines() if line.startswith(str(number) + " "))
    primitive = [vector(t) for t in tuples(row)[:3]]
    reps = {tuple(F(0) for _ in range(3))}
    pending = list(reps)
    while pending:
        current = pending.pop()
        for p in primitive:
            new = tuple((current[i] + p[i]) % 1 for i in range(3))
            if new not in reps:
                reps.add(new)
                pending.append(new)
            if len(reps) > 16:
                raise AssertionError("Unexpected translation quotient size.")
    ops = [re.sub(r"\s+", "", t) for t in tuples(raw) if re.search("[xyz]", t)]
    primitive_det = abs(det(mat(primitive)))
    assert len(reps) == 1 / primitive_det
    return {
        "case": case, "primitive_basis_columns": primitive,
        "primitive_basis_determinant_absolute": primitive_det,
        "translation_representatives_mod_conventional_cell": sorted(reps),
        "point_coset_count": len(ops),
        "conventional_cell_coset_count": len(ops) * len(reps),
        "validated": "Exact additive closure of translations and index=1/|det P|.",
        "not_validated": "No independent CIF symmetry loop was supplied for comparison.",
    }


def main():
    # Retain the small generated subgroup tables/logs that earlier probes saved
    # only in runs/. This is an evidence copy, never a compute-cache import.
    for case in sorted(BASE.iterdir()):
        if not case.is_dir() or not (case / "result.json").exists():
            continue
        result = json.loads((case / "result.json").read_text())
        run_dir = Path(result["run_directory"])
        if not run_dir.is_relative_to(ROOT / "runs"):
            raise ValueError("Result run directory is outside this project.")
        for source in (run_dir / "work").glob("*"):
            if source.is_file() and not source.is_symlink() and source.stat().st_size < 1024 * 1024:
                target = case / "generated-work" / source.name
                target.parent.mkdir(exist_ok=True)
                if target.exists():
                    assert sha(target) == sha(source)
                else:
                    shutil.copyfile(source, target)

    face = parent_group("parent-225-primitive-basis-elements", 225)
    rhom = parent_group("parent-160-primitive-basis-elements", 160)
    subgroup_row = next(line for line in read("subgroup-221-gm4minus-p3-centered").splitlines()
                        if line.startswith("221 "))
    # First tuple is OPD; the next three are the subgroup conventional basis.
    subgroup_basis = [vector(t) for t in tuples(subgroup_row)[1:4]]
    mapped = product(mat(subgroup_basis), mat(rhom["primitive_basis_columns"]))
    assert abs(det(mapped)) == 1
    assert all(x.denominator == 1 for row in mapped for x in row)

    alpha, beta, gamma = [math.radians(x) for x in (80, 85, 75)]
    cx = 6 * math.cos(beta)
    cy = 6 * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    lattice = ((4.0, 5 * math.cos(gamma), cx),
               (0.0, 5 * math.sin(gamma), cy),
               (0.0, 0.0, math.sqrt(36 - cx * cx - cy * cy)))
    expected_xyz = tuple(sum(lattice[i][j] * (0.123, 0.234, 0.345)[j]
                             for j in range(3)) for i in range(3))
    assert "(0.97523, 1.45392, 2.03654)" in read("gamma-1-free-sites-cartesian")
    printed_xyz = (0.97523, 1.45392, 2.03654)
    cartesian_residual = max(abs(a - b) for a, b in zip(expected_xyz, printed_xyz))
    assert cartesian_residual <= 0.5e-5

    # Each general-k real basis has cos(phi) and -sin(phi) columns. Check these
    # only for SG1's translation representation, not as a general crystal rule.
    blocks, current = [], None
    for line in read("gp-1-rational-two-sites-modes").splitlines():
        if line.startswith("1 P1"):
            current = []
            blocks.append(current)
        ts = tuples(line)
        if len(ts) == 3 and current is not None:
            current.append([tuple(float(v) for v in vector(t)) for t in ts])
    assert len(blocks) == 6 and all(len(b) == 12 for b in blocks)
    largest_residual = 0.0
    norms_squared = []
    for j, block in enumerate(blocks):
        axis = j % 3
        origin = block[0][0]
        norm = [0.0, 0.0]
        for point, first, second in block:
            shift = tuple(point[i] - origin[i] for i in range(3))
            assert max(abs(s - round(s)) for s in shift) < 1e-12
            phase = 2 * math.pi * (round(shift[0]) / 3 + round(shift[1]) / 4)
            for column, target in zip((first, second), (math.cos(phase), -math.sin(phase))):
                largest_residual = max(largest_residual, abs(column[axis] - target))
                assert all(column[k] == 0 for k in range(3) if k != axis)
            norm[0] += sum(x * x for x in first)
            norm[1] += sum(x * x for x in second)
        norms_squared.append(norm)
    assert largest_residual <= 0.5e-3  # Printed three-decimal rounding, not scientific solver tolerance.
    direction_rows = [line for line in read("direction-1-to1-twelve-cell").splitlines()
                      if line.startswith(("GM1 ", "Y1 ", "GP1GQ1 "))]
    assert len(direction_rows) == 7
    real_dimensions = [2 if row.startswith("GP1GQ1") else 1 for row in direction_rows]
    assert sum(real_dimensions) == 12

    audit = {
        "purpose": "Bounded independent checks of the retained official ISO outputs.",
        "scientific_status": "Finite contract probes only; not full workflow acceptance.",
        "sources": {
            "manual": "evidence/database/official/isotropy_doc.php.html",
            "manual_sha256": sha(ROOT / "evidence/database/official/isotropy_doc.php.html"),
            "target": "vendor/isobyu/iso",
            "target_sha256": sha(ROOT / "vendor/isobyu/iso"),
            "audit_script_sha256": sha(Path(__file__)),
        },
        "centering": [face, rhom],
        "centered_subgroup": {
            "case": "subgroup-221-gm4minus-p3-centered",
            "conventional_basis_determinant_absolute": abs(det(mat(subgroup_basis))),
            "mapped_primitive_basis_matrix": mapped,
            "mapped_primitive_determinant_absolute": abs(det(mapped)),
            "reported_size": 1,
        },
        "nonorthogonal_cartesian": {
            "case": "gamma-1-free-sites-cartesian",
            "lattice_columns_in_cartesian_rows": lattice,
            "fractional_point": [0.123, 0.234, 0.345],
            "expected_cartesian_point": expected_xyz,
            "printed_cartesian_point": printed_xyz,
            "max_abs_difference": cartesian_residual,
            "format_rounding_bound": 0.5e-5,
            "basis_norms": [4, 5, 6],
            "a_dot_b": 20 * math.cos(gamma),
            "conclusion": "Raw projections are neither unit-normalized nor necessarily orthogonal.",
        },
        "general_k_translation_basis": {
            "case": "gp-1-rational-two-sites-modes",
            "formula": "per-axis real columns cos(2*pi*k dot t), -sin(2*pi*k dot t)",
            "k_parent_conventional_reciprocal": ["1/3", "1/4", "0"],
            "blocks": len(blocks), "sites": 2, "points_per_block": 12,
            "columns_per_block": 2,
            "norm_squared_by_block": norms_squared,
            "max_abs_difference_from_formula": largest_residual,
            "format_rounding_bound": 0.5e-3,
            "applicability": "SG1 independent translation-orbit synthetic input only.",
        },
        "finite_k_enumeration": {
            "case": "direction-1-to1-twelve-cell", "rows": direction_rows,
            "real_dimensions": real_dimensions, "sum": sum(real_dimensions),
            "translation_quotient_order": 12,
            "two_sites_full_displacement_dimension": 2 * 3 * 12,
            "primary_irrep_only_dimension": 2 * 3 * 2,
        },
        "negative_or_limited_observations": [
            {"case": "builder-parent-225",
             "observation": "Documented SETTING INTERNATIONAL NEW bombs in this binary with exit zero. Initial builder probes preserve this failure; corrected default-setting probes and typed parser reject it."},
            {"case": "dt-221-quarter-subgroups",
             "observation": "Y declines database creation; Quit display; process exit zero."},
            {"case": "frequency-dt-221-quarter-p1",
             "observation": "Frequency column empty for this generated parameter-k subgroup."},
            {"case": "frequency-gp-1-rational-c1",
             "observation": "Frequency column empty despite nonempty fixed displacement space."},
            {"case": "negative-invalid-cell-gp",
             "observation": "Incompatible display window I returns only one parent cell without an error; insufficient for periodic export."},
            {"case": "gamma-221-forbidden-scalar-irrep",
             "observation": "GM1+ polar vector at 1a is empty without error; absence alone does not distinguish valid zero from a failed query."},
            {"case": "gamma-221-multisite",
             "observation": "SHOW CARTESIAN scales displacements but special-site Points remain exact fractional-looking values (b=(1/2,1/2,1/2) at lattice a=4). Do not assume all printed Point units from one toggle."},
        ],
    }
    (BASE / "audit.json").write_text(json.dumps(audit, indent=2, default=str) + "\n")
    index = []
    for case in sorted(BASE.iterdir()):
        if case.is_dir() and (case / "evidence.json").is_file():
            evidence = json.loads((case / "evidence.json").read_text())
            index.append({"case": case.name, "run_id": evidence["run_id"],
                          "status": evidence["status"],
                          "scientific_validation": "not_performed",
                          "files": {p.relative_to(case).as_posix(): sha(p)
                                    for p in sorted(case.rglob("*")) if p.is_file()}})
    (BASE / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(json.dumps({"cases_indexed": len(index), "checks": "passed",
                      "audit": str(BASE / "audit.json")}))


if __name__ == "__main__":
    main()
