"""Read-only independent comparison; imports no product calculation modules."""
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import html
import json
import math
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence/phase2/web-comparison/20261020_cubic_polar_axis"
LOCAL = ROOT / "evidence/phase2/workflow/20261020T145121147627Z/cubic_polar_axis"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_from_input(text):
    values = {}
    for kind, names in (("length", "abc"), ("angle", ("alpha", "beta", "gamma"))):
        for name in names:
            tag = f"_cell_{kind}_{name}"
            match = re.search(r"(?m)^" + tag + r"\s+(\S+)\s*$", text)
            if not match:
                raise ValueError(f"Missing unambiguous input cell field {tag}")
            values[name] = float(match[1])
    a, b, c = (values[x] for x in "abc")
    ca, cb, cg = (math.cos(math.radians(values[x])) for x in ("alpha", "beta", "gamma"))
    sg = math.sin(math.radians(values["gamma"]))
    cy = (ca - cb * cg) / sg
    return np.array([[a, b * cg, c * cb], [0, b * sg, c * cy],
                     [0, 0, c * math.sqrt(1 - cb * cb - cy * cy)]])


def rational_triplet(text):
    return tuple(F(value.strip()) for value in text.split(","))


def main():
    freeze = json.loads((OUT / "local_freeze_manifest.json").read_text())
    for name, expected in freeze["files"].items():
        path = ROOT / name
        assert sha(path) == expected["sha256"] and path.stat().st_size == expected["bytes"]
    local = json.loads((LOCAL / "exports/modes.json").read_text())
    raw = (OUT / "05_complete_modes.html").read_text(encoding="utf-8")
    pre = re.findall(r"<pre>(.*?)</pre>", raw, re.I | re.S)
    assert len(pre) == 1
    text = html.unescape(re.sub(r"<[^>]*>", "", pre[0]))
    subgroup = re.search(
        r"Subgroup details\s+(\d+)\s+(\S+), basis=\{([^}]+)\}, origin=\(([^)]+)\), s=(\d+), i=(\d+)",
        text,
    )
    assert subgroup
    columns = [rational_triplet(v) for v in re.findall(r"\(([^)]+)\)", subgroup[3])]
    basis = tuple(tuple(columns[j][i] for j in range(3)) for i in range(3))
    origin = rational_triplet(subgroup[4])
    website_subgroup = {
        "number": int(subgroup[1]), "symbol": subgroup[2],
        "basis": [[str(x) for x in row] for row in basis], "origin": list(map(str, origin)),
        "primitive_size": int(subgroup[5]), "index": int(subgroup[6]),
    }
    assert website_subgroup == local["subgroup"], (website_subgroup, local["subgroup"])
    lattice = matrix_from_input((LOCAL / "parent.cif").read_text()) @ np.array(basis, float)
    assert np.max(np.abs(lattice - local["lattice"])) < 1e-12

    positions = {
        tuple(F(value) % 1 for value in atom["fractional_exact"]): i
        for i, atom in enumerate(local["atoms"])
    }
    assert len(positions) == len(local["atoms"])
    section = text.split("Displacive mode definitions", 1)[1].split("Displacive mode amplitudes", 1)[0]
    modes, current = [], None
    for line in section.splitlines():
        header = re.match(r"(\S+)\s+normfactor\s*=\s*(\S+)\s*$", line)
        if header:
            current = {
                "label": header[1], "printed_normfactor": float(header[2]),
                "fractional_vector": np.zeros((len(positions), 3)),
                "rows": [], "seen_positions": set(),
            }
            modes.append(current)
            continue
        parts = line.split()
        if current is None or len(parts) not in (6, 7):
            continue
        data = parts[-6:]
        try:
            coordinates = tuple(F(value) % 1 for value in data[:3])
            vector = tuple(float(value) for value in data[3:])
        except ValueError:
            continue
        assert coordinates in positions, coordinates
        assert coordinates not in current["seen_positions"]
        current["seen_positions"].add(coordinates)
        atom_index = positions[coordinates]
        current["fractional_vector"][atom_index] = vector
        current["rows"].append({
            "printed_atom_label": parts[0] if len(parts) == 7 else None,
            "position": list(map(str, coordinates)), "local_atom_index": atom_index,
            "fractional_vector": vector,
        })
    assert modes and all(mode["rows"] for mode in modes)
    raw_cartesian = np.column_stack([(mode["fractional_vector"] @ lattice.T).ravel() for mode in modes])
    # QR changes only basis within the website's observed span. It never rewrites local output.
    q_web, triangular = np.linalg.qr(raw_cartesian, mode="reduced")
    bound = 2048 * np.finfo(float).eps * max(raw_cartesian.shape)
    assert min(abs(np.diag(triangular))) > bound
    local_basis = np.array(local["displacement_basis"], float)
    projector_error = float(np.max(np.abs(q_web @ q_web.T - local_basis @ local_basis.T)))
    singular_values = np.linalg.svd(local_basis.T @ q_web, compute_uv=False)
    dimension_agreement = local_basis.shape[1] == len(modes)
    norms = np.linalg.norm(raw_cartesian, axis=0)
    normalization = []
    for i, mode in enumerate(modes):
        calculated = 1 / norms[i]
        difference = abs(mode["printed_normfactor"] - calculated)
        normalization.append({
            "label": mode["label"], "printed_normfactor_angstrom_inverse": mode["printed_normfactor"],
            "calculated_from_printed_vector_and_input_cell": float(calculated),
            "absolute_difference": float(difference),
            "five_decimal_half_unit": 5e-6,
            "consistent_with_printed_precision": bool(difference <= 5e-6 + bound),
            "norm_using_printed_factor": float(norms[i] * mode["printed_normfactor"]),
        })
    outcome = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "same_input_finite_displacement_span_agrees" if (
            dimension_agreement and projector_error < bound) else "inconclusive_or_mismatch",
        "scope": "One newly uploaded synthetic ordinary 3D GM4-/P1 case; no all-case or strain-equivalence claim",
        "isolation": {
            "freeze_manifest": "local_freeze_manifest.json", "freeze_sha256": sha(OUT / "local_freeze_manifest.json"),
            "all_local_files_unchanged": True, "product_imports": False,
            "reference_role": "Read only after local output freeze; no production changes or reruns",
        },
        "inputs_read": {
            str(path.relative_to(ROOT)): sha(path) for path in (
                LOCAL / "parent.cif", LOCAL / "exports/modes.json", OUT / "05_complete_modes.html",
                OUT / "local_freeze_manifest.json",
            )
        },
        "comparator_sha256": sha(Path(__file__)), "numpy_version": np.__version__,
        "subgroup_exact_comparison": {"agrees": True, "website": website_subgroup},
        "coordinate_conversion": "Website complete-mode fractional rows -> input-parent A * website C -> Cartesian; rows matched by exact modulo-one coordinates",
        "website_mode_count": len(modes), "local_mode_count": local_basis.shape[1],
        "local_internal_count": sum(name.startswith("internal_") for name in local["mode_names"]),
        "local_translation_count": sum(name.startswith("translation_") for name in local["mode_names"]),
        "website_rows": [{"label": mode["label"], "rows": mode["rows"]} for mode in modes],
        "span_comparison": {
            "max_projector_component_difference": projector_error,
            "cross_basis_singular_values": singular_values.tolist(),
            "absolute_numerical_bound": float(bound),
            "bound_reason": "2048*binary64 epsilon*max(matrix shape) for this 15x4 well-conditioned example; printed displacement entries are exactly 0 or 1, so no rounded-vector uncertainty is inferred",
            "website_cartesian_raw_column_norms": norms.tolist(),
        },
        "normfactor_comparison": normalization,
        "limitations": [
            "No column order, signs or IR label equivalence claimed for the local basis.",
            "Official normfactors are rounded to five decimals; directly applying them need not give exact unit norm.",
            "Only displacement span compared; two printed strain modes retained as reference but not used to establish current official strain convention.",
            "HTTPS form replay sees request/response bodies, not server executable calls.",
            "The frozen local outputs precede the subsequent CIF-roundtrip export fix and are not presented as current repaired-export validation.",
        ],
    }
    target = OUT / "comparison.json"
    if target.exists():
        raise FileExistsError(target)
    target.write_text(json.dumps(outcome, indent=2) + "\n", encoding="utf-8")
    print(outcome["status"], "modes", len(modes), "projector max error", projector_error)


if __name__ == "__main__":
    main()
