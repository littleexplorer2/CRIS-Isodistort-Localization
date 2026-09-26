"""Read-only SMODES evidence for real quadrature partners at parametric k.

This diagnostic does not select or patch a production answer.  It records the
raw SMODES block degeneracy, whether ``-k`` is parent-reciprocal equivalent to
``k``, and whether the generic phase-translation construction yields a second
independent real Cartesian column.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from isocore.api import IsoDistort  # noqa: E402
from isocore.distortion.superspace import (  # noqa: E402
    _collect_smodes_blocks,
    _is_self_conjugate_parent_k,
    _map_smodes_block_to_supercell,
    _primary_candidate,
    _quadrature_partner,
    folding_k_candidates,
    parse_k_token,
    smodes_k_spec,
)
from isocore.structure import build_supercell  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def diagnose(args: argparse.Namespace) -> dict:
    api = IsoDistort(language="en")
    api.load_structure(args.cif.resolve())
    pool = api.list_subgroups_at_kpoint(
        args.k_point,
        k_parameters=list(args.k_parameter),
        generate_if_missing=False,
    )
    matches = [
        item
        for item in pool
        if item.irrep_label == args.irrep
        and item.opd_symbol == args.opd
        and (args.space_group is None or int(item.space_group_number) == args.space_group)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one selected subgroup, found {len(matches)}")
    subgroup = matches[0]
    parent_sg = int(api.symmetry_info["space_group_number"])
    basis = np.asarray(subgroup.basis_vectors, dtype=float)
    kpoints = api._iso.list_k_points(parent_sg)
    folded = folding_k_candidates(
        parent_sg,
        kpoints,
        basis,
        _primary_candidate(subgroup, parent_sg),
        0,
    )
    specs = []
    seen = set()
    for candidate in folded:
        active = candidate.active_coords if candidate.active_coords is not None else candidate.coords
        key = (
            candidate.label,
            tuple(candidate.params or ()),
            tuple(str(value) for value in active),
        )
        if key in seen:
            continue
        seen.add(key)
        specs.append((candidate, smodes_k_spec(parent_sg, candidate, kpoints)))
    blocks = _collect_smodes_blocks(
        api._smodes,
        api.structure,
        parent_sg,
        list(api.symmetry_info.get("wyckoff_sites") or []),
        specs,
    )
    supercell = build_supercell(api.structure, basis)

    def representative_parent_arrow(array: np.ndarray | None) -> list[float] | None:
        if array is None or not array.size:
            return None
        parent_arrows = np.asarray(array, dtype=float) @ basis
        arrow = parent_arrows[int(np.argmax(np.linalg.norm(parent_arrows, axis=1)))]
        scale = float(np.max(np.abs(arrow)))
        if scale <= 1.0e-12:
            return [0.0, 0.0, 0.0]
        arrow = arrow / scale
        first = next((value for value in arrow if abs(float(value)) > 1.0e-12), 1.0)
        if first < 0:
            arrow = -arrow
        return [round(float(value), 10) for value in arrow]

    rows = []
    for block in blocks:
        if args.block_irrep and block.irrep != args.block_irrep:
            continue
        mapped = _map_smodes_block_to_supercell(block, supercell)
        phase = block.active_k_coords or block.k_coords
        try:
            kvec = np.array([parse_k_token(value) for value in phase], dtype=float)
        except (ValueError, ZeroDivisionError):
            continue
        partner = (
            None
            if mapped is None
            else _quadrature_partner(
                mapped,
                supercell,
                kvec,
                basis,
                parent_sg=parent_sg,
            )
        )
        active_species = sorted(
            {
                atom.species
                for atom, vector in zip(
                    block.atoms,
                    block.dense_displacements(),
                    strict=True,
                )
                if float(np.linalg.norm(vector)) > 1.0e-10
            }
        )
        rank = 0
        if mapped is not None:
            columns = [(mapped @ supercell.lattice.matrix).ravel()]
            if partner is not None:
                columns.append((partner @ supercell.lattice.matrix).ravel())
            rank = int(np.linalg.matrix_rank(np.column_stack(columns), tol=1.0e-8))
        rows.append(
            {
                "irrep": block.irrep,
                "reported_degeneracy": int(block.degeneracy),
                "k_label": block.k_label,
                "k_coords": list(block.k_coords),
                "active_k_coords": list(block.active_k_coords or ()),
                "self_conjugate_mod_parent_reciprocal_lattice": (
                    _is_self_conjugate_parent_k(kvec, parent_sg)
                ),
                "active_species": active_species,
                "quadrature_partner_exists": partner is not None,
                "seed_plus_partner_cartesian_rank": rank,
                "seed_representative_parent_arrow": representative_parent_arrow(mapped),
                "partner_representative_parent_arrow": representative_parent_arrow(partner),
            }
        )
    signed_paths = [
        Path(__file__).resolve(),
        args.cif.resolve(),
        PROJECT_ROOT / "isocore" / "distortion" / "superspace.py",
        PROJECT_ROOT / "isocore" / "backend" / "smodes_wrapper.py",
    ]
    signed_files = [
        {"path": str(path), "sha256": _sha256(path), "bytes": path.stat().st_size}
        for path in signed_paths
    ]
    signature_material = json.dumps(
        {"files": signed_files, "rows": rows},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        "generated_at": datetime.now(UTC).isoformat(),
        "evidence_signature": hashlib.sha256(signature_material).hexdigest(),
        "signed_files": signed_files,
        "question": (
            "Can a non-self-conjugate 1D complex SMODES block supply two "
            "independent real cosine/sine columns?"
        ),
        "parent": str(args.cif.resolve()),
        "parent_space_group": parent_sg,
        "selected_subgroup": {
            "irrep": subgroup.irrep_label,
            "opd": subgroup.opd_symbol,
            "space_group_number": int(subgroup.space_group_number),
            "basis": np.asarray(subgroup.basis_vectors, dtype=float).tolist(),
        },
        "block_filter": args.block_irrep,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    parser.add_argument("k_point")
    parser.add_argument("irrep")
    parser.add_argument("opd")
    parser.add_argument("--k-parameter", action="append", default=[])
    parser.add_argument("--space-group", type=int)
    parser.add_argument("--block-irrep")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        payload = diagnose(args)
    except Exception as exc:  # noqa: BLE001 - CLI records backend diagnostic failure
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    rendered = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f"{args.output.name}.tmp")
        temporary.write_text(rendered, encoding="utf-8")
        temporary.replace(args.output)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
