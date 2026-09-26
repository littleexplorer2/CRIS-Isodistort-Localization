"""Print the current complete-mode basis for one Method 1/2 candidate.

This is a read-only diagnostic companion to ``validate_method_outputs.py``.
It deliberately selects by the complete candidate identity rather than by a
saved directory name, so equivalent settings and duplicate OPD labels remain
unambiguous.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SpaceGroup

PROJECT = Path(__file__).resolve().parents[2]
WORKSPACE = PROJECT.parent
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))

from isocore.api import IsoDistort  # noqa: E402
from isocore.distortion.superspace import (  # noqa: E402
    _child_origin_candidates,
    _displacement_actions,
    _map_smodes_block_to_supercell,
    _parent_map,
    _project_invariant,
)
from isocore.structure import build_supercell  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cif", type=Path)
    parser.add_argument("irrep")
    parser.add_argument("opd")
    parser.add_argument("--k-point")
    parser.add_argument("--k-parameter", action="append", default=[])
    parser.add_argument("--raw-smodes", action="store_true")
    parser.add_argument("--origin-scores", action="store_true")
    args = parser.parse_args()

    api = IsoDistort()
    api.load_structure(args.cif.expanduser().resolve())
    if args.k_point:
        candidates = api.list_subgroups_at_kpoint(
            args.k_point,
            k_parameters=list(args.k_parameter),
            generate_if_missing=False,
        )
    else:
        rows = api.search_method_1(distortion_types=["strain", "displacive"])
        candidates = [row.subgroup for row in rows]
    matches = [
        subgroup for subgroup in candidates
        if subgroup.irrep_label == args.irrep and subgroup.opd_symbol == args.opd
    ]
    if len(matches) != 1:
        print(f"expected one candidate, found {len(matches)}", file=sys.stderr)
        return 2
    subgroup = matches[0]
    if args.raw_smodes:
        payload = api._smodes._build_input(
            api.structure,
            int(api.symmetry_info["space_group_number"]),
            list(api.symmetry_info.get("wyckoff_sites") or []),
            [(subgroup.k_point_label, list(subgroup.k_parameters or []))],
        )
        print(api._smodes.run_stdin(api._smodes.binary, payload))
        return 0
    if args.origin_scores:
        basis = np.asarray(subgroup.basis_vectors, dtype=float)
        child = build_supercell(api.structure, basis)
        _parent_idx, translations = _parent_map(api.structure, child, basis)
        origins = _child_origin_candidates(subgroup.origin, basis, translations)
        blocks = api._smodes.compute_mode_blocks(
            api.structure,
            int(api.symmetry_info["space_group_number"]),
            list(api.symmetry_info.get("wyckoff_sites") or []),
            [(subgroup.k_point_label, list(subgroup.k_parameters or []))],
        ) or []
        primary = [block for block in blocks if block.irrep == subgroup.irrep_label]
        standard = list(
            SpaceGroup.from_int_number(int(subgroup.space_group_number)).symmetry_ops
        )
        for origin_i, origin in enumerate(origins):
            for sign in (1.0, -1.0, 0.0):
                shift = sign * origin
                ops = []
                for op in standard:
                    rotation = np.asarray(op.rotation_matrix, dtype=float)
                    translation = np.asarray(op.translation_vector, dtype=float)
                    ops.append(SymmOp.from_rotation_and_translation(
                        rotation,
                        translation + shift - rotation @ shift,
                    ))
                actions = _displacement_actions(child, ops)
                projected = []
                if actions is not None:
                    for block in primary:
                        mapped = _map_smodes_block_to_supercell(block, child)
                        result = None if mapped is None else _project_invariant(
                            mapped, child, ops, actions=actions,
                        )
                        projected.append(
                            0.0 if result is None else float(np.linalg.norm(result))
                        )
                print(
                    f"origin[{origin_i}]={origin.tolist()} sign={sign:+.0f} "
                    f"complete={actions is not None} primary={projected}"
                )
        return 0
    api.search_method_2(
        subgroup.index,
        distortion_type=["displacive"],
        number_of_independent_modulations=0,
        candidates=candidates,
    )
    print(f"{args.irrep} {args.opd}: {len(api.mode_displacements)} modes")
    for key, label in api._mode_labels_now().items():
        print(f"{key}\t{label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
