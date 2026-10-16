"""Orchestrate file validation, parsing, semantic comparison, and hashes."""

from __future__ import annotations

from pathlib import Path
from math import isfinite

from ...backend.cif import parse_structure, read_block, read_text, sha256
from ...backend.comparison import compare_structures
from ...backend.results import ComparisonResult
from ...backend.config import get_config


def compare_cif(
    local_path: str | Path,
    reference_path: str | Path,
    lattice_tolerance: float | None = None,
    coordinate_tolerance: float | None = None,
    scalar_tolerance: float | None = None,
    ignore_atom_order: bool | None = None,
    reference_sha256: str | None = None,
) -> ComparisonResult:
    """Compare two CIF files at byte, text, metadata, and structure levels."""
    cfg = get_config()
    lattice_tolerance = cfg.lattice_tolerance if lattice_tolerance is None else lattice_tolerance
    coordinate_tolerance = cfg.coordinate_tolerance if coordinate_tolerance is None else coordinate_tolerance
    scalar_tolerance = cfg.scalar_tolerance if scalar_tolerance is None else scalar_tolerance
    ignore_atom_order = cfg.ignore_atom_order if ignore_atom_order is None else ignore_atom_order
    local = Path(local_path).expanduser().resolve()
    reference = Path(reference_path).expanduser().resolve()
    for path in (local, reference):
        if not path.is_file():
            raise FileNotFoundError(f"CIF file not found: {path}")
    if any(not isfinite(value) or value < 0 for value in
           (lattice_tolerance, coordinate_tolerance, scalar_tolerance)):
        raise ValueError("tolerances must be finite and non-negative")

    local_bytes = local.read_bytes()
    reference_bytes = reference.read_bytes()
    local_sha = sha256(local)
    reference_sha = sha256(reference)
    local_structure = parse_structure(local)
    reference_structure = parse_structure(reference)
    local_block = read_block(local)
    reference_block = read_block(reference)
    structure_equal, details = compare_structures(
        local_structure,
        reference_structure,
        local_block,
        reference_block,
        lattice_tolerance,
        coordinate_tolerance,
        scalar_tolerance,
        ignore_atom_order,
    )
    hash_valid = None if reference_sha256 is None else reference_sha.lower() == reference_sha256.lower()
    if hash_valid is False:
        details["issues"].append("reference SHA-256 does not match expected hash")

    return ComparisonResult(
        local_path=str(local),
        reference_path=str(reference),
        byte_exact=local_bytes == reference_bytes,
        text_exact=read_text(local) == read_text(reference),
        structure_equal=structure_equal,
        reference_hash_valid=hash_valid,
        local_sha256=local_sha,
        reference_sha256=reference_sha,
        details=details,
    )
