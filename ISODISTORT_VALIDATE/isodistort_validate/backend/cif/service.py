"""Read-only CIF parsing and file fingerprints."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from pymatgen.core import Structure
from pymatgen.io.cif import CifFile, CifParser


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig")


def parse_structure(path: Path) -> Structure:
    parser = CifParser(str(path), occupancy_tolerance=100)
    structures = parser.parse_structures(primitive=False)
    if not structures:
        raise ValueError(f"CIF does not contain a parseable structure: {path}")
    if len(structures) > 1:
        raise ValueError(f"CIF contains {len(structures)} structures; expected one: {path}")
    return structures[0]


def read_block(path: Path) -> dict[str, Any]:
    cif = CifFile.from_file(path)
    if len(cif.data) != 1:
        raise ValueError(f"CIF must contain exactly one data block: {path}")
    return next(iter(cif.data.values())).data
