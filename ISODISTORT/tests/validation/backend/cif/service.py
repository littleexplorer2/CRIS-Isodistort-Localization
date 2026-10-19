"""Read-only CIF parsing and file fingerprints."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pymatgen.core import Structure
from pymatgen.io.cif import CifFile, CifParser

from ..provenance import FileSnapshot, sha256_file


def sha256(path: Path) -> str:
    return sha256_file(path)


def _snapshot(source: Path | FileSnapshot) -> FileSnapshot:
    return source if isinstance(source, FileSnapshot) else FileSnapshot.capture(source)


def read_text(source: Path | FileSnapshot) -> str:
    return _snapshot(source).text


def parse_structure(source: Path | FileSnapshot) -> Structure:
    snapshot = _snapshot(source)
    parser = CifParser.from_str(snapshot.text, occupancy_tolerance=100)
    structures = parser.parse_structures(primitive=False)
    if not structures:
        raise ValueError(f"CIF does not contain a parseable structure: {snapshot.path}")
    if len(structures) > 1:
        raise ValueError(f"CIF contains {len(structures)} structures; expected one: {snapshot.path}")
    return structures[0]


def read_block(source: Path | FileSnapshot) -> dict[str, Any]:
    snapshot = _snapshot(source)
    cif = CifFile.from_str(snapshot.text)
    if len(cif.data) != 1:
        raise ValueError(f"CIF must contain exactly one data block: {snapshot.path}")
    return next(iter(cif.data.values())).data
