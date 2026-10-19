"""Tests for support_validate_method_outputs; migrated without changing scientific expectations."""
from __future__ import annotations

from pathlib import Path


# Original coverage: test_method_validation.py.

def _complete_modes_html(
    mode_count: int = 1,
    *,
    subgroup_details: str | None = None,
) -> str:
    modes = "".join(
        f"<pre>:dsp] mode-{index} normfactor = 1</pre>"
        for index in range(mode_count)
    )
    return (
        "<html><head><title>ISODISTORT: complete modes details</title></head>"
        "<body><h1>ISODISTORT: complete modes details</h1>"
        "<strong>Subgroup details</strong>"
        f"{subgroup_details or ''}"
        "<strong>Undistorted superstructure</strong>"
        "<strong>Distorted superstructure</strong>"
        f"{modes}</body></html>"
    )

def _write_setting_cif(
    path: Path,
    *,
    basis: str = "(1,0,0),(0,1,0),(0,0,1)",
    origin: str = "(0,0,0)",
    cell: tuple[str, str, str, str, str, str] = (
        "1.00000",
        "1.00000",
        "1.00000",
        "90.00000",
        "90.00000",
        "90.00000",
    ),
    species: str = "C",
    coordinates: tuple[str, str, str] = ("0.10000", "0.20000", "0.30000"),
    occupancy: str = "1.00000",
    atom_rows: tuple[tuple[str, str, str, str, str, str], ...] | None = None,
    symmetry_operations: tuple[str, ...] | None = ("x,y,z",),
    magnetic_moments: tuple[tuple[str, str, str], ...] | None = None,
    space_group_symbol: str = "P 1",
    space_group_number: str | None = "1",
) -> None:
    a, b, c, alpha, beta, gamma = cell
    x, y, z = coordinates
    number_line = (
        f"_symmetry_Int_Tables_number {space_group_number}\n"
        if space_group_number is not None
        else ""
    )
    if atom_rows is None:
        atom_rows = ((f"{species}1", species, x, y, z, occupancy),)
    if magnetic_moments is not None and len(magnetic_moments) != len(atom_rows):
        raise ValueError("magnetic moment rows must align with atom rows")
    symmetry_loop = ""
    if symmetry_operations is not None:
        symmetry_loop = (
            "loop_\n_symmetry_equiv_pos_as_xyz\n"
            + "\n".join(f"'{operation}'" for operation in symmetry_operations)
            + "\n"
        )
    moment_columns = ""
    if magnetic_moments is not None:
        moment_columns = (
            "_atom_site_moment.crystalaxis_x\n"
            "_atom_site_moment.crystalaxis_y\n"
            "_atom_site_moment.crystalaxis_z\n"
        )
    atom_lines = []
    for index, row in enumerate(atom_rows):
        fields = list(row)
        if magnetic_moments is not None:
            fields.extend(magnetic_moments[index])
        atom_lines.append(" ".join(fields))
    path.write_text(
        f"""# IR: GM1+
# P1 (a) 1 P1, basis={{{basis}}}, origin={origin}, s=1, i=1, k-active= (0,0,0)
data_test
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha {alpha}
_cell_angle_beta {beta}
_cell_angle_gamma {gamma}
_symmetry_space_group_name_H-M '{space_group_symbol}'
{number_line}{symmetry_loop}loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
{moment_columns}{chr(10).join(atom_lines)}
""",
        encoding="utf-8",
    )
