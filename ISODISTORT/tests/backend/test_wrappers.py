"""Tests for backend.wrappers; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.models.iso_mode_models import InvariantDirection
from backend.models.iso_mode_models import MicroscopicVectorBlock
from backend.models.iso_mode_models import MicroscopicVectorRow
from backend.utils import OutputParseError as _formats_OutputParseError
from backend.utils import WrapperRunError as _basics_WrapperRunError
from backend.utils import WrapperRunError as _native_wrapper_WrapperRunError
from backend.utils import WrapperRunError as _wrapper_runtime_WrapperRunError
from backend.utils import WrapperTimeoutError
from backend.utils import format_opd_line
from backend.utils import parse_subgroup_table
from backend.utils.exceptions import OutputParseError as _microscopic_column_provenance_OutputParseError
from backend.utils.exceptions import WrapperRunError as _site_irreps_WrapperRunError
from backend.utils.text_parser import parse_domain_table
from backend.wrappers import FindsymResult
from backend.wrappers import FindsymWrapper
from backend.wrappers import IsoWrapper as _search_methods_IsoWrapper
from backend.wrappers import IsoWrapper as _wsl_IsoWrapper
from backend.wrappers import SubgroupInfo as _formats_SubgroupInfo
from backend.wrappers import SubgroupInfo as _search_methods_SubgroupInfo
from backend.wrappers import SubgroupInfo as _wsl_SubgroupInfo
from backend.wrappers import base_wrapper as module
from backend.wrappers.base_wrapper import BaseWrapper
from backend.wrappers.base_wrapper import _decode_process_output
from backend.wrappers.base_wrapper import _is_wsl_access_denied
from backend.wrappers.iso_wrapper import IsoWrapper as _isotropy_cache_IsoWrapper
from backend.wrappers.iso_wrapper import IsoWrapper as _microscopic_column_provenance_IsoWrapper
from backend.wrappers.iso_wrapper import IsoWrapper as _site_irreps_IsoWrapper
from backend.wrappers.iso_wrapper import KPointInfo
from backend.wrappers.iso_wrapper import SubgroupInfo as _microscopic_column_provenance_SubgroupInfo
from backend.wrappers.iso_wrapper import SubgroupInfo as _site_irreps_SubgroupInfo
from backend.wrappers.isotropy_cache import _safe_name
from backend.wrappers.isotropy_cache import delete_isotropy_cache
from backend.wrappers.isotropy_cache import list_isotropy_cache
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from pathlib import PurePosixPath
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.support_application_session import _wsl_required
from tests.support_backend_export import _formats_REPO
from tests.support_backend_method2 import _block
from tests.support_backend_method2 import _exact_fixed_space_subgroup
from tests.support_backend_method2 import _exact_fixed_space_wrapper
from tests.support_backend_method2 import _microscopic_column_provenance_subgroup
from tests.support_backend_method2 import _provenance
from types import SimpleNamespace
import gc
import numpy as np
import os
import pytest
import re
import shlex
import subprocess
import sys
import tempfile
import time


# Original coverage: test_basics.py.

def test_binary_launch_requires_the_staging_directory(monkeypatch):
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._stage_dir = "/var/tmp/id-12345678"
    monkeypatch.setattr(
        wrapper,
        "_stage_text",
        lambda _prefix, _text: "/var/tmp/id-12345678/iso_input.in",
    )
    monkeypatch.setattr(wrapper, "_wsl_bin_path", lambda _path: "/opt/isobyu/iso")
    monkeypatch.setattr(wrapper, "_isodata_path", lambda: "/var/tmp/id-12345678/data/")
    commands: list[str] = []

    def capture(command: str, timeout: float | None = None):
        _ = timeout
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(wrapper, "_wsl", capture)

    assert wrapper._run_program(Path("iso"), "QUIT\n", True) == "ok"
    assert len(commands) == 1
    assert (
        "export ISODATA=/var/tmp/id-12345678/data/ && "
        "cd /var/tmp/id-12345678 && /opt/isobyu/iso"
    ) in commands[0]
    assert "; cd " not in commands[0]

def test_wsl_stream_decoder_separates_linux_output_from_windows_diagnostics():
    assert _decode_process_output(b"/home/test\n") == "/home/test"
    diagnostic = "Access is denied.\r\nError code: Wsl/E_ACCESSDENIED\r\n"
    decoded = _decode_process_output(diagnostic.encode("utf-16-le"))
    assert decoded == "Access is denied.\nError code: Wsl/E_ACCESSDENIED"
    assert _is_wsl_access_denied(decoded)

def test_wsl_workspace_reports_onedrive_process_restriction(monkeypatch):
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    denied = "Access is denied.\nError code: Wsl/E_ACCESSDENIED"
    monkeypatch.setattr(
        wrapper,
        "_wsl",
        lambda _command: SimpleNamespace(returncode=0xFFFFFFFF, stdout=denied, stderr=""),
    )

    with pytest.raises(_basics_WrapperRunError) as excinfo:
        wrapper._init_wsl_environment()

    message = str(excinfo.value)
    assert "Wsl/E_ACCESSDENIED" in message
    assert ".\\run_cris.ps1 <原 Python 参数>" in message
    assert "不是缺少管理员权限" in message

# Original coverage: test_formats.py.

_OFFICIAL_HTML = (
    _formats_REPO / "webpage_info" / "a. ISODISTORT_ order parameter direction.html"
)

def test_findsym_parser_retains_standard_wyckoff_form_and_signed_parameter():
    raw = """\
Space Group: 139  D4h-17    I4/mmm
Origin at    0.00000   0.00000   0.00000
Vectors a,b,c:
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000
   0.00000   0.00000   1.00000
Atomic positions and occupancies in terms of a,b,c:
Wyckoff position e (11), z = -0.13885
    1   0.00000   0.00000   0.86115     1.00000
------------------------------------------
# CIF file created by FINDSYM, version 7.1.3
data_findsym-output
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_symmetry_multiplicity
_atom_site_Wyckoff_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_fract_symmform
11 1 4 e 0.0000000000 0.0000000000 0.8611500000 1.0000000000 0,0,Dz
"""
    result = FindsymWrapper._parse_output(raw)
    assert result.origin == (0.0, 0.0, 0.0)
    assert np.allclose(result.basis_vectors, np.eye(3))
    assert result.wyckoff_sites[0]["parameters"] == {"z": -0.13885}
    assert result.standardized_sites == [{
        "atom_type": "1",
        "multiplicity": 4,
        "wyckoff_letter": "e",
        "representative_frac_coords": [0.0, 0.0, 0.86115],
        "symmform": "0,0,z",
        "parameters": {"z": -0.13885},
        "input_atom_indices": [0],
    }]

def test_findsym_orbit_metadata_rejects_a_missing_setting_transform():
    structure = Structure(Lattice.cubic(4), ["H"], [[0, 0, 0]])
    wrapper = object.__new__(FindsymWrapper)
    wrapper.identify = lambda *args, **kwargs: FindsymResult(  # type: ignore[method-assign]
        space_group_number=1,
        standardized_sites=[{
            "atom_type": "1",
            "multiplicity": 1,
            "wyckoff_letter": "a",
            "representative_frac_coords": [0.0, 0.0, 0.0],
            "symmform": "x,y,z",
            "parameters": {"x": 0.0, "y": 0.0, "z": 0.0},
        }],
    )

    with pytest.raises(_formats_OutputParseError, match="basis transformation"):
        wrapper.standardize_wyckoff_orbits(
            structure,
            [{
                "wyckoff_letter": "a",
                "multiplicity": 1,
                "equivalent_indices": [0],
                "orbit_id": "p1:a:h",
            }],
            expected_space_group=1,
        )

def test_findsym_symmform_supports_implicit_crystallographic_products():
    assert FindsymWrapper._absolute_symmform(
        [0.2, 0.05, 0.3],
        "2Dx,1/2Dx,Dz",
        {"x": 0.1, "z": 0.3},
    ) == "2x,1/2x,z"

def test_findsym_parser_matches_reordered_repeated_letters_by_output_label():
    raw = """\
Space Group: 139  D4h-17    I4/mmm
Origin at    0.00000   0.00000   0.00000
Vectors a,b,c:
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000
   0.00000   0.00000   1.00000
Atomic positions and occupancies in terms of a,b,c:
Wyckoff position e (11), z = 0.125
    1   0.00000   0.00000   0.12500     1.00000
Wyckoff position e (21), z = 0.375
    2   0.00000   0.00000   0.37500     1.00000
------------------------------------------
# CIF file created by FINDSYM, version 7.1.3
data_findsym-output
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_symmetry_multiplicity
_atom_site_Wyckoff_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_fract_symmform
21 2 4 e 0.0000000000 0.0000000000 0.3750000000 1.0000000000 0,0,Dz
11 1 4 e 0.0000000000 0.0000000000 0.1250000000 1.0000000000 0,0,Dz
"""
    result = FindsymWrapper._parse_output(raw)
    by_type = {site["atom_type"]: site for site in result.standardized_sites}
    assert by_type["1"]["parameters"] == {"z": 0.125}
    assert by_type["2"]["parameters"] == {"z": 0.375}

def test_parse_subgroup_table_keeps_letter_opd_and_raw_basis():
    text = (
        "Subgroup  Max Index Size Dir         Basis Vectors                Origin\n"
        "139 I4/mmm yes 1     1    P1  (a) (1,0,0),(0,1,0),(0,0,1) (0,0,0)\n"
        "123 P4/mmm yes 4     4    P1  (a;a) (1,-1,0),(1,1,0),(0,0,1) (0,0,0)\n"
        " 12 C2/m  no  4     1    P1  (a,0) (0,1,1),(-1,0,0),(0,-1,0) (1/2,-1/2,0)\n"
    )
    rows = parse_subgroup_table(text)
    assert len(rows) == 3
    assert rows[0]["opd_dir_raw"] == "(a)"
    assert rows[0]["basis_raw"] == "(1,0,0),(0,1,0),(0,0,1)"
    assert rows[0]["origin_raw"] == "(0,0,0)"
    assert rows[1]["opd_dir_raw"] == "(a;a)"
    assert rows[2]["origin_raw"] == "(1/2,-1/2,0)"
    assert rows[2]["origin"] == [0.5, -0.5, 0.0]

def test_subgroup_from_row_uses_iso_basis_raw():
    """Display tokens come from iso output, not a memorized official OPD table."""
    from backend.wrappers.iso_wrapper import IsoWrapper

    row = {
        "space_group_number": 12,
        "space_group_symbol": "C2/m",
        "is_maximal": True,
        "subgroup_index": 4,
        "size": 1,
        "opd_symbol": "P1",
        "opd_vector": [1.0, 0.0],
        "opd_dir_raw": "(a,0)",
        "basis_vectors": [[0, 1, 1], [1, 0, 0], [0, 0, -1]],
        "basis_raw": "(0,1,1),(1,0,0),(0,0,-1)",
        "origin": [0.0, 0.0, 0.0],
        "origin_raw": "(0,0,0)",
    }
    sg = IsoWrapper._subgroup_from_row(
        row, index=0, k_point="GM", irrep_label="GM5+",
        k_coordinates=["0", "0", "0"], parent_sg=139,
    )
    assert sg.basis_raw == "(0,1,1),(1,0,0),(0,0,-1)"
    assert "basis={(0,1,1),(1,0,0),(0,0,-1)}" in sg.opd_line()

def test_subgroupinfo_opd_line_no_maximal_asterisk():
    sg = _formats_SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        subgroup_index=1,
        size=1,
        is_maximal=True,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        irrep_label="GM1+",
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    line = sg.opd_line()
    assert "*" not in line
    assert "k-active= (0,0,0)" in line

def test_official_html_radio_lines_parseable():
    if not _OFFICIAL_HTML.is_file():
        return
    html = _OFFICIAL_HTML.read_text(encoding="utf-8", errors="replace")
    lines = re.findall(r'name="orderparam"[^>]*>([^<]+)<br>', html)
    assert len(lines) >= 50
    assert lines[0].startswith("GM1+     P1   (a) 139 I4/mmm")
    assert "basis={(1,0,0),(0,1,0),(0,0,1)}" in lines[0]
    assert "k-active=" in lines[0]
    formatted = format_opd_line(
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=1,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert formatted == lines[0].strip()
    x_p1 = next(ln for ln in lines if ln.startswith("X1+") and "(a;a)" in ln)
    assert "k-active= (1/2,1/2,0),(1/2,1/2,1)" in x_p1
    gm4 = next(ln for ln in lines if ln.startswith("GM4+"))
    assert "69 Fmmm" in gm4
    assert "basis={(1,1,0),(-1,1,0),(0,0,1)}" in gm4
    gm5_bases = [
        ln for ln in lines if ln.startswith("GM5+")
    ]
    assert any("{(0,1,1),(-1,0,0),(0,-1,0)}" in ln for ln in gm5_bases)
    assert any("{(1,-1,0),(1,1,0),(-1/2,1/2,1/2)}" in ln for ln in gm5_bases)
    assert any("{(1,0,0),(0,1,0),(-1/2,-1/2,1/2)}" in ln for ln in gm5_bases)
    m1 = next(ln for ln in lines if ln.startswith("M1+"))
    assert "k-active= (1,1,1)" in m1

def test_official_method1_fields_match_radio_tokens():
    from backend.utils.opd_format import official_method1_fields

    fields = official_method1_fields(
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=1,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert fields["irrep"] == "GM1+"
    assert fields["opd"] == "P1"
    assert fields["dir"] == "(a)"
    assert fields["sg"] == "139 I4/mmm"
    assert fields["basis"] == "{(1,0,0),(0,1,0),(0,0,1)}"
    assert fields["origin"] == "(0,0,0)"
    assert fields["s"] == 1
    assert fields["i"] == 1
    assert fields["k_active"] == "(0,0,0)"
    sg = _formats_SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        subgroup_index=1,
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        irrep_label="GM1+",
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert sg.official_fields()["k_active"] == "(0,0,0)"
    assert "k-active= (0,0,0)" in sg.opd_line()

# Original coverage: test_isotropy_cache.py.

def test_safe_name_strips_path():
    assert _safe_name(r"a/b/i0123400.iso") == "i0123400.iso"
    assert _safe_name(r"C:\tmp\i0123400.iso") == "i0123400.iso"

def test_list_isotropy_cache_smoke():
    try:
        wrapper = _isotropy_cache_IsoWrapper()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"WSL/iso unavailable: {exc}")
    entries = list_isotropy_cache(wrapper)
    assert isinstance(entries, list)
    for e in entries:
        assert re.match(r"^i\d+\.iso$", e.name, re.I)
        assert e.size >= 0
        d = e.to_dict()
        assert d["name"] == e.name

def test_delete_rejects_non_iso_names():
    try:
        wrapper = _isotropy_cache_IsoWrapper()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"WSL/iso unavailable: {exc}")
    result = delete_isotropy_cache(wrapper, ["../etc/passwd", "data_isotropy.txt"])
    assert result["deleted"] == []
    assert len(result["skipped"]) == 2

# Original coverage: test_microscopic_column_provenance.py.

def _one_column_microscopic_table(irrep: str, vector: str) -> str:
    return (
        "********Irrep (ML) Wyckoff Irrep Point Projected Vectors\n"
        f"{irrep} a A (0,0,0) ({vector})\n"
    )

def test_parametric_microscopic_tables_preserve_repeated_irrep_boundaries():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._direction_transport_for_invariant = lambda _parent, _invariant: (
        SimpleNamespace(symbol="P1", domain_number=1)
    )
    wrapper._run_session = lambda _commands: (
        _one_column_microscopic_table("LD5", "1,0,0")
        + _one_column_microscopic_table("LD5", "0,1,0")
    )
    subgroup = _microscopic_column_provenance_subgroup(irrep="LD5")
    subgroup.k_coordinates = ["0", "0", "1/6"]
    invariants = [
        InvariantDirection(
            "LD5", "(a,0,0,0)", 1, "P1", 1,
            k_parameters=("1/12",), k_coordinates=("0", "0", "1/6"),
        ),
        InvariantDirection(
            "LD3", "(a,0)", 1, "P1", 1,
            k_parameters=("1/6",), k_coordinates=("0", "0", "1/3"),
        ),
        InvariantDirection(
            "LD5", "(0,a,0,0)", 1, "P1", 1,
            k_parameters=("1/4",), k_coordinates=("0", "0", "1/2"),
        ),
    ]

    queries, commands = _microscopic_column_provenance_IsoWrapper._parametric_microscopic_query_results(
        wrapper, 139, subgroup, invariants, ["a"]
    )

    assert commands.count("VALUE DIRECTION VECTOR,A,0,0,0") == 1
    assert commands.count("VALUE DIRECTION VECTOR,A,0") == 1
    assert commands.count("VALUE DIRECTION VECTOR,0,A,0,0") == 1
    assert commands.count("VALUE DOMAIN 1") == 3
    assert [query.transport_direction for query in queries] == [
        "VECTOR,A,0,0,0",
        "VECTOR,A,0",
        "VECTOR,0,A,0,0",
    ]
    assert [query.transport_domain_number for query in queries] == [1, 1, 1]
    assert [len(query.blocks) for query in queries] == [1, 0, 1]
    assert queries[0].blocks[0].direction_symbol == "VECTOR,A,0,0,0"
    assert queries[2].blocks[0].direction_symbol == "VECTOR,0,A,0,0"
    assert queries[0].blocks[0].rows[0].displacements[0] == (
        Fraction(1), Fraction(0), Fraction(0)
    )
    assert queries[2].blocks[0].rows[0].displacements[0] == (
        Fraction(0), Fraction(1), Fraction(0)
    )

def test_parametric_microscopic_tables_reject_ambiguous_empty_repeated_irrep():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._direction_transport_for_invariant = lambda _parent, _invariant: (
        SimpleNamespace(symbol="P1", domain_number=1)
    )
    wrapper._run_session = lambda _commands: _one_column_microscopic_table(
        "LD5", "1,0,0"
    )
    subgroup = _microscopic_column_provenance_subgroup(irrep="LD5")
    subgroup.k_coordinates = ["0", "0", "1/6"]
    invariants = [
        InvariantDirection(
            "LD5", "(a,0,0,0)", 1, "P1", 1,
            k_parameters=("1/12",), k_coordinates=("0", "0", "1/6"),
        ),
        InvariantDirection(
            "LD5", "(0,a,0,0)", 1, "P1", 1,
            k_parameters=("1/4",), k_coordinates=("0", "0", "1/2"),
        ),
    ]

    with pytest.raises(_microscopic_column_provenance_OutputParseError, match="无法唯一归属"):
        _microscopic_column_provenance_IsoWrapper._parametric_microscopic_query_results(
            wrapper, 139, subgroup, invariants, ["a"]
        )

def test_exact_fixed_space_wrapper_has_no_synthetic_primary_route():
    modes = _microscopic_column_provenance_IsoWrapper.calc_embedding_microscopic_distortion_modes(
        _exact_fixed_space_wrapper(),
        1,
        _exact_fixed_space_subgroup(),
        ["a"],
    )
    assert len(modes) == 3
    assert all(
        mode.microscopic_provenance.source_subgroup_context_kind
        == "exact_fixed_space"
        and mode.microscopic_provenance.source_subgroup_irrep_label is None
        and mode.microscopic_provenance.source_subgroup_opd_symbol is None
        and mode.microscopic_provenance.source_subgroup_primary_direction_selector
        is None
        for mode in modes
    )

@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("k_point_label", "X"),
        ("k_coordinates", ["1/2", "0", "0"]),
        ("k_parameters", ["1/2"]),
        ("k_active_raw", "(1/2,0,0)"),
        ("irrep_label", "X1"),
        ("opd_symbol", "P1"),
        ("opd_dir_raw", "(a)"),
        ("opd_vector", [1]),
    ),
)
def test_exact_fixed_space_wrapper_rejects_any_synthetic_primary_route(
    field,
    value,
):
    subgroup = _exact_fixed_space_subgroup()
    setattr(subgroup, field, value)

    with pytest.raises(ValueError, match="cannot claim one primary"):
        _microscopic_column_provenance_IsoWrapper.calc_embedding_microscopic_distortion_modes(
            _exact_fixed_space_wrapper(),
            1,
            subgroup,
            ["a"],
        )

def test_unverified_opd_alias_does_not_enter_scientific_source_digest():
    provenance = _provenance(0)
    bogus_alias = replace(provenance, source_subgroup_opd_symbol="BOGUS")

    assert bogus_alias.exact_source_token == provenance.exact_source_token
    assert bogus_alias.stable_token == provenance.stable_token

def test_direction_resolved_microscopic_query_binds_exact_route_and_opd():
    subgroup = _microscopic_column_provenance_subgroup(irrep="X1+")
    irreps, letters, commands = _microscopic_column_provenance_IsoWrapper._microscopic_query_commands(
        62,
        subgroup,
        ["X1+", "GM1+"],
        ["a"],
        direction_symbols={
            "X1+": "VECTOR,0,A",
            "GM1+": "VECTOR,A",
        },
    )

    assert irreps == ("X1+", "GM1+")
    assert letters == ("a",)
    text = "\n".join(commands)
    assert "VALUE KPOINT X\nVALUE IRREP X1+\nVALUE DIRECTION VECTOR,0,A" in text
    assert "VALUE KPOINT GM\nVALUE IRREP GM1+\nVALUE DIRECTION VECTOR,A" in text
    assert text.count("VALUE SUBGROUP 1") == 2
    assert text.count("DISPLAY DISTORTION") == 2

def test_microscopic_vector_query_rejects_syntax_error_with_parseable_full_block():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: """
        ******Syntax error: VECTOR,A/2,0
        Irrep (ML) Wyckoff Irrep Point Projected Vectors
        X1+ a A (0,0,0) (1,0,0)
        *
    """

    with pytest.raises(_microscopic_column_provenance_OutputParseError, match="VECTOR"):
        wrapper.get_microscopic_vector_blocks(
            62,
            _microscopic_column_provenance_subgroup(irrep="X1+"),
            ["X1+"],
            ["a"],
            direction_symbols={"X1+": "VECTOR,A"},
        )

def test_microscopic_vector_query_rejects_unconditioned_column_width():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: """
        Irrep (ML) Wyckoff Irrep Point Projected Vectors
        X1+ a A (0,0,0) (1,0,0), (0,1,0)
        *
    """

    with pytest.raises(_microscopic_column_provenance_OutputParseError, match="错误列数"):
        wrapper.get_microscopic_vector_blocks(
            62,
            _microscopic_column_provenance_subgroup(irrep="X1+"),
            ["X1+"],
            ["a"],
            direction_symbols={"X1+": "VECTOR,A"},
        )

@pytest.mark.parametrize(
    "unexpected_row",
    (
        "Y1+ a A (0,0,0) (1,0,0)",
        "X1+ b A (0,0,0) (1,0,0)",
    ),
)
def test_microscopic_vector_query_rejects_unrequested_output_blocks(
    unexpected_row,
):
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: f"""
        Irrep (ML) Wyckoff Irrep Point Projected Vectors
        {unexpected_row}
        *
    """

    with pytest.raises(_microscopic_column_provenance_OutputParseError, match="未请求"):
        wrapper.get_microscopic_vector_blocks(
            62,
            _microscopic_column_provenance_subgroup(irrep="X1+"),
            ["X1+"],
            ["a"],
            direction_symbols={"X1+": "VECTOR,A"},
        )

def test_exact_invariant_vector_resolves_one_cached_direction_symbol():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    calls = []

    def list_subgroups(parent_sg, k_point, irrep_label):
        calls.append((parent_sg, k_point, irrep_label))
        return [
            SimpleNamespace(
                opd_dir_raw="(a, 0)", opd_symbol="P1",
                space_group_number=2, size=2,
            ),
            SimpleNamespace(
                opd_dir_raw="(a; a)", opd_symbol="P3",
                space_group_number=1, size=1,
            ),
        ]

    wrapper.list_subgroups = list_subgroups
    invariant = InvariantDirection("X1+", "(a,a)", 1, "P1", 1)

    assert wrapper._direction_symbol_for_invariant(62, invariant) == "P3"
    assert wrapper._direction_symbol_for_invariant(62, invariant) == "P3"
    assert calls == [(62, "X", "X1+")]

    wrapper._direction_transport_cache = {}
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(a,a)", opd_symbol="P3",
            space_group_number=1, size=1,
        ),
        SimpleNamespace(
            opd_dir_raw="(b,b)", opd_symbol="C1",
            space_group_number=1, size=1,
        ),
    ]
    with pytest.raises(ValueError, match="one ISO OPD/domain"):
        wrapper._direction_symbol_for_invariant(62, invariant)

def test_exact_invariant_vector_resolves_noncanonical_named_domain():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    row = SimpleNamespace(
        opd_dir_raw="(a,0;0,0)",
        opd_symbol="P1",
        space_group_number=47,
        size=3,
    )
    wrapper.list_subgroups = lambda *_args, **_kwargs: [row]
    wrapper.get_domains = lambda *_args, **_kwargs: [
        SimpleNamespace(domain_number=2, opd_dir_raw="(0,0;a,0)"),
        SimpleNamespace(
            domain_number=4,
            opd_dir_raw="(0,0;-0.500a,0.866a)",
        ),
    ]
    invariant = InvariantDirection(
        "DT1",
        "(0,0,a,-1.732a)",
        47,
        "Pmmm",
        3,
        k_parameters=("1/3",),
        k_coordinates=("0", "1/3", "0"),
    )

    transport = wrapper._direction_transport_for_invariant(123, invariant)

    assert transport.symbol == "P1"
    assert transport.domain_number == 4
    assert transport.direction_raw == "(0,0;-0.500a,0.866a)"

def test_microscopic_transport_uses_only_legacy_safe_exact_vectors_directly():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._direction_transport_for_invariant = lambda *_args: (_ for _ in ()).throw(
        AssertionError("safe exact VECTOR must not require isotropy-subgroup lookup")
    )
    invariant = InvariantDirection(
        "SM1",
        "(a,0,a,0)",
        123,
        "P4/mmm",
        1,
        k_parameters=("1/6",),
        k_coordinates=("1/6", "1/6", "0"),
    )

    transport = wrapper._microscopic_direction_transport(123, invariant)

    assert transport.symbol == "VECTOR,A,0,A,0"
    assert transport.domain_number == 1
    assert transport.direction_raw == "(a,0,a,0)"

def test_microscopic_transport_falls_back_for_scaled_exact_vector():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    fallback = SimpleNamespace(
        symbol="P2",
        domain_number=1,
        direction_raw="(a,0.577a;0,0)",
    )
    wrapper._direction_transport_for_invariant = lambda *_args: fallback
    invariant = InvariantDirection(
        "Y1",
        "(a,0.577a;0,0)",
        51,
        "Pmma",
        3,
        k_parameters=("1/3",),
        k_coordinates=("0", "1/3", "0"),
    )

    transport = wrapper._microscopic_direction_transport(123, invariant)

    assert _microscopic_column_provenance_IsoWrapper._direction_vector_selector(invariant.direction_raw) == (
        "VECTOR,1000A,577A,0,0"
    )
    assert transport is fallback

def test_microscopic_transport_reuses_exact_selected_primary_opd():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._direction_transport_for_invariant = lambda *_args: (_ for _ in ()).throw(
        AssertionError("the selected primary OPD is already exact evidence")
    )
    subgroup = _microscopic_column_provenance_SubgroupInfo(
        index=3,
        space_group_number=127,
        opd_symbol="P4",
        opd_dir_raw="(a,0.577a;0,1.155a)",
        irrep_label="Y1",
        parent_sg=123,
        size=36,
        k_parameters=["1/3"],
    )
    invariant = InvariantDirection(
        "Y1",
        "(a,0.577a,0,1.155a)",
        127,
        "P4/mbm",
        36,
        k_parameters=("1/3",),
        k_coordinates=("1/3", "1/2", "0"),
    )

    transport = wrapper._microscopic_direction_transport(
        123, invariant, subgroup
    )

    assert transport.symbol == "P4"
    assert transport.domain_number == 1
    assert transport.direction_raw == "(a,0.577a;0,1.155a)"

def test_parametric_microscopic_query_selects_resolved_domain():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper._direction_transport_for_invariant = lambda *_args: SimpleNamespace(
        symbol="P1",
        domain_number=4,
    )
    wrapper._run_session = lambda commands: (
        "********Irrep (ML) Domain Wyckoff Irrep Point Projected Vectors\n"
        "DT1 4 a A (0,0,0) (0,1,0)\n"
    )
    subgroup = _microscopic_column_provenance_subgroup(irrep="Y1")
    invariant = InvariantDirection(
        "DT1",
        "(0,0,a,-1.732a)",
        47,
        "Pmmm",
        3,
        k_parameters=("1/3",),
        k_coordinates=("0", "1/3", "0"),
    )

    queries, commands = _microscopic_column_provenance_IsoWrapper._parametric_microscopic_query_results(
        wrapper,
        123,
        subgroup,
        [invariant],
        ["a"],
    )

    assert "VALUE DIRECTION P1" in commands
    assert "VALUE DOMAIN 4" in commands
    assert queries[0].transport_direction == "P1"
    assert queries[0].transport_domain_number == 4
    assert queries[0].blocks[0].column_count == 1
    assert queries[0].blocks[0].direction_symbol == "VECTOR,0,0,250A,-433A"

def test_domain_parser_retains_multiarm_direction_and_all_components():
    rows = parse_domain_table(
        """
        4 (C2b|-1,0,0) 47 Pmmm 6 P1 (0,0;-0.500a,0.866a) (0,1,0),(0,0,-1),(-3,0,0) (-1,0,0)
        """
    )

    assert rows[0]["opd_dir_raw"] == "(0,0;-0.500a,0.866a)"
    assert rows[0]["opd_vector"] == [0.0, 0.0, -0.5, 0.866]

def test_direction_comparison_accepts_identical_decimal_iso_tokens_only():
    assert _microscopic_column_provenance_IsoWrapper._direction_subspaces_same_component_order(
        "(a,0.577a,0,0)",
        "(a,0.577a;0,0)",
    )
    assert not _microscopic_column_provenance_IsoWrapper._direction_subspaces_same_component_order(
        "(a,0.577a,0,0)",
        "(a,0.578a;0,0)",
    )
    assert not _microscopic_column_provenance_IsoWrapper._direction_subspaces_same_component_order(
        "(0,a)",
        "(a;0)",
    )

def test_exact_direction_vector_serializer_preserves_component_order():
    assert _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(0,a)") == "VECTOR,0,A"
    assert (
        _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(a;b;b;-a)")
        == "VECTOR,A,B,B,-A"
    )
    assert (
        _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(2*a+b;4*a-3*b)")
        == "VECTOR,A+B,2A-3B"
    )
    assert _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(3*a,a)") == "VECTOR,3A,A"
    assert _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(a/2,a/3)") == "VECTOR,3A,2A"
    assert _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(0.5*a,0)") == "VECTOR,A,0"
    assert (
        _microscopic_column_provenance_IsoWrapper._direction_vector_selector("(a,0.577a;0,0)")
        == "VECTOR,1000A,577A,0,0"
    )

@pytest.mark.parametrize(
    "value",
    [
        "(0,0)",
        "(a*a,0)",
        "(sqrt(2)*a,0)",
        "(a+1,0)",
        "(a+b,2*a+2*b)",
        "(a+A,0)",
        "a,0",
    ],
)
def test_exact_direction_vector_serializer_rejects_unrepresentable_forms(value):
    with pytest.raises(ValueError):
        _microscopic_column_provenance_IsoWrapper._direction_vector_selector(value)

def test_primary_direction_symbol_rejects_permuted_star_arm():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.get_domains = lambda *_args, **_kwargs: []
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(a;a)", opd_symbol="P1",
            space_group_number=127, size=4,
        ),
        SimpleNamespace(
            opd_dir_raw="(2*b;0)", opd_symbol="P3",
            space_group_number=64, size=2,
        ),
        SimpleNamespace(
            opd_dir_raw="(a;b)", opd_symbol="C1",
            space_group_number=55, size=4,
        ),
    ]

    invariant = InvariantDirection("X2+", "(0,a)", 64, "Cmce", 2)

    with pytest.raises(ValueError, match=r"transports=\[\]"):
        wrapper._direction_symbol_for_invariant(139, invariant)

def test_direction_symbol_rejects_wrong_isotropy_identity_and_ambiguity():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.get_domains = lambda *_args, **_kwargs: []
    invariant = InvariantDirection("X2+", "(0,a)", 64, "Cmce", 2)
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(a;0)", opd_symbol="P3",
            space_group_number=65, size=2,
        )
    ]
    with pytest.raises(ValueError, match=r"transports=\[\]"):
        wrapper._direction_symbol_for_invariant(139, invariant)

    wrapper._direction_transport_cache = {}
    wrapper.list_subgroups = lambda *_args: [
        SimpleNamespace(
            opd_dir_raw="(0;b)", opd_symbol="P3",
            space_group_number=64, size=2,
        ),
        SimpleNamespace(
            opd_dir_raw="(0;c)", opd_symbol="C9",
            space_group_number=64, size=2,
        ),
    ]
    with pytest.raises(ValueError, match="one ISO OPD/domain"):
        wrapper._direction_symbol_for_invariant(139, invariant)

def test_parser_modes_retain_exact_column_and_unresolved_source_facts():
    block = replace(
        _block(irrep="X1+", source_order=7),
        direction_symbol="VECTOR,A,B",
    )
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.get_microscopic_vector_blocks = lambda *_args, **_kwargs: [block]
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1),
        InvariantDirection("X1+", "(a,a)", 1, "P1", 1),
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P1"
    wrapper.list_k_points = lambda _parent: [
        KPointInfo("X", ["0", "1/2", "0"], [], True)
    ]

    subgroup = _microscopic_column_provenance_subgroup()
    subgroup.basis_vectors = [[1, 0, 0], [0, 2, 0], [0, 0, 1]]
    subgroup.basis_raw = "(1,0,0),(0,2,0),(0,0,1)"
    modes = _microscopic_column_provenance_IsoWrapper.calc_microscopic_distortion_modes(
        wrapper,
        62,
        subgroup,
        ["a"],
        irrep_labels=["X1+"],
    )

    assert [mode.microscopic_provenance.source_block_order for mode in modes] == [7, 7]
    assert [mode.microscopic_provenance.source_column_index for mode in modes] == [0, 1]
    assert all(mode.microscopic_provenance.query_order == 0 for mode in modes)
    assert all(mode.microscopic_provenance.source_row_count == 2 for mode in modes)
    assert modes[0].microscopic_provenance.exact_row_digest == modes[1].microscopic_provenance.exact_row_digest
    assert modes[0].microscopic_provenance.exact_vector_digest != modes[1].microscopic_provenance.exact_vector_digest
    assert all(
        mode.microscopic_provenance.unresolved_fields
        == ("source_frame_id", "source_atom_order_id")
        for mode in modes
    )
    assert [mode.mode_identity.k_coordinates for mode in modes] == [
        ("0", "1/2", "0"),
        ("0", "1/2", "0"),
    ]
    assert all(mode.opd_dir_raw == "(a,a)" for mode in modes)
    assert all(
        mode.microscopic_provenance.source_direction_symbol == "VECTOR,A,B"
        for mode in modes
    )

@pytest.mark.parametrize(
    "listed",
    [
        [KPointInfo("X", ["a", "0", "0"], ["a"], False)],
        [
            KPointInfo("X", ["0", "1/2", "0"], [], True),
            KPointInfo("X", ["1/2", "0", "0"], [], True),
        ],
    ],
)
def test_parametric_or_ambiguous_listed_secondary_k_remains_unresolved(listed):
    block = _block(irrep="X1+")
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.get_microscopic_vector_blocks = lambda *_args, **_kwargs: [block]
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1),
        InvariantDirection("X1+", "(a,a)", 1, "P1", 1),
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P1"
    wrapper.list_k_points = lambda _parent: listed
    modes = _microscopic_column_provenance_IsoWrapper.calc_microscopic_distortion_modes(
        wrapper, 62, _microscopic_column_provenance_subgroup(), ["a"], irrep_labels=["X1+"]
    )
    assert all(mode.mode_identity.status == "unresolved" for mode in modes)
    assert all(
        mode.mode_identity.reason == "exact_special_k_identity_unresolved"
        for mode in modes
    )
    assert all(mode.mode_identity.k_coordinates == () for mode in modes)

def test_microscopic_calc_requires_exact_direction_for_explicit_irreps():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1)
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P1"

    with pytest.raises(ValueError, match="requested irrep"):
        _microscopic_column_provenance_IsoWrapper.calc_microscopic_distortion_modes(
            wrapper,
            62,
            _microscopic_column_provenance_subgroup(),
            ["a"],
            irrep_labels=["X1+"],
        )

def test_microscopic_calc_rejects_primary_opd_direction_mismatch():
    wrapper = object.__new__(_microscopic_column_provenance_IsoWrapper)
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM1+", "(a)", 1, "P1", 1)
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "P2"

    with pytest.raises(ValueError, match="selected subgroup OPD"):
        _microscopic_column_provenance_IsoWrapper.calc_microscopic_distortion_modes(
            wrapper,
            62,
            _microscopic_column_provenance_subgroup(),
            ["a"],
            irrep_labels=["GM1+"],
        )

# Original coverage: test_native_wrapper.py.

def _native_wrapper(tmp_path: Path, monkeypatch) -> BaseWrapper:
    temporary_directory = tempfile.TemporaryDirectory
    monkeypatch.setattr(
        tempfile, "TemporaryDirectory",
        lambda *, prefix, dir: temporary_directory(prefix=prefix, dir=tmp_path),
    )
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._mode = "native"
    wrapper._init_native_environment()
    return wrapper

def _remove_source_execute_permission(source: Path, monkeypatch) -> None:
    source.chmod(0o600)
    actual_access = os.access
    # Windows does not implement POSIX execute bits; emulate only the source
    # permission result while exercising real copying and file metadata.
    monkeypatch.setattr(
        os, "access",
        lambda path, mode: False if Path(path) == source and mode == os.X_OK else actual_access(path, mode),
    )

def test_native_nonexecutable_binary_is_privately_copied_without_source_changes(
    tmp_path, monkeypatch,
) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"binary fixture")
    _remove_source_execute_permission(source, monkeypatch)
    before = (source.read_bytes(), source.stat().st_mode)
    first = _native_wrapper(tmp_path, monkeypatch)
    second = _native_wrapper(tmp_path, monkeypatch)

    executable = Path(first._wsl_bin_path(source))

    assert executable != source
    assert executable.parent == Path(first._stage_dir)
    assert first._stage_dir != second._stage_dir
    assert executable.read_bytes() == source.read_bytes()
    assert (source.read_bytes(), source.stat().st_mode) == before
    assert first._wsl_bin_path(source) == str(executable)
    if os.name == "posix":
        assert executable.stat().st_mode & 0o777 == 0o700
        assert executable.parent.stat().st_mode & 0o777 == 0o700

def test_native_executable_binary_runs_directly(tmp_path, monkeypatch) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"binary fixture")
    monkeypatch.setattr(os, "access", lambda _path, _mode: True)
    wrapper = _native_wrapper(tmp_path, monkeypatch)

    assert wrapper._wsl_bin_path(source) == str(source.resolve())
    assert not list(Path(wrapper._stage_dir).iterdir())

def test_wsl_binary_path_keeps_existing_conversion(monkeypatch) -> None:
    wrapper = object.__new__(BaseWrapper)
    wrapper._mode = "wsl"
    monkeypatch.setattr(wrapper, "_win_to_wsl", lambda _path: "/mnt/c/isobyu/iso")

    assert wrapper._wsl_bin_path(Path("iso")) == "/mnt/c/isobyu/iso"
    assert wrapper._wsl_bin_path(PurePosixPath("/home/user/iso")) == "/home/user/iso"

def test_native_staging_failure_is_reported(monkeypatch) -> None:
    def deny_staging(**_kwargs):
        raise PermissionError("private directory denied")

    monkeypatch.setattr(tempfile, "TemporaryDirectory", deny_staging)
    wrapper = object.__new__(BaseWrapper)

    with pytest.raises(_native_wrapper_WrapperRunError, match="Cannot create private staging directory"):
        wrapper._init_native_environment()

def test_native_binary_copy_is_invalidated_when_source_is_replaced(
    tmp_path, monkeypatch,
) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"old")
    _remove_source_execute_permission(source, monkeypatch)
    wrapper = _native_wrapper(tmp_path, monkeypatch)
    previous = Path(wrapper._wsl_bin_path(source))

    source.write_bytes(b"new binary")
    current = Path(wrapper._wsl_bin_path(source))

    assert current != previous
    assert current.read_bytes() == b"new binary"
    assert source.read_bytes() == b"new binary"

def test_native_binary_copy_failure_has_actionable_error(tmp_path, monkeypatch) -> None:
    source = tmp_path / "iso"
    source.write_bytes(b"binary fixture")
    _remove_source_execute_permission(source, monkeypatch)
    wrapper = _native_wrapper(tmp_path, monkeypatch)
    def deny_chmod(_path, _mode):
        raise PermissionError("blocked")

    monkeypatch.setattr(Path, "chmod", deny_chmod)

    with pytest.raises(_native_wrapper_WrapperRunError, match=r"Cannot prepare native executable.*blocked"):
        wrapper._wsl_bin_path(source)

    assert source.read_bytes() == b"binary fixture"
    assert not list(Path(wrapper._stage_dir).iterdir())

def test_generated_isodata_uses_private_disposable_copy(tmp_path) -> None:
    source = tmp_path / "reference-isodata"
    source.mkdir()
    reference = source / "data_isotropy.txt"
    reference.write_text("reference\n", encoding="utf-8")
    runtime = tmp_path / "runtime"
    runtime.mkdir()
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._mode = "native"
    wrapper._runtime_isodata = None
    wrapper.cfg = SimpleNamespace(
        temp_dir=runtime,
        _cfg={"isobyu": {"data_dir": str(source)}},
        resolve_path=Path,
    )

    wrapper._ensure_writable_isodata()

    private = Path(wrapper._isodata_path())
    assert private != source
    assert (private / reference.name).read_text(encoding="utf-8") == "reference\n"
    (private / reference.name).write_text("generated\n", encoding="utf-8")
    assert reference.read_text(encoding="utf-8") == "reference\n"

def test_runtime_clone_survives_source_close_and_does_not_share_mutations(tmp_path, monkeypatch):
    source_data = tmp_path / "reference"
    source_data.mkdir()
    (source_data / "database").write_text("reference", encoding="utf-8")
    cfg = SimpleNamespace(
        temp_dir=tmp_path, _cfg={"isobyu": {"data_dir": str(source_data)}},
        resolve_path=Path,
    )
    source = _native_wrapper(tmp_path, monkeypatch)
    target = _native_wrapper(tmp_path, monkeypatch)
    source.cfg = target.cfg = cfg
    source._ensure_writable_isodata()
    source_private = Path(source._runtime_isodata.name)
    (source_private / "database").write_text("generated", encoding="utf-8")
    (Path(source._stage_dir) / "i123.iso").write_text("cache", encoding="utf-8")
    (Path(source._stage_dir) / "input.in").write_text("not a cache", encoding="utf-8")

    target.clone_runtime_from(source)
    target_private = Path(target._runtime_isodata.name)
    assert target_private != source_private
    assert (target_private / "database").read_text(encoding="utf-8") == "generated"
    assert (Path(target._stage_dir) / "i123.iso").read_text(encoding="utf-8") == "cache"
    assert not (Path(target._stage_dir) / "input.in").exists()
    (target_private / "database").write_text("changed", encoding="utf-8")
    assert (source_private / "database").read_text(encoding="utf-8") == "generated"
    source.close()
    assert not source_private.exists()
    assert (target_private / "database").read_text(encoding="utf-8") == "changed"
    assert (Path(target._stage_dir) / "i123.iso").is_file()
    target.close()
    assert not target_private.exists()
    assert (source_data / "database").read_text(encoding="utf-8") == "reference"

def test_runtime_clone_without_generated_data_does_not_copy_source_database(tmp_path, monkeypatch):
    source = _native_wrapper(tmp_path, monkeypatch)
    target = _native_wrapper(tmp_path, monkeypatch)
    target.clone_runtime_from(source)
    assert target._runtime_isodata is None
    assert list(Path(target._stage_dir).iterdir()) == []
    source.close()
    target.close()

def test_failed_runtime_clone_closes_destination_and_keeps_source(tmp_path, monkeypatch):
    source = _native_wrapper(tmp_path, monkeypatch)
    target = _native_wrapper(tmp_path, monkeypatch)
    cache = Path(source._stage_dir) / "i123.iso"
    cache.write_text("source cache", encoding="utf-8")
    target_dir = Path(target._stage_dir)

    def fail_copy(*_args, **_kwargs):
        raise OSError("copy interrupted")

    monkeypatch.setattr("backend.wrappers.base_wrapper.shutil.copy2", fail_copy)
    with pytest.raises(OSError, match="copy interrupted"):
        target.clone_runtime_from(source)
    assert target._closed
    assert not target_dir.exists()
    assert not source._closed
    assert cache.read_text(encoding="utf-8") == "source cache"
    source.close()

@pytest.mark.skipif(os.name != "posix", reason="requires native POSIX shell execution")
def test_native_run_executes_private_copy_with_isodata_and_stdin(tmp_path) -> None:
    source = tmp_path / "iso"
    source.write_text('#!/bin/sh\nprintf "%s\\n" "$ISODATA"\ncat\n', encoding="utf-8")
    source.chmod(0o600)
    data = tmp_path / "data"
    data.mkdir()
    temporary_inputs = tmp_path / "inputs"
    temporary_inputs.mkdir()
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._mode = "native"
    wrapper.cfg = SimpleNamespace(
        timeout=5, temp_dir=temporary_inputs,
        _cfg={"isobyu": {"data_dir": str(data)}},
        resolve_path=Path,
    )
    wrapper._init_native_environment()
    before = source.stat().st_mode

    assert wrapper.run_stdin(source, "input text\n").splitlines() == [str(data) + "/", "input text"]
    assert source.stat().st_mode == before

# Original coverage: test_search_methods.py.

def test_iso_bush_queries_every_wyckoff_orbit_before_merging(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wrapper = object.__new__(_search_methods_IsoWrapper)
    seen: list[str] = []

    def fake_session(commands: str, timeout=None) -> str:
        del timeout
        seen.append(commands)
        return """
Irrep (ML) Dir(dom) Wyckoff Point            Displacement
GM5+       C1(1)    d       (0,1/2,1/4)      (1,0,0), (0,1,0)
Irrep (ML) Dir(dom) Wyckoff Point            Displacement
GM5+       C1(1)    e       (0,0,z)          (1,0,0)
"""

    monkeypatch.setattr(wrapper, "_run_session", fake_session)
    subgroup = _search_methods_SubgroupInfo(
        index=0,
        space_group_number=2,
        k_point_label="GM",
        irrep_label="GM5+",
        opd_symbol="C1",
    )

    modes = wrapper.calc_distortion_modes(139, subgroup, ["d", "e", "d"])

    assert len(seen) == 1
    assert seen[0].count("VALUE WYCKOFF d\nDISPLAY BUSH") == 1
    assert seen[0].count("VALUE WYCKOFF e\nDISPLAY BUSH") == 1
    assert len(modes) == 3
    assert [mode.wyckoff_site for mode in modes] == ["d", "d", "e"]
    assert [mode.opd_component for mode in modes] == ["a", "b", "a"]
    assert [mode.dimension for mode in modes] == [2, 2, 1]
    assert [mode.amplitude_key for mode in modes] == [
        "GM5+__d__C1(1)__a",
        "GM5+__d__C1(1)__b",
        "GM5+__e__C1(1)__a",
    ]
    assert all(
        len(row.displacements) == 1
        for mode in modes
        for row in mode.bush_modes
    )
    assert modes[0].bush_modes[0].displacements == [[1.0, 0.0, 0.0]]
    assert modes[1].bush_modes[0].displacements == [[0.0, 1.0, 0.0]]

# Original coverage: test_site_irreps.py.

def test_iso_microscopic_block_order_defines_copy_and_component_identity():
    row = MicroscopicVectorRow(
        point_raw=("0", "0", "z"),
        displacements=(
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
        ),
    )
    blocks = [
        MicroscopicVectorBlock(
            "GM5+", "g", "A1", (row,), source_order=index,
            direction_symbol="VECTOR,A,B",
        )
        for index in range(2)
    ]
    wrapper = object.__new__(_site_irreps_IsoWrapper)
    wrapper.get_microscopic_vector_blocks = lambda *_args, **_kwargs: blocks
    wrapper.list_invariant_directions = lambda *_args: [
        InvariantDirection("GM5+", "(a;b)", 2, "P-1", 8)
    ]
    wrapper._direction_symbol_for_invariant = lambda *_args: "4D1"
    subgroup = _site_irreps_SubgroupInfo(
        index=0,
        space_group_number=2,
        opd_symbol="4D1",
        irrep_label="GM5+",
        parent_sg=139,
        k_coordinates=["0", "0", "0"],
        basis_vectors=np.eye(3).tolist(),
        origin=[0.0, 0.0, 0.0],
    )
    modes = _site_irreps_IsoWrapper.calc_microscopic_distortion_modes(
        wrapper,
        139,
        subgroup,
        ["g"],
        irrep_labels=["GM5+"],
    )
    assert [mode.site_irrep for mode in modes] == ["A1_1", "A1_1", "A1_2", "A1_2"]
    assert [mode.opd_component for mode in modes] == ["a", "b", "a", "b"]
    assert [mode.mode_identity.copy_index for mode in modes] == [1, 1, 2, 2]
    assert [mode.mode_identity.component_index for mode in modes] == [0, 1, 0, 1]
    assert all(mode.mode_identity.status == "partial" for mode in modes)
    assert all(
        mode.mode_identity.reason == "awaiting_bush_or_smodes_subspace_validation"
        for mode in modes
    )
    assert len({mode.mode_identity.stable_token for mode in modes}) == 4

def test_macroscopic_wrapper_assigns_one_based_copy_indices_and_fails_closed():
    output = """
    Irrep (ML) Basis Functions
    GM1+       xx+yy
               zz
    GM2+       xx-yy
    GM4+       xy
    GM5+       yz,-xz
    *
    """
    wrapper = object.__new__(_site_irreps_IsoWrapper)
    wrapper._run_session = lambda *_args, **_kwargs: output
    basis = _site_irreps_IsoWrapper.get_macroscopic_tensor_basis(wrapper, 139)
    assert basis.status == "verified"
    assert [(block.global_irrep, block.copy_index) for block in basis.blocks] == [
        ("GM1+", 1), ("GM1+", 2), ("GM2+", None),
        ("GM4+", None), ("GM5+", None),
    ]

    wrapper._run_session = lambda *_args, **_kwargs: "no table"
    missing = _site_irreps_IsoWrapper.get_macroscopic_tensor_basis(wrapper, 139)
    assert missing.status == "unresolved"
    assert missing.blocks == ()
    assert missing.reason == "iso_macroscopic_tensor_output_empty"

    def _unavailable(*_args, **_kwargs):
        raise _site_irreps_WrapperRunError("iso", 1, "unavailable")

    wrapper._run_session = _unavailable
    unavailable = _site_irreps_IsoWrapper.get_macroscopic_tensor_basis(wrapper, 139)
    assert unavailable.status == "unresolved"
    assert unavailable.blocks == ()
    assert unavailable.reason == "iso_macroscopic_tensor_output_unavailable"

# Original coverage: test_wrapper_runtime.py.

@pytest.mark.skipif(os.name != "nt", reason="regresses Windows inherited pipe handles")
@pytest.mark.parametrize("parent_waits", [True, False])
def test_wsl_capture_does_not_wait_for_descendant_output_handles(tmp_path, monkeypatch, parent_waits):
    """A finite-lived grandchild makes the old run(PIPE) hang measurable, not infinite."""
    ready = tmp_path / "grandchild-ready"
    finished = tmp_path / "grandchild-finished"
    grandchild = (
        "from pathlib import Path; import time; "
        f"Path({str(ready)!r}).write_text('ready'); "
        "time.sleep(5); "
        f"Path({str(finished)!r}).write_text('finished')"
    )
    parent = (
        "import subprocess, sys, time\nfrom pathlib import Path\n"
        f"subprocess.Popen([sys.executable, '-c', {grandchild!r}], "
        "stdout=sys.stdout, stderr=sys.stderr)\n"
        "deadline = time.monotonic() + 3\n"
        f"while not Path({str(ready)!r}).exists():\n"
        "    assert time.monotonic() < deadline, 'grandchild did not start'\n"
        "    time.sleep(0.01)\n"
        "print('captured stdout', flush=True)\n"
        "print('captured stderr', file=sys.stderr, flush=True)\n"
        + ("time.sleep(30)" if parent_waits else "")
    )
    real_popen = subprocess.Popen

    def launch_owned_python(args, **kwargs):
        assert args[:4] == ["wsl.exe", "-e", "sh", "-c"]
        return real_popen([sys.executable, "-c", parent], **kwargs)

    monkeypatch.setattr(module.subprocess, "Popen", launch_owned_python)
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._mode = "wsl"
    wrapper.cfg = SimpleNamespace(timeout=1)
    started = time.monotonic()
    try:
        if parent_waits:
            with pytest.raises(WrapperTimeoutError):
                wrapper._wsl("owned fixture", timeout=1)
        else:
            result = wrapper._wsl("owned fixture", timeout=1)
            assert result.returncode == 0
            assert result.stdout == "captured stdout"
            assert result.stderr == "captured stderr"
        elapsed = time.monotonic() - started
        assert ready.is_file(), "Fixture never established inherited output handles"
        assert elapsed < 3, "Capture waited for a descendant after the parent ended"
        assert not finished.exists(), "Descendant exited before the bound was checked"
    finally:
        # Only this test's finite-lived descendant is awaited, outside the
        # measured wrapper call; no process lookup or unrelated WSL kill occurs.
        deadline = started + 10
        while ready.exists() and not finished.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        assert not ready.exists() or finished.exists(), "Fixture descendant failed to exit"

@pytest.mark.parametrize("cleanup_failure", ["kill", "wait"])
def test_wsl_timeout_cleanup_is_bounded_and_preserves_primary_error(monkeypatch, cleanup_failure):
    calls = []
    streams = []
    primary = subprocess.TimeoutExpired("owned fixture", 0.01)
    cleanup = OSError("kill denied") if cleanup_failure == "kill" else subprocess.TimeoutExpired("reap", 1)

    class Process:
        def __init__(self, args, **kwargs):
            self.args = args
            streams.extend((kwargs["stdout"], kwargs["stderr"]))
            assert all(stream != subprocess.PIPE for stream in streams)

        def wait(self, timeout):
            calls.append(("wait", timeout))
            if len(calls) == 1:
                raise primary
            raise cleanup

        def kill(self):
            calls.append(("kill",))
            if cleanup_failure == "kill":
                raise cleanup

        def communicate(self, *args, **kwargs):
            pytest.fail("Pipe communication can wait forever for inherited handles")

    monkeypatch.setattr(module.subprocess, "Popen", Process)
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._mode = "wsl"
    wrapper.cfg = SimpleNamespace(timeout=0.01)
    with pytest.raises(WrapperTimeoutError) as captured:
        wrapper._wsl("owned fixture")
    assert captured.value.__cause__ is primary
    assert primary.__cause__ is cleanup
    assert calls[:2] == [("wait", 0.01), ("kill",)]
    if cleanup_failure == "wait":
        assert calls[2] == ("wait", 1.0)
    assert all(stream.closed for stream in streams)

def test_wsl_timeout_retires_runtime_before_execution_or_cache_reuse(monkeypatch):
    wrapper = object.__new__(BaseWrapper)
    wrapper._init_runtime_state()
    wrapper._mode = "wsl"
    wrapper.cfg = SimpleNamespace(timeout=0.01)
    cleaned = []
    wrapper._wsl_cleanup = lambda: cleaned.append("private workspace")

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired("owned fixture", 0.01)

    monkeypatch.setattr(module, "_run_captured_command", timeout)
    with pytest.raises(WrapperTimeoutError):
        wrapper._wsl("cache mutation")
    assert wrapper._closed
    assert cleaned == ["private workspace"]
    with pytest.raises(RuntimeError, match="closed"):
        wrapper._wsl("any later command")
    with pytest.raises(RuntimeError, match="closed"):
        wrapper.run_stdin(Path("iso"), "QUIT")
    with pytest.raises(RuntimeError, match="closed"):
        list_isotropy_cache(wrapper)
    with pytest.raises(RuntimeError, match="closed"):
        wrapper._ensure_writable_isodata()

@pytest.fixture
def isolated_transport(tmp_path, monkeypatch):
    """Simulate only external allocation/linking; exercise real local copies."""
    source = tmp_path / "source"
    source.mkdir()
    (source / "database").write_text("reference", encoding="utf-8")
    temporary = tmp_path / "runtime"
    temporary.mkdir()
    cfg = SimpleNamespace(
        timeout=5,
        temp_dir=temporary,
        _cfg={"isobyu": {"data_dir": str(source)}},
        resolve_path=Path,
    )
    paths = []
    links = {}
    removed = []

    def transport(_self, command, timeout=None):
        if command.startswith("umask 077; mktemp"):
            path = f"/var/tmp/id-{len(paths):08d}"
            paths.append(path)
            return subprocess.CompletedProcess(command, 0, path + "\n", "")
        tokens = shlex.split(command)
        if tokens[:2] == ["ln", "-sfn"]:
            links[tokens[3]] = Path(tokens[2])
            return subprocess.CompletedProcess(command, 0, "", "")
        raise AssertionError(f"Unexpected shell command: {command}")

    def cleanup(path, _timeout):
        removed.append(path)
        links.pop(path + "/data", None)

    monkeypatch.setattr(module, "get_config", lambda: cfg)
    monkeypatch.setattr(BaseWrapper, "_detect_mode", staticmethod(lambda: "wsl"))
    monkeypatch.setattr(BaseWrapper, "_wsl", transport)
    monkeypatch.setattr(BaseWrapper, "_win_to_wsl", lambda _self, path: str(path))
    monkeypatch.setattr(module, "_remove_wsl_workspace", cleanup)
    return SimpleNamespace(source=source, paths=paths, links=links, removed=removed)

def test_interleaved_wrappers_cannot_retarget_generated_database(isolated_transport):
    transport = isolated_transport
    with BaseWrapper() as first:
        first._ensure_writable_isodata()
        private = transport.links[first._isodata_path().rstrip("/")]
        (private / "database").write_text("generated A", encoding="utf-8")
        with BaseWrapper() as second:
            assert first._stage_dir != second._stage_dir
            assert first._isodata_path() != second._isodata_path()
            first._ensure_writable_isodata()
            assert transport.links[first._isodata_path().rstrip("/")] == private
            assert transport.links[second._isodata_path().rstrip("/")] == transport.source
            second._ensure_writable_isodata()
            private_b = transport.links[second._isodata_path().rstrip("/")]
            assert private_b != private
            assert (private_b / "database").read_text(encoding="utf-8") == "reference"
        assert private.exists()
        assert first._stage_dir not in transport.removed
    assert not private.exists()
    assert not private_b.exists()
    assert (transport.source / "database").read_text(encoding="utf-8") == "reference"
    assert len(transport.removed) == 2

def test_close_is_idempotent_and_never_reopens_cache(isolated_transport):
    wrapper = BaseWrapper()
    stage = wrapper._stage_dir
    wrapper.close()
    wrapper.close()
    assert isolated_transport.removed == [stage]
    with pytest.raises(RuntimeError, match="closed"):
        wrapper.run_stdin(Path("iso"), "QUIT")
    with pytest.raises(RuntimeError, match="closed"):
        wrapper._ensure_writable_isodata()
    with pytest.raises(RuntimeError, match="closed"):
        list_isotropy_cache(wrapper)

def test_abandoned_wrapper_releases_only_its_own_workspace(isolated_transport):
    first = BaseWrapper()
    second = BaseWrapper()
    stage = first._stage_dir
    del first
    gc.collect()
    assert isolated_transport.removed == [stage]
    assert second._stage_dir not in isolated_transport.removed
    second.close()

def test_failed_database_binding_cleans_partial_workspace(isolated_transport, monkeypatch):
    def fail_binding(_self, data_dir_win=None):
        raise _wrapper_runtime_WrapperRunError("wsl", 1, "link denied")

    monkeypatch.setattr(BaseWrapper, "_wsl_link_isodata", fail_binding)
    with pytest.raises(_wrapper_runtime_WrapperRunError, match="link denied"):
        BaseWrapper()
    assert isolated_transport.removed == isolated_transport.paths

def test_initialization_preserves_primary_error_when_cleanup_fails(isolated_transport, monkeypatch):
    primary_error = _wrapper_runtime_WrapperRunError("wsl", 1, "link denied")
    cleanup_error = OSError("workspace cleanup denied")
    original_cleanup = module._remove_wsl_workspace

    def fail_binding(_self, data_dir_win=None):
        raise primary_error

    def fail_cleanup(path, timeout):
        original_cleanup(path, timeout)
        raise cleanup_error

    monkeypatch.setattr(BaseWrapper, "_wsl_link_isodata", fail_binding)
    monkeypatch.setattr(module, "_remove_wsl_workspace", fail_cleanup)
    with pytest.raises(_wrapper_runtime_WrapperRunError, match="link denied") as captured:
        BaseWrapper()
    assert captured.value is primary_error
    assert captured.value.__cause__ is cleanup_error
    assert isolated_transport.removed == isolated_transport.paths

def test_runtime_clone_preserves_primary_error_when_cleanup_fails(isolated_transport, monkeypatch):
    primary_error = ValueError("cache cannot be copied")
    cleanup_error = OSError("workspace cleanup denied")
    with BaseWrapper() as source:
        target = BaseWrapper()
        original_cleanup = target._wsl_cleanup

        def fail_copy(_source):
            raise primary_error

        def fail_cleanup():
            original_cleanup()
            raise cleanup_error

        monkeypatch.setattr(target, "_copy_runtime_cache_from", fail_copy)
        monkeypatch.setattr(target, "_wsl_cleanup", fail_cleanup)
        with pytest.raises(ValueError, match="cache cannot be copied") as captured:
            target.clone_runtime_from(source)
        assert captured.value is primary_error
        assert captured.value.__cause__ is cleanup_error
        assert target._closed
        assert not source._closed
        assert target._stage_dir in isolated_transport.removed

def test_cleanup_rejects_nonowned_paths_without_launching_process(monkeypatch):
    def unexpected(*_args, **_kwargs):
        pytest.fail("Unsafe cleanup launched a process")

    monkeypatch.setattr(subprocess, "run", unexpected)
    for path in ("/tmp", "/var/tmp", "/", "/var/tmp/id-12345678/..", "/home/user/.id", "E:\\CRIS"):
        with pytest.raises(ValueError, match="unowned"):
            module._remove_wsl_workspace(path, 1)

@pytest.mark.skipif(os.name != "nt", reason="real WSL transport requires Windows")
def test_real_wsl_interleaved_and_parallel_workspaces(tmp_path, monkeypatch):
    """Use an independent tiny database, never generate into packaged ISODATA."""
    source = tmp_path / "source"
    source.mkdir()
    reference = source / "database"
    reference.write_bytes(b"reference\n")
    runtime = tmp_path / "inputs"
    runtime.mkdir()
    cfg = SimpleNamespace(
        timeout=20,
        temp_dir=runtime,
        _cfg={"isobyu": {"data_dir": str(source)}},
        resolve_path=Path,
    )
    monkeypatch.setattr(module, "get_config", lambda: cfg)
    before = reference.read_bytes()
    # Construction failure is a test failure, so environment errors cannot be
    # mistaken for acceptance of the isolation contract on the supported host.
    with BaseWrapper() as first:
        first._ensure_writable_isodata()
        with BaseWrapper() as second:
            first._ensure_writable_isodata()
            second._ensure_writable_isodata()
            stages = [first._stage_dir, second._stage_dir]
            assert len(set(stages)) == 2
            assert all(len(stage.encode("ascii")) == 20 for stage in stages)
            for wrapper in (first, second):
                result = wrapper._wsl(f"stat -c %a {shlex.quote(wrapper._stage_dir)}")
                assert result.returncode == 0 and result.stdout == "700"

            def invoke(pair):
                wrapper, content = pair
                script = (
                    'printf "%s\\n" ' + shlex.quote(content) + ' > "$ISODATA/database"\ncat "$ISODATA/database"\npwd\n'
                )
                output = wrapper.run_stdin(PurePosixPath("/bin/sh"), script).splitlines()
                assert output == [content, wrapper._stage_dir]
                staged = wrapper._stage_text("fs", "fixture\n")
                assert len(staged.encode("ascii")) < 48
                return output[0]

            with ThreadPoolExecutor(max_workers=2) as pool:
                assert list(pool.map(invoke, [(first, "first"), (second, "second")])) == ["first", "second"]
            assert reference.read_bytes() == before
        assert first._wsl(f"test -d {shlex.quote(stages[0])}").returncode == 0
        assert first._wsl(f"test ! -e {shlex.quote(stages[1])}").returncode == 0
    assert reference.read_bytes() == before
    assert not list(runtime.glob("isodata-*"))

@pytest.mark.skipif(os.name != "nt", reason="real WSL transport requires Windows")
def test_real_wsl_staging_survives_natural_distro_idle_lifecycle(tmp_path, monkeypatch):
    """The Windows owner can outlive WSL's automatic idle distro shutdown."""
    source = tmp_path / "source"
    source.mkdir()
    (source / "database").write_bytes(b"reference\n")
    cfg = SimpleNamespace(
        timeout=20, temp_dir=tmp_path,
        _cfg={"isobyu": {"data_dir": str(source)}}, resolve_path=Path,
    )
    monkeypatch.setattr(module, "get_config", lambda: cfg)
    with BaseWrapper() as wrapper:
        wrapper._ensure_writable_isodata()
        staged = wrapper._stage_text("fs", "owned input\n")
        assert len(staged.encode("ascii")) < 48
        # No Linux process is kept alive, and no shutdown/restart command is
        # issued. On this host /tmp is tmpfs and disappears during this gap.
        time.sleep(20)
        assert wrapper._wsl(f"cat {shlex.quote(staged)}").stdout == "owned input"
        assert wrapper.run_stdin(
            PurePosixPath("/bin/sh"), 'cat "$ISODATA/database"\n'
        ) == "reference"
        assert wrapper._wsl_cleanup.alive
        assert not wrapper._closed
    assert (source / "database").read_bytes() == b"reference\n"

@pytest.mark.skipif(os.name != "nt", reason="real WSL transport requires Windows")
def test_real_wsl_transaction_clone_keeps_generated_tables_after_parent_close(tmp_path, monkeypatch):
    source_data = tmp_path / "reference"
    source_data.mkdir()
    reference = source_data / "database"
    reference.write_bytes(b"reference\n")
    cfg = SimpleNamespace(
        timeout=20, temp_dir=tmp_path,
        _cfg={"isobyu": {"data_dir": str(source_data)}}, resolve_path=Path,
    )
    monkeypatch.setattr(module, "get_config", lambda: cfg)
    with BaseWrapper() as parent, BaseWrapper() as transaction:
        parent._ensure_writable_isodata()
        parent.run_stdin(
            PurePosixPath("/bin/sh"),
            'printf "generated\\n" > "$ISODATA/database"\n'
            'printf "Space Group 139, Irrep LD1, kparam 1/6\\n" > i123.iso\n',
        )
        assert [entry.name for entry in list_isotropy_cache(parent)] == ["i123.iso"]
        transaction.clone_runtime_from(parent)
        parent.close()
        # This step performs no generation or source recreation: it must use
        # the transaction's independent copy of the previous step's state.
        assert transaction.run_stdin(
            PurePosixPath("/bin/sh"), 'cat "$ISODATA/database"\ncat i123.iso\n'
        ).splitlines() == ["generated", "Space Group 139, Irrep LD1, kparam 1/6"]
        entries = list_isotropy_cache(transaction)
        assert [(entry.name, entry.parent_sg, entry.irrep, entry.kparam) for entry in entries] == [
            ("i123.iso", 139, "LD1", "1/6")
        ]
    assert reference.read_bytes() == b"reference\n"
    assert not list(tmp_path.glob("isodata-*"))

# Original coverage: test_wsl.py.

@_wsl_required
# --- from test_real_binaries.py ---

def test_findsym_identifies_nacl():
    """findsym 识别 NaCl（F 心）应为 Fm-3m #225。"""
    fs = FindsymWrapper()
    result = fs.identify(
        lattice_params=[5.63, 5.63, 5.63, 90, 90, 90],
        atom_types=["Na", "Cl"],
        atom_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        centering="F",
    )
    assert result.space_group_number == 225
    assert result.space_group_symbol == "Fm-3m"
    assert {s["wyckoff_letter"] for s in result.wyckoff_sites} == {"a", "b"}

@_wsl_required
def test_iso_kpoints_and_subgroups():
    """iso 枚举 SG 225 的 k 点与 GM5- 子群。"""
    iso = _wsl_IsoWrapper()
    kpoints = iso.list_k_points(225)
    labels = {kp.label for kp in kpoints}
    assert {"GM", "L", "X", "W"}.issubset(labels)

    subgroups = iso.list_subgroups(225, "GM", "GM5-")
    assert len(subgroups) >= 1
    # GM5- P1 对应子群 I-42m (#121)，指数 6
    p1 = next(sg for sg in subgroups if sg.opd_symbol == "P1")
    assert p1.space_group_number == 121
    assert p1.subgroup_index == 6

@_wsl_required
def test_iso_modes_and_domains():
    """BUSH 模式基矢与畴列表。"""
    iso = _wsl_IsoWrapper()

    target = _wsl_SubgroupInfo(
        index=0, space_group_number=107, space_group_symbol="I4mm",
        subgroup_index=6, size=1, is_maximal=True,
        opd_symbol="P1", opd_vector=[1.0, 0.0, 0.0],
        basis_vectors=[[0, 0.5, -0.5], [0, 0.5, 0.5], [1, 0, 0]],
        origin=[0, 0, 0], k_point_label="GM", irrep_label="GM4-",
    )
    modes = iso.calc_distortion_modes(225, target, wyckoff_letters=["a"])
    assert len(modes) >= 1
    assert modes[0].bush_modes, "GM4- P1 在 4a 位点应存在位移模式"

    domains = iso.get_domains(225, target)
    assert len(domains) == 6
    assert domains[0].domain_number == 1
