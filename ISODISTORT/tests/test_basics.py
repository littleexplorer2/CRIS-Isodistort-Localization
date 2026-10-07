"""Basics: config, i18n, parent header, official k-params."""
from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from pymatgen.core import Lattice, Structure

from backend.api import IsoDistort
from backend.tables.kpoints_official import official_kparams_to_iso
from backend.utils import IsodistortError, WrapperRunError, get_config
from backend.utils.parent_header import (
    format_fixed_coord,
    format_wyckoff_site,
    format_wyckoff_sites,
    format_wyckoff_sites_from_cif,
)
from backend.wrappers import KPointInfo
from backend.wrappers.base_wrapper import (
    BaseWrapper,
    _decode_process_output,
    _is_wsl_access_denied,
)
from features.input_cif import SymmetryValidator
from frontend.i18n import MESSAGES, t

# --- from test_config.py ---

def test_config_load():
    cfg = get_config()
    # isobyu 二进制路径可解析（存在性由部署保证；WSL 下为 Windows 侧路径）
    assert cfg.iso_bin.name == "iso"
    assert cfg.findsym_bin.name == "findsym"
    assert cfg.temp_dir.exists(), "临时目录已创建"
    assert cfg.output_dir.exists(), "输出目录已创建"
    assert "ISODATA" in os.environ, "ISODATA 环境变量已设置"
    print("配置加载测试通过")
    print(f"   ISODATA = {os.environ['ISODATA']}")


def test_binary_launch_requires_the_staging_directory(monkeypatch):
    wrapper = object.__new__(BaseWrapper)
    wrapper._stage_dir = "/home/test/.id/tmp"
    monkeypatch.setattr(
        wrapper,
        "_stage_text",
        lambda _prefix, _text: "/home/test/.id/tmp/iso_input.in",
    )
    monkeypatch.setattr(wrapper, "_wsl_bin_path", lambda _path: "/opt/isobyu/iso")
    monkeypatch.setattr(wrapper, "_isodata_path", lambda: "/home/test/.id/data/")
    commands: list[str] = []

    def capture(command: str, timeout: float | None = None):
        _ = timeout
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(wrapper, "_wsl", capture)

    assert wrapper._run_program(Path("iso"), "QUIT\n", True) == "ok"
    assert len(commands) == 1
    assert (
        "export ISODATA=/home/test/.id/data/ && "
        "cd /home/test/.id/tmp && /opt/isobyu/iso"
    ) in commands[0]
    assert "; cd " not in commands[0]


def test_wsl_stream_decoder_separates_linux_output_from_windows_diagnostics():
    assert _decode_process_output(b"/home/test\n") == "/home/test"
    diagnostic = "Access is denied.\r\nError code: Wsl/E_ACCESSDENIED\r\n"
    decoded = _decode_process_output(diagnostic.encode("utf-16-le"))
    assert decoded == "Access is denied.\nError code: Wsl/E_ACCESSDENIED"
    assert _is_wsl_access_denied(decoded)


def test_wsl_get_home_reports_onedrive_process_restriction(monkeypatch):
    wrapper = object.__new__(BaseWrapper)
    denied = "Access is denied.\nError code: Wsl/E_ACCESSDENIED"
    monkeypatch.setattr(
        wrapper,
        "_wsl",
        lambda _command: SimpleNamespace(returncode=0xFFFFFFFF, stdout=denied, stderr=""),
    )

    with pytest.raises(WrapperRunError) as excinfo:
        wrapper._wsl_get_home()

    message = str(excinfo.value)
    assert "Wsl/E_ACCESSDENIED" in message
    assert ".\\run_cris.ps1 <原 Python 参数>" in message
    assert "不是缺少管理员权限" in message


def test_crystallographic_tolerances_have_explicit_units_and_roles():
    cfg = get_config()

    assert cfg.symmetry_cartesian_tolerance_angstrom == 1e-3
    assert cfg.position_tolerance == cfg.symmetry_cartesian_tolerance_angstrom
    assert cfg.symmetry_angle_tolerance_degrees == 5.0
    assert cfg.affine_exact_cartesian_tolerance_angstrom == 1e-7
    assert cfg.fractional_coordinate_tolerance == 1e-7
    assert cfg.lattice_tolerance == 1e-5
    # Exact synthetic affine data have a tighter Cartesian error budget than
    # uploaded/refined structures; matrix and fractional residuals are
    # dimensionless and are intentionally separate values.
    assert (
        cfg.affine_exact_cartesian_tolerance_angstrom
        < cfg.symmetry_cartesian_tolerance_angstrom
    )


@pytest.mark.parametrize("generate_if_missing", [False, True])
def test_method2_all_ir_inventory_rejects_any_partial_backend_failure(
    generate_if_missing: bool,
) -> None:
    """A successful IR must not hide a failed IR in the advertised full table."""

    old_candidate = SimpleNamespace(index=7)
    first_candidate = SimpleNamespace(index=0)

    class PartialBackend:
        @staticmethod
        def list_irreps(*_args):
            return [SimpleNamespace(label="A1"), SimpleNamespace(label="A2")]

        @staticmethod
        def list_subgroups(*_args, **_kwargs):
            if _args[2] == "A1":
                return [first_candidate]
            raise IsodistortError("second IR parse failed")

    api = object.__new__(IsoDistort)
    api.structure = object()
    api.symmetry_info = {"space_group_number": 123}
    api.subgroups = [old_candidate]
    api._iso = PartialBackend()
    api._resolve_iso_kparams = lambda *_args: []
    api._filter_subgroups_for_search = lambda rows, *_args: rows
    api._tag_official_kparams = lambda *_args: None

    with pytest.raises(IsodistortError, match="second IR parse failed"):
        api.list_subgroups_at_kpoint(
            "X",
            generate_if_missing=generate_if_missing,
        )

    assert api.subgroups == [old_candidate]


def test_method2_all_ir_inventory_allows_a_genuine_empty_ir() -> None:
    candidate = SimpleNamespace(index=99)

    class CompleteBackend:
        @staticmethod
        def list_irreps(*_args):
            return [SimpleNamespace(label="A1"), SimpleNamespace(label="A2")]

        @staticmethod
        def list_subgroups(*_args, **_kwargs):
            return [] if _args[2] == "A1" else [candidate]

    api = object.__new__(IsoDistort)
    api.structure = object()
    api.symmetry_info = {"space_group_number": 123}
    api.subgroups = []
    api._iso = CompleteBackend()
    api._resolve_iso_kparams = lambda *_args: []
    api._filter_subgroups_for_search = lambda rows, *_args: rows
    api._tag_official_kparams = lambda *_args: None

    result = api.list_subgroups_at_kpoint("X")

    assert result == [candidate]
    assert candidate.index == 0
    assert api.subgroups == [candidate]


def test_orbit_identity_is_periodic_at_fractional_cell_boundary():
    lattice = Lattice.cubic(4.0)
    exact = Structure(lattice, ["He"], [[0.0, 0.25, 0.5]])
    wrapped_noise = Structure(lattice, ["He"], [[-1e-12, 0.25, 0.5]])
    validator = SymmetryValidator()

    exact_id = validator.validate(exact)["wyckoff_sites"][0]["orbit_id"]
    wrapped_id = validator.validate(wrapped_noise)["wyckoff_sites"][0]["orbit_id"]
    assert exact_id == wrapped_id




# --- from test_i18n.py ---

def test_messages_are_flat_english():
    assert isinstance(MESSAGES, dict)
    assert "zh" not in MESSAGES
    assert "en" not in MESSAGES
    assert t("load.done", sg=139, sym="I4/mmm", n=10).startswith("[Loaded]")
    assert "Space group" in t("load.done", sg=139, sym="I4/mmm", n=10)
    assert "139" in t("load.done", sg=139, sym="I4/mmm", n=10)


def test_unknown_key_returns_key():
    assert t("no.such.key") == "no.such.key"


def test_no_language_switch_keys():
    assert "ui.menu.language" not in MESSAGES
    assert "ui.lang.current" not in MESSAGES
    assert "dist.gen" not in MESSAGES
    assert "dist.domainsBtn" not in MESSAGES


def test_no_terminal_ui_message_keys():
    for key in (
        "ui.banner.title",
        "ui.menu.search",
        "prefs.terminalBlock",
        "m2.nmodRemoved",
        "dist.computeModesAsk",
    ):
        assert key not in MESSAGES


def test_distortion_and_method2_help_keys():
    assert "m2.genDbHelp" in MESSAGES
    assert "isotropy-subgroup" in MESSAGES["m2.genDbHelp"].lower() or (
        "isotropy" in MESSAGES["m2.genDbHelp"].lower()
        and "subgroup" in MESSAGES["m2.genDbHelp"].lower()
    )
    assert "cache" in MESSAGES["m2.genDbHelp"].lower()
    assert "\n" not in MESSAGES["m2.genDbHelp"]
    assert "m2.genDbManage" in MESSAGES
    assert "m2.genDbDelete" in MESSAGES
    assert "dist.method4" in MESSAGES
    assert "dist.tableLabel" in MESSAGES
    assert "not limited by the filter" not in MESSAGES["dist.tableNote"]
    assert "dist.zipNotM4" in MESSAGES
    assert "dist.zipWait" in MESSAGES
    assert "dist.zipDone" in MESSAGES
    assert "st.elapsed" in MESSAGES
    assert "st.busyHint" in MESSAGES
    assert "independent incommensurate" in MESSAGES["l.nmod"].lower()
    assert "ok.nmod" in MESSAGES
    assert "ss.title" not in MESSAGES
    assert "ui.menu.superspace" not in MESSAGES


def test_numbered_web_placeholders():
    text = t("ok.nsup", 3)
    assert "3" in text


# --- from test_parent_header.py ---

def test_format_fixed_coord_fractions():
    assert format_fixed_coord(0.0) == "0"
    assert format_fixed_coord(0.5) == "1/2"
    assert format_fixed_coord(0.25) == "1/4"


def test_format_wyckoff_site_free_z():
    line = format_wyckoff_site("Al", 2, 4, "e", [0.0, 0.0, 0.38])
    assert line == "Al2 4e (0,0,z), z= 0.38000"


def test_format_wyckoff_sites_eual4_like():
    lattice = Lattice.tetragonal(4.402, 11.163)
    structure = Structure(
        lattice,
        ["Eu", "Eu", "Al", "Al", "Al", "Al", "Al", "Al", "Al", "Al"],
        [
            [0, 0, 0], [0.5, 0.5, 0.5],
            [0, 0.5, 0.25], [0.5, 0, 0.25], [0.5, 0, 0.75], [0, 0.5, 0.75],
            [0, 0, 0.38], [0, 0, 0.62], [0.5, 0.5, 0.12], [0.5, 0.5, 0.88],
        ],
    )
    sites = [
        {"species": "Eu", "multiplicity": 2, "wyckoff_letter": "a",
         "representative_index": 0},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "d",
         "representative_index": 2},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "e",
         "representative_index": 6},
    ]
    lines = format_wyckoff_sites(structure, sites)
    assert lines == [
        "Eu1 2a (0,0,0)",
        "Al1 4d (0,1/2,1/4)",
        "Al2 4e (0,0,z), z= 0.38000",
    ]


def test_cif_label_matches_nonrepresentative_member_of_orbit(tmp_path):
    structure = Structure(
        Lattice.cubic(4.0),
        ["H", "H"],
        [[0.1, 0, 0], [0.9, 0, 0]],
    )
    sites = [{
        "species": "H",
        "multiplicity": 2,
        "wyckoff_letter": "a",
        "representative_index": 0,
        "equivalent_indices": [0, 1],
    }]
    cif = tmp_path / "equivalent_member.cif"
    cif.write_text(
        """data_test
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
H_second H 0.9 0 0
""",
        encoding="utf-8",
    )

    lines = format_wyckoff_sites_from_cif(cif, structure, sites, tol=0.01)

    assert lines is not None
    assert lines[0].startswith("H_second 2a ")
    assert sites[0]["display_label"] == "H_second"


def test_cif_label_matching_checks_species_and_repeated_orbit_distance(tmp_path):
    structure = Structure(
        Lattice.cubic(4.0),
        ["Na", "Na", "K"],
        [[0.1, 0, 0], [0.2, 0, 0], [0.201, 0, 0]],
    )
    sites = [
        {
            "species": "Na", "multiplicity": 1, "wyckoff_letter": "e",
            "representative_index": 0, "equivalent_indices": [0],
        },
        {
            "species": "Na", "multiplicity": 1, "wyckoff_letter": "e",
            "representative_index": 1, "equivalent_indices": [1],
        },
        {
            "species": "K", "multiplicity": 1, "wyckoff_letter": "e",
            "representative_index": 2, "equivalent_indices": [2],
        },
    ]
    cif = tmp_path / "nearby_orbits.cif"
    cif.write_text(
        """data_test
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
NaB Na 0.2009 0 0
K1 K 0.2010 0 0
NaA Na 0.1000 0 0
""",
        encoding="utf-8",
    )

    lines = format_wyckoff_sites_from_cif(cif, structure, sites, tol=0.01)

    assert lines is not None
    assert [site["display_label"] for site in sites] == ["NaA", "NaB", "K1"]
    assert [line.split()[0] for line in lines] == ["NaB", "K1", "NaA"]


def test_parent_wyckoff_from_cif_ndnio2():
    """CIF atom_site labels/order drive the web header (not memorized EuAl4)."""
    from pathlib import Path

    from backend.api import IsoDistort

    cris = Path(__file__).resolve().parents[2]
    cif = cris / "experiment_data" / "NdNiO2 own.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("NdNiO2 own.cif not in experiment_data")
    iso = IsoDistort()
    iso.load_structure(cif)
    labels_after_load = {
        site.get("display_label")
        for site in iso.symmetry_info["wyckoff_sites"]
    }
    assert {"O", "ND", "NI"} <= labels_after_load
    lines = iso.parent_wyckoff_display()
    assert lines[0].startswith("O 2f")
    assert any(ln.startswith("ND 1d") for ln in lines)
    assert any(ln.startswith("NI 1a") for ln in lines)


def test_parent_wyckoff_from_cif_4310_uses_standard_orbit_representatives():
    """Body-centering-equivalent CIF rows must use FINDSYM standard points."""
    from pathlib import Path

    from backend.api import IsoDistort

    cris = Path(__file__).resolve().parents[2]
    cif = cris / "experiment_data" / "4310_tetra.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("4310_tetra.cif not in experiment_data")

    iso = IsoDistort()
    iso.load_structure(cif)
    sites = iso.symmetry_info["wyckoff_sites"]
    assert len({site["orbit_id"] for site in sites}) == 8
    assert len({site["orbit_id"] for site in sites if site["wyckoff_letter"] == "e"}) == 5
    assert iso.parent_wyckoff_display() == [
        "La1 4e (0,0,z), z= 0.43204",
        "La2 4e (0,0,z), z= 0.30148",
        "Ni1 2a (0,0,0)",
        "Ni2 4e (0,0,z), z=-0.13885",
        "O1 8g (0,1/2,z), z= 0.36070",
        "O3 4c (0,1/2,0)",
        "O4 4e (0,0,z), z=-0.21680",
        "O2 4e (0,0,z), z=-0.06780",
    ]


def test_export_cif_uses_cif_wyckoff_lines_ndnio2():
    """Distortion CIF comments must reuse the same CIF-ordered parent header."""
    from pathlib import Path

    from backend.api import IsoDistort
    from backend.wrappers import SubgroupInfo
    from features.export.distortion_formats import render_cif

    cris = Path(__file__).resolve().parents[2]
    cif = cris / "experiment_data" / "NdNiO2 own.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("NdNiO2 own.cif not in experiment_data")
    iso = IsoDistort()
    iso.load_structure(cif)
    sg = SubgroupInfo(
        index=0,
        space_group_number=139,
        space_group_symbol="I4/mmm",
        size=2,
        subgroup_index=2,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        irrep_label="A1+",
        k_point_label="A",
        parent_sg=123,
        k_coordinates=["1/2", "1/2", "1/2"],
        basis_vectors=[[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
        origin=[0.0, 0.0, 0.0],
    )
    spec = iso._spec_for_subgroup(
        sg, use_current_modes=False, use_generated_structure=False
    )
    assert spec.parent_wyckoff_lines
    assert spec.parent_wyckoff_lines[0].startswith("O 2f")
    text = render_cif(spec.structure, spec)
    assert "# O 2f (0,1/2,0)" in text
    assert "# ND 1d (1/2,1/2,1/2)" in text
    assert "# NI 1a (0,0,0)" in text
    assert "# Nd1" not in text



# --- from test_kparams_official.py ---

def test_ld_g_to_iso():
    iso_kp = KPointInfo(
        label="LD", coordinates=["0", "0", "2a"], parameters=["2a"], is_special=False
    )
    out = official_kparams_to_iso(139, "LD", ["1/6"], iso_kp)
    assert out == ["1/12"]


def test_sm_a_to_iso():
    iso_kp = KPointInfo(
        label="SM", coordinates=["2a", "0", "0"], parameters=["2a"], is_special=False
    )
    out = official_kparams_to_iso(139, "SM", ["1/4"], iso_kp)
    assert out == ["1/8"]


def test_dt_unchanged():
    iso_kp = KPointInfo(
        label="DT", coordinates=["a", "a", "0"], parameters=["a"], is_special=False
    )
    out = official_kparams_to_iso(139, "DT", ["1/4"], iso_kp)
    assert out == ["1/4"]


def test_no_override_passthrough():
    iso_kp = KPointInfo(
        label="Z", coordinates=["0", "0", "1/2"], parameters=[], is_special=True
    )
    out = official_kparams_to_iso(221, "Z", ["1/6"], iso_kp)
    assert out == ["1/6"]
