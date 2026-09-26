from __future__ import annotations

from pathlib import Path

from isoviz_input.amplitudes import apply_amplitudes, patch_isoviz_file, read_amplitude_csv
from isoviz_input.paths import ensure_input_content

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def test_read_gd_csv():
    modes = read_amplitude_csv(FIXTURES / "sample.csv")
    assert len(modes) == 2
    assert modes[0].alias == "a1"
    assert modes[0].name.startswith("[0,0,1/6]")
    assert abs(modes[0].amplitude - 0.12345) < 1e-9
    assert modes[0].max_amplitude is not None
    assert abs(modes[0].max_amplitude - 2.44949) < 1e-9


def test_apply_by_mode_name():
    text = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    modes = read_amplitude_csv(FIXTURES / "sample.csv")
    patched, report = apply_amplitudes(text, modes)
    assert len(report.matched) == 2
    assert report.unmatched_csv == []
    assert "1    1   0.12345   2.44949" in patched
    assert "1    2   0.50000   2.82843" in patched
    assert report.unmatched_isoviz == ["GM1+strain_1(a)"]
    assert "  1    0.00000   0.10000" in patched


def test_apply_by_alias_when_names_missing():
    text = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    from isoviz_input.amplitudes import ModeAmplitude

    modes = [
        ModeAmplitude(name="", amplitude=0.01, alias="a1"),
        ModeAmplitude(name="", amplitude=0.02, alias="a2"),
        ModeAmplitude(name="", amplitude=0.03, alias="a3"),
    ]
    patched, report = apply_amplitudes(text, modes)
    assert len(report.matched) == 3
    assert "1    0.01000   0.10000" in patched
    assert "1    1   0.02000   2.44949" in patched
    assert "1    2   0.03000   2.82843" in patched


def test_patch_isoviz_file(tmp_path):
    dest = tmp_path / "out.isoviz"
    report = patch_isoviz_file(
        FIXTURES / "sample.isoviz",
        FIXTURES / "sample.csv",
        dest,
    )
    assert dest.is_file()
    assert len(report.matched) == 2
    body = dest.read_text(encoding="utf-8")
    assert "0.12345" in body
    assert "0.50000" in body


def test_ensure_input_content_creates_named_folders(tmp_path, monkeypatch):
    monkeypatch.setattr("isoviz_input.paths.ISOVIZ_ROOT", tmp_path)
    monkeypatch.setattr("isoviz_input.paths.INPUT_ROOT", tmp_path / "input_content")
    monkeypatch.setattr("isoviz_input.paths.DATA_DIR", tmp_path / "input_content" / "data.csv")
    monkeypatch.setattr(
        "isoviz_input.paths.STRUCTURE_DIR",
        tmp_path / "input_content" / "subgroup.isoviz",
    )
    data_dir, structure_dir = ensure_input_content()
    assert data_dir.is_dir()
    assert structure_dir.is_dir()
    assert data_dir.name == "data.csv"
    assert structure_dir.name == "subgroup.isoviz"


def test_gd_style_names_match_isoviz_labels():
    text = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    from isoviz_input.amplitudes import ModeAmplitude

    modes = [
        ModeAmplitude(
            name="I4/mmm[0,0,1/6]LD1(a,b)[Eu0:a:dsp]A2u(a)",
            amplitude=0.12345,
            alias="a1",
        ),
        ModeAmplitude(
            name="I4/mmm[0,0,1/3]LD1(a,b)[Eu0:a:dsp]A2u(a)",
            amplitude=0.5,
            alias="a2",
        ),
    ]
    patched, report = apply_amplitudes(text, modes)
    assert len(report.matched) == 2
    assert report.unmatched_csv == []
    assert "1    1   0.12345   2.44949" in patched
    assert "1    2   0.50000   2.82843" in patched
    assert report.unmatched_isoviz == ["GM1+strain_1(a)"]


def test_displacive_alias_fallback_skips_strain():
    text = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    from isoviz_input.amplitudes import ModeAmplitude

    modes = [
        ModeAmplitude(name="unmatched-name-a", amplitude=0.11, alias="a1"),
        ModeAmplitude(name="unmatched-name-b", amplitude=0.22, alias="a2"),
    ]
    patched, report = apply_amplitudes(text, modes)
    assert len(report.matched) == 2
    assert "1    1   0.11000   2.44949" in patched
    assert "1    2   0.22000   2.82843" in patched
    assert report.unmatched_isoviz == ["GM1+strain_1(a)"]


def test_main_uses_quoted_structure_without_launching(monkeypatch):
    import main as isoviz_main

    launched: dict[str, str] = {}

    def fake_open(path, launcher=None):
        launched["text"] = Path(path).read_text(encoding="utf-8")

    monkeypatch.setattr(isoviz_main, "open_isoviz", fake_open)
    monkeypatch.setattr(isoviz_main, "find_isoviz_launcher", lambda: None)
    rc = isoviz_main.main(
        [
            "--data",
            str(FIXTURES / "sample.csv"),
            "--structure",
            f'"{FIXTURES / "sample.isoviz"}"',
        ]
    )
    assert rc == 0
    assert "0.12345" in launched["text"]
    assert "0.50000" in launched["text"]
