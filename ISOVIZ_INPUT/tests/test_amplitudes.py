from __future__ import annotations

from pathlib import Path

from isoviz_input.amplitudes import apply_amplitudes, patch_isoviz_file, read_amplitude_csv

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
