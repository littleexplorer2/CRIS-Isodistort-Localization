from pathlib import Path

from isoviz_input.amplitudes import apply_amplitudes, read_amplitude_csv
from isoviz_input.validation import section_scalar, validate_patched_text

FIXTURES = Path(__file__).parent / "fixtures"


def test_static_validation_proves_mode_identity_and_amplitudes():
    original = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    patched, patch_report = apply_amplitudes(
        original,
        read_amplitude_csv(FIXTURES / "sample.csv"),
    )

    result = validate_patched_text(original, patched, patch_report)

    assert result["status"] == "pass_with_warnings"
    assert result["identity_preserved"] is True
    assert result["non_amplitude_content_preserved"] is True
    assert result["mode_count"] == 3
    assert result["matched_count"] == 2
    assert result["missing_expected"] == []
    assert result["metadata"]["isoversion"] == "6.12"
    assert result["expected_visible_amplitudes"] == [
        {"label": "[0,0,1/6]LD1[Eu1:a:dsp]A2u(a)", "amp": 0.12345},
        {"label": "[0,0,1/3]LD1[Eu1:a:dsp]A2u(a)", "amp": 0.5},
    ]


def test_static_validation_fails_when_patched_value_is_not_present():
    original = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    _patched, patch_report = apply_amplitudes(
        original,
        read_amplitude_csv(FIXTURES / "sample.csv"),
    )

    result = validate_patched_text(original, original, patch_report)

    assert result["status"] == "fail"
    assert result["missing_expected"]


def test_static_validation_fails_if_non_amplitude_content_changes():
    original = (FIXTURES / "sample.isoviz").read_text(encoding="utf-8")
    patched, patch_report = apply_amplitudes(
        original,
        read_amplitude_csv(FIXTURES / "sample.csv"),
    )
    patched = patched.replace("!isoversion 6.12", "!isoversion 9.99")

    result = validate_patched_text(original, patched, patch_report)

    assert result["status"] == "fail"
    assert result["non_amplitude_content_preserved"] is False


def test_section_scalar_returns_none_for_absent_metadata():
    assert section_scalar("!isoversion\n6.12\n", "child_space_group_number") is None
