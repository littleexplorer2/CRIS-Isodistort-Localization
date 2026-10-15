from __future__ import annotations

from pathlib import Path

import pytest

from isoviz_input.application.pipeline import service as pipeline
from isoviz_input.application.pipeline import PipelineRequest, dispatch_to_isoviz, prepare_isoviz
from isoviz_input.backend.amplitudes import parse_amplitude_csv
from isoviz_input.backend.inputs import document_from_text, load_text_document
from isoviz_input.frontend.gui import format_preparation_report

FIXTURES = Path(__file__).parent / "fixtures"


def _pasted_request() -> PipelineRequest:
    return PipelineRequest(
        amplitudes=document_from_text(
            (FIXTURES / "sample.csv").read_text(encoding="utf-8"),
            source="pasted CSV",
        ),
        structure=document_from_text(
            (FIXTURES / "sample.isoviz").read_text(encoding="utf-8"),
            source="pasted structure.isoviz",
        ),
    )


def test_pasted_text_uses_complete_pipeline(tmp_path):
    prepared = prepare_isoviz(_pasted_request(), output_dir=tmp_path)

    assert prepared.launchable is True
    assert prepared.validation["status"] == "pass_with_warnings"
    assert prepared.output_path is not None
    assert prepared.output_path.parent == tmp_path
    assert "0.12345" in prepared.output_path.read_text(encoding="utf-8")
    assert "0.50000" in prepared.patched_text
    report = format_preparation_report(prepared)
    assert "匹配模式：2" in report
    assert "pasted CSV" in report


def test_pipeline_blocks_launch_when_no_mode_matches(tmp_path):
    request = PipelineRequest(
        amplitudes=document_from_text(
            "Mode,Mode Name,Best Model Parameter\nz9,not-a-real-mode,0.5\n",
            source="bad CSV",
        ),
        structure=_pasted_request().structure,
    )

    prepared = prepare_isoviz(request, output_dir=tmp_path)

    assert prepared.validation["status"] == "fail"
    assert prepared.launchable is False
    assert prepared.output_path is None
    with pytest.raises(ValueError, match="not launchable"):
        dispatch_to_isoviz(prepared)


def test_dispatch_receipt_is_limited_to_process_dispatch(tmp_path, monkeypatch):
    prepared = prepare_isoviz(_pasted_request(), output_dir=tmp_path)
    launcher = tmp_path / "IsoVIZ.exe"
    launcher.write_bytes(b"exe")
    calls: list[tuple[Path, Path | None]] = []
    monkeypatch.setattr(
        pipeline,
        "open_isoviz",
        lambda path, launcher=None: calls.append((path, launcher)),
    )

    receipt = dispatch_to_isoviz(prepared, launcher=launcher)

    assert calls == [(prepared.output_path, launcher)]
    assert receipt.prepared_path == prepared.output_path
    assert receipt.launcher == launcher


def test_file_upload_reader_enforces_expected_suffix(tmp_path):
    csv_path = tmp_path / "values.csv"
    csv_path.write_text("Mode,Best Model Parameter\na1,0.1\n", encoding="utf-8")

    document = load_text_document(csv_path, expected_suffix=".csv")

    assert document.path == csv_path.resolve()
    assert document.source == str(csv_path.resolve())
    with pytest.raises(ValueError, match=r"Expected a \.isoviz"):
        load_text_document(csv_path, expected_suffix=".isoviz")


@pytest.mark.parametrize("value", ["NaN", "inf", "-inf"])
def test_parser_rejects_non_finite_amplitudes(value):
    text = f"Mode,Best Model Parameter\na1,{value}\n"

    with pytest.raises(ValueError, match="must be finite"):
        parse_amplitude_csv(text, source="pasted CSV")
