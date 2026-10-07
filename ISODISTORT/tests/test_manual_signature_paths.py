"""Guard manual provenance signatures against the retired package layout."""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANUAL_ROOT = PROJECT_ROOT / "tests" / "manual"


def test_manual_signature_targets_exist_in_the_four_part_layout() -> None:
    expected = (
        PROJECT_ROOT / "backend" / "api" / "core_api.py",
        PROJECT_ROOT / "backend" / "wrappers" / "iso_wrapper.py",
        PROJECT_ROOT / "backend" / "wrappers" / "smodes_wrapper.py",
        PROJECT_ROOT / "features" / "input_cif" / "__init__.py",
        PROJECT_ROOT / "features" / "input_cif" / "cif_io.py",
        PROJECT_ROOT / "features" / "method1" / "affine_embeddings.py",
        PROJECT_ROOT / "features" / "method1" / "search_methods.py",
        PROJECT_ROOT / "features" / "method2" / "distortion_mapper.py",
        PROJECT_ROOT / "features" / "method2" / "superspace.py",
        PROJECT_ROOT / "features" / "method3" / "inverse_landau.py",
        PROJECT_ROOT / "features" / "method3" / "inverse_landau_adapter.py",
        PROJECT_ROOT / "features" / "method4" / "distortion_engine.py",
        PROJECT_ROOT / "features" / "method4" / "strain.py",
        PROJECT_ROOT / "features" / "export" / "isodistort_cif.py",
        MANUAL_ROOT / "official_html_resolver.py",
        MANUAL_ROOT / "validate_method_outputs.py",
    )

    missing = [str(path) for path in expected if not path.is_file()]
    assert not missing, f"manual signature targets are missing: {missing}"


def test_manual_audits_do_not_sign_retired_package_paths() -> None:
    retired_fragments = (
        '/ "backend" / "distortion" /',
        '/ "backend" / "backend" /',
        '/ "backend" / "structure" /',
        '/ "backend" / "io" /',
        "ISODISTORT/isocore/",
        "ISODISTORT/tests_dev/",
    )
    offenders: list[str] = []
    for path in sorted(MANUAL_ROOT.glob("*.py")):
        source = path.read_text(encoding="utf-8-sig")
        for fragment in retired_fragments:
            if fragment in source:
                offenders.append(f"{path.name}: {fragment}")

    assert not offenders, "retired manual signature paths remain: " + "; ".join(offenders)
