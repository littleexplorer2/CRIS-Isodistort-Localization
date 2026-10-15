"""Guard manual provenance signatures against the retired package layout."""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANUAL_ROOT = PROJECT_ROOT / "tests" / "manual"


def test_manual_signature_targets_exist_in_the_three_layer_layout() -> None:
    expected = (
        PROJECT_ROOT / "application" / "session" / "core_api.py",
        PROJECT_ROOT / "backend" / "wrappers" / "iso_wrapper.py",
        PROJECT_ROOT / "backend" / "wrappers" / "smodes_wrapper.py",
        PROJECT_ROOT / "backend" / "input_cif" / "__init__.py",
        PROJECT_ROOT / "backend" / "input_cif" / "cif_io.py",
        PROJECT_ROOT / "backend" / "method1" / "affine_embeddings.py",
        PROJECT_ROOT / "backend" / "method1" / "search_methods.py",
        PROJECT_ROOT / "backend" / "method2" / "distortion_mapper.py",
        PROJECT_ROOT / "backend" / "method2" / "superspace.py",
        PROJECT_ROOT / "backend" / "method3" / "inverse_landau.py",
        PROJECT_ROOT / "backend" / "method3" / "inverse_landau_adapter.py",
        PROJECT_ROOT / "backend" / "method4" / "distortion_engine.py",
        PROJECT_ROOT / "backend" / "method4" / "strain.py",
        PROJECT_ROOT / "backend" / "export" / "isodistort_cif.py",
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


def test_comparison_signatures_cover_the_relocated_scientific_implementation():
    from tests.manual.audit_method1_4310_live import _validator_source_paths

    paths = _validator_source_paths()
    assert "cif-validator:application/single/service.py" in paths
    assert "cif-validator:backend/cif/service.py" in paths
    assert "cif-validator:backend/comparison/service.py" in paths
    assert "cif-validator:backend/results/__init__.py" in paths
    assert all(path.is_file() for path in paths.values())
    assert not any(path.name == "compare_cif.py" for path in paths.values())
