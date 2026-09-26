"""Tests for the read-only Method output audit helpers."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import pytest

import isocore.api
import tests_dev.manual.validate_method_outputs as validation_module
from tests_dev.manual.validate_method_outputs import (
    _count_local_displacive_modes,
    _live12_checkpoint_policy,
    _live12_signature,
    audit_output_root,
    bases_generate_same_lattice,
    parse_candidate_header,
    run_live_method12_audit,
)


def _basis(rows: tuple[tuple[int, int, int], ...]):
    return tuple(tuple(Fraction(value) for value in row) for row in rows)


def test_equivalent_basis_accepts_integer_unimodular_change() -> None:
    reference = _basis(((1, 0, 0), (0, 1, 0), (0, 0, 2)))
    equivalent = _basis(((1, 1, 0), (0, 1, 0), (0, 0, 2)))
    assert bases_generate_same_lattice(reference, equivalent)


def test_equivalent_basis_rejects_different_sublattice() -> None:
    reference = _basis(((1, 0, 0), (0, 1, 0), (0, 0, 2)))
    different = _basis(((2, 0, 0), (0, 1, 0), (0, 0, 2)))
    assert not bases_generate_same_lattice(reference, different)


def test_candidate_identity_reduces_origin_modulo_integer_translation(tmp_path: Path) -> None:
    template = """# IR: GM1+
# P1 (a) 1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=ORIGIN, s=1, i=1, k-active= (0,0,0)
data_test
"""
    first = tmp_path / "first.cif"
    second = tmp_path / "second.cif"
    first.write_text(template.replace("ORIGIN", "(5/4,-3/2,2)"), encoding="utf-8")
    second.write_text(template.replace("ORIGIN", "(1/4,1/2,0)"), encoding="utf-8")

    assert parse_candidate_header(first).identity == parse_candidate_header(second).identity


def test_local_mode_count_excludes_prose_headings() -> None:
    text = """Mode definitions (every atom)
Mode vectors are given in unitless coordinates
Mode LD1[0,0,1/6]__a__A2u(a)  pretty label
Mode GM5-[0,0,0]__e__Eu(b)  another label
Mode amplitudes (As = supercell-normalized)
"""
    assert _count_local_displacive_modes(text) == 2


def test_saved_method12_audit_delegates_method3_affine_rows(tmp_path: Path) -> None:
    for side in ("官网", "现有网页版交互"):
        candidate = tmp_path / "parent.cif" / side / "Method3" / "affine candidate"
        candidate.mkdir(parents=True)
        # An affine/coupled Method 3 CIF need not carry a single-IR header.
        (candidate / "subgroup.cif").write_text(
            "# Subgroup:  1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, "
            "origin=(0,0,0), s=1, i=1\n",
            encoding="utf-8",
        )

    summary = audit_output_root(tmp_path)

    assert summary.failures == 0
    assert summary.groups == {}


def test_live12_signature_covers_scientific_inputs_but_not_generated_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workspace = tmp_path / "workspace"
    package_root = workspace / "ISODISTORT"
    config_dir = package_root / "config"
    runtime_dir = package_root / "runtime"
    isocore_dir = package_root / "isocore"
    input_dir = workspace / "experiment_data"
    output_root = workspace / "output_compare"
    for directory in (config_dir, runtime_dir, isocore_dir, input_dir):
        directory.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / "settings.yaml"
    config_path.write_text(
        """isobyu:
  bin_dir: ../runtime
  data_dir: ../runtime
  iso_bin: iso
  smodes_bin: smodes
""",
        encoding="utf-8",
    )
    (isocore_dir / "engine.py").write_text("VALUE = 1\n", encoding="utf-8")
    for parent in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        (input_dir / parent).write_text(f"parent {parent}\n", encoding="utf-8")
    iso_path = runtime_dir / "iso"
    smodes_path = runtime_dir / "smodes"
    data_path = runtime_dir / "data_isotropy.txt"
    iso_path.write_bytes(b"iso-v1")
    smodes_path.write_bytes(b"smodes-v1")
    data_path.write_bytes(b"data-v1")
    reference_dir = output_root / "EuAl4 Parent.cif" / "官网" / "Method1" / "candidate"
    reference_dir.mkdir(parents=True)
    (reference_dir / "subgroup.cif").write_text("reference-v1\n", encoding="utf-8")
    details_path = reference_dir / "ISODISTORT complete modes details.html"
    details_path.write_text("details-v1\n", encoding="utf-8")

    runtime = {"python": {"version": "runtime-v1"}}
    monkeypatch.setattr(validation_module, "_live12_runtime_versions", lambda: dict(runtime))

    def signature() -> str:
        return _live12_signature(
            output_root,
            workspace=workspace,
            package_root=package_root,
        )

    initial = signature()
    generated_cache = runtime_dir / "i123.iso"
    generated_cache.write_bytes(b"generated-v1")
    assert signature() == initial
    generated_cache.write_bytes(b"generated-v2")
    assert signature() == initial
    unrelated_method3 = output_root / "EuAl4 Parent.cif" / "官网" / "Method3" / "candidate"
    unrelated_method3.mkdir(parents=True)
    (unrelated_method3 / "subgroup.cif").write_text("method3-only\n", encoding="utf-8")
    assert signature() == initial

    signatures = [initial]
    config_path.write_text(config_path.read_text(encoding="utf-8") + "# changed\n", encoding="utf-8")
    signatures.append(signature())
    iso_path.write_bytes(b"iso-v2")
    signatures.append(signature())
    smodes_path.write_bytes(b"smodes-v2")
    signatures.append(signature())
    data_path.write_bytes(b"data-v2")
    signatures.append(signature())
    runtime["python"] = {"version": "runtime-v2"}
    signatures.append(signature())
    details_path.write_text("details-v2\n", encoding="utf-8")
    signatures.append(signature())
    assert len(set(signatures)) == len(signatures)
    assert not _live12_checkpoint_policy()["generated_wsl_i_star_iso_cache"][
        "included_in_global_signature"
    ]


class _Live12Subgroup:
    irrep_label = "GM1+"
    origin = (0, 0, 0)
    index = 0
    basis_scale = 1

    @classmethod
    def opd_line_body(cls) -> str:
        return (
            f"P1 (a) 1 P1, basis={{({cls.basis_scale},0,0),(0,1,0),(0,0,1)}}, "
            "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
        )


def _write_live12_reference(output_root: Path, case_name: str, method: str) -> None:
    folder = output_root / case_name / "官网" / method / "GM1+ P1"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "subgroup.cif").write_text(
        """# IR: GM1+
# P1 (a) 1 P1, basis={(1,0,0),(0,1,0),(0,0,1)}, origin=(0,0,0), s=1, i=1, k-active= (0,0,0)
data_test
""",
        encoding="utf-8",
    )
    (folder / "ISODISTORT complete modes details.html").write_text(
        ":dsp] test normfactor = 1\n",
        encoding="utf-8",
    )


def test_live12_checkpoint_marks_interruption_and_safely_reuses_identity_counts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workspace = tmp_path / "workspace"
    output_root = workspace / "output_compare"
    input_root = workspace / "experiment_data"
    input_root.mkdir(parents=True)
    for case_name in ("EuAl4 Parent.cif", "NdNiO2 own.cif"):
        (input_root / case_name).write_text("parent\n", encoding="utf-8")
        for method in ("Method1", "Method2"):
            _write_live12_reference(output_root, case_name, method)

    state = {"mode_calls": 0, "fail_after": 2}
    monkeypatch.setattr(_Live12Subgroup, "basis_scale", 1)

    class FakeIsoDistort:
        def __init__(self) -> None:
            self.mode_displacements: list[object] = []

        def load_structure(self, _path: Path) -> None:
            return None

        def search_method_1(self, **_kwargs: object) -> list[SimpleNamespace]:
            return [SimpleNamespace(subgroup=_Live12Subgroup())]

        def list_subgroups_at_kpoint(self, *_args: object, **_kwargs: object) -> list[_Live12Subgroup]:
            return [_Live12Subgroup()]

        def search_method_2(self, **_kwargs: object) -> None:
            if state["fail_after"] is not None and state["mode_calls"] >= state["fail_after"]:
                raise RuntimeError("simulated interruption")
            state["mode_calls"] += 1
            self.mode_displacements = [object()]

    monkeypatch.setattr(validation_module, "WORKSPACE", workspace)
    monkeypatch.setattr(isocore.api, "IsoDistort", FakeIsoDistort)
    checkpoint_path = workspace / "output" / "validation" / "checkpoint.json"
    report_path = workspace / "output" / "validation" / "report.json"

    with pytest.raises(RuntimeError, match="simulated interruption"):
        run_live_method12_audit(output_root, checkpoint_path, report_path)
    interrupted_checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    interrupted_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert interrupted_checkpoint["status"] == "incomplete_error"
    assert not interrupted_checkpoint["execution_complete"]
    assert interrupted_report["status"] == "incomplete_error"
    assert interrupted_report["progress"]["completed_groups"] == 2
    assert interrupted_report["last_invocation_elapsed_seconds"] >= 0

    state["fail_after"] = None
    completed = run_live_method12_audit(output_root, checkpoint_path, report_path)
    assert completed["status"] == "complete_passed"
    assert completed["execution_complete"]
    assert completed["validation_complete"]
    assert completed["passed"]
    assert completed["resumed_from_checkpoint"]
    assert sum(report["checkpoint_counts_reused"] for report in completed["reports"]) == 2
    assert sum(report["mode_counts_computed"] for report in completed["reports"]) == 2
    assert state["mode_calls"] == 4
    saved_report = json.loads(report_path.read_text(encoding="utf-8"))
    assert saved_report["status"] == "complete_passed"
    assert saved_report["progress"]["completed_groups"] == 4

    # The basis-independent official identity stays the same, but a changed
    # live basis must invalidate each saved mode count rather than mixing runs.
    _Live12Subgroup.basis_scale = 2
    changed_candidate = run_live_method12_audit(output_root, checkpoint_path, report_path)
    assert sum(report["checkpoint_counts_reused"] for report in changed_candidate["reports"]) == 0
    assert sum(report["mode_counts_computed"] for report in changed_candidate["reports"]) == 4
    assert state["mode_calls"] == 8
