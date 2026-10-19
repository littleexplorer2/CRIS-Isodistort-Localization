"""Tests for backend.models; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.export.cif_displacive_model import DisplaciveModelError
from backend.export.cif_displacive_model import build_displacive_cif_model
from backend.models.iso_mode_models import IsoPointOperationCatalog
from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import OrientedSiteCharacterTable
from backend.models.iso_mode_models import SiteIrrepCharacter
from backend.models.iso_mode_models import SiteOperation
from backend.models.iso_mode_models import build_site_irrep_projectors
from backend.models.iso_mode_models import validate_cartesian_subspaces
from backend.models.phase_path import DEFAULT_DISTORTION_TYPES
from backend.models.phase_path import normalize_distortion_types
from backend.models.phase_path import validate_modulation_count
from backend.utils.text_parser import parse_wyckoff_character_tables
from dataclasses import replace
from tests.support_backend_export import _anti_coupled_fixture
from tests.support_backend_method2 import ISODISTORT_ROOT
from tests.support_backend_method2 import _block
from tests.support_backend_method2 import _microscopic_column_provenance_QUERY_DIGEST
import numpy as np
import pytest


# Original coverage: test_cif_displacive_model.py.

@pytest.mark.parametrize(
    ("identity_changes", "message"),
    (
        ({"orbit_id": ""}, "physical orbit ID"),
        ({"k_coordinates": ("0", "g", "0")}, "exact rational k coordinates"),
    ),
)
def test_verified_mode_identity_requires_orbit_and_exact_rational_k(
    identity_changes: dict[str, object],
    message: str,
) -> None:
    frame, atoms, modes = _anti_coupled_fixture()
    identity = replace(modes[0].mode_identity, **identity_changes)
    invalid = replace(modes[0], mode_id=identity.stable_token, mode_identity=identity)

    assert message in invalid.source_binding_issue()
    with pytest.raises(DisplaciveModelError) as caught:
        build_displacive_cif_model(frame, atoms, ("A",), (invalid, modes[1]))
    assert caught.value.code == "unresolved_mode_provenance"

# Original coverage: test_distortion.py.

# --- from test_distortion.py ---

def test_default_distortion_modes():
    """默认同时启用 strain 和 displacive；旧别名保持相同物理类型。"""
    # 默认对齐官网 Types 面板：Strain 勾选 + Displacive 全物种勾选
    assert DEFAULT_DISTORTION_TYPES == ["strain", "displacive"]
    assert normalize_distortion_types(None) == ["strain", "displacive"]
    assert normalize_distortion_types(["strain", "displacement", "strain"]) == [
        "strain",
        "displacive",  # 旧名 displacement 自动映射为 displacive
    ]
    assert normalize_distortion_types("order") == ["occupational"]

# Original coverage: test_microscopic_column_provenance.py.

def test_active_arm_resolution_kind_is_validated_and_part_of_source_token():
    block = _block()
    exact = MicroscopicColumnProvenance.from_exact_column(
        query_digest=_microscopic_column_provenance_QUERY_DIGEST,
        query_order=0,
        query_irrep_label=block.global_irrep,
        source_parent_sg=62,
        source_k_coordinates=("0", "1/2", "0"),
        source_active_k_coordinates=(("0", "1/2", "0"),),
        source_active_k_resolution_kind="exact_column_support",
        source_subgroup_space_group_number=1,
        source_subgroup_basis=np.eye(3, dtype=int).tolist(),
        source_subgroup_origin=(0, 0, 0),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        block=block,
        column_index=0,
    )
    envelope = replace(
        exact,
        source_active_k_resolution_kind="compatible_envelope",
    )

    assert exact.exact_source_token != envelope.exact_source_token
    with pytest.raises(ValueError, match="recognized resolution kind"):
        replace(exact, source_active_k_resolution_kind="guessed")
    with pytest.raises(ValueError, match="requires resolved arms"):
        replace(exact, source_active_k_coordinates=())

# Original coverage: test_session_result_contract.py.

@pytest.mark.parametrize("value", [True, False, 1.5, -1, "1", None])
def test_modulation_count_does_not_coerce_a_different_scientific_request(value):
    with pytest.raises(ValueError, match="nmod"):
        validate_modulation_count(value)

def test_modulation_count_reports_unsupported_domain_instead_of_clamping():
    for value in range(4):
        assert validate_modulation_count(value) == value
    with pytest.raises(NotImplementedError, match="nmod=0..3"):
        validate_modulation_count(4)

# Original coverage: test_site_irreps.py.

def _character_table(row: dict) -> OrientedSiteCharacterTable:
    return OrientedSiteCharacterTable(
        wyckoff_letter=row["wyckoff_letter"],
        point_group=row["point_group"],
        operations=tuple(SiteOperation(**item) for item in row["operations"]),
        irreps=tuple(SiteIrrepCharacter(**item) for item in row["irreps"]),
    )

def _matrix(projectors, label: str) -> np.ndarray:
    item = next(projector for projector in projectors if projector.label == label)
    return np.asarray(item.matrix, dtype=float)

def test_oriented_c2v_and_d2h_polar_vector_projectors():
    text = (
        "c D2h (E|0,0,0) Ag 1, B3g 1, B1g 1, B2g 1, Au 1, B3u 1, "
        "B1u 1, B2u 1, (C2x|0,1,0) Ag 1, B3g 1, B1g -1, B2g -1, "
        "Au 1, B3u 1, B1u -1, B2u -1, (C2y|0,0,0) Ag 1, B3g -1, "
        "B1g -1, B2g\n"
        "             1, Au 1, B3u -1, B1u -1, B2u 1, (C2z|0,1,0) "
        "Ag 1, B3g -1, B1g 1, B2g -1, Au 1, B3u -1, B1u 1, B2u -1, "
        "(I|0,1,0) Ag 1, B3g 1, B1g 1, B2g 1, Au -1, B3u -1, B1u -1, "
        "B2u -1,\n"
        "             (SGx|0,0,0) Ag 1, B3g 1, B1g -1, B2g -1, Au -1, "
        "B3u -1, B1u 1, B2u 1, (SGy|0,1,0) Ag 1, B3g -1, B1g -1, "
        "B2g 1, Au -1, B3u 1, B1u 1, B2u -1, (SGz|0,0,0) Ag 1, "
        "B3g -1, B1g 1,\n"
        "             B2g -1, Au -1, B3u 1, B1u -1, B2u 1\n"
        "g C2v (E|0,0,0) A1 1, B2 1, A2 1, B1 1, (C2z|0,1,0) A1 1, "
        "B2 -1, A2 1, B1 -1, (SGx|0,0,0) A1 1, B2 1, A2 -1, B1 -1, "
        "(SGy|0,1,0) A1 1, B2 -1, A2 -1, B1 1\n*\n"
    )
    tables = parse_wyckoff_character_tables(text)
    assert [(row["wyckoff_letter"], row["point_group"]) for row in tables] == [
        ("c", "D2h"), ("g", "C2v")
    ]
    catalog = IsoPointOperationCatalog.from_data_space(
        ISODISTORT_ROOT / "resources" / "isobyu" / "data_space.txt"
    )
    d2h = build_site_irrep_projectors(_character_table(tables[0]), catalog)
    c2v = build_site_irrep_projectors(_character_table(tables[1]), catalog)
    assert _matrix(c2v, "A1") == pytest.approx(np.diag([0.0, 0.0, 1.0]))
    assert _matrix(c2v, "B1") == pytest.approx(np.diag([1.0, 0.0, 0.0]))
    assert _matrix(c2v, "B2") == pytest.approx(np.diag([0.0, 1.0, 0.0]))
    assert _matrix(c2v, "A2") == pytest.approx(np.zeros((3, 3)))
    assert _matrix(d2h, "B3u") == pytest.approx(np.diag([1.0, 0.0, 0.0]))
    assert _matrix(d2h, "B2u") == pytest.approx(np.diag([0.0, 1.0, 0.0]))
    assert _matrix(d2h, "B1u") == pytest.approx(np.diag([0.0, 0.0, 1.0]))

def test_cartesian_subspace_validation_accepts_basis_rotation_and_rejects_mismatch():
    lattice = np.asarray([[2.0, 0.3, 0.0], [0.0, 3.0, 0.2], [0.0, 0.0, 4.0]])
    x = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    y = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])
    rotated = [x + y, x - y]
    matched = validate_cartesian_subspaces([x, y], rotated, lattice)
    assert matched.matched
    assert matched.reference_rank == matched.candidate_rank == 2
    mismatch = validate_cartesian_subspaces([x, y], [x], lattice)
    assert not mismatch.matched
    assert mismatch.reason == "rank_mismatch"
