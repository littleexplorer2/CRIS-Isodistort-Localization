"""Tests for support_backend_export; migrated without changing scientific expectations."""
from __future__ import annotations

from application.session.core_api import IsoDistort
from application.session.core_api import IsoDistort as _formats_IsoDistort
from backend.export.cif_displacive_model import ChildAtom
from backend.export.cif_displacive_model import ExactChildFrame
from backend.export.cif_displacive_model import ExactSeitzOperation
from backend.export.cif_displacive_model import ModeScaleProvenance
from backend.export.cif_displacive_model import RawModeColumn
from backend.export.distortion_formats import SubgroupExportSpec
from backend.models.iso_mode_models import InvariantDirection
from backend.models.iso_mode_models import MacroscopicTensorBasis
from backend.models.iso_mode_models import MacroscopicTensorBlock
from backend.models.iso_mode_models import MicroscopicColumnProvenance
from backend.models.iso_mode_models import ModeIdentity
from backend.models.iso_mode_models import SymbolicTensorComponent
from backend.models.iso_mode_models import microscopic_atom_order_id
from backend.utils import get_config
from backend.wrappers import SubgroupInfo
from backend.wrappers import SubgroupInfo as _formats_SubgroupInfo
from fractions import Fraction
from pymatgen.core import Structure
from tests.support_paths import TESTS_ROOT
import hashlib
import math
import numpy as np


# Original coverage: test_cif_displacive_model.py.

_cif_displacive_model_I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

_cif_displacive_model_ZERO = (0, 0, 0)

ANTI_DIAGONAL_MIRROR = ((0, -1, 0), (-1, 0, 0), (0, 0, 1))

# This skew metric has Gxx=Gyy and is invariant under ANTI_DIAGONAL_MIRROR.
SKEW_LATTICE = ((2.0, 0.0, 0.0), (0.5, math.sqrt(15.0) / 2.0, 0.0), (0.0, 0.0, 4.0))

def _op(rotation=_cif_displacive_model_I3, translation=_cif_displacive_model_ZERO) -> ExactSeitzOperation:
    return ExactSeitzOperation.from_values(rotation, translation)

def _canonical_source() -> ModeScaleProvenance:
    return ModeScaleProvenance(
        source="iso_microscopic_display_distortion",
        convention="authoritative_fractional_direction_arbitrary_positive_scale",
        direction_resolved=True,
        sign_resolved=True,
    )

def _mode(
    display_label: str,
    column: int,
    frame_id: str,
    atoms: tuple[ChildAtom, ...],
    vectors,
    *,
    provenance: ModeScaleProvenance | None = None,
) -> RawModeColumn:
    atom_ids = tuple(atom.atom_id for atom in atoms)
    global_irrep = "TEST1+"
    site_irrep = "A"
    identity = ModeIdentity(
        parent_sg=1,
        global_irrep=global_irrep,
        k_coordinates=("0", "0", "0"),
        wyckoff_letter="a",
        orbit_id=atom_ids[0] if atom_ids else "",
        site_irrep=site_irrep,
        component_index=column,
        component_label=display_label,
        source="iso_microscopic",
        status="verified",
    )
    row_digest = hashlib.sha256(
        repr((frame_id, atom_ids)).encode("utf-8")
    ).hexdigest()
    vector_digest = hashlib.sha256(
        repr((display_label, vectors)).encode("utf-8")
    ).hexdigest()
    microscopic_provenance = MicroscopicColumnProvenance(
        query_digest=hashlib.sha256(b"synthetic exact query").hexdigest(),
        query_order=0,
        query_irrep_label=global_irrep,
        source_parent_sg=1,
        source_k_coordinates=("0", "0", "0"),
        source_active_k_coordinates=(),
        source_subgroup_space_group_number=1,
        source_subgroup_basis=tuple(
            tuple(str(value) for value in row) for row in _cif_displacive_model_I3
        ),
        source_subgroup_origin=tuple(str(value) for value in _cif_displacive_model_ZERO),
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
        source_block_order=0,
        source_global_irrep=global_irrep,
        source_wyckoff_letter="a",
        source_site_irrep=site_irrep,
        source_direction_symbol="VECTOR,A",
        source_column_index=column,
        source_block_column_count=64,
        source_row_count=len(atom_ids),
        exact_row_digest=row_digest,
        exact_vector_digest=vector_digest,
        source_float_column_digest=vector_digest,
        source_frame_id=frame_id,
        source_atom_order_id=microscopic_atom_order_id(atom_ids),
        unresolved_fields=(),
        unresolved_diagnostics=(),
    )
    return RawModeColumn.from_values(
        identity.stable_token,
        column,
        frame_id,
        atom_ids,
        vectors,
        provenance or _canonical_source(),
        mode_identity=identity,
        microscopic_provenance=microscopic_provenance,
        label=display_label,
    )

def _anti_coupled_fixture(frame_id: str = "anti"):
    operations = (_op(), _op(ANTI_DIAGONAL_MIRROR))
    atoms = (ChildAtom.from_values("A", _cif_displacive_model_ZERO, "X"),)
    modes = (
        _mode("anti_xy", 0, frame_id, atoms, ((1.0, -1.0, 0.0),)),
        _mode("z", 1, frame_id, atoms, ((0.0, 0.0, 1.0),)),
    )
    frame = ExactChildFrame.from_values(frame_id, operations, SKEW_LATTICE)
    return frame, atoms, modes

# Original coverage: test_formats.py.

_formats_REPO = TESTS_ROOT.parent.parent

# --- from test_distortion_export.py ---

def _sg(irrep: str, opd: str, index: int = 0, number: int = 123,
        symbol: str = "P4/mmm") -> _formats_SubgroupInfo:
    return _formats_SubgroupInfo(
        index=index,
        space_group_number=number,
        space_group_symbol=symbol,
        subgroup_index=12,
        size=12,
        opd_symbol=opd,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 6]],
        origin=[0, 0, 0],
        k_point_label="LD",
        irrep_label=irrep,
        k_parameters=["1/6"],
    )

def _cubic() -> Structure:
    return Structure(
        [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
        ["Eu", "Al"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )

def _minimal_export_api(parent: Structure) -> _formats_IsoDistort:
    iso = object.__new__(_formats_IsoDistort)
    iso.structure = parent
    iso.symmetry_info = None
    iso.distortion_types = []
    iso._selected_subgroup = None
    iso.phase_path = None
    iso.distortion_modes = []
    iso.mode_displacements = {}
    iso.mode_occupancies = {}
    iso.mode_displacements_sc = {}
    iso._mode_label_overrides = {}
    iso.distorted_structure = None
    iso.number_of_independent_modulations = 0
    iso._spec_for_subgroup = lambda subgroup, **kwargs: SubgroupExportSpec(
        subgroup=subgroup,
        structure=parent,
        folder_name=kwargs["folder_name"],
    )
    return iso

# Original coverage: test_strain_export.py.

def _strain_export_subgroup(number: int, symbol: str, basis: list[list[float]] | None = None) -> SubgroupInfo:
    return SubgroupInfo(
        index=0,
        parent_sg=number,
        space_group_number=number,
        space_group_symbol=symbol,
        subgroup_index=1,
        size=1,
        opd_symbol="P1",
        opd_dir_raw="(a)",
        basis_vectors=basis or np.eye(3).tolist(),
        origin=[0.0, 0.0, 0.0],
        k_point_label="GM",
        k_coordinates=["0", "0", "0"],
        irrep_label="GM_TEST",
    )

def _bare_api(structure: Structure, distortion_types: list[str]) -> IsoDistort:
    """Construct only the pure export context; no WSL backend is needed."""

    api = object.__new__(IsoDistort)
    api.cfg = get_config()
    api.structure = structure
    api.distortion_types = distortion_types
    api._strain_representation_cache = None
    return api

def _macro_component(expression: str, values: list[int]) -> SymbolicTensorComponent:
    return SymbolicTensorComponent(
        expression_raw=expression,
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        coefficients=tuple(Fraction(value) for value in values),
        quality="exact_rational",
    )

def _tetragonal_macro_basis() -> MacroscopicTensorBasis:
    return MacroscopicTensorBasis(
        rank_signature="[12]",
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        blocks=(
            MacroscopicTensorBlock(
                "GM1+", (_macro_component("xx+yy", [1, 1, 0, 0, 0, 0]),), 1, 1
            ),
            MacroscopicTensorBlock(
                "GM1+", (_macro_component("zz", [0, 0, 1, 0, 0, 0]),), 2, 2
            ),
            MacroscopicTensorBlock(
                "GM2+", (_macro_component("xx-yy", [1, -1, 0, 0, 0, 0]),), 3
            ),
            MacroscopicTensorBlock(
                "GM4+", (_macro_component("xy", [0, 0, 0, 0, 0, 1]),), 4
            ),
            MacroscopicTensorBlock(
                "GM5+",
                (
                    _macro_component("yz", [0, 0, 0, 1, 0, 0]),
                    _macro_component("-xz", [0, 0, 0, 0, -1, 0]),
                ),
                5,
            ),
        ),
        status="verified",
    )

class _TetragonalMacroIso:
    def get_macroscopic_tensor_basis(self, _parent_sg: int) -> MacroscopicTensorBasis:
        return _tetragonal_macro_basis()

    def list_invariant_directions(
        self, _parent_sg: int, subgroup: SubgroupInfo
    ) -> list[InvariantDirection]:
        if int(subgroup.space_group_number) not in {1, 2}:
            raise AssertionError("test ISO stub is only valid for P1/P-1 embeddings")
        return [
            InvariantDirection(label, direction, 2, "P-1", subgroup.size)
            for label, direction in (
                ("GM1+", "(a)"),
                ("GM2+", "(a)"),
                ("GM4+", "(a)"),
                ("GM5+", "(a,b)"),
            )
        ]
