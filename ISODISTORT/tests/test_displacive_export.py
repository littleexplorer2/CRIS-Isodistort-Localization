from __future__ import annotations

import inspect
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction

import numpy as np
import pytest
import spglib
from pymatgen.core import IStructure, Lattice, Structure

from backend.models.iso_mode_models import (
    MicroscopicColumnProvenance,
    MicroscopicVectorBlock,
    MicroscopicVectorRow,
    ModeIdentity,
    microscopic_atom_order_id,
)
from backend.wrappers import SubgroupInfo
from features.export.cif_displacive_model import (
    ChildAtom,
    DisplaciveModelError,
    ExactChildFrame,
    ExactSeitzOperation,
    ModeScaleProvenance,
    RawModeColumn,
    build_displacive_cif_model,
    site_fixed_coordinate_model,
)
from features.export.displacive_export import (
    DisplaciveExportData,
    DisplaciveSubgroupIdentity,
    ExactParentChildEmbedding,
    ParentChildSiteMapping,
    ParentOrbitType,
    _periodic_cartesian_residual,
    _validate_emitted_space_group,
)
from features.export.distortion_formats import (
    SubgroupExportSpec,
    render_cif,
    render_complete_modes,
    render_isoviz,
    render_topas,
)
from features.export.isodistort_cif import (
    _apply_origin_choice,
    _coordinate_formula,
    _setting,
    _validate_displacive_writer_frame,
    _validated_subgroup_sites,
)
from features.export.isodistort_isoviz import render_isodistort_isoviz

I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ZERO = (0, 0, 0)


def _site_mapping(
    child_atom_id: str,
    parent_site_index: int,
    parent_orbit_id: str,
    translation=ZERO,
) -> ParentChildSiteMapping:
    return ParentChildSiteMapping(
        child_atom_id=child_atom_id,
        parent_site_index=parent_site_index,
        parent_orbit_id=parent_orbit_id,
        parent_cell_translation=translation,
    )


def _subgroup_identity(
    embedding: ExactParentChildEmbedding,
    *,
    target_space_group_number: int = 1,
    target_hall_number: int | None = None,
    irrep_label: str = "GM1+",
    opd_symbol: str = "P1",
    primary_direction_selector: str | None = "VECTOR,A",
    context_kind: str = "single_irrep",
) -> DisplaciveSubgroupIdentity:
    return DisplaciveSubgroupIdentity(
        embedding,
        target_space_group_number=target_space_group_number,
        target_hall_number=target_hall_number,
        irrep_label=(None if context_kind == "exact_fixed_space" else irrep_label),
        opd_symbol=(None if context_kind == "exact_fixed_space" else opd_symbol),
        primary_direction_selector=primary_direction_selector,
        context_kind=context_kind,
    )


def _resolved_mode(
    *,
    axis: int,
    frame_id: str,
    atom_ids: tuple[str, ...],
    displacements,
    global_irrep: str,
    parent_sg: int,
    wyckoff_letter: str,
    orbit_id: str | None = None,
    column_index: int | None = None,
    source_child_sg: int = 1,
    source_basis=I3,
    source_origin=ZERO,
    source_subgroup_irrep: str = "GM1+",
    source_subgroup_opd: str = "P1",
    context_kind: str = "single_irrep",
) -> RawModeColumn:
    source_order = axis if column_index is None else column_index
    site_irrep = f"A{axis + 1}"
    block = MicroscopicVectorBlock(
        global_irrep=global_irrep,
        wyckoff_letter=wyckoff_letter,
        site_irrep=site_irrep,
        rows=(
            MicroscopicVectorRow(
                point_raw=("0", "0", "0"),
                displacements=((Fraction(1), Fraction(0), Fraction(0)),),
            ),
        ),
        source_order=source_order,
        direction_symbol="VECTOR,A",
    )
    provenance = MicroscopicColumnProvenance.from_exact_column(
        query_digest=f"{source_order:064x}",
        query_order=source_order,
        query_irrep_label=global_irrep,
        source_parent_sg=parent_sg,
        source_k_coordinates=("0", "0", "0"),
        source_subgroup_space_group_number=source_child_sg,
        source_subgroup_basis=source_basis,
        source_subgroup_origin=source_origin,
        source_subgroup_irrep_label=(
            None if context_kind == "exact_fixed_space"
            else source_subgroup_irrep
        ),
        source_subgroup_opd_symbol=(
            None if context_kind == "exact_fixed_space"
            else source_subgroup_opd
        ),
        source_subgroup_primary_direction_selector=(
            None if context_kind == "exact_fixed_space" else "VECTOR,A"
        ),
        source_subgroup_context_kind=context_kind,
        block=block,
        column_index=0,
        source_frame_id=frame_id,
        source_atom_order_id=microscopic_atom_order_id(atom_ids),
    )
    identity = ModeIdentity(
        parent_sg=parent_sg,
        global_irrep=global_irrep,
        k_coordinates=("0", "0", "0"),
        wyckoff_letter=wyckoff_letter,
        orbit_id=orbit_id or f"orbit-{wyckoff_letter}",
        site_irrep=site_irrep,
        component_index=0,
        component_label=f"c{axis + 1}",
        source="iso_microscopic",
        status="verified",
    )
    return RawModeColumn.from_values(
        identity.stable_token,
        source_order,
        frame_id,
        atom_ids,
        displacements,
        ModeScaleProvenance(
            source="iso_microscopic_display_distortion",
            convention="unmixed_fractional_source_column",
            direction_resolved=True,
            sign_resolved=True,
        ),
        mode_identity=identity,
        microscopic_provenance=provenance,
        label=f"mode-{source_order}",
    )


def _export_data(
    *,
    final_x: float | None = None,
    global_irreps: tuple[str, ...] = ("GM1+", "GM1+", "GM1+"),
    context_kind: str = "single_irrep",
) -> DisplaciveExportData:
    operation = ExactSeitzOperation.from_values(I3, ZERO)
    frame = ExactChildFrame.from_values(
        "p1_exact",
        (operation,),
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atom = ChildAtom.from_values(
        "a0", (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)), "Fe"
    )
    modes = tuple(
        _resolved_mode(
            axis=axis,
            frame_id=frame.frame_id,
            atom_ids=(atom.atom_id,),
            displacements=(tuple(float(value) for value in np.eye(3)[axis]),),
            global_irrep=global_irreps[axis],
            parent_sg=1,
            wyckoff_letter="a",
            context_kind=context_kind,
        )
        for axis in range(3)
    )
    model = build_displacive_cif_model(
        frame,
        (atom,),
        (atom.atom_id,),
        modes,
        amplitudes=(0.4, 0.0, 0.0),
    )
    reference = Structure(
        Lattice.cubic(4.0),
        ["Fe"],
        [[0.25, 0.5, 0.75]],
        labels=["Fe1"],
    )
    final = None
    if final_x is not None:
        final = Structure(
            reference.lattice,
            ["Fe"],
            [[final_x, 0.5, 0.75]],
            labels=["Fe1"],
        )
    parent = Structure.from_sites(reference)
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        I3,
        ZERO,
        (_site_mapping("a0", 0, "orbit-a"),),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(
            embedding,
            primary_direction_selector=(
                None if context_kind == "exact_fixed_space" else "VECTOR,A"
            ),
            context_kind=context_kind,
        ),
        reference_structure=reference,
        final_structure=final,
        representative_labels={"a0": "Fe1_1"},
        parent_orbit_types={
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
        },
    )


def _spec(data: DisplaciveExportData) -> SubgroupExportSpec:
    exact_fixed_space = (
        data.subgroup_identity.context_kind == "exact_fixed_space"
    )
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=1,
        space_group_symbol="P1",
        size=1,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
        irrep_label="" if exact_fixed_space else "GM1+",
        opd_symbol="" if exact_fixed_space else "P1",
        opd_dir_raw="" if exact_fixed_space else "(a)",
        parent_sg=1,
    )
    if exact_fixed_space:
        subgroup._method3_route_resolution = "exact_fixed_space"
    subgroup._displacive_embedding_id = data.embedding.embedding_id
    return SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )


def _two_same_species_orbit_export_data() -> DisplaciveExportData:
    operation = ExactSeitzOperation.from_values(I3, ZERO)
    frame = ExactChildFrame.from_values(
        "p1_two_orbits",
        (operation,),
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atoms = (
        ChildAtom.from_values("a0", (Fraction(1, 4),) * 3, "Fe"),
        ChildAtom.from_values("b0", (Fraction(3, 4),) * 3, "Fe"),
    )
    atom_ids = tuple(atom.atom_id for atom in atoms)
    modes = []
    for orbit_position, (orbit_id, atom_position) in enumerate(
        (("orbit-a", 0), ("orbit-b", 1))
    ):
        for axis in range(3):
            rows = np.zeros((2, 3), dtype=float)
            rows[atom_position, axis] = 1.0
            modes.append(
                _resolved_mode(
                    axis=axis,
                    frame_id=frame.frame_id,
                    atom_ids=atom_ids,
                    displacements=rows,
                    global_irrep="GM1+",
                    parent_sg=1,
                    wyckoff_letter="a",
                    orbit_id=orbit_id,
                    column_index=3 * orbit_position + axis,
                )
            )
    model = build_displacive_cif_model(
        frame,
        atoms,
        atom_ids,
        tuple(modes),
        amplitudes=(0.0,) * 6,
    )
    reference = Structure(
        Lattice.cubic(4.0),
        ["Fe", "Fe"],
        [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
        labels=["Fe1", "Fe2"],
    )
    parent = Structure.from_sites(reference)
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        I3,
        ZERO,
        (
            _site_mapping("a0", 0, "orbit-a"),
            _site_mapping("b0", 1, "orbit-b"),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels={"a0": "Fe1_1", "b0": "Fe2_1"},
        parent_orbit_types={
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
            "orbit-b": ParentOrbitType(2, "Fe2", ("b0",)),
        },
    )


def _index_two_export_data() -> DisplaciveExportData:
    basis = ((2, 0, 0), (0, 1, 0), (0, 0, 1))
    operation = ExactSeitzOperation.from_values(I3, ZERO)
    frame = ExactChildFrame.from_values(
        "p1_index_two",
        (operation,),
        ((8.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atoms = (
        ChildAtom.from_values(
            "a0", (Fraction(1, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
        ChildAtom.from_values(
            "a1", (Fraction(5, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
    )
    atom_ids = tuple(atom.atom_id for atom in atoms)
    modes = []
    for atom_position, global_irrep in enumerate(("GM1+", "GM2+")):
        for axis in range(3):
            rows = np.zeros((2, 3), dtype=float)
            rows[atom_position, axis] = 1.0
            modes.append(
                _resolved_mode(
                    axis=axis,
                    frame_id=frame.frame_id,
                    atom_ids=atom_ids,
                    displacements=rows,
                    global_irrep=global_irrep,
                    parent_sg=1,
                    wyckoff_letter="a",
                    orbit_id="orbit-a",
                    column_index=3 * atom_position + axis,
                    source_basis=basis,
                )
            )
    model = build_displacive_cif_model(
        frame,
        atoms,
        atom_ids,
        tuple(modes),
        amplitudes=(0.4, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
    parent = Structure(
        Lattice.cubic(4.0),
        ["Fe"],
        [[0.25, 0.25, 0.25]],
        labels=["Fe1"],
    )
    reference = Structure(
        Lattice.orthorhombic(8.0, 4.0, 4.0),
        ["Fe", "Fe"],
        [[0.125, 0.25, 0.25], [0.625, 0.25, 0.25]],
        labels=["Fe1", "Fe1"],
    )
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        basis,
        ZERO,
        (
            _site_mapping("a0", 0, "orbit-a", (0, 0, 0)),
            _site_mapping("a1", 0, "orbit-a", (1, 0, 0)),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels={"a0": "Fe1_1", "a1": "Fe1_2"},
        parent_orbit_types={
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0", "a1")),
        },
    )


def _subgroup_for_data(data: DisplaciveExportData) -> SubgroupInfo:
    identity = data.subgroup_identity
    exact_fixed_space = identity.context_kind == "exact_fixed_space"
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=identity.target_space_group_number,
        space_group_symbol="P1",
        size=identity.size,
        basis_vectors=[
            [value for value in row] for row in identity.basis
        ],
        origin=[value for value in identity.origin],
        irrep_label=identity.irrep_label or "",
        opd_symbol=identity.opd_symbol or "",
        opd_dir_raw="" if exact_fixed_space else "(a)",
        parent_sg=identity.parent_space_group_number,
    )
    if exact_fixed_space:
        subgroup._method3_route_resolution = "exact_fixed_space"
    subgroup._displacive_embedding_id = identity.embedding_id
    return subgroup


def test_shared_export_contract_separates_reference_and_final_coordinates() -> None:
    data = _export_data()
    spec = _spec(data)

    assert spec.structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert spec.cif_structure is not None
    assert spec.cif_structure.frac_coords[0, 0] == pytest.approx(0.35)
    assert tuple(spec.mode_displacements_sc or {}) == data.mode_keys
    assert spec.amplitudes == dict(zip(data.mode_keys, (0.4, 0.0, 0.0), strict=True))
    assert spec.mode_global_irreps == {
        key: "GM1+" for key in data.mode_keys
    }
    assert spec.mode_keys == data.mode_keys
    assert spec.atom_ids == ("a0",)
    assert spec.parent_sg == 1
    assert spec.parent_symbol == "P1"
    assert spec.parent_orbit_types is data.parent_orbit_types
    rows = data.modes()
    assert {row.parent_orbit_id for row in rows} == {"orbit-a"}
    assert {row.parent_type_index for row in rows} == {1}
    assert {row.parent_type_label for row in rows} == {"Fe1"}

    with pytest.raises(DisplaciveModelError) as caught:
        _export_data(final_x=0.25)
    assert caught.value.code == "final_displacement_mismatch"


def test_parent_child_mapping_proves_complete_index_two_translation_quotient() -> None:
    data = _index_two_export_data()
    subgroup = _subgroup_for_data(data)

    assert data.embedding.conventional_translation_index == 2
    assert data.subgroup_identity.size == 2
    assert tuple(
        mapping.parent_cell_translation
        for mapping in data.embedding.atom_mappings
    ) == ((0, 0, 0), (1, 0, 0))
    with pytest.raises(FrozenInstanceError):
        data.embedding.atom_mappings[0].parent_site_index = 1  # type: ignore[misc]
    data.validate_subgroup(subgroup)

    spec = SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    complete = render_complete_modes(spec)
    assert "As=  0.400000" in complete
    assert "Ap=  0.282843" in complete


def test_parent_child_mapping_rejects_duplicate_translation_coset() -> None:
    source = _index_two_export_data()
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        ZERO,
        (
            _site_mapping("a0", 0, "orbit-a", (0, 0, 0)),
            _site_mapping("a1", 0, "orbit-a", (0, 0, 0)),
        ),
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=embedding,
            subgroup_identity=_subgroup_identity(embedding),
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "duplicate_parent_child_translation_coset"


def test_parent_child_mapping_rejects_identity_basis_zero_origin_global_shift() -> None:
    source = _export_data()
    shifted_atom = ChildAtom.from_values(
        "a0", (Fraction(3, 8), Fraction(1, 2), Fraction(3, 4)), "Fe"
    )
    shifted_model = build_displacive_cif_model(
        source.model.frame,
        (shifted_atom,),
        (shifted_atom.atom_id,),
        source.model.canonical_modes,
        amplitudes=source.model.amplitudes,
    )
    shifted_reference = Structure(
        source.reference_structure.lattice,
        ["Fe"],
        [[0.375, 0.5, 0.75]],
        labels=["Fe1"],
    )
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        I3,
        ZERO,
        source.embedding.atom_mappings,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=shifted_model,
            embedding=embedding,
            subgroup_identity=_subgroup_identity(embedding),
            reference_structure=shifted_reference,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "parent_child_coordinate_mismatch"


def test_parent_child_mapping_rejects_chemical_or_orbit_substitution() -> None:
    source = _export_data()
    wrong_chemistry_parent = Structure(
        source.parent_structure.lattice,
        ["Co"],
        source.parent_structure.frac_coords,
        labels=["Co1"],
    )
    chemical_embedding = ExactParentChildEmbedding(
        wrong_chemistry_parent,
        1,
        I3,
        ZERO,
        source.embedding.atom_mappings,
    )
    with pytest.raises(DisplaciveModelError) as chemical:
        DisplaciveExportData(
            model=source.model,
            embedding=chemical_embedding,
            subgroup_identity=_subgroup_identity(chemical_embedding),
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )
    assert chemical.value.code == "parent_child_chemistry_mismatch"

    orbit_embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        I3,
        ZERO,
        (_site_mapping("a0", 0, "another-orbit"),),
    )
    with pytest.raises(DisplaciveModelError) as orbit:
        DisplaciveExportData(
            model=source.model,
            embedding=orbit_embedding,
            subgroup_identity=_subgroup_identity(orbit_embedding),
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )
    assert orbit.value.code == "parent_child_orbit_mapping_mismatch"


def test_subgroup_identity_derives_size_and_rejects_every_mutable_identity_field() -> None:
    data = _index_two_export_data()
    assert "size" not in inspect.signature(DisplaciveSubgroupIdentity).parameters

    mutations = (
        ("irrep_label", "GM2+", "subgroup_identity_mismatch"),
        ("opd_symbol", "P2", "subgroup_identity_mismatch"),
        ("opd_dir_raw", "(0,a)", "subgroup_identity_mismatch"),
        ("size", 1, "subgroup_identity_mismatch"),
        ("parent_sg", 2, "subgroup_embedding_mismatch"),
        ("space_group_number", 2, "subgroup_identity_mismatch"),
    )
    for field, value, expected_code in mutations:
        subgroup = _subgroup_for_data(data)
        setattr(subgroup, field, value)
        with pytest.raises(DisplaciveModelError) as caught:
            data.validate_subgroup(subgroup)
        assert caught.value.code == expected_code

    subgroup = _subgroup_for_data(data)
    subgroup.basis_vectors[0][0] = 3
    with pytest.raises(DisplaciveModelError) as basis:
        data.validate_subgroup(subgroup)
    assert basis.value.code == "subgroup_embedding_mismatch"

    subgroup = _subgroup_for_data(data)
    subgroup.origin[0] = 1
    with pytest.raises(DisplaciveModelError) as origin:
        data.validate_subgroup(subgroup)
    assert origin.value.code == "subgroup_embedding_mismatch"

    subgroup = _subgroup_for_data(data)
    subgroup._displacive_embedding_id = "forged-embedding"
    with pytest.raises(DisplaciveModelError) as embedding_id:
        data.validate_subgroup(subgroup)
    assert embedding_id.value.code == "subgroup_embedding_mismatch"


def test_export_query_context_distinguishes_exact_embedding_and_primary_direction() -> None:
    source = _export_data()

    def _model_with_provenance_change(**changes):
        modes = tuple(
            replace(
                mode,
                microscopic_provenance=replace(
                    mode.microscopic_provenance,
                    **changes,
                ),
            )
            for mode in source.model.canonical_modes
        )
        return build_displacive_cif_model(
            source.model.frame,
            source.model.atoms,
            tuple(orbit.representative_atom_id for orbit in source.model.orbits),
            modes,
            amplitudes=source.model.amplitudes,
        )

    context_changes = (
        {"source_subgroup_space_group_number": 2},
        {"source_subgroup_basis": (("2", "0", "0"), ("0", "1", "0"), ("0", "0", "1"))},
        {"source_subgroup_origin": ("1", "0", "0")},
        {"source_subgroup_irrep_label": "GM2+"},
        {"source_subgroup_opd_symbol": "P2"},
        {"source_subgroup_primary_direction_selector": "VECTOR,0,A"},
    )
    for changes in context_changes:
        with pytest.raises(DisplaciveModelError) as caught:
            DisplaciveExportData(
                model=_model_with_provenance_change(**changes),
                embedding=source.embedding,
                subgroup_identity=source.subgroup_identity,
                reference_structure=source.reference_structure,
                final_structure=source.final_structure,
                representative_labels=source.representative_labels,
                parent_orbit_types=source.parent_orbit_types,
            )
        assert caught.value.code == "microscopic_query_context_mismatch"


def test_export_query_context_does_not_require_a_primary_displacement_block() -> None:
    data = _export_data(global_irreps=("GM2+", "GM2+", "GM2+"))
    assert data.subgroup_identity.primary_direction_selector == "VECTOR,A"
    assert {row.global_irrep_label for row in data.modes()} == {"GM2+"}


def test_export_contract_has_no_parallel_mode_or_atom_identity_inputs() -> None:
    parameters = inspect.signature(DisplaciveExportData).parameters
    assert not {
        "mode_keys",
        "mode_labels",
        "mode_global_irreps",
        "atom_ids",
    }.intersection(parameters)


def test_parent_orbit_type_contract_rejects_unmapped_mode_orbit() -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types={
                "another-orbit": ParentOrbitType(1, "Fe1", ("a0",)),
            },
        )

    assert caught.value.code == "missing_parent_orbit_type"


@pytest.mark.parametrize(
    "parent_orbit_types",
    [
        {
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
            "orbit-b": ParentOrbitType(1, "Fe2", ("other",)),
        },
        {
            "orbit-a": ParentOrbitType(1, "Fe1", ("a0",)),
            "orbit-b": ParentOrbitType(2, "Fe1", ("other",)),
        },
    ],
)
def test_parent_orbit_type_contract_requires_unique_indices_and_labels(
    parent_orbit_types: dict[str, ParentOrbitType],
) -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=parent_orbit_types,
        )

    assert caught.value.code == "invalid_parent_orbit_type"


def test_parent_orbit_type_contract_requires_exact_child_atom_partition() -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types={
                "orbit-a": ParentOrbitType(1, "Fe1", ("unknown",)),
            },
        )

    assert caught.value.code == "parent_orbit_atom_mismatch"


def test_mode_source_rejects_empty_global_irrep_before_export_contract() -> None:
    with pytest.raises(ValueError, match="query irrep cannot be empty"):
        _export_data(global_irreps=("GM1+", "", "GM1+"))


def test_production_writers_fail_closed_without_verified_displacive_data() -> None:
    reference = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
    spec = SubgroupExportSpec(
        subgroup=SubgroupInfo(index=0, space_group_number=1, space_group_symbol="P1"),
        structure=reference,
        mode_displacements_sc={"anonymous": np.array([[1.0, 0.0, 0.0]])},
        require_verified_displacive_data=True,
    )

    for writer in (
        lambda: render_cif(reference, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    ):
        with pytest.raises(ValueError, match="exact emitted-setting frame"):
            writer()


def test_four_writers_share_validated_modes_without_double_counting() -> None:
    spec = _spec(_export_data())

    cif = render_cif(spec.cif_structure, spec)
    complete = render_complete_modes(spec)
    isoviz = render_isoviz(spec)
    topas = render_topas(spec)

    assert "0.35000  0.50000  0.75000" in cif
    assert '"1/4 + Fe1_1_dx"' in cif
    assert "_iso_displacivemodenorm_value" in cif
    assert "_iso_displacivemodematrix_value" in cif
    assert "Undistorted superstructure:" in complete
    assert "Distorted superstructure:" in complete
    assert "  0.250000   0.500000   0.750000" in complete
    assert "  0.350000   0.500000   0.750000" in complete
    assert "0.25000   0.00000   0.00000" in isoviz
    assert "+  0.25000*a1" in topas


def test_exact_fixed_space_four_writers_use_real_irreps_without_fake_primary():
    data = _export_data(
        global_irreps=("GM1+", "GM2+", "GM3+"),
        context_kind="exact_fixed_space",
    )
    spec = _spec(data)
    spec.subgroup.basis_raw = "BOGUS-BASIS"
    spec.subgroup.origin_raw = "BOGUS-ORIGIN"

    cif = render_cif(spec.cif_structure, spec)
    complete = render_complete_modes(spec)
    isoviz = render_isoviz(spec)
    topas = render_topas(spec)

    assert "Exact fixed-space coupled embedding" in cif
    assert "basis=(1,0,0),(0,1,0),(0,0,1), origin=(0,0,0)" in cif
    assert "[0,0,0]" not in cif.split("# Order parameter values:", 1)[1]
    assert "exact fixed-space coupled embedding" in complete
    assert "IR=None" not in complete
    irrep_block = isoviz.split("!irreplist ", 1)[1].split("\n\n", 1)[0]
    assert "IR" not in irrep_block.split()
    assert all(irrep in irrep_block for irrep in ("GM1+", "GM2+", "GM3+"))
    assert "'{{{mode definitions" in topas
    assert all(
        "BOGUS-BASIS" not in payload and "BOGUS-ORIGIN" not in payload
        for payload in (cif, complete, isoviz, topas)
    )

    spec.subgroup.opd_symbol = "BOGUS"
    for writer in (
        lambda: render_cif(spec.cif_structure, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    ):
        with pytest.raises(
            DisplaciveModelError,
            match="incorrectly claims a primary route",
        ):
            writer()

    spec.subgroup.opd_symbol = ""
    spec.subgroup.k_point_label = "X"
    spec.subgroup.k_coordinates = ["1/2", "0", "0"]
    for writer in (
        lambda: render_cif(spec.cif_structure, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    ):
        with pytest.raises(
            DisplaciveModelError,
            match="incorrectly claims a primary route",
        ):
            writer()


def test_exact_fixed_space_export_rejects_mixed_primary_context_provenance():
    source = _export_data(context_kind="exact_fixed_space")
    first = source.model.canonical_modes[0]
    mixed_provenance = replace(
        first.microscopic_provenance,
        source_subgroup_context_kind="single_irrep",
        source_subgroup_irrep_label="GM1+",
        source_subgroup_opd_symbol="P1",
        source_subgroup_primary_direction_selector="VECTOR,A",
    )
    mixed_modes = (
        replace(first, microscopic_provenance=mixed_provenance),
        *source.model.canonical_modes[1:],
    )
    mixed_model = build_displacive_cif_model(
        source.model.frame,
        source.model.atoms,
        tuple(orbit.representative_atom_id for orbit in source.model.orbits),
        mixed_modes,
        amplitudes=source.model.amplitudes,
    )

    with pytest.raises(
        DisplaciveModelError,
        match=r"query context differs|mix or mismatch",
    ):
        DisplaciveExportData(
            model=mixed_model,
            embedding=source.embedding,
            subgroup_identity=source.subgroup_identity,
            reference_structure=source.reference_structure,
            final_structure=source.final_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )


def test_isoviz_assigns_stable_numbers_to_two_explicit_displacive_irreps() -> None:
    spec = _spec(
        _export_data(global_irreps=("GM1+", "GM2+", "GM1+"))
    )
    text = render_isoviz(spec)

    irrep_block = text.split("!irreplist ", 1)[1].split("\n\n", 1)[0]
    irrep_rows = [
        (int(parts[0]), parts[1])
        for line in irrep_block.splitlines()
        if len(parts := line.split()) == 2 and parts[0].isdigit()
    ]
    assert irrep_rows == [(1, "GM1+"), (2, "GM2+")]

    mode_block = text.split("!displacivemodelist ", 1)[1].split("\n\n", 1)[0]
    mode_rows = [
        parts
        for line in mode_block.splitlines()
        if len(parts := line.split()) >= 6
        and parts[0].isdigit()
        and parts[1].isdigit()
    ]
    assert [int(parts[4]) for parts in mode_rows] == [1, 2, 1]


def test_isoviz_routes_same_species_modes_by_parent_orbit_contract() -> None:
    spec = _spec(_two_same_species_orbit_export_data())

    text = render_isoviz(spec)

    atom_types = text.split("!atomtypelist ", 1)[1].split("\n\n", 1)[0]
    assert "1 Fe1 Fe" in atom_types
    assert "2 Fe2 Fe" in atom_types
    mode_block = text.split("!displacivemodelist ", 1)[1].split("\n\n", 1)[0]
    mode_rows = [
        parts
        for line in mode_block.splitlines()
        if len(parts := line.split()) >= 6
        and parts[0].isdigit()
        and parts[1].isdigit()
    ]
    assert [int(parts[0]) for parts in mode_rows] == [1, 1, 1, 2, 2, 2]


def test_export_data_rejects_final_lattice_outside_validated_frame() -> None:
    data = _export_data()
    assert data.final_structure is not None
    wrong_final = Structure(
        Lattice.tetragonal(9.0, 11.0),
        [site.species for site in data.final_structure],
        data.final_structure.frac_coords,
        labels=[site.label for site in data.final_structure],
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=data.model,
            embedding=data.embedding,
            subgroup_identity=data.subgroup_identity,
            reference_structure=data.reference_structure,
            final_structure=wrong_final,
            representative_labels=data.representative_labels,
            parent_orbit_types=data.parent_orbit_types,
        )

    assert caught.value.code == "final_frame_mismatch"


def _pminus_one_export_data() -> DisplaciveExportData:
    identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    inversion = ((-1, 0, 0), (0, -1, 0), (0, 0, -1))
    operations = (
        ExactSeitzOperation.from_values(identity, (0, 0, 0)),
        ExactSeitzOperation.from_values(inversion, (0, 0, 0)),
    )
    frame = ExactChildFrame.from_values(
        "pminus1_exact",
        operations,
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)),
    )
    atoms = (
        ChildAtom.from_values("first", (Fraction(1, 5),) * 3, "Fe"),
        ChildAtom.from_values("representative", (Fraction(4, 5),) * 3, "Fe"),
    )
    modes = tuple(
        _resolved_mode(
            axis=axis,
            frame_id=frame.frame_id,
            atom_ids=tuple(atom.atom_id for atom in atoms),
            displacements=(-np.eye(3)[axis], np.eye(3)[axis]),
            global_irrep="GM1+",
            parent_sg=2,
            wyckoff_letter="i",
            source_child_sg=2,
        )
        for axis in range(3)
    )
    model = build_displacive_cif_model(
        frame,
        atoms,
        ("representative",),
        modes,
        amplitudes=(0.1, 0.0, 0.0),
    )
    reference = Structure(
        Lattice.cubic(4.0),
        ["Fe", "Fe"],
        [[0.2, 0.2, 0.2], [0.8, 0.8, 0.8]],
        labels=["first", "representative"],
    )
    parent = Structure.from_sites(reference)
    embedding = ExactParentChildEmbedding(
        parent,
        2,
        I3,
        ZERO,
        (
            _site_mapping("first", 0, "orbit-i"),
            _site_mapping("representative", 1, "orbit-i"),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(
            embedding, target_space_group_number=2
        ),
        reference_structure=reference,
        representative_labels={"representative": "Fe1_1"},
        parent_orbit_types={
            "orbit-i": ParentOrbitType(
                1,
                "Fe1",
                ("first", "representative"),
            ),
        },
    )


def _pminus_one_nonzero_origin_export_data() -> DisplaciveExportData:
    """Return an exact child already expressed after a nonzero B,q embedding."""

    source = _pminus_one_export_data()
    source_origin = ("1/2", "0", "0")
    source_modes = tuple(
        replace(
            mode,
            microscopic_provenance=replace(
                mode.microscopic_provenance,
                source_subgroup_origin=source_origin,
            ),
        )
        for mode in source.model.canonical_modes
    )
    model = build_displacive_cif_model(
        source.model.frame,
        source.model.atoms,
        ("first",),
        source_modes,
        amplitudes=source.model.amplitudes,
    )
    parent = Structure(
        source.parent_structure.lattice,
        [site.species for site in source.parent_structure],
        [[0.7, 0.2, 0.2], [0.3, 0.8, 0.8]],
        labels=["first", "representative"],
    )
    embedding = ExactParentChildEmbedding(
        parent,
        2,
        I3,
        (Fraction(1, 2), Fraction(0), Fraction(0)),
        (
            _site_mapping("first", 0, "orbit-i"),
            _site_mapping("representative", 1, "orbit-i", (1, 0, 0)),
        ),
    )
    return DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(
            embedding, target_space_group_number=2
        ),
        reference_structure=source.reference_structure,
        representative_labels={"first": "Fe1_1"},
        parent_orbit_types=source.parent_orbit_types,
    )


def test_cif_uses_exact_contract_frame_without_reapplying_subgroup_origin() -> None:
    data = _pminus_one_nonzero_origin_export_data()
    spec = SubgroupExportSpec(
        subgroup=_subgroup_for_data(data),
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    setting = _setting(spec.subgroup, force_p1=False)
    _legacy_structure, legacy_shift = _apply_origin_choice(
        spec.cif_structure,
        setting,
        (-0.5, 0.0, 0.0),
    )

    assert legacy_shift == pytest.approx((0.5, 0.0, 0.0))
    _setting_value, shifted, _sites, writer_shift = _validated_subgroup_sites(
        spec,
        spec.cif_structure,
    )
    assert writer_shift == pytest.approx((0.0, 0.0, 0.0))
    assert shifted.frac_coords == pytest.approx(spec.cif_structure.frac_coords)

    cif = render_cif(spec.cif_structure, spec)
    assert "_iso_parent-to-child.transform_Pp_abc a,b,c;1/2,0,0" in cif


def test_cif_rejects_unrecorded_shift_after_exact_contract_frame() -> None:
    data = _pminus_one_nonzero_origin_export_data()
    spec = SubgroupExportSpec(
        subgroup=_subgroup_for_data(data),
        displacive_data=data,
        require_verified_displacive_data=True,
    )
    setting, _shifted, sites, _writer_shift = _validated_subgroup_sites(
        spec,
        spec.cif_structure,
    )

    with pytest.raises(ValueError, match="unrecorded writer origin shift"):
        _validate_displacive_writer_frame(
            spec,
            setting,
            sites,
            np.asarray((0.5, 0.0, 0.0)),
        )


def test_cif_fails_closed_when_writer_selects_an_equivalent_nonrepresentative_atom() -> None:
    data = _pminus_one_export_data()
    subgroup = SubgroupInfo(
        index=0,
        space_group_number=2,
        space_group_symbol="P-1",
        size=1,
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        origin=[0, 0, 0],
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        parent_sg=2,
    )
    subgroup._displacive_embedding_id = data.embedding.embedding_id
    spec = SubgroupExportSpec(
        subgroup=subgroup,
        displacive_data=data,
        require_verified_displacive_data=True,
    )

    with pytest.raises(ValueError, match="representative differs"):
        render_cif(spec.cif_structure, spec)


def test_cif_binds_contract_labels_when_coordinate_sort_reverses_orbit_order() -> None:
    source = _index_two_export_data()
    atoms = (
        ChildAtom.from_values(
            "a0", (Fraction(5, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
        ChildAtom.from_values(
            "a1", (Fraction(1, 8), Fraction(1, 4), Fraction(1, 4)), "Fe"
        ),
    )
    model = build_displacive_cif_model(
        source.model.frame,
        atoms,
        ("a0", "a1"),
        source.model.canonical_modes,
        amplitudes=source.model.amplitudes,
    )
    reference = Structure(
        source.reference_structure.lattice,
        ["Fe", "Fe"],
        [[0.625, 0.25, 0.25], [0.125, 0.25, 0.25]],
        labels=["Fe1", "Fe1"],
    )
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        ZERO,
        (
            _site_mapping("a0", 0, "orbit-a", (1, 0, 0)),
            _site_mapping("a1", 0, "orbit-a", (0, 0, 0)),
        ),
    )
    data = DisplaciveExportData(
        model=model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels={"a0": "Fe1_1", "a1": "Fe1_2"},
        parent_orbit_types=source.parent_orbit_types,
    )
    spec = SubgroupExportSpec(
        subgroup=_subgroup_for_data(data),
        displacive_data=data,
        require_verified_displacive_data=True,
    )

    cif = render_cif(spec.cif_structure, spec)
    atom_rows = [
        line.split()
        for line in cif.splitlines()
        if line.startswith(("Fe1_1 ", "Fe1_2 ")) and len(line.split()) >= 8
    ]

    assert [(row[0], float(row[4])) for row in atom_rows] == [
        ("Fe1_2", 0.125),
        ("Fe1_1", 0.675),
    ]
    assert 'Fe1_1_x                      "5/8 + Fe1_1_dx"' in cif
    assert 'Fe1_2_x                      "1/8 + Fe1_2_dx"' in cif


def test_all_writers_revalidate_and_reject_changed_subgroup_embedding() -> None:
    data = _export_data()
    spec = _spec(data)
    spec.subgroup.origin = [0.5, 0, 0]

    writers = (
        lambda: render_cif(spec.cif_structure, spec),
        lambda: render_isoviz(spec),
        lambda: render_complete_modes(spec),
        lambda: render_topas(spec),
    )
    for writer in writers:
        with pytest.raises(DisplaciveModelError) as caught:
            writer()
        assert caught.value.code == "subgroup_embedding_mismatch"


def test_embedding_proves_reference_lattice_is_basis_times_parent_lattice() -> None:
    data = _export_data()
    wrong_embedding = ExactParentChildEmbedding(
        data.parent_structure,
        data.parent_space_group_number,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        ZERO,
        data.embedding.atom_mappings,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=data.model,
            embedding=wrong_embedding,
            subgroup_identity=_subgroup_identity(wrong_embedding),
            reference_structure=data.reference_structure,
            representative_labels=data.representative_labels,
            parent_orbit_types=data.parent_orbit_types,
        )
    assert caught.value.code == "parent_child_lattice_mismatch"


def test_embedding_parent_space_group_must_match_canonical_mode_identity() -> None:
    data = _export_data()
    wrong_parent = ExactParentChildEmbedding(
        data.parent_structure,
        2,
        I3,
        ZERO,
        data.embedding.atom_mappings,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=data.model,
            embedding=wrong_parent,
            subgroup_identity=_subgroup_identity(wrong_parent),
            reference_structure=data.reference_structure,
            representative_labels=data.representative_labels,
            parent_orbit_types=data.parent_orbit_types,
        )
    assert caught.value.code == "mode_parent_space_group_mismatch"


def test_contract_snapshots_structures_mapping_and_rebinds_writer_inputs() -> None:
    source = _export_data()
    parent = Structure.from_sites(source.parent_structure)
    reference = Structure.from_sites(source.reference_structure)
    labels = dict(source.representative_labels)
    embedding = ExactParentChildEmbedding(
        parent,
        1,
        I3,
        ZERO,
        source.embedding.atom_mappings,
    )
    data = DisplaciveExportData(
        model=source.model,
        embedding=embedding,
        subgroup_identity=_subgroup_identity(embedding),
        reference_structure=reference,
        representative_labels=labels,
        parent_orbit_types=source.parent_orbit_types,
    )
    parent.translate_sites([0], [0.2, 0, 0], frac_coords=True)
    reference.translate_sites([0], [0.2, 0, 0], frac_coords=True)
    labels["a0"] = "changed"

    assert isinstance(data.parent_structure, IStructure)
    assert isinstance(data.reference_structure, IStructure)
    assert isinstance(data.final_structure, IStructure)
    assert data.parent_structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert data.reference_structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert data.representative_label("a0") == "Fe1_1"
    with pytest.raises(TypeError):
        data.representative_labels["a0"] = "forbidden"  # type: ignore[index]
    with pytest.raises(TypeError):
        data.parent_orbit_types["orbit-a"] = ParentOrbitType(  # type: ignore[index]
            1,
            "changed",
            ("a0",),
        )

    spec = _spec(data)
    assert spec.mode_displacements_sc is not None
    first_key = data.mode_keys[0]
    spec.structure = Structure(Lattice.cubic(9), ["Fe"], [[0, 0, 0]])
    spec.mode_displacements_sc[first_key][0, 0] = 99.0
    spec.mode_labels[first_key] = "changed"
    spec.parent_symbol = "forged"
    render_complete_modes(spec)
    assert isinstance(spec.structure, IStructure)
    assert spec.structure.frac_coords[0, 0] == pytest.approx(0.25)
    assert spec.mode_displacements_sc[first_key][0, 0] == pytest.approx(1.0)
    assert spec.mode_labels[first_key] == "mode-0"
    assert spec.parent_symbol == "P1"


@pytest.mark.parametrize("target", ["parent", "reference", "final"])
def test_contract_rejects_non_single_species_or_partial_occupancy(target: str) -> None:
    source = _export_data()
    disordered = Structure(
        source.reference_structure.lattice,
        [{"Fe": 0.5, "Co": 0.5}],
        source.reference_structure.frac_coords,
        labels=["Fe1"],
    )
    parent = disordered if target == "parent" else source.parent_structure
    if target == "parent":
        with pytest.raises(DisplaciveModelError) as caught:
            ExactParentChildEmbedding(
                parent,
                1,
                I3,
                ZERO,
                source.embedding.atom_mappings,
            )
    else:
        with pytest.raises(DisplaciveModelError) as caught:
            DisplaciveExportData(
                model=source.model,
                embedding=source.embedding,
                subgroup_identity=source.subgroup_identity,
                reference_structure=(
                    disordered if target == "reference" else source.reference_structure
                ),
                final_structure=(disordered if target == "final" else None),
                representative_labels=source.representative_labels,
                parent_orbit_types=source.parent_orbit_types,
            )
    assert caught.value.code == "unsupported_site_occupancy"


def test_spec_rejects_parallel_inputs_when_displacive_contract_is_present() -> None:
    data = _export_data()
    with pytest.raises(ValueError, match="only source for structure"):
        SubgroupExportSpec(
            subgroup=_spec(data).subgroup,
            structure=data.reference_structure,
            displacive_data=data,
        )

    with pytest.raises(ValueError, match="only source for parent_symbol"):
        SubgroupExportSpec(
            subgroup=_spec(data).subgroup,
            parent_symbol="forged",
            displacive_data=data,
        )


def test_direct_isoviz_entry_rejects_removed_contract_even_with_empty_views() -> None:
    spec = _spec(_export_data())
    spec.displacive_data = None
    spec.mode_displacements_sc = {}
    spec.mode_labels = {}
    spec.amplitudes = {}
    spec.require_verified_displacive_data = False

    with pytest.raises(ValueError, match="contract was removed"):
        render_isodistort_isoviz(spec)


def test_coordinate_formula_does_not_reduce_fixed_space_coefficients_modulo_one() -> None:
    identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    sheared_reflection = ((1, 0, 0), (4, -1, 0), (0, 0, 1))
    operations = (
        ExactSeitzOperation.from_values(identity, (0, 0, 0)),
        ExactSeitzOperation.from_values(sheared_reflection, (0, 0, 0)),
    )
    shear = np.array([[1.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    inverse = np.linalg.inv(shear)
    invariant_metric = inverse.T @ inverse
    lattice = np.linalg.cholesky(invariant_metric)
    frame = ExactChildFrame.from_values("sheared_order2", operations, lattice)
    site = site_fixed_coordinate_model("a", (0, 0, 0), frame.operations)

    assert site.symmform == "Dx,2*Dx,Dz"
    parameters = tuple(f"p{index}" for index in range(site.dimension))
    assert _coordinate_formula("0", site.fixed_basis[1], parameters) == "0 + 2*p0"


def _database_operations(hall_number: int) -> tuple[ExactSeitzOperation, ...]:
    database = spglib.get_symmetry_from_database(hall_number)
    assert database is not None
    return tuple(
        ExactSeitzOperation.from_values(rotation, translation)
        for rotation, translation in zip(
            database["rotations"],
            database["translations"],
            strict=True,
        )
    )


def _database_centering_index(operations: tuple[ExactSeitzOperation, ...]) -> int:
    identity = tuple(tuple(Fraction(value) for value in row) for row in I3)
    return sum(operation.rotation == identity for operation in operations)


@pytest.mark.parametrize(
    "hall_number",
    (1, 424, 523, 185, 186, 188, 458, 459),
)
def test_emitted_space_group_accepts_exact_declared_hall_setting(
    hall_number: int,
) -> None:
    operations = _database_operations(hall_number)
    target = spglib.get_spacegroup_type(hall_number)
    assert target is not None

    _validate_emitted_space_group(
        operations,
        model_centering_index=_database_centering_index(operations),
        target_space_group_number=target.number,
        target_hall_number=hall_number,
    )


@pytest.mark.parametrize(
    ("target_hall", "emitted_hall", "error_code"),
    (
        (1, 2, "target_space_group_type_mismatch"),
        (424, 1, "target_centering_mismatch"),
        (523, 1, "target_centering_mismatch"),
        (185, 186, "target_space_group_setting_mismatch"),
        (186, 188, "target_space_group_setting_mismatch"),
        (188, 185, "target_space_group_setting_mismatch"),
        (458, 459, "target_centering_mismatch"),
        (459, 458, "target_centering_mismatch"),
    ),
)
def test_emitted_space_group_rejects_wrong_centering_type_or_setting(
    target_hall: int,
    emitted_hall: int,
    error_code: str,
) -> None:
    operations = _database_operations(emitted_hall)
    target = spglib.get_spacegroup_type(target_hall)
    assert target is not None

    with pytest.raises(DisplaciveModelError) as caught:
        _validate_emitted_space_group(
            operations,
            model_centering_index=_database_centering_index(operations),
            target_space_group_number=target.number,
            target_hall_number=target_hall,
        )

    assert caught.value.code == error_code


def test_export_contract_rejects_model_with_wrong_target_space_group_type() -> None:
    source = _export_data()
    identity = _subgroup_identity(
        source.embedding,
        target_space_group_number=2,
        target_hall_number=2,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=source.embedding,
            subgroup_identity=identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "target_space_group_type_mismatch"


def test_subgroup_identity_rejects_hall_setting_from_another_space_group() -> None:
    source = _export_data()

    with pytest.raises(DisplaciveModelError) as caught:
        _subgroup_identity(
            source.embedding,
            target_space_group_number=1,
            target_hall_number=2,
        )

    assert caught.value.code == "invalid_target_hall_setting"


def test_export_contract_rejects_model_with_wrong_target_centering() -> None:
    source = _export_data()
    embedding = ExactParentChildEmbedding(
        source.parent_structure,
        1,
        ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
        ZERO,
        source.embedding.atom_mappings,
    )
    identity = _subgroup_identity(
        embedding,
        target_space_group_number=5,
        target_hall_number=9,
    )

    with pytest.raises(DisplaciveModelError) as caught:
        DisplaciveExportData(
            model=source.model,
            embedding=embedding,
            subgroup_identity=identity,
            reference_structure=source.reference_structure,
            representative_labels=source.representative_labels,
            parent_orbit_types=source.parent_orbit_types,
        )

    assert caught.value.code == "target_centering_mismatch"


def test_periodic_cartesian_residual_uses_skew_lattice_nearest_image() -> None:
    lattice = np.asarray(
        ((1.0, 0.0, 0.0), (0.9, 0.1, 0.0), (0.0, 0.0, 1.0)),
        dtype=float,
    )
    first = np.asarray(((0.49, 0.49, 0.0),), dtype=float)
    second = np.zeros((1, 3), dtype=float)

    residual = _periodic_cartesian_residual(first, second, lattice)

    assert residual == pytest.approx(np.hypot(0.031, 0.051))
    componentwise_residual = np.linalg.norm(
        (first[0] - np.round(first[0])) @ lattice
    )
    assert componentwise_residual > 15 * residual
