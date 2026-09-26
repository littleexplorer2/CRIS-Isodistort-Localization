from __future__ import annotations

from fractions import Fraction

import pytest

from isocore.distortion.inverse_landau import (
    FiniteGroup,
    FixedSpaceFeasibilityAnalyzer,
    RationalRepresentation,
    analyze_fixed_space,
    direct_sum_representations,
    homogeneous_strain_representation,
    polar_vector_representation,
    remove_uniform_site_vectors,
    site_displacement_representation,
    site_vector_character_representation,
)


def _cyclic_two() -> FiniteGroup:
    return FiniteGroup.from_operation(
        ("e", "a"),
        "e",
        lambda left, right: "e" if left == right else "a",
        name="C2",
    )


def _klein_four() -> FiniteGroup:
    elements = ((0, 0), (1, 0), (0, 1), (1, 1))
    return FiniteGroup.from_operation(
        elements,
        (0, 0),
        lambda left, right: (left[0] ^ right[0], left[1] ^ right[1]),
        name="C2 x C2",
    )


def test_trivial_representation_rejects_a_spurious_broken_subgroup() -> None:
    group = _cyclic_two()
    trivial = RationalRepresentation.from_matrices(
        group,
        {"e": ((1,),), "a": ((1,),)},
        name="trivial",
    )

    broken = analyze_fixed_space(trivial, ("e",))
    unbroken = analyze_fixed_space(trivial, ("e", "a"))

    assert broken.fixed_dimension == 1
    assert broken.pointwise_stabilizer == ("e", "a")
    assert not broken.is_reachable
    assert "outside H" in (broken.obstruction or "")
    assert unbroken.is_reachable


def test_arbitrary_matrices_still_require_full_exact_homomorphism() -> None:
    group = _cyclic_two()

    with pytest.raises(ValueError, match="do not form a representation"):
        RationalRepresentation.from_matrices(
            group,
            {"e": ((1,),), "a": ((2,),)},
            name="invalid arbitrary matrices",
        )


def test_homogeneous_strain_cannot_break_a_pure_translation() -> None:
    group = FiniteGroup.from_operation(
        ("identity", "translation"),
        "identity",
        lambda left, right: "identity" if left == right else "translation",
    )
    identity_rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    strain = homogeneous_strain_representation(
        group,
        {
            "identity": identity_rotation,
            "translation": identity_rotation,
        },
    )

    analysis = analyze_fixed_space(strain, ("identity",))

    assert strain.dimension == 6
    assert analysis.fixed_dimension == 6
    assert analysis.pointwise_stabilizer == ("identity", "translation")
    assert not analysis.is_reachable


def test_polar_vector_has_no_inversion_even_fixed_space() -> None:
    group = FiniteGroup.from_operation(
        ("e", "inversion"),
        "e",
        lambda left, right: "e" if left == right else "inversion",
    )
    polar = polar_vector_representation(
        group,
        {
            "e": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            "inversion": ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
        },
    )

    polar_subgroup = analyze_fixed_space(polar, ("e",))
    centrosymmetric = analyze_fixed_space(polar, ("e", "inversion"))

    assert polar_subgroup.is_reachable
    assert polar_subgroup.pointwise_stabilizer == ("e",)
    assert centrosymmetric.projector == (
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert centrosymmetric.fixed_dimension == 0
    assert not centrosymmetric.is_reachable


def test_coupled_direct_sum_intersects_component_stabilizers() -> None:
    group = _klein_four()
    first = RationalRepresentation.from_matrices(
        group,
        {
            element: ((-1 if element[1] else 1,),)
            for element in group.elements
        },
        name="first character",
    )
    second = RationalRepresentation.from_matrices(
        group,
        {
            element: ((-1 if element[0] else 1,),)
            for element in group.elements
        },
        name="second character",
    )

    first_only = analyze_fixed_space(first, ((0, 0),))
    second_only = analyze_fixed_space(second, ((0, 0),))
    coupled = analyze_fixed_space(
        direct_sum_representations((first, second)),
        ((0, 0),),
    )

    assert set(first_only.pointwise_stabilizer) == {(0, 0), (1, 0)}
    assert set(second_only.pointwise_stabilizer) == {(0, 0), (0, 1)}
    assert not first_only.is_reachable
    assert not second_only.is_reachable
    assert coupled.fixed_dimension == 2
    assert coupled.pointwise_stabilizer == ((0, 0),)
    assert coupled.is_reachable


def test_character_analyzer_is_exactly_equivalent_to_dense_projectors() -> None:
    group = _klein_four()
    first = RationalRepresentation.from_matrices(
        group,
        {element: ((-1 if element[1] else 1,),) for element in group.elements},
    )
    second = RationalRepresentation.from_matrices(
        group,
        {element: ((-1 if element[0] else 1,),) for element in group.elements},
    )
    representation = direct_sum_representations((first, second))
    analyzer = FixedSpaceFeasibilityAnalyzer(representation)

    for subgroup in (
        ((0, 0),),
        ((0, 0), (1, 0)),
        group.elements,
    ):
        dense = analyze_fixed_space(representation, subgroup)
        summary = analyzer.analyze(subgroup)
        assert summary.fixed_dimension == dense.fixed_dimension
        assert summary.pointwise_stabilizer == dense.pointwise_stabilizer
        assert summary.is_reachable == dense.is_reachable
        assert summary.obstruction == dense.obstruction


def test_commensurate_site_displacements_can_break_a_translation() -> None:
    group = FiniteGroup.from_operation(
        ("e", "translation"),
        "e",
        lambda left, right: "e" if left == right else "translation",
    )
    identity_rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    displacements = site_displacement_representation(
        group,
        {"e": (0, 1), "translation": (1, 0)},
        {"e": identity_rotation, "translation": identity_rotation},
    )

    broken = analyze_fixed_space(displacements, ("e",))
    translation_even = analyze_fixed_space(displacements, ("e", "translation"))

    assert broken.fixed_dimension == 6
    assert broken.pointwise_stabilizer == ("e",)
    assert broken.is_reachable
    assert translation_even.fixed_dimension == 3
    assert translation_even.projector[0][0] == Fraction(1, 2)
    assert translation_even.projector[0][3] == Fraction(1, 2)
    assert translation_even.is_reachable


def test_internal_site_character_matches_dense_quotient_exactly() -> None:
    group = FiniteGroup.from_operation(
        ("e", "translation"),
        "e",
        lambda left, right: "e" if left == right else "translation",
    )
    identity_rotation = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    permutations = {"e": (0, 1), "translation": (1, 0)}
    vector = polar_vector_representation(
        group,
        {"e": identity_rotation, "translation": identity_rotation},
    )
    full = site_displacement_representation(
        group,
        permutations,
        {"e": identity_rotation, "translation": identity_rotation},
    )
    dense_internal = remove_uniform_site_vectors(full, site_count=2)
    character_internal = site_vector_character_representation(
        group,
        permutations,
        vector,
        remove_uniform_vectors=True,
    )

    dense_analyzer = FixedSpaceFeasibilityAnalyzer(dense_internal)
    character_analyzer = FixedSpaceFeasibilityAnalyzer(character_internal)
    for subgroup in (("e",), group.elements):
        assert character_analyzer.analyze(subgroup) == dense_analyzer.analyze(subgroup)
