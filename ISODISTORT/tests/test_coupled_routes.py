from fractions import Fraction

import pytest

from backend.utils.lattice import identity_matrix, rational_matrix
from backend.wrappers import SubgroupInfo
from features.method1.affine_embeddings import (
    AffineEmbedding,
    AffineOperation,
    ParentAffineGroup,
    build_affine_quotient,
)
from features.method3.coupled_routes import (
    CoupledRouteCandidate,
    lift_embedding_subgroup,
    resolve_coupled_route_witness,
)

ZERO = (Fraction(0), Fraction(0), Fraction(0))


def _operation(rotation):
    return AffineOperation(rational_matrix(rotation), ZERO)


def _route(label: str, space_group_number: int) -> SubgroupInfo:
    return SubgroupInfo(
        index=0,
        space_group_number=space_group_number,
        irrep_label=label,
        opd_symbol="P1",
        basis_vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )


def _klein_parent():
    identity = _operation([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    first = _operation([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
    second = _operation([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    product = _operation([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    parent = ParentAffineGroup(
        lattice=identity_matrix(),
        operations=(identity, first, second, product),
        space_group_number=47,
        hall_number=227,
        symbol="Pmmm",
    )
    target = AffineEmbedding(identity_matrix(), (identity,))
    candidates = [
        CoupledRouteCandidate(
            _route("GM2+", 10),
            AffineEmbedding(identity_matrix(), (identity, first)),
        ),
        CoupledRouteCandidate(
            _route("GM4+", 10),
            AffineEmbedding(identity_matrix(), (identity, second)),
        ),
    ]
    return parent, target, candidates


def test_resolve_coupled_route_witness_uses_exact_stabilizer_intersection():
    parent, target, candidates = _klein_parent()

    witness = resolve_coupled_route_witness(parent, target, candidates)

    assert witness is not None
    assert [route.irrep_label for route in witness.components] == ["GM2+", "GM4+"]
    assert witness.target_subgroup_order == 1
    assert witness.quotient_order == 4
    assert witness.verification_method == "exact_affine_stabilizer_intersection"


def test_lift_embedding_subgroup_adds_translation_cosets():
    parent, _, _ = _klein_parent()
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
    quotient = build_affine_quotient(parent, child_lattice)
    identity_route = AffineEmbedding(identity_matrix(), (parent.operations[0],))

    lifted = lift_embedding_subgroup(identity_route, quotient)

    assert len(lifted) == 2
    assert quotient.identity_index in lifted


def test_coupled_route_state_budget_fails_without_truncation():
    parent, target, candidates = _klein_parent()

    with pytest.raises(RuntimeError, match="No candidate result was truncated"):
        resolve_coupled_route_witness(
            parent,
            target,
            candidates,
            max_intersection_states=1,
        )
