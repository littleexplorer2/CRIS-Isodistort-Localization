"""Tests for test_verify_parametric_k_grid; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.utils.lattice import centering_primitive_matrix
from backend.utils.lattice import determinant
from backend.utils.lattice import identity_matrix
from backend.utils.lattice import inverse
from backend.utils.lattice import multiply
from backend.utils.lattice import rational_matrix
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
from tests import measure_method3_closure as retired
from tests import verify_parametric_k_grid as audit
import json
import pytest


# Original coverage: test_parametric_k_grid.py.

@pytest.mark.parametrize('centering', list('PA BCIFR'.replace(' ', '')))
def test_reciprocal_selection_rules_against_conventional_translations(centering):
    # Independent conventional centering conditions (R: obverse hexagonal axes).
    predicates = {
        'P': lambda h, k, l: True, 'A': lambda h, k, l: (k+l) % 2 == 0,
        'B': lambda h, k, l: (h+l) % 2 == 0, 'C': lambda h, k, l: (h+k) % 2 == 0,
        'I': lambda h, k, l: (h+k+l) % 2 == 0,
        'F': lambda h, k, l: h % 2 == k % 2 == l % 2,
        'R': lambda h, k, l: (-h+k+l) % 3 == 0,
    }
    primitive = centering_primitive_matrix(centering)
    for delta in product(range(-2, 3), repeat=3):
        assert audit.same_point_mod_lattice(delta, (0, 0, 0), primitive) == predicates[centering](*delta)
    assert not audit.same_point_mod_lattice(('1/2', 0, 0), (0, 0, 0), primitive)

def test_body_centered_gamma_and_111_are_distinct():
    assert not audit.same_point_mod_lattice((1, 1, 1), (0, 0, 0), centering_primitive_matrix('I'))
    assert not audit.same_star((1, 1, 1), (0, 0, 0), centering_primitive_matrix('I'))
    assert audit.same_point_mod_lattice((1, 1, 0), (0, 0, 0), centering_primitive_matrix('I'))
    assert audit.same_point_mod_lattice((1, 1, 1), (0, 0, 0), centering_primitive_matrix('P'))

def test_primitive_generators_and_coordinate_frame_do_not_change_equivalence():
    primitive = centering_primitive_matrix('I')
    unimodular = rational_matrix(((1, 1, 0), (0, 1, 0), (0, 0, -1)))
    frame = rational_matrix(((2, 1, 0), (0, 1, 1), (0, 0, 1)))
    def mv(matrix, vector):
        return tuple(sum(a*b for a, b in zip(row, vector)) for row in matrix)
    for vector in product(range(-1, 2), repeat=3):
        expected = audit.same_point_mod_lattice(vector, (0, 0, 0), primitive)
        assert audit.same_point_mod_lattice(vector, (0, 0, 0), multiply(unimodular, primitive)) == expected
        assert audit.same_point_mod_lattice(mv(inverse(frame), vector), (0, 0, 0),
                                            multiply(primitive, frame)) == expected

def test_nonorthogonal_hexagonal_rotation_uses_inverse_transpose():
    rotation = rational_matrix(((0, -1, 0), (1, -1, 0), (0, 0, 1)))
    k = (F(1, 3), F(1, 5), F(1, 7))
    dual_k = audit.reciprocal_action(rotation, k)
    assert dual_k == (-k[0]-k[1], k[0], k[2])
    t = (2, -3, 1)
    rotated_t = tuple(sum(a*b for a, b in zip(row, t)) for row in rotation)
    assert sum(a*b for a, b in zip(dual_k, rotated_t)) == sum(a*b for a, b in zip(k, t))

def test_symbolic_generic_stabilizer_detects_distinct_parameter_exception():
    primitive = identity_matrix()
    template = ['a', 'b', '0']
    generic = audit.generic_little_group(template, ('a', 'b'), primitive)
    # a != b does not exclude a+b=1, where a swapped/reflected arm coincides.
    actual = audit.little_group(audit.coordinates(template, {'a': '1/3', 'b': '2/3'}), primitive)
    assert len(generic) == 2
    assert generic < actual

@pytest.mark.parametrize('template', [['a*b', '0', '0'], ['a*a', '0', '0'], ['unknown', '0', '0']])
def test_nonaffine_or_unknown_parameter_is_rejected(template):
    with pytest.raises(ValueError, match='affine'):
        audit.affine_template(template, ('a', 'b'))

def test_centering_changes_boundary_little_groups():
    for k, expected in [(('1/2', '1/4', '1/2'), 2), (('1/2', '1/2', '1/3'), 4)]:
        assert len(audit.little_group(k, centering_primitive_matrix('I'))) == expected
        assert len(audit.little_group(k, centering_primitive_matrix('P'))) == 2 * expected

def test_star_equivalence_includes_rotated_nonidentical_points():
    primitive = centering_primitive_matrix('P')
    left, right = ('1/5', 0, 0), (0, '1/5', 0)
    assert not audit.same_point_mod_lattice(left, right, primitive)
    assert audit.same_star(left, right, primitive)

@pytest.fixture(scope='module')
def manifest():
    return json.loads(audit.DEFAULT_MANIFEST.read_text(encoding='utf-8'))

def test_frozen_input_scope_and_derived_evidence(manifest):
    report = audit.build_report(manifest)
    assert report['verdict'] == 'pass', report['issues']
    assert report['query_count']['download_queries'] == 5
    assert report['query_count']['registered_queries'] == 3
    assert report['query_count']['blocked_queries'] == 1
    assert report['query_count']['parametric_queries_download'] == 5
    assert report['query_count']['parametric_queries_total'] == 8
    assert report['query_count']['scope_queries_total'] == 9
    assert report['query_count']['parametric_labels_total'] == 4
    assert report['label_aliases'] == report['same_star_pairs'] == []
    assert len(audit.POINT_GROUP) == 16
    for left, right in product(audit.POINT_GROUP, repeat=2):
        assert multiply(left, right) in audit.POINT_GROUP
    for query in manifest['method2_parametric']['queries']:
        entry = report['parents'][query['parent']]['parametric_labels'][query['k_label']]
        for key in ('little_group_order', 'generic_little_group_order', 'general_position'):
            assert query[key] == entry[key]
        download = query.get('download_nmod_values', query['nmod_values'])
        registered = query.get('registered_nmod_values', [])
        blocked = query.get('blocked_nmod_values', [])
        assert sorted(download + registered + blocked) == query['nmod_values']
        assert not set(download) & set(registered)
        assert entry['frozen_archived_value'] == bool(registered)

def test_registered_nmod0_queries_are_the_only_archive_reuse_queries(manifest):
    registered = {}
    for query in manifest['method2_fixed_k_queries']['queries'] + manifest['method2_parametric']['queries']:
        if query.get('registered_nmod_values'):
            registered[query['id']] = query['registered_nmod_values']
    assert registered == {'M2P-EU-LD': [0], 'M2P-ND-Y': [0], 'M2F-C4310-X': [0]}
    for query in manifest['method2_fixed_k_queries']['queries'] + manifest['method2_parametric']['queries']:
        if query['id'] == 'M2F-C4310-X':
            assert query['nmod_values'] == [0]
            assert query['download_nmod_values'] == []
        elif query['id'] in registered:
            assert query['download_nmod_values'] == [1]
        elif query['id'] == 'M2P-ND-GP':
            assert query['download_nmod_values'] == [0]
            assert query['blocked_nmod_values'] == [1]
        else:
            assert query.get('download_nmod_values', query['nmod_values']) == [0, 1]

def test_manifest_claims_and_expected_answers_cannot_affect_arithmetic(manifest):
    changed = deepcopy(manifest)
    for query in changed['method2_parametric']['queries']:
        query.update(little_group_order=999, general_position=False, expected_candidates_nmod0=123456)
    assert audit.input_contract(changed) == audit.input_contract(manifest)

@pytest.mark.parametrize('mutation', ['coordinate', 'missing', 'duplicate', 'hash', 'parameter', 'partition'])
def test_input_drift_fails_closed(manifest, mutation):
    changed = deepcopy(manifest)
    rows = changed['method2_parametric']['queries']
    if mutation == 'coordinate':
        rows[0]['coordinates'][0] = '7/13'
    elif mutation == 'missing':
        rows.pop()
    elif mutation == 'duplicate':
        rows.append(deepcopy(rows[0]))
    elif mutation == 'hash':
        next(iter(changed['parents'].values()))['sha256'] = '0' * 64
    elif mutation == 'partition':
        rows[0]['download_nmod_values'] = [0]
        rows[0]['registered_nmod_values'] = []
    else:
        rows[0]['values']['unknown'] = '1/5'
    assert audit.build_report(changed)['verdict'] == 'fail'

def test_unsupported_setting_and_float_parameters_rejected(manifest):
    changed = deepcopy(manifest)
    next(iter(changed['parents'].values()))['space_group_symbol'] = 'P4/mmm'
    with pytest.raises(ValueError, match='setting'):
        audit.build_report(changed)
    changed = deepcopy(manifest)
    values = changed['method2_parametric']['queries'][0]['values']
    values[next(iter(values))] = 0.1666666666666
    with pytest.raises(ValueError, match='exact rational'):
        audit.build_report(changed)


def test_subset_requires_explicit_scope_and_cannot_silently_drop_queries(manifest):
    changed = deepcopy(manifest)
    changed.pop('method2_query_scope')
    assert audit.build_report(changed)['verdict'] == 'fail'
    changed['method2_query_scope'] = {'mode': 'unknown'}
    assert audit.build_report(changed)['verdict'] == 'fail'
    changed = deepcopy(manifest)
    changed['method2_query_scope']['query_ids'].pop()
    assert audit.build_report(changed)['verdict'] == 'fail'


@pytest.mark.parametrize('nmods', [[], [0, 0], [False], [2]])
def test_invalid_selected_nmod_scope_is_rejected(manifest, nmods):
    changed = deepcopy(manifest)
    changed['method2_fixed_k_queries']['queries'][0]['nmod_values'] = nmods
    assert audit.build_report(changed)['verdict'] == 'fail'


@pytest.mark.parametrize('mutation', ['missing', 'overlap', 'duplicate', 'boolean'])
def test_blocked_queries_must_be_a_disjoint_explicit_partition(manifest, mutation):
    changed = deepcopy(manifest)
    query = next(q for q in changed['method2_parametric']['queries']
                 if q.get('blocked_nmod_values'))
    if mutation == 'missing':
        query.pop('blocked_nmod_values')
    elif mutation == 'overlap':
        query['download_nmod_values'].append(1)
    elif mutation == 'duplicate':
        query['blocked_nmod_values'] = [1, 1]
    else:
        query['blocked_nmod_values'] = [True]
    assert audit.build_report(changed)['verdict'] == 'fail'


def test_old_method3_command_cannot_resume_contaminated_checkpoints(capsys):
    assert retired.main(['--resume', '--dry-run']) == 2
    assert 'INVALIDATED' in capsys.readouterr().err

def test_conventional_determinant_is_not_translation_index():
    # Same conventional cell, I parent -> P child: index TWO, not det(B)=1.
    parent, child = centering_primitive_matrix('I'), identity_matrix()
    basis = identity_matrix()
    child_in_parent = multiply(multiply(child, basis), inverse(parent))
    assert all(value.denominator == 1 for row in child_in_parent for value in row)
    assert abs(determinant(child_in_parent)) == 2
    # Special k of order three exists on a hexagonal reciprocal lattice.
    rotation = rational_matrix(((0, -1, 0), (1, -1, 0), (0, 0, 1)))
    k = (F(1, 3), F(1, 3), F(0))
    assert audit.same_point_mod_lattice(audit.reciprocal_action(rotation, k), k, identity_matrix())
    assert any((2 * component).denominator != 1 for component in k)
