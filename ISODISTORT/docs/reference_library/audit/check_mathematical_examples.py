"""Independent examples for the written derivations, not CRIS engine tests."""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def exact_stabilizers():
    identity = sp.eye(2)
    generators = [sp.Matrix([[0, -1], [1, 0]]), sp.diag(1, -1)]
    group = {tuple(identity): identity}
    queue = [identity]
    while queue:
        item = queue.pop()
        for gen in generators:
            product = gen * item
            if tuple(product) not in group:
                group[tuple(product)] = product
                queue.append(product)
    keys, matrices = list(group), list(group.values())
    unit = keys.index(tuple(identity))
    multiply = [[keys.index(tuple(a * b)) for b in matrices] for a in matrices]
    subgroups = []
    for bits in range(1 << len(keys)):
        subgroup = {i for i in range(len(keys)) if bits & (1 << i)}
        if unit in subgroup and all(multiply[i][j] in subgroup for i in subgroup for j in subgroup):
            subgroups.append(subgroup)
    witnesses = 0
    for representation in (matrices, [sp.diag(1, m) for m in matrices]):
        dim = representation[0].rows
        for subgroup in subgroups:
            fixed = sp.Matrix.vstack(*(representation[i] - sp.eye(dim) for i in sorted(subgroup))).nullspace()
            if not fixed:
                continue
            basis = sp.Matrix.hstack(*fixed)
            pointwise = {i for i, matrix in enumerate(representation) if (matrix - sp.eye(dim)) * basis == sp.zeros(dim, basis.cols)}
            if pointwise != subgroup:
                continue
            bound = (len(keys) - len(subgroup)) * (len(fixed) - 1)
            found = False
            for parameter in range(bound + 1):
                vector = sum((parameter ** j * column for j, column in enumerate(fixed)), sp.zeros(dim, 1))
                actual = {i for i, matrix in enumerate(representation) if matrix * vector == vector}
                if actual == subgroup and vector != sp.zeros(dim, 1):
                    found = True
                    witnesses += 1
                    break
            assert found, 'The proved finite witness bound failed on an exact example.'
    return {'group_order': len(keys), 'subgroups_enumerated': len(subgroups), 'representations': 2, 'witnesses_checked': witnesses}


def positive_sqrt(matrix):
    values, vectors = np.linalg.eigh(matrix)
    assert np.min(values) > 0
    return (vectors * np.sqrt(values)) @ vectors.T


def main():
    results = {'scope': 'Illustrations of manuscript contracts only; not a general proof or production-code validation.'}
    results['exact_stabilizers'] = exact_stabilizers()
    lattice = np.array([[1., .9], [0., .1]])
    displacement = np.array([.49, .49])
    rounded = float(np.linalg.norm(lattice @ displacement))
    improved = float(np.linalg.norm(lattice @ (displacement - [1, 0])))
    radius = rounded / np.linalg.svd(lattice, compute_uv=False)[-1]
    ranges = [range(int(np.ceil(v-radius)), int(np.floor(v+radius))+1) for v in displacement]
    minimum = min(float(np.linalg.norm(lattice @ (displacement - np.array(n)))) for n in itertools.product(*ranges))
    assert improved < rounded and minimum <= improved + 1e-14
    results['oblique_nearest_image'] = {'rounding_distance': rounded, 'better_candidate_distance': improved, 'global_distance_with_proved_bound': minimum}
    metric = np.array([[3., .4, .2], [.4, 2., .3], [.2, .3, 1.]])
    expected = np.array([[1.2, .1, 0.], [.1, .8, .08], [0., .08, 1.1]])
    target = expected @ metric @ expected
    root = positive_sqrt(metric)
    invroot = np.linalg.inv(root)
    recovered = invroot @ positive_sqrt(root @ target @ root) @ invroot
    error = float(np.linalg.norm(recovered - expected))
    assert error < 1e-12
    results['spd_noncommuting_metric'] = {'recovered_multiplier_error': error, 'metric_residual': float(np.linalg.norm(recovered @ metric @ recovered - target))}
    r = sp.Matrix([[1, 1], [0, 1]])
    s = sp.Matrix([[1, 0], [1, 1]])
    tensor = sp.diag(2, 3)
    pullback = lambda m, x: m.T * x * m
    assert pullback(r, pullback(s, tensor)) == pullback(s*r, tensor)
    assert pullback(r, pullback(s, tensor)) != pullback(r*s, tensor)
    results['right_action_order'] = 'Exact noncommuting counterexample checked with rational matrices.'
    b = np.eye(2)
    c = np.array([[1., -1.], [1., 1.]]) / np.sqrt(2)
    assert np.allclose(b @ b.T, c @ c.T)
    assert not np.allclose(b.T @ [1, 0], c.T @ [1, 0])
    results['subspace_not_column_identity'] = 'Equal projectors and different coordinates verified.'
    results['status'] = 'passed_for_listed_examples'
    (ROOT / 'audit/mathematical_examples.json').write_text(json.dumps(results, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(results, ensure_ascii=True))


if __name__ == '__main__':
    main()
