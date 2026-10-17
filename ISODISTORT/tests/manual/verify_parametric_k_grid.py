#!/usr/bin/env python3
"""Exact audit of the declared conventional-tetragonal Method 2 input grid.

Reads only the input manifest, original parent hashes and general k tables.
No official result, candidate list, local search or sampled generic order is
used. See reference_library/latex/unsupported_algorithm_proofs.tex, grid-duality.
This checks generic little co-groups on affine templates, not physical modes,
incommensurability, Brillouin-zone domains or completeness of an infinite scope.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = PACKAGE_ROOT.parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from backend.tables.kpoints_official import KPOINT_OFFICIAL
from backend.utils.lattice import (
    centering_primitive_matrix, identity_matrix, inverse, multiply,
    rational_matrix, rational_vector, transpose, wavevector_translation_phases,
)
from backend.utils.text_parser import parse_affine_coordinate

DEFAULT_MANIFEST = PACKAGE_ROOT / 'docs/manifests/official_output_download_manifest.json'
DEFAULT_REPORT = PACKAGE_ROOT / 'output/validation/reciprocal_scope_20261017/grid-report.json'
SETTINGS = {123: 'P4/mmm', 139: 'I4/mmm'}


def point_group():
    """4/mmm in its conventional tetragonal fractional direct frame."""
    generators = [rational_matrix(r) for r in (
        ((0, -1, 0), (1, 0, 0), (0, 0, 1)),
        ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
        ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
    )]
    group, pending = {identity_matrix()}, [identity_matrix()]
    while pending:
        current = pending.pop()
        for generator in generators:
            candidate = multiply(generator, current)
            if candidate not in group:
                group.add(candidate)
                pending.append(candidate)
    return tuple(sorted(group))


POINT_GROUP = point_group()


def reciprocal_action(rotation, vector):
    dual = transpose(inverse(rational_matrix(rotation)))
    k = rational_vector(vector)
    return tuple(sum(a * b for a, b in zip(row, k, strict=True)) for row in dual)


def same_point_mod_lattice(left, right, primitive):
    """P delta integral, with primitive direct translations as rows of P."""
    delta = tuple(a - b for a, b in zip(
        rational_vector(left), rational_vector(right), strict=True))
    return all(phase.denominator == 1 for phase in wavevector_translation_phases(
        delta, rational_matrix(primitive)))


def little_group(k, primitive, rotations=POINT_GROUP):
    return frozenset(r for r in rotations if same_point_mod_lattice(
        reciprocal_action(r, k), k, primitive))


def affine_template(template, names):
    if len(template) != 3 or len(set(names)) != len(names):
        raise ValueError('Expected three affine coordinates and distinct parameter names')
    rows = [parse_affine_coordinate(component, names) for component in template]
    if any(row is None for row in rows):
        raise ValueError(f'Invalid affine template: {template}')
    return (tuple(row[1] for row in rows),
            tuple(tuple(row[0][i] for row in rows) for i in range(len(names))))


def coordinates(template, values):
    names = tuple(values)
    offset, columns = affine_template(template, names)
    return tuple(offset[i] + sum(Fraction(values[name]) * column[i]
                 for name, column in zip(names, columns, strict=True)) for i in range(3))


def generic_little_group(template, names, primitive, rotations=POINT_GROUP):
    """Affine map into a discrete lattice is constant on an open domain."""
    offset, columns = affine_template(template, names)
    return frozenset(r for r in rotations
                     if all(reciprocal_action(r, col) == col for col in columns)
                     and same_point_mod_lattice(reciprocal_action(r, offset), offset, primitive))


def same_star(left, right, primitive, rotations=POINT_GROUP):
    return any(same_point_mod_lattice(reciprocal_action(r, left), right, primitive)
               for r in rotations)


def input_contract(manifest):
    """Allowlist inputs; ignore all claimed orders, counts and pass statuses."""
    parents = {name: {key: row[key] for key in
               ('path', 'sha256', 'space_group_type', 'space_group_symbol')}
               for name, row in manifest['parents'].items()}
    rows = []
    for role, queries in (
        ('fixed', manifest['method2_fixed_k_queries']['queries']),
        ('parametric', manifest['method2_parametric']['queries']),
        ('registered', manifest['method2_parametric']['frozen_archived_values']),
    ):
        for query in queries:
            rows.append(dict(role=role, **{key: query[key] for key in
                ('parent', 'k_label', 'kovalev', 'coordinates')},
                template=query.get('template'), values=query.get('values', {}),
                nmod_values=query.get('nmod_values', [0, 1])))
    return {'parents': parents, 'queries': rows}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    separators=(',', ':')).encode()).hexdigest()


def source_signatures():
    paths = [Path(__file__), PACKAGE_ROOT / 'backend/tables/kpoints_official.py',
             PACKAGE_ROOT / 'backend/utils/lattice.py', PACKAGE_ROOT / 'backend/utils/text_parser.py']
    return {path.relative_to(REPO_ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def build_report(manifest, *, verify_parent_hashes=True):
    contract = input_contract(manifest)
    report = dict(schema='reciprocal-grid-audit/v2', created_at=datetime.now(timezone.utc).isoformat(),
                  input_contract_sha256=digest(contract), source_signatures=source_signatures(),
                  arithmetic='exact rational; symbolic generic stabilizer', parents={},
                  issues=[], label_aliases=[], same_star_pairs=[])
    issues = report['issues']
    if not contract['parents']:
        issues.append('Empty parent scope')
    for query in contract['queries']:
        if query['parent'] not in contract['parents']:
            issues.append(f"Unknown parent: {query['parent']}")
    for parent, metadata in contract['parents'].items():
        sg = metadata['space_group_type']
        if SETTINGS.get(sg) != metadata['space_group_symbol']:
            raise ValueError(f'Unsupported conventional setting: {parent}: {metadata}')
        if verify_parent_hashes:
            path = (REPO_ROOT / metadata['path']).resolve()
            if not path.is_relative_to((REPO_ROOT / 'experiment_data').resolve()):
                raise ValueError('Parent must be an independent experiment_data input')
            if hashlib.sha256(path.read_bytes()).hexdigest() != metadata['sha256']:
                issues.append(f'{parent}: parent SHA-256 mismatch')
        primitive = centering_primitive_matrix(SETTINGS[sg][0])
        table = KPOINT_OFFICIAL[sg]
        queries = [q for q in contract['queries'] if q['parent'] == parent]
        labels = Counter(q['k_label'] for q in queries)
        if set(labels) != set(table) or any(count != 1 for count in labels.values()):
            issues.append(f'{parent}: labels missing, unknown or duplicated')
        result = dict(centering=SETTINGS[sg][0], fixed_labels={}, parametric_labels={})
        points = {}
        for query in queries:
            label = query['k_label']
            if label not in table:
                continue
            kovalev, template, names = table[label]
            if (query['kovalev'] != kovalev or set(query['values']) != set(names)
                    or (query['template'] is not None and query['template'] != template)
                    or (query['role'] == 'fixed') != (not names)
                    or query['nmod_values'] != [0, 1]):
                issues.append(f'{parent}/{label}: input/table contract mismatch')
                continue
            # Only exact text/integer parameters are accepted, never binary floats.
            if any(type(v) not in (str, int) for v in query['values'].values()):
                raise ValueError('Parameters must be exact rational text or integers')
            k = coordinates(template, query['values'])
            if list(map(str, k)) != list(map(str, rational_vector(query['coordinates']))):
                issues.append(f'{parent}/{label}: stored coordinates disagree with parameters')
            actual = little_group(k, primitive)
            generic = generic_little_group(template, names, primitive)
            general = actual == generic
            if not general:
                issues.append(f'{parent}/{label}: enhanced symmetry at representative value')
            points[label] = k
            result['parametric_labels' if names else 'fixed_labels'][label] = dict(
                kovalev=kovalev, template=template, values=query['values'], coordinates=list(map(str, k)),
                little_group_order=len(actual), generic_little_group_order=len(generic),
                general_position=general, star_size=len(POINT_GROUP) // len(actual),
                frozen_archived_value=query['role'] == 'registered')
        for (left, k1), (right, k2) in combinations(points.items(), 2):
            pair = {'parent': parent, 'labels': [left, right]}
            if same_point_mod_lattice(k1, k2, primitive):
                report['label_aliases'].append(pair)
                issues.append(f'{parent}: reciprocal-equivalent labels {left}/{right}')
            if same_star(k1, k2, primitive):
                report['same_star_pairs'].append(pair)
                issues.append(f'{parent}: same-star labels {left}/{right}')
        report['parents'][parent] = result
    counts = Counter(q['role'] for q in contract['queries'])
    report['query_count'] = dict(fixed_labels=counts['fixed'], fixed_queries=2 * counts['fixed'],
        parametric_labels_total=counts['parametric'] + counts['registered'],
        parametric_labels_new=counts['parametric'], parametric_queries_new=2 * counts['parametric'],
        registered_queries=2 * counts['registered'],
        download_queries=2 * (counts['fixed'] + counts['parametric']))
    report['scope_limit'] = ('Generic little co-group on each stated affine template only; '
        'rational representatives, nmod 0/1. No physical-mode, nmod-equivalence, '
        'full Brillouin-zone or Method 3/4 completeness certification.')
    report['verdict'] = 'fail' if issues else 'pass'
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument('--report', type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args(argv)
    report = build_report(json.loads(args.manifest.read_text(encoding='utf-8')))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: report[key] for key in ('verdict', 'issues', 'query_count')}, ensure_ascii=False))
    return int(report['verdict'] != 'pass')


if __name__ == '__main__':
    raise SystemExit(main())
