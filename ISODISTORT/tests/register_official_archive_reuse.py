#!/usr/bin/env python3
"""Register the reusable part of the existing official archive (read-only).

The user decided on 2026-10-17 that new official downloads must live in the
existing ``output_compare`` tree and that already archived official files must
be reused instead of downloaded again.  This tool keeps that decision
auditable:

* it only READS ``output_compare/`` (never writes there);
* it locates every archived Method 2 candidate directory by its
  ``subgroup.cif`` file and verifies the four core artifacts are present and
  non-empty;
* it reads the submitted-input provenance that the official page wrote into
  each candidate (``# k point:`` header line and ``!numberOfModulations``);
* it fingerprints every official method subtree and every archived Method 2
  candidate, so the reuse registry can later prove the archive did not change;
* it writes one JSON registry OUTSIDE the archive (default:
  ``ISODISTORT/docs/manifests/official_archive_reuse_known_method2.json``).

The script refuses to run if the observed candidate count, k point or nmod
value disagrees with the registered contract: a changed archive must be
investigated, not silently re-registered.  No official result is copied into
any computation input and no candidate identity is inferred from the output;
the candidate directory names are recorded as the official page printed them.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parent
DEFAULT_ARCHIVE_ROOT = REPO_ROOT / 'output_compare'
DEFAULT_OUTPUT = PACKAGE_ROOT / 'docs/manifests/official_archive_reuse_known_method2.json'

CORE_FILES = ('subgroup.cif', 'data.isoviz', 'topas.str')
FINGERPRINT_ALGORITHM = (
    'For each file, sorted by POSIX relative path: sha256(relative_path UTF-8, NUL, '
    'decimal size, NUL, sha256(file bytes)). Candidate and subtree aggregates sort '
    'the same way; directory aggregates are not used.')

# Registered Method 2 reuse contract: query id, parent, archived layout, expected
# candidate count, official k-point header and archived nmod value.  These facts
# are checked against the archive; a mismatch aborts the registration.
METHOD2_REUSE_CONTRACT = (
    {
        'query_id': 'M2P-EU-LD',
        'parent': 'EuAl4 Parent.cif',
        'nmod': 0,
        'expected_candidates': 22,
        'archive_subdir': 'LD_g1-6/nmod0',
        'layout': '<K>/nmod<N>/<IR>_<OPD>',
        'expected_k_header': 'LD, k10 (0,0,g), g=1/6',
        'registered_input_evidence': (
            'independent input values reconfirmed in the 20261019 download manifest '
            '(k template (0,0,g), g=1/6); this archive predates it and is matched '
            'only by the official page header'),
        'submit_state_capture': None,
        'action': 'register_only_no_download',
    },
    {
        'query_id': 'M2P-ND-Y',
        'parent': 'NdNiO2 own.cif',
        'nmod': 0,
        'expected_candidates': 48,
        'archive_subdir': 'Y_a1-3/nmod0',
        'layout': '<K>/nmod<N>/<IR>_<OPD>',
        'expected_k_header': 'Y, k6 (a,1/2,0), a=1/3',
        'registered_input_evidence': (
            'independent input values reconfirmed in the 20261019 download manifest '
            '(k template (a,1/2,0), a=1/3); this archive predates it and is '
            'matched only by the official page header'),
        'submit_state_capture': None,
        'action': 'register_only_no_download',
    },
    {
        'query_id': 'M2F-C4310-X',
        'parent': '4310_tetra.cif',
        'nmod': 0,
        'expected_candidates': 24,
        'archive_subdir': 'X/nmod0',
        'layout': 'flat <IR>_<OPD> under Method2/X/nmod0 (built from archived result lists)',
        'expected_k_header': 'X, k13 (1/2,1/2,0)',
        'registered_input_evidence': (
            'submit-state screenshot in the archive shows k point X, k13 '
            '(1/2,1/2,0) and "# of independent incommensurate modulations" = 0'),
        'submit_state_capture': 'Screenshot 2026-10-05 222306.png',
        'action': 'register_built_structure',
        'list_evidence': (
            'webpage_info/4310_tetra.cif/3. X_nmod0 - ... plus the eight '
            '4.1-4.8 X_nmod0_<IR> result-list pages; 24 rows verified against '
            'these candidate folders in '
            'output/validation/rebuild_20261019/list_identity_checks.json'),
    },
)

# Labels whose Method 2 nmod=1 output was never archived: the registered archive
# stores only nmod=0 candidates (verified below from !numberOfModulations).
METHOD2_NMOD1_CONTRACT = (
    {
        'query_id': 'M2P-EU-LD',
        'parent': 'EuAl4 Parent.cif',
        'nmod': 1,
        'status': 'absent_conditional_candidate_reuse',
        'rule': (
            'download the nmod=1 official result-list pages first; candidate files may '
            'be reused only if every candidate identity and every core file byte matches '
            'the registered nmod=0 archive, otherwise save them'),
    },
    {
        'query_id': 'M2P-ND-Y',
        'parent': 'NdNiO2 own.cif',
        'nmod': 1,
        'status': 'absent_conditional_candidate_reuse',
        'rule': (
            'download the nmod=1 official result-list pages first; candidate files may '
            'be reused only if every candidate identity and every core file byte matches '
            'the registered nmod=0 archive, otherwise save them'),
    },
    {
        'query_id': 'M2F-C4310-X',
        'parent': '4310_tetra.cif',
        'nmod': 1,
        'status': 'absent_must_download',
        'rule': 'no archived nmod=1 candidate set exists for this label',
    },
)


class RegistrationError(RuntimeError):
    """The archive does not match the frozen reuse contract."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _posix(path: Path, base: Path) -> str:
    return path.relative_to(base).as_posix()


def fingerprint_files(paths, base: Path) -> dict:
    """Aggregate fingerprint over files given as paths below ``base``."""
    entries = {}
    for path in paths:
        rel = _posix(path, base)
        entries[rel] = (path.stat().st_size, sha256_file(path))
    digest = hashlib.sha256()
    for rel in sorted(entries):
        size, sha = entries[rel]
        digest.update(rel.encode('utf-8'))
        digest.update(b'\0')
        digest.update(str(size).encode('ascii'))
        digest.update(b'\0')
        digest.update(bytes.fromhex(sha))
    return {
        'file_count': len(entries),
        'total_bytes': sum(size for size, _sha in entries.values()),
        'aggregate_sha256': digest.hexdigest(),
    }


def fingerprint_tree(root: Path) -> dict:
    return fingerprint_files([p for p in root.rglob('*') if p.is_file()], root)


def find_candidates(archive_dir: Path) -> list[Path]:
    """Locate candidate directories directly below one registered archive dir.

    The contract names the exact archive location of each reused set (the legacy
    K/nmod query directory for each parent), so the scan
    never picks up candidates from other queries that are added later.
    """
    return sorted(
        (path for path in archive_dir.iterdir()
         if path.is_dir() and (path / 'subgroup.cif').is_file()),
        key=lambda path: path.as_posix())


def candidate_id(candidate_dir: Path, method2_root: Path) -> str:
    return candidate_dir.relative_to(method2_root).as_posix()


def read_k_header(candidate_dir: Path) -> str:
    for line in (candidate_dir / 'subgroup.cif').read_text(
            encoding='utf-8', errors='replace').splitlines()[:30]:
        if line.startswith('# k point:'):
            return line[len('# k point:'):].strip()
    return ''


def read_number_of_modulations(candidate_dir: Path) -> int:
    text = (candidate_dir / 'data.isoviz').read_text(
        encoding='utf-8', errors='replace')
    match = re.search(r'^!numberOfModulations\s+(\d+)\s*$', text, re.MULTILINE)
    if match is None:
        raise RegistrationError(f'{candidate_dir}: no !numberOfModulations value')
    return int(match.group(1))


def candidate_resource_dirs(candidate_dir: Path) -> list[dict]:
    resources = []
    for html in sorted(candidate_dir.glob('*.htm*')):
        companion = html.with_name(html.stem + '_files')
        if companion.is_dir():
            resources.append({'path': companion.name, **fingerprint_tree(companion)})
    return resources


def _recorded_path(candidate_dir: Path, method2_root: Path) -> str:
    """Prefer a repository-relative path; fall back for temporary test trees."""
    for base in (REPO_ROOT, *(method2_root.parents[:3])):
        try:
            return _posix(candidate_dir, base)
        except ValueError:
            continue
    return candidate_dir.as_posix()


def inventory_candidate(candidate_dir: Path, method2_root: Path,
                        nmod: int | None = None) -> dict:
    files = {}
    missing = []
    for name in CORE_FILES:
        path = candidate_dir / name
        if path.is_file() and path.stat().st_size > 0:
            files[name] = {'bytes': path.stat().st_size, 'sha256': sha256_file(path)}
        else:
            missing.append(name)
    detail_pages = sorted(candidate_dir.glob('*.htm*'))
    if not detail_pages:
        missing.append('complete_modes_details.html')
    for page in detail_pages:
        if page.is_file() and page.stat().st_size > 0:
            files[page.name] = {'bytes': page.stat().st_size, 'sha256': sha256_file(page)}
        else:
            missing.append(page.name)
    extra = [path for path in candidate_dir.rglob('*')
             if path.is_file() and path.name not in files
             and path.parent == candidate_dir]
    entry = {
        'id': candidate_id(candidate_dir, method2_root),
        'path': _recorded_path(candidate_dir, method2_root),
        'missing_files': missing,
        'files': files,
        'resource_dirs': candidate_resource_dirs(candidate_dir),
        'k_header': read_k_header(candidate_dir),
        'number_of_modulations': read_number_of_modulations(candidate_dir),
        'aggregate_sha256': fingerprint_files(
            [path for path in candidate_dir.rglob('*') if path.is_file()],
            candidate_dir)['aggregate_sha256'],
    }
    if extra:
        entry['extra_files'] = sorted(path.name for path in extra)
    if nmod is not None and entry['number_of_modulations'] != nmod:
        raise RegistrationError(
            f'{candidate_dir}: !numberOfModulations={entry["number_of_modulations"]} '
            f'but the registered query has nmod={nmod}')
    return entry


def method2_registration(contract: dict, archive_root: Path) -> dict:
    parent = contract['parent']
    method2_root = archive_root / parent / '官网' / 'Method2'
    subdir = contract.get('archive_subdir', '')
    archive_dir = method2_root / subdir if subdir else method2_root
    if not archive_dir.is_dir():
        raise RegistrationError(f'missing archive subtree: {archive_dir}')
    candidate_dirs = find_candidates(archive_dir)
    if len(candidate_dirs) != contract['expected_candidates']:
        raise RegistrationError(
            f'{parent}: found {len(candidate_dirs)} archived Method 2 candidates in '
            f'{archive_dir}, contract expects {contract["expected_candidates"]}')
    candidates = [inventory_candidate(path, archive_dir, nmod=contract['nmod'])
                  for path in candidate_dirs]
    k_headers = {entry['k_header'] for entry in candidates}
    if k_headers != {contract['expected_k_header']}:
        raise RegistrationError(
            f'{parent}: archived k headers {sorted(k_headers)} disagree with '
            f'the registered input {contract["expected_k_header"]!r}')
    missing = [(entry['id'], name) for entry in candidates
               for name in entry['missing_files']]
    if missing:
        raise RegistrationError(f'{parent}: incomplete reused candidates: {missing[:5]}')
    entry = dict(contract)
    entry['archive_path'] = f'{parent}/官网/Method2' + (f'/{subdir}' if subdir else '')
    entry['candidate_count'] = len(candidates)
    entry['candidate_nmod_values'] = {str(contract['nmod']): len(candidates)}
    entry['k_headers'] = sorted(k_headers)
    candidate_files = [path for directory in candidate_dirs
                       for path in directory.rglob('*') if path.is_file()]
    entry['candidates_fingerprint'] = fingerprint_files(candidate_files, archive_dir)
    entry['archive_top_level_files'] = fingerprint_files(
        [path for path in archive_dir.iterdir() if path.is_file()], archive_dir)
    capture = contract.get('submit_state_capture')
    if capture:
        capture_path = archive_dir / capture
        if not capture_path.is_file():
            raise RegistrationError(
                f'{parent}: registered submit-state capture {capture!r} is missing')
        entry['submit_state_capture'] = {
            'path': _posix(capture_path, REPO_ROOT),
            'bytes': capture_path.stat().st_size,
            'sha256': sha256_file(capture_path),
        }
    entry['candidates'] = candidates
    return entry


def official_method_inventory(archive_root: Path, parents: list[str]) -> dict:
    inventory = {}
    for parent in parents:
        for method in ('Method1', 'Method2', 'Method3', 'Method4'):
            root = archive_root / parent / '官网' / method
            if not root.is_dir():
                continue
            files = [path for path in root.rglob('*') if path.is_file()]
            inventory[f'{parent}/官网/{method}'] = {
                'candidate_or_case_dirs': len([path for path in root.rglob('subgroup.cif')]),
                'top_level_dirs': len([path for path in root.iterdir() if path.is_dir()]),
                'file_count': len(files),
                'total_bytes': sum(path.stat().st_size for path in files),
                'aggregate_sha256': fingerprint_tree(root)['aggregate_sha256'],
            }
    return inventory


def build_registry(archive_root: Path = DEFAULT_ARCHIVE_ROOT,
                   generated_at: str | None = None) -> dict:
    archive_root = Path(archive_root)
    if not archive_root.is_dir():
        raise RegistrationError(f'official archive root not found: {archive_root}')
    parents = sorted({contract['parent'] for contract in METHOD2_REUSE_CONTRACT})
    method2_entries = [method2_registration(contract, archive_root)
                       for contract in METHOD2_REUSE_CONTRACT]
    registry = {
        'schema': 'official-archive-reuse/v1',
        'status': 'registered',
        'generated_at': generated_at or datetime.now(timezone.utc).isoformat(),
        'generator': {
            'path': _posix(Path(__file__).resolve(), REPO_ROOT),
            'sha256': sha256_file(Path(__file__).resolve()),
            'command': (
                '.venv/Scripts/python.exe ISODISTORT/tests/'
                'register_official_archive_reuse.py'),
        },
        'archive_root': _posix(archive_root, REPO_ROOT),
        'fingerprint_algorithm': FINGERPRINT_ALGORITHM,
        'write_policy': {
            'archive_read_only': (
                'Existing official files under output_compare are verifiable '
                'reference material: never overwrite, never delete, never rerun a '
                'query to replace them, and never copy them into a computation input.'),
            'new_downloads': (
                'New official page downloads (result lists) go into the registered '
                'webpage_info/<parent>/ page archive using the documented naming; new '
                'official exports may only be written into Method2/<K>/nmod<N>/ folders '
                'that were generated from a parsed official result list. Unknown lists '
                'stay blank until downloaded; no local-web-version output may be used '
                'to infer them.'),
            'reuse': (
                'Before every download, look the query up in this registry and in the '
                'webpage_info page index. A registered candidate set or already archived '
                'official page is reused by reference; it is not downloaded again and not '
                'copied into a computation input.'),
            'fingerprints': (
                'An audit re-fingerprints the archive and fails closed if a '
                'registered subtree or candidate aggregate changes.'),
        },
        'method2_registered_reuse': method2_entries,
        'method2_nmod1_status': [dict(entry) for entry in METHOD2_NMOD1_CONTRACT],
        'other_methods': official_method_inventory(archive_root, parents),
        'scope_notes': [
            'Only Method 2 nmod=0 candidate sets are archived for the reusable k points; '
            'no nmod=1 candidate set exists in the archive.',
            'The 4310 X/nmod=0 candidate set was moved into the official query folder '
            'Method2/X/nmod0 after verifying all 24 rows against the archived official '
            'result-list pages; EuAl4 LD and NdNiO2 Y files are now organized under '
            'their own K/nmod0 folders, with missing result lists tracked separately.',
            'Method 1 (330 candidates in total) stays registered-only and is not '
            'downloaded again; its complete official result lists live in the read-only '
            'webpage_info page archive and are referenced, not copied.',
            'Method 3/4 archive cases are historical evidence for the superseded '
            'scope. Previous Method 3/4 freezes are revoked; independent inputs remain '
            'unresolved. This helper checks only the three known Method 2 sets; '
            'the current complete inventory is official_archive_reuse_20261019.json.',
            'The EuAl4 Method 4 G05 auto-origin basis supplement is already '
            'delivered under output_compare and needs no further download.',
        ],
    }
    registry['official_archive_snapshot'] = {
        'path': _posix(archive_root, REPO_ROOT) + '/<parent>/官网',
        'parents': parents,
        'subtree_count': len(registry['other_methods']),
        'file_count': sum(row['file_count'] for row in registry['other_methods'].values()),
        'total_bytes': sum(row['total_bytes'] for row in registry['other_methods'].values()),
        'aggregate_sha256': _combine_subtree_aggregates(registry['other_methods']),
        'note': (
            'Snapshot at registration time.  Method 2 (and later Method 3/4) '
            'subtrees may grow only through the pre-registered new download '
            'directories; the immutable part is method2_registered_reuse, whose '
            'candidate fingerprints must not change.'),
    }
    registry['verdict'] = 'pass'
    return registry


def _combine_subtree_aggregates(subtrees: dict) -> str:
    digest = hashlib.sha256()
    for key in sorted(subtrees):
        row = subtrees[key]
        digest.update(key.encode('utf-8'))
        digest.update(b'\0')
        digest.update(row['aggregate_sha256'].encode('ascii'))
    return digest.hexdigest()


def _ensure_outside_archive(output: Path, archive_root: Path) -> None:
    output = output.resolve()
    archive_root = archive_root.resolve()
    if output == archive_root or archive_root in output.parents:
        raise RegistrationError(
            'the reuse registry must be written outside the read-only official archive')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-root', type=Path, default=DEFAULT_ARCHIVE_ROOT)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    try:
        _ensure_outside_archive(args.output, args.archive_root)
        registry = build_registry(args.archive_root)
    except RegistrationError as error:
        print(f'INVALID reuse registration: {error}', file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(registry, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    summary = {
        'verdict': registry['verdict'],
        'reused_queries': [
            {'query_id': row['query_id'], 'candidates': row['candidate_count'],
             'nmod': row['nmod'], 'action': row['action']}
            for row in registry['method2_registered_reuse']],
        'archive_file_count': registry['official_archive_snapshot']['file_count'],
        'archive_total_bytes': registry['official_archive_snapshot']['total_bytes'],
        'registry': _posix(args.output.resolve(), REPO_ROOT),
    }
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
