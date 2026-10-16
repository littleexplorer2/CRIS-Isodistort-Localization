"""Validate library references, acquisition hashes, module pointers and compiled PDFs."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, value):
    (ROOT / 'audit' / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def main():
    catalogue = json.loads((ROOT / 'sources.json').read_text(encoding='utf-8'))['sources']
    ids = {row['id'] for row in catalogue}
    assert len(ids) == len(catalogue), 'Duplicate bibliography ID'
    bib = (ROOT / 'latex/references.bib').read_text(encoding='utf-8')
    bib_ids = set(re.findall(r'@\w+\{([^,]+),', bib))
    assert ids == bib_ids
    doi_values = [row['doi'].lower() for row in catalogue if row.get('doi')]
    assert len(doi_values) == len(set(doi_values)), 'Duplicate work DOI'
    texts = {p.relative_to(ROOT).as_posix(): p.read_text(encoding='utf-8') for p in (ROOT / 'latex').glob('*.tex')}
    citations = set()
    for content in texts.values():
        for match in re.findall(r'\\cite\w*(?:\[[^]]*\])?\{([^}]+)\}', content):
            citations.update(match.split(','))
    assert citations <= ids, sorted(citations - ids)
    coverage = json.loads((ROOT / 'coverage.json').read_text(encoding='utf-8'))['capabilities']
    code_files = set()
    for row in coverage:
        assert set(row['source_ids']) <= ids, row['id']
        for proof in row['proofs']:
            assert '\\label{' + proof['label'] + '}' in texts[proof['file']], proof
        for relative in row['repository_modules']:
            path = (REPO / relative).resolve()
            assert path.is_relative_to(REPO) and path.exists(), relative
            code_files.update(path.rglob('*.py') if path.is_dir() else [path])
    implementation = json.loads((ROOT / 'audit/implementation_audit.json').read_text(encoding='utf-8'))
    reviewed = implementation['capabilities']
    assert len({row['id'] for row in reviewed}) == len(reviewed), 'Duplicate implementation capability'
    assert {row['id'] for row in reviewed} == {row['id'] for row in coverage}
    assert (ROOT / 'audit' / implementation['report']).is_file()
    for row in reviewed:
        assert set(row['source_ids']) <= ids, row['id']
        assert row['status'] and row['conclusion'], row['id']
        for proof in row['proofs']:
            assert '\\label{' + proof['label'] + '}' in texts[proof['file']], proof
        for relative in row['repository_modules']:
            path = (REPO / relative).resolve()
            assert path.is_relative_to(REPO) and path.exists(), relative
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    write_json('code_context.json', {'head': head, 'meaning': 'Module-location context, NOT a claim that these source files were fully audited or tested.', 'files': {p.relative_to(REPO).as_posix(): digest(p) for p in sorted(code_files)}})
    acquired = json.loads((ROOT / 'audit/download_manifest.json').read_text(encoding='utf-8'))
    assert {row['id'] for row in acquired} == ids, 'Acquisition manifest must account for every source'
    pdf_count, scanned, available = 0, [], 0
    for row in acquired:
        if row['status'] == 'unavailable':
            continue
        path = ROOT / row['local_file']
        assert path.exists() and digest(path) == row['sha256'], row['id']
        available += 1
        if path.suffix == '.pdf':
            pdf_count += 1
            assert path.read_bytes().startswith(b'%PDF-'), path
            if row.get('text_characters', 0) < 100:
                scanned.append(row['id'])
    extra = []
    for path in sorted((ROOT / 'snapshots').glob('*.txt')):
        extra.append({'file': path.relative_to(ROOT).as_posix(), 'sha256': digest(path), 'bytes': path.stat().st_size, 'provenance': 'Exa primary-page extraction or query excerpts; filename records retrieval topic.'})
    write_json('supplementary_acquisition.json', extra)
    pdfs = []
    for stem in ('unsupported_algorithm_proofs', 'isodistort_theoretical_foundations'):
        path = ROOT / 'latex' / (stem + '.pdf')
        log = (ROOT / 'latex/build' / (stem + '.log')).read_text(encoding='utf-8', errors='replace')
        problems = re.findall(r'(?m)^!.*|.*undefined.*|.*Overfull.*|.*Missing character.*', log)
        assert not problems, (stem, problems)
        blg = (ROOT / 'latex/build' / (stem + '.blg')).read_text(encoding='utf-8', errors='replace')
        assert 'Warning--' not in blg, stem + ': bibliography warning'
        assert digest(path) == digest(ROOT / 'latex/build' / (stem + '.pdf'))
        output = ROOT / 'audit/extracted' / (stem + '_compiled.txt')
        output.parent.mkdir(exist_ok=True)
        subprocess.run(['pdftotext', '-layout', str(path), str(output)], check=True)
        extracted = output.read_text(encoding='utf-8')
        assert '参考文献' in extracted and '证明' in extracted
        info = subprocess.check_output(['pdfinfo', str(path)], text=True, encoding='utf-8', errors='replace')
        pages = int(re.search(r'Pages:\s+(\d+)', info).group(1))
        pdfs.append({'file': path.relative_to(ROOT).as_posix(), 'pages': pages, 'sha256': digest(path), 'text_characters': len(extracted), 'citations_resolved': True})
    math_results = json.loads((ROOT / 'audit/mathematical_examples.json').read_text(encoding='utf-8'))
    assert math_results['status'] == 'passed_for_listed_examples'
    result = {'checked_at': datetime.now(timezone.utc).isoformat(), 'status': 'passed_document_integrity_checks', 'source_records': len(ids), 'cited_records': len(citations), 'capability_topics': len(coverage), 'downloaded_primary_files': available, 'reference_pdfs': pdf_count, 'scanned_pdfs_not_fully_text_reviewed': scanned, 'compiled_documents': pdfs, 'meaning': 'Checks document integrity, hashes, citations and selected examples; NOT full scientific-engine validation or formal verification. Visual review is recorded separately in DELIVERY_REVIEW.md.'}
    write_json('validation.json', result)
    print(json.dumps(result, ensure_ascii=True, indent=2))


if __name__ == '__main__':
    main()
