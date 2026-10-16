"""Fetch only the explicit bibliography; keep all artifacts inside this library."""
from __future__ import annotations
import concurrent.futures
import argparse
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


class TextParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []
        self.skip = 0
    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style'): self.skip += 1
        if tag in ('p', 'br', 'div', 'h1', 'h2', 'h3', 'li', 'tr'): self.parts.append('\n')
    def handle_endtag(self, tag):
        if tag in ('script', 'style'): self.skip = max(0, self.skip - 1)
    def handle_data(self, data):
        if not self.skip: self.parts.append(data)


def fetch(row):
    target = ROOT / row['local_file']
    assert target.resolve().is_relative_to(ROOT)
    try:
        if target.exists():
            raw, final_url, status = target.read_bytes(), row['url'], 'cached_local'
        else:
            request = urllib.request.Request(row['url'], headers={'User-Agent': 'CRIS-reference-library/1.0'})
            with urllib.request.urlopen(request, timeout=35) as response:
                raw, final_url = response.read(), response.url
            if target.suffix == '.pdf' and not raw.startswith(b'%PDF-'):
                raise ValueError('not a PDF response; no file saved')
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
            status = 'downloaded'
        text_path = ROOT / 'audit/extracted' / (row['id'] + '.txt')
        text_path.parent.mkdir(parents=True, exist_ok=True)
        if target.suffix == '.pdf':
            subprocess.run(['pdftotext', '-layout', str(target), str(text_path)], check=True, capture_output=True, timeout=30)
            text = text_path.read_text(encoding='utf-8', errors='replace')
        else:
            parser = TextParser()
            parser.feed(raw.decode('utf-8', errors='replace'))
            text = '\n'.join(x.strip() for x in ''.join(parser.parts).splitlines() if x.strip())
            text_path.write_text(text, encoding='utf-8')
        return {'id': row['id'], 'status': status, 'url': row['url'], 'final_url': final_url,
                'local_file': row['local_file'], 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(),
                'checked_at': datetime.now(timezone.utc).isoformat(), 'text_characters': len(text),
                'doi_candidates': list(dict.fromkeys(re.findall(r'10\.\d{4,9}/[^\s<>]+', text[:14000], re.I)))}
    except Exception as error:
        return {'id': row['id'], 'status': 'unavailable', 'url': row['url'], 'error': str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ids', nargs='+', help='Only retrieve these catalogue IDs; preserve other acquisition records.')
    args = parser.parse_args()
    rows = json.loads((ROOT / 'sources.json').read_text(encoding='utf-8'))['sources']
    if args.ids:
        unknown = set(args.ids) - {row['id'] for row in rows}
        if unknown:
            parser.error('Unknown catalogue IDs: ' + ', '.join(sorted(unknown)))
        rows = [row for row in rows if row['id'] in args.ids]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        results = list(pool.map(fetch, rows))
    manifest = ROOT / 'audit/download_manifest.json'
    if args.ids and manifest.exists():
        previous = json.loads(manifest.read_text(encoding='utf-8'))
        replaced = {row['id'] for row in results}
        if 'Wagner2009' in replaced:
            replaced.add('Wagner2011')
        merged = [row for row in previous if row['id'] not in replaced] + results
    else:
        merged = results
    manifest.write_text(json.dumps(merged, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    for row in results: print(json.dumps(row, ensure_ascii=True))


if __name__ == '__main__':
    main()
