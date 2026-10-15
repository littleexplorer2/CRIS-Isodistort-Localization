"""Convert uploaded files or pasted text into one pipeline input type."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..paths import strip_user_path


@dataclass(frozen=True, slots=True)
class TextDocument:
    """Text plus its provenance; algorithms never need to know how it arrived."""

    text: str
    source: str
    path: Path | None = None

    @property
    def stem(self) -> str:
        return self.path.stem if self.path is not None else Path(self.source).stem


def document_from_text(text: str, *, source: str) -> TextDocument:
    """Create a document from pasted or otherwise in-memory text."""
    if not text.strip():
        raise ValueError(f"Input is empty: {source}")
    return TextDocument(text=text.lstrip("\ufeff"), source=source)


def load_text_document(
    raw_path: str | Path,
    *,
    expected_suffix: str,
    errors: str = "strict",
) -> TextDocument:
    """Read a user-selected text file and validate its filename contract."""
    path = Path(strip_user_path(str(raw_path))).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix.lower() != expected_suffix.lower():
        raise ValueError(f"Expected a {expected_suffix} file, got: {path}")
    text = path.read_text(encoding="utf-8-sig", errors=errors)
    if not text.strip():
        raise ValueError(f"Input file is empty: {path}")
    return TextDocument(text=text, source=str(path), path=path)
