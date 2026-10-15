"""Public API for backend.inputs."""

from .documents import (
    TextDocument,
    document_from_text,
    load_text_document,
)

__all__ = ['TextDocument', 'document_from_text', 'load_text_document']
