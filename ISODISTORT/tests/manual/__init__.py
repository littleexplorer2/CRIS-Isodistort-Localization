"""Long-running manual audit and validation scripts (not collected by pytest).

These scripts are importable as :mod:`tests.manual.*` so unit tests can exercise
their parsers and helpers directly.  They must never modify read-only official
evidence; new outputs go to ``output/validation/``.
"""
