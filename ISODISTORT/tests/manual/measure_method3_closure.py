#!/usr/bin/env python3
"""Retired unsafe scope measurement; retained only to reject old resume commands.

Original source and derived evidence: docs/history/ds-scope-before-correction-20261017.zip.
This command must not read official results or launch local computations.
"""
import sys


def main(argv=None):
    print("INVALIDATED: official-result-seeded closure and centering-blind categories "
          "cannot establish an independent download scope. Old checkpoints must not "
          "be resumed. See docs/DEVELOPMENT_PLAN.md for independent input enumeration.",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
