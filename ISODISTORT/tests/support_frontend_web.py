"""Tests for support_frontend_web; migrated without changing scientific expectations."""
from __future__ import annotations

import copy


# Original coverage: test_web.py.

class _SessionIsoStub:
    """Small mutable IsoDistort stand-in for revision/session tests."""

    def __init__(self) -> None:
        self.structure = None
        self.subgroups = [object()]
        self.mode_displacements = {"old": object()}
        self.mode_occupancies = {"old-occ": object()}
        self.distorted_structure = [object()]
        self.scope = None
        self.types = None
        self.selected_modes_cleared = False

    def snapshot_for_action(self):
        return copy.deepcopy(self)

    def close(self):
        pass

    @staticmethod
    def species():
        return []

    def set_distortion_scope(self, scope):
        self.scope = scope

    def set_distortion_types(self, types):
        self.types = list(types)

    def set_subgroup_candidates(self, subgroups):
        self.subgroups = list(subgroups)

    def clear_selected_modes(self):
        self.selected_modes_cleared = True
        self.mode_displacements = {}
        self.mode_occupancies = {}
        self.distorted_structure = None
