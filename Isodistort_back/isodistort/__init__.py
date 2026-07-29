"""isodistort public package.

This module keeps the new package name stable while re-exporting the
implementation that still lives under :mod:`isocore`.
"""

from importlib import import_module
import sys

from isocore import __version__


_ALIASES = {
    "isodistort.api": "isocore.api",
    "isodistort.api.core_api": "isocore.api.core_api",
    "isodistort.backend": "isocore.backend",
    "isodistort.backend.base_wrapper": "isocore.backend.base_wrapper",
    "isodistort.backend.findsym_wrapper": "isocore.backend.findsym_wrapper",
    "isodistort.backend.iso_wrapper": "isocore.backend.iso_wrapper",
    "isodistort.distortion": "isocore.distortion",
    "isodistort.distortion.phase_path": "isocore.distortion.phase_path",
    "isodistort.distortion.distortion_mapper": "isocore.distortion.distortion_mapper",
    "isodistort.distortion.distortion_engine": "isocore.distortion.distortion_engine",
    "isodistort.distortion.domain_generator": "isocore.distortion.domain_generator",
    "isodistort.io": "isocore.io",
    "isodistort.io.structure_exporter": "isocore.io.structure_exporter",
    "isodistort.io.result_serializer": "isocore.io.result_serializer",
    "isodistort.structure": "isocore.structure",
    "isodistort.structure.cif_io": "isocore.structure.cif_io",
    "isodistort.structure.coordinate_transform": "isocore.structure.coordinate_transform",
    "isodistort.structure.site_mapping": "isocore.structure.site_mapping",
    "isodistort.structure.symmetry_validator": "isocore.structure.symmetry_validator",
    "isodistort.utils": "isocore.utils",
    "isodistort.utils.config_loader": "isocore.utils.config_loader",
    "isodistort.utils.exceptions": "isocore.utils.exceptions",
    "isodistort.utils.text_parser": "isocore.utils.text_parser",
}

for alias, target in _ALIASES.items():
    sys.modules.setdefault(alias, import_module(target))
