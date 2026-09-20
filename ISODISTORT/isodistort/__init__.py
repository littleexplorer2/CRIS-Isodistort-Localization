"""isodistort public package.

This module keeps the new package name stable while re-exporting the
implementation that still lives under :mod:`isocore`.
"""

import sys
from importlib import import_module

from isocore import __version__
from isocore.api import IsoDistort

__all__ = ["IsoDistort", "__version__"]

_ALIASES = {
    "isodistort.api": "isocore.api",
    "isodistort.api.core_api": "isocore.api.core_api",
    "isodistort.backend": "isocore.backend",
    "isodistort.backend.base_wrapper": "isocore.backend.base_wrapper",
    "isodistort.backend.findsym_wrapper": "isocore.backend.findsym_wrapper",
    "isodistort.backend.isotropy_cache": "isocore.backend.isotropy_cache",
    "isodistort.backend.iso_wrapper": "isocore.backend.iso_wrapper",
    "isodistort.backend.smodes_wrapper": "isocore.backend.smodes_wrapper",
    "isodistort.data": "isocore.data",
    "isodistort.data.irreps_cdml": "isocore.data.irreps_cdml",
    "isodistort.data.kpoints_official": "isocore.data.kpoints_official",
    "isodistort.data.spacegroup_meta": "isocore.data.spacegroup_meta",
    "isodistort.distortion": "isocore.distortion",
    "isodistort.distortion.domain_generator": "isocore.distortion.domain_generator",
    "isodistort.distortion.distortion_engine": "isocore.distortion.distortion_engine",
    "isodistort.distortion.distortion_mapper": "isocore.distortion.distortion_mapper",
    "isodistort.distortion.occupational_modes": "isocore.distortion.occupational_modes",
    "isodistort.distortion.phase_path": "isocore.distortion.phase_path",
    "isodistort.distortion.search_methods": "isocore.distortion.search_methods",
    "isodistort.io": "isocore.io",
    "isodistort.io.distortion_formats": "isocore.io.distortion_formats",
    "isodistort.io.isodistort_cif": "isocore.io.isodistort_cif",
    "isodistort.io.isodistort_isoviz": "isocore.io.isodistort_isoviz",
    "isodistort.io.result_serializer": "isocore.io.result_serializer",
    "isodistort.io.structure_exporter": "isocore.io.structure_exporter",
    "isodistort.i18n": "isocore.i18n",
    "isodistort.structure": "isocore.structure",
    "isodistort.structure.cif_io": "isocore.structure.cif_io",
    "isodistort.structure.coordinate_transform": "isocore.structure.coordinate_transform",
    "isodistort.structure.symmetry_validator": "isocore.structure.symmetry_validator",
    "isodistort.utils": "isocore.utils",
    "isodistort.utils.config_loader": "isocore.utils.config_loader",
    "isodistort.utils.exceptions": "isocore.utils.exceptions",
    "isodistort.utils.opd_format": "isocore.utils.opd_format",
    "isodistort.utils.parent_header": "isocore.utils.parent_header",
    "isodistort.utils.schoenflies": "isocore.utils.schoenflies",
    "isodistort.utils.self_check": "isocore.utils.self_check",
    "isodistort.utils.text_parser": "isocore.utils.text_parser",
}

for alias, target in _ALIASES.items():
    sys.modules.setdefault(alias, import_module(target))
