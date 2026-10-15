"""Build-time resource copy from the authoritative config/settings.yaml."""
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildWithSettings(build_py):
    """Include settings in wheels without introducing a second source copy."""

    def run(self):
        super().run()
        source = Path(__file__).resolve().parent / "config" / "settings.yaml"
        destination = Path(self.build_lib) / "isodistort_validate" / "_resources" / "settings.yaml"
        self.mkpath(str(destination.parent))
        self.copy_file(str(source), str(destination))


setup(cmdclass={"build_py": BuildWithSettings})
