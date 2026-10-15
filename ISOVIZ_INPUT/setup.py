"""Build-time resource copy; editable settings remain in config/settings.yaml."""
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildWithSettings(build_py):
    """Include authoritative settings in wheels without a second source copy."""

    def run(self):
        super().run()
        source = Path(__file__).resolve().parent / "config" / "settings.yaml"
        destination = Path(self.build_lib) / "isoviz_input" / "_resources" / "settings.yaml"
        self.mkpath(str(destination.parent))
        self.copy_file(str(source), str(destination))


setup(cmdclass={"build_py": BuildWithSettings})
