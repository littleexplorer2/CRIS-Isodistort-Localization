from pathlib import Path
import tempfile
import unittest

from cris_rea.backend.runtime import _output_diagnostics


class NativeBombTests(unittest.TestCase):
    def test_prompt_prefixed_fatal_diagnostic(self):
        # Actual ISO 9.6.1 SETTING INTERNATIONAL NEW failure, despite exit code 0.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "stdout.txt").write_text(
                "***This program has bombed: isotropy_setting_commands 84\n")
            self.assertEqual(_output_diagnostics(path)[0]["code"], "program_bombed")


if __name__ == "__main__":
    unittest.main()
