import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class CorpusReleaseTest(unittest.TestCase):
    def test_release_validator(self) -> None:
        completed = subprocess.run(
            [str(ROOT / ".venv/bin/python"), str(ROOT / "scripts/validate_corpus.py")],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn("120 total records", completed.stdout)


if __name__ == "__main__":
    unittest.main()
