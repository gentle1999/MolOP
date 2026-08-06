from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def test_molop_import_loads_rdkit_before_openbabel() -> None:
    source_root = Path(__file__).resolve().parents[1] / "src"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (str(source_root), environment.get("PYTHONPATH")))
    )

    result = subprocess.run(
        [sys.executable, "-c", "import molop; import openbabel.pybel"],
        cwd=source_root.parent,
        env=environment,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
