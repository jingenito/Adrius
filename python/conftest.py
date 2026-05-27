# Copyright (c) 2025 InGenifold Research LLC. MIT License.
#
# Auto-discover the compiled _adrius extension in any build_*/python/ directory
# so pytest can be run from the repo root or python/ without setting PYTHONPATH.

import sys
from pathlib import Path


def _find_build_python() -> str | None:
    repo_root = Path(__file__).parent.parent
    for ext_file in sorted(repo_root.glob("build*/python/adrius/_adrius*")):
        if ext_file.suffix in (".pyd", ".so") or ".so." in ext_file.name:
            return str(ext_file.parent.parent)  # build_*/python/
    return None


_build_python = _find_build_python()
if _build_python and _build_python not in sys.path:
    sys.path.insert(0, _build_python)
