# file: src/__init__.py   (add this once)
import sys
from pathlib import Path

# Add the repository root to the Python path **once** when the package is first imported
REPO_ROOT = Path(__file__).resolve().parents[1]   # src → repo root
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))