"""Put playground/ on sys.path so tests can `import pcca`."""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "playground"))
sys.path.insert(0, str(Path(__file__).parent))
