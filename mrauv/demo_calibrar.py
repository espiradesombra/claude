"""Demo MRAUV filtrada — delega en mdc_lib si AntiPC está en el path."""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1] / "antipc" / "src" / "antipc"
sys.path.insert(0, str(_ROOT))

from mdc_lib.mrauv import calibrar, format_calibrar  # noqa: E402


def main() -> int:
    n0 = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
    print(format_calibrar(calibrar(n0)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())