"""Demo discriminante filtrada — delega en mdc_lib si AntiPC está en el path."""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1] / "antipc" / "src" / "antipc"
sys.path.insert(0, str(_ROOT))

from mdc_lib.discriminant import discriminant_factor, format_result  # noqa: E402


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 11449
    print(format_result(discriminant_factor(n)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())