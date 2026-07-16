"""Tests k3simil + k3search."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "src" / "antipc"
sys.path.insert(0, str(ROOT))

from k3_text_tools import jaccard, k3_search_index, k3_search_query, k3_simil_pairs


def test_jaccard_identical():
    a = [1, 2, 3]
    assert jaccard(a, a) == 1.0


def test_simil_and_search():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        t1 = d / "a.txt"
        t2 = d / "b.txt"
        t3 = d / "c.txt"
        base = "factorizacion de primos con metodo mdc y criba hibrida vma antipc"
        t1.write_text(base, encoding="utf-8")
        t2.write_text(base + " extra", encoding="utf-8")
        t3.write_text("completamente distinto sin palabras comunes", encoding="utf-8")

        pairs = k3_simil_pairs([t1, t2, t3], threshold=0.2)
        assert any(p["path_a"].endswith("a.txt") or p["path_b"].endswith("a.txt") for p in pairs)

        idx = d / "test.k3idx"
        meta = k3_search_index([t1, t2, t3], idx)
        assert meta["documents"] == 3
        hits = k3_search_query(idx, ["factorizacion", "primos"])
        assert len(hits) >= 1
        assert hits[0]["score"] > 0


if __name__ == "__main__":
    test_jaccard_identical()
    test_simil_and_search()
    print("OK")