"""Tests microkernel v0.1.2-alpha."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "antipc"))

from runtime.kop import KOPId, K3MicroKernel
from runtime.ledger import KnowledgeLedger
from runtime.microkernel import AntiPCMicroKernel


def test_boot():
    mk = AntiPCMicroKernel()
    mk.boot(verbose=False)
    assert mk.registry is not None


def test_create_and_reuse():
    mk = AntiPCMicroKernel()
    mk.boot(verbose=False)
    oid1 = mk.create_knowledge(b"test-payload", verbose=False)
    assert oid1 is not None
    oid2 = mk.create_knowledge(b"test-payload", verbose=False)
    assert mk.metrics.kop_reused >= 1 or oid1 == oid2


def test_ledger_k3mk():
    ledger = KnowledgeLedger()
    ledger.append("CREATE", knowledge_id="ko-1", producer="test")
    blobs = ledger.export_blobs()
    assert len(blobs) == 1
    mk = K3MicroKernel.unpack(blobs[0])
    assert KOPId.HISTORY in mk.blocks


def test_kop_registry():
    mk = AntiPCMicroKernel()
    mk.boot(verbose=False)
    desc = mk.executor.describe(KOPId.IDENTITY)
    assert desc["registered"] is True


def test_hash_kop_reuse():
    mk = AntiPCMicroKernel()
    mk.boot(verbose=False)
    d1 = mk.hash_payload(b"hash-test", verbose=False)
    d2 = mk.hash_payload(b"hash-test", verbose=False)
    assert d1 is not None and d1 == d2
    assert mk.metrics.kop_reused >= 1


def test_metrics_json_export(tmp_path=None):
    import tempfile
    from pathlib import Path

    mk = AntiPCMicroKernel()
    mk.boot(verbose=False)
    mk.hash_payload(b"x", verbose=False)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "metrics.json"
        mk.export_metrics(path)
        data = path.read_text(encoding="utf-8")
        assert "kop_breakdown" in data
        assert "p_util" in data


if __name__ == "__main__":
    test_boot()
    test_create_and_reuse()
    test_ledger_k3mk()
    test_kop_registry()
    test_hash_kop_reuse()
    test_metrics_json_export()
    print("OK: 6 tests passed")