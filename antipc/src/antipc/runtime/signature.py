"""S = f(Operation, Inputs) — transformation signatures."""

from __future__ import annotations

import hashlib
import json
from typing import Any


def make_signature(operation: str, inputs: list[str], version: int = 1) -> str:
    """Deterministic signature independent of node, time, or host."""
    body = json.dumps(
        {"op": operation, "inputs": sorted(inputs), "v": version},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(body.encode()).hexdigest()[:32]