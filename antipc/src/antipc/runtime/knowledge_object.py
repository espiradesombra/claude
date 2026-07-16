"""KnowledgeObject — primer objeto del microkernel con identidad K3."""

from __future__ import annotations

import uuid
from dataclasses import dataclass, field

from .kop import build_knowledge_blob
from .reference import Reference
from .registry import Descriptor, IRegistryObject, ObjectID
from .signature import make_signature


@dataclass
class KnowledgeObject(IRegistryObject):
    """Payload + Identity + DNA empaquetados en K3MK."""

    ref: Reference
    producer: str = "antipc"
    trust: float = 1.0
    _blob: bytes | None = field(default=None, repr=False)

    @classmethod
    def create(cls, payload: bytes, label: str = "RAW", producer: str = "antipc") -> KnowledgeObject:
        sig = make_signature("IDENTITY", [payload.hex(), label])
        ref_id = f"ko-{uuid.uuid4().hex[:12]}"
        ref = Reference(
            ref_id=ref_id,
            signature=sig,
            payload=payload,
            metadata={"label": label, "producer": producer},
        )
        return cls(ref=ref, producer=producer)

    def id(self) -> ObjectID:
        return self.ref.ref_id

    def descriptor(self) -> Descriptor:
        return Descriptor(
            kind="KnowledgeObject",
            name=self.ref.metadata.get("label", "object"),
            capabilities=("identity", "signature", "reuse"),
        )

    def k3_blob(self) -> bytes:
        if self._blob is None:
            self._blob = build_knowledge_blob(
                self.ref.ref_id,
                self.ref.signature,
                self.ref.payload or b"",
                producer=self.producer,
                trust=self.trust,
            )
        return self._blob