"""MDC_PHASE — análisis de fase modular (Regla Mecánica / analisi_modular)."""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mdc_lib.analisi_modular import curvatura_modular, fase_modular, patron_bits
from runtime.plugin import Capability, Cost, IPlugin, PluginContext, PluginInfo, PluginType
from runtime.reference import Reference
from runtime.signature import make_signature


class MdcPhasePlugin(IPlugin):
    PLUGIN_ID = "MDC_PHASE"

    def info(self) -> PluginInfo:
        return PluginInfo(
            name=self.PLUGIN_ID,
            author="Víctor Manzanares Alberola",
            version=1,
            api_version=1,
            plugin_type=PluginType.INFER,
        )

    def capability(self) -> Capability:
        return Capability(deterministic=True, cacheable=True, distributable=True)

    def signature(self, ctx: PluginContext) -> str:
        n = int.from_bytes(ctx.input_payloads[0][:8], "little")
        m = int.from_bytes(ctx.input_payloads[0][8:16], "little")
        return make_signature("MDC_PHASE", [str(n), str(m)])

    def estimate(self, ctx: PluginContext) -> Cost:
        return Cost(alu_units=3)

    def validate(self, ctx: PluginContext) -> bool:
        return len(ctx.input_payloads) == 1 and len(ctx.input_payloads[0]) >= 16

    def execute(self, ctx: PluginContext) -> Reference:
        raw = ctx.input_payloads[0]
        n = int.from_bytes(raw[:8], "little")
        m = int.from_bytes(raw[8:16], "little")
        result = {
            "n": n,
            "m": m,
            "fase": fase_modular(n, m),
            "patron_bits": patron_bits(n, 8),
            "curvatura": curvatura_modular(n, m),
        }
        payload = json.dumps(result).encode("utf-8")
        return Reference.create(
            signature=self.signature(ctx),
            payload=payload,
            metadata={"plugin": self.PLUGIN_ID, **result},
        )