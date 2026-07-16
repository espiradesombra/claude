"""MDC_REGLA — Regla de cálculo analógica (planificador Battery Mode)."""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mdc_lib.regla_calculo import ReglaCalculo
from runtime.plugin import Capability, Cost, IPlugin, PluginContext, PluginInfo, PluginType
from runtime.reference import Reference
from runtime.signature import make_signature

_REGLA = ReglaCalculo(precision=500, tamano_escala=2000)


class MdcReglaPlugin(IPlugin):
    PLUGIN_ID = "MDC_REGLA"

    def info(self) -> PluginInfo:
        return PluginInfo(
            name=self.PLUGIN_ID,
            author="Víctor Manzanares Alberola",
            version=1,
            api_version=1,
            plugin_type=PluginType.TRANSFORM,
        )

    def capability(self) -> Capability:
        return Capability(deterministic=True, cacheable=True)

    def signature(self, ctx: PluginContext) -> str:
        a = int.from_bytes(ctx.input_payloads[0][:8], "little")
        b = int.from_bytes(ctx.input_payloads[0][8:16], "little")
        return make_signature("MDC_REGLA", [str(a), str(b)])

    def estimate(self, ctx: PluginContext) -> Cost:
        return Cost(alu_units=1)

    def validate(self, ctx: PluginContext) -> bool:
        return len(ctx.input_payloads) == 1 and len(ctx.input_payloads[0]) >= 16

    def execute(self, ctx: PluginContext) -> Reference:
        raw = ctx.input_payloads[0]
        a = max(1, int.from_bytes(raw[:8], "little")) / 1_000_000.0
        b = max(1, int.from_bytes(raw[8:16], "little")) / 1_000_000.0
        producto = _REGLA.multiplicar_valores(a, b)
        result = {"a": a, "b": b, "producto_regla": producto}
        return Reference.create(
            signature=self.signature(ctx),
            payload=json.dumps(result).encode(),
            metadata={"plugin": self.PLUGIN_ID, **result},
        )