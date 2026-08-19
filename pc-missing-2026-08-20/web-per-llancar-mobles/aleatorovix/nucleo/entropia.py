"""Fuentes de entropía física — nanosegundos, pila, basura, red opcional."""
from __future__ import annotations

import secrets
import socket
import time
from dataclasses import dataclass


@dataclass
class MuestraEntropia:
    nanos: int
    inercia_a: int
    red_x: int
    bit_pila: int
    basura: int
    ping_us: int | None


class EntropiaFisica:
    """Roba señal del entorno como en organismo_lila_v99.c (port Windows)."""

    def __init__(
        self,
        ping_host: str = "8.8.8.8",
        ping_timeout_ms: int = 80,
        usar_ping: bool = True,
    ) -> None:
        self.ping_host = ping_host
        self.ping_timeout_ms = ping_timeout_ms
        self.usar_ping = usar_ping

    @staticmethod
    def get_nanos() -> int:
        return time.perf_counter_ns() % 1_000_000_000

    @staticmethod
    def robar_bit_pila() -> int:
        variable_local = 0
        sp = id(variable_local)
        return (sp >> 6) & 0x01

    def robar_basura_memoria(self) -> int:
        basura = secrets.randbits(32)
        puntero_sucio = id(basura) ^ self.get_nanos()
        return puntero_sucio & 0xFF

    def medir_ping_us(self) -> int | None:
        try:
            t0 = time.perf_counter_ns()
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(self.ping_timeout_ms / 1000.0)
            sock.connect((self.ping_host, 53))
            sock.close()
            return (time.perf_counter_ns() - t0) // 1000
        except OSError:
            return None

    def jitter(self, inercia_a: int) -> None:
        micro = (inercia_a % 97) / 100_000.0
        if micro > 0:
            time.sleep(micro)

    def muestrear(self) -> MuestraEntropia:
        nanos = self.get_nanos()
        inercia_a = nanos % 1000
        self.jitter(inercia_a)
        ping = self.medir_ping_us() if self.usar_ping else None
        x_fuente = ping if ping and ping > 0 else self.get_nanos()
        red_x = x_fuente % 1000
        return MuestraEntropia(
            nanos=nanos,
            inercia_a=inercia_a,
            red_x=red_x,
            bit_pila=self.robar_bit_pila(),
            basura=self.robar_basura_memoria(),
            ping_us=ping,
        )