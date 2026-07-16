#!/usr/bin/env python3
"""
Aleatorovix v1.0 — organismo unificado Lila + MDC + criba desmemoriada.
Víctor Manzanares Alberola · VMA / 33×1

Ciclo: entropía → máscara → intérprete mutante → acción → olvido.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from nucleo.acciones import ACCION_NOMBRES, AccionesLila
from nucleo.entropia import EntropiaFisica
from nucleo.mascara_lila import MascaraLila


class AleatorovixOrganismo:
    """Organismo Lila v1.0 — respira, decide, actúa y olvida."""

    def __init__(
        self,
        n_mdc: int | None = None,
        usar_ping: bool = True,
        direccion: str = "rtl",
    ) -> None:
        self.entropia = EntropiaFisica(usar_ping=usar_ping)
        self.mascara = MascaraLila()
        self.acciones = AccionesLila(n_mdc=n_mdc, direccion=direccion)
        self.usar_ping = usar_ping
        self.direccion = direccion
        self.ciclos_completados = 0
        self.primos_encontrados: list[int] = []
        self.candidatos_vistos: set[int] = set()
        self.k_vistos: set[int] = set()
        self.k_historial: list[int] = []
        self.decisiones: list[int] = []

    def ciclo(self, verbose: bool = True) -> dict:
        from nucleo.acciones import EstadoVolatil

        estado = EstadoVolatil()
        muestra = self.entropia.muestrear()
        decision = self.mascara.decidir(muestra)
        # Mezcla entropía del ciclo; RTL = escritura derecha→izquierda en el espacio k
        k_raw = (
            muestra.nanos
            ^ (muestra.basura << 7)
            ^ (muestra.bit_pila << 15)
            ^ (muestra.inercia_a << 22)
            ^ self.ciclos_completados
        ) % 10_000
        k_semilla = (9_999 - k_raw) if self.direccion == "rtl" else k_raw
        self.k_vistos.add(k_semilla)
        self.k_historial.append(k_semilla)
        self.decisiones.append(decision.decision)

        resultado = self.acciones.ejecutar(decision, k_semilla)
        if resultado.candidato is not None:
            self.candidatos_vistos.add(resultado.candidato)
        estado.ultimo_resultado = resultado
        estado.traza.append(resultado.mensaje)

        if resultado.es_primo and resultado.candidato is not None:
            self.primos_encontrados.append(resultado.candidato)

        if verbose:
            ping_txt = f"{muestra.ping_us} us" if muestra.ping_us else "sin red"
            print(f"Inercia(a): {decision.inercia_a} | Red(x): {decision.red_x} | ping: {ping_txt}")
            print(f"Medida: {decision.medida} | Curva: {decision.curva:.4f}")
            print(f"Eventos (0,5,9) -> r0={decision.r0} r5={decision.r5} r9={decision.r9}")
            print(f"Select Case mutante: {decision.decision} -> {ACCION_NOMBRES[decision.decision]}")
            print(f">>> {resultado.mensaje}")
            print(">>> CRIBA DESMEMORIADA: rastro del codigo de operar borrado de RAM.")

        estado.desmemoriar()
        self.ciclos_completados += 1

        return {
            "decision": decision.decision,
            "k": k_semilla,
            "candidato": resultado.candidato,
            "primo": resultado.es_primo,
            "m_mdc": resultado.m_mdc,
            "d_mdc": resultado.d_mdc,
        }

    def respirar(self, n: int = 10, verbose: bool = True) -> None:
        if verbose:
            print("=== ALEATOROVIX v1.0 — ORGANISMO DESMEMORIADO ===\n")
        for i in range(n):
            if verbose:
                print(f"--- Ciclo {i + 1}/{n} ---")
            self.ciclo(verbose=verbose)
            if verbose:
                print("----------------------------")
        if verbose and self.primos_encontrados:
            print(f"\nPrimos sintonizados en esta sesión: {sorted(set(self.primos_encontrados))}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Aleatorovix v1.0 — organismo Lila")
    parser.add_argument("-n", "--ciclos", type=int, default=10, help="Ciclos a ejecutar")
    parser.add_argument("--n-mdc", type=int, default=None, help="N para perfil d(m) en resonancia")
    parser.add_argument("--sin-ping", action="store_true", help="Solo entropía local (sin red)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Sin salida detallada")
    parser.add_argument(
        "--direccion",
        choices=("rtl", "ltr"),
        default="rtl",
        help="rtl=derecha→izquierda (pinza desde arriba); ltr=izquierda→derecha",
    )
    args = parser.parse_args()

    org = AleatorovixOrganismo(
        n_mdc=args.n_mdc,
        usar_ping=not args.sin_ping,
        direccion=args.direccion,
    )
    org.respirar(n=args.ciclos, verbose=not args.quiet)


if __name__ == "__main__":
    main()