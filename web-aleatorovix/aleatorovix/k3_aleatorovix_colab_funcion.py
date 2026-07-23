#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K3 + Aleatorovix — FUNCIÓN COLAB (100% headless) ❤️🌹
====================================================
Copia este archivo a Colab o pégalo en una celda.

Uso en Colab:
  !python k3_aleatorovix_colab_funcion.py
  # o en celdas:
  from k3_aleatorovix_colab_funcion import demo_todo, mural_industrial, ciclo, ventanas
  demo_todo()

NO usa tkinter. pandas es opcional.
GUI local: usa k3_aleatorovix_gui.py en tu PC (Windows).
"""

from __future__ import annotations

import math
import random
import time
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from typing import Dict, List, Optional, Tuple

getcontext().prec = 100

# ---------------------------------------------------------------------------
# Constantes (¡con salto de línea correcto!)
# ---------------------------------------------------------------------------
K_HEART = "❤️"  # bit 1 / factor_msl / K
L_ROSE = "🌹"  # bit 0 / factor_lsl / L
PAIR = f"{K_HEART}{L_ROSE}"
SEP = "·"


# ===========================================================================
# Notación floral
# ===========================================================================

def bits_a_flores(bits: str) -> str:
    return "".join(K_HEART if b == "1" else L_ROSE for b in bits if b in "01")


def digito_a_flores(d: int) -> str:
    return bits_a_flores(format(abs(int(d)) % 10, "04b"))


def numero_floral(n: int) -> str:
    if n < 0:
        return "−" + numero_floral(-n)
    if n == 0:
        return digito_a_flores(0)
    return SEP.join(digito_a_flores(int(c)) for c in str(int(n)))


def entero_a_flores(n: int, ancho: Optional[int] = None) -> str:
    if n < 0:
        return "−" + entero_a_flores(-n, ancho)
    bits = bin(int(n))[2:] if n else "0"
    if ancho:
        bits = bits.zfill(ancho)
    return bits_a_flores(bits)


def render_acordeon(bits: str, fase: int = 0) -> List[str]:
    flores = [K_HEART if b == "1" else L_ROSE for b in bits if b in "01"]
    if not flores:
        return [PAIR]
    gap = "  " if fase % 2 == 0 else " "
    return [
        "".join(flores),
        " ".join(flores),
        gap.join(flores),
        " ".join(flores),
        "".join(flores),
    ]


def formatear_valor_floral(etiqueta: str, n: Optional[int]) -> str:
    if n is None:
        return f"  {etiqueta}: {L_ROSE * 4}  (∅)"
    return f"  {etiqueta}: {numero_floral(n)}  ⟦bin {entero_a_flores(n)}⟧"


def contar_flores(texto: str) -> Tuple[int, int]:
    return texto.count(K_HEART), texto.count(L_ROSE)


def _mostrar_tabla(rows: List[dict]) -> None:
    try:
        import pandas as pd

        df = pd.DataFrame(rows)
        try:
            from IPython.display import display  # type: ignore

            display(df)
        except Exception:
            print(df.to_string(index=False))
    except ImportError:
        for r in rows:
            print(r)


# ===========================================================================
# Núcleo Lila + Aleatorovix
# ===========================================================================

@dataclass
class MuestraEntropia:
    nanos: int
    bit_pila: int
    basura: int
    inercia_a: int
    ping_us: Optional[int] = None


@dataclass
class DecisionLila:
    inercia_a: int
    red_x: int
    medida: int
    curva: float
    r0: int
    r5: int
    r9: int
    decision: int


@dataclass
class ResultAccion:
    candidato: Optional[int]
    es_primo: bool
    mensaje: str
    m_mdc: Optional[int]
    d_mdc: Optional[int]


class EntropiaFisica:
    def __init__(self, usar_ping: bool = False):
        self.usar_ping = usar_ping
        self.contador_pila = 0
        self.basura_heap = [random.randint(0, 255) for _ in range(128)]

    def muestrear(self) -> MuestraEntropia:
        nanos = int(time.time_ns()) % 10000
        self.contador_pila += 1
        bit_pila = (self.contador_pila >> 6) & 0x01
        basura = self.basura_heap[nanos % len(self.basura_heap)]
        inercia_a = nanos % 1000
        ping_us = random.randint(100, 500) if self.usar_ping else None
        return MuestraEntropia(nanos, bit_pila, basura, inercia_a, ping_us)


class MascaraLila:
    def mascara_sagrada(self, x: int, a: int) -> float:
        if x <= 0:
            x = 1
        dx, da = float(x), float(max(a, 0))
        return (-math.pow(10.0, 1.0 / dx) + math.pow(10.0, 1.0 / (dx + da)) - 1.0) * dx

    def pulso_gauss(self, v: int, p: int) -> int:
        return 1 if abs(v - p) <= 1 else 0

    def interprete_mutante(self, b_ext: int, b_pila: int, ruido: int) -> int:
        par = (b_ext << 1) | b_pila
        rot = ruido % 3
        if par in (0b11, 0b00):
            base = 0
        elif par == 0b01:
            base = 1
        else:
            base = 2
        return (base + rot) % 3

    def decidir(self, muestra: MuestraEntropia) -> DecisionLila:
        x = muestra.ping_us if muestra.ping_us is not None else (muestra.nanos % 500 + 100)
        a = max(muestra.inercia_a, 1)
        curva = self.mascara_sagrada(x, a)
        medida = abs(int(curva)) % 10
        ruido = muestra.basura & 0xFF
        r0 = self.interprete_mutante(self.pulso_gauss(medida, 0), muestra.bit_pila, ruido)
        r5 = self.interprete_mutante(self.pulso_gauss(medida, 5), muestra.bit_pila, ruido ^ 0x55)
        r9 = self.interprete_mutante(self.pulso_gauss(medida, 9), muestra.bit_pila, ruido ^ 0xAA)
        return DecisionLila(a, x, medida, curva, r0, r5, r9, (r0 + r5 + r9) % 4)


class AccionesLila:
    NOMBRES = {
        0: f"[OLVIDO] Criba Desmemoriada {PAIR}",
        1: f"[SALTA] 6k+1 {K_HEART}",
        2: f"[BAILA] 6k-1 {L_ROSE}",
        3: f"[ZYPYZAPE] Resonancia {PAIR}",
    }

    @staticmethod
    def _es_primo(n: int) -> bool:
        if n < 2:
            return False
        if n % 2 == 0:
            return n == 2
        r = int(math.isqrt(n))
        for i in range(3, r + 1, 2):
            if n % i == 0:
                return False
        return True

    def ejecutar(self, decision: DecisionLila, k: int) -> ResultAccion:
        d = decision.decision
        if d == 0:
            return ResultAccion(None, False, self.NOMBRES[0], None, None)
        if d == 1:
            c = 6 * k + 1
            return ResultAccion(c, self._es_primo(c), f"{self.NOMBRES[1]}: {c}", None, None)
        if d == 2:
            c = 6 * k - 1
            return ResultAccion(c, self._es_primo(c), f"{self.NOMBRES[2]}: {c}", None, None)
        resto = k % 121
        return ResultAccion(None, False, f"{self.NOMBRES[3]} resto {resto}", resto, 10403)


class AleatorovixAcordeon:
    def __init__(self, semilla: int = 33):
        self.estado = Decimal(str(semilla or 33))
        self.factor_msl = Decimal("42.123456789")  # ❤️
        self.factor_lsl = Decimal("0.000123456")  # 🌹
        self.fase = 0
        self.corazones = 0
        self.rosas = 0

    def bit(self) -> str:
        term = Decimal(str(math.cos(float(self.estado * self.factor_lsl))))
        fase = (self.estado * self.factor_msl) + term
        self.estado = fase - fase.to_integral_value(rounding="ROUND_FLOOR")
        b = "1" if self.estado > Decimal("0.5") else "0"
        if b == "1":
            self.corazones += 1
        else:
            self.rosas += 1
        return b

    def flor(self) -> str:
        self.fase += 1
        return K_HEART if self.bit() == "1" else L_ROSE

    def chorro(self, n: int) -> str:
        return "".join(self.bit() for _ in range(n))

    def chorro_flores(self, n: int) -> str:
        return "".join(self.flor() for _ in range(n))

    def numero_posible(self, bits: int = 8) -> Tuple[int, str]:
        b = self.chorro(bits)
        return int(b, 2), bits_a_flores(b)

    def ventanas_moviles(
        self,
        total_bits: int = 512,
        num_ventanas: int = 12,
        min_ancho: int = 8,
        max_ancho: int = 64,
    ) -> List[dict]:
        total_bits = max(1, total_bits)
        min_ancho = max(1, min(min_ancho, total_bits))
        max_ancho = max(min_ancho, min(max_ancho, total_bits))
        base = self.chorro(total_bits)
        out = []
        for i in range(num_ventanas):
            ancho = random.randint(min_ancho, max_ancho)
            start = random.randint(0, max(0, total_bits - ancho))
            seg = base[start : start + ancho]
            if not seg:
                continue
            valor = int(seg, 2)
            out.append(
                {
                    "ventana": i + 1,
                    "start": start,
                    "end": start + ancho,
                    "ancho": ancho,
                    "flores": bits_a_flores(seg),
                    "floral": numero_floral(valor % 10**12),
                    "valor": valor,
                    "corazones": seg.count("1"),
                    "rosas": seg.count("0"),
                }
            )
        return out


class K3Industrial:
    """Acordeón tipo k3_core.c (base=33, rel=1) → bytes → flores."""

    def __init__(self, v: int = 33, L: int = 1):
        self.v = v & 0xFFFFFFFFFFFFFFFF
        self.L = L & 0xFFFFFFFFFFFFFFFF
        self.corazones = 0
        self.rosas = 0

    def step_byte(self) -> int:
        self.L = (self.L + (self.v + 1)) & 0xFFFFFFFFFFFFFFFF
        self.v = (self.v + 2) & 0xFFFFFFFFFFFFFFFF
        if (5 * self.L) <= (2 * self.v + 1):
            self.L = (self.L + self.v * 2) & 0xFFFFFFFFFFFFFFFF
        stream = ((self.L ^ self.v) * 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        return stream & 0xFF

    def chorro_bits(self, n: int) -> str:
        bits = []
        while len(bits) < n:
            for ch in f"{self.step_byte():08b}":
                bits.append(ch)
                if ch == "1":
                    self.corazones += 1
                else:
                    self.rosas += 1
                if len(bits) >= n:
                    break
        return "".join(bits)

    def chorro_flores(self, n: int) -> str:
        return bits_a_flores(self.chorro_bits(n))


class AleatorovixOrganismo:
    def __init__(self, usar_ping: bool = False):
        self.entropia = EntropiaFisica(usar_ping=usar_ping)
        self.mascara = MascaraLila()
        self.acciones = AccionesLila()
        self.acordeon = AleatorovixAcordeon(33)
        self.ciclos = 0
        self.primos: List[int] = []
        self.hist: List[int] = []

    @staticmethod
    def k3_transform(n: int) -> int:
        if n == 12:
            return n + 12
        if n == 1:
            return n + 10
        if n == 11:
            return n + 2
        return n

    def ciclo_completo(self) -> dict:
        m = self.entropia.muestrear()
        d = self.mascara.decidir(m)
        k_raw = (m.nanos ^ (m.basura << 7) ^ (m.bit_pila << 15)) % 10000
        k = max(1, 9999 - k_raw)
        res = self.acciones.ejecutar(d, k)
        if res.candidato is not None:
            self.primos.append(res.candidato)
        self.hist.append(d.decision)
        self.ciclos += 1
        bits = self.acordeon.chorro(16)
        num, flores_p = self.acordeon.numero_posible(8)
        return {
            "decision": d.decision,
            "k_semilla": k,
            "candidato": res.candidato,
            "primo": res.es_primo,
            "mensaje": res.mensaje,
            "inercia_a": d.inercia_a,
            "red_x": d.red_x,
            "medida": d.medida,
            "curva": d.curva,
            "r0": d.r0,
            "r5": d.r5,
            "r9": d.r9,
            "flores_acordeon": bits_a_flores(bits),
            "flores_posible": flores_p,
            "num_posible": num,
            "k_floral": numero_floral(k),
            "decision_floral": entero_a_flores(d.decision, 2),
            "candidato_floral": numero_floral(res.candidato) if res.candidato else None,
            "acordeon_vista": "\n".join(render_acordeon(bits, self.ciclos)),
        }


# ===========================================================================
# FUNCIONES PÚBLICAS (para celdas de Colab)
# ===========================================================================

def ciclo(n: int = 1, verbose: bool = True) -> List[dict]:
    """Ejecuta n ciclos Aleatorovix (sin GUI)."""
    org = AleatorovixOrganismo(usar_ping=False)
    out = []
    for _ in range(n):
        r = org.ciclo_completo()
        out.append(r)
        if verbose:
            print(f"Ciclo {numero_floral(org.ciclos)} {PAIR}")
            print(formatear_valor_floral("inercia", r["inercia_a"]))
            print(formatear_valor_floral("red_x", r["red_x"]))
            print(f"  decisión: {r['decision_floral']}  {r['mensaje']}")
            if r["candidato_floral"]:
                tag = f" {L_ROSE}PRIMO{L_ROSE}" if r["primo"] else ""
                print(f"  candidato: {r['candidato_floral']}{tag}")
            print(f"  chorro:  {r['flores_acordeon']}")
            print(f"  posible: {r['flores_posible']}")
            print(r["acordeon_vista"])
            print("-" * 40)
    return out


def ventanas(
    total_bits: int = 512,
    num: int = 12,
    min_ancho: int = 8,
    max_ancho: int = 64,
) -> List[dict]:
    """Ventanas móviles → muchos números en ❤️🌹."""
    print(f"=== VENTANAS MÓVILES {PAIR} ===")
    ac = AleatorovixAcordeon(33)
    data = ac.ventanas_moviles(total_bits, num, min_ancho, max_ancho)
    # tabla reducida
    tabla = [
        {
            "v": w["ventana"],
            "start": w["start"],
            "ancho": w["ancho"],
            f"{K_HEART}": w["corazones"],
            f"{L_ROSE}": w["rosas"],
            "floral": w["floral"][:40] + ("…" if len(w["floral"]) > 40 else ""),
            "flores_preview": w["flores"][:48] + ("…" if len(w["flores"]) > 48 else ""),
        }
        for w in data
    ]
    _mostrar_tabla(tabla)
    for w in data[:5]:
        print(f"\nv{w['ventana']} [{w['start']}:{w['end']}] ancho={w['ancho']}")
        print(f"  {w['flores'][:96]}{'…' if len(w['flores']) > 96 else ''}")
        print(f"  floral: {w['floral']}")
        # pulso acordeón compacto (3 fases)
        bits = "".join(
            "1" if ch == "1" else "0"
            for ch in bin(w["valor"])[2:].zfill(min(w["ancho"], 32))
        )
        for line in render_acordeon(bits, w["ventana"])[:3]:
            print(f"    {line}")
    print(
        f"\nTotal en ventanas: {K_HEART}{sum(w['corazones'] for w in data)} "
        f"{L_ROSE}{sum(w['rosas'] for w in data)}"
    )
    return data


def mural_industrial(
    lineas: int = 24,
    ancho: int = 64,
    dual: bool = True,
) -> Dict:
    """Mural denso de ❤️🌹 (modo industrial Colab)."""
    print(f"{PAIR * 12}")
    print(f"  PLANTA FLORAL COLAB 33×1  {PAIR}")
    print(f"{PAIR * 12}\n")
    ale = AleatorovixAcordeon(33)
    k3 = K3Industrial(33, 1)
    filas = []
    half = lineas // 2 if dual else lineas
    for _ in range(half):
        filas.append(ale.chorro_flores(ancho))
    if dual:
        for _ in range(lineas - half):
            filas.append(k3.chorro_flores(ancho))
    for row in filas:
        print(row)
    texto = "\n".join(filas)
    h, r = contar_flores(texto)
    info = {
        "lineas": len(filas),
        "ancho": ancho,
        "corazones": h,
        "rosas": r,
        "total": h + r,
        "ratio_heart": h / (h + r) if (h + r) else 0,
    }
    print(
        f"\n{K_HEART}={h}  {L_ROSE}={r}  total={h + r}  "
        f"ratio❤️={info['ratio_heart']:.3f}"
    )
    print(f"floral total: {numero_floral(h + r)}")
    return info


def lluvia(n_flores: int = 5000, por_linea: int = 80) -> Dict:
    """Lluvia de miles de flores."""
    print(f"=== LLUVIA {n_flores} flores {PAIR} ===\n")
    ale = AleatorovixAcordeon(33)
    k3 = K3Industrial(33, 1)
    filas = []
    left = n_flores
    use_ale = True
    while left > 0:
        take = min(por_linea, left)
        filas.append(ale.chorro_flores(take) if use_ale else k3.chorro_flores(take))
        use_ale = not use_ale
        left -= take
    for row in filas[:20]:
        print(row)
    if len(filas) > 20:
        print(f"... +{len(filas) - 20} filas")
    texto = "".join(filas)
    h, r = contar_flores(texto)
    print(f"\n{K_HEART}={h}  {L_ROSE}={r}  total={h + r}")
    return {"corazones": h, "rosas": r, "total": h + r, "filas": len(filas)}


def secuencia_k3(length: int = 10) -> None:
    print(f"=== SECUENCIA K3 {PAIR} (floral) ===")
    base = [7, 1, 12, 4, 11, 3, 0, 8, 1, 12]
    for i in range(length):
        n = base[i] if i < len(base) else random.randint(0, 12)
        t = AleatorovixOrganismo.k3_transform(n)
        print(
            f"N={entero_a_flores(i, 4)}: {numero_floral(n)} → {numero_floral(t)} → "
            f"{bits_a_flores(bin(t)[2:])}"
        )


def alfabeto() -> None:
    print(f"Alfabeto nibble 0–9 {PAIR}:")
    for d in range(10):
        print(f"  {d} → {digito_a_flores(d)}")


def demo_todo() -> None:
    """Una sola función: todo lo útil en Colab."""
    print(f"\n{'=' * 60}")
    print(f"  K3 + ALEATOROVIX — DEMO COLAB {PAIR}")
    print(f"{'=' * 60}\n")
    alfabeto()
    print()
    secuencia_k3(8)
    print()
    ciclo(3)
    print()
    ventanas(256, 8, 8, 32)
    print()
    mural_industrial(16, 48, dual=True)
    print()
    lluvia(3000, 64)
    print(f"\n{PAIR} fin demo — en local con GUI: python k3_aleatorovix_gui.py")
    print("En Colab no hay pantalla gráfica; este script es el sustituto.\n")


# ===========================================================================
if __name__ == "__main__":
    import sys

    arg = (sys.argv[1] if len(sys.argv) > 1 else "todo").lower()
    if arg in ("ciclo", "sim"):
        ciclo(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
    elif arg in ("ventanas", "windows"):
        ventanas()
    elif arg in ("mural", "industrial"):
        lineas = int(sys.argv[2]) if len(sys.argv) > 2 else 24
        mural_industrial(lineas, 64)
    elif arg in ("lluvia", "rain"):
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
        lluvia(n)
    elif arg in ("k3",):
        secuencia_k3()
    else:
        demo_todo()
