#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K3 + Aleatorovix GUI v2.0 ❤️🌹
Víctor Manzanares Alberola — uso civil / educativo

Basado en:
  - organismo_lila_v99.c (repo/just run/aleatorovix/)
  - gemini-code-1784158392232.py (factor_msl / factor_lsl)
  - Reglas K3: 12→+12, 1→+10, 11→+2

Símbolos:
  ❤️ = K = factor_msl (lazo significativo mayor)
  🌹 = L = factor_lsl (lazo significativo menor)

Ejecutar:
  python k3_aleatorovix_gui.py
"""

from __future__ import annotations

import math
import random
import time
import tkinter as tk
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from typing import List, Optional, Tuple

getcontext().prec = 80

# ---------------------------------------------------------------------------
# Constantes simbólicas
# ---------------------------------------------------------------------------
K_HEART = "❤️"   # factor_msl / bit 1
L_ROSE = "🌹"    # factor_lsl / bit 0
PAIR = f"{K_HEART}{L_ROSE}"
SEP_DIGITO = "·"  # separador entre dígitos en notación floral


# ===========================================================================
# Notación acordeón ❤️🌹 (números NO estándar)
# ===========================================================================
# bit 1 → ❤️ (K / msl)   bit 0 → 🌹 (L / lsl)
# Un entero se lee como serie de flores, no como dígitos arábigos.


def bits_a_flores(bits: str) -> str:
    """Chorrito binario → serie ❤️/🌹."""
    return "".join(K_HEART if b == "1" else L_ROSE for b in bits if b in "01")


def flores_a_bits(flores: str) -> str:
    """Serie ❤️/🌹 → bits (ignora separadores)."""
    out = []
    i = 0
    s = flores
    while i < len(s):
        if s.startswith(K_HEART, i):
            out.append("1")
            i += len(K_HEART)
        elif s.startswith(L_ROSE, i):
            out.append("0")
            i += len(L_ROSE)
        else:
            i += 1  # ·, espacios, etc.
    return "".join(out)


def entero_a_flores(n: int, ancho: Optional[int] = None) -> str:
    """Entero → binario en flores (lectura no estándar del número)."""
    if n < 0:
        return "−" + entero_a_flores(-n, ancho)
    bits = bin(int(n))[2:] if n != 0 else "0"
    if ancho is not None and ancho > 0:
        bits = bits.zfill(ancho)
    return bits_a_flores(bits)


def digito_a_flores(d: int) -> str:
    """Un dígito 0–9 como nibble de 4 flores (alfabeto floral fijo)."""
    d = abs(int(d)) % 10
    return bits_a_flores(format(d, "04b"))


def numero_floral(n: int) -> str:
    """
    Número en notación floral por dígitos (no arábigo).
    Ej: 12 → 🌹🌹❤️🌹·🌹🌹❤️❤️  (1=0001, 2=0010 en nibble)
    """
    if n < 0:
        return "−" + numero_floral(-n)
    if n == 0:
        return digito_a_flores(0)
    digitos = [int(c) for c in str(int(n))]
    return SEP_DIGITO.join(digito_a_flores(d) for d in digitos)


def flores_a_entero(flores: str) -> int:
    """Recupera el entero desde binario floral (sin nibbles)."""
    bits = flores_a_bits(flores.replace(SEP_DIGITO, ""))
    return int(bits, 2) if bits else 0


def render_acordeon(bits: str, fase: int = 0) -> List[str]:
    """
    Visual modo acordeón: estira / encoge la serie ❤️🌹.
    fase par = expandido (espacios), impar = comprimido.
    """
    flores = [K_HEART if b == "1" else L_ROSE for b in bits if b in "01"]
    if not flores:
        return [PAIR]
    lineas = []
    # contracción (cerrado)
    lineas.append("".join(flores))
    # semi-abierto
    lineas.append(" ".join(flores))
    # abierto (estirado msl)
    gap = "  " if fase % 2 == 0 else " "
    lineas.append(gap.join(flores))
    # semi otra vez
    lineas.append(" ".join(flores))
    # cerrado
    lineas.append("".join(flores))
    return lineas


def formatear_valor_floral(etiqueta: str, n: Optional[int], bits_vivo: Optional[str] = None) -> str:
    """Una línea: etiqueta + valor solo en flores (+ binario floral compacto)."""
    if n is None:
        return f"  {etiqueta}: {L_ROSE}{L_ROSE}{L_ROSE}{L_ROSE}  (∅)"
    floral_dig = numero_floral(n)
    floral_bin = entero_a_flores(n)
    extra = f"  | chorro:{bits_a_flores(bits_vivo)}" if bits_vivo else ""
    return f"  {etiqueta}: {floral_dig}  ⟦bin {floral_bin}⟧{extra}"


# ===========================================================================
# Núcleo matemático (motor real, sin bloquear la GUI)
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
    """Muestreo de entropía (nanos / pila / basura / ping simulado)."""

    def __init__(self, usar_ping: bool = True):
        self.usar_ping = usar_ping
        self.contador_pila = 0
        self.basura_heap = [random.randint(0, 255) for _ in range(128)]

    def muestrear(self) -> MuestraEntropia:
        nanos = int(time.time_ns()) % 10000
        self.contador_pila += 1
        bit_pila = (self.contador_pila >> 6) & 0x01
        basura = self.basura_heap[nanos % len(self.basura_heap)]
        inercia_a = nanos % 1000  # como a_propio en C (get_nanos % 1000)
        ping_us = random.randint(100, 500) if self.usar_ping else None
        return MuestraEntropia(nanos, bit_pila, basura, inercia_a, ping_us)


class MascaraLila:
    """Máscara de curvatura Lila + campanas de Gauss + intérprete mutante."""

    def mascara_sagrada(self, x: int, a: int) -> float:
        # C: ( -10^{1/x} + 10^{1/(x+a)} - 1 ) * x
        if x <= 0:
            x = 1
        dx = float(x)
        da = float(max(a, 0))
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
        decision = (r0 + r5 + r9) % 4

        return DecisionLila(
            inercia_a=a,
            red_x=x,
            medida=medida,
            curva=curva,
            r0=r0,
            r5=r5,
            r9=r9,
            decision=decision,
        )


class AccionesLila:
    ACCION_NOMBRES = {
        0: f"[OLVIDO] Criba Desmemoriada {PAIR}",
        1: f"[SALTA] 6k+1 {K_HEART}",
        2: f"[BAILA] 6k-1 {L_ROSE}",
        3: f"[ZYPYZAPE] Resonancia {PAIR}",
    }

    def __init__(self, n_mdc: Optional[int] = 10403):
        self.n_mdc = n_mdc

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

    def ejecutar(self, decision: DecisionLila, k_semilla: int) -> ResultAccion:
        d = decision.decision
        if d == 0:
            return ResultAccion(None, False, self.ACCION_NOMBRES[0], None, None)
        if d == 1:
            c = 6 * k_semilla + 1
            return ResultAccion(c, self._es_primo(c), f"{self.ACCION_NOMBRES[1]}: {c}", None, None)
        if d == 2:
            c = 6 * k_semilla - 1
            return ResultAccion(c, self._es_primo(c), f"{self.ACCION_NOMBRES[2]}: {c}", None, None)
        resto = k_semilla % 121
        return ResultAccion(
            None,
            False,
            f"{self.ACCION_NOMBRES[3]} con resto {resto}",
            resto,
            self.n_mdc,
        )


class AleatorovixAcordeon:
    """
    Chorro caótico K3 (gemini): ❤️ factor_msl / 🌹 factor_lsl.
    Cada bit es un pliegue del acordeón: 1=❤️ (estirar), 0=🌹 (plegar).
    """

    def __init__(self, semilla: int = 123456789):
        self.estado = Decimal(str(semilla if semilla != 0 else 123456789))
        self.factor_msl = Decimal("42.123456789")   # ❤️ K
        self.factor_lsl = Decimal("0.000123456")    # 🌹 L
        self.fase_acordeon = 0
        self.historial_flores: List[str] = []

    def bit(self) -> str:
        # Decimal no expone cos(); usamos float solo para el término armónico
        term_cos = Decimal(str(math.cos(float(self.estado * self.factor_lsl))))
        fase = (self.estado * self.factor_msl) + term_cos
        self.estado = fase - fase.to_integral_value(rounding="ROUND_FLOOR")
        return "1" if self.estado > Decimal("0.5") else "0"

    def flor(self) -> str:
        """Un pliegue: ❤️ o 🌹."""
        b = self.bit()
        f = K_HEART if b == "1" else L_ROSE
        self.historial_flores.append(f)
        self.fase_acordeon += 1
        return f

    def chorro(self, n: int) -> str:
        return "".join(self.bit() for _ in range(n))

    def chorro_flores(self, n: int) -> str:
        """Serie de rosas y corazones generada por el acordeón vivo."""
        return "".join(self.flor() for _ in range(n))

    def numero_posible(self, bits: int = 8) -> Tuple[int, str]:
        """
        Extrae un número 'posible' del chorro (no se escribe en arábigo).
        Devuelve (valor_interno, serie_floral).
        """
        b = self.chorro(bits)
        return int(b, 2), bits_a_flores(b)

    def mostrar_serie(self, n_pliegues: int = 16) -> str:
        """Texto multi-línea del acordeón abriéndose y cerrándose."""
        bits = self.chorro(n_pliegues)
        lineas = render_acordeon(bits, self.fase_acordeon)
        valor = int(bits, 2)
        cab = (
            f"acordeón {PAIR}  pliegues={n_pliegues}  "
            f"lectura_floral={bits_a_flores(bits)}"
        )
        # valor solo en notación floral (sin dígitos arábigos)
        pie = f"valor⟦floral⟧: {numero_floral(valor)}  |  bin⟦{entero_a_flores(valor)}⟧"
        self.fase_acordeon += 1
        return "\n".join([cab, *lineas, pie])

    def chorro_ventanas_moviles(
        self,
        total_bits: int = 128,
        num_ventanas: int = 5,
        min_ancho: int = 8,
        max_ancho: int = 32,
    ) -> List[dict]:
        """
        Varios números a partir de un chorro grande: ventanas móviles
        con posición y ancho variables (acordeón que se mueve).
        """
        total_bits = max(1, int(total_bits))
        min_ancho = max(1, min(min_ancho, total_bits))
        max_ancho = max(min_ancho, min(max_ancho, total_bits))
        bits_base = self.chorro(total_bits)
        resultados: List[dict] = []

        for i in range(num_ventanas):
            ancho = random.randint(min_ancho, max_ancho)
            max_start = total_bits - ancho
            start_idx = random.randint(0, max_start) if max_start > 0 else 0
            end_idx = start_idx + ancho
            segmento = bits_base[start_idx:end_idx]
            if not segmento:
                continue
            valor = int(segmento, 2)
            resultados.append(
                {
                    "ventana": i + 1,
                    "start_idx": start_idx,
                    "end_idx": end_idx,
                    "ancho": ancho,
                    "bits": segmento,
                    "flores": bits_a_flores(segmento),
                    "floral_digitos": numero_floral(valor),
                    "valor": valor,
                    "vista": "\n".join(render_acordeon(segmento, i)),
                }
            )
        self.fase_acordeon += 1
        return resultados


class AleatorovixOrganismo:
    """Organismo Lila + reglas K3 + acordeón ❤️🌹."""

    def __init__(self, usar_ping: bool = True):
        self.entropia = EntropiaFisica(usar_ping=usar_ping)
        self.mascara = MascaraLila()
        self.acciones = AccionesLila()
        self.acordeon = AleatorovixAcordeon()
        self.ciclos_completados = 0
        self.primos_encontrados: List[int] = []
        self.historial_decisiones: List[int] = []

    @staticmethod
    def k3_transform(n: int) -> int:
        """Reglas K3 del usuario: 12→+12, 1→+10, 11→+2."""
        if n == 12:
            return n + 12
        if n == 1:
            return n + 10
        if n == 11:
            return n + 2
        return n

    def generar_secuencia_k3(self, length: int = 10) -> List[Tuple[int, int, str]]:
        seq = []
        # mezcla determinista + aleatoria; incluye casos especiales
        base = [7, 1, 12, 4, 11, 3, 0, 8, 1, 12]
        for i in range(length):
            n = base[i] if i < len(base) else random.randint(0, 12)
            t = self.k3_transform(n)
            seq.append((n, t, bin(t)[2:]))
        return seq

    def ciclo_completo(self) -> dict:
        muestra = self.entropia.muestrear()
        decision = self.mascara.decidir(muestra)
        k_raw = (muestra.nanos ^ (muestra.basura << 7) ^ (muestra.bit_pila << 15)) % 10000
        k_semilla = max(1, 9999 - k_raw)
        resultado = self.acciones.ejecutar(decision, k_semilla)

        if resultado.candidato is not None:
            self.primos_encontrados.append(resultado.candidato)
        self.historial_decisiones.append(decision.decision)
        self.ciclos_completados += 1

        bits = self.acordeon.chorro(16)
        num_posible, flores_posible = self.acordeon.numero_posible(8)
        flores_vivo = bits_a_flores(bits)

        return {
            "decision": decision.decision,
            "k_semilla": k_semilla,
            "candidato": resultado.candidato,
            "primo": resultado.es_primo,
            "mensaje": resultado.mensaje,
            "inercia_a": decision.inercia_a,
            "red_x": decision.red_x,
            "medida": decision.medida,
            "curva": decision.curva,
            "r0": decision.r0,
            "r5": decision.r5,
            "r9": decision.r9,
            "ping_us": muestra.ping_us,
            "bits_acordeon": bits,
            "flores_acordeon": flores_vivo,
            "num_posible": num_posible,
            "flores_posible": flores_posible,
            "k_floral": numero_floral(k_semilla),
            "k_bin_floral": entero_a_flores(k_semilla),
            "decision_floral": entero_a_flores(decision.decision, 2),
            "medida_floral": entero_a_flores(decision.medida, 4),
            "candidato_floral": (
                numero_floral(resultado.candidato) if resultado.candidato is not None else None
            ),
            "acordeon_vista": "\n".join(render_acordeon(bits, self.ciclos_completados)),
            "msl": str(self.acordeon.factor_msl),
            "lsl": str(self.acordeon.factor_lsl),
        }

    def serie_acordeon_flores(self, pliegues: int = 24, cuantos: int = 8) -> List[dict]:
        """
        Genera una serie de números posibles solo como ❤️🌹 (modo acordeón).
        Cada entrada es un pliegue que 'muestra' un número sin dígitos estándar.
        """
        serie = []
        for i in range(cuantos):
            valor, flores = self.acordeon.numero_posible(pliegues // max(cuantos // 2, 1) or 8)
            # ancho variable (estirar/plegar): 4..12 bits
            ancho = 4 + (i % 9)
            valor2, flores2 = self.acordeon.numero_posible(ancho)
            serie.append(
                {
                    "pliegue": i + 1,
                    "flores": flores2,
                    "floral_digitos": numero_floral(valor2),
                    "ancho": ancho,
                    "vista": "\n".join(render_acordeon(flores_a_bits(flores2), i)),
                    # valor interno solo para cálculo; la lectura pública es floral
                    "_valor": valor2,
                }
            )
        return serie

    def reiniciar(self) -> None:
        self.ciclos_completados = 0
        self.primos_encontrados.clear()
        self.historial_decisiones.clear()
        self.acordeon = AleatorovixAcordeon(random.randint(1, 10**9))


# ===========================================================================
# GUI
# ===========================================================================

class K3AleatorovixGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title(f"{K_HEART} K3 + Aleatorovix GUI v2.0 {L_ROSE}")
        self.root.geometry("960x720")
        self.root.minsize(800, 600)

        self.organismo = AleatorovixOrganismo(usar_ping=True)
        self.running = False
        self.velocidad_ms = 1000
        self.detallado = tk.BooleanVar(value=True)
        self.usar_ping = tk.BooleanVar(value=True)

        self._build()

    def _build(self) -> None:
        style_font = ("Segoe UI Emoji", 10)
        title_font = ("Segoe UI Emoji", 16, "bold")

        main = tk.Frame(self.root, padx=10, pady=10)
        main.pack(fill=tk.BOTH, expand=True)

        tk.Label(
            main,
            text=f"{K_HEART} K3 + Aleatorovix {L_ROSE} — Panel de Control",
            font=title_font,
        ).pack(pady=(0, 4))

        tk.Label(
            main,
            text=f"K = {K_HEART} factor_msl  |  L = {L_ROSE} factor_lsl",
            font=style_font,
        ).pack(pady=(0, 8))

        # Notebook-like tabs with frames
        notebook = tk.Frame(main)
        notebook.pack(fill=tk.BOTH, expand=True)

        self.tab_buttons = tk.Frame(notebook)
        self.tab_buttons.pack(fill=tk.X)

        self.pages = {}
        self.page_frames = tk.Frame(notebook)
        self.page_frames.pack(fill=tk.BOTH, expand=True)

        for key, label in [
            ("control", f"🎛️ Control {PAIR}"),
            ("viz", f"📊 Visualización {PAIR}"),
            ("k3", f"🔺 K3 {PAIR}"),
            ("stats", f"📈 Estadísticas {PAIR}"),
            ("info", f"ℹ️ Información {PAIR}"),
        ]:
            btn = tk.Button(
                self.tab_buttons,
                text=label,
                command=lambda k=key: self._show_page(k),
                font=("Segoe UI Emoji", 9),
            )
            btn.pack(side=tk.LEFT, padx=2, pady=2)
            frame = tk.Frame(self.page_frames)
            self.pages[key] = frame

        self._build_control(self.pages["control"])
        self._build_viz(self.pages["viz"])
        self._build_k3(self.pages["k3"])
        self._build_stats(self.pages["stats"])
        self._build_info(self.pages["info"])
        self._show_page("control")

        self.status = tk.Label(
            main,
            text=f"Listo {PAIR}  |  Ciclos: 0 {K_HEART}  |  Primos: 0 {L_ROSE}",
            font=style_font,
            anchor="w",
        )
        self.status.pack(fill=tk.X, pady=(8, 0))

    def _show_page(self, key: str) -> None:
        for f in self.pages.values():
            f.pack_forget()
        self.pages[key].pack(fill=tk.BOTH, expand=True)

    def _build_control(self, parent: tk.Frame) -> None:
        box = tk.LabelFrame(parent, text=f"Acciones {PAIR}", padx=10, pady=10)
        box.pack(fill=tk.X, pady=5)

        row = tk.Frame(box)
        row.pack(fill=tk.X)

        tk.Button(row, text=f"🎲 Generar K3 {PAIR}", command=self.generate_k3, width=18).pack(
            side=tk.LEFT, padx=4, pady=4
        )
        tk.Button(
            row, text=f"🌊 Ciclo Aleatorovix {PAIR}", command=self.run_cycle, width=22
        ).pack(side=tk.LEFT, padx=4, pady=4)
        tk.Button(
            row, text=f"🪗 Acordeón {PAIR}", command=self.run_acordeon_serie, width=16
        ).pack(side=tk.LEFT, padx=4, pady=4)
        tk.Button(
            row, text=f"🪟 Ventanas {PAIR}", command=self.run_ventanas_moviles, width=16
        ).pack(side=tk.LEFT, padx=4, pady=4)
        tk.Button(
            row, text=f"🏭 Industrial {PAIR}", command=self.run_industrial, width=16
        ).pack(side=tk.LEFT, padx=4, pady=4)
        tk.Button(row, text=f"🔄 Auto Run {PAIR}", command=self.toggle_auto, width=14).pack(
            side=tk.LEFT, padx=4, pady=4
        )
        tk.Button(row, text=f"🧹 Reiniciar {PAIR}", command=self.reset_all, width=14).pack(
            side=tk.LEFT, padx=4, pady=4
        )

        opts = tk.Frame(box)
        opts.pack(fill=tk.X, pady=6)
        tk.Checkbutton(
            opts,
            text=f"Usar ping simulado {K_HEART}",
            variable=self.usar_ping,
            command=self._toggle_ping,
        ).pack(side=tk.LEFT, padx=6)
        tk.Checkbutton(
            opts, text=f"Salida detallada {L_ROSE}", variable=self.detallado
        ).pack(side=tk.LEFT, padx=6)

        speed = tk.Frame(box)
        speed.pack(fill=tk.X, pady=4)
        tk.Label(speed, text=f"⚡ Velocidad (ms entre ciclos) {PAIR}:").pack(side=tk.LEFT)
        self.speed_scale = tk.Scale(
            speed, from_=100, to=2000, orient=tk.HORIZONTAL, length=280, resolution=50
        )
        self.speed_scale.set(1000)
        self.speed_scale.pack(side=tk.LEFT, padx=8)

        help_txt = (
            f"Reglas K3: 12→+12, 1→+10, 11→+2\n"
            f"{K_HEART}=bit1/msl  |  {L_ROSE}=bit0/lsl  |  números en flores (no arábigos)\n"
            f"Acciones: OLVIDO · SALTA 6k+1 · BAILA 6k-1 · ZYPYZAPE\n"
            f"🪗 Acordeón: serie de {PAIR} que codifican números posibles"
        )
        tk.Label(parent, text=help_txt, justify=tk.LEFT, font=("Consolas", 9)).pack(
            anchor="w", pady=10
        )

    def _build_viz(self, parent: tk.Frame) -> None:
        box = tk.LabelFrame(
            parent, text=f"📊 Visualización en Tiempo Real {PAIR}", padx=6, pady=6
        )
        box.pack(fill=tk.BOTH, expand=True)
        self.out = tk.Text(
            box,
            wrap=tk.WORD,
            font=("Consolas", 10),
            bg="#111111",
            fg="#e8e8e8",
            insertbackground="white",
        )
        sb = tk.Scrollbar(box, command=self.out.yview)
        self.out.config(yscrollcommand=sb.set)
        self.out.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._log(f"=== Sistema listo {PAIR} ===\n")

    def _build_k3(self, parent: tk.Frame) -> None:
        box = tk.LabelFrame(
            parent, text=f"🔺 K3: Encriptación Geométrica {PAIR}", padx=8, pady=8
        )
        box.pack(fill=tk.BOTH, expand=True)

        top = tk.Frame(box)
        top.pack(fill=tk.X)
        tk.Label(top, text="Longitud secuencia:").pack(side=tk.LEFT)
        self.k3_len = tk.Spinbox(top, from_=1, to=50, width=5)
        self.k3_len.delete(0, tk.END)
        self.k3_len.insert(0, "10")
        self.k3_len.pack(side=tk.LEFT, padx=6)
        tk.Button(top, text=f"Generar secuencia K3 {PAIR}", command=self.generate_k3).pack(
            side=tk.LEFT, padx=6
        )

        self.k3_out = tk.Text(box, wrap=tk.WORD, font=("Consolas", 10), height=20)
        self.k3_out.pack(fill=tk.BOTH, expand=True, pady=6)
        self.k3_out.insert(
            tk.END,
            f"Reglas: 12→+12, 1→+10, 11→+2\n"
            f"K={K_HEART} factor_msl, L={L_ROSE} factor_lsl\n\n"
            f"Pulsa el botón para generar la secuencia {PAIR}.\n",
        )

    def _build_stats(self, parent: tk.Frame) -> None:
        box = tk.LabelFrame(parent, text=f"📈 Estadísticas del Sistema {PAIR}", padx=8, pady=8)
        box.pack(fill=tk.BOTH, expand=True)
        tk.Button(box, text=f"🔄 Actualizar {PAIR}", command=self.refresh_stats).pack(
            anchor="w", pady=4
        )
        self.stats_out = tk.Text(box, wrap=tk.WORD, font=("Consolas", 10), height=22)
        self.stats_out.pack(fill=tk.BOTH, expand=True)
        self.refresh_stats()

    def _build_info(self, parent: tk.Frame) -> None:
        box = tk.LabelFrame(parent, text=f"ℹ️ Documentación {PAIR}", padx=8, pady=8)
        box.pack(fill=tk.BOTH, expand=True)
        info = tk.Text(box, wrap=tk.WORD, font=("Segoe UI Emoji", 10))
        info.pack(fill=tk.BOTH, expand=True)
        info.insert(
            tk.END,
            f"""K3 + Aleatorovix {PAIR}
========================

Origen
------
• organismo_lila_v99.c — Criba Desmemoriada, máscara Lila, intérprete mutante
• gemini-code-1784158392232.py — acordeón ❤️ factor_msl / 🌹 factor_lsl
• Uso civil y educativo (proyecto 33×1 / VMA)

Símbolos
--------
• {K_HEART} = K = factor_msl (Most Significant Level / amplitud)
• {L_ROSE} = L = factor_lsl (Least Significant Level / frecuencia)

Reglas K3 (secuencia de demostración)
-------------------------------------
• Si n == 12 → n + 12
• Si n == 1  → n + 10
• Si n == 11 → n + 2
• En otro caso → n sin cambio

Acciones Lila (decisión 0..3)
-----------------------------
0 OLVIDO   — Criba Desmemoriada (olvido / wipe lógico)
1 SALTA    — candidato 6k+1
2 BAILA    — candidato 6k-1
3 ZYPYZAPE — resonancia con resto k % 121

Ejecutar
--------
  python k3_aleatorovix_gui.py

No requiere dependencias externas (solo Python 3 + tkinter).
""",
        )
        info.config(state=tk.DISABLED)

    # ----- actions -----

    def _toggle_ping(self) -> None:
        self.organismo.entropia.usar_ping = self.usar_ping.get()

    def _log(self, text: str) -> None:
        self.out.insert(tk.END, text)
        self.out.see(tk.END)

    def _update_status(self) -> None:
        n = self.organismo.ciclos_completados
        p = len([x for x in self.organismo.primos_encontrados if AccionesLila._es_primo(x)])
        self.status.config(
            text=f"Ciclos: {n} {K_HEART}  |  Candidatos: {len(self.organismo.primos_encontrados)}  |  "
            f"Primos únicos: {p} {L_ROSE}  |  Auto={'ON' if self.running else 'OFF'} {PAIR}"
        )

    def generate_k3(self) -> None:
        try:
            length = int(self.k3_len.get())
        except Exception:
            length = 10
        length = max(1, min(50, length))
        seq = self.organismo.generar_secuencia_k3(length)

        lines = [
            f"=== SECUENCIA K3 {PAIR} (lectura floral, no arábiga) ===",
            "",
            "Reglas: 12→+12, 1→+10, 11→+2",
            f"K={K_HEART} factor_msl, L={L_ROSE} factor_lsl",
            f"Alfabeto: {K_HEART}=1  {L_ROSE}=0  |  dígitos en nibbles de 4 flores",
            "",
        ]
        for i, (orig, transformed, binary) in enumerate(seq):
            lines.append(
                f"N={entero_a_flores(i, 4)}: "
                f"orig {numero_floral(orig)} → K3 {numero_floral(transformed)} "
                f"→ {bits_a_flores(binary)}"
            )
        text = "\n".join(lines) + "\n"

        self.k3_out.delete("1.0", tk.END)
        self.k3_out.insert(tk.END, text)
        self._log(text + "\n")
        self._show_page("k3")
        self._update_status()

    def run_cycle(self) -> None:
        r = self.organismo.ciclo_completo()
        n = self.organismo.ciclos_completados

        lines = [
            f"Ciclo {numero_floral(n)} {PAIR}",
            formatear_valor_floral("inercia(a)", r["inercia_a"]),
            formatear_valor_floral("red(x)", r["red_x"]),
            formatear_valor_floral("ping_us", r["ping_us"]),
        ]
        if self.detallado.get():
            lines += [
                formatear_valor_floral("medida", r["medida"]),
                f"  curva(fase): {entero_a_flores(abs(int(r['curva'])) % 1024, 10)}",
                f"  eventos: r0={r['decision_floral'] and entero_a_flores(r['r0'], 2)} "
                f"r5={entero_a_flores(r['r5'], 2)} r9={entero_a_flores(r['r9'], 2)}",
                f"  k_semilla floral: {r['k_floral']}",
                f"  k_semilla bin:    {r['k_bin_floral']}",
                f"  chorro acordeón:  {r['flores_acordeon']}",
                f"  número posible:   {r['flores_posible']}  (= {r['num_posible']} interno)",
                f"  {K_HEART} msl  /  {L_ROSE} lsl",
                "  --- pliegue acordeón ---",
                r["acordeon_vista"],
            ]
        if r["candidato_floral"]:
            tag = f" {L_ROSE}PRIMO{L_ROSE}" if r["primo"] else ""
            lines.append(f"  candidato: {r['candidato_floral']}{tag}")
        lines.append(f"  decisión: {r['decision_floral']}  {r['mensaje']}")
        lines.append(f">>> CRIBA DESMEMORIADA {PAIR}")
        lines.append("-" * 50)

        self._log("\n".join(lines) + "\n")
        self._show_page("viz")
        self._update_status()
        self.refresh_stats(silent=True)

    def run_acordeon_serie(self) -> None:
        """Serie de rosas/corazones: cada pliegue muestra un número posible (no estándar)."""
        self._show_page("viz")
        self._log(f"\n=== MODO ACORDEÓN {PAIR} — números solo en flores ===\n")
        self._log(f"Leyenda: {K_HEART}=1 (msl/K)   {L_ROSE}=0 (lsl/L)\n")
        self._log("Cada fila es un número posible leído en pliegues, no en dígitos arábigos.\n\n")
        serie = self.organismo.serie_acordeon_flores(pliegues=24, cuantos=10)
        for item in serie:
            self._log(
                f"pliegue {entero_a_flores(item['pliegue'], 4)}  "
                f"ancho={entero_a_flores(item['ancho'], 4)}\n"
            )
            self._log(f"  serie:  {item['flores']}\n")
            self._log(f"  dígitos-flor: {item['floral_digitos']}\n")
            self._log(item["vista"] + "\n")
            self._log("-" * 40 + "\n")
        self._log(f"=== fin serie acordeón {PAIR} ===\n")
        self._update_status()

    def run_ventanas_moviles(self) -> None:
        """Ventanas móviles sobre un chorro largo → varios números ❤️🌹."""
        self._show_page("viz")
        total_bits, num_ventanas = 128, 8
        self._log(f"\n=== VENTANAS MÓVILES {PAIR} ===\n")
        self._log(
            f"Chorro base: {entero_a_flores(total_bits)} bits | "
            f"ventanas: {entero_a_flores(num_ventanas)}\n\n"
        )
        data = self.organismo.acordeon.chorro_ventanas_moviles(
            total_bits=total_bits,
            num_ventanas=num_ventanas,
            min_ancho=8,
            max_ancho=32,
        )
        for w in data:
            self._log(
                f"ventana {entero_a_flores(w['ventana'], 4)}  "
                f"[{entero_a_flores(w['start_idx'], 8)}…{entero_a_flores(w['end_idx'], 8)}]  "
                f"ancho={entero_a_flores(w['ancho'], 6)}\n"
            )
            self._log(f"  flores: {w['flores']}\n")
            self._log(f"  dígitos-flor: {w['floral_digitos']}\n")
            self._log(w["vista"] + "\n")
            self._log("-" * 40 + "\n")
        self._log(f"=== fin ventanas {PAIR} ===\n")
        self._update_status()

    def run_industrial(self) -> None:
        """Modo industrial: mural denso de miles de ❤️🌹 (escala std en GUI)."""
        self._show_page("viz")
        self._log(f"\n=== 🏭 PLANTA INDUSTRIAL {PAIR} ===\n")
        try:
            from k3_aleatorovix_industrial import PlantaFloral, banner
        except ImportError:
            self._log("No se encontró k3_aleatorovix_industrial.py en el mismo directorio.\n")
            return
        self._log(banner() + "\n\n")
        planta = PlantaFloral(semilla=33)
        # std: ~32×64 = 2048 flores — denso pero usable en Text widget
        filas, inf = planta.producir_mural(escala="std", motor="dual")
        for row in filas:
            self._log(row + "\n")
        self._log("\n" + inf.resumen_floral() + "\n")
        # mini lluvia extra
        lluvia, inf2 = planta.producir_lluvia(n_flores=2000, por_linea=64)
        self._log(f"\n--- lluvia extra {PAIR} ---\n")
        for row in lluvia[:20]:
            self._log(row + "\n")
        self._log(inf2.resumen_floral() + "\n")
        total_h = inf.corazones + inf2.corazones
        total_r = inf.rosas + inf2.rosas
        self._log(
            f"\nTOTAL GUI: {K_HEART}{total_h} + {L_ROSE}{total_r} = {total_h + total_r} flores\n"
        )
        self._log(f"=== fin industrial {PAIR} ===\n")
        self._update_status()

    def toggle_auto(self) -> None:
        self.running = not self.running
        if self.running:
            self._auto_tick()
        self._update_status()

    def _auto_tick(self) -> None:
        if not self.running:
            return
        self.run_cycle()
        delay = int(self.speed_scale.get())
        self.root.after(delay, self._auto_tick)

    def reset_all(self) -> None:
        self.running = False
        self.organismo.reiniciar()
        self.out.delete("1.0", tk.END)
        self._log(f"=== Reiniciado {PAIR} ===\n")
        self.refresh_stats()
        self._update_status()

    def refresh_stats(self, silent: bool = False) -> None:
        o = self.organismo
        counts = Counter(o.historial_decisiones)
        primos_ok = sorted({c for c in o.primos_encontrados if AccionesLila._es_primo(c)})
        lines = [
            f"=== ESTADÍSTICAS {PAIR} ===",
            "",
            f"General {PAIR}",
            f"  Ciclos completados: {o.ciclos_completados} {K_HEART}",
            f"  Candidatos generados: {len(o.primos_encontrados)}",
            f"  Primos únicos: {len(primos_ok)} {L_ROSE}",
            "",
            "Frecuencia de decisiones:",
        ]
        for d in range(4):
            name = AccionesLila.ACCION_NOMBRES[d]
            lines.append(f"  {name}: {counts.get(d, 0)}")
        lines += [
            "",
            f"Primos encontrados {L_ROSE}:",
            f"  {primos_ok if primos_ok else '(ninguno aún)'}",
            "",
            f"{K_HEART} factor_msl = {o.acordeon.factor_msl}",
            f"{L_ROSE} factor_lsl = {o.acordeon.factor_lsl}",
        ]
        text = "\n".join(lines) + "\n"
        self.stats_out.delete("1.0", tk.END)
        self.stats_out.insert(tk.END, text)
        if not silent:
            pass


def main() -> None:
    root = tk.Tk()
    # Fuente amigable con emoji en Windows
    try:
        root.option_add("*Font", "Segoe UI Emoji 10")
    except Exception:
        pass
    K3AleatorovixGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
