#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
═══════════════════════════════════════════════════════════════
  ALEATOROVIX UNO  ❤️🌹
  Un solo archivo · fácil de leer · funciona en Colab y en PC
═══════════════════════════════════════════════════════════════

  ¿Qué es?
    Un generador de azar “vivo”: cada bit sale como flor.
      ❤️ = 1  (corazón · fuerza · K)
      🌹 = 0  (rosa · matiz · L)

  ¿Por qué no tirar una moneda “yo solo en la cabeza”?
    Si inventas números de memoria, repites hábitos (fechas, 7, 3…).
    Aquí el azar se mezcla de reloj + pila de memoria + una fórmula
    que “respira” (acordeón). Dos personas con la misma semilla
    pueden REPETIR el mismo chorro (útil para compartir un juego,
    una clave de sesión, un sorteo verificable). Eso es mejor que
    “me lo inventé yo” — es azar + rastro + posible acuerdo.

  Uso:
    python aleatorovix_uno.py              # modo INTERACTIVO (ON + Generar)
    python aleatorovix_uno.py demo         # explicación + demo de texto
    python aleatorovix_uno.py corto
    python aleatorovix_uno.py mural
    python aleatorovix_uno.py lluvia 2000
    python aleatorovix_uno.py fijo 33

  Interactivo:
    · Botón ON  → el ramo crece solo, poco a poco (más ❤️ y 🌹)
    · Botón OFF → se detiene el crecimiento
    · GENERAR   → los índices parpadean rápido varias veces y se fijan
    · Reiniciar → ramo vacío + nueva semilla

  En Colab (sin ventana):  !python aleatorovix_uno.py consola
  Uso civil / educativo. No es “seguridad bancaria”.
"""

from __future__ import annotations

import math
import random
import sys
import time
from decimal import Decimal, getcontext
from typing import List, Optional

getcontext().prec = 60

# ── Símbolos (cámbialos si quieres otro “alfabeto”) ───────────
CORAZON = "❤️"   # bit 1
ROSA = "🌹"      # bit 0
PAREJA = CORAZON + ROSA


# ═════════════════════════════════════════════════════════════
#  Motor mínimo: un pliegue del acordeón = una flor
# ═════════════════════════════════════════════════════════════

class Acordeon:
    """
    Estado que se estira y se pliega (como un acordeón).
    factor_msl (K/❤️) = amplitud grande
    factor_lsl (L/🌹) = temblor fino
    """

    def __init__(self, semilla: int = 0):
        # Si no das semilla: reloj + basura de memoria → distinto cada vez
        if semilla == 0:
            basura = id([]) ^ int(time.time_ns())
            semilla = (basura % 999_983) or 33
        self.estado = Decimal(str(semilla))
        self.msl = Decimal("42.123456789")   # ❤️ K
        self.lsl = Decimal("0.000123456")    # 🌹 L
        self.n_corazones = 0
        self.n_rosas = 0
        self.semilla_usada = semilla

    def un_bit(self) -> str:
        # Evolución: fase caótica determinista (misma semilla → mismo chorro)
        temblor = Decimal(str(math.cos(float(self.estado * self.lsl))))
        fase = (self.estado * self.msl) + temblor
        self.estado = fase - fase.to_integral_value(rounding="ROUND_FLOOR")
        return "1" if self.estado > Decimal("0.5") else "0"

    def una_flor(self) -> str:
        b = self.un_bit()
        if b == "1":
            self.n_corazones += 1
            return CORAZON
        self.n_rosas += 1
        return ROSA

    def ramo(self, cuantas: int) -> str:
        """Devuelve una ristra de flores."""
        return "".join(self.una_flor() for _ in range(cuantas))

    def numero_desde_flores(self, bits: int = 8) -> tuple[int, str]:
        """Lee N bits del acordeón y los muestra como flores + valor."""
        cadena_bits = "".join(self.un_bit() for _ in range(bits))
        for b in cadena_bits:
            if b == "1":
                self.n_corazones += 1
            else:
                self.n_rosas += 1
        flores = "".join(CORAZON if b == "1" else ROSA for b in cadena_bits)
        return int(cadena_bits, 2), flores


# ═════════════════════════════════════════════════════════════
#  Entropía simple: “por qué no soy yo solo inventando”
# ═════════════════════════════════════════════════════════════

def semilla_del_mundo() -> int:
    """
    Mezcla tres fuentes que TÚ no controlas del todo:
      1) nanosegundos del reloj
      2) un id de objeto en memoria (dirección “sucia”)
      3) un random del sistema
    Resultado: número de arranque distinto casi siempre.
    Si quieres REPETIR el experimento, pasa semilla fija a Acordeon(123).
    """
    t = int(time.time_ns()) % 10_000_000
    sucio = id(object()) % 10_000_000
    r = random.randint(1, 10_000_000)
    return (t ^ sucio ^ r) or 33


# ═════════════════════════════════════════════════════════════
#  Pantallas (TODOS los print van con comentario para que los toques)
# ═════════════════════════════════════════════════════════════

def pantalla_titulo() -> None:
    # [TOCAR] cabecera del programa
    print(PAREJA * 16)
    # [TOCAR] nombre que quieras darle al proyecto
    print(f"  ALEATOROVIX UNO  {PAREJA}")
    # [TOCAR] subtítulo
    print("  Azar con rastro · compartible · no solo 'me lo inventé'")
    # [TOCAR] cierre del banner
    print(PAREJA * 16)
    # [TOCAR] línea en blanco después del título
    print()


def pantalla_por_que() -> None:
    # [TOCAR] bloque explicativo — reescribe con tus palabras
    print("── ¿Por qué esto y no un número de mi cabeza? ──")
    # [TOCAR]
    print("  • Si eliges tú solo, repites manías (cumpleaños, 7, 3…).")
    # [TOCAR]
    print("  • Aquí el pliegue del acordeón + el reloj hacen un chorro")
    # [TOCAR]
    print("    que no 'piensa' como tú: es más justo en un sorteo.")
    # [TOCAR]
    print("  • Si dos personas usan la MISMA semilla, ven las MISMAS flores:")
    # [TOCAR]
    print("    se puede comprobar. Eso es mejor para un 'alguien' que confía")
    # [TOCAR]
    print("    en ti: no es magia opaca, es un rastro que se puede repetir.")
    # [TOCAR]
    print("  • ❤️ y 🌹 no son adorno: son el número en otro alfabeto")
    # [TOCAR]
    print("    (1 y 0). Leer flores = leer azar sin dígitos arábigos.")
    # [TOCAR]
    print()


def pantalla_leyenda() -> None:
    # [TOCAR] leyenda de símbolos
    print("── Leyenda ──")
    # [TOCAR]
    print(f"  {CORAZON}  = bit 1  = K  = empuje grande (msl)")
    # [TOCAR]
    print(f"  {ROSA}  = bit 0  = L  = detalle fino (lsl)")
    # [TOCAR]
    print()


def pantalla_semilla(semilla: int, modo: str) -> None:
    # [TOCAR] cómo se eligió el arranque
    print("── Semilla de arranque ──")
    # [TOCAR] muestra el número (interno; las flores son la cara pública)
    print(f"  modo: {modo}")
    # [TOCAR]
    print(f"  semilla interna: {semilla}")
    # [TOCAR] misma semilla en flores (binario corto)
    bits = bin(semilla)[2:][-16:].zfill(16)
    flores = "".join(CORAZON if b == "1" else ROSA for b in bits)
    print(f"  semilla en flores (16 bits): {flores}")
    # [TOCAR]
    print()


def pantalla_ramo(flores: str, titulo: str = "Ramo") -> None:
    # [TOCAR] título del ramo
    print(f"── {titulo} ──")
    # [TOCAR] la serie de corazones y rosas (esto es lo bonito)
    print(f"  {flores}")
    # [TOCAR]
    print()


def pantalla_numeros_posibles(acordeon: Acordeon, cuantos: int = 5) -> None:
    # [TOCAR] cabecera de números “posibles”
    print("── Números posibles (cada uno es un pliegue del acordeón) ──")
    # [TOCAR]
    print("  (se leen en flores; el número entre paréntesis es solo para ti)")
    # [TOCAR]
    print()
    for i in range(1, cuantos + 1):
        valor, flores = acordeon.numero_desde_flores(bits=8)
        # [TOCAR] una línea por número — cambia el texto a tu gusto
        print(f"  pliegue {i}: {flores}   ←→  ({valor})")
    # [TOCAR]
    print()


def pantalla_mural(acordeon: Acordeon, lineas: int = 12, ancho: int = 40) -> None:
    # [TOCAR] cabecera del mural
    print("── Mural (pared de azar) ──")
    # [TOCAR]
    print()
    for _ in range(lineas):
        # [TOCAR] cada fila del mural (no hace falta tocar el print: es la fila)
        print(acordeon.ramo(ancho))
    # [TOCAR]
    print()


def pantalla_lluvia(acordeon: Acordeon, total: int = 500, por_fila: int = 40) -> None:
    # [TOCAR] cabecera lluvia
    print(f"── Lluvia de {total} flores ──")
    # [TOCAR]
    print()
    hechas = 0
    while hechas < total:
        n = min(por_fila, total - hechas)
        # [TOCAR] cada racha de la lluvia
        print(acordeon.ramo(n))
        hechas += n
    # [TOCAR]
    print()


def pantalla_pulso(flores: str) -> None:
    # [TOCAR] el acordeón se abre y se cierra (visual)
    print("── Pulso del acordeón (cierra → abre → cierra) ──")
    # Separar emojis bien (❤️ y 🌹 no son 1 carácter “simple”)
    lista = []
    i = 0
    while i < len(flores):
        if flores.startswith(CORAZON, i):
            lista.append(CORAZON)
            i += len(CORAZON)
        elif flores.startswith(ROSA, i):
            lista.append(ROSA)
            i += len(ROSA)
        else:
            i += 1
    # [TOCAR] cerrado
    print(f"  cerrado:  {''.join(lista)}")
    # [TOCAR] semi abierto
    print(f"  semi:     {' '.join(lista)}")
    # [TOCAR] bien abierto
    print(f"  abierto:  {'  '.join(lista)}")
    # [TOCAR] vuelve a cerrar
    print(f"  cerrado:  {''.join(lista)}")
    # [TOCAR]
    print()


def pantalla_conteo(acordeon: Acordeon) -> None:
    total = acordeon.n_corazones + acordeon.n_rosas
    # [TOCAR] resumen final
    print("── Conteo ──")
    # [TOCAR]
    print(f"  {CORAZON} corazones: {acordeon.n_corazones}")
    # [TOCAR]
    print(f"  {ROSA} rosas:     {acordeon.n_rosas}")
    # [TOCAR]
    print(f"  total flores:  {total}")
    if total:
        # [TOCAR]
        print(f"  proporción ❤️: {acordeon.n_corazones / total:.3f}")
    # [TOCAR]
    print()


def pantalla_cierre() -> None:
    # [TOCAR] despedida
    print("── Fin ──")
    # [TOCAR] invita a tocar los prints
    print("  Tip: busca en este archivo la marca  [TOCAR]  y cambia los textos.")
    # [TOCAR]
    print(f"  Interactivo:        python aleatorovix_uno.py")
    # [TOCAR]
    print(f"  Demo texto:         python aleatorovix_uno.py demo")
    # [TOCAR]
    print(f"  Mural grande:       python aleatorovix_uno.py mural")
    # [TOCAR]
    print(f"  Muchas flores:      python aleatorovix_uno.py lluvia 3000")
    # [TOCAR]
    print(f"  Semilla fija:       python aleatorovix_uno.py fijo 33")
    # [TOCAR]
    print()
    # [TOCAR]
    print(PAREJA * 8)


# ═════════════════════════════════════════════════════════════
#  MODO INTERACTIVO (ventana) — ON crece el ramo · GENERAR agita índices
# ═════════════════════════════════════════════════════════════

class InteractivoRamo:
    """
    Ventana simple:
      · ON  → cada poco se añade 1–3 flores al ramo (crece solo)
      · GENERAR → los índices (números posibles) cambian rápido N veces
    Todos los textos de la UI llevan marca [TOCAR] en el código abajo.
    """

    # Ritmo: cada cuántos ms crece el ramo estando ON
    MS_CRECER = 400
    # Cuántas flores nuevas por tick (mín, máx)
    FLORES_POR_TICK = (1, 3)
    # Tope del ramo (para no colgar la ventana)
    MAX_FLORES_RAMO = 400
    # Al generar: cuántas veces parpadean los índices
    PARPADEOS = 12
    # ms entre parpadeos
    MS_PARPADEO = 55
    # Cuántos índices se muestran
    N_INDICES = 6

    def __init__(self, semilla: Optional[int] = None):
        import tkinter as tk
        from tkinter import font as tkfont

        self.tk = tk
        self.root = tk.Tk()
        # [TOCAR] título de la ventana
        self.root.title(f"Aleatorovix UNO  {PAREJA}  — interactivo")
        self.root.geometry("720x560")
        self.root.minsize(560, 420)

        self.semilla = semilla if semilla is not None else semilla_del_mundo()
        self.ac = Acordeon(self.semilla)
        self.ramo_lista: List[str] = []  # flores del ramo (crece con ON)
        self.encendido = False
        self.generando = False
        self.indices: List[tuple[int, str]] = []  # (valor, flores)
        self._job_crecer = None
        self._job_gen = None
        self._parpadeo_restante = 0

        # —— UI ——
        pad = {"padx": 10, "pady": 6}
        # [TOCAR] cabecera
        self.lbl_titulo = tk.Label(
            self.root,
            text=f"ALEATOROVIX UNO  {PAREJA}",
            font=("Segoe UI Emoji", 16, "bold"),
        )
        self.lbl_titulo.pack(pady=(12, 2))

        # [TOCAR] instrucción
        self.lbl_ayuda = tk.Label(
            self.root,
            text="ON = el ramo crece solo  ·  GENERAR = los índices saltan rápido  ·  OFF = pausa",
            font=("Segoe UI", 9),
            fg="#444",
        )
        self.lbl_ayuda.pack()

        # [TOCAR] estado ON/OFF
        self.lbl_estado = tk.Label(
            self.root,
            text="Estado: OFF  (pulsa ON para crecer)",
            font=("Segoe UI", 11, "bold"),
            fg="#a00",
        )
        self.lbl_estado.pack(pady=4)

        # Botones
        fila = tk.Frame(self.root)
        fila.pack(pady=8)

        # [TOCAR] textos de botones
        self.btn_on = tk.Button(fila, text="▶  ON", width=12, command=self.poner_on, bg="#d4edda")
        self.btn_on.pack(side=tk.LEFT, padx=4)
        self.btn_off = tk.Button(fila, text="■  OFF", width=12, command=self.poner_off, bg="#f8d7da")
        self.btn_off.pack(side=tk.LEFT, padx=4)
        self.btn_gen = tk.Button(
            fila, text=f"⚡ GENERAR {PAREJA}", width=16, command=self.generar_indices, bg="#fff3cd"
        )
        self.btn_gen.pack(side=tk.LEFT, padx=4)
        self.btn_reset = tk.Button(fila, text="↺ Reiniciar", width=12, command=self.reiniciar)
        self.btn_reset.pack(side=tk.LEFT, padx=4)

        # Ramo (crece)
        # [TOCAR] marco del ramo
        marco_ramo = tk.LabelFrame(self.root, text=f"Ramo vivo {PAREJA}  (crece con ON)", padx=8, pady=8)
        marco_ramo.pack(fill=tk.BOTH, expand=True, padx=12, pady=6)

        self.txt_ramo = tk.Text(
            marco_ramo,
            height=8,
            wrap=tk.WORD,
            font=("Segoe UI Emoji", 14),
            bg="#111",
            fg="#f5f5f5",
        )
        self.txt_ramo.pack(fill=tk.BOTH, expand=True)
        # [TOCAR] mensaje inicial del ramo
        self.txt_ramo.insert("1.0", "(vacío — pulsa ON y verás nacer ❤️ y 🌹…)\n")

        self.lbl_tamano = tk.Label(marco_ramo, text="Flores en el ramo: 0", font=("Segoe UI", 9))
        self.lbl_tamano.pack(anchor="w")

        # Índices
        # [TOCAR] marco de índices
        marco_idx = tk.LabelFrame(
            self.root, text="Índices (números posibles) — GENERAR los agita", padx=8, pady=8
        )
        marco_idx.pack(fill=tk.BOTH, expand=False, padx=12, pady=6)

        self.txt_idx = tk.Text(
            marco_idx,
            height=8,
            wrap=tk.WORD,
            font=("Consolas", 11),
            bg="#1a1a2e",
            fg="#eaeaea",
        )
        self.txt_idx.pack(fill=tk.BOTH, expand=True)
        # [TOCAR]
        self.txt_idx.insert("1.0", "Pulsa GENERAR para ver índices saltar…\n")

        self.lbl_semilla = tk.Label(
            self.root,
            text=f"Semilla: {self.semilla}",
            font=("Segoe UI", 8),
            fg="#666",
        )
        self.lbl_semilla.pack(pady=(0, 8))

        # Primera tira de índices en reposo
        self._rellenar_indices(animando=False)

        self.root.protocol("WM_DELETE_WINDOW", self._cerrar)

    def _cerrar(self) -> None:
        self.poner_off()
        if self._job_gen:
            try:
                self.root.after_cancel(self._job_gen)
            except Exception:
                pass
        self.root.destroy()

    def poner_on(self) -> None:
        if self.encendido:
            return
        self.encendido = True
        # [TOCAR] texto cuando está ON
        self.lbl_estado.config(text="Estado: ON  ●  el ramo crece…", fg="#080")
        self._tick_crecer()

    def poner_off(self) -> None:
        self.encendido = False
        if self._job_crecer is not None:
            try:
                self.root.after_cancel(self._job_crecer)
            except Exception:
                pass
            self._job_crecer = None
        # [TOCAR] texto cuando está OFF
        self.lbl_estado.config(text="Estado: OFF  (pausado)", fg="#a00")

    def _tick_crecer(self) -> None:
        if not self.encendido:
            return
        if len(self.ramo_lista) >= self.MAX_FLORES_RAMO:
            # [TOCAR] tope alcanzado
            self.lbl_estado.config(
                text=f"Estado: ON · ramo al máximo ({self.MAX_FLORES_RAMO}) — GENERAR sigue activo",
                fg="#060",
            )
        else:
            n = random.randint(*self.FLORES_POR_TICK)
            n = min(n, self.MAX_FLORES_RAMO - len(self.ramo_lista))
            for _ in range(n):
                self.ramo_lista.append(self.ac.una_flor())
            self._pintar_ramo()
        self._job_crecer = self.root.after(self.MS_CRECER, self._tick_crecer)

    def _pintar_ramo(self) -> None:
        texto = "".join(self.ramo_lista)
        self.txt_ramo.delete("1.0", self.tk.END)
        self.txt_ramo.insert("1.0", texto if texto else "(vacío)\n")
        # [TOCAR] contador de tamaño
        h = self.ramo_lista.count(CORAZON)
        r = self.ramo_lista.count(ROSA)
        self.lbl_tamano.config(
            text=f"Flores en el ramo: {len(self.ramo_lista)}   "
            f"{CORAZON} {h}   {ROSA} {r}"
        )

    def _rellenar_indices(self, animando: bool = False) -> None:
        self.indices = []
        lineas = []
        for i in range(1, self.N_INDICES + 1):
            valor, flores = self.ac.numero_desde_flores(8)
            self.indices.append((valor, flores))
            # [TOCAR] formato de cada índice en pantalla
            marca = "  …" if animando else ""
            lineas.append(f"  índice {i}:  {flores}   ←→  ({valor}){marca}")
        self.txt_idx.delete("1.0", self.tk.END)
        if animando:
            # [TOCAR] cabecera mientras agita
            self.txt_idx.insert("1.0", f"⚡ agitando índices… {PAREJA}\n\n")
        else:
            # [TOCAR] cabecera en reposo
            self.txt_idx.insert("1.0", f"Índices fijados {PAREJA}\n\n")
        self.txt_idx.insert(self.tk.END, "\n".join(lineas) + "\n")

    def generar_indices(self) -> None:
        """Los índices cambian rápido unas cuantas veces y se quedan."""
        if self.generando:
            return
        self.generando = True
        self.btn_gen.config(state=self.tk.DISABLED)
        self._parpadeo_restante = self.PARPADEOS
        self._tick_generar()

    def _tick_generar(self) -> None:
        if self._parpadeo_restante <= 0:
            self._rellenar_indices(animando=False)
            self.generando = False
            self.btn_gen.config(state=self.tk.NORMAL)
            # [TOCAR] aviso al terminar de generar
            self.lbl_estado.config(
                text=("Estado: ON  ●  índices listos" if self.encendido else "Estado: OFF  ·  índices listos"),
                fg="#080" if self.encendido else "#a00",
            )
            return
        self._rellenar_indices(animando=True)
        self._parpadeo_restante -= 1
        self._job_gen = self.root.after(self.MS_PARPADEO, self._tick_generar)

    def reiniciar(self) -> None:
        self.poner_off()
        self.semilla = semilla_del_mundo()
        self.ac = Acordeon(self.semilla)
        self.ramo_lista.clear()
        self._pintar_ramo()
        # [TOCAR] ramo vacío tras reinicio
        self.txt_ramo.delete("1.0", self.tk.END)
        self.txt_ramo.insert("1.0", "(vacío — pulsa ON…)\n")
        self.lbl_semilla.config(text=f"Semilla: {self.semilla}")
        self._rellenar_indices(animando=False)
        # [TOCAR]
        self.lbl_estado.config(text="Estado: OFF  (reiniciado)", fg="#a00")

    def run(self) -> None:
        self.root.mainloop()


def modo_interactivo(semilla: Optional[int] = None) -> None:
    try:
        app = InteractivoRamo(semilla=semilla)
    except Exception as e:
        # [TOCAR] si no hay ventana (Colab, etc.)
        print(f"No se pudo abrir la ventana interactiva: {e}")
        print("Uso modo consola en su lugar…\n")
        modo_consola_interactivo()
        return
    app.run()


def modo_consola_interactivo() -> None:
    """
    Sin GUI (Colab / SSH):
      on     → crece el ramo en pasos (Enter entre pasos o auto)
      gen    → agita índices varias veces
      off / q
    """
    # [TOCAR]
    print(PAREJA * 12)
    # [TOCAR]
    print("  MODO CONSOLA INTERACTIVO")
    # [TOCAR]
    print("  comandos:  on  |  paso  |  gen  |  ver  |  reset  |  q")
    # [TOCAR]
    print(PAREJA * 12)
    ac = Acordeon(semilla_del_mundo())
    ramo: List[str] = []
    # [TOCAR]
    print(f"  Semilla: {ac.semilla_usada}")
    # [TOCAR]
    print()

    def ver() -> None:
        # [TOCAR]
        print(f"  Ramo ({len(ramo)}): {''.join(ramo) if ramo else '(vacío)'}")
        # [TOCAR]
        print(f"  {CORAZON}{ramo.count(CORAZON)}  {ROSA}{ramo.count(ROSA)}")

    def gen() -> None:
        # [TOCAR]
        print("  ⚡ agitando índices…")
        for n in range(8):
            lineas = []
            for i in range(1, 5):
                v, f = ac.numero_desde_flores(8)
                lineas.append(f"    índice {i}: {f} ({v})")
            # [TOCAR] parpadeo en consola
            print(f"  · ronda {n + 1}")
            print("\n".join(lineas))
            time.sleep(0.08)
        # [TOCAR]
        print("  → índices fijados:")
        for i in range(1, 5):
            v, f = ac.numero_desde_flores(8)
            # [TOCAR]
            print(f"    índice {i}: {f}  ←→  ({v})")

    while True:
        try:
            # [TOCAR] prompt
            cmd = input("  > ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            # [TOCAR]
            print("\n  adiós")
            break
        if cmd in ("q", "quit", "salir", "exit"):
            # [TOCAR]
            print("  adiós")
            break
        if cmd == "on":
            # [TOCAR] crece en ráfaga de varios pasos
            print("  ON: creciendo 15 pasos (poco a poco)…")
            for _ in range(15):
                for __ in range(random.randint(1, 3)):
                    ramo.append(ac.una_flor())
                ver()
                time.sleep(0.25)
            # [TOCAR]
            print("  (fin de ráfaga ON — escribe on otra vez o gen)")
        elif cmd in ("paso", "p", "+"):
            for _ in range(random.randint(1, 3)):
                ramo.append(ac.una_flor())
            ver()
        elif cmd in ("gen", "generar", "g"):
            gen()
        elif cmd in ("ver", "v", "ramo"):
            ver()
        elif cmd in ("reset", "reiniciar", "r"):
            ac = Acordeon(semilla_del_mundo())
            ramo.clear()
            # [TOCAR]
            print(f"  reiniciado · semilla {ac.semilla_usada}")
        else:
            # [TOCAR]
            print("  on | paso | gen | ver | reset | q")


# ═════════════════════════════════════════════════════════════
#  Flujos (qué se muestra en cada modo)
# ═════════════════════════════════════════════════════════════

def demo_normal(semilla: int | None = None) -> None:
    """Demo fácil de entender: título → por qué → ramo → números → pulso → conteo."""
    pantalla_titulo()
    pantalla_por_que()
    pantalla_leyenda()

    if semilla is None:
        s = semilla_del_mundo()
        modo = "mundo (reloj + memoria + azar del sistema)"
    else:
        s = semilla
        modo = "fija (repetible: misma semilla → mismas flores)"

    pantalla_semilla(s, modo)
    ac = Acordeon(s)

    ramo = ac.ramo(24)
    pantalla_ramo(ramo, titulo="Primer ramo (24 flores de azar)")
    pantalla_pulso(ramo)
    pantalla_numeros_posibles(ac, cuantos=6)
    pantalla_conteo(ac)
    pantalla_cierre()


def demo_corto() -> None:
    """Muy corto: solo un ramo y un número."""
    pantalla_titulo()
    # [TOCAR] modo corto
    print("── Modo corto ──")
    # [TOCAR]
    print()
    ac = Acordeon(semilla_del_mundo())
    # [TOCAR]
    print(f"  {ac.ramo(32)}")
    # [TOCAR]
    print()
    v, f = ac.numero_desde_flores(8)
    # [TOCAR]
    print(f"  un número posible: {f}  ({v})")
    # [TOCAR]
    print()
    pantalla_conteo(ac)


def demo_mural() -> None:
    pantalla_titulo()
    pantalla_por_que()
    ac = Acordeon(semilla_del_mundo())
    pantalla_mural(ac, lineas=16, ancho=48)
    pantalla_conteo(ac)
    pantalla_cierre()


def demo_lluvia(n: int = 2000) -> None:
    pantalla_titulo()
    ac = Acordeon(semilla_del_mundo())
    pantalla_lluvia(ac, total=n, por_fila=50)
    pantalla_conteo(ac)
    pantalla_cierre()


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    # Por defecto: interactivo (ventana ON / GENERAR)
    if not argv:
        modo_interactivo()
        return
    cmd = argv[0].lower()
    if cmd in ("on", "ui", "gui", "interactivo", "ventana"):
        s = int(argv[1]) if len(argv) > 1 else None
        modo_interactivo(semilla=s)
    elif cmd in ("consola", "colab", "cli"):
        modo_consola_interactivo()
    elif cmd in ("demo", "texto", "explicar"):
        demo_normal()
    elif cmd in ("corto", "short", "mini"):
        demo_corto()
    elif cmd in ("mural", "pared"):
        demo_mural()
    elif cmd in ("lluvia", "rain"):
        n = int(argv[1]) if len(argv) > 1 else 2000
        demo_lluvia(n)
    elif cmd in ("fijo", "seed", "semilla"):
        s = int(argv[1]) if len(argv) > 1 else 33
        # semilla fija + ventana interactiva
        modo_interactivo(semilla=s)
    elif cmd in ("help", "-h", "--help", "ayuda"):
        # [TOCAR] ayuda
        print("aleatorovix_uno.py")
        print("  (sin args)     ventana interactiva: ON crece · GENERAR agita índices")
        print("  on | gui       igual")
        print("  consola        interactivo por texto (Colab)")
        print("  demo           explicación en texto")
        print("  corto | mural | lluvia [N]")
        print("  fijo [semilla] interactivo con semilla fija")
    else:
        # [TOCAR] comando desconocido
        print(f"No entiendo '{cmd}'. Prueba: on | consola | demo | corto | mural")
        modo_consola_interactivo()


if __name__ == "__main__":
    main()
