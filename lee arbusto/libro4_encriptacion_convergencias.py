#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Libro 4 — Encriptación por Convergencias
==========================================
Método gráfico de restos + 4 semillas + Energía Verde + secuencias de figuras.

Secuencias geométricas (Números i numeritos, pág. 15-16):
  - Figuras: cuadrado → octógono → 16 → 32 → 64 lados (inscritas, r=1)
  - Lados:   4, 8, 16, 32, 64, …
  - Alturas: apotema cos(π/n) por figura

Tras calcular diversas aproximaciones a π, la contraseña guía la elección de
referencias (método π + índice de figura/lados/altura) necesarias para cifrar
y recuperar el mensaje.

Autor: Víctor Manzanares Alberola — Proyecto 33×1
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import hmac
import math
import struct
import sys
from dataclasses import dataclass
from typing import Callable, Iterable

# ---------------------------------------------------------------------------
# Constantes Libro 4
# ---------------------------------------------------------------------------

SEMILLAS_FUNDAMENTALES = (2, 3, 5, 7)
E_VERDE = (math.e - 1) / 2
CONVERGENCIA_S = 0.5
ITER_CONVERGENCIA = 16
MAGICO_K3 = 0x9E3779B9
PI_OBJETIVO = math.pi
HEADER_MAGIC = "L4v3"


# ---------------------------------------------------------------------------
# Secuencias de figuras, lados y alturas (Pitágoras / polígonos inscritos)
# ---------------------------------------------------------------------------

NOMBRES_FIGURAS = {
    4: "cuadrado",
    8: "octógono",
    16: "hexadecágono",
    32: "triacontádigono",
    64: "hexacontatetrágono",
    128: "hectogon",
    256: "256-ágono",
    512: "512-ágono",
}


def generar_secuencias_geometricas(pasos: int = 8) -> dict[str, list]:
    """
    Cuadrado inscrito en círculo r=1; se duplican lados (4→8→16→…).
    Altura = apotema = cos(π/n).  Lado = 2·sin(π/n).
    """
    lados: list[int] = []
    figuras: list[str] = []
    alturas: list[float] = []
    lados_lado: list[float] = []
    perimetros: list[float] = []

    for i in range(pasos):
        n = 4 * (2**i)
        lado = 2.0 * math.sin(math.pi / n)
        altura = math.cos(math.pi / n)
        perimetro = n * lado
        lados.append(n)
        figuras.append(NOMBRES_FIGURAS.get(n, f"{n}-ágono"))
        alturas.append(altura)
        lados_lado.append(lado)
        perimetros.append(perimetro)

    return {
        "figuras": figuras,
        "lados": lados,
        "alturas": alturas,
        "lado": lados_lado,
        "perimetros": perimetros,
    }


# ---------------------------------------------------------------------------
# Métodos diversos de aproximación a π
# ---------------------------------------------------------------------------

def pi_poligono_inscrito(n_lados: int) -> float:
    """π ≈ n·sin(π/n) — perímetro inscrito / 2 con r=1."""
    n = max(4, n_lados)
    return n * math.sin(math.pi / n)


def pi_poligono_circunscrito(n_lados: int) -> float:
    """Cota superior: n·tan(π/n)."""
    n = max(4, n_lados)
    return n * math.tan(math.pi / n)


def pi_archimedes_medio(pasos: int) -> float:
    """Promedio inscrito/circunscrito duplicando lados (estilo Arquímedes)."""
    n = 6
    p_in = p_out = 0.0
    for _ in range(max(1, pasos)):
        p_in = n * math.sin(math.pi / n)
        p_out = n * math.tan(math.pi / n)
        n *= 2
    return (p_in + p_out) / 2.0


def pi_pitagoras_victor(pasos: int) -> float:
    """
    Victor (pág. 15): cuadrado inscrito → octógono → 16… vía Pitágoras.
    Iteración: duplicar lados; lado_nuevo desde geometría del polígono inscrito.
    """
    n = 4
    for _ in range(max(1, pasos)):
        n *= 2
    return n * math.sin(math.pi / n)


def pi_leibniz(terminos: int) -> float:
    """π/4 = 1 − 1/3 + 1/5 − …"""
    s = 0.0
    for k in range(max(1, terminos)):
        s += (-1) ** k / (2 * k + 1)
    return 4.0 * s


def pi_euler(terminos: int) -> float:
    """π²/6 = Σ 1/n²  ⇒  π = √(6·Σ)."""
    s = sum(1.0 / (k * k) for k in range(1, max(2, terminos) + 1))
    return math.sqrt(6.0 * s)


def pi_convergencia_criba(pasos: int, objetivo: float = PI_OBJETIVO) -> float:
    """Tendencia 1/2 hacia π: D_{k+1} = (D_k + objetivo) / 2."""
    d = objetivo / 4.0
    for _ in range(max(1, pasos)):
        d = (d + objetivo) / 2.0
    return d


def pi_densidad_primos(x: int) -> float:
    """Heurística C(x)≈π(x)/x escalada (convergencia modular del TXT)."""
    if x < 10:
        x = 97
    # conteo simple de primos ≤ x
    criba = [True] * (x + 1)
    criba[0] = criba[1] = False
    for i in range(2, int(math.isqrt(x)) + 1):
        if criba[i]:
            for j in range(i * i, x + 1, i):
                criba[j] = False
    cuenta = sum(1 for i in range(2, x + 1) if criba[i])
    c = cuenta / x
    # escalar hacia π usando factor de densidad asintótica ln(x)
    return c * math.log(max(x, 3)) * 1.15


@dataclass(frozen=True)
class MetodoPi:
    id: str
    nombre: str
    calcular: Callable[[int], float]


METODOS_PI: tuple[MetodoPi, ...] = (
    MetodoPi("inscrito", "Polígono inscrito n·sin(π/n)", pi_poligono_inscrito),
    MetodoPi("circunscrito", "Polígono circunscrito n·tan(π/n)", pi_poligono_circunscrito),
    MetodoPi("archimedes", "Arquímedes (media inscrito/circunscrito)", pi_archimedes_medio),
    MetodoPi("pitagoras", "Pitágoras Victor (4→8→16→…)", pi_pitagoras_victor),
    MetodoPi("leibniz", "Serie de Leibniz", pi_leibniz),
    MetodoPi("euler", "Serie de Euler Σ1/n²", pi_euler),
    MetodoPi("criba12", "Convergencia criba 1/2", pi_convergencia_criba),
    MetodoPi("densidad", "Densidad primos C(x)·ln(x)", pi_densidad_primos),
)


def calcular_todas_pi(seq: dict[str, list], terminos: int = 12) -> list[dict]:
    """Calcula todas las aproximaciones; cada método usa lados/índice propio."""
    resultados: list[dict] = []
    for i, metodo in enumerate(METODOS_PI):
        idx = min(i % len(seq["lados"]), len(seq["lados"]) - 1)
        n_lados = seq["lados"][idx]
        if metodo.id == "densidad":
            valor = metodo.calcular(n_lados * 37 + 100)
        else:
            valor = metodo.calcular(max(4, n_lados if metodo.id in ("inscrito", "circunscrito") else terminos))
        error = abs(valor - PI_OBJETIVO)
        resultados.append({
            "id": metodo.id,
            "nombre": metodo.nombre,
            "valor": valor,
            "error": error,
            "idx_figura": idx,
            "lados": n_lados,
            "figura": seq["figuras"][idx],
            "altura": seq["alturas"][idx],
        })
    return resultados


# ---------------------------------------------------------------------------
# Referencias elegidas por la contraseña
# ---------------------------------------------------------------------------

@dataclass
class ReferenciasCifrado:
    metodo: MetodoPi
    idx_figura: int
    figura: str
    lados: int
    altura: float
    lado: float
    perimetro: float
    pi_elegida: float
    todas_pi: list[dict]
    semillas: tuple[int, ...]


def elegir_referencias_por_clave(
    clave: str,
    semillas: tuple[int, ...] = SEMILLAS_FUNDAMENTALES,
    pasos_geom: int = 8,
) -> ReferenciasCifrado:
    """
    La contraseña guía la elección de referencias:
      - qué método π usar
      - qué figura/lados/altura de la secuencia
    """
    seq = generar_secuencias_geometricas(pasos_geom)
    digest = hashlib.sha256(clave.encode("utf-8")).digest()

    idx_metodo = digest[0] % len(METODOS_PI)
    idx_figura = digest[1] % len(seq["lados"])
    # refinar con segundo byte y longitud de clave
    idx_figura = (idx_figura + digest[2] + len(clave)) % len(seq["lados"])

    metodo = METODOS_PI[idx_metodo]
    n_lados = seq["lados"][idx_figura]
    terminos = 8 + (digest[3] % 16)

    if metodo.id == "densidad":
        pi_val = metodo.calcular(n_lados * 37 + 100 + digest[4])
    elif metodo.id in ("inscrito", "circunscrito"):
        pi_val = metodo.calcular(n_lados)
    else:
        pi_val = metodo.calcular(terminos)

    todas = calcular_todas_pi(seq, terminos)

    return ReferenciasCifrado(
        metodo=metodo,
        idx_figura=idx_figura,
        figura=seq["figuras"][idx_figura],
        lados=n_lados,
        altura=seq["alturas"][idx_figura],
        lado=seq["lado"][idx_figura],
        perimetro=seq["perimetros"][idx_figura],
        pi_elegida=pi_val,
        todas_pi=todas,
        semillas=semillas,
    )


def huella_clave(clave: str, refs: ReferenciasCifrado) -> str:
    """Huella para validar contraseña + referencias en descifrado."""
    msg = f"{clave}|{refs.metodo.id}|{refs.idx_figura}|{refs.lados}"
    return hmac.new(clave.encode(), msg.encode(), hashlib.sha256).hexdigest()[:8]


# ---------------------------------------------------------------------------
# Cuatro lenguajes Sofi + método gráfico (sin cambios sustanciales)
# ---------------------------------------------------------------------------

def _clase_l1(n: int) -> int:
    return ((n * 6 + 1) & 0xFFFFFFFF)


def _clase_l2(n: int) -> int:
    return ((n * 6 + 5) & 0xFFFFFFFF)


def _clase_l3(n: int) -> int:
    return ((n * n + 7) & 0xFFFFFFFF)


def _clase_l4(n: int) -> int:
    return ((2 * n + 1) & 0xFFFFFFFF)


_CLASES_SOFI = (_clase_l1, _clase_l2, _clase_l3, _clase_l4)


def expandir_cuatro_semillas(semillas: tuple[int, ...], clave: str) -> list[int]:
    digest = hashlib.sha256(clave.encode("utf-8")).digest()
    registros: list[int] = []
    for i, semilla in enumerate(semillas[:4]):
        mezcla = struct.unpack(">I", digest[i * 4 : i * 4 + 4])[0]
        reg = _CLASES_SOFI[i](semilla) ^ mezcla ^ (MAGICO_K3 * (i + 1))
        registros.append(reg & 0xFFFFFFFF)
    while len(registros) < 4:
        registros.append(MAGICO_K3 ^ len(registros))
    return registros


def resto_grafico(x: int, a: int) -> int:
    if a < 1 or a >= x:
        return 0
    return x % (x - a)


def secuencia_restos(x: int, limite: int = 48) -> list[int]:
    tope = min(limite, max(2, x - 1))
    return [resto_grafico(x, a) for a in range(1, tope + 1)]


def pendiente_en(i: int, restos: list[int]) -> int:
    if i < 1 or i >= len(restos):
        return 1
    dy = restos[i] - restos[i - 1]
    if dy == 0:
        return 1
    return max(1, min(12, abs(dy)))


def corte_hipotenusa(a: int, b: int, n: int, pendiente: int) -> tuple[float, float]:
    if pendiente == 0:
        return float(a), float(n - a)
    m = float(pendiente)
    cx = (n - b + m * a) / (1.0 + m)
    cy = n - cx
    return cx, cy


def energia_grafica(x: int, altura: float, lados: int) -> float:
    if x < 5:
        x = x + 37
    restos = secuencia_restos(x)
    if not restos:
        return E_VERDE * altura

    acum = 0.0
    n = len(restos) + 2
    for i in range(1, len(restos)):
        m = pendiente_en(i, restos)
        cx, cy = corte_hipotenusa(i, restos[i - 1], n, m)
        dist = math.hypot(cx / n, cy / n)
        delta = abs(restos[i] - restos[i - 1]) / max(1, x)
        acum += dist * 0.5 + delta * 0.3 + altura * 0.2

    factor_lados = math.log2(max(4, lados)) / 10.0
    return min(1.0, acum / max(1, len(restos) - 1) + factor_lados * 0.05)


def convergencia_iterativa(objetivo: float, pasos: int = ITER_CONVERGENCIA) -> float:
    d = 0.0
    for _ in range(pasos):
        d = (d + objetivo) / 2.0
    return d


def convergencia_suma_primos(x: int) -> float:
    if x < 4:
        x += 41
    ant = 1.0
    i2 = 2.0
    suma = 0.0
    lim = min(x // 2, 64)
    for i in range(2, lim + 1):
        j = x % i
        entry = (x - j) // i if i else 0
        if entry > 0 and (x % entry) == (j % entry) and x % entry == 0:
            return 0.0
        if suma <= 1.0:
            ant2 = ant * i2
            suma += (1.0 / i2) / ant + 1.0 / x
            ant = ant2
            i2 += 1.0
    return min(1.0, suma)


def criterio_raiz_cuarta(prod_a: int, prod_A: int, talla: int) -> bool:
    if prod_a <= 0:
        return False
    ratio_inv = (prod_A + prod_a) / prod_a
    umbral = (max(4, talla) / 4.0) ** 0.25
    return ratio_inv >= umbral


@dataclass
class EstadoVerde:
    registros: list[int]
    refs: ReferenciasCifrado
    fase: str = "A"

    def alternar(self) -> None:
        self.fase = "V" if self.fase == "A" else "A"

    def x_en_posicion(self, pos: int) -> int:
        r = self.registros
        base = r[0] ^ r[1] ^ r[2] ^ r[3]
        pi_i = int(self.refs.pi_elegida * 1e6) & 0xFFFFFFFF
        alt_i = int(self.refs.altura * 1e6) & 0xFFFF
        if self.fase == "A":
            return (base + pos * r[0] + (pos ** 2) * r[2] + pi_i + alt_i) & 0xFFFFFFFF
        return (base - pos * r[3] - (pos ** 2) * r[1] - pi_i) & 0xFFFFFFFF


def byte_keystream(estado: EstadoVerde, pos: int) -> int:
    refs = estado.refs
    x = estado.x_en_posicion(pos)
    eg = energia_grafica(x, refs.altura, refs.lados)
    sp = convergencia_suma_primos(x)
    pi_frac = refs.pi_elegida % 1.0
    objetivo = (eg * E_VERDE + sp * (1.0 - E_VERDE) * 0.5 + pi_frac * 0.5) % 1.0
    conv = convergencia_iterativa(objetivo)

    prod_a = estado.registros[0] * estado.registros[2] + 1
    prod_A = estado.registros[1] * estado.registros[3] + 1
    talla = max(x.bit_length(), refs.lados)
    if not criterio_raiz_cuarta(prod_a, prod_A, talla):
        conv = (conv + E_VERDE) / 2.0

    conv = conv * (1.0 - CONVERGENCIA_S) + objetivo * CONVERGENCIA_S
    kb = int(conv * 255.0) & 0xFF

    if estado.fase == "A":
        kb ^= (pos * 17 + estado.registros[pos % 4] + refs.lados) & 0xFF
    else:
        kb ^= ((255 - pos) * 13 + estado.registros[(3 - pos) % 4] + refs.idx_figura) & 0xFF

    return kb & 0xFF


# ---------------------------------------------------------------------------
# Cifrado con cabecera de referencias
# ---------------------------------------------------------------------------

def _empaquetar_cabecera(refs: ReferenciasCifrado, huella: str) -> str:
    """Referencias necesarias para recuperar (método + figura + huella)."""
    return (
        f"{HEADER_MAGIC}|{refs.metodo.id}|{refs.idx_figura}|"
        f"{refs.lados}|{refs.figura}|{huella}"
    )


def _parsear_cabecera(token: str) -> tuple[str, bytes]:
    if "|" not in token:
        raise ValueError("Token sin cabecera de referencias")
    partes = token.split("|")
    if len(partes) < 7 or partes[0] != HEADER_MAGIC:
        raise ValueError("Cabecera L4v3 inválida")
    payload_b64 = partes[-1]
    cabecera = "|".join(partes[:-1])
    return cabecera, base64.b64decode(payload_b64)


def _validar_referencias(cabecera: str, clave: str, semillas: tuple[int, ...]) -> ReferenciasCifrado:
    partes = cabecera.split("|")
    metodo_id = partes[1]
    idx_figura = int(partes[2])
    lados = int(partes[3])
    huella_esperada = partes[5]

    refs = elegir_referencias_por_clave(clave, semillas)
    if refs.metodo.id != metodo_id or refs.idx_figura != idx_figura or refs.lados != lados:
        raise ValueError(
            "Contraseña incorrecta o referencias no coinciden. "
            f"La clave selecciona: {refs.metodo.id}/{refs.figura}({refs.lados} lados)"
        )
    if huella_clave(clave, refs) != huella_esperada:
        raise ValueError("Huella de contraseña no válida")
    return refs


def procesar_bytes(
    datos: bytes,
    clave: str,
    refs: ReferenciasCifrado,
) -> bytes:
    estado = EstadoVerde(expandir_cuatro_semillas(refs.semillas, clave), refs)
    salida = bytearray(len(datos))
    for i, b in enumerate(datos):
        salida[i] = b ^ byte_keystream(estado, i)
        estado.alternar()
    return bytes(salida)


def encriptar_texto(texto: str, clave: str, semillas: tuple[int, ...]) -> str:
    refs = elegir_referencias_por_clave(clave, semillas)
    cifrado = procesar_bytes(texto.encode("utf-8"), clave, refs)
    huella = huella_clave(clave, refs)
    cabecera = _empaquetar_cabecera(refs, huella)
    return cabecera + "|" + base64.b64encode(cifrado).decode("ascii")


def desencriptar_texto(token: str, clave: str, semillas: tuple[int, ...]) -> str:
    cabecera, datos = _parsear_cabecera(token)
    refs = _validar_referencias(cabecera, clave, semillas)
    plano = procesar_bytes(datos, clave, refs)
    return plano.decode("utf-8")


# ---------------------------------------------------------------------------
# Demos
# ---------------------------------------------------------------------------

def demo_secuencias() -> dict[str, list]:
    seq = generar_secuencias_geometricas()
    print("\n=== Secuencias geométricas (r=1, inscritas) ===")
    print(f"{'Figura':<22} {'Lados':>6} {'Altura':>10} {'Lado':>10} {'π≈':>12}")
    print("-" * 64)
    for i in range(len(seq["lados"])):
        pi_a = seq["perimetros"][i] / 2.0
        print(
            f"{seq['figuras'][i]:<22} {seq['lados'][i]:>6} "
            f"{seq['alturas'][i]:>10.6f} {seq['lado'][i]:>10.6f} {pi_a:>12.8f}"
        )
    return seq


def demo_pi(seq: dict[str, list]) -> list[dict]:
    print("\n=== Aproximaciones diversas a π ===")
    print(f"π real = {PI_OBJETIVO:.12f}\n")
    todas = calcular_todas_pi(seq)
    print(f"{'Método':<28} {'Valor':>14} {'Error':>12} {'Figura':>12}")
    print("-" * 70)
    for r in todas:
        print(
            f"{r['nombre']:<28} {r['valor']:>14.10f} "
            f"{r['error']:>12.2e} {r['figura']:>12}"
        )
    return todas


def demo_grafico(n: int = 13) -> None:
    seq = demo_secuencias()
    demo_pi(seq)

    clave = "33x1-verde"
    refs = elegir_referencias_por_clave(clave)
    print(f"\n=== Referencias elegidas por clave '{clave}' ===")
    print(f"  Método π:   {refs.metodo.nombre} [{refs.metodo.id}]")
    print(f"  Figura:     {refs.figura}")
    print(f"  Lados:      {refs.lados}")
    print(f"  Altura:     {refs.altura:.8f}")
    print(f"  π elegida:  {refs.pi_elegida:.12f}")
    print(f"  Huella:     {huella_clave(clave, refs)}")

    restos = secuencia_restos(n)
    print(f"\n=== Restos gráficos (n={n}) ===")
    for a, b in enumerate(restos[:8], start=1):
        print(f"  a={a}  b={b}")

    msg = "Proyecto 33×1 — Energía Verde"
    token = encriptar_texto(msg, clave, SEMILLAS_FUNDAMENTALES)
    rec = desencriptar_texto(token, clave, SEMILLAS_FUNDAMENTALES)
    print(f"\n--- Round-trip ---")
    print(f"Original:   {msg}")
    print(f"Token:      {token[:80]}...")
    print(f"Descifrado: {rec}")
    print(f"OK:         {msg == rec}")


def modo_refs(clave: str, semillas: tuple[int, ...]) -> None:
    """Muestra qué referencias selecciona una contraseña (sin cifrar)."""
    seq = generar_secuencias_geometricas()
    demo_pi(seq)
    refs = elegir_referencias_por_clave(clave, semillas)
    print(f"\n=== Guía de referencias para clave '{clave}' ===")
    print(f"  Método:  {refs.metodo.id} — {refs.metodo.nombre}")
    print(f"  Figura:  {refs.figura} ({refs.lados} lados)")
    print(f"  Altura:  {refs.altura:.8f}")
    print(f"  π:       {refs.pi_elegida:.12f}")
    print(f"  Huella:  {huella_clave(clave, refs)}")
    print("\n  Para descifrar necesitas esta clave Y el token con estas referencias.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Libro 4 — Encriptación por convergencias y figuras π"
    )
    parser.add_argument(
        "modo",
        nargs="?",
        choices=["enc", "dec", "demo", "refs", "pi"],
        default="demo",
        help="demo|pi|refs|enc|dec",
    )
    parser.add_argument("texto", nargs="?", default="", help="Texto o token")
    parser.add_argument("--clave", "-k", default="33x1-verde", help="Contraseña guía")
    parser.add_argument("--semillas", default="2,3,5,7", help="4 semillas")
    args = parser.parse_args(list(argv) if argv is not None else None)

    partes = [int(s.strip()) for s in args.semillas.split(",")]
    if len(partes) != 4:
        print("Error: se requieren exactamente 4 semillas.", file=sys.stderr)
        return 1
    semillas = tuple(partes)

    if args.modo == "demo":
        demo_grafico()
        return 0

    if args.modo == "pi":
        demo_pi(generar_secuencias_geometricas())
        return 0

    if args.modo == "refs":
        modo_refs(args.clave, semillas)
        return 0

    if not args.texto:
        print("Error: falta texto/token.", file=sys.stderr)
        return 1

    try:
        if args.modo == "enc":
            print(encriptar_texto(args.texto, args.clave, semillas))
        else:
            print(desencriptar_texto(args.texto, args.clave, semillas))
    except (ValueError, UnicodeDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())