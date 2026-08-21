# -*- coding: utf-8 -*-
"""
mdcMaestro — núcleo unificado MDC
==================================
AUTOR: Víctor Manzanares Alberola (EPSA/UPV Alcoi)
ESCRIPTURA: Claude (Anthropic)

Fusiona en una sola función las piezas que quedaron validadas por separado
en el corpus L5:

  F0  Trivial          -> primos pequeños, N cuadrado perfecto.
  F0b Criba p210        -> militar_test.py, barrido rápido de m pequeño
                            (factores muy desbalanceados).
  F1  SALTO PREDICTIVO  -> "el truco de SP":
                            pinza doble 4+4 (avance + retroceso, filtrados
                            por L1 igual que mdc_benchmark.py) -> V/A/J
                            (velocidad/aceleración discretas de d(m)) ->
                            predicción cuadrática entera del vértice
                            (militar_test.k_sweep_predictiu) -> verificación
                            en radio ±RADI_PRED -> si falla, contraste de
                            "misma parábola" (mdc_libro5.misma_parabola)
                            antes de reducir el paso a la mitad.

No se ha tocado la fórmula del vértice (a_r·p² + 2v_r·p + 2(r0-D0) = 0);
solo se ha integrado con el filtro L1 y la doble dirección de la pinza,
que en el benchmark reducían evaluaciones respecto al barrido crudo.
"""

import math
import time

# ══════════════════════════════════════════════════════════════════════
#  CONSTANTES / RUEDA p210
# ══════════════════════════════════════════════════════════════════════

PRIMOS_TRIVIALES = [2, 3, 5, 7, 11, 13, 17, 19, 23]
RADI_PRED = 4          # radio de verificación tras la predicción
PINZA_PASOS = 4         # puntos por dirección (mínimo 3, aquí 4 -> V,A,J)

_R210 = [r for r in range(1, 211, 2) if math.gcd(r, 210) == 1]
_S210 = tuple(((_R210[(i + 1) % 48] - _R210[i]) % 210 or 210) // 2 for i in range(48))
_R210S = frozenset(_R210)


def es_L1(m):
    """m válido para D=2m+3 si D no es múltiplo de 3 (filtro barato, igual que benchmark)."""
    return m >= 1 and (2 * m + 3) % 3 != 0


def sig_L1(m):
    m += 1
    while not es_L1(m):
        m += 1
    return m


def ant_L1(m):
    m -= 1
    while m >= 1 and not es_L1(m):
        m -= 1
    return m


# ══════════════════════════════════════════════════════════════════════
#  F0b — CRIBA p210 (barrido rápido de m pequeño, desbalanceados)
# ══════════════════════════════════════════════════════════════════════

def _criba_p210(N, m_ini, m_max_c):
    if m_ini < 1:
        m_ini = 1
    m = m_ini
    while m <= m_max_c:
        if (2 * m + 3) % 210 in _R210S:
            break
        m += 1
    if m > m_max_c:
        return None, 0
    idx = _R210.index((2 * m + 3) % 210)
    ev = 0
    while m <= m_max_c:
        D = 2 * m + 3
        ev += 1
        if N % D == 0 and D < N:
            return D, ev
        m += _S210[idx % 48]
        idx += 1
    return None, ev


# ══════════════════════════════════════════════════════════════════════
#  F1 — EL TRUCO DE SP (SALTO PREDICTIVO)
#  Pinza doble 4+4 -> V/A/J -> predicción cuadrática entera -> verificación
# ══════════════════════════════════════════════════════════════════════

def _pinza_4(N, m0, avance):
    """
    Toma 4 puntos L1 en una dirección desde m0 (m0 incluido si avance,
    excluido si retroceso ya viene "gastado" desde fuera).
    Devuelve listas ms, r (restos N mod 2D), D (denominadores).
    """
    ms, rs, Ds = [], [], []
    m = m0
    for _ in range(PINZA_PASOS):
        if m < 1:
            break
        D = 2 * m + 3
        ms.append(m)
        rs.append(N % (2 * D))
        Ds.append(D)
        m = sig_L1(m) if avance else ant_L1(m)
    return ms, rs, Ds


def _vaj(rs):
    """Velocidad (1ª diferencia) y aceleración (2ª diferencia) discretas."""
    if len(rs) < 3:
        return [], 0
    v = [rs[i + 1] - rs[i] for i in range(len(rs) - 1)]
    a = sum(v[i + 1] - v[i] for i in range(len(v) - 1)) // max(1, len(v) - 1)
    return v, a


def _misma_parabola(a_fwd, a_bwd, D0, tol_rel=0.35):
    """Contraste 'misma parábola' (mdc_libro5): coherencia de curvatura
    entre las dos ramas de la pinza. Si coinciden, la predicción es fiable."""
    if a_fwd == 0 and a_bwd == 0:
        return True
    ref = max(1, D0 // 3)
    return abs(a_fwd - a_bwd) <= tol_rel * ref + 1


def _prediccion_cuadratica(r0, D0, v_r, a_r):
    """
    Resuelve a_r·p² + 2v_r·p + 2(r0-D0) = 0 (vértice/raíz entera),
    igual que militar_test.k_sweep_predictiu. Devuelve lista de candidatos
    enteros de desplazamiento p (en pasos de m) a partir del punto base.
    """
    dist = r0 - D0
    candidatos = []
    if a_r == 0:
        if v_r != 0 and dist * v_r < 0:
            p = -dist / v_r
            if p > 0:
                candidatos.append(int(round(p)))
    else:
        disc = 4 * v_r * v_r - 8 * a_r * dist
        if disc >= 0:
            sq = math.isqrt(disc)
            for sq_t in (sq, sq + 1):
                if sq_t * sq_t <= disc:
                    for sgn in (1, -1):
                        num = -2 * v_r + sgn * sq_t
                        if num * a_r > 0:
                            pf = num / (2 * a_r)
                            if pf > 0:
                                candidatos.append(int(round(pf)))
    return candidatos


def _verificar_radio(N, m_centro, m_ini, m_fi):
    """
    Comprueba m_centro y un radio ±RADI_PRED de puntos L1 alrededor,
    en AMBAS direcciones (bug anterior: solo subía con sig_L1, nunca
    bajaba, así que un objetivo por debajo de m_centro nunca se veía).
    """
    ev = 0
    visitados = set()

    # dirección ascendente: m_centro, m_centro+1L1, ...
    cursor = m_centro
    while cursor >= 1 and not es_L1(cursor):
        cursor += 1
    for _ in range(RADI_PRED + 1):
        if m_ini <= cursor <= m_fi and cursor not in visitados:
            visitados.add(cursor)
            D = 2 * cursor + 3
            ev += 1
            if N % D == 0 and 1 < D < N:
                return D, ev
        cursor = sig_L1(cursor)

    # dirección descendente: m_centro-1L1, m_centro-2L1, ...
    cursor = ant_L1(m_centro if es_L1(m_centro) else m_centro + 1)
    for _ in range(RADI_PRED):
        if cursor < 1:
            break
        if m_ini <= cursor <= m_fi and cursor not in visitados:
            visitados.add(cursor)
            D = 2 * cursor + 3
            ev += 1
            if N % D == 0 and 1 < D < N:
                return D, ev
        cursor = ant_L1(cursor)

    return None, ev


def _salto_predictivo(N, m_ini, m_fi, max_ev=200_000):
    """
    Núcleo F1: recorre [m_ini, m_fi] de arriba hacia abajo (factores
    equilibrados primero) usando pinza doble 4+4 + predicción cuadrática.
    Si la predicción no es fiable ('misma_parabola' falla) o no converge,
    retrocede con paso fino L1 en vez de saltar a ciegas.

    max_ev: presupuesto de evaluaciones. Si |p-q| es grande (RSA genuí,
    fuera del alcance del salto predictivo sin F2/bisección), abandonamos
    con honestidad en vez de arrastrarnos m_ini..m_fi entero.
    """
    if m_ini < 1:
        m_ini = 1
    if m_fi < m_ini:
        return None, 0

    ev = 0
    m = m_fi
    while m >= m_ini:
        if ev >= max_ev:
            break
        if not es_L1(m):
            m = ant_L1(m)
            if m < m_ini:
                break
            continue

        D0 = 2 * m + 3
        ev += 1
        if N % D0 == 0 and 1 < D0 < N:
            return D0, ev

        # pinza doble: retrocedemos (hacia m_ini) y avanzamos (hacia m_fi)
        ms_b, rs_b, Ds_b = _pinza_4(N, ant_L1(m), avance=False)
        ms_f, rs_f, Ds_f = _pinza_4(N, sig_L1(m), avance=True)
        ev += len(rs_b) + len(rs_f)

        if len(rs_b) < 3 and len(rs_f) < 3:
            m = ant_L1(m)
            continue

        v_b, a_b = _vaj(rs_b)
        v_f, a_f = _vaj(rs_f)
        v_r = v_b[0] if v_b else (v_f[0] if v_f else 0)
        a_r = a_b if rs_b else a_f

        fiable = _misma_parabola(a_f if rs_f else 0, a_b if rs_b else 0, D0)

        r0 = N % (2 * D0)
        candidatos = _prediccion_cuadratica(r0, D0, v_r, a_r) if fiable else []

        saltado = False
        for pc in set(candidatos):
            m_pred = m - pc
            if m_pred < m_ini - RADI_PRED:
                continue
            Dt, ev2 = _verificar_radio(N, m_pred, m_ini, m_fi)
            ev += ev2
            if Dt:
                return Dt, ev
            # el salto no dio factor pero es válido: avanzamos el cursor
            m = min(m - 1, m_pred - RADI_PRED - 1)
            saltado = True
            break

        if not saltado:
            m = ant_L1(m)

    return None, ev


# ══════════════════════════════════════════════════════════════════════
#  mdcMaestro — ENTRY POINT
# ══════════════════════════════════════════════════════════════════════

def mdcMaestro(N, lim_criba=500_000, max_ev_sp=200_000, verbose=False):
    """
    Factoriza N (devuelve un factor no trivial p tal que p*q=N, o None si
    N es primo / no se encontró dentro del alcance).

    Fases:
      F0  trivial (primos pequeños, cuadrado perfecto)
      F0b criba p210 hasta lim_criba (desbalanceados)
      F1  EL TRUCO DE SP: salto predictivo (pinza 4+4 + cuadrática)
          sobre el resto del rango [lim_criba+1, m_max] (equilibrados)

    Devuelve: (p, q, evaluaciones, fase, tiempo_ms)
    """
    t0 = time.perf_counter()

    for p in PRIMOS_TRIVIALES:
        if N % p == 0 and p < N:
            return p, N // p, 1, "F0 trivial", (time.perf_counter() - t0) * 1000

    r = math.isqrt(N)
    if r * r == N:
        return r, r, 1, "F0 cuadrado", (time.perf_counter() - t0) * 1000

    m_max = (math.isqrt(N) - 3) // 2
    lim_c = min(m_max, lim_criba)

    f0, ev0 = _criba_p210(N, 1, lim_c)
    if f0:
        return f0, N // f0, ev0, "F0b criba p210", (time.perf_counter() - t0) * 1000

    if lim_c >= m_max:
        return None, None, ev0, "primo", (time.perf_counter() - t0) * 1000

    f1, ev1 = _salto_predictivo(N, lim_c + 1, m_max, max_ev=max_ev_sp)
    t_ms = (time.perf_counter() - t0) * 1000
    if f1:
        return f1, N // f1, ev0 + ev1, "F1 salto predictivo (SP)", t_ms

    fase = "no encontrado" if N % 2 else "no encontrado"
    if ev0 + ev1 >= max_ev_sp:
        fase = "fuera d'abast (|p-q| gran, requereix F2)"
    return None, None, ev0 + ev1, fase, t_ms


# ══════════════════════════════════════════════════════════════════════
#  TESTS
# ══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sympy

    CASOS = [
        (143, [11, 13]),
        (221, [13, 17]),
        (10403, [101, 103]),
        (15251, [101, 151]),
        (999983, None),          # primo
        (104729 * 104723, [104723, 104729]),  # equilibrado
    ]

    print("=" * 70)
    print("mdcMaestro — tests fijos")
    print("=" * 70)
    for N, esperado in CASOS:
        p, q, ev, fase, t = mdcMaestro(N)
        if esperado is None:
            ok = (p is None) and sympy.isprime(N)
        else:
            ok = (p is not None) and sorted([p, q]) == sorted(esperado)
        print(f"  {'✓' if ok else '✗'}  N={N:>15,}  -> {(p, q) if p else 'primo'}  "
              f"ev={ev:<5} fase={fase:<25} {t:.3f}ms")

    print()
    print("=" * 70)
    print("mdcMaestro — RÉGIM SP: |p−q| petit, escalada 10-100 dígits")
    print("(aquí es on ha de brillar el truco de salto predictiu)")
    print("=" * 70)
    for exp in [10, 15, 20, 25, 30, 40, 50, 70, 100]:
        p = sympy.nextprime(10 ** exp + 37)
        q = sympy.nextprime(p + 1000)
        N = p * q
        p_r, q_r, ev, fase, t = mdcMaestro(N)
        ok = p_r is not None and (p_r == p or p_r == q)
        print(f"  {'✓' if ok else '✗'}  N={len(str(N)):>3}d  p={len(str(p)):>3}d  "
              f"ev={ev:<7} fase={fase:<32} {t:9.3f}ms")

    print()
    print("=" * 70)
    print("mdcMaestro — desbalanceados (p pequeño, vía F0b criba p210)")
    print("=" * 70)
    for nd in [8, 12, 16, 20]:
        p = sympy.randprime(5, 200)
        hi = 10 ** nd // p
        q = sympy.randprime(hi // 10 + 2, hi)
        N = p * q
        p_r, q_r, ev, fase, t = mdcMaestro(N)
        ok = p_r is not None and N % p_r == 0
        print(f"  {'✓' if ok else '✗'}  N={N:>22,} ({len(str(N))}d)  "
              f"ev={ev:<6} fase={fase:<32} {t:8.3f}ms")

    print()
    print("=" * 70)
    print("mdcMaestro — equilibrado GENUINO (|p−q| grande, fuera d'abast)")
    print("=" * 70)
    for nd in [10, 12, 14]:
        lo = 10 ** ((nd - 1) // 2)
        hi = 10 ** (nd // 2 + 1)
        p = sympy.randprime(lo, hi)
        q = sympy.randprime(lo, hi)
        N = p * q
        p_r, q_r, ev, fase, t = mdcMaestro(N, max_ev_sp=50_000)
        ok = p_r is not None and N % p_r == 0
        print(f"  {'✓' if ok else '⚠️ '}  N={N:>16,} ({len(str(N))}d)  "
              f"ev={ev:<6} fase={fase:<32} {t:8.3f}ms")
