# -*- coding: utf-8 -*-
"""
mdcMaestro con trampa N * siguiente_primo(N)
============================================
AUTOR: Víctor Manzanares Alberola
ESCRIPTURA: Claude (Anthropic)

Fases:
  F0  Trivial (primos pequeños, cuadrado perfecto)
  F0b Criba p210 (hasta lim_criba)
  F1  Salto predictivo (pinza 4+4 + cuadrática)
  F2  TRAMPA: N * siguiente_primo(N)  (para casos donde F1 no llega)
"""

import math
import time
import sys

# ══════════════════════════════════════════════════════════════════════
#  CONSTANTES
# ══════════════════════════════════════════════════════════════════════

PRIMOS_TRIVIALES = [2, 3, 5, 7, 11, 13, 17, 19, 23]
RADI_PRED = 4
PINZA_PASOS = 4
LIM_CRIBA = 500_000
MAX_EV_SP = 200_000

# Rueda p210 (para criba F0b)
_R210 = [r for r in range(1, 211, 2) if math.gcd(r, 210) == 1]
_S210 = tuple(((_R210[(i + 1) % 48] - _R210[i]) % 210 or 210) // 2 for i in range(48))
_R210S = frozenset(_R210)


def es_L1(m):
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
#  SIGUIENTE_PRIMO (Libro 2 con 9 módulos + Karnaugh)
# ══════════════════════════════════════════════════════════════════════

def siguiente_primo(inicio: int) -> int:
    """
    Algoritmo del Libro 2 (capítulo final) para encontrar el siguiente primo
    después de 'inicio'. Funciona incluso si 'inicio' no es primo.
    """
    if inicio < 2:
        return 2

    n = inicio
    m = n + 1
    ny = n - 1
    t = ny * ny
    tt = n * n
    nt = m * m

    antp1 = False
    ant2p1 = False
    antp2 = False

    # Límite de seguridad: no más de 10000 iteraciones
    for _ in range(10000):
        # 9 módulos
        t1 = t % ny if ny != 0 else 0
        t2 = t % n  if n  != 0 else 0
        t3 = t % m  if m  != 0 else 0
        tt1 = tt % ny if ny != 0 else 0
        tt2 = tt % n  if n  != 0 else 0
        tt3 = tt % m  if m  != 0 else 0
        nt1 = nt % ny if ny != 0 else 0
        nt2 = nt % n  if n  != 0 else 0
        nt3 = nt % m  if m  != 0 else 0

        # PASO 1 (base)
        paso1_new = (t3 > 0 and nt2 == 0)
        paso1 = paso1_new or antp1 or ant2p1
        ant2p1_next = paso1_new or antp1
        antp1_next = paso1_new

        # PASO 2 (intermedio)
        paso2_new = paso1 and t2 > 0 and tt2 > 0 and t3 == 0
        paso2 = paso2_new or antp2
        antp2_next = paso2_new

        # PASO 3 (decisión final)
        if antp2 and t1 > 0 and (nt1 + nt2 == 0):
            return ny

        # Actualizar acumuladores y contadores
        t *= ny
        tt *= n
        nt *= m
        antp1, ant2p1, antp2 = antp1_next, ant2p1_next, antp2_next
        ny += 1
        n += 1
        m += 1

    # Fallback: método tradicional por si el algoritmo del libro falla
    cand = inicio + 1
    while True:
        es_primo = True
        for p in range(2, int(math.isqrt(cand)) + 1):
            if cand % p == 0:
                es_primo = False
                break
        if es_primo:
            return cand
        cand += 1


# ══════════════════════════════════════════════════════════════════════
#  CRIBA p210
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
#  SALTO PREDICTIVO (igual que en mdc_maestro)
# ══════════════════════════════════════════════════════════════════════

def _pinza_4(N, m0, avance):
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
    if len(rs) < 3:
        return [], 0
    v = [rs[i + 1] - rs[i] for i in range(len(rs) - 1)]
    a = sum(v[i + 1] - v[i] for i in range(len(v) - 1)) // max(1, len(v) - 1)
    return v, a


def _misma_parabola(a_fwd, a_bwd, D0, tol_rel=0.35):
    if a_fwd == 0 and a_bwd == 0:
        return True
    ref = max(1, D0 // 3)
    return abs(a_fwd - a_bwd) <= tol_rel * ref + 1


def _prediccion_cuadratica(r0, D0, v_r, a_r):
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
    ev = 0
    m = m_centro
    while m >= 1 and not es_L1(m):
        m -= 1
    if m < 1:
        m = sig_L1(0)
    visitados = set()
    cursor = m
    for _ in range(2 * RADI_PRED + 1):
        if cursor in visitados or cursor < m_ini or cursor > m_fi:
            cursor = sig_L1(cursor)
            continue
        visitados.add(cursor)
        D = 2 * cursor + 3
        ev += 1
        if N % D == 0 and 1 < D < N:
            return D, ev
        cursor = sig_L1(cursor)
    return None, ev


def _salto_predictivo(N, m_ini, m_fi, max_ev=MAX_EV_SP):
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
            m = min(m - 1, m_pred - RADI_PRED - 1)
            saltado = True
            break

        if not saltado:
            m = ant_L1(m)

    return None, ev


# ══════════════════════════════════════════════════════════════════════
#  MDC MAESTRO CON TRAMPA N * siguiente_primo(N)
# ══════════════════════════════════════════════════════════════════════

def mdc_maestro_con_trampa(N, lim_criba=LIM_CRIBA, max_ev_sp=MAX_EV_SP, verbose=False):
    """
    Factoriza N usando F0 (trivial), F0b (criba p210), F1 (salto predictivo)
    y, si todo falla, F2 (trampa N * siguiente_primo(N)).
    """
    t0 = time.perf_counter()

    # ── F0: TRIVIAL ──────────────────────────────────────────────
    for p in PRIMOS_TRIVIALES:
        if N % p == 0 and p < N:
            t = (time.perf_counter() - t0) * 1000
            return p, N // p, 1, "F0 trivial", t

    r = math.isqrt(N)
    if r * r == N:
        t = (time.perf_counter() - t0) * 1000
        return r, r, 1, "F0 cuadrado", t

    m_max = (math.isqrt(N) - 3) // 2
    lim_c = min(m_max, lim_criba)

    # ── F0b: CRIBA p210 ──────────────────────────────────────────
    f0, ev0 = _criba_p210(N, 1, lim_c)
    if f0:
        t = (time.perf_counter() - t0) * 1000
        return f0, N // f0, ev0, "F0b criba p210", t

    if lim_c >= m_max:
        t = (time.perf_counter() - t0) * 1000
        return None, None, ev0, "primo", t

    # ── F1: SALTO PREDICTIVO ─────────────────────────────────────
    f1, ev1 = _salto_predictivo(N, lim_c + 1, m_max, max_ev=max_ev_sp)
    ev_total = ev0 + ev1
    if f1:
        t = (time.perf_counter() - t0) * 1000
        return f1, N // f1, ev_total, "F1 salto predictivo", t

    # ── F2: TRAMPA N * siguiente_primo(N) ────────────────────────
    if verbose:
        print(f"  [F2] Intentando trampa: sp = siguiente_primo({N})")

    sp = siguiente_primo(N)
    N2 = N * sp

    # Calculamos el rango de búsqueda para N2 (de lim_c+1 a m_max2)
    m_max2 = (math.isqrt(N2) - 3) // 2
    lim_c2 = min(m_max2, lim_c)  # mantenemos el mismo límite inferior

    f2, ev2 = _salto_predictivo(N2, lim_c2 + 1, m_max2, max_ev=max_ev_sp)
    ev_total += ev2

    if f2 is not None:
        # ¡Filtro crítico! Verificamos que f2 divida a N original
        if math.gcd(N, f2) > 1:
            t = (time.perf_counter() - t0) * 1000
            # f2 puede ser p o q, o algún divisor de N
            return f2, N // f2, ev_total, "F2 trampa (N*sp)", t
        else:
            if verbose:
                print(f"  [F2] factor {f2} no divide a N (descartado)")

    t = (time.perf_counter() - t0) * 1000
    return None, None, ev_total, "no encontrado", t


# ══════════════════════════════════════════════════════════════════════
#  PRUEBA
# ══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sympy

    # Ejemplo: semiprimo con gap pequeño pero grandes (> 2e6 para que F0b no lo pille)
    p = sympy.nextprime(10**15 + 123456789)
    q = sympy.nextprime(p + 1000)
    N = p * q

    print(f"N = {N}  ({len(str(N))} dígitos)")
    print(f"p = {p}  ({len(str(p))} dígitos)")
    print(f"q = {q}  ({len(str(q))} dígitos)")
    print(f"gap = {q - p}")
    print()

    p_res, q_res, ev, fase, t = mdc_maestro_con_trampa(N, verbose=True)

    if p_res:
        print(f"✅ FACTOR ENCONTRADO: {p_res} × {q_res}")
        print(f"   fase: {fase}, evaluaciones: {ev}, tiempo: {t:.2f} ms")
    else:
        print("❌ No se encontró factor")