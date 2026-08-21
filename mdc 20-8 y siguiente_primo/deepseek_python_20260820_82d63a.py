import math
import time

def siguiente_primo(inicio: int) -> int:
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

    MAX_ITER = 50000

    for _ in range(MAX_ITER):
        # ── 9 MÓDULOS DE FASE (TODOS DEBEN CALCULARSE) ──
        t1 = t % ny if ny != 0 else 0
        t2 = t % n  if n  != 0 else 0
        t3 = t % m  if m  != 0 else 0

        tt1 = tt % ny if ny != 0 else 0
        tt2 = tt % n  if n  != 0 else 0
        tt3 = tt % m  if m  != 0 else 0

        nt1 = nt % ny if ny != 0 else 0
        nt2 = nt % n  if n  != 0 else 0
        nt3 = nt % m  if m  != 0 else 0

        # ── PASO 1 ──────────────────────────────────────
        paso1_new = (t3 > 0 and nt2 == 0)
        paso1 = paso1_new or antp1 or ant2p1

        # ── PASO 2 ──────────────────────────────────────
        paso2_new = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)

        # ── PASO 3 ──────────────────────────────────────
        if antp2 and (t1 > 0) and (nt1 + nt2 == 0):
            return ny

        # ── ACTUALIZAR ACUMULADORES (¡IMPRESCINDIBLE!) ──
        ant2p1_next = paso1_new or antp1
        antp1_next = paso1_new
        antp2_next = paso2_new

        t *= ny
        tt *= n
        nt *= m

        antp1, ant2p1, antp2 = antp1_next, ant2p1_next, antp2_next
        ny += 1
        n += 1
        m += 1

    # ── FALLBACK (casi nunca se ejecuta) ──────────────
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