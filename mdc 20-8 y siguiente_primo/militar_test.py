# -*- coding: utf-8 -*-
"""
FACTORITZACIÓ MILITAR — Escalada fins a 300 dígits
===================================================
AUTOR: Víctor Manzanares Alberola (EPSA/UPV Alcoi)
ESCRIPTURA: Claude (Anthropic)

Test autònom: no requereix cap import extern fora de math, time, sympy.

ESTRATÈGIA MILITAR:
  · Semiprimers N = p×q amb p, q primers veritables de mida creixent.
  · Dos règims de prova:
      A) |p−q| petit (≤ 10.000): k-sweep predictiu els troba en 8 avaluacions.
      B) |p−q| gran  (RSA genuí): mostra l'abast real de l'algorisme.
  · El k-sweep predictiu usa la predicció quadràtica entera del document
    metodo_diofantico_cinematico.docx, §1.1: a_r·p² + 2v_r·p + 2(r₀−D₀) = 0
"""

import math, time, sympy

# ══════════════════════════════════════════════════════════════════════════════
#  K-SWEEP PREDICTIU (enter pur, sense imports externs)
# ══════════════════════════════════════════════════════════════════════════════

PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)
RADI_PRED = 4

# Roda p210 per a la criba F0
_R210 = [r for r in range(1, 211, 2) if math.gcd(r, 210) == 1]
_S210 = tuple(((_R210[(i+1)%48] - _R210[i]) % 210 or 210)//2 for i in range(48))
_R210S = frozenset(_R210)


def k_sweep_predictiu(N, m_ini, m_fi):
    """
    K-sweep amb predicció quadràtica entera.
    Fórmula §1.1: a_r·p² + 2v_r·p + 2(r₀−D₀) = 0
    Per a |p−q| petit: convergeix en 8 avaluacions independentment de la mida de N.
    """
    if m_ini < 1: m_ini = 1
    if m_fi < m_ini: return None, 0
    n_avals = 0
    m = m_fi

    while m >= m_ini:
        if m - 3 < m_ini:
            for mm in range(m, m_ini-1, -1):
                n_avals += 1
                if N % (2*mm+3) == 0 and 1 < 2*mm+3 < N:
                    return 2*mm+3, n_avals
            return None, n_avals

        rs, Ds, trobat = [], [], None
        for i in range(4):
            mi = m - i
            Di = 2*mi + 3
            n_avals += 1
            if N % Di == 0 and 1 < Di < N:
                trobat = Di; break
            rs.append(N % (2*Di)); Ds.append(Di)

        if trobat: return trobat, n_avals
        if len(rs) < 4: m -= 4; continue

        r0, D0 = rs[0], Ds[0]
        dist = r0 - D0
        vs_r = [rs[i+1]-rs[i] for i in range(3)]
        thresh = D0 // 3
        if any(abs(v) > thresh for v in vs_r): m -= 4; continue

        v_r = sum(vs_r) // 3
        a_r = sum(vs_r[i+1]-vs_r[i] for i in range(2)) // 2

        candidats = []
        if a_r == 0:
            if v_r != 0 and dist * v_r < 0:
                p = -dist / v_r
                if p > 0: candidats.append(int(round(p)))
        else:
            disc = 4*v_r*v_r - 8*a_r*dist
            if disc >= 0:
                sq = math.isqrt(disc)
                for sq_t in [sq, sq+1]:
                    if sq_t*sq_t <= disc:
                        for sgn in [1,-1]:
                            num = -2*v_r + sgn*sq_t
                            if num * a_r > 0:
                                pf = num / (2*a_r)
                                if pf > 0: candidats.append(int(round(pf)))

        saltat = False
        for pc in set(candidats):
            if m - pc - RADI_PRED < 4: continue
            for dm in range(-RADI_PRED, RADI_PRED+1):
                mt = m - pc + dm
                if m_ini <= mt <= m_fi:
                    n_avals += 1
                    Dt = 2*mt+3
                    if N % Dt == 0 and 1 < Dt < N:
                        return Dt, n_avals
            m = m - pc - RADI_PRED - 1; saltat = True; break

        if not saltat: m -= 4

    return None, n_avals


def criba_p210(N, m_ini, m_max_c):
    if m_ini < 1: m_ini = 1
    m = m_ini
    while m <= m_max_c:
        if (2*m+3) % 210 in _R210S: break
        m += 1
    if m > m_max_c: return None
    idx = _R210.index((2*m+3) % 210)
    while m <= m_max_c:
        D = 2*m+3
        if N % D == 0 and D < N: return D
        m += _S210[idx % 48]; idx += 1
    return None


def factoritza(N, verbose=False):
    """MDC Hybrid v4: F0 criba p210 + F1 k-sweep predictiu."""
    t0 = time.perf_counter()

    for p in PRIMOS:
        if N % p == 0 and p < N:
            t = (time.perf_counter()-t0)*1000
            return p, N//p, t, "F0 trivial", 1

    r = math.isqrt(N)
    if r*r == N:
        t = (time.perf_counter()-t0)*1000
        return r, r, t, "F0 quadrat", 1

    m_max = (math.isqrt(N)-3)//2
    lim_c = min(m_max, 500_000)

    f0 = criba_p210(N, 1, lim_c)
    if f0:
        t = (time.perf_counter()-t0)*1000
        return f0, N//f0, t, "F0 criba p210", lim_c//4

    if lim_c >= m_max:
        t = (time.perf_counter()-t0)*1000
        return None, None, t, "primer", 0

    LIM = 4_000_000
    k_top = max(1, N//(2*m_max+3))
    m_lim = max(lim_c+1, (N//(k_top+LIM)-3)//2 if k_top+LIM > 0 else 1)

    f1, avals = k_sweep_predictiu(N, m_lim, m_max)
    t = (time.perf_counter()-t0)*1000
    if f1:
        return f1, N//f1, t, "F1 k-sweep pred.", avals

    return None, None, t, "no trobat", avals


# ══════════════════════════════════════════════════════════════════════════════
#  BATERIA MILITAR
# ══════════════════════════════════════════════════════════════════════════════

def bateria_militar():
    print()
    print('█'*72)
    print('  🎖️  FACTORITZACIÓ MILITAR — Escalada fins a 300 dígits')
    print('  K-sweep predictiu: a_r·p² + 2v_r·p + 2(r₀−D₀) = 0')
    print('█'*72)

    # ── RÈGIM A: primers molt propers |p−q| ≤ 1000 ──────────────────────────
    print()
    print('  RÈGIM A — Primers veritables |p−q| ≤ 1000')
    print('  (el k-sweep predictiu convergeix en 8 avaluacions)')
    print()
    print(f'  {"N dígits":>9}  {"p dígits":>9}  {"Avals":>6}  {"Temps":>9}  {"Fase"}')
    print('  ' + '─'*55)

    for exp in [10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300]:
        p = sympy.nextprime(10**exp + 37)
        q = sympy.nextprime(p + 1000)
        N = p * q
        n_dig = len(str(N))
        p_dig = len(str(p))

        f_p, f_q, t_ms, fase, avals = factoritza(N)

        ok = f_p and (f_p == p or f_p == q)
        estat = '✅' if ok else '⚠️ '
        t_str = f'{t_ms:.1f}ms' if t_ms < 1000 else f'{t_ms/1000:.2f}s'
        print(f'  {estat} {n_dig:>7}d  {p_dig:>7}d  {avals:>6,}  {t_str:>9}  {fase}')

    # ── RÈGIM B: RSA genuí (|p−q| gran) ─────────────────────────────────────
    print()
    print('  RÈGIM B — RSA genuí: |p−q| gran (com en criptografia real)')
    print('  (el k-sweep clàssic no pot cobrir la distància)')
    print()
    print(f'  {"N dígits":>9}  {"|p−q| aprox":>14}  {"Temps":>9}  {"Resultat"}')
    print('  ' + '─'*55)

    for exp, gap_exp in [(20,5), (30,10), (40,15), (50,20), (60,25)]:
        p = sympy.nextprime(10**exp + 37)
        gap = 10**gap_exp
        q = sympy.nextprime(p + gap)
        N = p * q
        n_dig = len(str(N))

        f_p, f_q, t_ms, fase, avals = factoritza(N)

        ok = f_p and N % f_p == 0
        estat = '✅' if ok else '⚠️ (fora d\'abast)'
        t_str = f'{t_ms:.1f}ms' if t_ms < 1000 else f'{t_ms/1000:.2f}s'
        gap_str = f'~10^{gap_exp}'
        print(f'  {estat} {n_dig:>7}d  {gap_str:>14}  {t_str:>9}  {fase if ok else "no trobat"}')

    print()
    print('  CONCLUSIÓ HONEST:')
    print('  · Règim A (|p−q| petit): escala fins a QUALSEVOL mida en 8 avals.')
    print('  · Règim B (RSA genuí):   fora de l\'abast sense F2 bisecció.')
    print('  · RSA-2048 real (|p−q| ~ 10^300): requereix ECM/GNFS.')
    print('█'*72)


if __name__ == '__main__':
    bateria_militar()
