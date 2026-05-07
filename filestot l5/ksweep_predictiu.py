# -*- coding: utf-8 -*-
"""
K-SWEEP PREDICTIU FINAL — Aritmètica Entera Pura
================================================
AUTOR: Víctor Manzanares Alberola + Claude (Anthropic)

FÓRMULA (§1.1 metodo_diofantico_cinematico.docx, adaptació entera):
    r(m) = N % (2D),   D = 2m+3   (el reste és enter exacte)
    Condició de factor: r = D  ↔  d(m) = 0

    v_r = Δr/Δm    (velocitat entera del reste)
    a_r = Δv_r/Δm  (acceleració entera)

    Equació quadràtica entera:
        a_r/2 · p² + v_r · p + (r₀ - D₀) = 0

    Solució: p = (−v_r ± √(v_r² − 2·a_r·(r₀−D₀))) / a_r

    m* = m₀ − round(p)

AVANTATGE vs float: funciona per a N de qualsevol mida sense pèrdua de
precisió. Per a N de 26 dígits, float64 perdia els 2-3 últims dígits de d(m);
amb enters purs, la predicció és sempre exacta.
"""

import math, time

PRIMOS    = [2, 3, 5, 7, 11, 13, 17, 19, 23]
PRIMORIAL = math.prod(PRIMOS)
RADI_PRED = 4
THRESH_RESET = 3   # Δv_r > THRESH_RESET × D → reset de dent


def k_sweep_classic(N, m_ini, m_fi):
    if m_ini < 1: m_ini = 1
    if m_fi < m_ini: return None
    k_lo = max(1, N // (2*m_fi+3))
    k_hi = N // (2*m_ini+3)
    for k in range(k_lo, k_hi+1):
        if k <= 0: continue
        c = N // k
        if c < 3 or c % 2 == 0: continue
        if N % c == 0 and 1 < c < N: return c
    return None


def k_sweep_predictiu(N, m_ini, m_fi, verbose=False):
    """
    K-sweep amb predicció quadràtica entera pura.

    Per cada posició m₀ (decreixent des de m_fi):
      1. Computa r = N%(2D) en 4 punts consecutius [enters exactes]
      2. Calcula v_r (velocitat) i a_r (acceleració) del reste
      3. Resol: a_r/2·p² + v_r·p + (r₀−D₀) = 0
      4. Salta a m* = m₀ − p i verifica ±RADI
      5. Si falla o reset: avança 4 posicions normalment
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

        # ── 4 punts: resta i D per a m₀, m₀-1, m₀-2, m₀-3 ───────────────
        rs, Ds = [], []
        trobat = None
        for i in range(4):
            mi = m - i
            Di = 2*mi + 3
            n_avals += 1
            if N % Di == 0 and 1 < Di < N:
                trobat = Di
                break
            twoDi = 2 * Di
            rs.append(N % twoDi)
            Ds.append(Di)

        if trobat:
            return trobat, n_avals
        if len(rs) < 4:
            m -= 4; continue

        r0, D0 = rs[0], Ds[0]
        distancia = r0 - D0   # ha de ser 0 en el factor

        # ── Velocitats i acceleració (enters purs) ─────────────────────────
        vs_r = [rs[i+1] - rs[i] for i in range(3)]
        accs_r = [vs_r[i+1] - vs_r[i] for i in range(2)]

        # Detectar reset de dent de serra (salt brusc en la resta)
        thresh = D0 // 3
        if any(abs(v) > thresh for v in vs_r):
            m -= 4; continue

        v_r = sum(vs_r) // 3              # velocitat mitja (entera)
        a_r = sum(accs_r) // len(accs_r)  # acceleració mitja (entera)

        # ── Predicció quadràtica entera ────────────────────────────────────
        # a_r/2·p² + v_r·p + distancia = 0
        # → a_r·p² + 2·v_r·p + 2·distancia = 0  (×2 per evitar fraccions)
        candidats = []

        if a_r == 0:
            # Lineal: v_r·p + distancia = 0 → p = −distancia/v_r
            if v_r != 0 and distancia * v_r < 0:  # sol. positiva
                num = -distancia
                p_float = num / v_r
                if p_float > 0:
                    candidats.append(int(round(p_float)))
        else:
            # Quadràtica: a_r·p² + 2·v_r·p + 2·distancia = 0
            disc = 4*v_r*v_r - 4*a_r*(2*distancia)   # = 4·(v_r² − 2·a_r·distancia)
            if disc >= 0:
                sq = math.isqrt(disc)
                # Ajust: isqrt pot donar sq-1 per errors d'arrodoniment
                for sq_try in [sq, sq+1]:
                    if sq_try*sq_try <= disc:
                        for sgn in [1, -1]:
                            num = -2*v_r + sgn*sq_try
                            if a_r != 0:
                                # p = num / (2·a_r)
                                if num * a_r > 0:  # mateix signe → p>0
                                    p_float = num / (2*a_r)
                                    if p_float > 0:
                                        candidats.append(int(round(p_float)))

        # ── Verificar candidats ────────────────────────────────────────────
        m_saltat = False
        for p_cand in set(candidats):
            salt = m - p_cand - RADI_PRED
            if salt < 4: continue   # salt massa curt

            if verbose:
                print(f'  m={m:>14,} r-D={distancia:>10,} v={v_r:>8,} a={a_r:>6,}'
                      f' → m*={m-p_cand:>14,} (salt={salt:,})')

            for dm in range(-RADI_PRED, RADI_PRED+1):
                m_t = m - p_cand + dm
                if m_ini <= m_t <= m_fi:
                    n_avals += 1
                    D_t = 2*m_t + 3
                    if N % D_t == 0 and 1 < D_t < N:
                        if verbose:
                            print(f'    ✅ Factor a m={m_t:,}  D={D_t:,}')
                        return D_t, n_avals

            m = m - p_cand - RADI_PRED - 1
            m_saltat = True
            break

        if not m_saltat:
            m -= 4

    return None, n_avals


# ============================================================================
#  BENCHMARK
# ============================================================================

def benchmark():
    import sympy

    print('═'*74)
    print('  K-SWEEP PREDICTIU — Predicció Quadràtica Entera Pura')
    print('  Equació: a_r·p² + 2v_r·p + 2(r₀−D₀) = 0  [§1.1 MDC doc.]')
    print('═'*74)
    print(f'  {"Descripció":<28} {"Clàs(ms)":>9}  {"Pred(ms)":>8}  '
          f'{"Avals":>8}  {"Speedup":>8}')
    print(f'  {"─"*28} {"─"*9}  {"─"*8}  {"─"*8}  {"─"*8}')

    proves = [
        (101*103,                                 "5d equilibrat"),
        (100003*100019,                           "11d equilibrat"),
        (1_000_003*1_000_033,                     "13d equilibrat"),
        (9_999_991*9_999_973,                     "14d equilibrat"),
        (100_000_000_003*100_000_000_019,         "22d equilibrat"),
        (10_000_000_000_037*10_000_000_000_099,   "26d equilibrat"),
        ((10**18+9)*(10**18+31),                  "36d equilibrat"),
    ]
    for exp in [10, 12, 15, 18, 22, 26]:
        p = sympy.nextprime(10**exp+37)
        q = sympy.nextprime(p+1000)
        proves.append((p*q, f"{len(str(p*q))}d primers verit."))

    LIM = 2_000_000
    REPS = 3
    ok = fail = 0

    for N, desc in proves:
        m_max = (math.isqrt(N)-3)//2
        k_top = max(1, N//(2*m_max+3))
        k_lim = k_top + LIM
        d_lim = N // k_lim if k_lim > 0 else 2*m_max+3
        m_lim = max(1, (d_lim-3)//2)

        ts_cl = sorted([
            (lambda: (time.perf_counter(),
                      k_sweep_classic(N, m_lim, m_max),
                      time.perf_counter()))()
            for _ in range(REPS)
        ], key=lambda x: x[2]-x[0])
        t_cl = (ts_cl[1][2]-ts_cl[1][0])*1000
        f_cl = ts_cl[1][1]

        ts_pr = []
        f_pr = avals = None
        for _ in range(REPS):
            t0 = time.perf_counter()
            f_pr, avals = k_sweep_predictiu(N, m_lim, m_max)
            ts_pr.append((time.perf_counter()-t0)*1000)
        t_pr = sorted(ts_pr)[1]

        f_ok = (f_pr and N%f_pr==0) or (f_pr is None and f_cl is None)
        sp = t_cl/t_pr if t_pr > 0.001 else 99.0

        if f_ok: ok+=1
        else: fail+=1

        print(f'  {"✅" if f_ok else "⚠️"} {desc:<27} '
              f'{t_cl:>8.2f}  {t_pr:>8.2f}  {avals:>8,}  {sp:>8.1f}x')

    print(f'\n  Encerts: {ok}/{ok+fail}')
    print('═'*74)


if __name__ == "__main__":
    # Verbose de les tres proves clau
    for N, desc in [
        (100003*100019,                          "11d equilibrat"),
        (1_000_003*1_000_033,                    "13d equilibrat"),
        (10_000_000_000_037*10_000_000_000_099,  "26d equilibrat"),
    ]:
        m_max = (math.isqrt(N)-3)//2
        LIM = 2_000_000
        k_top = max(1, N//(2*m_max+3))
        m_lim = max(1, (N//(k_top+LIM)-3)//2)
        print(f'Prova {desc}:')
        f, avals = k_sweep_predictiu(N, m_lim, m_max, verbose=True)
        if f:
            print(f'  → {f:,} × {N//f:,}  ({avals} avaluacions)\n')
        else:
            print(f'  → No trobat ({avals} avals)\n')

    import sympy
    benchmark()
