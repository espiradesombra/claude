
import math

def is_candidate(v):
    return v % 6 == 1 or v % 6 == 5

def value(i):
    return 2 * i + 3

def criba_modular(L):
    Nmax = L * L
    p = [0] * (Nmax + 1)
    for v in range(Nmax + 1):
        if is_candidate(v):
            p[v] = 1
    for seed in [2, 3, 5, 7]:
        if seed <= Nmax:
            p[seed] = 1

    for n in range((L - 3) // 2 + 1):
        pn = value(n)
        if pn > Nmax or p[pn] != 1:
            continue
        m = n
        siguiente_m_multiplo3 = 0
        while True:
            qm = value(m)
            if qm > Nmax // pn:
                break
            if m == siguiente_m_multiplo3:
                siguiente_m_multiplo3 += 3
                m += 1
                continue
            if not is_candidate(qm):
                m += 1
                continue
            if pn <= qm:
                comp = pn * qm
                salto1 = 2 * pn
                salto2 = 4 * pn
                j = comp
                # Fase correcta: p≡5 (mod 6) → 2p primer; p≡1 → 4p primer
                usar_salto1 = (pn % 6 == 5)
                while j <= Nmax:
                    if j <= Nmax and p[j] == 1:
                        p[j] = 3
                    j += salto1 if usar_salto1 else salto2
                    usar_salto1 = not usar_salto1
            m += 1
    # 2 i 3 són prims; 0,1 i compostos marcat 3
    if Nmax >= 2:
        p[2] = 1
    if Nmax >= 3:
        p[3] = 1
    # treure 1 com a candidat residual
    if Nmax >= 1:
        p[1] = 0
    return p
