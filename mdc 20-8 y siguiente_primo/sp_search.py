import sympy, itertools

def run(inicio, grp_paso1, update_before_mod, max_iter=200):
    n = inicio
    m = n + 1
    ny = n - 1
    t = ny*ny
    tt = n*n
    nt = m*m
    n += 1; m += 1; ny += 1

    antp1 = False
    ant2p1 = False
    antp2 = False

    for it in range(max_iter):
        if update_before_mod:
            t2_,tt2_,nt2_ = t*ny, tt*n, nt*m
        else:
            t2_,tt2_,nt2_ = t, tt, nt

        t1 = t2_ % ny
        t2 = t2_ % n
        t3 = t2_ % m
        tt1 = tt2_ % ny
        tt2 = tt2_ % n
        tt3 = tt2_ % m
        nt1 = nt2_ % ny
        nt2 = nt2_ % n
        nt3 = nt2_ % m

        base1 = (t3 > 0 and nt2 == 0)
        if grp_paso1 == 'A':   # (A and B) or C or D
            paso1 = base1 or antp1 or ant2p1
        else:                  # A and (B or C or D)
            paso1 = (t3 > 0) and (nt2 == 0 or antp1 or ant2p1)

        ant2p1_next = base1 or antp1
        antp1_next = base1

        paso2_new = paso1 and t2 > 0 and tt2 > 0 and t3 == 0
        antp2_next = paso2_new

        paso3 = antp2 and t1 > 0 and (nt1 + nt2 == 0)

        if paso3:
            return ny

        if not update_before_mod:
            t *= ny; tt *= n; nt *= m
        antp1, ant2p1, antp2 = antp1_next, ant2p1_next, antp2_next
        ny += 1; n += 1; m += 1

    return None

configs = list(itertools.product(['A','B'], [False, True]))
primos = list(sympy.primerange(2, 300))

for grp, upd in configs:
    fails = 0
    for inicio in primos:
        esperado = sympy.nextprime(inicio)
        got = run(inicio, grp, upd, max_iter=100)
        if got != esperado:
            fails += 1
    print(f"grp_paso1={grp}  update_before_mod={upd}  -> fallos {fails}/{len(primos)}")
