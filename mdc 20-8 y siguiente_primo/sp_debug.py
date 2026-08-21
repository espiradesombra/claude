import sympy

def siguiente_primo_trace(inicio, max_iter=200):
    """Transliteración literal de libro2 (pág 73-75) con traza de depuración."""
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

    trace = []
    for it in range(max_iter):
        t1 = t % ny
        t2 = t % n
        t3 = t % m
        tt1 = tt % ny
        tt2 = tt % n
        tt3 = tt % m
        nt1 = nt % ny
        nt2 = nt % n
        nt3 = nt % m

        paso1_new = (t3 > 0 and nt2 == 0)
        paso1 = paso1_new or antp1 or ant2p1
        ant2p1_next = paso1_new or antp1
        antp1_next = paso1_new

        paso2_new = paso1 and t2 > 0 and tt2 > 0 and t3 == 0
        paso2 = paso2_new or antp2
        antp2_next = paso2_new

        paso3 = antp2 and t1 > 0 and (nt1 + nt2 == 0)

        trace.append(dict(it=it, ny=ny, n=n, m=m, t1=t1,t2=t2,t3=t3,
                           tt1=tt1,tt2=tt2,tt3=tt3,nt1=nt1,nt2=nt2,nt3=nt3,
                           paso1=paso1, paso2=paso2, paso3=paso3,
                           antp1=antp1, antp2=antp2, ant2p1=ant2p1))

        if paso3:
            return ny, n, m, it, trace

        # avanzar
        t *= ny; tt *= n; nt *= m
        antp1, ant2p1, antp2 = antp1_next, ant2p1_next, antp2_next
        ny += 1; n += 1; m += 1

    return None, None, None, None, trace

for inicio in [2,3,5,7,11]:
    esperado = sympy.nextprime(inicio)
    ny,n,m,it,trace = siguiente_primo_trace(inicio, max_iter=30)
    print(f"inicio={inicio} esperado={esperado}  -> ny={ny} n={n} m={m} it={it}")
