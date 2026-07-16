# -*- coding: utf-8 -*-
"""
BENCHMARK DE FRONTERA G*(D) — MDC Híbrid
=========================================
AUTOR: Víctor Manzanares Alberola + Claude (Anthropic)

Objectiu: caracteritzar empíricament la frontera G*(D) del mètode,
és a dir, per a cada mida D de N, quin és el gap màxim G=log10(|p-q|)
on el mètode encerta amb probabilitat ≥ 50%.

Variables:
  D = log10(N) ≈ nombre de dígits de N
  G = log10(|p-q|)  — el gap entre factors
  Resultat: P(èxit) vs (D, G) → heatmap + corba G*(D)

Mètode: multiprocessing per a executar ~1000 casos en minuts.
"""

import math, time, sympy, random, csv, os, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

# ──────────────────────────────────────────────────────────────────────────────
#  MOTOR MDC (autònom, sense imports externs)
# ──────────────────────────────────────────────────────────────────────────────

_R210  = [r for r in range(1,211,2) if math.gcd(r,210)==1]
_S210  = tuple(((_R210[(i+1)%48]-_R210[i])%210 or 210)//2 for i in range(48))
_R210S = frozenset(_R210)
_RADI  = 4

def _factoritza_mdc(N, timeout_ms=2000):
    """
    MDC Híbrid v4 complet: F0 criba p210 + F1 k-sweep predictiu.
    Retorna (factor, avals, temps_ms) o (None, avals, temps_ms).
    """
    t0 = time.perf_counter()

    def elapsed():
        return (time.perf_counter()-t0)*1000

    # Trivials
    for pp in [2,3,5,7,11,13,17,19,23]:
        if N%pp==0 and pp<N: return pp, 1, elapsed()

    r = math.isqrt(N)
    if r*r==N: return r, 1, elapsed()

    mm  = (math.isqrt(N)-3)//2
    lc  = min(mm, 500_000)

    # F0 criba p210
    m = 1
    while m <= lc:
        if (2*m+3)%210 in _R210S: break
        m += 1
    if m <= lc:
        idx = _R210.index((2*m+3)%210)
        while m <= lc:
            if elapsed() > timeout_ms: return None, -1, elapsed()
            D = 2*m+3
            if N%D==0 and D<N: return D, lc//4, elapsed()
            m += _S210[idx%48]; idx += 1

    if lc >= mm: return None, 0, elapsed()

    # F1 k-sweep predictiu (rang complet)
    n = 0; m = mm
    while m >= lc:
        if elapsed() > timeout_ms: return None, n, elapsed()
        if m-3 < lc:
            for x in range(m, lc-1, -1):
                n += 1
                if N%(2*x+3)==0 and 1<2*x+3<N:
                    return 2*x+3, n, elapsed()
            return None, n, elapsed()
        rs, Ds, t2 = [], [], None
        for i in range(4):
            mi=m-i; Di=2*mi+3; n+=1
            if N%Di==0 and 1<Di<N: t2=Di; break
            rs.append(N%(2*Di)); Ds.append(Di)
        if t2: return t2, n, elapsed()
        if len(rs)<4: m-=4; continue
        r0,D0=rs[0],Ds[0]; dist=r0-D0
        vs=[rs[i+1]-rs[i] for i in range(3)]
        if any(abs(v)>D0//3 for v in vs): m-=4; continue
        vr=sum(vs)//3; ar=sum(vs[i+1]-vs[i] for i in range(2))//2
        cc=[]
        if ar==0:
            if vr!=0 and dist*vr<0:
                pi=(-dist)//vr
                for pt in [pi,pi+1,pi-1]:
                    if pt>0: cc.append(pt)
        else:
            d=4*vr*vr-8*ar*dist
            if d>=0:
                sq=math.isqrt(d)
                for st in [sq,sq+1]:
                    if st*st<=d:
                        for s in [1,-1]:
                            num=-2*vr+s*st
                            if num*ar>0:
                                pi=num//(2*ar)
                                for pt in [pi,pi+1,pi-1]:
                                    if pt>0: cc.append(pt)
        sl=False
        for pc in set(cc):
            if m-pc-_RADI<4: continue
            for dm in range(-_RADI,_RADI+1):
                mt=m-pc+dm
                if lc<=mt<=mm:
                    n+=1
                    if N%(2*mt+3)==0 and 1<2*mt+3<N:
                        return 2*mt+3, n, elapsed()
            m=m-pc-_RADI-1; sl=True; break
        if not sl: m-=4
    return None, n, elapsed()


# ──────────────────────────────────────────────────────────────────────────────
#  WORKER
# ──────────────────────────────────────────────────────────────────────────────

def worker(task):
    """
    Executa un cas (D, G, seed) i retorna el resultat.
    D = dígits de N, G = log10(|p-q|), seed = reproductibilitat.
    """
    D, G, seed = task
    rng = random.Random(seed)

    try:
        offset = rng.randint(0, 10**6)
        p = sympy.nextprime(10**D + offset)
        gap = int(10**G)
        jitter = rng.randint(-gap//10, gap//10) if gap > 10 else 0
        q = sympy.nextprime(p + gap + jitter)
        N = p * q
        n_dig = len(str(N))

        f, avals, t_ms = _factoritza_mdc(N, timeout_ms=3000)

        ok = 1 if (f and N%f==0 and 1<f<N) else 0
        gap_real = abs(p-q)
        G_real = math.log10(gap_real) if gap_real > 0 else 0

        return {
            'D': D, 'G': G, 'seed': seed,
            'N_digits': n_dig, 'p_digits': len(str(p)),
            'success': ok, 'time_ms': round(t_ms,2),
            'avals': avals,
            'G_real': round(G_real, 2),
            'gap_real': gap_real,
        }
    except Exception as e:
        return {
            'D': D, 'G': G, 'seed': seed,
            'N_digits': -1, 'p_digits': -1,
            'success': 0, 'time_ms': -1,
            'avals': -1, 'G_real': -1, 'gap_real': -1,
        }


# ──────────────────────────────────────────────────────────────────────────────
#  CONFIGURACIÓ DE L'EXPERIMENT
# ──────────────────────────────────────────────────────────────────────────────

# D: dígits de N (log10 aproximat)
D_VALUES = [20, 30, 40, 50, 60, 80, 100]

# G: log10(|p-q|)
G_VALUES = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20]

# Mostres per cada (D, G)
N_MOSTRES = 10  # ← reduir per a test ràpid, augmentar per a paper (30-50)

TIMEOUT_S  = 3.0
CSV_OUTPUT = 'benchmark_frontera.csv'


# ──────────────────────────────────────────────────────────────────────────────
#  RUNNER
# ──────────────────────────────────────────────────────────────────────────────

def run_benchmark(n_workers=None, quick=False):
    if quick:
        # Versió ràpida per a testejar
        D_vals = [20, 30, 40, 50, 60]
        G_vals = [0, 2, 4, 6, 8, 10, 12, 15]
        n_m    = 5
    else:
        D_vals = D_VALUES
        G_vals = G_VALUES
        n_m    = N_MOSTRES

    tasks = []
    for D in D_vals:
        for G in G_vals:
            if G >= D:          # G >= D → factor > N → impossible
                continue
            for seed in range(n_m):
                tasks.append((D, G, seed * 1000 + D * 100 + G))

    # Barrejar per evitar sesgo temporal
    random.shuffle(tasks)

    n_tot = len(tasks)
    n_workers = n_workers or max(1, (os.cpu_count() or 4) - 1)

    print(f'  {n_tot} casos × {n_workers} workers')
    print(f'  D={D_vals}')
    print(f'  G={G_vals}')
    print()

    resultats = []
    n_done = 0
    t_inici = time.perf_counter()

    # Escriptura CSV
    with open(CSV_OUTPUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'D','G','seed','N_digits','p_digits',
            'success','time_ms','avals','G_real','gap_real'
        ])
        writer.writeheader()

        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            chunk = 5
            futures = {
                ex.submit(worker, t): t
                for t in tasks
            }
            for fut in as_completed(futures):
                r = fut.result()
                writer.writerow(r)
                f.flush()
                resultats.append(r)
                n_done += 1
                if n_done % 20 == 0 or n_done == n_tot:
                    elapsed = time.perf_counter() - t_inici
                    eta = elapsed / n_done * (n_tot - n_done)
                    ok = sum(1 for x in resultats if x['success'])
                    print(f'  [{n_done:>4}/{n_tot}]  '
                          f'{ok/n_done*100:.0f}% èxit  '
                          f'ETA {eta:.0f}s', end='\r', flush=True)

    print()
    return resultats


# ──────────────────────────────────────────────────────────────────────────────
#  ANÀLISI: FRONTERA G*(D)
# ──────────────────────────────────────────────────────────────────────────────

def analitzar_frontera(resultats):
    """
    Per a cada (D, G): calcula la probabilitat d'èxit.
    Troba G*(D) = màxim G on P(èxit) ≥ 0.5.
    """
    from collections import defaultdict
    data = defaultdict(list)
    for r in resultats:
        data[(r['D'], r['G'])].append(r['success'])

    print()
    print('FRONTERA G*(D) — Taxa d\'èxit per (D, G):')
    print()

    # Capçalera
    D_vals = sorted(set(r['D'] for r in resultats))
    G_vals = sorted(set(r['G'] for r in resultats))

    # Taula
    header = f'  {"D\\G":>6}' + ''.join(f'{g:>7}' for g in G_vals)
    print(header)
    print('  ' + '─' * (6 + 7*len(G_vals)))

    frontera = {}  # D → G* màxim on P≥0.5
    for D in D_vals:
        row = f'  {D:>6}'
        g_star = 0
        for G in G_vals:
            casos = data.get((D,G), [])
            if not casos:
                row += f'  {"─":>5}'
                continue
            p_exit = sum(casos)/len(casos)
            simbol = f'{p_exit:.0%}'
            row += f'{simbol:>7}'
            if p_exit >= 0.5:
                g_star = G
        print(row)
        frontera[D] = g_star

    print()
    print('FRONTERA G*(D) (màxim G amb P≥50%):')
    for D, G_star in frontera.items():
        print(f'  D={D:>4}d → G*={G_star}  (|p-q| ≤ 10^{G_star})')

    # Ajust lineal log(G*) ~ a*log(D) + b
    if len(frontera) >= 3:
        xs = [math.log10(D) for D in frontera.keys() if frontera[D] > 0]
        ys = [frontera[D]   for D in frontera.keys() if frontera[D] > 0]
        if len(xs) >= 2:
            n = len(xs)
            sx  = sum(xs);  sy  = sum(ys)
            sxy = sum(x*y for x,y in zip(xs,ys))
            sx2 = sum(x*x for x in xs)
            a = (n*sxy - sx*sy) / (n*sx2 - sx*sx)
            b = (sy - a*sx) / n
            print()
            print(f'Ajust lineal: G*(D) ≈ {a:.2f}·log₁₀(D) + {b:.2f}')
            print(f'(millor ajust amb les dades actuals)')

    return frontera


# ──────────────────────────────────────────────────────────────────────────────
#  VISUALITZACIÓ ASCII
# ──────────────────────────────────────────────────────────────────────────────

def heatmap_ascii(resultats):
    from collections import defaultdict
    data = defaultdict(list)
    for r in resultats:
        data[(r['D'], r['G'])].append(r['success'])

    D_vals = sorted(set(r['D'] for r in resultats))
    G_vals = sorted(set(r['G'] for r in resultats))

    print()
    print('HEATMAP ASCII (P(èxit)):')
    print('  █ = 100%  ▓ = 75%  ▒ = 50%  ░ = 25%  · = 0%')
    print()
    print(f'  D\\G  ' + ' '.join(f'{g:>3}' for g in G_vals))
    print('  ' + '─'*(5+4*len(G_vals)))
    for D in D_vals:
        row = f'  {D:>4} '
        for G in G_vals:
            casos = data.get((D,G),[])
            if not casos:
                row += '  · '
                continue
            p = sum(casos)/len(casos)
            if   p >= 0.95: row += '  █ '
            elif p >= 0.75: row += '  ▓ '
            elif p >= 0.50: row += '  ▒ '
            elif p >= 0.25: row += '  ░ '
            else:           row += '  · '
        print(row)
    print()


# ──────────────────────────────────────────────────────────────────────────────
#  PUNT D'ENTRADA
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    quick = '--quick' in sys.argv or len(sys.argv) == 1

    print('█'*66)
    print('  BENCHMARK FRONTERA G*(D) — MDC Híbrid')
    print('  Qui funciona el predictiu per a factors DESEQUILIBRATS?')
    print('█'*66)
    print()
    print(f'  Mode: {"ràpid (5 mostres)" if quick else "complet (10 mostres)"}')
    print(f'  Sortida CSV: {CSV_OUTPUT}')
    print()

    resultats = run_benchmark(quick=quick)

    print()
    heatmap_ascii(resultats)
    frontera = analitzar_frontera(resultats)

    print()
    print(f'  Resultats guardats a: {CSV_OUTPUT}')
    print(f'  Total casos: {len(resultats)}')
    ok = sum(1 for r in resultats if r['success'])
    print(f'  Encerts globals: {ok}/{len(resultats)} ({ok/len(resultats)*100:.1f}%)')
    print('█'*66)
