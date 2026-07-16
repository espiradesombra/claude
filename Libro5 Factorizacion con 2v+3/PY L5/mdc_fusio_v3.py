"""
MDC FUSIÓ v3 — Método Diofántico Cinemático
============================================
Correcció clau v3:
  - Teorema 2 és LOCAL (finestra de max_local passos)
  - Quan Teorema 2 dispara → salt parabòlic (E4), NO parada global
  - Cerca bidireccional: salt parabòlic recorre k creixent
    (cada k genera un S candidat molt allunyat de l'inicial)

Autor: Víctor Manzanares Alberola (EPSA / UPV Alcoi)

TEOREMA 1: Δ(S)=k² → N=(2v+3)(2b+3), v=(S+k)/2, b=(S-k)/2
TEOREMA 2: 4S+16 > 2k+1 → no hi ha solució en la finestra local actual
           → acció: saltar a nova paràbola (salt E4)
"""

import time, math, random
from dataclasses import dataclass, field
from typing import Optional, Tuple, List

try:
    from gmpy2 import mpz, isqrt as gisqrt
    HAS_GMPY2 = True
except ImportError:
    HAS_GMPY2 = False
    mpz = int
    def gisqrt(n): return int(math.isqrt(int(n)))

try:
    import sympy as sp
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False


# ─────────────────────────────────────────────
# UTILITATS
# ─────────────────────────────────────────────

def isqrt_f(n):    return int(math.isqrt(int(n)))
def isqrt_c(n):
    r = isqrt_f(n)
    return r+1 if r*r < int(n) else r
def is_sq(n):
    if n < 0: return False, 0
    r = isqrt_f(n); return r*r == n, r
def bitlen(n):     return int(n).bit_length()
def now_ms():      return time.perf_counter() * 1000.0

def is_prime_test(n):
    if HAS_SYMPY: return sp.isprime(int(n))
    n = int(n)
    if n < 2: return False
    if n < 4: return True
    if n%2==0 or n%3==0: return False
    i = 5
    while i*i <= n:
        if n%i==0 or n%(i+2)==0: return False
        i += 6
    return True

def ecm_rescue(N):
    if not HAS_SYMPY: return None
    try:
        fs = sp.factorint(int(N))
        if len(fs) > 1:
            fl = []
            for p,e in fs.items():
                for _ in range(e): fl.append(int(p))
            if len(fl) >= 2:
                A = fl[0]; B = int(N)//A
                if A*B == int(N): return (min(A,B), max(A,B))
    except: pass
    return None


# ─────────────────────────────────────────────
# FILTRES
# ─────────────────────────────────────────────

_QR16 = {0,1,4,9}
_QRT  = {p:{(a*a)%p for a in range(p)} for p in (3,5,7,11,13)}

def quick_filter(delta):
    if delta < 0: return False
    if (int(delta)&15) not in _QR16: return False
    for p,qr in _QRT.items():
        if int(delta)%p not in qr: return False
    return True

def valid_mod4(S, M):
    return (int(M) - 6*int(S)) % 4 == 0


# ─────────────────────────────────────────────
# SENYAL φ + VECTOR DIRECTOR
# ─────────────────────────────────────────────

def compute(S, N):
    S=int(S); N=int(N)
    M = N-9
    Delta = S*S + 6*S - M
    if Delta < 0: return None
    k = isqrt_f(Delta)
    sq_d = math.sqrt(max(0, Delta))
    v    = (S + sq_d) / 2.0
    phi  = v - math.floor(v)
    dv   = 0.5 + (S+3)/(2.0*sq_d) if sq_d > 0 else 1.0
    lam  = 1.0/dv if dv > 1e-9 else 1e9
    return {'Delta':Delta,'k':k,'phi':phi,'lam':lam,'sq':k*k==Delta}


# ─────────────────────────────────────────────
# SALT ADAPTATIU (pinces 1-4)
# ─────────────────────────────────────────────

def salt_phi(phi, lam):
    """
    Pinça 1 (φ<0.04):   validar ara
    Pinça 2 (φ<0.5):    avanç fi
    Pinça 3 (φ>0.85):   prop del màxim → saltar dent
    Zona neutra:         pas=2
    """
    if phi < 0.04:              return 0, 'validar'
    elif phi < 0.5:             return max(2, int(lam*phi*0.8)), 'avanc'
    elif phi > 0.85:            return max(2, int(lam*1.2)), 'salt_dent'
    else:                       return 2, 'avanc'


# ─────────────────────────────────────────────
# SALT PARABÒLIC  S = -3 + ⌈√(k²+N)⌉
# ─────────────────────────────────────────────

def salt_par(k, N):
    T = int(k)*int(k) + int(N)
    return -3 + isqrt_c(T)


# ─────────────────────────────────────────────
# VALIDACIÓ (Teorema 1)
# ─────────────────────────────────────────────

def validar(S, k, N):
    S=int(S); k=int(k); N=int(N)
    Sp=S+k; Sm=S-k
    if Sp%2!=0 or Sm%2!=0: return None
    v=Sp//2; b=Sm//2
    if b < 0: return None          # factor trivial: B=2b+3 < 3
    A=2*v+3; B=2*b+3
    if A <= 1 or B <= 1: return None
    if A == N or B == N: return None  # factor trivial 1*N
    return (min(A,B),max(A,B)) if A*B==N else None


# ─────────────────────────────────────────────
# GRAFCET MDC FUSIÓ v3
# ─────────────────────────────────────────────

@dataclass
class MDCResult:
    N: int; bits: int
    factors: Optional[Tuple[int,int]]
    resultat: str; metode_final: str
    iteracions: int; temps_ms: float
    log: List[str] = field(default_factory=list)


def mdc_fusio(N, max_iter=20000, max_local=128, verbose=False):
    """
    GRAFCET MDC Fusió v3

    Arquitectura:
      E0: inicialitza S0 = ⌈√(N+64)⌉-3
      Bucle principal:
        E1: filtra + evalua φ + vector director
        E2: avanç fi (segueix el dent actual)
        E3: salt de dent (φ alt → salta dent complet)
        E4: Teorema 2 dispara (local) → salt parabòlic:
               k_iter++ → S_nou = -3+⌈√(k_iter²+N)⌉
        E5: Δ quadrat perfecte → Teorema 1 → validar → FINAL
      Pipeline ECM si s'esgoten tots els salts parabòlics
    """
    t0 = now_ms()
    N = int(N)
    log = []
    def vlog(m):
        if verbose: log.append(m)

    # ── E0 ───────────────────────────────────────────
    if N <= 8:
        return MDCResult(N=N,bits=bitlen(N),factors=None,
                         resultat='fail',metode_final='trivial',
                         iteracions=0,temps_ms=now_ms()-t0,log=log)
    if is_prime_test(N):
        return MDCResult(N=N,bits=bitlen(N),factors=None,
                         resultat='no_solutions',metode_final='primalitat',
                         iteracions=0,temps_ms=now_ms()-t0,log=log)

    M = N - 9
    S  = isqrt_c(N+64) - 3
    if S < 1: S = 1
    for d in range(8):
        if valid_mod4(S+d, M): S=S+d; break

    # k_iter per al salt parabòlic (comença a 0, s'incrementa ABANS d'usar)
    # k=1 → S = -3+ceil(√(1+N)) ≈ √N-2  (factor equilibrat)
    # k gran → S gran → factor desigual
    k_iter   = 0
    local_n  = 0      # passos des de l'últim salt parabòlic
    it       = 0
    phi_prv  = None

    vlog(f"E0: N={N}, S0={S}")

    # Calculem rang màxim de k: fins k ≈ (N-9)/6 (factor mínim p=3)
    k_max = int(N) // 6 + 10

    while it < max_iter and k_iter <= k_max:
        it += 1

        # Filtre mod4
        if not valid_mod4(S, M):
            S += 1; local_n += 1; continue

        # Càlcul senyal
        info = compute(S, N)
        if info is None:
            S += 2; local_n += 1; continue

        Delta = info['Delta']
        k     = info['k']
        phi   = info['phi']
        lam   = info['lam']

        # Filtre Δ ràpid
        if not quick_filter(Delta):
            S += 2; local_n += 1; continue

        # ── E5: Δ quadrat perfecte → Teorema 1 ──────────────
        if info['sq']:
            facs = validar(S, k, N)
            if facs:
                vlog(f"E5✔ S={S},k={k} → {facs}")
                return MDCResult(N=N,bits=bitlen(N),factors=facs,
                                 resultat='factor_found',metode_final='mdc_fusio',
                                 iteracions=it,temps_ms=now_ms()-t0,log=log)
            S += 2; local_n += 1; continue

        # ── TEOREMA 2 (LOCAL): si dispara → E4 salt parabòlic ──
        if (4*S + 16) > (2*k + 1) or local_n >= max_local:
            # Salt parabòlic: nova k → nou S directament
            k_iter += 1
            S_nou = salt_par(k_iter, N)
            if S_nou < 1: S_nou = 1
            # Ajust mod4
            for d in range(8):
                if valid_mod4(S_nou+d, M):
                    S = S_nou+d; break
            else:
                S = S_nou
            local_n = 0
            phi_prv = None
            vlog(f"E4: k_iter={k_iter} → S={S}")
            continue

        # ── E1: AVALUACIÓ PINCES ─────────────────────────────
        vlog(f"E1: S={S},φ={phi:.3f},λ={lam:.2f},k={k}")

        # ── E2/E3: SALT ADAPTATIU ────────────────────────────
        pas, accio = salt_phi(phi, lam)
        if pas == 0:
            # φ≈0 però Δ no era quadrat → avanç
            S += 2; local_n += 1; phi_prv = phi; continue

        # Ajust paritat
        if pas % 2 != 0: pas += 1
        S       += pas
        local_n += pas
        phi_prv  = phi

    # ── ECM rescat ───────────────────────────────────────────
    vlog("Rescat ECM")
    facs = ecm_rescue(N)
    if facs:
        return MDCResult(N=N,bits=bitlen(N),factors=facs,
                         resultat='ecm_rescue',metode_final='ecm',
                         iteracions=it,temps_ms=now_ms()-t0,log=log)

    return MDCResult(N=N,bits=bitlen(N),factors=None,
                     resultat='fail',metode_final='none',
                     iteracions=it,temps_ms=now_ms()-t0,log=log)


# ─────────────────────────────────────────────
# BENCHMARK
# ─────────────────────────────────────────────

def gen_semi(bits):
    if not HAS_SYMPY: return None,None,None
    h = bits//2
    p = int(sp.nextprime(random.getrandbits(h)|(1<<(h-1))))
    q = int(sp.nextprime(random.getrandbits(bits-h)|(1<<(bits-h-1))))
    return p*q, p, q

def gen_semi_unbal(bits, small_bits=16):
    if not HAS_SYMPY: return None,None,None
    p = int(sp.nextprime(random.getrandbits(small_bits)|(1<<(small_bits-1))))
    qb = bits - p.bit_length()
    if qb < 8: qb = 8
    q = int(sp.nextprime(random.getrandbits(qb)|(1<<(qb-1))))
    return p*q, p, q


def benchmark():
    print("\n" + "="*70)
    print("  BENCHMARK MDC FUSIÓ v3")
    print("="*70)
    print(f"{'Tipus':>10} {'Bits':>5} {'N':>22} {'Resultat':>15} {'It':>7} {'ms':>8} {'Mètode':>10}")
    print("-"*70)
    contadors = {'mdc_fusio':0,'ecm':0,'fail':0}

    for bits in (20, 28, 36, 44):
        # Equilibrats
        for _ in range(4):
            N,p,q = gen_semi(bits)
            if N is None: continue
            res = mdc_fusio(N, max_iter=30000)
            nm = str(N)[:20]
            ok = "✔" if res.resultat in ('factor_found','ecm_rescue') else "✗"
            m  = res.metode_final
            contadors[m if m in contadors else 'fail'] += 1
            print(f"{'equil':>10} {bits:>5} {nm:>22} {ok+' '+res.resultat:>15} "
                  f"{res.iteracions:>7} {res.temps_ms:>8.2f} {m:>10}")

        # Desequilibrats
        for _ in range(4):
            N,p,q = gen_semi_unbal(bits, small_bits=min(16,bits//2))
            if N is None: continue
            res = mdc_fusio(N, max_iter=30000)
            nm = str(N)[:20]
            ok = "✔" if res.resultat in ('factor_found','ecm_rescue') else "✗"
            m  = res.metode_final
            contadors[m if m in contadors else 'fail'] += 1
            print(f"{'desig':>10} {bits:>5} {nm:>22} {ok+' '+res.resultat:>15} "
                  f"{res.iteracions:>7} {res.temps_ms:>8.2f} {m:>10}")

    print("="*70)
    print(f"  MDC directe: {contadors['mdc_fusio']}  ECM: {contadors['ecm']}  Fallits: {contadors['fail']}")
    print("="*70)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("MDC FUSIÓ v3")
    print("="*65)

    casos = [
        # (N, p, q, etiqueta)
        (5*7,           5,    7,    "equil-petit"),
        (11*13,        11,   13,    "equil-petit"),
        (29*31,        29,   31,    "equil"),
        (101*103,     101,  103,    "equil"),
        (9973*9967,  9973, 9967,    "equil-gran"),
        (7*1009,        7, 1009,    "desig"),
        (11*9973,      11, 9973,    "desig"),
        (13*99991,     13,99991,    "desig"),
        (17*999983,    17,999983,   "desig-gran"),
        (5*1000003,     5,1000003,  "desig-5"),
    ]

    for N, p, q, etiq in casos:
        res = mdc_fusio(N, max_iter=50000, verbose=False)
        mdc_ok = res.metode_final == 'mdc_fusio'
        tag = "✔ MDC" if mdc_ok else ("✔ ECM" if res.resultat=='ecm_rescue' else "✗ FAIL")
        print(f"N={N:>12} ({p}×{q}) [{etiq:>12}] | {tag} | "
              f"iters={res.iteracions:>6} | {res.temps_ms:>6.2f}ms | {res.factors}")

    if HAS_SYMPY:
        benchmark()
    else:
        print("\n(Instal·la sympy per al benchmark)")
