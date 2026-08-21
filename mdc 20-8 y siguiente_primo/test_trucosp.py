import sympy, math, time
import sys
sys.path.insert(0, '/home/claude')
from mdc_maestro import mdcMaestro, _salto_predictivo, _criba_p210

random_seed = 7
import random
random.seed(random_seed)

# Caso desbalanceado real "hueco": p no lo coge F0b, no lo coge F1
p = sympy.nextprime(10**14 + 7)      # ~15 dígitos
q = sympy.nextprime(10**24 + 91)     # ~25 dígitos
N = p * q
print(f"N = p*q  ->  {len(str(N))} dígitos  (p {len(str(p))}d, q {len(str(q))}d)")
print(f"p-q = {q-p}  ¿(p-q) > sqrt(N)? {abs(p-q) > math.isqrt(N)}")
print()

# --- 1) intento directo sobre N (referencia, sabemos que falla / es caro) ---
t0=time.time()
res_directo = mdcMaestro(N, max_ev_sp=300_000)
print("Directo sobre N:", res_directo, f"  ({time.time()-t0:.2f}s)")
print()

# --- 2) truco: N2 = N * randprime(tamaño N) ---
s = sympy.randprime(10**38, 10**39)   # primo aleatorio del mismo orden que N (39 dígitos)
N2 = N * s
print(f"s (random prime): {len(str(s))} dígitos")
print(f"N2 = N*s: {len(str(N2))} dígitos, sqrt(N2) ~ {len(str(math.isqrt(N2)))} dígitos")
print()

t0=time.time()
res_N2 = mdcMaestro(N2, max_ev_sp=300_000)
print("mdcMaestro sobre N2 = N*s:", res_N2, f"  ({time.time()-t0:.2f}s)")

if res_N2[0]:
    f = res_N2[0]
    print("factor encontrado:", f)
    print("  ¿es p?", f == p, "  ¿es q?", f == q, "  ¿es s?", f == s)
