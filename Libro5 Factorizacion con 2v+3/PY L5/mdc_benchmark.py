"""
MDC Benchmark Masivo — Víctor Manzanares Alberola, 2026
-------------------------------------------------------
Compara MDC vs División de Prueba (L1) en tres categorías:
  A) Semiprimos desbalanceados: p pequeño (5–199), q ≈ N/p
  B) Primos:                    peor caso para MDC
  C) Factor en posición aleatoria: p ≈ N^alpha, alpha ∈ (0.1, 0.9)

Rango: 4–18 dígitos.
Salida: benchmark_mdc.csv + benchmark_mdc.png
"""

import math
import time
import csv
import random
import sympy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

random.seed(42)

# ══════════════════════════════════════════════════════════════════════════════
# Motor MDC (inline, autocontenido)
# ══════════════════════════════════════════════════════════════════════════════

def d(m, N):
    denom = 2 * (2 * m + 3)
    return (N % denom) / denom

def es_L1(m):
    return m >= 1 and (2 * m + 3) % 3 != 0

def sig_L1(m):
    m += 1
    while not es_L1(m):
        m += 1
    return m

def mdc_factorizar(N):
    if N < 4:      return N, 1, 1, 0
    if N % 2 == 0: return 2, N//2, 1, 0
    if N % 3 == 0: return 3, N//3, 1, 0

    m_max = (int(math.isqrt(N)) - 3) // 2 + 2
    m = 1
    while not es_L1(m):
        m += 1

    ev = 0
    saltos = 0

    while m <= m_max:
        ev += 1
        if d(m, N) == 0.5:
            f = 2*m+3
            if N % f == 0:
                return f, N//f, ev, saltos

        # 4 puntos para predicción
        pts = [m]
        for _ in range(3):
            nxt = sig_L1(pts[-1])
            if nxt > m_max: break
            pts.append(nxt)

        if len(pts) == 4:
            ds = [d(mi, N) for mi in pts]
            ev += 3  # ya contamos el primero
            vs = [ds[i+1]-ds[i] for i in range(3)]

            if not any(v < -0.15 for v in vs):  # sin reset
                vs_pos = [v for v in vs if v > 1e-10]
                if vs_pos:
                    v_med = sum(vs_pos)/len(vs_pos)
                    delta = 0.5 - ds[0]
                    paso = (pts[-1]-pts[0])/3
                    m_pred = int(round(pts[0] + delta/v_med * paso))
                    if m_pred < pts[-1]: m_pred = pts[-1]
                    while not es_L1(m_pred): m_pred += 1

                    if m_pred <= m_max:
                        # verificar zona ±8
                        mc = max(1, m_pred - 8)
                        while not es_L1(mc): mc += 1
                        while mc <= min(m_pred+8, m_max):
                            ev += 1
                            if d(mc, N) == 0.5:
                                f = 2*mc+3
                                if N % f == 0:
                                    return f, N//f, ev, saltos+1
                            mc = sig_L1(mc)
                        saltos += 1
                        m = m_pred if m_pred > pts[-1] else sig_L1(pts[-1])
                        continue

        m = sig_L1(m)

    return None, None, ev, saltos


def td_factorizar(N):
    if N % 2 == 0: return 2, N//2, 1
    if N % 3 == 0: return 3, N//3, 1
    m_max = (int(math.isqrt(N))-3)//2+2
    ev = 0
    m = 1
    while m <= m_max:
        if es_L1(m):
            ev += 1
            f = 2*m+3
            if N % f == 0:
                return f, N//f, ev
        m += 1
    return None, None, ev


# ══════════════════════════════════════════════════════════════════════════════
# Generadores de casos
# ══════════════════════════════════════════════════════════════════════════════

PRIMOS_PEQUENOS = [p for p in sympy.primerange(5, 200) if p % 6 in (1, 5)]

def gen_desbalanceado(n_digitos):
    """p pequeño (5–199), q tal que N tiene ~n_digitos."""
    p = random.choice(PRIMOS_PEQUENOS)
    lo = 10**(n_digitos-1) // p
    hi = (10**n_digitos - 1) // p
    if lo > hi or lo < 2: return None
    try:
        q = sympy.randprime(max(lo, p+1), hi+1)
    except Exception:
        return None
    if q is None: return None
    N = p * q
    if len(str(N)) != n_digitos: return None
    return N, p, q, 'desbal'

def gen_primo(n_digitos):
    """Primo de n_digitos exactos."""
    lo, hi = 10**(n_digitos-1), 10**n_digitos - 1
    N = sympy.randprime(lo, hi)
    return N, None, None, 'primo'

def gen_aleatorio(n_digitos):
    """p ≈ N^alpha con alpha aleatorio en (0.15, 0.85)."""
    alpha = random.uniform(0.15, 0.85)
    lo = 10**(n_digitos-1)
    hi = 10**n_digitos - 1
    p_target = int(lo**alpha)
    if p_target < 5: p_target = 5
    p = sympy.nextprime(p_target)
    if p % 6 not in (1, 5): p = sympy.nextprime(p)
    q_lo = lo // p
    q_hi = hi // p
    if q_lo > q_hi or q_lo < p: return None
    try:
        q = sympy.randprime(max(q_lo, p), q_hi+1)
    except Exception:
        return None
    if q is None: return None
    N = p * q
    if len(str(N)) != n_digitos: return None
    return N, p, q, f'aleatorio(α={alpha:.2f})'


# ══════════════════════════════════════════════════════════════════════════════
# Ejecución del benchmark
# ══════════════════════════════════════════════════════════════════════════════

DIGITOS   = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
REPS      = 3   # repeticiones por (categoría, n_digitos)

generadores = [
    ('Desbalanceado', gen_desbalanceado),
    ('Primo',         gen_primo),
    ('Aleatorio',     gen_aleatorio),
]

filas = []
print(f"\n{'N_dig':>5} {'Categoría':>14} {'N':>20} {'f':>12} "
      f"{'ev_MDC':>8} {'ev_TD':>7} {'ratio':>7} {'saltos':>7} {'ok':>4}")
print('─'*90)

for nd in DIGITOS:
    for cat_nombre, gen_fn in generadores:
        exitos = 0
        for _ in range(REPS * 3):          # intentos extras por si gen falla
            if exitos >= REPS: break
            res = gen_fn(nd)
            if res is None: continue
            N, p_true, q_true, subtipo = res

            t0 = time.perf_counter()
            f_mdc, g_mdc, ev_mdc, saltos = mdc_factorizar(N)
            t_mdc = time.perf_counter() - t0

            t0 = time.perf_counter()
            f_td, g_td, ev_td = td_factorizar(N)
            t_td = time.perf_counter() - t0

            ok_mdc = (f_mdc is not None and f_mdc * g_mdc == N) or \
                     (f_mdc is None and sympy.isprime(N))
            ok_td  = (f_td  is not None and f_td  * g_td  == N) or \
                     (f_td  is None and sympy.isprime(N))

            ratio = ev_td / ev_mdc if ev_mdc > 0 else float('nan')
            ok_str = '✓' if ok_mdc else '✗'

            fila = {
                'n_digitos':  nd,
                'categoria':  cat_nombre,
                'subtipo':    subtipo,
                'N':          N,
                'factor_real': p_true,
                'factor_MDC': f_mdc,
                'ev_MDC':     ev_mdc,
                'ev_TD':      ev_td,
                'ratio_TD_MDC': round(ratio, 3),
                'saltos_MDC': saltos,
                't_MDC_ms':   round(t_mdc*1000, 3),
                't_TD_ms':    round(t_td*1000, 3),
                'MDC_correcto': ok_mdc,
                'TD_correcto':  ok_td,
            }
            filas.append(fila)
            exitos += 1

            N_str = str(N)[:18] + ('…' if len(str(N))>18 else '')
            f_str = str(f_mdc)[:10] if f_mdc else 'primo'
            print(f"{nd:>5} {cat_nombre:>14} {N_str:>20} {f_str:>12} "
                  f"{ev_mdc:>8} {ev_td:>7} {ratio:>7.2f} {saltos:>7} {ok_str:>4}")


# ══════════════════════════════════════════════════════════════════════════════
# Exportar CSV
# ══════════════════════════════════════════════════════════════════════════════

csv_path = 'benchmark_mdc.csv'
campos = list(filas[0].keys())
with open(csv_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=campos)
    w.writeheader()
    w.writerows(filas)
print(f'\n✓ CSV guardado: {csv_path}  ({len(filas)} filas)')


# ══════════════════════════════════════════════════════════════════════════════
# Gráfica
# ══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(16, 6))
fig.patch.set_facecolor('#0f1117')

CAT_COLOR = {
    'Desbalanceado': '#3fb950',
    'Primo':         '#f78166',
    'Aleatorio':     '#d29922',
}

for ax, cat in zip(axes, ['Desbalanceado', 'Primo', 'Aleatorio']):
    ax.set_facecolor('#161b22')
    for spine in ax.spines.values():
        spine.set_edgecolor('#30363d')
    ax.tick_params(colors='#8b949e')

    sub = [r for r in filas if r['categoria'] == cat]
    # Agrupar por n_digitos: mediana del ratio
    nd_vals = sorted(set(r['n_digitos'] for r in sub))
    ratios_med = []
    ratios_p25 = []
    ratios_p75 = []
    ev_mdc_med = []
    ev_td_med  = []

    for nd in nd_vals:
        rs = [r['ratio_TD_MDC'] for r in sub if r['n_digitos']==nd and not math.isnan(r['ratio_TD_MDC'])]
        em = [r['ev_MDC'] for r in sub if r['n_digitos']==nd]
        et = [r['ev_TD']  for r in sub if r['n_digitos']==nd]
        if rs:
            ratios_med.append(float(np.median(rs)))
            ratios_p25.append(float(np.percentile(rs, 25)))
            ratios_p75.append(float(np.percentile(rs, 75)))
        else:
            ratios_med.append(float('nan'))
            ratios_p25.append(float('nan'))
            ratios_p75.append(float('nan'))
        ev_mdc_med.append(float(np.median(em)) if em else float('nan'))
        ev_td_med.append(float(np.median(et))  if et else float('nan'))

    col = CAT_COLOR[cat]
    ax.plot(nd_vals, ev_mdc_med, color=col, linewidth=2, marker='o',
            markersize=5, label='MDC', zorder=3)
    ax.plot(nd_vals, ev_td_med, color='#388bfd', linewidth=2,
            marker='s', markersize=5, linestyle='--', label='TD (L1)', zorder=3)
    ax.fill_between(nd_vals,
                    [a*b for a,b in zip(ev_mdc_med, ratios_p25)],
                    [a*b for a,b in zip(ev_mdc_med, ratios_p75)],
                    alpha=0.12, color=col)

    # Ratio como anotación
    ax2 = ax.twinx()
    ax2.set_facecolor('#161b22')
    ax2.tick_params(colors='#8b949e')
    ax2.plot(nd_vals, ratios_med, color='#bc8cff', linewidth=1.2,
             linestyle=':', marker='^', markersize=4, label='ratio TD/MDC')
    ax2.axhline(1.0, color='#bc8cff', linewidth=0.6, alpha=0.4)
    ax2.set_ylabel('ratio TD/MDC', color='#bc8cff', fontsize=9)
    ax2.tick_params(axis='y', colors='#bc8cff')
    ax2.set_ylim(0, max(max(r for r in ratios_med if not math.isnan(r)), 1.5) * 1.2)

    ax.set_yscale('log')
    ax.set_xlabel('dígitos de N', color='#c9d1d9', fontsize=10)
    ax.set_ylabel('evaluaciones (mediana)', color='#c9d1d9', fontsize=10)
    ax.set_title(f'Categoría: {cat}', color='#c9d1d9', fontsize=12, pad=8)
    ax.grid(True, color='#21262d', linewidth=0.5, which='both')
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2))

    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labs1+labs2, loc='upper left',
              facecolor='#161b22', edgecolor='#30363d',
              labelcolor='#c9d1d9', fontsize=9)

fig.suptitle('MDC vs División de Prueba — Benchmark honesto (4–18 dígitos)',
             color='#c9d1d9', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('benchmark_mdc.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.close()
print('✓ Gráfica guardada: benchmark_mdc.png')
