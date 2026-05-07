"""
MDC Factorizador — Método Diofántico Cinemático
Víctor Manzanares Alberola, 2026
------------------------------------------------
Busca factores de N en el espacio L1 = {6k±1}
usando la condición frac(N / (2*(2m+3))) = 0.5
con predicción cinemática de 4 puntos.

Salida: traza detallada + gráfica del perfil d(m).
"""

import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import time


# ─── Aritmética exacta ────────────────────────────────────────────────────────

def d(m, N):
    """frac(N / (2*(2m+3))) usando aritmética entera — sin errores float."""
    denom = 2 * (2 * m + 3)
    return (N % denom) / denom          # exactamente frac(N/denom)

def es_L1(m):
    return m >= 1 and (2 * m + 3) % 3 != 0

def sig_L1(m):
    m += 1
    while not es_L1(m):
        m += 1
    return m

def prev_L1(m):
    m -= 1
    while m >= 1 and not es_L1(m):
        m -= 1
    return m if m >= 1 else None


# ─── Predicción cinemática ────────────────────────────────────────────────────

def predecir(m0, N, m_max):
    """
    Toma 4 puntos consecutivos en L1 desde m0.
    Calcula velocidades y predice m* donde d(m*)=0.5.

    Devuelve: (m_pred, d0, vs, tipo_salto)
    tipo_salto: 'ok' | 'reset' | 'plano'
    """
    pts = [m0]
    for _ in range(3):
        nxt = sig_L1(pts[-1])
        if nxt > m_max:
            break
        pts.append(nxt)

    if len(pts) < 4:
        return None, None, None, 'borde'

    ds = [d(mi, N) for mi in pts]
    vs = [ds[i+1] - ds[i] for i in range(3)]

    # Detectar reset de diente de sierra
    if any(v < -0.15 for v in vs):
        return None, ds[0], vs, 'reset'

    vs_pos = [v for v in vs if v > 1e-10]
    if not vs_pos:
        return None, ds[0], vs, 'plano'

    v_media = sum(vs_pos) / len(vs_pos)
    delta = 0.5 - ds[0]
    paso_m = (pts[-1] - pts[0]) / (len(pts) - 1)   # paso medio en m
    m_pred = pts[0] + round(delta / v_media * paso_m)

    # Snap al candidato L1 más cercano
    if m_pred < pts[-1]:
        m_pred = pts[-1]
    while not es_L1(m_pred):
        m_pred += 1

    return m_pred, ds[0], vs, 'ok'


# ─── Verificación local ───────────────────────────────────────────────────────

def verificar_zona(m_centro, N, m_max, radio=6):
    """
    Evalúa todos los candidatos L1 en [m_centro-radio, m_centro+radio].
    Devuelve (m_factor, evaluaciones) o (None, evaluaciones).
    """
    ev = 0
    m = max(1, m_centro - radio)
    while not es_L1(m) and m <= m_centro + radio:
        m += 1
    while m <= m_centro + radio and m <= m_max:
        ev += 1
        if d(m, N) == 0.5:
            f = 2 * m + 3
            if N % f == 0:
                return m, ev
        m = sig_L1(m)
    return None, ev


# ─── Motor principal MDC ──────────────────────────────────────────────────────

def mdc_factorizar(N, verbose=True):
    """
    Factoriza N mediante el MDC.

    Devuelve:
        factor, cofactor, estadísticas, registro de saltos
    """
    if N < 4:
        return N, 1, {}, []
    if N % 2 == 0:
        return 2, N // 2, {'evaluaciones': 1, 'saltos': 0}, []
    if N % 3 == 0:
        return 3, N // 3, {'evaluaciones': 1, 'saltos': 0}, []

    m_max = (int(math.isqrt(N)) - 3) // 2 + 2

    m = 1
    while not es_L1(m):
        m += 1

    evaluaciones = 0
    saltos = 0
    registro = []          # lista de dicts con info de cada salto

    if verbose:
        print(f"\n{'═'*60}")
        print(f"  MDC Factorizador — N = {N:,}")
        print(f"  m_max = {m_max}  (factor mínimo posible: {2*m_max+3})")
        print(f"{'═'*60}")
        print(f"  {'Salto':>5}  {'m0':>6}  {'f0':>6}  {'d0':>8}  "
              f"{'v_media':>9}  {'m_pred':>8}  {'error':>7}  tipo")
        print(f"  {'-'*56}")

    while m <= m_max:
        # ── Verificación directa antes de saltar ──────────────────────────
        evaluaciones += 1
        di = d(m, N)
        if di == 0.5:
            f = 2 * m + 3
            if N % f == 0:
                registro.append({
                    'salto': saltos, 'm0': m, 'f0': 2*m+3,
                    'd0': di, 'tipo': 'DIRECTO', 'm_pred': m, 'error': 0
                })
                break

        # ── Predicción cinemática ─────────────────────────────────────────
        m_pred, d0, vs, tipo = predecir(m, N, m_max)
        evaluaciones += 4   # 4 puntos evaluados en predecir()

        if tipo == 'reset':
            # Avanzar al siguiente nivel
            for _ in range(4):
                m = sig_L1(m)
            saltos += 1
            registro.append({
                'salto': saltos, 'm0': m, 'f0': 2*m+3,
                'd0': d0 if d0 is not None else 0,
                'tipo': 'RESET', 'm_pred': m, 'error': 0
            })
            continue

        if tipo in ('plano', 'borde'):
            m = sig_L1(m)
            continue

        # ── Verificar zona alrededor de la predicción ─────────────────────
        m_factor, ev_local = verificar_zona(m_pred, N, m_max, radio=8)
        evaluaciones += ev_local
        saltos += 1

        error_pred = abs(m_pred - (m_factor if m_factor else m_pred))
        v_media = sum(v for v in vs if v > 0) / max(1, sum(1 for v in vs if v > 0))

        entrada_log = {
            'salto': saltos, 'm0': m, 'f0': 2*m+3,
            'd0': d0, 'v_media': v_media,
            'm_pred': m_pred,
            'error': error_pred,
            'tipo': 'ENCONTRADO' if m_factor else 'fallido'
        }
        registro.append(entrada_log)

        if verbose:
            tipo_str = '✓ ENCONTRADO' if m_factor else '  →'
            print(f"  {saltos:>5}  {m:>6}  {2*m+3:>6}  {d0:>8.4f}  "
                  f"{v_media:>9.5f}  {m_pred:>8}  {error_pred:>7}  {tipo_str}")

        if m_factor:
            f = 2 * m_factor + 3
            m = m_factor
            break

        # Avanzar desde la predicción
        m = m_pred if m_pred > m else sig_L1(m)

    else:
        # No encontrado (N es primo)
        stats = {
            'evaluaciones': evaluaciones,
            'saltos': saltos,
            'factor': None,
            'cofactor': None,
            'primo': True
        }
        if verbose:
            print(f"\n  → N = {N:,} es primo (o sin factores en L1)")
            print(f"  Evaluaciones: {evaluaciones}  |  Saltos: {saltos}")
        return None, None, stats, registro

    f = 2 * m + 3
    g = N // f
    stats = {
        'evaluaciones': evaluaciones,
        'saltos': saltos,
        'factor': f,
        'cofactor': g,
        'primo': False
    }

    if verbose:
        print(f"\n  ✓ Factor encontrado: {f} × {g} = {N:,}")
        print(f"  Evaluaciones: {evaluaciones}  |  Saltos: {saltos}")
        print(f"{'═'*60}\n")

    return f, g, stats, registro


# ─── Gráfica del perfil d(m) ─────────────────────────────────────────────────

def graficar(N, registro, ruta_salida="mdc_perfil.png"):
    """
    Dibuja el perfil completo d(m) para todos los candidatos L1 hasta √N,
    marcando el factor encontrado, la predicción y la trayectoria de saltos.
    """
    m_max = (int(math.isqrt(N)) - 3) // 2 + 3

    # Recoger todos los puntos L1
    ms_todo, ds_todo = [], []
    m = 1
    while m <= m_max:
        if es_L1(m):
            ms_todo.append(m)
            ds_todo.append(d(m, N))
        m += 1

    fig, axes = plt.subplots(2, 1, figsize=(13, 9),
                              gridspec_kw={'height_ratios': [3, 1]})
    ax, ax2 = axes
    fig.patch.set_facecolor('#0f1117')
    for a in axes:
        a.set_facecolor('#161b22')
        a.tick_params(colors='#8b949e')
        for spine in a.spines.values():
            spine.set_edgecolor('#30363d')

    # ── Panel superior: perfil completo ──────────────────────────────────────
    ax.plot(ms_todo, ds_todo,
            color='#388bfd', linewidth=0.6, alpha=0.5, zorder=1, label='d(m)')
    ax.scatter(ms_todo, ds_todo,
               c='#388bfd', s=8, alpha=0.6, zorder=2)

    # Línea objetivo δ=0.5
    ax.axhline(0.5, color='#f78166', linewidth=1.2,
               linestyle='--', alpha=0.8, label='objetivo δ=0.5')

    # Marcar factor(es) exactos
    for mi, di in zip(ms_todo, ds_todo):
        if di == 0.5:
            ax.scatter([mi], [di], c='#3fb950', s=120, zorder=5, marker='*')
            ax.annotate(f'  f={2*mi+3}',
                        xy=(mi, di), fontsize=9, color='#3fb950',
                        va='bottom')

    # Marcar trayectoria de saltos
    colores_tipo = {'RESET': '#d29922', 'fallido': '#f78166', 'DIRECTO': '#3fb950'}
    for r in registro:
        if r['tipo'] in ('RESET', 'fallido'):
            ax.axvline(r['m0'], color=colores_tipo[r['tipo']],
                       linewidth=0.8, alpha=0.4, linestyle=':')
        if 'm_pred' in r and r['tipo'] not in ('RESET', 'DIRECTO'):
            ax.axvline(r['m_pred'], color='#bc8cff',
                       linewidth=0.8, alpha=0.35, linestyle='-.')

    ax.set_ylabel('d(m)', color='#c9d1d9', fontsize=11)
    ax.set_title(f'MDC — Perfil d(m) para N = {N:,}', color='#c9d1d9', fontsize=13, pad=10)
    ax.set_xlim(0, m_max + 1)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, color='#21262d', linewidth=0.5)

    legend_elements = [
        mpatches.Patch(facecolor='#388bfd', label='d(m) en L1'),
        mpatches.Patch(facecolor='#3fb950', label='factor exacto'),
        mpatches.Patch(facecolor='#f78166', label='objetivo δ=0.5'),
        mpatches.Patch(facecolor='#d29922', alpha=0.6, label='avance RESET'),
        mpatches.Patch(facecolor='#bc8cff', alpha=0.6, label='predicción m_pred'),
    ]
    ax.legend(handles=legend_elements, loc='upper right',
              facecolor='#161b22', edgecolor='#30363d',
              labelcolor='#c9d1d9', fontsize=9)

    # ── Panel inferior: error de predicción por salto ─────────────────────────
    saltos_ok = [r for r in registro if r['tipo'] not in ('RESET', 'DIRECTO')]
    if saltos_ok:
        xs = [r['salto'] for r in saltos_ok]
        errs = [r['error'] for r in saltos_ok]
        colores_err = ['#3fb950' if e == 0 else '#f78166' if e > 5 else '#d29922'
                       for e in errs]
        ax2.bar(xs, errs, color=colores_err, alpha=0.85, width=0.6)
        ax2.set_xlabel('Número de salto', color='#c9d1d9', fontsize=10)
        ax2.set_ylabel('Error (pasos)', color='#c9d1d9', fontsize=10)
        ax2.set_title('Error de predicción por salto cinemático',
                      color='#8b949e', fontsize=10, pad=6)
        ax2.grid(True, color='#21262d', linewidth=0.5, axis='y')
        if max(errs) > 0:
            ax2.set_ylim(0, max(errs) * 1.3)
    else:
        ax2.text(0.5, 0.5, 'Sin saltos registrados', color='#8b949e',
                 ha='center', va='center', transform=ax2.transAxes)

    plt.tight_layout(pad=1.5)
    plt.savefig(ruta_salida, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Gráfica guardada: {ruta_salida}")


# ─── Comparación con división de prueba ──────────────────────────────────────

def trial_division_L1(N):
    """Cuenta evaluaciones de división de prueba en espacio L1."""
    m_max = (int(math.isqrt(N)) - 3) // 2 + 2
    ev = 0
    m = 1
    while m <= m_max:
        if es_L1(m):
            ev += 1
            if N % (2 * m + 3) == 0:
                f = 2 * m + 3
                return f, N // f, ev
        m += 1
    return None, None, ev


# ─── Demo principal ───────────────────────────────────────────────────────────

def demo(N):
    t0 = time.perf_counter()
    f, g, stats, registro = mdc_factorizar(N, verbose=True)
    t_mdc = time.perf_counter() - t0

    _, _, ev_td = trial_division_L1(N)

    print(f"  División de prueba (L1): {ev_td} evaluaciones")
    print(f"  MDC:                     {stats['evaluaciones']} evaluaciones")
    ratio = ev_td / stats['evaluaciones'] if stats['evaluaciones'] else 0
    print(f"  Ratio TD/MDC:            {ratio:.2f}×")
    print(f"  Tiempo MDC:              {t_mdc*1000:.2f} ms\n")

    nombre = f"mdc_N{N}.png"
    graficar(N, registro, ruta_salida=nombre)
    return nombre


if __name__ == "__main__":
    # Casos de prueba
    casos = [
        143,        # 11 × 13
        3127,       # 53 × 59
        10403,      # 101 × 103
        47053,      # 211 × 223
        1656119,    # 1283 × 1291
    ]

    archivos = []
    for N in casos:
        arch = demo(N)
        archivos.append(arch)

    print("\nArchivos generados:")
    for a in archivos:
        print(f"  {a}")
