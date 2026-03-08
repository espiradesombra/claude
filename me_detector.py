"""
MEcuation Factor Detector
=========================
Dado un conjunto de números generados por una regla, detecta si existe
una MEcuation local (colapso SVD) y extrae la ecuación en los factores.

Familias soportadas:
  - cuadrados: E = p^2
  - cubos: E = p^3
  - semiprimos: E = p * q
  - k_primo: E = k * p  (k fijo)
  - mersenne: E = 2^n - 1
  - sophie_germain: E = p donde 2p+1 también es primo
  - custom: cualquier lista
"""

import numpy as np
from itertools import combinations
import math

# ──────────────────────────────────────────────
# 1. GENERADORES DE FAMILIAS
# ──────────────────────────────────────────────

def es_primo(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

def primos_hasta(n):
    return [p for p in range(2, n+1) if es_primo(p)]

def generar_familia(regla, n_elementos=15, **kwargs):
    """
    Genera números según la regla dada.
    Devuelve lista de (E, metadata) donde metadata describe los factores reales.
    """
    primos = primos_hasta(500)
    datos = []

    if regla == "cuadrados":
        # E = p^2
        for p in primos[:n_elementos]:
            datos.append((p**2, {"tipo": "p^2", "p": p}))

    elif regla == "cubos":
        # E = p^3
        for p in primos[:n_elementos]:
            datos.append((p**3, {"tipo": "p^3", "p": p}))

    elif regla == "semiprimos":
        # E = p * q, p < q ambos primos
        count = 0
        for i, p in enumerate(primos):
            for q in primos[i+1:]:
                if count >= n_elementos: break
                datos.append((p*q, {"tipo": "p*q", "p": p, "q": q}))
                count += 1
            if count >= n_elementos: break

    elif regla == "semiprimos_cercanos":
        # E = p * q donde q = siguiente primo después de p
        for i in range(min(n_elementos, len(primos)-1)):
            p, q = primos[i], primos[i+1]
            datos.append((p*q, {"tipo": "p*q_cercanos", "p": p, "q": q}))

    elif regla == "k_primo":
        # E = k * p, k fijo
        k = kwargs.get("k", 6)
        for p in primos[:n_elementos]:
            datos.append((k*p, {"tipo": f"{k}*p", "k": k, "p": p}))

    elif regla == "mersenne":
        # E = 2^n - 1
        count = 0
        for n in range(2, 60):
            if count >= n_elementos: break
            val = 2**n - 1
            datos.append((val, {"tipo": "2^n-1", "n": n}))
            count += 1

    elif regla == "potencia_mixta":
        # E = p^2 * q
        count = 0
        for i, p in enumerate(primos):
            for q in primos:
                if q == p: continue
                if count >= n_elementos: break
                datos.append((p**2 * q, {"tipo": "p^2*q", "p": p, "q": q}))
                count += 1
            if count >= n_elementos: break

    elif regla == "sophie_germain":
        # p primo de Sophie Germain: p y 2p+1 ambos primos
        count = 0
        for p in primos:
            if count >= n_elementos: break
            if es_primo(2*p + 1):
                datos.append((p, {"tipo": "sophie", "p": p, "2p+1": 2*p+1}))
                count += 1

    elif regla == "custom":
        lista = kwargs.get("lista", [])
        for e in lista:
            datos.append((e, {"tipo": "custom"}))

    return datos

# ──────────────────────────────────────────────
# 2. FEATURE MAP
# ──────────────────────────────────────────────

def feature_map(E):
    """
    Transforma E en vector de features logarítmicas.
    Diseñado para capturar relaciones entre factores.
    """
    if E <= 0: return None
    lE = math.log(E)
    features = [
        lE,                          # f1: log(E)
        lE / 2,                      # f2: log(E)/2  → log(√E)
        lE / 3,                      # f3: log(E)/3  → log(∛E)
        lE ** 2,                     # f4: log(E)^2
        math.sqrt(lE),               # f5: √log(E)
        lE / math.log(2),            # f6: log_2(E)
        lE / math.log(3),            # f7: log_3(E)
        lE / math.log(5),            # f8: log_5(E)
        lE / math.log(6),            # f9: log_6(E)
        lE / math.log(10),           # f10: log_10(E)
        math.log(lE) if lE > 1 else 0,  # f11: log(log(E))
        E % 6,                        # f12: mod 6  (clase modular)
        E % 30,                       # f13: mod 30
    ]
    return np.array(features, dtype=float)

FEATURE_NAMES = [
    "log(E)",
    "log(E)/2",
    "log(E)/3",
    "log(E)²",
    "√log(E)",
    "log₂(E)",
    "log₃(E)",
    "log₅(E)",
    "log₆(E)",
    "log₁₀(E)",
    "log(log(E))",
    "E mod 6",
    "E mod 30",
]

# ──────────────────────────────────────────────
# 3. DETECTOR SVD + BOOTSTRAP
# ──────────────────────────────────────────────

def detectar_mecuation(datos, tau=1e-4, bootstrap_N=300, consensus=0.90, verbose=True):
    """
    Dado un conjunto de (E, metadata), detecta si existe MEcuation local.
    Devuelve dict con resultados.
    """
    valores = [d[0] for d in datos]
    Y = np.array([feature_map(e) for e in valores])
    Y_centrado = Y - Y.mean(axis=0)

    # SVD principal
    U, sigma, Vt = np.linalg.svd(Y_centrado, full_matrices=False)
    kappa = sigma[-1] / sigma[0] if sigma[0] > 0 else 1.0
    var_explicada_r1 = sigma[0]**2 / (sigma**2).sum()

    # Bootstrap
    deg_count = 0
    epsilons = [0.001, 0.005, 0.01]  # perturbaciones relativas
    for _ in range(bootstrap_N):
        eps = np.random.choice(epsilons)
        ruido = 1 + np.random.uniform(-eps, eps, size=len(valores))
        vals_pert = [v * r for v, r in zip(valores, ruido)]
        Yp = np.array([feature_map(e) for e in vals_pert])
        Yp_c = Yp - Yp.mean(axis=0)
        _, sp, _ = np.linalg.svd(Yp_c, full_matrices=False)
        kp = sp[-1] / sp[0] if sp[0] > 0 else 1.0
        if kp < tau:
            deg_count += 1
    deg_frac = deg_count / bootstrap_N

    # Vector principal → MEcuation candidata
    v1 = Vt[0]  # dirección principal

    existe = kappa < tau and deg_frac >= consensus

    resultado = {
        "existe_mecuation": existe,
        "kappa": kappa,
        "var_explicada_r1": var_explicada_r1,
        "sigma": sigma,
        "deg_frac": deg_frac,
        "v1": v1,
        "n_datos": len(datos),
    }

    if verbose:
        print(f"\n{'='*55}")
        print(f"  DETECTOR MEcuation — {len(datos)} observaciones")
        print(f"{'='*55}")
        print(f"  κ (colapso)       : {kappa:.2e}  {'✅ COLAPSO' if kappa < tau else '❌ no colapsa'}")
        print(f"  Var explicada r=1 : {var_explicada_r1*100:.1f}%")
        print(f"  Bootstrap deg_frac: {deg_frac:.2f}  {'✅ ROBUSTO' if deg_frac >= consensus else '❌ frágil'}")
        print(f"  Valores singulares: {sigma[:5].round(4)}")
        print()
        if existe:
            print("  🏆 MEcuation LOCAL DETECTADA")
            print("  Ecuación: Σ cᵢ·fᵢ(E) ≈ 0")
            print("  Coeficientes principales:")
            idx_sorted = np.argsort(np.abs(v1))[::-1]
            for idx in idx_sorted[:5]:
                if abs(v1[idx]) > 0.05:
                    print(f"    [{idx:2d}] {FEATURE_NAMES[idx]:20s}  c = {v1[idx]:+.6f}")
        else:
            print("  ⚪ No se detecta MEcuation global para esta familia.")
            print("  → Buscar subfamilias o ajustar features.")

    return resultado

# ──────────────────────────────────────────────
# 4. INTERPRETADOR DE MEcuation → FACTOR
# ──────────────────────────────────────────────

def interpretar_mecuation(v1, datos, verbose=True):
    """
    Intenta interpretar el vector principal como una relación entre factores.
    Si c[f1] ≈ 2*c[f2]: E ≈ p^2  → factor = √E
    Si c[f1] ≈ 3*c[f3]: E ≈ p^3  → factor = ∛E
    etc.
    """
    c = v1
    interpretaciones = []

    # Detectar ratio f1/f2
    if abs(c[1]) > 1e-6:
        ratio = c[0] / c[1]
        if abs(ratio - 2.0) < 0.15:
            interpretaciones.append(("cuadrado", "E = p²  →  factor = √E = exp(log(E)/2)"))
        if abs(ratio - 3.0) < 0.2:
            interpretaciones.append(("cubo", "E = p³  →  factor = ∛E = exp(log(E)/3)"))

    # Detectar si f2 domina (log(E)/2 grande)
    if abs(c[1]) > 0.5 and abs(c[0]) > 0.5:
        interpretaciones.append(("potencia_par", f"posible potencia par, ratio c[0]/c[1]={c[0]/c[1]:.3f}"))

    # Probar extracción de factor directamente
    factores_extraidos = []
    for E, meta in datos:
        lE = math.log(E)
        # Factor primario estimado por el vector v1
        # Proyección: α = y · v1 (escalar), factor ≈ exp(α * v1[0])
        y = feature_map(E)
        y_c = y - y.mean()
        alpha = float(np.dot(y_c, v1))
        factor_est = math.exp(abs(alpha))
        factores_extraidos.append({
            "E": E,
            "meta": meta,
            "factor_estimado": factor_est,
            "sqrt_E": math.sqrt(E),
            "cbrt_E": E**(1/3),
            "log_E": math.log(E),
        })

    if verbose and interpretaciones:
        print("\n  📐 INTERPRETACIÓN DE LA MEcuation:")
        for nombre, desc in interpretaciones:
            print(f"    → [{nombre}] {desc}")

        print("\n  🔍 MUESTRA DE FACTORES EXTRAÍDOS:")
        print(f"  {'E':>12}  {'factor_real':>14}  {'√E':>10}  {'∛E':>10}")
        print(f"  {'-'*52}")
        for item in factores_extraidos[:8]:
            meta = item["meta"]
            factor_real = meta.get("p", "?")
            print(f"  {item['E']:>12}  {str(factor_real):>14}  "
                  f"{item['sqrt_E']:>10.4f}  {item['cbrt_E']:>10.4f}")

    return interpretaciones, factores_extraidos

# ──────────────────────────────────────────────
# 5. SCAN DE FAMILIAS
# ──────────────────────────────────────────────

def scan_familias(familias=None, n=12):
    """
    Escanea múltiples familias y reporta cuáles tienen MEcuation.
    """
    if familias is None:
        familias = [
            ("cuadrados",          {}),
            ("cubos",              {}),
            ("semiprimos_cercanos",{}),
            ("semiprimos",         {}),
            ("k_primo",            {"k": 2}),
            ("k_primo",            {"k": 6}),
            ("mersenne",           {}),
            ("sophie_germain",     {}),
            ("potencia_mixta",     {}),
        ]

    resumen = []
    print("\n" + "█"*55)
    print("  SCAN DE FAMILIAS — búsqueda de MEcuations")
    print("█"*55)

    for regla, kwargs in familias:
        nombre = regla + (f"(k={kwargs['k']})" if 'k' in kwargs else "")
        datos = generar_familia(regla, n_elementos=n, **kwargs)
        if len(datos) < 3:
            continue

        # Detección silenciosa
        res = detectar_mecuation(datos, verbose=False)
        estado = "✅ MEcuation" if res["existe_mecuation"] else "⚪ sin ME"
        kappa_str = f"{res['kappa']:.1e}"
        deg_str = f"{res['deg_frac']:.2f}"

        print(f"  {nombre:30s}  κ={kappa_str:>10}  boot={deg_str}  {estado}")

        if res["existe_mecuation"]:
            interp, _ = interpretar_mecuation(res["v1"], datos, verbose=False)
            for _, desc in interp:
                print(f"    ↳ {desc}")

        resumen.append((nombre, res, datos))

    return resumen

# ──────────────────────────────────────────────
# 6. MAIN
# ──────────────────────────────────────────────

if __name__ == "__main__":
    np.random.seed(42)

    # ── A) Scan rápido de todas las familias
    resumen = scan_familias(n=12)

    # ── B) Análisis detallado de familias con MEcuation
    print("\n\n" + "="*55)
    print("  ANÁLISIS DETALLADO — familias con MEcuation")
    print("="*55)

    for nombre, res, datos in resumen:
        if res["existe_mecuation"]:
            print(f"\n▶ Familia: {nombre}")
            detectar_mecuation(datos, verbose=True)
            interpretar_mecuation(res["v1"], datos, verbose=True)

    # ── C) Caso especial: cuadrados con análisis completo
    print("\n\n" + "="*55)
    print("  CASO ESPECIAL: cuadrados — Newton Rápido + MEcuation")
    print("="*55)
    datos_cuad = generar_familia("cuadrados", n_elementos=15)
    res_cuad = detectar_mecuation(datos_cuad, verbose=True)
    if res_cuad["existe_mecuation"]:
        interp, factores = interpretar_mecuation(res_cuad["v1"], datos_cuad, verbose=True)
        print("\n  💡 El oráculo para Newton Rápido en cuadrados:")
        print("     j_inicial = log(E)/2  (en lugar de j=1.0)")
        print("     → El factor p = exp(log(E)/2) = √E")
        print("     → MEcuation actúa como guess perfecto: 1 iteración")
        print()
        print("  Verificación directa:")
        for E, meta in datos_cuad[:6]:
            p_real = meta["p"]
            p_est = math.exp(math.log(E)/2)
            error = abs(p_est - p_real)
            print(f"    E={E:6d}  p_real={p_real:4d}  "
                  f"p_estimado={p_est:.4f}  error={error:.2e}")
