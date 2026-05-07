"""
metodo_restos.py
================
Método Gráfico de Restos — VMA (Víctor Manzanares Alberola)

Concepto original:
    Dado un entero x y un parámetro a < x, se define
        b(a, x) = x mod (x - a)
    Al incrementar a en 1, b sigue secuencias de pendiente creciente
    (1, 2, 3, ...) separadas por "caídas" a 1.  Los primos presentan
    una estructura regular en estas secuencias que permite inferir
    primalidad con pocos módulos.

Estructura del módulo
---------------------
1. calcular_restos(x)       → lista de pares (a, b) para a=1..x-2
2. detectar_secuencias(x)   → segmenta los restos en tramos de pendiente constante
3. checker_primo(x)         → heurística VMA basada en coherencia de secuencias
4. tabla_restos(limite)     → imprime tabla n × divisores (como en el libro)
5. visualizar_restos(x)     → visualizador Tkinter + Matplotlib (interactivo)

Autoría
-------
Investigación original: Víctor Manzanares Alberola
Asistencia de escritura/código: IA
"""

import math


# ─────────────────────────────────────────────────────────────
# 1. FUNCIÓN BASE
# ─────────────────────────────────────────────────────────────

def calcular_restos(x: int) -> list[tuple[int, int]]:
    """
    Devuelve [(a, b)] donde b = x mod (x - a), para a = 1 ... x-2.

    Puntualidades (VMA):
      · Cuando x - a = (x+1)/2  →  b = (x-1)/2
      · Cuando x - a = (x+1)/2 - 1  →  b = 1
      · Desde a=1 hasta x - (x+1)//2, b incrementa en 1 con a
        (primera secuencia, propia de los impares).
    """
    restos = []
    for a in range(1, x - 1):
        divisor = x - a
        if divisor > 0:
            restos.append((a, x % divisor))
    return restos


# ─────────────────────────────────────────────────────────────
# 2. DETECCIÓN DE SECUENCIAS DE PENDIENTE CONSTANTE
# ─────────────────────────────────────────────────────────────

def detectar_secuencias(x: int) -> list[dict]:
    """
    Segmenta los restos en tramos donde Δb = constante.
    Devuelve lista de {'pendiente': m, 'tramo': [(a,b), ...]}
    """
    datos = calcular_restos(x)
    if len(datos) < 2:
        return []

    secuencias = []
    tramo_actual = [datos[0]]
    pendiente_actual = None

    for i in range(1, len(datos)):
        a_prev, b_prev = datos[i - 1]
        a_curr, b_curr = datos[i]
        delta = b_curr - b_prev

        if pendiente_actual is None:
            pendiente_actual = delta
            tramo_actual.append(datos[i])
        elif delta == pendiente_actual:
            tramo_actual.append(datos[i])
        else:
            secuencias.append({'pendiente': pendiente_actual, 'tramo': tramo_actual})
            tramo_actual = [datos[i - 1], datos[i]]
            pendiente_actual = delta

    secuencias.append({'pendiente': pendiente_actual, 'tramo': tramo_actual})
    return secuencias


def contar_secuencias(x: int) -> int:
    """Número de tramos de pendiente constante para x."""
    return len(detectar_secuencias(x))


# ─────────────────────────────────────────────────────────────
# 3. CHECKER DE PRIMALIDAD POR COHERENCIA DE SECUENCIAS (VMA)
# ─────────────────────────────────────────────────────────────

def checker_primo(n: int) -> bool:
    """
    Heurística VMA para primalidad basada en coherencia de restos.

    Algoritmo (del libro):
        Para i desde 2 hasta (n-1)//2:
            j = n % i
            x = (n - j) // i       # entrada de f(y) = n mod (n-x)
            Si n % x ≠ j % x → retorna False   (nunca debería ocurrir en primo)
            Si n % x == 0    → retorna False   (x es divisor)
        Retorna True

    Nota: es una condición necesaria derivada de la estructura de
    secuencias, no una prueba de primalidad estándar.  Complementar
    con test de divisibilidad para certeza.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False

    for i in range(2, (n - 1) // 2 + 1):
        j = n % i
        if i == 0:
            continue
        x = (n - j) // i
        if x == 0:
            continue
        if n % x == 0:
            return False
        if (n % x) != (j % x):
            return False   # coherencia de secuencia rota → compuesto
    return True


def verificar_checker(limite: int = 200) -> dict:
    """
    Compara checker_primo con isprime de sympy hasta 'limite'.
    Devuelve estadísticas de aciertos/falsos.
    """
    try:
        from sympy import isprime
    except ImportError:
        print("sympy no disponible; usando criba básica para referencia.")
        def isprime(k):
            if k < 2: return False
            for d in range(2, int(k**0.5) + 1):
                if k % d == 0: return False
            return True

    resultados = {'verdaderos_positivos': 0, 'verdaderos_negativos': 0,
                  'falsos_positivos': 0, 'falsos_negativos': 0}

    for k in range(2, limite + 1):
        pred = checker_primo(k)
        real = isprime(k)
        if pred and real:
            resultados['verdaderos_positivos'] += 1
        elif not pred and not real:
            resultados['verdaderos_negativos'] += 1
        elif pred and not real:
            resultados['falsos_positivos'] += 1
        else:
            resultados['falsos_negativos'] += 1

    total = sum(resultados.values())
    aciertos = resultados['verdaderos_positivos'] + resultados['verdaderos_negativos']
    resultados['precision'] = aciertos / total if total else 0
    return resultados


# ─────────────────────────────────────────────────────────────
# 4. TABLA DE RESTOS n × DIVISORES
# ─────────────────────────────────────────────────────────────

def tabla_restos(n_min: int = 2, n_max: int = 31,
                 div_min: int = 2, div_max: int = 17) -> None:
    """
    Imprime la tabla de n mod d, como en el prefacio de 'Números otra vez'.
    Filas: n desde n_min hasta n_max.
    Columnas: divisores desde div_min hasta div_max.
    """
    encabezado = f"{'n':>4} |" + "".join(f"{d:>4}" for d in range(div_min, div_max + 1))
    print(encabezado)
    print("-" * len(encabezado))
    for n in range(n_min, n_max + 1):
        fila = f"{n:>4} |"
        for d in range(div_min, div_max + 1):
            if d > n:
                fila += f"{d:>4}"   # marca el divisor igual a n (diagonal)
            else:
                fila += f"{n % d:>4}"
        print(fila)


# ─────────────────────────────────────────────────────────────
# 5. VISUALIZADOR TKINTER + MATPLOTLIB
# ─────────────────────────────────────────────────────────────

def visualizar_restos(n_inicial: int = 13) -> None:
    """
    Visualizador interactivo de restos y pendientes (VMA).

    Dibuja:
      · Puntos (i, n mod i) para i=2..n-1
      · Hipotenusa del cuadrado de lado n
      · Líneas rojas de referencia
      · Rectas de decremento/incremento extendidas hasta la hipotenusa
      · Línea imaginaria desde (1,0) con pendiente 1
      · Botones +/- para cambiar n
      · Zoom con rueda del ratón
    """
    try:
        import tkinter as tk
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    except ImportError as e:
        print(f"Dependencias no disponibles: {e}")
        print("Instala: pip install matplotlib")
        return

    def calcular(n):
        x_vals, y_vals = [], []
        for i in range(2, n):
            x_vals.append(i)
            y_vals.append(n % i)
        return x_vals, y_vals

    def encontrar_corte(xi, yi, pendiente, n):
        intercepto = yi - pendiente * xi
        if abs(pendiente) < 1e-10:
            return (n, yi)
        corte_x = (n - intercepto) / pendiente
        corte_y = intercepto + pendiente * corte_x
        return corte_x, corte_y

    def trazar(n, ax):
        x, y = calcular(n)
        ax.clear()
        ax.plot(x, y, 'bo-', markersize=4, label='Restos n mod i')
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
        ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
        # Hipotenusa
        ax.plot([0, n], [n, 0], 'g-', linewidth=1.2, label='Hipotenusa')
        # Marcas rojas
        for xi, yi in zip(x, y):
            ax.plot([xi, xi], [0, yi], 'r--', linewidth=0.5, alpha=0.4)
            ax.plot([0, xi], [yi, yi], 'r--', linewidth=0.5, alpha=0.4)
            ax.text(xi, yi, f'{xi},{yi}', fontsize=6, ha='right', alpha=0.7)
        # Rectas de decremento/incremento extendidas
        for i in range(1, len(y)):
            if y[i] != y[i - 1]:
                pendiente = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
                cx, cy = encontrar_corte(x[i], y[i], pendiente, n)
                ax.plot([x[i], cx], [y[i], cy], 'b-', linewidth=0.8, alpha=0.5)
        # Línea imaginaria pendiente 1
        ax.plot([1, n], [0, n - 1], 'b:', linewidth=1, label='Pendiente 1 desde (1,0)')
        ax.set_xlim(0, n)
        ax.set_ylim(0, n)
        ax.set_xlabel('Divisores (i)')
        ax.set_ylabel('Restos (n mod i)')
        ax.set_title(f'Trazado de Restos y Pendientes  —  n = {n}  (VMA)')
        ax.legend(fontsize=7)
        canvas.draw()

    root = tk.Tk()
    root.title("Método Gráfico de Restos — VMA")

    fig, ax = plt.subplots(figsize=(8, 6))
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    n_var = tk.IntVar(value=n_inicial)
    n = [n_inicial]

    def actualizar(delta):
        n[0] = max(5, n[0] + delta)
        n_var.set(n[0])
        trazar(n[0], ax)

    def zoom(event):
        base = 1.1
        xlim = ax.get_xlim(); ylim = ax.get_ylim()
        xd, yd = event.xdata or 0, event.ydata or 0
        sf = 1 / base if event.button == 'up' else base
        nw = (xlim[1] - xlim[0]) * sf; nh = (ylim[1] - ylim[0]) * sf
        rx = (xlim[1] - xd) / (xlim[1] - xlim[0])
        ry = (ylim[1] - yd) / (ylim[1] - ylim[0])
        ax.set_xlim([xd - nw * (1 - rx), xd + nw * rx])
        ax.set_ylim([yd - nh * (1 - ry), yd + nh * ry])
        canvas.draw()

    canvas.mpl_connect('scroll_event', zoom)

    frame = tk.Frame(root)
    frame.pack()
    tk.Label(frame, text="n =", font=("Arial", 12)).pack(side=tk.LEFT)
    tk.Label(frame, textvariable=n_var, font=("Arial", 14, "bold")).pack(side=tk.LEFT, padx=5)
    tk.Button(frame, text="−", font=("Arial", 14),
              command=lambda: actualizar(-1)).pack(side=tk.LEFT, padx=5)
    tk.Button(frame, text="+", font=("Arial", 14),
              command=lambda: actualizar(+1)).pack(side=tk.LEFT, padx=5)
    tk.Button(frame, text="+10", font=("Arial", 12),
              command=lambda: actualizar(+10)).pack(side=tk.LEFT, padx=5)

    trazar(n[0], ax)
    root.mainloop()


# ─────────────────────────────────────────────────────────────
# 6. DEMO EN CONSOLA
# ─────────────────────────────────────────────────────────────

def demo():
    print("=" * 60)
    print("MÉTODO GRÁFICO DE RESTOS — VMA")
    print("=" * 60)

    # Tabla de restos (prefacio de 'Números otra vez')
    print("\n--- Tabla de restos (n mod d) ---")
    tabla_restos(n_min=2, n_max=20, div_min=2, div_max=10)

    # Secuencias para varios valores
    print("\n--- Secuencias de pendiente constante ---")
    for x in [7, 11, 13, 17, 25, 29, 31]:
        seqs = detectar_secuencias(x)
        pendientes = [s['pendiente'] for s in seqs]
        print(f"  x={x:3d}  →  {len(seqs)} secuencias, pendientes: {pendientes}")

    # Checker de primalidad
    print("\n--- Checker VMA de primalidad ---")
    candidatos = list(range(2, 50))
    marcados = [n for n in candidatos if checker_primo(n)]
    print(f"  Checker marca como primo: {marcados}")

    stats = verificar_checker(200)
    print(f"\n  Estadísticas hasta n=200:")
    print(f"    Verdaderos positivos : {stats['verdaderos_positivos']}")
    print(f"    Verdaderos negativos : {stats['verdaderos_negativos']}")
    print(f"    Falsos positivos     : {stats['falsos_positivos']}")
    print(f"    Falsos negativos     : {stats['falsos_negativos']}")
    print(f"    Precisión            : {stats['precision']:.2%}")

    # Coherencia de secuencias: primos vs compuestos
    print("\n--- Estructura de secuencias: primos vs compuestos ---")
    print(f"  {'n':>5}  {'#sec':>5}  {'primo?':>8}")
    try:
        from sympy import isprime
    except ImportError:
        def isprime(k):
            if k < 2: return False
            for d in range(2, int(k**0.5)+1):
                if k % d == 0: return False
            return True

    for n in range(5, 40):
        ns = contar_secuencias(n)
        print(f"  {n:>5}  {ns:>5}  {'SÍ' if isprime(n) else 'no':>8}")

    print("\nPara el visualizador gráfico interactivo ejecuta: visualizar_restos(13)")


if __name__ == "__main__":
    demo()
