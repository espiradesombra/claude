import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk

# Función para calcular los restos
def calcular_restos(n):
    x = []
    y = []
    for i in range(2, n):
        x.append(i)
        y.append(n % i)
    return x, y

# Función para encontrar el corte con la hipotenusa
def encontrar_corte(xi, yi, pendiente, n):
    intercepto = yi - pendiente * xi
    corte_x = (n - intercepto) / pendiente
    corte_y = intercepto + pendiente * corte_x
    return corte_x, corte_y

# Función para trazar los restos y pendientes
def trazar_pendientes(n, ax):
    x, y = calcular_restos(n)
    
    ax.clear()
    ax.plot(x, y, 'bo-', label='Restos')
    ax.axhline(y=0, color='black', linestyle='--', label='Eje x')
    ax.axvline(x=0, color='black', linestyle='--', label='Eje y')

    # Dibujar hipotenusa con catetos n
    ax.plot([0, n], [n, 0], 'g-', label='Hipotenusa')

    # Añadir marcas y líneas rojas para n
    for xi, yi in zip(x, y):
        ax.plot([xi, xi], [0, yi], 'r--')
        ax.plot([0, xi], [yi, yi], 'r--')
        ax.text(xi, yi, f'{xi},{yi}', fontsize=8, ha='right')

    # Dibujar las rectas en los decrementos extendidas hasta la hipotenusa
    for i in range(1, len(y)):
        if y[i] < y[i-1]:  # Hay un decremento
            pendiente = (y[i] - y[i-1]) / (x[i] - x[i-1])
            corte_x, corte_y = encontrar_corte(x[i], y[i], pendiente, n)
            ax.plot([x[i], corte_x], [y[i], corte_y], 'b-', label=f'Decremento desde ({x[i]},{y[i]})')
    
    # Dibujar una línea desde la hipotenusa a la recta de pendiente 2
    for i in range(1, len(y)-1):
        if y[i] < y[i-1] and y[i+1] < y[i]:  # Decremento seguido de incremento
            pendiente_hipotenusa = -1
            pendiente_2 = (y[i+1] - y[i]) / (x[i+1] - x[i])
            corte_x, corte_y = encontrar_corte(x[i], y[i], pendiente_hipotenusa, n)
            ax.plot([x[i], corte_x], [y[i], corte_y], 'm-', linestyle=':', label=f'Cruce hipotenusa a pendiente 2 desde ({x[i]},{y[i]})')
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], 'c-', label=f'Decremento a Incremento desde ({x[i]},{y[i]})')
        if y[i] > y[i-1] :  # Incremento seguido de Decremento
            pendiente_hipotenusa = -1
            pendiente_2 = -y[i]+y[i-1]
            corte_x, corte_y = encontrar_corte(x[i], y[i], pendiente_hipotenusa, n)
            ax.plot([x[i-1], corte_x], [y[i-1], corte_y], 'm-', linestyle=':', label=f'Cruce hipotenusa a pendiente 2 desde ({x[i-1]},{y[i-1]})')
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], 'c-', label=f'Incremento a Decremento desde ({x[i]},{y[i]})')


    # Línea imaginaria desde (1,0) con pendiente 1
    ax.plot([1, n], [0, n-1], 'b-', linestyle=':', label='Línea Imaginaria')

    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_xlabel('Divisores (i)')
    ax.set_ylabel('Restos (n % i)')
    ax.set_title(f'Trazado de Restos y Pendientes para n = {n}')
    ax.legend()
    canvas.draw()

# Función para actualizar el valor de n
def actualizar_n(incremento):
    global n
    n += incremento
    n_var.set(n)
    trazar_pendientes(n, ax)

# Configuración de la interfaz gráfica
root = tk.Tk()
root.title("Trazado Dinámico de Restos y Pendientes")

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Variable global para n
n = 13
n_var = tk.IntVar(value=n)

# Cuadro de texto para mostrar el valor de n
n_label = tk.Label(root, textvariable=n_var, font=("Arial", 14))
n_label.pack()

# Botones para aumentar y disminuir el valor de n
btn_mas = tk.Button(root, text="+", command=lambda: actualizar_n(1), font=("Arial", 14))
btn_mas.pack(side=tk.LEFT, padx=10, pady=10)

btn_menos = tk.Button(root, text="-", command=lambda: actualizar_n(-1), font=("Arial", 14))
btn_menos.pack(side=tk.RIGHT, padx=10, pady=10)

# Trazar el gráfico inicial
trazar_pendientes(n, ax)

root.mainloop()
