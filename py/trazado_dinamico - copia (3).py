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

# Función para trazar los restos y pendientes
def trazar_pendientes(n, ax):
    x, y = calcular_restos(n)
    
    ax.clear()
    ax.plot(x, y, 'bo-', label='Restos')
    ax.axhline(y=0, color='black', linestyle='--', label='Eje x')
    ax.axvline(x=0, color='black', linestyle='--', label='Eje y')

    # Añadir marcas y líneas rojas para n
    for xi, yi in zip(x, y):
        ax.plot([xi, xi], [0, yi], 'r--')
        ax.plot([0, xi], [yi, yi], 'r--')
        ax.text(xi, yi, f'{xi},{yi}', fontsize=8, ha='right')

    # Dibujar hipotenusa con catetos n
    ax.plot([0, n], [n, 0], 'g-', label='Hipotenusa')

    # Alargar las rectas de las secuencias de los restos
    for xi, yi in zip(x, y):
        if yi != 0:
            ax.plot([xi, (n + xi - yi)], [yi, yi], 'r--')  # Línea horizontal hasta la hipotenusa
            ax.plot([xi, xi], [yi, n], 'r--')  # Línea vertical

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
