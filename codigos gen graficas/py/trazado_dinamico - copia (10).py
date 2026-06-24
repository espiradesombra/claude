import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import filedialog  # Para guardar la imagen

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

    # Línea imaginaria desde (1,0) con pendiente 1
    ax.plot([1, n], [0, n-1], 'b-', linestyle=':', label='Línea Imaginaria')

    # Asegurarse de que la gráfica sea siempre un cuadrado
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_aspect('equal', adjustable='box')  # Mantener proporción cuadrada
    ax.set_xlabel('Divisores (i)')
    ax.set_ylabel('Restos (n % i)')
    ax.set_title(f'Trazado de Restos y Pendientes para n = {n}')
    ax.legend()

# Función para guardar las gráficas según parámetros
def generar_graficas(cuantas_graficas, inicio, fin, salto, solo_seis_k):
    # Crear un directorio para las gráficas guardadas
    save_dir = filedialog.askdirectory(title="Selecciona el directorio donde guardar las gráficas")
    if not save_dir:
        print("No se seleccionó directorio. Operación cancelada.")
        return

    fig, ax = plt.subplots()
    n_actual = inicio
    
    while n_actual <= fin:
        # Generar la cantidad de gráficas por número
        for i in range(cuantas_graficas):
            trazar_pendientes(n_actual, ax)
            file_path = f"{save_dir}/grafica_n_{n_actual}_iter_{i+1}.png"
            fig.savefig(file_path)
            print(f"Gráfica guardada: {file_path}")
        
        # Avanzar al siguiente número según el tipo de salto
        if solo_seis_k:
            # Considerar solo números de la forma 6k ± 1
            if (n_actual - 1) % 6 == 0:
                n_actual += 4  # Pasar al siguiente 6k + 1
            elif (n_actual + 1) % 6 == 0:
                n_actual += 2  # Pasar al siguiente 6k - 1
            else:
                n_actual += 1  # Asegurarse de avanzar correctamente
        else:
            n_actual += salto

    plt.close(fig)
    print("Todas las gráficas han sido generadas y guardadas.")

# Configuración de la interfaz gráfica
root = tk.Tk()
root.title("Trazado Dinámico de Restos y Pendientes")

# Configurar el área principal (Frame) para colocar elementos
main_frame = tk.Frame(root)
main_frame.pack(side=tk.TOP, padx=10, pady=10)

# Botón para generar las gráficas automáticamente
btn_generar_graficas = tk.Button(main_frame, text="Generar Gráficas Automáticas", font=("Arial", 14), 
                                  command=lambda: generar_graficas(
                                      cuantas_graficas=2,    # Cuántas gráficas por número
                                      inicio=10,             # Número inicial
                                      fin=100,              # Número final
                                      salto=5,               # Salto entre números
                                      solo_seis_k=True       # Activar 6k ± 1
                                  ))
btn_generar_graficas.pack(side=tk.TOP, pady=5)

root.mainloop()
