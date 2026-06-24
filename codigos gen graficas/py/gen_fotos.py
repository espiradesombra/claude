import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
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

# Nueva función trazar_visor sin ax.clear()
def trazar_visor(ax, n):
    x, y = calcular_restos(n)
    
    # Limpiar el eje al inicio
    ax.clear()
    
    ax.plot(x, y, 'bo-', label='Restos')
    ax.axhline(y=0, color='black', linestyle='--')
    ax.axvline(x=0, color='black', linestyle='--')
    
    # Dibujar hipotenusa
    ax.plot([0, n], [n, 0], 'g-', label='Hipotenusa')
    
    # Añadir marcas y líneas
    for xi, yi in zip(x, y):
        ax.plot([xi, xi], [0, yi], 'r--')
        ax.plot([0, xi], [yi, yi], 'r--')
        ax.text(xi, yi, f'({xi},{yi})', fontsize=8, ha='right')
    
        # Dibujar las rectas en los decrementos extendidas hasta la hipotenusa
    for i in range(1, len(y)):
        if y[i] < y[i-1]:  # Hay un decremento
            pendiente = (y[i] - y[i-1]) / (x[i] - x[i-1])
            corte_x, corte_y = encontrar_corte(x[i], y[i], pendiente, n)
            ax.plot([x[i], corte_x], [y[i], corte_y], 'b-', label=f'Decremento desde ({x[i]},{y[i]})')
            # Actualizar límites
            min_x = min(min_x, corte_x)
            max_x = max(max_x, corte_x)
            min_y = min(min_y, corte_y)
            max_y = max(max_y, corte_y)

    # Dibujar las rectas en los aumentos extendidas hasta la hipotenusa
    for i in range(1, len(y)):
        if y[i] > y[i-1]:  # Hay un aumento
            pendiente = (y[i] - y[i-1]) / (x[i] - x[i-1])
            corte_x, corte_y = encontrar_corte(x[i], y[i], pendiente, n)
            ax.plot([x[i], corte_x], [y[i], corte_y], 'm-', label=f'Aumento desde ({x[i]},{y[i]})')
            # Actualizar límites
            min_x = min(min_x, corte_x)
            max_x = max(max_x, corte_x)
            min_y = min(min_y, corte_y)
            max_y = max(max_y, corte_y)
    # Línea imaginaria
    ax.plot([1, n], [0, n-1], 'b:', label='Línea Imaginaria')
    
    # Ajustar límites
    ax.set_xlim(0, n + 10)
    ax.set_ylim(0, n + 10)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('Divisores (i)')
    ax.set_ylabel('Restos (n % i)')
    ax.set_title(f'Trazado para n = {n}')
    ax.legend()

# Función para iniciar el visor dinámico con zoom independiente
def iniciar_visor():
    def graficar():
        try:
            # Obtener números n1 y n2 desde las entradas
            n1 = int(entry_n1.get())
            n2 = int(entry_n2.get()) if entry_n2.get() else 1
        except ValueError:
            print("Introduce valores válidos.")
            return

        # Graficar en la primera figura
        trazar_pendientes(n1,ax1)
        canvas1.draw()

        # Graficar en la segunda figura si n2 está definido
       
        trazar_pendientes(n2,ax2)
        canvas2.draw()
        

    # Crear una ventana separada para el visor dinámico
    visor_root = tk.Toplevel()
    visor_root.title("Visor Dinámico de Restos y Pendientes")

    # Frame para los parámetros de entrada
    frame_parametros = tk.Frame(visor_root)
    frame_parametros.pack(pady=10)

    tk.Label(frame_parametros, text="Número 1:", font=("Arial", 12)).grid(row=0, column=0, padx=5)
    entry_n1 = tk.Entry(frame_parametros, width=10, font=("Arial", 12))
    entry_n1.grid(row=0, column=1, padx=5)

    tk.Label(frame_parametros, text="Número 2 (opcional):", font=("Arial", 12)).grid(row=1, column=0, padx=5)
    entry_n2 = tk.Entry(frame_parametros, width=10, font=("Arial", 12))
    entry_n2.grid(row=1, column=1, padx=5)

    btn_graficar = tk.Button(frame_parametros, text="Graficar", font=("Arial", 12), command=graficar)
    btn_graficar.grid(row=2, column=0, columnspan=2, pady=10)

    # Crear frames para las gráficas
    frame_graficas = tk.Frame(visor_root)
    frame_graficas.pack()

    # Gráfica 1
    frame_grafica1 = tk.Frame(frame_graficas)
    frame_grafica1.pack(side=tk.LEFT, padx=5)

    fig1, ax1 = plt.subplots(figsize=(6, 5))  # Figura y eje para gráfica 1
    canvas1 = FigureCanvasTkAgg(fig1, master=frame_grafica1)
    canvas1.draw()
    canvas1.get_tk_widget().pack()

    toolbar1 = NavigationToolbar2Tk(canvas1, frame_grafica1)
    toolbar1.update()
    canvas1.get_tk_widget().pack()

    # Gráfica 2
    frame_grafica2 = tk.Frame(frame_graficas)
    frame_grafica2.pack(side=tk.RIGHT, padx=5)

    fig2, ax2 = plt.subplots(figsize=(6, 5))  # Figura y eje para gráfica 2
    canvas2 = FigureCanvasTkAgg(fig2, master=frame_grafica2)
    canvas2.draw()
    canvas2.get_tk_widget().pack()

    toolbar2 = NavigationToolbar2Tk(canvas2, frame_grafica2)
    toolbar2.update()
    canvas2.get_tk_widget().pack()

# ... Aquí seguiría el resto de su código, como las funciones para el generador de gráficas y la configuración de la ventana principal.

# Función para guardar las gráficas según parámetros

# Función para guardar las gráficas según parámetros
# Generador principal de gráficas
def generar_graficas():
    try:
        cuantas_graficas = int(entry_cuantas_graficas.get() or 1)
        inicio = int(entry_inicio.get() or 10)
        fin = int(entry_fin.get() or 50)
        salto = int(entry_salto.get() or 1)
    except ValueError:
        print("Introduce valores numéricos válidos.")
        return

    save_dir = filedialog.askdirectory(title="Selecciona el directorio donde guardar las gráficas")
    if not save_dir:
        print("No se seleccionó directorio.")
        return

    fig, ax = plt.subplots()
    for n_actual in range(inicio, fin + 1, salto):
        for _ in range(cuantas_graficas):
            trazar_pendientes(n_actual, ax)
            file_path = f"{save_dir}/grafica_n_{n_actual}.png"
            fig.savefig(file_path)
            print(f"Gráfica guardada: {file_path}")

    plt.close(fig)
    print("Todas las gráficas han sido generadas.")
# Función trazar_pendientes original para el generador de gráficas

def trazar_pendientes(n, ax):
    x, y = calcular_restos(n)
    
    ax.clear()
    ax.plot(x, y, 'bo-', label='Restos')
    ax.axhline(y=0, color='black', linestyle='--', label='Eje x')
    ax.axvline(x=0, color='black', linestyle='--', label='Eje y')

    # Dibujar hipotenusa con catetos n
    ax.plot([0, n], [n, 0], 'g-', label='Hipotenusa')

    # Variables para calcular los límites de los ejes
    min_x, max_x = 0, n
    min_y, max_y = 0, n

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
            # Actualizar límites
            min_x = min(min_x, corte_x)
            max_x = max(max_x, corte_x)
            min_y = min(min_y, corte_y)
            max_y = max(max_y, corte_y)

    # Dibujar las rectas en los aumentos extendidas hasta la hipotenusa
    for i in range(1, len(y)):
        if y[i] > y[i-1]:  # Hay un aumento
            pendiente = (y[i] - y[i-1]) / (x[i] - x[i-1])
            corte_x, corte_y = encontrar_corte(x[i], y[i], pendiente, n)
            ax.plot([x[i], corte_x], [y[i], corte_y], 'm-', label=f'Aumento desde ({x[i]},{y[i]})')
            # Actualizar límites
            min_x = min(min_x, corte_x)
            max_x = max(max_x, corte_x)
            min_y = min(min_y, corte_y)
            max_y = max(max_y, corte_y)

    # Línea imaginaria desde (1,0) con pendiente 1
    ax.plot([1, n], [0, n-1], 'b-', linestyle=':', label='Línea Imaginaria')

    # Ajustar los límites para incluir todos los cruces
    ax.set_xlim(min_x - 1, max_x + 1)
    ax.set_ylim(min_y - 1, max_y + 1)
    ax.set_aspect('equal', adjustable='box')  # Mantener proporción cuadrada

    # Configurar la leyenda fuera de la gráfica
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=10)

    # Etiquetas y título
    ax.set_xlabel('Divisores (i)')
    ax.set_ylabel('Restos (n % i)')
    ax.set_title(f'Trazado de Restos y Pendientes para n = {n}')

# Configuración de la interfaz gráfica principal
root = tk.Tk()
root.title("Trazado Dinámico de Restos y Pendientes")

# Frame para los parámetros de entrada
parametros_frame = tk.Frame(root)
parametros_frame.pack(side=tk.TOP, padx=10, pady=10)

# Campo para "Cuántas gráficas"
tk.Label(parametros_frame, text="Cuántas gráficas:", font=("Arial", 12)).grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
entry_cuantas_graficas = tk.Entry(parametros_frame, width=10, font=("Arial", 12))
entry_cuantas_graficas.grid(row=0, column=1, padx=5, pady=5)

# Campo para "Número inicial"
tk.Label(parametros_frame, text="Número inicial:", font=("Arial", 12)).grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
entry_inicio = tk.Entry(parametros_frame, width=10, font=("Arial", 12))
entry_inicio.grid(row=1, column=1, padx=5, pady=5)

# Campo para "Número final"
tk.Label(parametros_frame, text="Número final:", font=("Arial", 12)).grid(row=2, column=0, sticky=tk.W, padx=5, pady=5)
entry_fin = tk.Entry(parametros_frame, width=10, font=("Arial", 12))
entry_fin.grid(row=2, column=1, padx=5, pady=5)

# Campo para "Salto"
tk.Label(parametros_frame, text="Salto:", font=("Arial", 12)).grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)
entry_salto = tk.Entry(parametros_frame, width=10, font=("Arial", 12))
entry_salto.grid(row=3, column=1, padx=5, pady=5)

# Checkbutton para "Solo 6k ± 1" (si aún lo utilizas)
# var_seis_k = tk.IntVar()
# check_seis_k = tk.Checkbutton(parametros_frame, text="Solo números de la forma 6k ± 1", variable=var_seis_k, font=("Arial", 12))
# check_seis_k.grid(row=4, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

# Botón para generar las gráficas
btn_generar_graficas = tk.Button(root, text="Generar Gráficas", font=("Arial", 14), command=generar_graficas)
btn_generar_graficas.pack(side=tk.BOTTOM, pady=10)

# Botón para abrir el visor dinámico
btn_visor = tk.Button(root, text="Abrir Visor Dinámico", font=("Arial", 14), command=iniciar_visor)
btn_visor.pack(pady=10)

root.mainloop()
