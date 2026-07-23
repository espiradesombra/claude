
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np

def plot_restos(n):
    i_vals = np.arange(1, n)
    y_vals = n % i_vals
    fig, ax = plt.subplots()
    ax.plot(i_vals, y_vals, label=f"n mod i")
    ax.plot([0, n], [0, n], 'k--', label="Hipotenusa")
    for slope in range(1, 6):
        ax.plot(i_vals, slope * i_vals, '--', label=f"Pendiente {slope}")
    ax.set_title(f"Restos de n={n}")
    ax.legend()
    return fig

def launch_gui(n):
    root = tk.Tk()
    root.title("Visualización de restos")
    fig = plot_restos(n)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack()
    root.mainloop()

if __name__ == "__main__":
    launch_gui(97)
