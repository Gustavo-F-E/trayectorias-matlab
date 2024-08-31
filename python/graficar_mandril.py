import numpy as np
import matplotlib.pyplot as plt

def graficar_mandril(a, b, factor, nombre):
    # Generar la hélice usando la función generador_helice
    X_grafico = a['X']
    Y_grafico = a['Y']
    Z_grafico = a['Z']

    x_min = np.min(a['Z']) / factor
    x_max = np.max(a['Z']) / factor

    # Graficar la hélice
    plt.figure(figsize=(15, 10))
    plt.plot(Z_grafico, Y_grafico, X_grafico, linewidth=2)
    plt.grid(True)
    plt.xlabel('X', color='k', fontsize=14)
    plt.ylabel('Y', color='k', fontsize=14)
    plt.gca().set_zlabel('Z', color='k', fontsize=14)  # Se necesita Axes3D para el etiquetado del eje Z
    plt.title(b, color='k', fontsize=20)
    plt.gca().set_aspect('equal')
    plt.gcf().set_facecolor((0.92, 0.96, 1))

    # Cambiar el color de los ejes a negro
    ax = plt.gca()
    ax.xaxis.label.set_color('k')
    ax.yaxis.label.set_color('k')

    # Cambiar el color de los números de los ejes a negro
    ax.tick_params(axis='x', colors='k')
    ax.tick_params(axis='y', colors='k')

    # Establecer límites para el eje X
    plt.xlim([x_min, x_max])
    plt.savefig(nombre)
    plt.show()