import numpy as np
import matplotlib.pyplot as plt

def grafico_mandril(D, L):
    #Primer gráfico: el mandril
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Crear los datos para el cilindro
    theta = np.linspace(0, 2*np.pi, 100)
    x = np.linspace(0, L, 100)
    theta, x = np.meshgrid(theta, x)
    r = D / 2
    y = r * np.cos(theta)
    z = r * np.sin(theta)

    # Graficar el cilindro
    ax.plot_surface(x, y, z, color='cyan', alpha=0.8)

    # Configurar los límites de los ejes
    ax.set_xlim(0, L)
    ax.set_ylim(-r, r)
    ax.set_zlim(-r, r)

    # Configurar etiquetas de los ejes
    ax.set_xlabel('X (Longitud)')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Añadir título y leyenda
    plt.title('Cilindro (Mandril)')
    ax.text2D(0.05, 0.95, f'Diámetro (D) = {D}\nLongitud (L) = {L}', transform=ax.transAxes)

    # Mostrar el gráfico
    plt.show()