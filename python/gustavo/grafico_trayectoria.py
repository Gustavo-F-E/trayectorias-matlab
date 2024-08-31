import numpy as np
import matplotlib.pyplot as plt

def evaluar_xy(y_origen, X_lim, Y_lim, b, alfa):
    b_eff_y = b / np.cos(alfa)
    b_eff_x = b / np.sin(alfa)
    n = 0
    x_ini = 0
    paso = X_lim / 1000
    x_values = []
    y_values = []
    contador = 0

    while x_ini < X_lim:
        contador = contador + 1
        #print(f"\nPaso Nº {contador}")
        x = x_ini
        y = x * np.tan(alfa) - n * Y_lim + y_origen
        if y <= Y_lim and x < X_lim:
            #print(f"Caso 1: y <= Y_lim and x < X_lim")
            x_ini = x_ini + paso
            xx = x
            yy = y
            x_values.append(xx)
            y_values.append(yy)
            #print(f"{xx} {yy} {n}")
        elif y <= Y_lim and x >= X_lim:
            x_ini = X_lim
            xx = x_ini
            yy = x * np.tan(alfa) - n * Y_lim + y_origen
            x_values.append(xx)
            y_values.append(yy)
            #print(f"Caso 2: {xx} {yy} {n} final")
        elif y > Y_lim and x < X_lim:
            n = n +1
            x_ini = (n) * Y_lim / np.tan(alfa)
            xx = x_ini
            yy = xx * np.tan(alfa) - n * Y_lim + y_origen
            x_values.append(xx)
            y_values.append(yy)
            x_ini = x_ini + paso
            #print(f"Caso 3: y > Y_lim and x < X_lim")
            #print(f"{xx} {yy} {n}")
        else:
            print(f"Caso 4: problemas")
    #print(f"{x_values} {y_values}")
    xx = X_lim
    yy = X_lim * np.tan(alfa) - n * Y_lim + y_origen
    x_values.append(xx)
    y_values.append(yy)
    return (x_values, y_values)

def grafico_trayectoria(i, y_origen, b, b_eff, alfa, D, L, X_lim,Y_lim, numero_de_divisiones, vector_y):
    x_values, y_values = evaluar_xy(y_origen[i], X_lim, Y_lim, b, alfa)
    x = x_values

    # Calcular Y para la línea principal
    y = y_values

    # Crear la figura y los ejes
    fig, ax = plt.subplots(figsize=(10, 6))

    # Graficar la línea principal
    ax.plot(x, y, 'b-', linewidth = 0.1, label='Línea principal')

    # Calcular y graficar la línea paralela para representar el ancho
    y_paralela = y + b_eff

    ax.plot(x, y_paralela, 'b--', linewidth=0.1, label='Línea paralela')

    # Rellenar el área entre las líneas
    ax.fill_between(x, y, y_paralela, alpha=0.5, color='blue', label='Área entre líneas')

    # Configurar los límites de los ejes
    ax.set_xlim(0, X_lim)
    ax.set_ylim(0, Y_lim)

    # Configurar etiquetas de los ejes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Añadir título y leyenda
    plt.title(f'Línea con ancho y ángulo específicos, trayectoria {i}')
    ax.legend()

    # Añadir texto con información sobre los parámetros
    info_text = f'Ángulo (alfa) = {np.degrees(alfa)}°\nAncho efectivo (b_eff) = {b_eff}\nDiámetro = {D}\nLongitud = {L}\nLongitud del cilindro extendida = {Y_lim}\nNúmero_de_divisiones = {numero_de_divisiones}'

    ax.text(0.05, 0.95, info_text, transform=ax.transAxes, verticalalignment='top')

    # Mostrar el gráfico

    plt.grid(True)

    plt.show()