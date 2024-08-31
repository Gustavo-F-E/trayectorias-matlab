import numpy as np

def primeros_calculos(ancho_de_la_fibra, diametro_del_mandril, longitud_a_bobinar, angulo_de_bobinado):
    b = ancho_de_la_fibra
    D = diametro_del_mandril
    L = longitud_a_bobinar
    X_lim = L

    calculo_1 = 2 * np.pi * diametro_del_mandril / (ancho_de_la_fibra / np.cos(np.radians(angulo_de_bobinado)))
    calculo_2 = 2 * np.pi * diametro_del_mandril // (ancho_de_la_fibra / np.cos(np.radians(angulo_de_bobinado)))
    calculo_3 = calculo_1 - calculo_2

    # Comprobación condicional
    if 0 < calculo_3 < 0.5:
        calculo_4 = (calculo_2 * ancho_de_la_fibra) / (2 * np.pi * diametro_del_mandril)
        alfa = np.arccos(calculo_4)
        print("Se ha aplicado el caso donde 0 < calculo_3 < 0.5")
    elif calculo_3 >= 0.5:
        calculo_2 = np.ceil(calculo_1)
        calculo_4 = (calculo_2 * ancho_de_la_fibra) / (2 * np.pi * diametro_del_mandril)
        alfa = np.arccos(calculo_4)
        print("Se ha aplicado el caso donde calculo_3 >= 0.5")
    elif calculo_3 == 0:
        alfa = np.radians(angulo_de_bobinado)  # Convertimos a radianes
        print("Se ha aplicado el caso donde calculo_3 == 0")
    else:
        print("No se ha aplicado ningún caso")
    b_eff = b / np.cos(alfa)

    print(f"Ancho de la fibra: {ancho_de_la_fibra}")
    print(f"Diámetro del mandril: {diametro_del_mandril}")
    print(f"Longitud a bobinar: {longitud_a_bobinar}")
    print(f"calculo_1: {calculo_1}")
    print(f"calculo_2: {calculo_2}")
    print(f"calculo_3: {calculo_3}")
    print(f"calculo_4: {calculo_4}")
    print(f"alfa: {np.degrees(alfa)} grados")
    print(f"b_eff: {b_eff}")

    Y_lim = 2 * np.pi * D
    numero_de_divisiones = int(Y_lim // b_eff)
    paso_en_y = Y_lim / numero_de_divisiones
    print(f"Paso en Y: {paso_en_y}")
    vector_y = np.arange(0, Y_lim, paso_en_y)
    print(f"Vector Y: {vector_y}")
    print(f"Número de divisiones: {numero_de_divisiones}")

    return alfa, b_eff, numero_de_divisiones, X_lim, Y_lim, D, L, b, vector_y