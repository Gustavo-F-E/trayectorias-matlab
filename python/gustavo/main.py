#Módulos externos
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Módulos propios
from variables_entrada import obtener_entero_positivo, obtener_flotante_positivo, obtener_angulo_bobinado
from grafico_mandril import grafico_mandril
from grafico_trayectoria import grafico_trayectoria
from primeros_calculos import primeros_calculos

###############################################################################

# Solicitar el ancho de la fibra con validación
ancho_de_la_fibra = obtener_entero_positivo("Por favor, ingrese un valor numérico entero positivo para el ancho de la fibra: ")

angulo_de_bobinado = obtener_angulo_bobinado("Por favor, ingrese un valor numérico entero positivo para el ángulo del bobinado: ")

diametro_del_mandril = obtener_flotante_positivo("Por favor, ingrese un valor numérico positivo para el diametro del mandril: ")

longitud_a_bobinar = obtener_flotante_positivo("Por favor, ingrese un valor numérico positivo para la longitud a bobinar: ")

###############################################################################
#Primeros cálculos
alfa, b_eff, numero_de_divisiones, X_lim, Y_lim, D, L, b, vector_y = primeros_calculos(ancho_de_la_fibra, diametro_del_mandril, longitud_a_bobinar, angulo_de_bobinado)

###############################################################################
#Graficar el mandril
grafico_mandril(D, L)

###############################################################################
# Gráfico en los ejes X e Y

#grafico_trayectoria(y_origen, b, b_eff, alfa, D, L, X_lim,Y_lim, numero_de_divisiones)

for i in range(numero_de_divisiones):
    grafico_trayectoria(i, vector_y, b, b_eff, alfa, D, L, X_lim,Y_lim, numero_de_divisiones, vector_y)