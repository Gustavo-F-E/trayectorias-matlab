import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from openpyxl import Workbook, load_workbook
from write_to_excel import write_to_excel
from update_excel_from_dict import update_excel_from_dict
from variables_de_salida import variables_de_salida
from variables_de_inicio import variables_de_inicio
from variables_de_la_maquina import variables_de_la_maquina
from calculos_iniciales import calculos_iniciales
from generador_superficie import generador_superficie
from funcion_de_Christoffel import funcion_de_Christoffel
from generacion_del_mandril import generacion_del_mandril
from graficar_mandril import graficar_mandril

# Nombre del archivo Excel de salida
#script_directory = os.path.dirname(os.path.abspath(__file__))
#excel = os.path.join(script_directory, 'trayectorias_FW_1.xlsx')
excel = 'trayectorias_FW_python_1.xlsx'  # Cambiar el nombre del archivo de salida

# Crear un diccionario para los tiempos de cálculo
tiempos_de_calculo = {}

# Guardar el tiempo de inicio
tiempos_de_calculo['tiempo_de_inicio'] = datetime.now()

# Opcional: imprimir el tiempo de inicio
print(f"Tiempo de inicio: {tiempos_de_calculo['tiempo_de_inicio']}")

var_sal = variables_de_salida(excel)

#Invocamos a la funcion con los valores de inicio:
var_ini = variables_de_inicio(excel) #variables de seteo del programa

#Invocamos a la funcion con la configuración de la máquina:
var_maq = variables_de_la_maquina(excel)

# Datos específicos del Estudio Nº 1:
var_ini.update({
    'b': 6,
    'rim': 75,
    'rfm': 75,
    'zim': -500,
    'zi': 0,
    'zfg': 300,
    'zfm': 1000,
    'alfainicio': 45,
    'lambda1': 0.14
})

# Llamada a la función para actualizar en el excel los valores ingresados anteriormente
update_excel_from_dict(excel, 'Datos de entrada', var_ini)

# Calculos iniciales:
vector_calculos_iniciales = calculos_iniciales(var_ini, excel)

# Generación del mandril y su gráfico:
vector_mandril = generacion_del_mandril(var_ini, vector_calculos_iniciales, excel)
graficar_mandril(vector_mandril, 'Gráfico del mandril', 1, 'mandril.png')