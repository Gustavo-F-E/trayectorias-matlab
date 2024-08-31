import math
import pandas as pd
import numpy as np
from openpyxl import load_workbook

def calculos_iniciales(var_ini, excel_path):
    # Definición de variables de entrada
    alfainicio = var_ini['alfainicio']
    b = var_ini['b']
    fibdens = var_ini['fibdens']
    fmc = var_ini['fmc']
    resdens = var_ini['resdens']
    rfm = var_ini['rfm']
    rim = var_ini['rim']
    RW = var_ini['RW']
    TEX = var_ini['TEX']
    zfg = var_ini['zfg']
    zfm = var_ini['zfm']
    zi = var_ini['zi']
    zig = var_ini['zig']
    zim = var_ini['zim']
    
    # Inicializamos las variables de salida de la función
    vector_salida_calculos_iniciales = {}
    celda_salida_calculos_iniciales = []

    # Cálculos iniciales
    alfai = alfainicio * math.pi / 180
    lm = zfm - zim
    lg = zfg - zig
    dm = rfm - rim
    rig = ((zig * dm / lm) + rim - (dm * zim / lm))
    ri = ((zi * dm / lm) + rim - (dm * zim / lm))
    rfg = ((zfg * dm / lm) + rim - (dm * zim / lm))
    tau = math.atan(dm / lm)
    zeroteta = (zfg + zig) / 2
    rc = rim

    if rig == rfg:
        beff = b / math.cos(alfai)
        nreal = 2 * math.pi * rc / beff
        n1 = math.ceil(nreal)
        vector_salida_calculos_iniciales['beff'] = beff
        vector_salida_calculos_iniciales['nreal'] = nreal
        vector_salida_calculos_iniciales['n1'] = n1
        celda_salida_calculos_iniciales.append(['beff', vector_salida_calculos_iniciales['beff']])
        celda_salida_calculos_iniciales.append(['nreal', vector_salida_calculos_iniciales['nreal']])
        celda_salida_calculos_iniciales.append(['n1', vector_salida_calculos_iniciales['n1']])

    # Espesor
    fd = fmc * fibdens + (1 - fmc) * resdens
    esp = TEX / (1000 * RW * fibdens)

    # Vector de variables de salida
    vector_salida_calculos_iniciales['alfai'] = alfai
    vector_salida_calculos_iniciales['dm'] = dm
    vector_salida_calculos_iniciales['esp'] = esp
    vector_salida_calculos_iniciales['fd'] = fd
    vector_salida_calculos_iniciales['lg'] = lg
    vector_salida_calculos_iniciales['lm'] = lm
    vector_salida_calculos_iniciales['rc'] = rc
    vector_salida_calculos_iniciales['rfg'] = rfg
    vector_salida_calculos_iniciales['ri'] = ri
    vector_salida_calculos_iniciales['rig'] = rig
    vector_salida_calculos_iniciales['tau'] = tau
    vector_salida_calculos_iniciales['zeroteta'] = zeroteta

    # Agregar a celda_salida_calculos_iniciales
    celda_salida_calculos_iniciales.append(['alfai', vector_salida_calculos_iniciales['alfai']])
    celda_salida_calculos_iniciales.append(['dm', vector_salida_calculos_iniciales['dm']])
    celda_salida_calculos_iniciales.append(['esp', vector_salida_calculos_iniciales['esp']])
    celda_salida_calculos_iniciales.append(['fd', vector_salida_calculos_iniciales['fd']])
    celda_salida_calculos_iniciales.append(['lg', vector_salida_calculos_iniciales['lg']])
    celda_salida_calculos_iniciales.append(['lm', vector_salida_calculos_iniciales['lm']])
    celda_salida_calculos_iniciales.append(['rc', vector_salida_calculos_iniciales['rc']])
    celda_salida_calculos_iniciales.append(['rfg', vector_salida_calculos_iniciales['rfg']])
    celda_salida_calculos_iniciales.append(['ri', vector_salida_calculos_iniciales['ri']])
    celda_salida_calculos_iniciales.append(['rig', vector_salida_calculos_iniciales['rig']])
    celda_salida_calculos_iniciales.append(['tau', vector_salida_calculos_iniciales['tau']])
    celda_salida_calculos_iniciales.append(['zeroteta', vector_salida_calculos_iniciales['zeroteta']])

    # Convertir a DataFrame
    df_celda_salida_calculos_iniciales = pd.DataFrame(celda_salida_calculos_iniciales, columns=['Variable', 'Valor'])

    # Escribir al archivo Excel
    try:
        with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a') as writer:
            df_celda_salida_calculos_iniciales.to_excel(writer, sheet_name='Cálculos iniciales', index=False)
        print(f'Se escribió la hoja "Cálculos iniciales" en el archivo Excel.')
    except Exception as e:
        print(f'Error al escribir en el archivo Excel: {e}')

    return vector_salida_calculos_iniciales