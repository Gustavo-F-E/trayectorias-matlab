import numpy as np
import sympy as sp
import pandas as pd
from generador_superficie import generador_superficie
from funcion_de_Christoffel import funcion_de_Christoffel

# Supongamos que estas funciones ya están definidas
# from tu_modulo import generador_superficie, funcion_de_Christoffel

def generacion_del_mandril(var_ini, vector_calculos_iniciales, excel):
    # Variables de entrada
    dm = vector_calculos_iniciales['dm']
    lm = vector_calculos_iniciales['lm']
    rc = vector_calculos_iniciales['rc']
    rfm = var_ini['rfm']
    rim = var_ini['rim']
    zfm = var_ini['zfm']
    zim = var_ini['zim']
    
    # Inicializamos las variables de salida de la función
    celda_mandril = []
    mandril = {}

    # Únicamente define la curva para un cilindro o un cono
    if rim == rfm:
        X, Y, Z, tc = generador_superficie(rc, rc, zim, zfm)
        # Superficie parametrizada del mandril en coordenadas esféricas
        teta, phi = sp.symbols('teta phi')
        ro = rc / sp.sin(teta)
        Se = [ro * sp.sin(teta) * sp.cos(phi), ro * sp.sin(teta) * sp.sin(phi), ro * sp.cos(teta)]
        
        # Propiedades de la superficie - coordenadas esféricas
        coeficientes_de_Christoffel = funcion_de_Christoffel(Se, phi, teta)
        
        # Escribimos las variables de salida para un cilindro (estructuras y celdas)
        mandril['eteta'] = coeficientes_de_Christoffel['ev']
        mandril['gteta'] = coeficientes_de_Christoffel['gv']
        celda_mandril.append(['eteta', str(mandril['eteta'])])
        celda_mandril.append(['gteta', str(mandril['gteta'])])

    else:
        X, Y, Z, tc = generador_superficie(rim, rfm, zim, zfm)
        # Superficie parametrizada del mandril en coordenadas polares
        phi, ro = sp.symbols('phi ro')
        z = ro * lm / dm - rim * lm / dm + zim
        S = [ro * sp.cos(phi), ro * sp.sin(phi), ro * lm / dm - rim * lm / dm + zim]
        # Propiedades de la superficie - coordenadas polares
        coeficientes_de_Christoffel = funcion_de_Christoffel(S, phi, ro)
        zlinha = sp.diff(z, ro)
        z2linha = sp.diff(zlinha, ro)

        # Escribimos la variable de salida para un cono (estructuras y celdas)
        mandril['ero'] = coeficientes_de_Christoffel['ev']
        mandril['gro'] = coeficientes_de_Christoffel['gv']
        mandril['zlinha'] = zlinha
        mandril['z2linha'] = z2linha
        celda_mandril.append(['ero', str(mandril['ero'])])
        celda_mandril.append(['gro', str(mandril['gro'])])
        celda_mandril.append(['zlinha', str(mandril['zlinha'])])
        celda_mandril.append(['z2linha', str(mandril['z2linha'])])

    mandril['ro'] = ro
    mandril['chr'] = coeficientes_de_Christoffel['chr']
    mandril['E'] = coeficientes_de_Christoffel['E']
    mandril['F'] = coeficientes_de_Christoffel['F']
    mandril['G'] = coeficientes_de_Christoffel['G']
    mandril['gg1'] = coeficientes_de_Christoffel['gg1']
    mandril['L'] = coeficientes_de_Christoffel['L']
    mandril['M'] = coeficientes_de_Christoffel['M']
    mandril['N'] = coeficientes_de_Christoffel['N']
    mandril['gg2'] = coeficientes_de_Christoffel['gg2']
    mandril['e'] = coeficientes_de_Christoffel['e']
    mandril['g'] = coeficientes_de_Christoffel['g']
    mandril['efi'] = coeficientes_de_Christoffel['eu']
    mandril['gfi'] = coeficientes_de_Christoffel['gu']
    mandril['X'] = X
    mandril['Y'] = Y
    mandril['Z'] = Z
    mandril['tc'] = tc
    
    # Escribir una celda, ya que es más fácil exportar un archivo con formato de celdas
    celda_mandril.append(['ro', str(mandril['ro'])])
    celda_mandril.append(['E', str(mandril['E'])])
    celda_mandril.append(['F', str(mandril['F'])])
    celda_mandril.append(['G', str(mandril['G'])])
    celda_mandril.append(['L', str(mandril['L'])])
    celda_mandril.append(['M', str(mandril['M'])])
    celda_mandril.append(['N', str(mandril['N'])])
    celda_mandril.append(['efi', str(mandril['efi'])])
    celda_mandril.append(['gfi', str(mandril['gfi'])])
    celda_mandril.append(['tc', mandril['tc']])
    celda_mandril.append(['gg1', str(mandril['gg1'])])
    celda_mandril.append(['gg2', str(mandril['gg2'])])
    celda_mandril.append(['e', str(mandril['e'])])
    celda_mandril.append(['g', str(mandril['g'])])

    # Crear un DataFrame de pandas y escribirlo en Excel
    df = pd.DataFrame(celda_mandril, columns=['Variable', 'Valor'])
    with pd.ExcelWriter(excel, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Variables del mandril', startrow=0, index=False)
        pd.DataFrame({'X': mandril['X'].flatten()}).to_excel(writer, sheet_name='Variables del mandril', startrow=20, startcol=1, index=False)
        pd.DataFrame({'Y': mandril['Y'].flatten()}).to_excel(writer, sheet_name='Variables del mandril', startrow=22, startcol=1, index=False)
        pd.DataFrame({'Z': mandril['Z'].flatten()}).to_excel(writer, sheet_name='Variables del mandril', startrow=24, startcol=1, index=False)
        pd.DataFrame({'chr1': [str(mandril['chr'][:, :, 0])]}).to_excel(writer, sheet_name='Variables del mandril', startrow=26, startcol=1, index=False)
        pd.DataFrame({'chr2': [str(mandril['chr'][:, :, 1])]}).to_excel(writer, sheet_name='Variables del mandril', startrow=27, startcol=1, index=False)

    return mandril
