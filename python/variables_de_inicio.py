import pandas as pd
from openpyxl import Workbook, load_workbook
from write_to_excel import write_to_excel

def variables_de_inicio(excel):
    vector_variables_de_inicio = {
        'alfainicio': 20,  # alfainicio: ángulo de bobinado en la coordenada inicial de bobinado (zi) en grados α ϵ [-90; 90].
        'b': 10,           # b: ancho de hilo de fibra en mm.
        'ds': 2,           # ds: incremento para el cálculo a lo largo de los caminos geodésicos en mm.
        'dsta': 2,         # dsta: incremento para el cálculo a lo largo de las rutas de retorno en mm.
        'fibdens': 1.77,   # fibdens: densidad de la fibra en g/cm3.
        'fmc': 60,         # fmc: contenido en masa de fibra.
        'fvc': 0,          # fvc: contenido en volumen de fibra.
        'niu': 0.14,       # niu: Coeficiente de fricción estática µ.
        'NP': 6,           # NP: Número de caminos por hilo (Number of rovings). El número real de caminos será 4NP.
        'NR': 10,          # NR: Número de mechas (Number of rovings).
        'resdens': 1.2,    # resdens: densidad de la resina en g/cm3.
        'rfm': 15,         # rfm: Radio del mandril en el extremo derecho en mm.
        'rim': 15,         # rim: Radio del mandril en el extremo izquierdo en mm.
        'TEX': 800,        # TEX: valor TEX en g/km
        'vel': 1.5,        # vel: Velocidad de la máquina (multiplica el tiempo que tarda cada línea de comando en realizar el movimiento de los ejes).
        'zfg': 300,        # zfg: Coordenada del extremo geodésico derecho del mandril en mm.
        'zfm': 1000,       # zfm: Coordenada del extremo derecho del mandril en mm.
        'zi': 0,           # zi: Coordenada del inicio del devanado en mm.
        'zig': 0,          # zig: Coordenada del extremo geodésico izquierdo del mandril en mm.
        'zim': -500,       # zim: Coordenada del extremo izquierdo del mandril en mm.
    }
    vector_variables_de_inicio['lambda1'] = vector_variables_de_inicio['niu']                           # lambda1: parámetro λ donde λ = cµ y c ϵ [0; 1].
    vector_variables_de_inicio['RW'] = vector_variables_de_inicio['b'] / vector_variables_de_inicio['NR']  # RW: b/NR - ancho de mecha en mm.

    # Convertir el diccionario a una lista de listas para facilitar la exportación a Excel
    celda_variables_de_inicio = [
        ['alfainicio', vector_variables_de_inicio['alfainicio']],
        ['b', vector_variables_de_inicio['b']],
        ['ds', vector_variables_de_inicio['ds']],
        ['dsta', vector_variables_de_inicio['dsta']],
        ['fibdens', vector_variables_de_inicio['fibdens']],
        ['fmc', vector_variables_de_inicio['fmc']],
        ['fvc', vector_variables_de_inicio['fvc']],
        ['niu', vector_variables_de_inicio['niu']],
        ['NP', vector_variables_de_inicio['NP']],
        ['NR', vector_variables_de_inicio['NR']],
        ['resdens', vector_variables_de_inicio['resdens']],
        ['rfm', vector_variables_de_inicio['rfm']],
        ['rim', vector_variables_de_inicio['rim']],
        ['TEX', vector_variables_de_inicio['TEX']],
        ['vel', vector_variables_de_inicio['vel']],
        ['zfg', vector_variables_de_inicio['zfg']],
        ['zfm', vector_variables_de_inicio['zfm']],
        ['zi', vector_variables_de_inicio['zi']],
        ['zig', vector_variables_de_inicio['zig']],
        ['zim', vector_variables_de_inicio['zim']],
        ['lambda1', vector_variables_de_inicio['lambda1']],
        ['RW', vector_variables_de_inicio['RW']]
    ]
    # Convertir la lista de listas a un DataFrame de pandas
    # Escribir el DataFrame a un archivo Excel
    df = pd.DataFrame(celda_variables_de_inicio, columns=['Variable', 'Valor'])
    write_to_excel(excel, 'Datos de entrada', df)

    return vector_variables_de_inicio

    print('Se ejecutó la función "variables_de_inicio"')