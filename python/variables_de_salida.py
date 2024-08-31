import pandas as pd
from openpyxl import Workbook, load_workbook
from write_to_excel import write_to_excel

def variables_de_salida(excel):
    vector_variables_de_salida = {
        'nomeoutput': '04.03.2016_rc75_b10_-500,0,300,1000,alfa20',
        'filepath': r'C:\Tesis de maestria\Archivos 0\0-Filament winding\fundamentos de patrones\nuevo_trayectoria_filament_winding',
        'gravartrajectorias': 's',
        'cnc': 's',
        'gravarcnc': 's',
        'D': 50,
        'seccoes': 'rectangular',
        'belip': 1,
        'sobre': 'sobreposiçao2'
    }

    # Convertir el diccionario a una lista de listas para facilitar la exportación a Excel
    celda_variables_de_salida = [
        ['nomeoutput', vector_variables_de_salida['nomeoutput']],
        ['filepath', vector_variables_de_salida['filepath']],
        ['gravartrajectorias', vector_variables_de_salida['gravartrajectorias']],
        ['cnc', vector_variables_de_salida['cnc']],
        ['gravarcnc', vector_variables_de_salida['gravarcnc']],
        ['D', vector_variables_de_salida['D']],
        ['seccoes', vector_variables_de_salida['seccoes']],
        ['belip', vector_variables_de_salida['belip']],
        ['sobre', vector_variables_de_salida['sobre']]
    ]
    # Convertir la lista de listas a un DataFrame de pandas
    # Escribir el DataFrame a un archivo Excel
    df = pd.DataFrame(celda_variables_de_salida, columns=['Variable', 'Valor'])
    write_to_excel(excel, 'Salida de archivos', df)

    return vector_variables_de_salida

    print('Se ejecutó la función "variables_de_salida"')
