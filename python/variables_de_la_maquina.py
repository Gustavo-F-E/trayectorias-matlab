import pandas as pd
from openpyxl import Workbook, load_workbook
from write_to_excel import write_to_excel

def variables_de_la_maquina(excel):
    vector_variables_de_la_maquina = {
        'RefA': 0,
        'RefX': 0,
        'RefY': 260,
        'RefB': 0,
        'RefZ': 0,
        'RefC': 0,
        'Xmin': 165.15,  # mm
        'Xmax': 2644.7,  # mm
        'Ymin': -120,    # mm
        'Ymax': 311,     # mm
        'Zmin': 0,       # mm
        'Zmax': 364,     # mm
        'Cmin': -50,     # º
        'Cmax': 50,      # º
        'POew': 50,      # mm
        'POe': 0,        # mm
        'YawR': 260,     # mm
        'VA': 150,       # 360º/min
        'VX': 1,         # m/s
        'VY': 0.5,       # m/s
        'VB': 115,       # 360º/min
        'VZ': 0.5,       # m/s
        'VC': 100,       # 360º/min
        'AA': 8,         # 360º/min
        'AX': 4,         # m/s
        'AY': 30,        # m/s
        'AB': 60,        # 360º/min
        'AZ': 1,         # m/s
        'AC': 45,        # 360º/min
        'RA': 360,       # inc/360º
        'RX': 1,         # inc/mm
        'RY': -1,        # inc/mm
        'RB': 360,       # inc/360º
        'RZ': 1,         # inc/mm
        'RC': 360,       # inc/360º
        'minproc': 0.03
    }

    # Velocidades Máximas - inc/s
    vector_variables_de_la_maquina['VAmax'] = abs(vector_variables_de_la_maquina['VA'] * vector_variables_de_la_maquina['RA'] / 60)
    vector_variables_de_la_maquina['VXmax'] = abs(vector_variables_de_la_maquina['VX'] * 1000 * vector_variables_de_la_maquina['RX'])
    vector_variables_de_la_maquina['VYmax'] = abs(vector_variables_de_la_maquina['VY'] * 1000 * vector_variables_de_la_maquina['RY'])
    vector_variables_de_la_maquina['VBmax'] = abs(vector_variables_de_la_maquina['VB'] * vector_variables_de_la_maquina['RB'] / 60)
    vector_variables_de_la_maquina['VZmax'] = abs(vector_variables_de_la_maquina['VZ'] * 1000 * vector_variables_de_la_maquina['RZ'])
    vector_variables_de_la_maquina['VCmax'] = abs(vector_variables_de_la_maquina['VC'] * vector_variables_de_la_maquina['RC'] / 60)

    # Aceleraciones máximas
    vector_variables_de_la_maquina['AAmax'] = abs(vector_variables_de_la_maquina['AA'] * vector_variables_de_la_maquina['RA'])
    vector_variables_de_la_maquina['AXmax'] = abs(vector_variables_de_la_maquina['AX'] * 60 * 1000 * vector_variables_de_la_maquina['RX'])
    vector_variables_de_la_maquina['AYmax'] = abs(vector_variables_de_la_maquina['AY'] * 60 * 1000 * vector_variables_de_la_maquina['RY'])
    vector_variables_de_la_maquina['ABmax'] = abs(vector_variables_de_la_maquina['AB'] * vector_variables_de_la_maquina['RB'])
    vector_variables_de_la_maquina['AZmax'] = abs(vector_variables_de_la_maquina['AZ'] * 60 * 1000 * vector_variables_de_la_maquina['RZ'])
    vector_variables_de_la_maquina['ACmax'] = abs(vector_variables_de_la_maquina['AC'] * vector_variables_de_la_maquina['RC'])

    # Convertir el diccionario a una lista de listas para facilitar la exportación a Excel
    celda_variables_de_la_maquina = [
        ['RefA', vector_variables_de_la_maquina['RefA']],
        ['RefX', vector_variables_de_la_maquina['RefX']],
        ['RefY', vector_variables_de_la_maquina['RefY']],
        ['RefB', vector_variables_de_la_maquina['RefB']],
        ['RefZ', vector_variables_de_la_maquina['RefZ']],
        ['RefC', vector_variables_de_la_maquina['RefC']],
        ['Xmin', vector_variables_de_la_maquina['Xmin']],
        ['Xmax', vector_variables_de_la_maquina['Xmax']],
        ['Ymin', vector_variables_de_la_maquina['Ymin']],
        ['Ymax', vector_variables_de_la_maquina['Ymax']],
        ['Zmin', vector_variables_de_la_maquina['Zmin']],
        ['Zmax', vector_variables_de_la_maquina['Zmax']],
        ['Cmin', vector_variables_de_la_maquina['Cmin']],
        ['Cmax', vector_variables_de_la_maquina['Cmax']],
        ['POew', vector_variables_de_la_maquina['POew']],
        ['POe', vector_variables_de_la_maquina['POe']],
        ['YawR', vector_variables_de_la_maquina['YawR']],
        ['VA', vector_variables_de_la_maquina['VA']],
        ['VX', vector_variables_de_la_maquina['VX']],
        ['VY', vector_variables_de_la_maquina['VY']],
        ['VB', vector_variables_de_la_maquina['VB']],
        ['VZ', vector_variables_de_la_maquina['VZ']],
        ['VC', vector_variables_de_la_maquina['VC']],
        ['AA', vector_variables_de_la_maquina['AA']],
        ['AX', vector_variables_de_la_maquina['AX']],
        ['AY', vector_variables_de_la_maquina['AY']],
        ['AB', vector_variables_de_la_maquina['AB']],
        ['AZ', vector_variables_de_la_maquina['AZ']],
        ['AC', vector_variables_de_la_maquina['AC']],
        ['RA', vector_variables_de_la_maquina['RA']],
        ['RX', vector_variables_de_la_maquina['RX']],
        ['RY', vector_variables_de_la_maquina['RY']],
        ['RB', vector_variables_de_la_maquina['RB']],
        ['RZ', vector_variables_de_la_maquina['RZ']],
        ['RC', vector_variables_de_la_maquina['RC']],
        ['minproc', vector_variables_de_la_maquina['minproc']],
        ['VAmax', vector_variables_de_la_maquina['VAmax']],
        ['VXmax', vector_variables_de_la_maquina['VXmax']],
        ['VYmax', vector_variables_de_la_maquina['VYmax']],
        ['VBmax', vector_variables_de_la_maquina['VBmax']],
        ['VZmax', vector_variables_de_la_maquina['VZmax']],
        ['VCmax', vector_variables_de_la_maquina['VCmax']],
        ['AAmax', vector_variables_de_la_maquina['AAmax']],
        ['AXmax', vector_variables_de_la_maquina['AXmax']],
        ['AYmax', vector_variables_de_la_maquina['AYmax']],
        ['ABmax', vector_variables_de_la_maquina['ABmax']],
        ['AZmax', vector_variables_de_la_maquina['AZmax']],
        ['ACmax', vector_variables_de_la_maquina['ACmax']]
    ]

    # Convertir la lista de listas a un DataFrame de pandas
    # Escribir el DataFrame a un archivo Excel
    df = pd.DataFrame(celda_variables_de_la_maquina, columns=['Variable', 'Valor'])
    write_to_excel(excel, 'Variables de la maquina', df)

    return vector_variables_de_la_maquina

    print('Se ejecutó la función "variables_de_la_maquina"')