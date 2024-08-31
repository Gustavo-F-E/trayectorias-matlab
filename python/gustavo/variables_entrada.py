#Parametros de entrada del usuario
def obtener_entero_positivo(mensaje):
    while True:
        try:
            valor = int(input(mensaje))
            if valor > 0:
                return valor
            else:
                print("Error: Por favor, ingrese un número entero positivo.")
        except ValueError:
            print("Error: Por favor, ingrese un número entero válido.")

def obtener_flotante_positivo(mensaje):
    while True:
        try:
            valor = float(input(mensaje))
            if valor > 0:
                return valor
            else:
                print("Error: Por favor, ingrese un número positivo.")
        except ValueError:
            print("Error: Por favor, ingrese un número válido.")

def obtener_angulo_bobinado(mensaje):
    while True:
        try:
            valor = int(input(mensaje))
            if valor <= 0:
                print("Error: Por favor, ingrese un número entero positivo.")
            elif valor < 5:
                print("No es posible realizar bobinados de menos de 5º")
            elif valor > 89:
                print("No es posible realizar bobinados mayores a 89º")
            else:
                return valor
        except ValueError:
            print("Error: Por favor, ingrese un número entero válido.")