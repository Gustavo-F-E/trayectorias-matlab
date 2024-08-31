import numpy as np

def generador_superficie(a, b, c, d):
    tc = np.arange(0, 2 * np.pi + np.pi / 20, np.pi / 20)
    X = np.zeros((2, len(tc)))
    Y = np.zeros((2, len(tc)))
    Z = np.zeros((2, len(tc)))

    for k, t in enumerate(tc):
        X[0, k] = a * np.cos(t)
        X[1, k] = b * np.cos(t)
        Y[0, k] = a * np.sin(t)
        Y[1, k] = b * np.sin(t)
        Z[0, k] = c
        Z[1, k] = d

    return X, Y, Z, tc