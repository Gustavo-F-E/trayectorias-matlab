import sympy as sp

def funcion_de_Christoffel(S, u, v):
    # S = Ecuación de la forma de la superficie
    # u = Primera variable independiente
    # v = Segunda variable independiente
    u, v = sp.symbols('u v')
    phi, teta = sp.symbols('phi teta')
    Se = sp.Matrix([75*sp.cos(phi), 75*sp.sin(phi), 75*sp.cos(teta)/sp.sin(teta)])


    # Primera fórmula fundamental
        # Verificar si S es una función o una expresión
    if isinstance(S, sp.Matrix):
        e = sp.Matrix([sp.diff(S[i], u) for i in range(S.shape[0])])
        g = sp.Matrix([sp.diff(S[i], v) for i in range(S.shape[0])])
    else:
        e = sp.diff(S, u)
        g = sp.diff(S, v)

    E1 = e.dot(e)
    E = sp.simplify(e.dot(e))
    g = sp.diff(S, v)
    G1 = g.dot(g)
    G = sp.simplify(g.dot(g))
    F1 = g.dot(e)
    F = sp.simplify(e.dot(g))

    gg1 = sp.Matrix([[E, F], [F, G]])

    n = sp.simplify(sp.cross(e, g) / sp.norm(sp.cross(e, g)))
    # Segunda fórmula fundamental
    L = sp.simplify(n.dot(sp.diff(e, u)))
    M = sp.simplify(n.dot(sp.diff(e, v)))
    N = sp.simplify(n.dot(sp.diff(g, v)))

    gg2 = sp.Matrix([[L, M], [M, N]])

    # Derivadas
    ev = sp.simplify(sp.diff(E, v))
    eu = sp.simplify(sp.diff(E, u))
    gv = sp.simplify(sp.diff(G, v))
    gu = sp.simplify(sp.diff(G, u))

    # Símbolos de Christoffel
    chr = sp.MutableDenseNDimArray.zeros(2, 2, 2)
    chr[0, 0, 0] = sp.simplify((G * sp.diff(E, u) - 2 * F * sp.diff(F, u) + F * sp.diff(E, v)) / (2 * (E * G - F**2)))
    chr[0, 0, 1] = sp.simplify((2 * E * sp.diff(F, u) - E * sp.diff(E, v) - F * sp.diff(E, u)) / (2 * (E * G - F**2)))
    chr[0, 1, 0] = sp.simplify((G * sp.diff(E, v) - F * sp.diff(G, u)) / (2 * (E * G - F**2)))
    chr[1, 0, 0] = chr[0, 1, 0]
    chr[0, 1, 1] = sp.simplify((E * sp.diff(G, u) - F * sp.diff(E, v)) / (2 * (E * G - F**2)))
    chr[1, 0, 1] = chr[0, 1, 1]
    chr[1, 1, 0] = sp.simplify((2 * G * sp.diff(F, v) - G * sp.diff(G, u) - F * sp.diff(G, v)) / (2 * (E * G - F**2)))
    chr[1, 1, 1] = sp.simplify((E * sp.diff(G, v) - 2 * F * sp.diff(F, v) - F * sp.diff(G, u)) / (2 * (E * G - F**2)))

    coeficientes_de_Christoffel = {
        'chr': chr,
        'E': E,
        'F': F,
        'G': G,
        'gg1': gg1,
        'L': L,
        'M': M,
        'N': N,
        'gg2': gg2,
        'e': e,
        'g': g,
        'ev': ev,
        'eu': eu,
        'gv': gv,
        'gu': gu
    }

    return coeficientes_de_Christoffel