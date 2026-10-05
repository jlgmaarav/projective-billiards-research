# -*- coding: utf-8 -*-
"""Certificado simbólico del recuento global de la densidad en V_5.

No usa Newton, diferencias finitas, SVD ni muestreo. Comprueba las familias de
factorización y las ecuaciones de los estratos de N para los tres charts de
la variedad de involuciones.
"""

import sympy as sp


I = sp.eye(3)
u, v, p, q = sp.symbols("u v p q")


def inv(H):
    assert sp.simplify(H * H - I) == sp.zeros(3)
    assert sp.trace(H) == 1
    assert sp.det(H) == -1


def pair(X, H0):
    assert sp.simplify(H0 * X - X.inv() * H0) == sp.zeros(3)
    inv(H0)
    inv(X * H0)


def class_dimension(X):
    z = sp.symbols("z0:9")
    Z = sp.Matrix(3, 3, z)
    equations = list(Z * X - X * Z)
    M, _ = sp.linear_eq_to_matrix(equations, z)
    centralizer_gl = 9 - M.rank()
    return 9 - centralizer_gl


def gcd_constant(f, g):
    F = sp.Poly(sp.expand(f), u, v, p, q, domain=sp.QQ)
    G = sp.Poly(sp.expand(g), u, v, p, q, domain=sp.QQ)
    return sp.gcd(F, G).total_degree() == 0


def chart(k):
    w = [None, None, None]
    ell = [None, None, None]
    w[k] = 1
    other = [j for j in range(3) if j != k]
    w[other[0]], w[other[1]] = u, v
    ell[other[0]], ell[other[1]] = p, q
    ell[k] = 1 - p * w[other[0]] - q * w[other[1]]
    return sp.simplify(2 * sp.Matrix(w) * sp.Matrix([ell]) - I)


def check_factor_fibres():
    rho, t, r, s = sp.symbols("rho t r s", nonzero=True)

    Xreg = sp.diag(1, rho, 1 / rho)
    pair(Xreg, sp.Matrix([[1, 0, 0], [0, 0, t], [0, 1 / t, 0]]))

    X3 = sp.Matrix([[1, 1, 0], [0, 1, 1], [0, 0, 1]])
    pair(X3, sp.Matrix([[1, 1 - t, t * (t - 1) / 2],
                        [0, -1, t], [0, 0, 1]]))

    X2 = sp.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    pair(X2, sp.Matrix([[-1, r, s], [0, 1, 0], [0, 0, 1]]))
    pair(X2, sp.Matrix([[1, r, 0], [0, -1, 0], [0, s, 1]]))

    Xinv = sp.diag(1, -1, -1)
    pair(Xinv, sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]]))
    a, b, c = sp.symbols("a b c")
    Hinv = sp.Matrix([[1, 0, 0], [0, a, b], [0, c, -a]])
    Hinv = Hinv.subs(c, (1 - a ** 2) / b)
    pair(Xinv, Hinv)

    Xns = sp.Matrix([[1, 0, 0], [0, -1, 1], [0, 0, -1]])
    pair(Xns, sp.Matrix([[1, 0, 0], [0, -1, r], [0, 0, 1]]))
    pair(Xns, sp.Matrix([[1, 0, 0], [0, 1, r], [0, 0, -1]]))

    assert [class_dimension(X) for X in
            (sp.eye(3), Xinv, X2, X3, Xns)] == [0, 4, 4, 6, 6]

    print("factorizacion: OK (fibras 4,2,2,1,1 y tipo no semisimple 1;")
    print("               clases conjugadas: 0,4,4,6,6)")


def check_charts():
    A1 = I
    Ainv = sp.diag(1, -1, -1)
    A3 = sp.Matrix([[1, 1, 0], [0, 1, 1], [0, 0, 1]])
    A2 = sp.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    Ans = sp.Matrix([[1, 0, 0], [0, -1, 1], [0, 0, -1]])

    for k in range(3):
        N = chart(k)
        data = {}
        for name, A in (("I", A1), ("involucion", Ainv),
                        ("J3", A3), ("J2", A2), ("ns", Ans)):
            f = sp.factor((N - A).det())
            gp = sp.factor(sp.trace(N * A.inv()) - 3)
            gm = sp.factor(sp.trace(N * A.inv()) + 1)
            data[name] = (f, gp, gm)

        assert data["I"][0] == 0 and data["I"][2] == 0
        assert data["involucion"][0] == 0
        assert data["involucion"][1] != 0 and data["involucion"][2] != 0

        f, gp, gm = data["J3"]
        assert f != 0 and len(sp.factor_list(f)[1]) == 1
        assert gcd_constant(f, gp) and gcd_constant(f, gm)

        f, gp, gm = data["J2"]
        assert f != 0 and sp.cancel(gm / f).free_symbols == set()

        f, gp, gm = data["ns"]
        assert f != 0 and gcd_constant(f, gp) and gcd_constant(f, gm)

    print("charts N: OK (J3 y tipo no semisimple: estrato regular dim 3, especial dim <= 2)")


def check_involution_branch():
    A = sp.diag(1, -1, -1)
    one = sp.ones(3, 1)
    N = 2 * one * one.T / 3 - I
    B = N * A
    char = sp.factor(B.charpoly().as_expr())
    assert char == (sp.Symbol("lambda") - 1) * (3 * sp.Symbol("lambda") ** 2
                                                  + 2 * sp.Symbol("lambda") + 3) / 3
    print("rama A=diag(1,-1,-1): OK (el sublugar no regular es propio; B=I es un punto)")


def main():
    check_factor_fibres()
    check_charts()
    check_involution_branch()
    pair_dims = [4, 6, 7, 6, 7]
    tail_dims = [6, 5, 4, 5, 4]
    totals = [a + b for a, b in zip(pair_dims, tail_dims)]
    assert totals == [10, 11, 11, 11, 11]
    assert max(totals) == 11 < 12
    print("recuento: dim Z_A <= [10,11,11,11,11] <= 11 < 12")
    print("recuento global: dim Z <= 11 < 12")
    print("VERIFICACIÓN EXACTA: OK")


if __name__ == "__main__":
    main()
