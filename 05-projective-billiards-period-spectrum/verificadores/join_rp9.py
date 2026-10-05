# -*- coding: utf-8 -*-
"""Certificado exacto de la construcción join (n,k)=(9,14).

El factor es el billar B_5 de ``reflectivo_rp4_b5.py``.  Dos copias se
sincronizan sobre una misma recta proyectiva mediante la homografía racional

    A = [[-200, 40000], [-1, 100]].

Todas las decisiones (acotación, primeros impactos, cierre y orientación) se
hacen sobre Q.  El proceso termina con código 0 sólo si todo pasa.
"""

import sys
from fractions import Fraction as F


FAIL = []


def claim(label, ok, extra=""):
    print(f"[{'OK  ' if ok else 'FAIL'}] {label}" + (f"  {extra}" if extra else ""))
    if not ok:
        FAIL.append(label)


def mm(a, b):
    return [[sum(a[i][r] * b[r][j] for r in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def mv(a, v):
    return [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a))]


def eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def inverse(a):
    n = len(a)
    b = [row[:] + eye(n)[i] for i, row in enumerate(a)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if b[i][col])
        b[col], b[pivot] = b[pivot], b[col]
        z = b[col][col]
        b[col] = [x / z for x in b[col]]
        for i in range(n):
            if i != col:
                z = b[i][col]
                b[i] = [b[i][j] - z * b[col][j] for j in range(2 * n)]
    return [row[n:] for row in b]


def reflection(normal):
    d = sum(F(x * x) for x in normal)
    return [[F(int(i == j)) - F(2 * normal[i] * normal[j], d)
             for j in range(5)] for i in range(5)]


def block_diag(a, b):
    na, nb = len(a), len(b)
    out = [[F(0) for _ in range(na + nb)] for _ in range(na + nb)]
    for i in range(na):
        for j in range(na):
            out[i][j] = a[i][j]
    for i in range(nb):
        for j in range(nb):
            out[na + i][na + j] = b[i][j]
    return out


ROOTS = [
    (0, 1, 0, 0, -1), (0, 0, 1, 0, -1), (0, 0, 0, 1, 0),
    (1, -1, 0, 0, 0), (0, 1, 1, 0, 0), (1, 0, 0, 0, 1),
    (1, 0, -1, 0, 0),
]
EPS = (-1, -1, 1, -1, -1, -1, -1)
R = [reflection(n) for n in ROOTS]

# X=T(x,y,z,w,1), igual que en el certificado del factor.
T = [
    [F(-1, 3), F(-1, 3), F(1, 3), F(1, 3), F(1, 3)],
    [F(1), 0, 0, 0, 0], [0, F(1), 0, 0, 0],
    [0, 0, F(1), 0, 0], [0, 0, 0, F(1), 0],
]
TI = inverse(T)
ELL = TI[4]
L = [[F(EPS[i] * ROOTS[i][j]) for j in range(5)] for i in range(7)]


def dot(covector, vector):
    return sum(covector[i] * vector[i] for i in range(len(vector)))


def factor_lift(point, tangent=False):
    return mv(T, point + [F(0 if tangent else 1)])


def fline(e, x0, u0):
    return [e[0] * x0[i] + e[1] * u0[i] for i in range(5)]


print("C1. Dos factores B_5 y monodromía")
p0 = [F(0), F(0), F(-4, 7), F(-2, 21)]
u0 = [F(1, 5), F(1), F(-1, 4), F(-1, 5)]
x0, v0 = factor_lift(p0), factor_lift(u0, tangent=True)

factor_product = eye(5)
for i in [5, 6, 0, 1, 2, 3, 4]:
    factor_product = mm(R[i], factor_product)
claim("producto de cada factor = -I5", factor_product == [[-x for x in row] for row in eye(5)])

H = [block_diag(r, eye(5)) for r in R] + [block_diag(eye(5), r) for r in R]
total = eye(10)
itinerary = [5, 6, 0, 1, 2, 3, 4, 12, 13, 7, 8, 9, 10, 11]
for i in itinerary:
    total = mm(H[i], total)
claim("producto de las catorce reflexiones = -I10",
      total == [[-x for x in row] for row in eye(10)])

print("\nC2. Politopo join: no vacío, acotado y con catorce facetas")
facets = [row + [F(0)] * 5 for row in L] + [[F(0)] * 5 + row for row in L]
ell_join = ELL + ELL
xint_factor = factor_lift([F(1, 10), F(1, 10), F(-1, 10), F(0)])
xint = xint_factor + xint_factor
xint = [z / dot(ell_join, xint) for z in xint]
claim("punto interior racional", all(dot(row, xint) < 0 for row in facets))
claim("funcional de corte positivo en el punto interior", dot(ell_join, xint) == 1)

# En cada bloque, la dependencia positiva de los siete normales espaciales
# certificada en el factor hace compacto el corte; el funcional ell_1+ell_2
# da el corte compacto del cono producto.
spatial = [[sum(L[i][a] * T[a][j] for a in range(5)) for j in range(4)]
           for i in range(7)]
claim("dependencia positiva del primer factor",
      all(sum(spatial[i][j] for i in range(7)) == 0 for j in range(4)))
claim("dependencia positiva del segundo factor",
      all(sum(spatial[i][j] for i in range(7)) == 0 for j in range(4)))

print("\nC3. Homografía racional de sincronización")
A = [[F(-200), F(40000)], [F(-1), F(100)]]
det_a = A[0][0] * A[1][1] - A[0][1] * A[1][0]
claim("det(A)>0", det_a > 0, f"det={det_a}")
e_start = [F(1), F(1, 100)]
ae = mv(A, e_start)
claim("el segundo factor empieza en su cara 4", ae[1] == 0 and ae[0] > 0, f"Ae={ae}")

raw_x = fline(e_start, x0, v0) + fline(ae, x0, v0)
raw_v = fline([F(0), F(1)], x0, v0) + fline([A[0][1], A[1][1]], x0, v0)
s = dot(ell_join, raw_x)
X = [z / s for z in raw_x]
U = [(raw_v[i] - X[i] * dot(ell_join, raw_v)) / s for i in range(10)]
X_start, U_start = X[:], U[:]
claim("inicio en la cara 11 (= segunda copia de la cara 4)", dot(facets[11], X) == 0)
claim("inicio interior relativo", all(dot(facets[i], X) < 0 for i in range(14) if i != 11))
claim("tangente afín", dot(ell_join, U) == 0)

print("\nC4. Catorce primeros impactos, todos estrictos")
current = 11
orbit = [X[:]]
for step, wanted in enumerate(itinerary, 1):
    hits = []
    for i, row in enumerate(facets):
        if i == current:
            continue
        den = dot(row, U)
        if den:
            t = -dot(row, X) / den
            if t > 0:
                hits.append((t, i))
    hits.sort()
    claim(f"paso {step}: existe impacto", bool(hits))
    if not hits:
        break
    t, face = hits[0]
    claim(f"paso {step}: cara prevista {wanted}", face == wanted,
          f"cara={face}, t={t}")
    claim(f"paso {step}: perdedores estrictamente posteriores",
          all(other_t > t for other_t, _ in hits[1:]))
    X = [X[i] + t * U[i] for i in range(10)]
    claim(f"paso {step}: impacto interior relativo",
          all(dot(facets[i], X) < 0 for i in range(14) if i != face))
    W = mv(H[face], U)
    U = [W[i] - dot(ell_join, W) * X[i] for i in range(10)]
    current = face
    orbit.append(X[:])

print("\nC5. Cierre, orientación y primitivad")
claim("p14=p0", X == X_start)
nonzero = next(i for i, z in enumerate(U_start) if z)
kappa = U[nonzero] / U_start[nonzero]
claim("u14 proporcional a u0",
      all(U[i] * U_start[j] == U[j] * U_start[i]
          for i in range(10) for j in range(10)))
claim("factor de retorno positivo", kappa > 0, f"kappa={kappa}")
claim("catorce facetas distintas", len(set(itinerary)) == 14)

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodo pasa exactamente sobre Q: billar 14-reflectivo en R^9.")
