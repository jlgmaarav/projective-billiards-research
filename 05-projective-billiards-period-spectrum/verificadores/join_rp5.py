# -*- coding: utf-8 -*-
"""Certificado exacto del join de dos pentágonos: (n,k)=(5,10).

El factor es el ejemplo racional ``thm:second`` de ``verifier_rational.py``.
La sincronización usa A=[[-4,16],[-1,2]] y e_*=(1,1/2). Todas las
decisiones se hacen en Q; código 0 significa que todas las comprobaciones pasan.
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


def reflection(normal):
    d = sum(F(x * x) for x in normal)
    return [[F(int(i == j)) - F(2 * normal[i] * normal[j], d)
             for j in range(3)] for i in range(3)]


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


def cross3(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


ROOTS = [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)]
R = [reflection(n) for n in ROOTS]
VH = [cross3(ROOTS[(i - 1) % 5], ROOTS[i]) for i in range(5)]
V = [(F(v[0], v[2]), F(v[1], v[2])) for v in VH]


def point(i, t):
    return tuple(V[i][j] + t * (V[(i + 1) % 5][j] - V[i][j]) for j in range(2))


# Orientaciones del pentágono en la carta z=1.
EPS = (1, -1, 1, 1, -1)
L = [[F(EPS[i] * ROOTS[i][j]) for j in range(3)] for i in range(5)]
ELL = [F(0), F(0), F(1)]


print("C1. Factor racional y monodromía")
factor_product = eye(3)
for i in [1, 2, 3, 4, 0]:
    factor_product = mm(R[i], factor_product)
claim("producto de cada pentágono = -I3",
      factor_product == [[-x for x in row] for row in eye(3)])

H = [block_diag(r, eye(3)) for r in R] + [block_diag(eye(3), r) for r in R]
itinerary = [1, 2, 3, 4, 0, 6, 7, 8, 9, 5]
total = eye(6)
for i in itinerary:
    total = mm(H[i], total)
claim("producto de las diez reflexiones = -I6",
      total == [[-x for x in row] for row in eye(6)])

print("\nC2. El join es un politopo acotado de diez facetas")
facets = [row + [F(0)] * 3 for row in L] + [[F(0)] * 3 + row for row in L]
ell_join = ELL + ELL
xint2 = (sum(p[0] for p in V) / 5, sum(p[1] for p in V) / 5)
xint_raw = [xint2[0], xint2[1], F(1)] * 2
xint = [z / dot(ell_join, xint_raw) for z in xint_raw]
claim("punto interior racional", all(dot(row, xint) < 0 for row in facets))
mu = [F(47), F(10), F(1), F(2), F(1)]
claim("dependencia positiva de normales espaciales",
      all(sum(mu[i] * L[i][j] for i in range(5)) == 0 for j in range(2)))

print("\nC3. Sincronización racional")
p0 = point(0, F(1, 14))
p1 = point(1, F(1, 14))
x0 = [p0[0], p0[1], F(1)]
u0 = [p1[0] - p0[0], p1[1] - p0[1], F(0)]
A = [[F(-4), F(16)], [F(-1), F(2)]]
det_a = A[0][0] * A[1][1] - A[0][1] * A[1][0]
e_start = [F(1), F(1, 2)]
ae = mv(A, e_start)
claim("det(A)>0", det_a > 0, f"det={det_a}")
claim("Ae_* es un múltiplo positivo del inicio del segundo factor",
      ae[1] == 0 and ae[0] > 0, f"Ae={ae}")


def fline(e):
    return [e[0] * x0[i] + e[1] * u0[i] for i in range(3)]


raw_x = fline(e_start) + fline(ae)
raw_v = fline([F(0), F(1)]) + fline([A[0][1], A[1][1]])
s = dot(ell_join, raw_x)
X = [z / s for z in raw_x]
U = [(raw_v[i] - X[i] * dot(ell_join, raw_v)) / s for i in range(6)]
X_start, U_start = X[:], U[:]
claim("inicio en la cara 5 (= segunda copia de la cara 0)", dot(facets[5], X) == 0)
claim("inicio interior relativo", all(dot(facets[i], X) < 0 for i in range(10) if i != 5))
claim("tangente afín", dot(ell_join, U) == 0)

print("\nC4. Diez primeros impactos estrictos")
current = 5
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
    X = [X[i] + t * U[i] for i in range(6)]
    claim(f"paso {step}: impacto interior relativo",
          all(dot(facets[i], X) < 0 for i in range(10) if i != face))
    W = mv(H[face], U)
    U = [W[i] - dot(ell_join, W) * X[i] for i in range(6)]
    current = face

print("\nC5. Cierre, orientación y periodo")
claim("p10=p0", X == X_start)
nonzero = next(i for i, z in enumerate(U_start) if z)
kappa = U[nonzero] / U_start[nonzero]
claim("u10 proporcional a u0",
      all(U[i] * U_start[j] == U[j] * U_start[i]
          for i in range(6) for j in range(6)))
claim("factor de retorno positivo", kappa > 0, f"kappa={kappa}")
claim("diez facetas distintas", len(set(itinerary)) == 10)

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodo pasa exactamente sobre Q: billar 10-reflectivo en R^5.")
