# -*- coding: utf-8 -*-
"""Certificado exacto de un join mixto (n,k)=(6,11).

Une el pentágono racional ``thm:second`` con el caleidoscopio D4 de
``reflectivo_rp3_d4.py``. La homografía racional de sincronización es
A=[[-70,280],[-1,2]]. Todas las decisiones se hacen sobre Q.
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
    n = len(normal)
    return [[F(int(i == j)) - F(2 * normal[i] * normal[j], d)
             for j in range(n)] for i in range(n)]


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


# Factor P: pentágono racional en RP2.
NP = [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)]
EP = (1, -1, 1, 1, -1)
RP = [reflection(n) for n in NP]
LP = [[F(EP[i] * NP[i][j]) for j in range(3)] for i in range(5)]
VHP = [cross3(NP[(i - 1) % 5], NP[i]) for i in range(5)]
VP = [(F(v[0], v[2]), F(v[1], v[2])) for v in VHP]


def ppoint(i, t):
    return tuple(VP[i][j] + t * (VP[(i + 1) % 5][j] - VP[i][j]) for j in range(2))


pp0, pp1 = ppoint(0, F(1, 14)), ppoint(1, F(1, 14))
XP = [pp0[0], pp0[1], F(1)]
UP = [pp1[0] - pp0[0], pp1[1] - pp0[1], F(0)]

# Factor D: ejemplo D4 en RP3.
ND = [(1, 0, 1, 0), (1, 0, -1, 0), (1, 0, 0, -1),
      (0, 1, 1, 0), (0, 1, -1, 0), (0, 1, 0, -1)]
ED = (-1, -1, 1, -1, -1, 1)
RD = [reflection(n) for n in ND]
LD = [[F(ED[i] * ND[i][j]) for j in range(4)] for i in range(6)]
XD = [F(4, 7), F(4, 5), F(-4, 7), F(1)]
UD = [F(15), F(5), F(3), F(0)]

print("C1. Monodromías de los dos factores")
prod_p = eye(3)
for i in [1, 2, 3, 4, 0]:
    prod_p = mm(RP[i], prod_p)
prod_d = eye(4)
for i in [2, 1, 3, 5, 4, 0]:
    prod_d = mm(RD[i], prod_d)
claim("pentágono: producto = -I3", prod_p == [[-x for x in row] for row in eye(3)])
claim("D4: producto = -I4", prod_d == [[-x for x in row] for row in eye(4)])

H = [block_diag(r, eye(4)) for r in RP] + [block_diag(eye(3), r) for r in RD]
itinerary = [1, 2, 3, 4, 0, 7, 6, 8, 10, 9, 5]
total = eye(7)
for i in itinerary:
    total = mm(H[i], total)
claim("join: producto = -I7", total == [[-x for x in row] for row in eye(7)])

print("\nC2. Politopo join acotado de once facetas")
facets = [row + [F(0)] * 4 for row in LP] + [[F(0)] * 3 + row for row in LD]
ell_join = [F(0), F(0), F(1), F(0), F(0), F(0), F(1)]
xint_p2 = (sum(v[0] for v in VP) / 5, sum(v[1] for v in VP) / 5)
xint_raw = [xint_p2[0], xint_p2[1], F(1), F(1, 2), F(1, 2), F(-1, 4), F(1)]
xint = [z / dot(ell_join, xint_raw) for z in xint_raw]
claim("punto interior racional", all(dot(row, xint) < 0 for row in facets))
mu_p = [F(47), F(10), F(1), F(2), F(1)]
mu_d = [F(1, 2), F(1, 2), F(1), F(1), F(1), F(2)]
claim("dependencia positiva en el pentágono",
      all(sum(mu_p[i] * LP[i][j] for i in range(5)) == 0 for j in range(2)))
claim("dependencia positiva en D4",
      all(sum(mu_d[i] * LD[i][j] for i in range(6)) == 0 for j in range(3)))

print("\nC3. Sincronización racional")
A = [[F(-70), F(280)], [F(-1), F(2)]]
e_start = [F(1), F(1, 2)]
ae = mv(A, e_start)
det_a = A[0][0] * A[1][1] - A[0][1] * A[1][0]
claim("det(A)>0", det_a > 0, f"det={det_a}")
claim("D4 empieza en su cara 0", ae[1] == 0 and ae[0] > 0, f"Ae={ae}")


def fline(e, x, u):
    return [e[0] * x[i] + e[1] * u[i] for i in range(len(x))]


raw_x = fline(e_start, XP, UP) + fline(ae, XD, UD)
raw_v = fline([F(0), F(1)], XP, UP) + fline([A[0][1], A[1][1]], XD, UD)
s = dot(ell_join, raw_x)
X = [z / s for z in raw_x]
U = [(raw_v[i] - X[i] * dot(ell_join, raw_v)) / s for i in range(7)]
X_start, U_start = X[:], U[:]
claim("inicio en la cara 5 (= cara 0 de D4)", dot(facets[5], X) == 0)
claim("inicio interior relativo", all(dot(facets[i], X) < 0 for i in range(11) if i != 5))
claim("tangente afín", dot(ell_join, U) == 0)

print("\nC4. Once primeros impactos estrictos")
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
    X = [X[i] + t * U[i] for i in range(7)]
    claim(f"paso {step}: impacto interior relativo",
          all(dot(facets[i], X) < 0 for i in range(11) if i != face))
    W = mv(H[face], U)
    U = [W[i] - dot(ell_join, W) * X[i] for i in range(7)]
    current = face

print("\nC5. Cierre, orientación y periodo")
claim("p11=p0", X == X_start)
nonzero = next(i for i, z in enumerate(U_start) if z)
kappa = U[nonzero] / U_start[nonzero]
claim("u11 proporcional a u0",
      all(U[i] * U_start[j] == U[j] * U_start[i]
          for i in range(7) for j in range(7)))
claim("factor de retorno positivo", kappa > 0, f"kappa={kappa}")
claim("once facetas distintas", len(set(itinerary)) == 11)

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodo pasa exactamente sobre Q: billar 11-reflectivo en R^6.")
