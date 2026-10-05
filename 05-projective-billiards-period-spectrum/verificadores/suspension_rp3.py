# -*- coding: utf-8 -*-
"""Certificado exacto de la suspensión de rango uno: (2,5) -> (3,6).

Parte del pentágono racional ``thm:second`` y usa el funcional
    h(a,b) = -9 a + 19 b
sobre su recta desarrollada. Código 0 significa que todos los cálculos y todas
las desigualdades se han comprobado exactamente sobre Q.
"""
import sys
from fractions import Fraction as F


FAIL = []


def claim(label, ok, extra=""):
    print(f"[{'OK  ' if ok else 'FAIL'}] {label}" + (f"  {extra}" if extra else ""))
    if not ok:
        FAIL.append(label)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


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


def cross3(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


ROOTS = [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)]
EPS = (1, -1, 1, 1, -1)
L = [[F(EPS[i] * ROOTS[i][j]) for j in range(3)] + [F(0)] for i in range(5)]
L.append([F(0), F(0), F(0), F(-1)])       # interior: t > 0
ELL = [F(0), F(0), F(1), F(1)]            # z+t=1

R3 = [reflection(n) for n in ROOTS]
R = [[row + [F(0)] for row in a] + [[F(0), F(0), F(0), F(1)]] for a in R3]
R.append([[F(1), F(0), F(0), F(0)],
          [F(0), F(1), F(0), F(0)],
          [F(0), F(0), F(1), F(0)],
          [F(0), F(0), F(0), F(-1)]])

VH = [cross3(ROOTS[(i - 1) % 5], ROOTS[i]) for i in range(5)]
V = [(F(v[0], v[2]), F(v[1], v[2])) for v in VH]


def point(i, t):
    return tuple(V[i][j] + t * (V[(i + 1) % 5][j] - V[i][j]) for j in range(2))


p0 = point(0, F(1, 14))
p1 = point(1, F(1, 14))
x0 = [p0[0], p0[1], F(1)]
u0 = [p1[0] - p0[0], p1[1] - p0[1], F(0)]


def fline(e):
    return [e[0] * x0[i] + e[1] * u0[i] for i in range(3)]


print("C1. Monodromía de la pirámide")
word = [1, 2, 3, 4, 0, 5]
total = eye(4)
for i in word:
    total = mm(R[i], total)
claim("producto de las seis reflexiones = -I4",
      total == [[-x for x in row] for row in eye(4)])
for i, a in enumerate(R):
    claim(f"R_{i}^2=I", mm(a, a) == eye(4))

print("\nC2. Cono piramidal acotado y seis facetas genuinas")
mu = [F(47), F(10), F(1), F(2), F(1)]
old_sum = [sum(mu[i] * L[i][j] for i in range(5)) for j in range(4)]
claim("dependencia positiva antigua: sum mu_i L_i = -152 z",
      old_sum == [F(0), F(0), F(-152), F(0)])
claim("la nueva desigualdad fuerza t>0", L[5] == [F(0), F(0), F(0), F(-1)])
claim("por tanto ell=z+t es positiva en el cierre del cono salvo en 0",
      old_sum[2] < 0 and L[5][3] < 0)

e_start = [F(1), F(1, 2)]
h = [F(-9), F(19)]
raw_x = fline(e_start) + [dot(h, e_start)]
raw_v = fline([F(0), F(1)]) + [h[1]]
s = dot(ELL, raw_x)
X = [z / s for z in raw_x]
U = [(raw_v[i] - X[i] * dot(ELL, raw_v)) / s for i in range(4)]
claim("punto de lanzamiento estrictamente interior", all(dot(row, X) < 0 for row in L))
claim("tangente en la carta z+t=1", dot(ELL, U) == 0)

print("\nC3. Seis primeros impactos estrictos")
hits, points, directions = [], [], []
current = -1
for step in range(7):
    candidates = []
    for i, row in enumerate(L):
        if i == current:
            continue
        den = dot(row, U)
        if den:
            t = -dot(row, X) / den
            if t > 0:
                candidates.append((t, i))
    candidates.sort()
    claim(f"paso {step + 1}: existe primer impacto", bool(candidates))
    if not candidates:
        break
    t, face = candidates[0]
    claim(f"paso {step + 1}: ganador único", all(tt > t for tt, _ in candidates[1:]))
    X = [X[i] + t * U[i] for i in range(4)]
    claim(f"paso {step + 1}: interior relativo",
          all(dot(L[i], X) < 0 for i in range(6) if i != face))
    W = mv(R[face], U)
    U = [W[i] - dot(ELL, W) * X[i] for i in range(4)]
    hits.append(face)
    points.append(X[:])
    directions.append(U[:])
    current = face

claim("itinerario esperado", hits == [1, 2, 3, 4, 0, 5, 1], f"hits={hits}")
claim("seis facetas distintas en el periodo", len(set(hits[:6])) == 6)

print("\nC4. Cierre exacto y grupo infinito")
if len(points) == 7:
    claim("p_6=p_0 exactamente", points[6] == points[0])
    nz = next(i for i, z in enumerate(directions[0]) if z)
    kappa = directions[6][nz] / directions[0][nz]
    claim("u_6=kappa u_0", all(directions[6][i] == kappa * directions[0][i]
                                  for i in range(4)))
    claim("kappa=1>0", kappa == 1, f"kappa={kappa}")

pair = mm(R3[0], R3[1])
pair_trace = sum(pair[i][i] for i in range(3))
claim("tr(R0 R1)=-67/85, fuera de los valores racionales de orden finito",
      pair_trace == F(-67, 85), f"tr={pair_trace}")
claim("el grupo viejo se inyecta por g -> diag(g,1)",
      all(R[i][3] == [F(0), F(0), F(0), F(1)] for i in range(5)))

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodo pasa exactamente sobre Q: suspensión 6-reflectiva en R^3.")
