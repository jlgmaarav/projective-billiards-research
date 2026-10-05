# -*- coding: utf-8 -*-
"""Certificado exacto de un hipercubo 8-reflectivo en RP^4 con monodromía +I.

El dominio es |x|,|y|,|t|,|s|<z. Código 0 significa que involuciones,
monodromía, acotación, trazado, cierre, una vuelta e infinitud pasan sobre Q.
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


def homology(L, O):
    q = dot(L, O)
    return [[F(int(i == j)) - F(2) * O[i] * L[j] / q
             for j in range(5)] for i in range(5)]


# Coordinates are (x,y,z,t,s); the interior convention is L_i X < 0.
L = [[F(x) for x in row] for row in
     [(1,0,-1,0,0), (0,1,-1,0,0), (-1,0,-1,0,0),
      (0,-1,-1,0,0), (0,0,-1,1,0), (0,0,-1,-1,0),
      (0,0,-1,0,1), (0,0,-1,0,-1)]]

c = [F(-6, 5), F(-1, 5), F(4, 5), F(9, 5)]
d = [F(-3), F(-17, 10), F(-2, 5), F(9, 10), F(3, 10), F(3, 10)]
O = [[F(0), F(0), F(1), c[i], d[i]] for i in range(4)]
O += [[F(0), F(0), F(0), F(1), d[i]] for i in range(4, 6)]
O += [[F(0), F(0), F(0), F(0), F(1)],
      [F(0), F(0), F(0), F(0), F(1)]]
R = [homology(L[i], O[i]) for i in range(8)]


print("C1. Ocho homologías y monodromía positiva")
for i, a in enumerate(R):
    claim(f"R_{i}^2=I", mm(a, a) == eye(5))
    claim(f"el centro O_{i} es transversal a la cara", dot(L[i], O[i]) != 0)

word = [1, 2, 3, 4, 5, 6, 7, 0]
total = eye(5)
for i in word:
    total = mm(R[i], total)
claim("producto en el orden del itinerario = +I5", total == eye(5))


print("\nC2. El dominio es un hipercubo proyectivo acotado")
for a, b in ((0, 2), (1, 3), (4, 5), (6, 7)):
    claim(f"el par de caras {a},{b} fuerza z>=0",
          [L[a][j] + L[b][j] for j in range(5)] ==
          [F(0), F(0), F(-2), F(0), F(0)])
claim("las ocho formas de cara son distintas", len({tuple(row) for row in L}) == 8)


print("\nC3. Nueve primeros impactos exactos")
ELL = [F(0), F(0), F(1), F(0), F(0)]
X_START = [F(4, 11), F(8, 11), F(1), F(2, 3), F(-8, 15)]
U_START = [F(-4, 5), F(1, 3), F(0), F(-2, 15), F(1, 2)]
X, U = X_START[:], U_START[:]
claim("punto inicial estrictamente interior", all(dot(row, X) < 0 for row in L))
claim("dirección tangente a z=1", dot(ELL, U) == 0)

hits, points, directions = [], [], []
current = -1
for step in range(9):
    candidates = []
    for i, row in enumerate(L):
        if i == current:
            continue
        den = dot(row, U)
        if den:
            tau = -dot(row, X) / den
            if tau > 0:
                candidates.append((tau, i))
    candidates.sort()
    claim(f"paso {step + 1}: existe primer impacto", bool(candidates))
    if not candidates:
        break
    tau, face = candidates[0]
    claim(f"paso {step + 1}: ganador único",
          all(other_tau > tau for other_tau, _ in candidates[1:]))
    X = [X[i] + tau * U[i] for i in range(5)]
    claim(f"paso {step + 1}: interior relativo",
          all(dot(L[i], X) < 0 for i in range(8) if i != face))
    W = mv(R[face], U)
    U = [W[i] - dot(ELL, W) * X[i] for i in range(5)]
    hits.append(face)
    points.append(X[:])
    directions.append(U[:])
    current = face

claim("itinerario esperado", hits == [1,2,3,4,5,6,7,0,1], f"hits={hits}")
claim("las ocho caras se visitan una vez", len(set(hits[:8])) == 8)


print("\nC4. Cierre y una sola vuelta vectorial")
if len(points) == 9:
    claim("p_8=p_0 exactamente", points[8] == points[0])
    nz = next(i for i, z in enumerate(directions[0]) if z)
    kappa = directions[8][nz] / directions[0][nz]
    claim("u_8=kappa u_0",
          all(directions[8][i] == kappa * directions[0][i] for i in range(5)))
    claim("kappa=1>0", kappa == 1, f"kappa={kappa}")

    pinv = eye(5)
    developed = []
    for j, point in enumerate(points):
        dev = mv(pinv, point)
        a = dev[2]
        b = (dev[0] - a * X_START[0]) / U_START[0]
        claim(f"impacto desarrollado {j} está en el plano inicial",
              dev == [a * X_START[q] + b * U_START[q] for q in range(5)])
        developed.append((a, b))
        if j < 8:
            pinv = mm(pinv, R[hits[j]])
    dets = [developed[j][0] * developed[j + 1][1]
            - developed[j][1] * developed[j + 1][0] for j in range(8)]
    claim("todos los arcos giran positivamente menos de pi", all(q > 0 for q in dets))
    signs = [(1 if a > 0 else -1, 1 if b > 0 else -1) for a, b in developed]
    claim("el patrón de cuadrantes da exactamente una vuelta",
          signs == [(1,1),(1,1),(-1,-1),(-1,-1),(1,-1),
                    (1,-1),(1,-1),(1,-1),(1,1)], f"signs={signs}")


print("\nC5. El grupo de espejos es infinito")
translation = mm(R[7], R[6])
N = [[translation[i][j] - eye(5)[i][j] for j in range(5)] for i in range(5)]
zero = [[F(0)] * 5 for _ in range(5)]
claim("las tapas producen un unipotente no trivial", N != zero)
claim("(R_7 R_6-I)^2=0", mm(N, N) == zero)

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodo pasa sobre Q: hipercubo 8-reflectivo en RP^4 con producto +I5.")
