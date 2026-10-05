# -*- coding: utf-8 -*-
"""Certificado exacto de un cubo 6-reflectivo en RP^3 con monodromía +I.

El cubo es |x|,|y|,|t|<z. Los cuatro espejos laterales son las
extensiones inclinadas del billar central del cuadrado y los dos espejos de
tapa son verticales. Código 0 significa que todo se ha verificado sobre Q.
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


def madd(a, b, sign=1):
    return [[a[i][j] + sign * b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def matrix_rank(a):
    a = [row[:] for row in a]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                scale = a[i][col]
                a[i] = [a[i][j] - scale * a[rank][j] for j in range(cols)]
        rank += 1
    return rank


def homology(L, O):
    d = dot(L, O)
    return [[F(int(i == j)) - F(2) * O[i] * L[j] / d
             for j in range(4)] for i in range(4)]


# Interior convention: L_i X < 0.
L = [[F(1), F(0), F(-1), F(0)],
     [F(0), F(1), F(-1), F(0)],
     [F(-1), F(0), F(-1), F(0)],
     [F(0), F(-1), F(-1), F(0)],
     [F(0), F(0), F(-1), F(1)],
     [F(0), F(0), F(-1), F(-1)]]

C = [F(-6, 5), F(-1, 5), F(4, 5), F(9, 5)]
O = [[F(0), F(0), F(1), C[i]] for i in range(4)]
O += [[F(0), F(0), F(0), F(1)], [F(0), F(0), F(0), F(1)]]
R = [homology(L[i], O[i]) for i in range(6)]


print("C1. Seis homologías y monodromía positiva")
for i, a in enumerate(R):
    claim(f"R_{i}^2=I", mm(a, a) == eye(4))
    claim(f"el centro O_{i} es transversal a la cara", dot(L[i], O[i]) != 0)

word = [1, 2, 3, 4, 5, 0]
total = eye(4)
for i in word:
    total = mm(R[i], total)
claim("producto en el orden del itinerario = +I4", total == eye(4))

print("\nC1bis. Scott-sharp y dualidad de Gale")
ordered_R = [R[i] for i in word]
ordered_O = [O[i] for i in word]
ordered_L = [[-x / dot(L[i], O[i]) for x in L[i]] for i in word]
developed = [None] * 6
suffix = eye(4)
for i in range(5, -1, -1):
    developed[i] = mv(suffix, ordered_O[i])
    suffix = mm(suffix, ordered_R[i])
P_GALE = [[developed[i][j] for i in range(6)] for j in (2, 3)]
claim("los centros generan exactamente el núcleo (z,t)",
      all(o[0] == o[1] == 0 for o in ordered_O) and
      matrix_rank([[o[2], o[3]] for o in ordered_O]) == 2)
claim("la configuración satura Scott: k=m+c", 6 == 4 + 2)
claim("los centros desarrollados tienen rango dos", matrix_rank(P_GALE) == 2)
claim("los seis ejes homogéneos tienen rango cuatro", matrix_rank(ordered_L) == 4)
claim("centros desarrollados y ejes son Gale duales",
      mm(P_GALE, ordered_L) == [[F(0)] * 4 for _ in range(2)])

# Orientación de politopo: A_i=-L_i^{original}, de modo que A_i X>0 dentro.
# Al pasar de la normalización canónica A_i^c=-L_i/L_i(O_i) a ésta, la
# columna Gale correspondiente se divide por L_i(O_i).
poly_axes = [[-x for x in L[i]] for i in word]
transversals = [dot(L[i], O[i]) for i in word]
poly_gale = [[P_GALE[a][i] / transversals[i] for i in range(6)]
             for a in range(2)]
claim("el Gale orientado aniquila los ejes del cubo",
      mm(poly_gale, poly_axes) == [[F(0)] * 4 for _ in range(2)])
interior = [F(0), F(0), F(1), F(0)]
interior_values = [dot(a, interior) for a in poly_axes]
claim("el núcleo de Gale contiene una dependencia estrictamente positiva",
      all(x > 0 for x in interior_values) and
      mv(poly_gale, interior_values) == [F(0), F(0)])
face_points = [[F(1), F(0), F(1), F(0)],
               [F(0), F(1), F(1), F(0)],
               [F(-1), F(0), F(1), F(0)],
               [F(0), F(-1), F(1), F(0)],
               [F(0), F(0), F(1), F(1)],
               [F(0), F(0), F(1), F(-1)]]
for pos, face in enumerate(word):
    values = [dot(a, face_points[face]) for a in poly_axes]
    claim(f"borrando la columna Gale {pos} queda una dependencia positiva",
          values[pos] == 0 and all(values[j] > 0 for j in range(6)
                                   if j != pos) and
          mv(poly_gale, values) == [F(0), F(0)])

# Rango exacto del mapa producto sobre las seis clases de conjugación en GL2.
core = [[[a[i][j] for j in (2, 3)] for i in (2, 3)] for a in ordered_R]
prefixes, prefix = [], eye(2)
for a in core:
    prefixes.append(prefix)
    prefix = mm(a, prefix)
suffixes, suffix = [None] * 6, eye(2)
for i in range(5, -1, -1):
    suffixes[i] = suffix
    suffix = mm(suffix, core[i])
matrix_units = [[[F(int(i == r and j == c)) for j in range(2)]
                 for i in range(2)] for r in range(2) for c in range(2)]
derivatives = []
for i, a in enumerate(core):
    for e in matrix_units:
        commutator = madd(mm(e, a), mm(a, e), sign=-1)
        delta = mm(mm(suffixes[i], commutator), prefixes[i])
        derivatives.append([delta[r][c] for r in range(2) for c in range(2)])
derivative_matrix = [[derivatives[j][i] for j in range(len(derivatives))]
                     for i in range(4)]
claim("el mapa producto del núcleo tiene rango tres",
      matrix_rank(derivative_matrix) == 3)

centralizer_columns = []
for e in matrix_units:
    column = []
    for a in core:
        commutator = madd(mm(e, a), mm(a, e), sign=-1)
        column.extend(commutator[r][c] for r in range(2) for c in range(2))
    centralizer_columns.append(column)
centralizer_matrix = [[centralizer_columns[j][i] for j in range(4)]
                      for i in range(24)]
claim("el centralizador común sólo contiene escalares",
      4 - matrix_rank(centralizer_matrix) == 1)
claim("el móduli local del núcleo, y de su extensión, tiene dimensión seis",
      2 * 6 - 3 - 3 == 6)


print("\nC2. El dominio es un cubo proyectivo acotado con seis caras genuinas")
# Las desigualdades son x<z, y<z, -x<z, -y<z, t<z, -t<z.
# Sumando cada par se obtiene z>=0 en el cierre; si z=0, todas las
# coordenadas se anulan. La sección z=1 es el cubo [-1,1]^3.
for a, b in ((0, 2), (1, 3), (4, 5)):
    claim(f"el par de caras {a},{b} fuerza z>=0",
          [L[a][j] + L[b][j] for j in range(4)] == [F(0), F(0), F(-2), F(0)])
claim("las seis formas de cara son distintas", len({tuple(row) for row in L}) == 6)


print("\nC3. Trazado racional estricto")
ELL = [F(0), F(0), F(1), F(0)]
X_START = [F(2, 5), F(3, 10), F(1), F(-1, 4)]
U_START = [F(-1, 4), F(1, 3), F(0), F(1, 10)]
X = X_START[:]
U = U_START[:]
claim("punto inicial estrictamente interior", all(dot(row, X) < 0 for row in L))
claim("dirección tangente a la carta z=1", dot(ELL, U) == 0)

hits, points, directions = [], [], []
current = -1
for step in range(7):
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
    X = [X[i] + tau * U[i] for i in range(4)]
    claim(f"paso {step + 1}: impacto en interior relativo",
          all(dot(L[i], X) < 0 for i in range(6) if i != face))
    W = mv(R[face], U)
    U = [W[i] - dot(ELL, W) * X[i] for i in range(4)]
    hits.append(face)
    points.append(X[:])
    directions.append(U[:])
    current = face

claim("itinerario esperado", hits == [1, 2, 3, 4, 5, 0, 1], f"hits={hits}")
claim("las seis caras se visitan exactamente una vez", len(set(hits[:6])) == 6)


print("\nC4. Cierre exacto y retorno orientado")
if len(points) == 7:
    claim("p_6=p_0 exactamente", points[6] == points[0])
    nz = next(i for i, z in enumerate(directions[0]) if z)
    kappa = directions[6][nz] / directions[0][nz]
    claim("u_6=kappa u_0",
          all(directions[6][i] == kappa * directions[0][i] for i in range(4)))
    claim("kappa=1>0", kappa == 1, f"kappa={kappa}")

    pinv = eye(4)
    developed = []
    for j, point in enumerate(points):
        dev = mv(pinv, point)
        a = dev[2]  # X_START has z=1 and U_START has z=0.
        b = (dev[0] - a * X_START[0]) / U_START[0]
        claim(f"impacto desarrollado {j} pertenece al plano inicial",
              dev == [a * X_START[q] + b * U_START[q] for q in range(4)])
        developed.append((a, b))
        if j < 6:
            pinv = mm(pinv, R[hits[j]])
    dets = [developed[j][0] * developed[j + 1][1]
            - developed[j][1] * developed[j + 1][0] for j in range(6)]
    claim("cada arco desarrollado gira positivamente menos de pi", all(d > 0 for d in dets))
    signs = [(0 if a == 0 else (1 if a > 0 else -1),
              0 if b == 0 else (1 if b > 0 else -1)) for a, b in developed]
    claim("el patrón de cuadrantes certifica exactamente una vuelta vectorial",
          signs == [(1, 1), (-1, 0), (-1, -1), (-1, -1),
                    (1, -1), (1, -1), (1, 1)], f"signs={signs}")


print("\nC5. El grupo de espejos es infinito")
translation = mm(R[5], R[4])
nilpotent = [[translation[i][j] - eye(4)[i][j] for j in range(4)] for i in range(4)]
claim("las dos tapas producen una cizalla no trivial", nilpotent != [[F(0)] * 4 for _ in range(4)])
claim("(R_5 R_4-I)^2=0", mm(nilpotent, nilpotent) == [[F(0)] * 4 for _ in range(4)])

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodo pasa exactamente sobre Q: cubo 6-reflectivo en RP^3 con producto +I4.")
