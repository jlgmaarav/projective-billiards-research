# -*- coding: utf-8 -*-
"""UN BILLAR 9-REFLECTIVO EN R^6.  TODO EXACTO SOBRE Q.

    python reflectivo_rp6_b7.py     ->  codigo de salida 0 solo si TODO pasa

QUE SE AFIRMA

El politopo convexo ACOTADO de R^6 con las nueve caras de normales

    n_0 = (0,0,0,0,0,1,-1)   n_1 = (0,1,0,1,0,0,0)   n_2 = (0,0,1,0,0,0,0)
    n_3 = (0,0,0,1,0,1, 0)   n_4 = (1,0,0,0,0,0,0)   n_5 = (0,1,0,0,0,0,-1)
    n_6 = (0,0,0,0,1,0, 0)   n_7 = (0,0,0,1,0,0,-1)  n_8 = (0,1,0,0,0,1, 0)

y signos  eps = (-,+,+,+,+,-,-,+,-),  en la carta afin  ell.X = 1  con
ell = (-1,1,-1,-3,1,1,-1)  despejando la CUARTA coordenada homogenea, y con el
campo transversal en la cara i dado por el haz que pasa por el polo O_i = n_i
respecto de Q = x_1^2 + ... + x_7^2, es un billar proyectivo 9-REFLECTIVO: hay un
abierto de orbitas 9-periodicas.

POR QUE IMPORTA

  - Es el primer ejemplo en R^6, y sube en dos la dimension maxima alcanzada.
  - k = 9 es IMPAR, con n = 6 par: exactamente lo que permite el teorema de
    paridad `thm:parity`, que para n impar obligaria a k par.
  - Cae en el umbral k = n+3 = 9, como los cuatro ejemplos anteriores.
  - Y confirma el Teorema B de `DECIMA-RONDA-RIGIDEZ.md` sobre datos que se
    obtuvieron ANTES de demostrarlo: la palabra tiene que ser un bloque
    irreducible de rango <= 4 mas cambios de signo, y lo es -- seis raices de tipo
    e_i +- e_j sobre las coordenadas {1,3,5,6}, que es un bloque D_4, mas los tres
    cambios de signo e_0, e_2, e_4.

PROCEDENCIA DE LOS DATOS

Las nueve raices se perdieron el 8 de agosto de 2026 -- se generaron dentro de la
salida de un script y no se escribieron -- y bloquearon el proyecto durante dos
dias; ver `OCTAVA-RONDA-CANDIDATO-RP6.md`.  Se recuperaron el 8 de agosto por la
tarde.  El corte afin, el punto y la direccion si se habian conservado, y encajan.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  las nueve son reflexiones de O(Q) y el producto en el orden del itinerario
      es -I_7  (hipotesis H1 del criterio de apertura)
  C2  el politopo es ACOTADO: las nueve normales exteriores SUMAN CERO, o sea
      mu = (1,...,1), y tienen rango 6
  C3  el punto de partida esta sobre la cara 6 y es estrictamente interior a las
      otras ocho
  C4  la orbita testigo racional: nueve impactos con itinerario 7,8,0,1,...,6,
      cada uno estrictamente dentro de su cara, cada cara perdedora excluida por
      desigualdad ESTRICTA, y todos los parametros de paso positivos  (H2)
  C5  el cierre: p_9 = p_0 EXACTAMENTE y u_9 = kappa u_0 con kappa > 0  (H3)
  C6  las nueve caras son distintas y las normales tienen rango 7
  C7  la estructura que predice el Teorema B, y el corango 3 que le corresponde
  C8  el bloque de rango 4 es el MISMO nucleo D_4 de `thm:rp4b5`, `thm:rp5b6` y
      `reflectivo_rp3_d4b.py`
"""
import sys
from fractions import Fraction as F

M, K = 7, 9
FAILURES = []


def claim(name, cond, detail=""):
    ok = bool(cond)
    print(f"   [{'OK  ' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)
    return ok


def header(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# ------------------------------------------------------------------- los datos
N = [(0, 0, 0, 0, 0, 1, -1),
     (0, 1, 0, 1, 0, 0, 0),
     (0, 0, 1, 0, 0, 0, 0),
     (0, 0, 0, 1, 0, 1, 0),
     (1, 0, 0, 0, 0, 0, 0),
     (0, 1, 0, 0, 0, 0, -1),
     (0, 0, 0, 0, 1, 0, 0),
     (0, 0, 0, 1, 0, 0, -1),
     (0, 1, 0, 0, 0, 1, 0)]
EPS = [-1, 1, 1, 1, 1, -1, -1, 1, -1]
ELL = [F(-1), F(1), F(-1), F(-3), F(1), F(1), F(-1)]
LIBRE = 3                                   # se despeja X_4 (indice 3)
POS = [c for c in range(M) if c != LIBRE]

P0A = [F(-6, 19), F(7, 228), F(-6, 19), F(0), F(7, 228), F(-7, 228)]
U0A = [F(x, 10**6) for x in (-1003425, -270995, 31421, 1824019, -819948, -545098)]
ITIN = [7, 8, 0, 1, 2, 3, 4, 5, 6]
CARA0 = 6                                   # p_0 esta sobre esta cara

NS = [[F(x) for x in n] for n in N]


def refl(nv):
    d = sum((nv[a] * nv[a] for a in range(M)), F(0))
    assert d != 0
    return [[F(int(a == b)) - F(2) * nv[a] * nv[b] / d for b in range(M)]
            for a in range(M)]


R = [refl(n) for n in NS]


def mul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(M)), F(0)) for j in range(M)]
            for i in range(M)]


def eye():
    return [[F(int(i == j)) for j in range(M)] for i in range(M)]


def rank(rows):
    if not rows:
        return 0
    A = [r[:] for r in rows]
    nr, nc = len(A), len(A[0])
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = F(1) / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == nr:
            break
    return r


def carta(n):
    n0, l0 = n[LIBRE], ELL[LIBRE]
    return [n[c] - n0 * ELL[c] / l0 for c in POS], n0 / l0


CH = [carta(n) for n in NS]
A6 = [a for a, d in CH]
D6 = [d for a, d in CH]


def lado(i, x):
    return sum((A6[i][j] * x[j] for j in range(M - 1)), F(0)) + D6[i]


def lift(v, punto):
    X = [F(0)] * M
    for a, c in enumerate(POS):
        X[c] = v[a]
    X[LIBRE] = ((F(1) if punto else F(0))
                - sum((ELL[c] * X[c] for c in POS), F(0))) / ELL[LIBRE]
    return X


# ---------------------------------------------------------------------- C1
def C1():
    header("C1.  NUEVE REFLEXIONES DE O(Q), Y EL PRODUCTO ES -I_7")
    for i in range(K):
        claim(f"R_{i}^2 = I", mul(R[i], R[i]) == eye())
        claim(f"tr R_{i} = {M-2}  (reflexion)",
              sum((R[i][a][a] for a in range(M)), F(0)) == F(M - 2))
    P = eye()
    for i in range(K):
        P = mul(R[i], P)
    menos = [[F(-1) if a == b else F(0) for b in range(M)] for a in range(M)]
    claim("R_8 R_7 ... R_0 = -I_7", P == menos)


# ---------------------------------------------------------------------- C2
def C2():
    header("C2.  EL POLITOPO ES ACOTADO: LAS NORMALES EXTERIORES SUMAN CERO")
    out = [[F(EPS[i]) * A6[i][j] for j in range(M - 1)] for i in range(K)]
    claim("las normales exteriores tienen rango 6", rank(out) == M - 1)
    for j in range(M - 1):
        s = sum((out[i][j] for i in range(K)), F(0))
        claim(f"suma de las nueve normales exteriores, coordenada {j}: = 0", s == 0,
              f"= {s}")
    print("   Es decir mu = (1,1,1,1,1,1,1,1,1), estrictamente positivo.")


# ---------------------------------------------------------------------- C3
def C3():
    header("C3.  EL PUNTO DE PARTIDA")
    print(f"   p_0 = {[str(x) for x in P0A]}")
    for i in range(K):
        v = F(EPS[i]) * lado(i, P0A)
        if i == CARA0:
            claim(f"cara {i}: p_0 esta SOBRE ella", v == 0, f"= {v}")
        else:
            claim(f"cara {i}: p_0 estrictamente dentro", v < 0, f"= {v}")


# ------------------------------------------------------------------ C4 y C5
def paso(P, D, ult):
    cand = []
    for i in range(K):
        if i == ult:
            continue
        den = sum((NS[i][a] * D[a] for a in range(M)), F(0))
        if den == 0:
            continue
        cand.append((i, -sum((NS[i][a] * P[a] for a in range(M)), F(0)) / den))
    pos = [(i, t) for i, t in cand if t > 0]
    if not pos:
        return None
    best = min(pos, key=lambda c: c[1])
    Q = [P[a] + best[1] * D[a] for a in range(M)]
    Dn = [sum((R[best[0]][a][b] * D[b] for b in range(M)), F(0)) for a in range(M)]
    s = sum((ELL[a] * Dn[a] for a in range(M)), F(0))
    return best[0], best[1], Q, [Dn[a] - s * Q[a] for a in range(M)], cand


def C4C5():
    header("C4 y C5.  LA ORBITA TESTIGO, NUEVE IMPACTOS, Y EL CIERRE")
    print(f"   u_0 = {[str(x) for x in U0A]}")
    P0, D0 = lift(P0A, True), lift(U0A, False)
    P, D, ult = P0, D0, CARA0
    itin, ts = [], []
    for j in range(K):
        r = paso(P, D, ult)
        claim(f"paso {j+1}: hay impacto", r is not None)
        if r is None:
            return
        c, t, P, D, cand = r
        itin.append(c)
        ts.append(t)
        claim(f"paso {j+1}: la cara ganadora es {ITIN[j]}", c == ITIN[j], f"t = {t}")
        claim(f"paso {j+1}: parametro de paso ESTRICTAMENTE positivo", t > 0)
        claim(f"paso {j+1}: las demas caras, excluidas por desigualdad ESTRICTA",
              all(tt <= 0 or tt > t for i, tt in cand if i != c))
        h = sum((ELL[a] * P[a] for a in range(M)), F(0))
        claim(f"paso {j+1}: el impacto no esta en el infinito", h != 0)
        if h != 0:
            xaf = [P[c2] / h for c2 in POS]
            claim(f"paso {j+1}: el impacto es interior a su cara",
                  all(F(EPS[i]) * lado(i, xaf) < 0 for i in range(K) if i != c))
        ult = c
    claim("el itinerario recorre las nueve caras sin repetir",
          len(set(itin)) == K, f"itinerario = {itin}")
    claim("CIERRE: p_9 = p_0 EXACTAMENTE", P == P0)
    kap = next((D[a] / D0[a] for a in range(M) if D0[a] != 0), None)
    claim("u_9 = kappa u_0", kap is not None
          and all(D[a] == kap * D0[a] for a in range(M)))
    claim("y kappa > 0: misma orientacion", kap is not None and kap > 0,
          f"kappa = {kap}")


# ---------------------------------------------------------------------- C6
def C6():
    header("C6.  NO DEGENERADA")
    claim("las nueve caras son distintas", len({tuple(n) for n in N}) == K)
    claim("las normales tienen rango 7", rank([list(n) for n in NS]) == M)


# ---------------------------------------------------------------------- C7
def C7():
    header("C7.  LA ESTRUCTURA QUE PREDICE EL TEOREMA B")
    signos, transp = [], []
    for i, n in enumerate(N):
        nz = [c for c in range(M) if n[c] != 0]
        (signos if len(nz) == 1 else transp).append(i)
    coords = sorted({c for i in transp for c in range(M) if N[i][c] != 0})
    r = len(coords)
    claim("el bloque no trivial tiene rango <= 4", r <= 4, f"rango = {r}")
    claim(f"y usa exactamente r+2 = {r+2} raices", len(transp) == r + 2,
          f"tiene {len(transp)}")
    claim(f"mas 7-r = {M-r} cambios de signo", len(signos) == M - r,
          f"tiene {len(signos)} en coordenadas "
          f"{[c for i in signos for c in range(M) if N[i][c] != 0]}")
    claim("todas las raices del bloque son de tipo e_i +- e_j (o sea, de D_4)",
          all(sorted(abs(N[i][c]) for c in range(M) if N[i][c] != 0) == [1, 1]
              for i in transp))
    # corango: cuadricas invariantes de mas
    basis = [(i, j) for i in range(M) for j in range(i, M)]

    def comp(i, j, a, rr):
        v = F(0)
        if i == rr:
            v += a[j]
        if j == rr and i != j:
            v += a[i]
        return v

    rows = []
    for a in NS:
        for p in range(M):
            for q in range(p + 1, M):
                rows.append([a[p] * comp(i, j, a, q) - a[q] * comp(i, j, a, p)
                             for i, j in basis])
    rows.append([F(1) if i == j else F(0) for i, j in basis])
    co = len(basis) - rank(rows)
    claim("el corango es 3, uno por cada coordenada aislada", co == 3,
          f"corango = {co}")
    print("   Luego la configuracion es REDUCIBLE, como obliga el Teorema B, y esta")
    print("   en el estrato degenerado descrito en NOVENA-RONDA-TOLERANCIAS.md.")


# ---------------------------------------------------------------------- C8
def C8():
    header("C8.  EL BLOQUE ES EL MISMO NUCLEO D_4 DE B_5 Y B_6")
    coords = [1, 3, 5, 6]
    bloque = [tuple(N[i][c] for c in coords) for i in range(K)
              if sum(1 for c in range(M) if N[i][c] != 0) == 2]
    NUCLEO = [(0, 1, 0, -1), (0, 0, 1, -1), (1, -1, 0, 0),
              (0, 1, 1, 0), (1, 0, 0, 1), (1, 0, -1, 0)]

    def gram(S):
        out = []
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                a, b = S[i], S[j]
                num = sum(F(a[t] * b[t]) for t in range(4)) ** 2
                den = (sum(F(a[t] * a[t]) for t in range(4))
                       * sum(F(b[t] * b[t]) for t in range(4)))
                out.append(num / den)
        return sorted(out)

    print(f"   bloque sobre las coordenadas {coords}:")
    for v in bloque:
        print("     ", v)
    g1, g2 = gram(bloque), gram(NUCLEO)
    claim("mismo multiconjunto de cos^2 que el nucleo de B_5 y B_6", g1 == g2,
          f"pares ortogonales: {g1.count(F(0))}")
    print("   Es el mismo nucleo que realiza `reflectivo_rp3_d4b.py` en R^3 por si")
    print("   solo.  Aqui aparece con tres coordenadas de relleno, en R^6.")


def C9():
    header("C9.  LA ELEVACION INGENUA DE B_5 A B_7 ES ESTA MISMA, Y SI ES REFLECTIVA")
    print("   `OCTAVA-RONDA-CANDIDATO-RP6.md` afirmaba que la elevacion natural de la")
    print("   palabra de B_5 a B_7 -- las siete raices rellenadas con ceros, mas e_6 y")
    print("   e_7 -- NO se realiza dinamicamente, tras un barrido de 256 camaras.")
    print("   Es FALSO, y aqui se demuestra: es esta misma configuracion, y su testigo")
    print("   se obtiene transportando el de arriba.")
    B5 = [(0, 1, 0, 0, -1), (0, 0, 1, 0, -1), (0, 0, 0, 1, 0), (1, -1, 0, 0, 0),
          (0, 1, 1, 0, 0), (1, 0, 0, 0, 1), (1, 0, -1, 0, 0)]
    ING = [tuple(list(r) + [0, 0]) for r in B5] + [(0, 0, 0, 0, 0, 1, 0),
                                                   (0, 0, 0, 0, 0, 0, 1)]
    PERM, SGN = (3, 0, 5, 1, 6, 2, 4), (1, 1, 1, -1, 1, -1, -1)
    T = [[F(SGN[k]) if j == PERM[k] else F(0) for j in range(M)] for k in range(M)]
    TT = [[T[k][j] for k in range(M)] for j in range(M)]

    def apl(A, v):
        return [sum((A[i][j] * v[j] for j in range(M)), F(0)) for i in range(M)]

    def canon(v):
        for x in v:
            if x > 0:
                return tuple(v)
            if x < 0:
                return tuple(-y for y in v)
        return tuple(v)

    claim("T es ortogonal (permutacion con signo)",
          all(sum((T[a][c] * T[b][c] for c in range(M)), F(0)) == F(int(a == b))
              for a in range(M) for b in range(M)))
    mapa = {}
    for i, c in enumerate(N):
        for j, g in enumerate(ING):
            if canon(tuple(apl(T, [F(x) for x in g]))) == canon(tuple(F(x) for x in c)):
                mapa[i] = j
    claim("T lleva la palabra ingenua sobre esta, raiz a raiz",
          len(mapa) == K and len(set(mapa.values())) == K, f"biyeccion = {mapa}")

    ELLI = apl(TT, ELL)
    P0X = apl(TT, lift(P0A, True))
    U0X = apl(TT, lift(U0A, False))
    print(f"   carta transportada: ell' = {[str(x) for x in ELLI]}")
    claim("ell'.p = 1", sum((ELLI[a] * P0X[a] for a in range(M)), F(0)) == 1)
    claim("ell'.u = 0", sum((ELLI[a] * U0X[a] for a in range(M)), F(0)) == 0)

    NSI = [[F(x) for x in v] for v in ING]
    RI = [refl(n) for n in NSI]
    P = eye()
    for i in range(K):
        P = mul(RI[mapa[i]], P)
    claim("el producto de la palabra ingenua, en el orden transportado, es -I_7",
          P == [[F(-1) if a == b else F(0) for b in range(M)] for a in range(M)])

    def pasoI(Pp, D, ult):
        cand = []
        for i in range(K):
            if i == ult:
                continue
            den = sum((NSI[i][a] * D[a] for a in range(M)), F(0))
            if den == 0:
                continue
            cand.append((i, -sum((NSI[i][a] * Pp[a] for a in range(M)), F(0)) / den))
        pos = [(i, t) for i, t in cand if t > 0]
        if not pos:
            return None
        b = min(pos, key=lambda c: c[1])
        Q = [Pp[a] + b[1] * D[a] for a in range(M)]
        Dn = [sum((RI[b[0]][a][c] * D[c] for c in range(M)), F(0)) for a in range(M)]
        s = sum((ELLI[a] * Dn[a] for a in range(M)), F(0))
        return b[0], b[1], Q, [Dn[a] - s * Q[a] for a in range(M)]

    Pp, D, ult = P0X, U0X, mapa[CARA0]
    itin, ts = [], []
    for _ in range(K):
        r = pasoI(Pp, D, ult)
        if r is None:
            break
        c, t, Pp, D = r
        itin.append(c)
        ts.append(t)
        ult = c
    esperado = [mapa[i] for i in ITIN]
    claim("la orbita transportada recorre el itinerario correspondiente",
          itin == esperado, f"itinerario = {itin}")
    claim("todos los parametros de paso positivos", all(t > 0 for t in ts))
    claim("y CIERRA exactamente", Pp == P0X)
    kap = next((D[a] / U0X[a] for a in range(M) if U0X[a] != 0), None)
    claim("con kappa > 0", kap is not None and kap > 0, f"kappa = {kap}")
    print()
    print("   CONCLUSION: la elevacion ingenua SI es 9-reflectiva, con la carta")
    print("   ell' de arriba.  El barrido de la octava ronda fijaba un signo -- 256")
    print("   camaras para nueve espejos son 2^8, no 2^9 -- y eso es exactamente el")
    print("   atajo que ya habia ocultado un ejemplo antes (ver README, errores).")
    print("   No habia DOS palabras: habia UNA, vista en dos sistemas de coordenadas.")


def main():
    C1(); C2(); C3(); C4C5(); C6(); C7(); C8(); C9()
    print()
    print("=" * 78)
    if FAILURES:
        print("FALLOS:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("TODO CORRECTO.")
    print()
    print("  Con C1-C6 el criterio de apertura da la 9-REFLECTIVIDAD en R^6.")
    print("  Es el primer ejemplo en dimension 6, con k = 9 IMPAR y n = 6 par,")
    print("  en el umbral k = n+3.  Y C7 y C8 dicen que su algebra es la misma")
    print("  configuracion D_4 de rango 4 de todos los demas, con relleno.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
