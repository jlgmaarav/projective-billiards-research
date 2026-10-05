# -*- coding: utf-8 -*-
"""UN TERCER BILLAR 6-REFLECTIVO EN R^3.  TODO EXACTO SOBRE Q.

    python reflectivo_rp3_d4b.py     ->  codigo de salida 0 solo si TODO pasa

QUE SE AFIRMA

El politopo convexo ACOTADO de R^3 obtenido cortando las seis caras

    n_0 = ( 0, 1, 0,-1)   n_1 = ( 0, 0, 1,-1)   n_2 = ( 1,-1, 0, 0)
    n_3 = ( 0, 1, 1, 0)   n_4 = ( 1, 0, 0, 1)   n_5 = ( 1, 0,-1, 0)

con el corte afin  ell = (-2,-2,-2,-1)  y el campo transversal dado por el polo
O_i = n_i respecto de Q = x0^2+x1^2+x2^2+x3^2, es un billar proyectivo
6-REFLECTIVO: hay un abierto de orbitas 6-periodicas.

DE DONDE SALE, Y POR QUE IMPORTA

Las seis raices son de tipo e_i +- e_j, o sea raices de D_4, y son EXACTAMENTE el
nucleo de rango 4 al que se reducen las palabras de B_5 y B_6 de `thm:rp4b5` y
`thm:rp5b6` al quitarles las coordenadas aisladas.  Es decir: la "familia
algebraica uniforme" de `SEPTIMA-RONDA.md` es esta configuracion mas relleno, en
todos los rangos.  Aqui se demuestra que ese nucleo no es solo algebra: **se
realiza dinamicamente por si mismo**, en R^3 y en el umbral k = n+3 = 6.

Y NO es el caleidoscopio D_4 ya certificado en `reflectivo_rp3_d4.py`.  Se
comprueba abajo con un invariante ortogonal: el multiconjunto de los cos^2 de los
angulos entre normales.  Aquel tiene seis pares ortogonales; este, tres.  Luego
ninguna transformacion ortogonal lleva una configuracion en la otra.

Con esto, los ejemplos certificados en R^3 pasan a ser TRES:
  - `thm:rp3example`, grupo de espejos infinito;
  - `thm:rp3d4`, caleidoscopio D_4;
  - este, tambien D_4 y no equivalente al anterior.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  las seis son reflexiones de O(Q) y el producto en el orden del itinerario
      es -I  (hipotesis H1 del criterio de apertura)
  C2  el politopo es ACOTADO: coeficientes racionales mu_i > 0 con suma cero de
      las normales exteriores, y rango 3
  C3  es NO VACIO: un punto interior racional explicito, con margen estricto
  C4  la orbita testigo racional: seis impactos, cada uno estrictamente dentro de
      su cara, cada cara perdedora excluida por desigualdad ESTRICTA, y todos los
      parametros de paso positivos  (H2)
  C5  el cierre: p_6 = p_0 EXACTAMENTE, y u_6 = kappa u_0 con kappa > 0  (H3)
  C6  las seis caras son distintas, las normales tienen rango 4, y la
      configuracion es IRREDUCIBLE (corango 0), como pide el Corolario B'
  C7  no es ortogonalmente equivalente al caleidoscopio de `thm:rp3d4`
"""
import sys
from fractions import Fraction as F

M, K = 4, 6
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
N = [(0, 1, 0, -1), (0, 0, 1, -1), (1, -1, 0, 0),
     (0, 1, 1, 0), (1, 0, 0, 1), (1, 0, -1, 0)]
ELL = [F(-2), F(-2), F(-2), F(-1)]
LIBRE = 0
EPS = [1] * K                      # el semiespacio interior es eps_i (a_i.x + d_i) < 0
POS = [c for c in range(M) if c != LIBRE]

NS = [[F(x) for x in n] for n in N]


def refl(nv):
    d = sum((nv[a] * nv[a] for a in range(M)), F(0))
    assert d != 0
    return [[F(int(a == b)) - F(2) * nv[a] * nv[b] / d for b in range(M)] for a in range(M)]


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


# ---- la carta afin:  n.X = a.x + d,  donde x son las coordenadas != LIBRE
def carta(n):
    n0, l0 = n[LIBRE], ELL[LIBRE]
    return [n[c] - n0 * ELL[c] / l0 for c in POS], n0 / l0


CH = [carta(n) for n in NS]
A3 = [a for a, d in CH]
D3 = [d for a, d in CH]


def lado(i, x):
    return sum((A3[i][j] * x[j] for j in range(3)), F(0)) + D3[i]


def lift(v, punto):
    X = [F(0)] * M
    for a, c in enumerate(POS):
        X[c] = v[a]
    X[LIBRE] = ((F(1) if punto else F(0))
                - sum((ELL[c] * X[c] for c in POS), F(0))) / ELL[LIBRE]
    return X


# ---------------------------------------------------------------------- C1
def C1():
    header("C1.  SEIS REFLEXIONES DE O(Q), Y EL PRODUCTO ES -I")
    for i in range(K):
        claim(f"R_{i}^2 = I", mul(R[i], R[i]) == eye())
        claim(f"tr R_{i} = 2  (reflexion, no -I)",
              sum((R[i][a][a] for a in range(M)), F(0)) == F(2))
    P = eye()
    for i in range(K):
        P = mul(R[i], P)
    claim("R_5 R_4 R_3 R_2 R_1 R_0 = -I_4", P == [[F(-1) if a == b else F(0)
                                                   for b in range(M)] for a in range(M)])
    claim("las seis raices son de tipo e_i +- e_j  (o sea, de D_4)",
          all(sorted(abs(x) for x in n if x != 0) == [1, 1] for n in N))


# ---------------------------------------------------------------------- C2
def C2():
    header("C2.  EL POLITOPO ES ACOTADO, CON CERTIFICADO RACIONAL")
    out = [[F(EPS[i]) * A3[i][j] for j in range(3)] for i in range(K)]
    claim("las normales exteriores tienen rango 3", rank(out) == 3)
    # buscar mu > 0 racional con sum mu_i out_i = 0, por barrido entero pequeno
    hallado = None
    rng = range(1, 9)
    for m0 in rng:
        for m1 in rng:
            for m2 in rng:
                for m3 in rng:
                    # las dos ultimas se despejan de dos ecuaciones y se comprueba la tercera
                    a = [[out[4][j], out[5][j]] for j in range(3)]
                    b = [-(F(m0) * out[0][j] + F(m1) * out[1][j]
                           + F(m2) * out[2][j] + F(m3) * out[3][j]) for j in range(3)]
                    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
                    if det == 0:
                        continue
                    m4 = (b[0] * a[1][1] - b[1] * a[0][1]) / det
                    m5 = (a[0][0] * b[1] - a[1][0] * b[0]) / det
                    if m4 <= 0 or m5 <= 0:
                        continue
                    if a[2][0] * m4 + a[2][1] * m5 == b[2]:
                        hallado = [F(m0), F(m1), F(m2), F(m3), m4, m5]
                        break
                if hallado:
                    break
            if hallado:
                break
        if hallado:
            break
    claim("existe mu racional ESTRICTAMENTE positivo con sum mu_i n_i^ext = 0",
          hallado is not None, f"mu = {[str(x) for x in hallado]}" if hallado else "")
    if hallado:
        for j in range(3):
            claim(f"   suma nula, coordenada {j}",
                  sum((hallado[i] * out[i][j] for i in range(K)), F(0)) == 0)


# ---------------------------------------------------------------------- C3
P0 = [F(-1, 9), F(-1, 9), F(1, 9)]


def C3():
    header("C3.  EL POLITOPO ES NO VACIO: PUNTO INTERIOR RACIONAL")
    print(f"   p = {[str(x) for x in P0]}")
    for i in range(K):
        v = F(EPS[i]) * lado(i, P0)
        claim(f"cara {i}: eps_i (a_i . p + d_i) < 0 ESTRICTO", v < 0, f"= {v}")


# ------------------------------------------------------------------ C4 y C5
U0 = [F(5, 4), F(-1, 4), F(2)]
ITIN = [3, 4, 5, 0, 1, 2]


def paso(P, D, ult, verbose=False):
    cand = []
    for i in range(K):
        if i == ult:
            continue
        den = sum((NS[i][a] * D[a] for a in range(M)), F(0))
        if den == 0:
            continue
        t = -sum((NS[i][a] * P[a] for a in range(M)), F(0)) / den
        cand.append((i, t))
    pos = [(i, t) for i, t in cand if t > 0]
    if not pos:
        return None
    best = min(pos, key=lambda c: c[1])
    if verbose:
        print("      t de cada cara: " + ", ".join(f"{i}:{t}" for i, t in cand))
    Q = [P[a] + best[1] * D[a] for a in range(M)]
    Dn = [sum((R[best[0]][a][b] * D[b] for b in range(M)), F(0)) for a in range(M)]
    s = sum((ELL[a] * Dn[a] for a in range(M)), F(0))
    Dn = [Dn[a] - s * Q[a] for a in range(M)]
    return best[0], best[1], Q, Dn, cand


def C4C5():
    header("C4 y C5.  ORBITA TESTIGO RACIONAL, SEIS IMPACTOS, Y EL CIERRE")
    print(f"   se lanza desde el interior p = {[str(x) for x in P0]}")
    print(f"   con direccion u = {[str(x) for x in U0]};  el primer impacto define")
    print("   el punto de partida de la orbita cerrada.")
    P, D = lift(P0, True), lift(U0, False)
    r = paso(P, D, -1)
    claim("hay primer impacto", r is not None)
    if r is None:
        return
    cara0, t0, Q0, D1, _ = r
    claim("el primer impacto cae en la cara 3", cara0 == 3, f"t = {t0}")
    P, D, ult = Q0, D1, cara0
    itin, ts = [], []
    for paso_n in range(K):
        rr = paso(P, D, ult, verbose=False)
        claim(f"paso {paso_n+1}: hay impacto", rr is not None)
        if rr is None:
            return
        c, t, P, D, cand = rr
        itin.append(c)
        ts.append(t)
        claim(f"paso {paso_n+1}: la cara ganadora es {ITIN[(paso_n+1) % K]}",
              c == ITIN[(paso_n + 1) % K], f"t = {t}")
        claim(f"paso {paso_n+1}: parametro de paso ESTRICTAMENTE positivo", t > 0)
        otras = [(i, tt) for i, tt in cand if i != c]
        claim(f"paso {paso_n+1}: las demas caras quedan excluidas por desigualdad "
              f"ESTRICTA", all(tt <= 0 or tt > t for i, tt in otras))
        # el impacto tiene que estar estrictamente dentro de su cara
        xaf = [P[c2] / sum((ELL[a] * P[a] for a in range(M)), F(0)) for c2 in POS]
        dentro = all(F(EPS[i]) * lado(i, xaf) < 0 for i in range(K) if i != c)
        claim(f"paso {paso_n+1}: el impacto es interior a su cara", dentro)
        ult = c
    claim("el itinerario recorre las seis caras, sin repetir",
          len(set(itin)) == K, f"itinerario = {itin}")
    claim("CIERRE: p_6 = p_0 EXACTAMENTE", P == Q0)
    kap = None
    for a in range(M):
        if D1[a] != 0:
            kap = D[a] / D1[a]
            break
    claim("u_6 = kappa u_0", all(D[a] == kap * D1[a] for a in range(M)))
    claim("y kappa > 0: misma orientacion", kap is not None and kap > 0, f"kappa = {kap}")


# ---------------------------------------------------------------------- C6
def C6():
    header("C6.  NO DEGENERADA, Y IRREDUCIBLE")
    claim("las seis caras son distintas",
          len({tuple(n) for n in N}) == K)
    claim("las normales tienen rango 4", rank([list(n) for n in NS]) == M)
    # irreducibilidad: {S simetrica de traza nula : S n_i // n_i} = 0
    basis = [(i, j) for i in range(M) for j in range(i, M)]

    def comp(i, j, a, r):
        v = F(0)
        if i == r:
            v += a[j]
        if j == r and i != j:
            v += a[i]
        return v

    rows = []
    for a in NS:
        for p in range(M):
            for q in range(p + 1, M):
                rows.append([a[p] * comp(i, j, a, q) - a[q] * comp(i, j, a, p)
                             for i, j in basis])
    rows.append([F(1) if i == j else F(0) for i, j in basis])
    claim("corango 0: la configuracion es IRREDUCIBLE  [Corolario B']",
          len(basis) - rank(rows) == 0, f"dim del conucleo = {len(basis) - rank(rows)}")


# ---------------------------------------------------------------------- C7
def C7():
    header("C7.  NO ES EL CALEIDOSCOPIO D_4 YA CERTIFICADO")
    OTRO = [(1, 0, 1, 0), (1, 0, -1, 0), (1, 0, 0, -1),
            (0, 1, 1, 0), (0, 1, -1, 0), (0, 1, 0, -1)]

    def gram(S):
        out = []
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                a, b = S[i], S[j]
                num = sum(F(a[t] * b[t]) for t in range(M)) ** 2
                den = (sum(F(a[t] * a[t]) for t in range(M))
                       * sum(F(b[t] * b[t]) for t in range(M)))
                out.append(num / den)
        return sorted(out)

    g1, g2 = gram(N), gram(OTRO)
    print(f"   este:   {[str(x) for x in g1]}")
    print(f"   aquel:  {[str(x) for x in g2]}")
    claim("los multiconjuntos de cos^2 entre normales DIFIEREN", g1 != g2)
    claim("luego ninguna transformacion ortogonal lleva uno en el otro",
          g1.count(F(0)) != g2.count(F(0)),
          f"pares ortogonales: {g1.count(F(0))} frente a {g2.count(F(0))}")
    print()
    print("   Y de ahi la inequivalencia PROYECTIVA, que es lo que hace falta:")
    print("   las dos configuraciones son irreducibles (C6 aqui, y tolerancias.py")
    print("   para la otra), luego cada una preserva UNA SOLA cuadrica salvo escala.")
    print("   Una equivalencia proyectiva entre ellas tendria que llevar la cuadrica")
    print("   invariante de una en la de la otra, o sea Q en un multiplo de Q, luego")
    print("   seria ortogonal salvo escala, y conservaria el multiconjunto de cos^2.")
    print("   Como no lo conserva, no hay tal equivalencia.")


def C8():
    header("C8.  EL GRUPO DE ESPEJOS: ES UN CALEIDOSCOPIO, Y CONVIENE DECIRLO")

    def orden(S):
        gens = [refl([F(x) for x in n]) for n in S]
        I = tuple(tuple(row) for row in eye())
        G, frontera = {I}, [I]
        while frontera:
            nueva = []
            for g in frontera:
                for s in gens:
                    h = tuple(tuple(r) for r in mul([list(r) for r in s],
                                                    [list(r) for r in g]))
                    if h not in G:
                        G.add(h)
                        nueva.append(h)
            frontera = nueva
            if len(G) > 1000:
                return None
        return len(G)

    o = orden(N)
    claim("el grupo de espejos es FINITO", o is not None, f"orden = {o}")
    claim("y es todo W(D_4), de orden 192", o == 192)
    print()
    print("   Luego este ejemplo es CALEIDOSCOPICO, igual que el de `thm:rp3d4`,")
    print("   que tambien genera un grupo de orden 192.  No es un ejemplo de grupo")
    print("   infinito: el unico de esos en R^3 sigue siendo `thm:rp3example`.")
    print("   Lo nuevo aqui es que son DOS caleidoscopios inequivalentes, y que")
    print("   este es exactamente el nucleo de las palabras de B_5 y B_6.")
    print()
    print("   Consecuencia espectral, y va en contra del interes de este ejemplo:")
    print("   al desdoblarse a un grupo finito, TODA orbita es periodica, asi que")
    print("   -- como el pentagono simetrico del plano -- no sirve para el")
    print("   experimento del espectro de longitudes.  Ver PREDICCION-ESPECTRAL.md.")


def main():
    C1(); C2(); C3(); C4C5(); C6(); C7(); C8()
    print()
    print("=" * 78)
    if FAILURES:
        print("FALLOS:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("TODO CORRECTO.")
    print()
    print("  Con C1-C6 el criterio de apertura da la 6-reflectividad, y C7 dice que")
    print("  es un ejemplo NUEVO: el tercero certificado en R^3, y el primero que")
    print("  demuestra que el nucleo D_4 de las palabras de B_5 y B_6 se realiza")
    print("  dinamicamente por si solo, sin el relleno de coordenadas aisladas.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
