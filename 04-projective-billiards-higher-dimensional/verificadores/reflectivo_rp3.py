# -*- coding: utf-8 -*-
"""UN BILLAR PROYECTIVO 6-REFLECTIVO EN R^3.  TODO EXACTO SOBRE Q.

    python reflectivo_rp3.py     ->  codigo de salida 0 solo si TODO pasa

QUE SE AFIRMA

El politopo convexo ACOTADO de R^3 con las seis caras

    cara i  =  { x : e_i (a_i . x + d_i) = 0 },     n_i = (a_i, d_i),

    n_0 = (-1,  3,  -4,  -1)      e_0 = +1
    n_1 = (-1,  0,   4,  +1)      e_1 = +1
    n_2 = ( 3,  4,  -1,  -2)      e_2 = +1
    n_3 = ( 0, -8,  91,   0)      e_3 = +1
    n_4 = ( 0, 2198, 2119, 8345)  e_4 = -1
    n_5 = (-1577, 276, 1117, -146) e_5 = +1

con el campo transversal en la cara i dado por el haz que pasa por el polo O_i =
n_i de esa cara respecto de Q = x0^2+x1^2+x2^2+x3^2, es un billar proyectivo
6-REFLECTIVO: hay un abierto de orbitas 6-periodicas.

POR QUE IMPORTA

Todos los ejemplos conocidos de billares reflectivos proyectivos son PLANOS.
Este es el primero en dimension 3, y cae exactamente en k = n+3 = 6, el umbral
de la tricotomia:  rigido k=n+2=5 (imposible por paridad), umbral k=n+3=6,
generico k>=7.  Ademas su k es PAR, como obliga el teorema de paridad para n
impar.

COMO SE OBTUVO

n_0, n_1, n_2 al azar sobre Z, y n_3,n_4,n_5 por la construccion 3+3 de la
seccion "Higher dimensions": todo elemento de O(4) de determinante -1 tiene
autovalor 1, luego es producto de a lo sumo tres reflexiones, y la factorizacion
es racional por Cartan-Dieudonne.  La orbita testigo se encontro por Monte Carlo
y aqui se re-certifica con datos RACIONALES.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  las seis son reflexiones de O(Q), y el producto en el orden del itinerario
      es -I  (hipotesis H1 del criterio de apertura)
  C2  el politopo es ACOTADO: certificado exacto, coeficientes positivos mu_i con
      suma cero de las normales interiores, y rango 3
  C3  el politopo es NO VACIO: un punto interior racional explicito
  C4  la orbita testigo racional: itinerario 5,4,3,2,1,0 con seis impactos
      interiores a sus caras, cada cara perdedora excluida por desigualdad
      ESTRICTA, transversalidad no nula, y parametros de paso positivos  (H2)
  C5  el cierre:  p_6 = p_0  exactamente, y u_6 = kappa u_0 con kappa > 0  (H3)
  C6  las seis caras son distintas y las normales tienen rango 4 (no degenerado)

Con C1-C6 el criterio de apertura (Teorema "Openness criterion", cuya
demostracion vale verbatim para un politopo convexo de R^n) da la
6-reflectividad.
"""
import sys
from fractions import Fraction as F

M = 4          # dimension homogenea
K = 6          # numero de espejos
FAILURES = []

def claim(name, cond, detail=""):
    ok = bool(cond)
    print(f"   [{'OK  ' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)
    return ok

def header(s):
    print(); print("=" * 78); print(s); print("=" * 78)

# ------------------------------------------------------------ algebra exacta 4x4
def eye(m=M): return [[F(int(i == j)) for j in range(m)] for i in range(m)]
def mul(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))), F(0))
             for j in range(len(B[0]))] for i in range(len(A))]
def smul(c, A): return [[c*x for x in r] for r in A]
def eqm(A, B): return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def mv(A, v): return [sum((A[i][j]*v[j] for j in range(len(v))), F(0)) for i in range(len(A))]

def rank(rows):
    if not rows: return 0
    A = [r[:] for r in rows]
    cols = len(A[0]); r = 0
    for c in range(cols):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = F(1)/A[r][c]; A[r] = [x*inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [A[i][j]-f*A[r][j] for j in range(cols)]
        r += 1
    return r

def solve(rows, rhs, cols):
    """Una solucion particular, o None."""
    A = [rows[i][:] + [rhs[i]] for i in range(len(rows))]
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = F(1)/A[r][c]; A[r] = [x*inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [A[i][j]-f*A[r][j] for j in range(cols+1)]
        piv.append(c); r += 1
    for i in range(r, len(A)):
        if A[i][cols] != 0: return None
    x = [F(0)]*cols
    for i, c in enumerate(piv): x[c] = A[i][cols]
    return x

# ------------------------------------------------------------------- los datos
N = [(-1, 3, -4, -1),
     (-1, 0, 4, 1),
     (3, 4, -1, -2),
     (0, -8, 91, 0),
     (0, 2198, 2119, 8345),
     (-1577, 276, 1117, -146)]
EPS = (1, 1, 1, 1, -1, 1)
NS = [[F(x) for x in n] for n in N]

def refl(nv):
    d = sum((nv[a]*nv[a] for a in range(M)), F(0))
    assert d != 0
    return [[F(int(a == b)) - F(2)*nv[a]*nv[b]/d for b in range(M)] for a in range(M)]

R = [refl(n) for n in NS]
# el semiespacio interior es  e_i (a_i . x + d_i) < 0
A3 = [[F(EPS[i])*NS[i][j] for j in range(3)] for i in range(K)]
D3 = [F(EPS[i])*NS[i][3] for i in range(K)]

def lado(i, x):
    return sum((A3[i][j]*x[j] for j in range(3)), F(0)) + D3[i]

# ---------------------------------------------------------------------- C1
def C1():
    header("C1.  SEIS REFLEXIONES DE O(Q), Y EL PRODUCTO ES -I")
    for i in range(K):
        claim(f"R_{i}^2 = I", eqm(mul(R[i], R[i]), eye()))
        claim(f"tr R_{i} = 2  (reflexion, no -I)",
              sum((R[i][a][a] for a in range(M)), F(0)) == F(2))
        Rt = [list(c) for c in zip(*R[i])]
        claim(f"R_{i} en O(Q), Q = I", eqm(mul(Rt, R[i]), eye()))
    print("   El itinerario de la orbita testigo es 5,4,3,2,1,0, luego el producto")
    print("   relevante, en la convencion de puntos, es R_0 R_1 R_2 R_3 R_4 R_5.")
    P = eye()
    for i in range(K):
        P = mul(P, R[i])                      # R_0 R_1 ... R_5
    claim("R_0 R_1 R_2 R_3 R_4 R_5 = -I", eqm(P, smul(F(-1), eye())))
    P2 = eye()
    for i in range(K-1, -1, -1):
        P2 = mul(P2, R[i])                    # R_5 R_4 ... R_0
    claim("R_5 R_4 R_3 R_2 R_1 R_0 = -I  (el otro orden, calculado aparte)",
          eqm(P2, smul(F(-1), eye())))
    claim("las seis normales tienen rango 4: la configuracion no es degenerada",
          rank([n[:] for n in NS]) == 4, f"rango = {rank([n[:] for n in NS])}")
    caras = set()
    for i in range(K):
        v = NS[i]
        g = None
        for x in v:
            if x != 0:
                g = x; break
        caras.add(tuple(x/g for x in v))
    claim("las seis caras son distintas", len(caras) == K, f"{len(caras)} distintas")

# ---------------------------------------------------------------------- C2, C3
def C2_C3():
    header("C2.  EL POLITOPO ES ACOTADO: CERTIFICADO EXACTO")
    print("   El cono de recesion es { y : a_i . y <= 0 para todo i }.  Es {0} si y")
    print("   solo si las a_i generan POSITIVAMENTE R^3, y eso equivale a que")
    print("   existan mu_i > 0 con suma mu_i a_i = 0 y que las a_i tengan rango 3.")
    print("   Se busca ese mu resolviendo el sistema lineal sobre Q.")
    claim("las a_i tienen rango 3", rank([r[:] for r in A3]) == 3,
          f"rango = {rank([r[:] for r in A3])}")
    # buscar mu > 0 con sum mu_i a_i = 0.  El nucleo de la 3x6 tiene dim 3.
    rows = [[A3[i][j] for i in range(K)] for j in range(3)]   # 3 x 6
    # base del nucleo
    A = [r[:] for r in rows]
    piv, r = [], 0
    for c in range(K):
        p = next((i for i in range(r, 3) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        inv = F(1)/A[r][c]; A[r] = [x*inv for x in A[r]]
        for i in range(3):
            if i != r and A[i][c] != 0:
                f = A[i][c]; A[i] = [A[i][j]-f*A[r][j] for j in range(K)]
        piv.append(c); r += 1
    free = [c for c in range(K) if c not in piv]
    base = []
    for fc in free:
        v = [F(0)]*K; v[fc] = F(1)
        for i, c in enumerate(piv): v[c] = -A[i][fc]
        base.append(v)
    print(f"   dim del nucleo = {len(base)}")
    mu = None
    rng = range(-6, 7)
    import itertools
    for coef in itertools.product(rng, repeat=len(base)):
        if all(c == 0 for c in coef): continue
        v = [sum((F(coef[t])*base[t][i] for t in range(len(base))), F(0)) for i in range(K)]
        if all(x > 0 for x in v):
            mu = v; break
    claim("existe mu > 0 con suma mu_i a_i = 0", mu is not None,
          f"mu = {[str(x) for x in mu]}" if mu else "")
    if mu is not None:
        chk = [sum((mu[i]*A3[i][j] for i in range(K)), F(0)) for j in range(3)]
        claim("y en efecto suma mu_i a_i = 0 exactamente", all(x == 0 for x in chk))
        print("   Con rango 3 y mu > 0, el cono de recesion es {0}: ACOTADO.")

    header("C3.  EL POLITOPO ES NO VACIO: UN PUNTO INTERIOR RACIONAL")
    print("   Se toma el punto p_0 de la orbita testigo, que esta EN la cara 0 y es")
    print("   estrictamente interior respecto de las otras cinco, y se desplaza")
    print("   hacia dentro a lo largo de -a_0:  x0 = p_0 - t a_0  con t = 1/100.")
    p0 = [F(3)*F(-3, 2) - F(4)*F(-5, 4) - 1, F(-3, 2), F(-5, 4)]
    t = F(1, 100)
    x0 = [p0[j] - t*A3[0][j] for j in range(3)]
    print(f"   p_0 = ({p0[0]}, {p0[1]}, {p0[2]})")
    print(f"   x0  = ({x0[0]}, {x0[1]}, {x0[2]})")
    for i in range(K):
        v = lado(i, x0)
        claim(f"lado {i}: e_i(a_i.x0 + d_i) < 0", v < 0, f"= {v}")

# ---------------------------------------------------------------------- C4, C5
def C4_C5():
    header("C4.  LA ORBITA TESTIGO, RACIONAL Y EXACTA")
    print("   Se parte de un punto racional de la cara 0 y una direccion racional.")
    print("   La cara 0 es  -x0 + 3x1 - 4x2 - 1 = 0,  luego x0 = 3x1 - 4x2 - 1.")
    x1, x2 = F(-3, 2), F(-5, 4)
    p = [3*x1 - 4*x2 - 1, x1, x2]
    claim("p_0 esta en la cara 0", lado(0, p) == 0, f"lado 0 = {lado(0,p)}")
    for i in range(1, K):
        v = lado(i, p)
        claim(f"p_0 es interior respecto de la cara {i} (estricto)", v < 0, f"= {v}")
    u = [F(-4, 9), F(-7, 10), F(5, 9)]
    print(f"   p_0 = ({p[0]}, {p[1]}, {p[2]})")
    print(f"   u_0 = ({u[0]}, {u[1]}, {u[2]})")
    claim("u_0 apunta hacia dentro del semiespacio de la cara 0",
          sum((A3[0][j]*u[j] for j in range(3)), F(0)) < 0,
          f"a_0.u_0 = {sum((A3[0][j]*u[j] for j in range(3)), F(0))}")

    itin_esperado = [5, 4, 3, 2, 1, 0]
    cara = 0
    p_act, u_act = p[:], u[:]
    u_ini = u[:]
    for step in range(K):
        cands = []
        for i in range(K):
            if i == cara: continue
            den = sum((A3[i][j]*u_act[j] for j in range(3)), F(0))
            if den == 0: continue
            s = -lado(i, p_act)/den
            if s > 0:
                q = [p_act[j] + s*u_act[j] for j in range(3)]
                cands.append((s, i, q))
        claim(f"paso {step+1}: toda cara no actual es no paralela al rayo",
              len(cands) + 1 >= 1 and all(
                  sum((A3[i][j]*u_act[j] for j in range(3)), F(0)) != 0
                  for i in range(K) if i != cara))
        claim(f"paso {step+1}: hay al menos un impacto con s > 0", len(cands) > 0)
        if not cands: return
        cands.sort(key=lambda t: t[0])
        s, i, q = cands[0]
        esperada = itin_esperado[step]
        claim(f"paso {step+1}: el primer impacto (s minimo) es la cara {esperada}",
              i == esperada, f"cara {i}, s = {s}")
        claim(f"paso {step+1}: s > 0 estricto", s > 0)
        estricto = True
        for (s2, j2, q2) in cands[1:]:
            if s2 <= s: estricto = False
        claim(f"paso {step+1}: las caras perdedoras quedan excluidas por"
              f" desigualdad ESTRICTA", estricto)
        interior = all(lado(j, q) < 0 for j in range(K) if j != i)
        claim(f"paso {step+1}: el impacto es INTERIOR a la cara {i}", interior)
        claim(f"paso {step+1}: y esta en la cara {i}", lado(i, q) == 0)
        # transversalidad: el centro O_i no esta en la cara i
        tr = sum((NS[i][j]*NS[i][j] for j in range(M)), F(0))
        claim(f"paso {step+1}: transversalidad no nula (L_i O_i != 0)", tr != 0)
        # reflexion:  (w, w4) = H_i (u, 0);  u' = w - w4 q
        ue = [u_act[0], u_act[1], u_act[2], F(0)]
        wv = mv(R[i], ue)
        u_new = [wv[j] - wv[3]*q[j] for j in range(3)]
        claim(f"paso {step+1}: la direccion reflejada no es nula",
              any(x != 0 for x in u_new))
        p_act, u_act, cara = q, u_new, i

    header("C5.  EL CIERRE")
    claim("p_6 = p_0 EXACTAMENTE", p_act == p, f"p_6 = {[str(x) for x in p_act]}")
    # u_6 = kappa u_0 con kappa > 0
    kappa = None
    for j in range(3):
        if u_ini[j] != 0:
            kappa = u_act[j]/u_ini[j]; break
    prop = all(u_act[j]*u_ini[jj] - u_act[jj]*u_ini[j] == 0
               for j in range(3) for jj in range(3))
    claim("u_6 es proporcional a u_0", prop)
    claim("kappa > 0  (la orientacion se conserva)", kappa is not None and kappa > 0,
          f"kappa = {kappa}")
    claim("los seis impactos estan en seis caras distintas: periodo primitivo 6",
          len(set(itin_esperado)) == K)

def main():
    print(__doc__)
    C1()
    C2_C3()
    C4_C5()
    header("VEREDICTO")
    if FAILURES:
        print(f"   {len(FAILURES)} COMPROBACIONES FALLIDAS:")
        for f in FAILURES:
            print(f"     - {f}")
        sys.exit(1)
    print("   Todas las comprobaciones pasan, en aritmetica exacta sobre Q.")
    print("   Sin coma flotante en ninguna decision.")
    print()
    print("   Con H1 (C1), el politopo convexo acotado (C2, C3) y la orbita")
    print("   testigo con todas sus desigualdades estrictas (C4, C5), el criterio")
    print("   de apertura da un ABIERTO de orbitas 6-periodicas:")
    print()
    print("      EXISTE UN BILLAR PROYECTIVO 6-REFLECTIVO EN R^3.")
    sys.exit(0)

if __name__ == "__main__":
    main()
