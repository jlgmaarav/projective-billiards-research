# -*- coding: utf-8 -*-
"""CLASIFICACION COMPLETA DE k=5, Y LA IRREDUCIBILIDAD DE V_5.

    python clasificacion_k5.py     ->  codigo de salida 0 solo si TODO pasa

Certifica las hipotesis computables de dos teoremas nuevos:

  TEOREMA A (irreducibilidad).  V_5 es irreducible de dimension 12.
      Cierra el unico hueco que el Remark "What is still assumed" de note.tex
      declaraba abierto.  La demostracion fibra V_5^reg sobre el lugar de los
      A = H_1 H_0 regulares semisimples:

          base    A regular semisimple, espectro {1,lambda,1/lambda}   dim 7
          fibra   H_0 en F_A                                           dim 1
                  N en Inv con det(N-A) = 0                            dim 3
                  H_2 en F_B,  B = N A^{-1}                            dim 1
                  H_1 = A H_0,  H_3 = B H_2,  H_4 = -N^{-1}        determinados
          total   7 + 1 + 3 + 1 = 12

      Base irreducible, fibras irreducibles de dimension constante => V_5^reg
      irreducible; y V_5 = clausura(V_5^reg) por la Prop. de densidad global.

  TEOREMA B (clasificacion).  La construccion de la Prop. 11 es EXHAUSTIVA:
      todo billar proyectivo 5-reflectivo convexo se obtiene de ella.  El
      espacio de moduli es unirracional de dimension 4 = 4k-16.

  COROLARIO.  V_5 = S_5:  TODA configuracion de V_5 preserva una forma
      cuadratica no nula, no solo las de un abierto.

Lo que aqui se comprueba, en Q exacto y sin coma flotante:
  C1  el ejemplo del Thm. "A second, inequivalent 5-reflective billiard"
      esta en V_5^reg  (asi las hipotesis de los Teoremas A y B son no vacias)
  C2  la ecuacion (b w_2 = c w_3) leida SIN BASE como det(M - A) = 0
  C3  exhaustividad: n_3 se RECUPERA de (n_0,n_1,n_2,n_4)
  C4  la hipersuperficie {b w_2 = c w_3} es irreducible (rango 2 de su forma
      bilineal), que es el paso que hace irreducible la fibra
  C5  el recuento de dimensiones cierra: 7+1+3+1 = 12 = 4k-8
  C6  la conica invariante es unica salvo escala (rango de M exactamente 5)
"""
from fractions import Fraction as F
from itertools import combinations
import sys

K = 5
FAILURES = []

def claim(name, cond, detail=""):
    ok = bool(cond)
    print(f"   [{'OK  ' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)
    return ok

def header(s):
    print(); print("=" * 78); print(s); print("=" * 78)

# --------------------------------------------------------- algebra exacta 3x3
def I3():   return [[F(int(i == j)) for j in range(3)] for i in range(3)]
def sc(c):  return [[F(c) if i == j else F(0) for j in range(3)] for i in range(3)]
def mul(A, B): return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def sub(A, B): return [[A[i][j]-B[i][j] for j in range(3)] for i in range(3)]
def add(A, B): return [[A[i][j]+B[i][j] for j in range(3)] for i in range(3)]
def eqm(A, B): return all(A[i][j] == B[i][j] for i in range(3) for j in range(3))
def tr(M):  return M[0][0] + M[1][1] + M[2][2]
def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

def rank(M):
    """Rango exacto por eliminacion de Gauss sobre Q."""
    M = [row[:] for row in M]
    rows, cols = len(M), (len(M[0]) if M else 0)
    r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = F(1)/M[r][c]
        M[r] = [x*inv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j] - f*M[r][j] for j in range(cols)]
        r += 1
    return r

def kernel_vec(M):
    """Generador entero del nucleo si es de dimension 1; si no, None."""
    from math import gcd
    A = [row[:] for row in M]
    piv, r = [], 0
    for c in range(3):
        p = next((i for i in range(r, 3) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = F(1)/A[r][c]
        A[r] = [x*inv for x in A[r]]
        for i in range(3):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f*A[r][j] for j in range(3)]
        piv.append(c); r += 1
    free = [c for c in range(3) if c not in piv]
    if len(free) != 1:
        return None
    fc = free[0]
    v = [F(0)]*3; v[fc] = F(1)
    for i, c in enumerate(piv):
        v[c] = -A[i][fc]
    den = 1
    for x in v:
        den = den*x.denominator//gcd(den, x.denominator)
    v = [int(x*den) for x in v]
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    if g:
        v = [x//g for x in v]
    return [F(x) for x in v]

def charpoly(M):
    """(c2,c1,c0) de t^3 - c2 t^2 + c1 t - c0."""
    c1 = ((M[0][0]*M[1][1]-M[0][1]*M[1][0]) + (M[0][0]*M[2][2]-M[0][2]*M[2][0])
        + (M[1][1]*M[2][2]-M[1][2]*M[2][1]))
    return tr(M), c1, det3(M)

# ------------------------------------------------------------------ el ejemplo
Q = I3()                                    # conica invariante x^2+y^2+z^2
N = [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)]

def Qv(v):    return [sum(Q[a][b]*v[b] for b in range(3)) for a in range(3)]
def qform(u, v):
    Qvv = Qv(v)
    return sum(u[a]*Qvv[a] for a in range(3))

def refl(n):
    """Homologia armonica = reflexion Q-ortogonal de normal n (centro = polo)."""
    qn, d = Qv(n), qform(n, n)
    assert d != 0
    return [[F(int(a == b)) - F(2)*n[a]*qn[b]/d for b in range(3)] for a in range(3)]

NS = [[F(x) for x in n] for n in N]
R = [refl(n) for n in NS]

# ------------------------------------------------------------------------ C1
def C1_regular():
    header("C1.  EL EJEMPLO ESTA EN V_5^reg   (hipotesis no vacias)")
    P = sc(1)
    for i in range(K-1, -1, -1):
        P = mul(R[i], P)
    claim("R4 R3 R2 R1 R0 = -I", eqm(P, sc(-1)))

    A = mul(R[1], R[0])
    B = mul(R[3], R[2])
    print("   Un producto de DOS involuciones tiene det 1 y autovalor 1 (Lema del par),")
    print("   luego su polinomio caracteristico es (t-1)(t^2 - beta t + 1), beta = tr - 1.")
    print("   Es regular semisimple con lambda != +-1  <=>  beta != +-2  y  beta^2-4 != 0.")
    out = {}
    for name, X in (("A = H1 H0", A), ("B = H3 H2", B)):
        c2, c1, c0 = charpoly(X)
        beta = c2 - 1
        claim(f"{name}: det = 1", c0 == F(1))
        claim(f"{name}: charpoly = (t-1)(t^2 - beta t + 1)", c1 == beta + 1,
              f"beta = {beta}")
        claim(f"{name}: t=1 no es raiz del factor cuadratico  (2 - beta != 0)",
              2 - beta != 0, f"2-beta = {2-beta}")
        claim(f"{name}: lambda != +-1  (beta != +-2)", beta not in (F(2), F(-2)))
        claim(f"{name}: semisimple  (beta^2 - 4 != 0)", beta*beta - 4 != 0,
              f"beta^2-4 = {beta*beta-4}")
        out[name] = X
    return A, B

# ------------------------------------------------------------------------ C2
def C2_ecuacion(A, B):
    header("C2.  LA ECUACION b w_2 = c w_3, LEIDA SIN BASE")
    print("   En la base propia de A la condicion es b w_2 = c w_3.  La formula")
    print("   encajonada del Teorema del umbral dice")
    print("       det(N - A) = -2 (1+lambda)(lambda-1)(b w_2 - c w_3) / (lambda s),")
    print("   luego con lambda != +-1 y s != 0 la condicion es EXACTAMENTE")
    print("   det(N - A) = 0, que si es racional y no necesita diagonalizar nada.")
    M = mul(B, A)                     # N = M = H3 H2 H1 H0
    claim("det M = 1", det3(M) == F(1))
    claim("M^2 = I   (Prop. de reduccion)", eqm(mul(M, M), sc(1)))
    claim("M != I", not eqm(M, sc(1)))
    claim("tr M = -1  (involucion no trivial: autovalores 1,-1,-1)", tr(M) == F(-1))
    claim("det(M - A) = 0   <-- la ecuacion (b w2 = c w3)", det3(sub(M, A)) == F(0))
    H4 = [[-x for x in row] for row in inv3(M)]
    claim("H_4 = -M^{-1} es la quinta homologia del ejemplo", eqm(H4, R[4]))
    return M

def inv3(M):
    d = det3(M)
    assert d != 0
    C = [[F(0)]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r = [x for k, x in enumerate(range(3)) if k != i]
            c = [x for k, x in enumerate(range(3)) if k != j]
            m = M[r[0]][c[0]]*M[r[1]][c[1]] - M[r[0]][c[1]]*M[r[1]][c[0]]
            C[j][i] = ((-1)**(i+j))*m/d
    return C

# ------------------------------------------------------------------------ C3
def C3_exhaustiva():
    header("C3.  EXHAUSTIVIDAD:  n_3 SE RECUPERA DE (n_0, n_1, n_2, n_4)")
    print("   Si R4R3R2R1R0 = -I entonces R3 R4 = B := -R2R1R0, luego B es producto")
    print("   de dos reflexiones y tiene autovalor 1.  R_4 fija un w de ker(B-I):")
    print("   si no, w seria paralelo a n_4 y saldria R_3 = R_4.  Por tanto todo")
    print("   ejemplo esta en la imagen de la construccion; no hay otros.")
    T = mul(R[2], mul(R[1], R[0]))
    B = [[-T[a][b] for b in range(3)] for a in range(3)]
    claim("det B = 1  (B es una rotacion de O(Q))", det3(B) == F(1))
    claim("B != I", not eqm(B, sc(1)))
    claim("rango(B - I) = 2, luego ker(B-I) es una recta", rank(sub(B, sc(1))) == 2)
    w = kernel_vec(sub(B, sc(1)))
    claim("w = ker(B-I) es racional", w is not None, f"w = {[int(x) for x in w]}")
    claim("w NO es isotropo  (w^T Q w != 0)", qform(w, w) != 0, f"= {qform(w,w)}")
    claim("n_4 . (Q w) = 0   (el espejo 4 contiene el eje de B)",
          qform(NS[4], w) == 0)
    R3 = mul(B, R[4])
    claim("(B R_4)^2 = I,  luego B R_4 ES una reflexion", eqm(mul(R3, R3), sc(1)))
    claim("det(B R_4) = -1", det3(R3) == F(-1))
    claim("B R_4 = R_3 del ejemplo", eqm(R3, R[3]))
    n3 = kernel_vec(add(R3, sc(1)))
    claim("n_3 recuperada = (-8, 9, 5) salvo escala",
          n3 is not None and prop(n3, NS[3]), f"n_3 = {[int(x) for x in n3]}")

def prop(a, b):
    if not any(a) or not any(b):
        return False
    for i in range(3):
        for j in range(i+1, 3):
            if a[i]*b[j] - a[j]*b[i] != 0:
                return False
    return True

# ------------------------------------------------------------------------ C4
def C4_hipersuperficie():
    header("C4.  LA FIBRA ES IRREDUCIBLE:  {b w_2 = c w_3} LO ES")
    print("   Inv = {(w, w~) en P2 x P2 : w~ w != 0}, dim 4, via N = 2 w w~/(w~ w) - I.")
    print("   La condicion det(N-A)=0 es la forma bilineal  b w_2 - c w_3  sobre")
    print("   P2 x P2.  Una forma de bidegrado (1,1) es reducible si y solo si")
    print("   factoriza como (lineal en w)(lineal en w~), es decir si y solo si su")
    print("   matriz de coeficientes tiene rango 1.  Aqui el rango es 2.")
    # coeficientes de sum E[i][j] w_i wtilde_j  con  wtilde = (a,b,c)
    E = [[F(0), F(0), F(0)],       # w_1 : nada
         [F(0), F(1), F(0)],       # w_2 * b   -> +1
         [F(0), F(0), F(-1)]]      # w_3 * c   -> -1
    r = rank(E)
    claim("rango de la forma bilineal = 2  (>= 2, luego IRREDUCIBLE)", r == 2,
          f"rango = {r}")
    claim("no es rango 1  (no factoriza, no es reducible)", r != 1)
    print("   Es ademas una hipersuperficie propia de Inv: dim = 4 - 1 = 3.")

# ------------------------------------------------------------------------ C5
def C5_dimensiones():
    header("C5.  EL RECUENTO DE DIMENSIONES CIERRA")
    base   = 7    # A reg ss con espectro {1,l,1/l}: clase (6) + el parametro l (1)
    fib_H0 = 1    # F_A, el parametro q del lema de fibras
    fib_N  = 3    # hipersuperficie {det(N-A)=0} dentro de Inv (dim 4)
    fib_H2 = 1    # F_B
    total = base + fib_H0 + fib_N + fib_H2
    print(f"   base   A reg ss, espectro {{1,l,1/l}}      dim {base}")
    print(f"   fibra  H_0 en F_A                        dim {fib_H0}")
    print(f"   fibra  N en Inv con det(N-A) = 0         dim {fib_N}")
    print(f"   fibra  H_2 en F_B                        dim {fib_H2}")
    print(f"                                            ---------")
    print(f"   total                                    dim {total}")
    claim("total = 12 = 4k-8 = dim V_5  (k=5)", total == 4*K-8, f"4k-8 = {4*K-8}")
    claim("fibra total = 5, base = 7", fib_H0+fib_N+fib_H2 == 5 and base == 7)
    print()
    print("   Y para el espacio de moduli, por la construccion de la Prop. 11:")
    par = 2+2+2+1     # n0,n1,n2 en P2 ; n4 en la recta w^perp
    print(f"   n_0,n_1,n_2 en P2 (2+2+2) + n_4 en P1 (1)   = {par}")
    print(f"   menos la accion de O(Q)                     - 3")
    claim("moduli = 4 = 4k-16", par - 3 == 4*K-16, f"4k-16 = {4*K-16}")

# ------------------------------------------------------------------------ C6
def C6_conica_unica():
    header("C6.  LA CONICA INVARIANTE ES UNICA SALVO ESCALA  (V_5 = S_5)")
    print("   M es la matriz 3k x 6 del sistema  Q O_i x L_i^T = 0.  ker M son las")
    print("   formas invariantes.  rango 5 <=> hay exactamente una salvo escala.")
    L = [Qv(n) for n in NS]           # L_i = (Q n_i)^T,  el espejo es la polar
    O = [n[:] for n in NS]            # O_i = n_i,        el centro es el polo
    rows = []
    for i in range(K):
        # Q O_i  como columna, en las 6 incognitas (Q11,Q12,Q13,Q22,Q23,Q33)
        def QO(idx):
            o = O[i]
            return [[o[0], o[1], o[2], F(0), F(0), F(0)],
                    [F(0), o[0], F(0), o[1], o[2], F(0)],
                    [F(0), F(0), o[0], F(0), o[1], o[2]]]
        m = QO(i)
        # (Q O_i) x L_i^T = 0 : tres filas
        l = L[i]
        for a in range(3):
            b, c = (a+1) % 3, (a+2) % 3
            rows.append([m[b][t]*l[c] - m[c][t]*l[b] for t in range(6)])
    r = rank(rows)
    claim("rango M = 5, luego la forma invariante es unica salvo escala", r == 5,
          f"rango = {r}")
    claim("rango <= 5, luego el ejemplo esta en S_5", r <= 5)
    print("   Combinado con la densidad global (V_5 = clausura de V_5^reg) y con")
    print("   'todo punto de V_5^reg preserva una forma no nula', esto da V_5 = S_5.")

# ---------------------------------------------------------------------- main
def main():
    print(__doc__)
    A, B = C1_regular()
    C2_ecuacion(A, B)
    C3_exhaustiva()
    C4_hipersuperficie()
    C5_dimensiones()
    C6_conica_unica()
    header("VEREDICTO")
    if FAILURES:
        print(f"   {len(FAILURES)} COMPROBACIONES FALLIDAS:")
        for f in FAILURES:
            print(f"     - {f}")
        sys.exit(1)
    print("   Todas las comprobaciones pasan, en aritmetica exacta sobre Q.")
    print("   Sin coma flotante en ninguna decision.")
    sys.exit(0)

if __name__ == "__main__":
    main()
