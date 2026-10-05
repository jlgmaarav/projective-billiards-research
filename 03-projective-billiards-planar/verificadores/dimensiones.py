# -*- coding: utf-8 -*-
"""Dimension del espacio de moduli, EN ARITMETICA EXACTA sobre Q.

Convierte en teorema lo que antes era "numerico, 30 soluciones aleatorias".

El argumento:
  - El rango del jacobiano es SEMICONTINUO INFERIORMENTE, luego si en UN punto
    racional vale 8, vale 8 en un abierto de Zariski de esa componente.
  - La dimension del estabilizador es SEMICONTINUA SUPERIORMENTE, luego si en UN
    punto es 0, es 0 genericamente.
  - Por tanto basta UNA cuenta exacta en un punto racional explicito.

    dim V = 6k - (2k escalados) - rango dF = 4k - 8
    dim moduli = (4k - 8) - (8 - dim Stab) = 4k - 16   si Stab es trivial.

Y para el umbral: el estrato con conica tiene dimension 2k+2.
    k = 5:  2k+2 = 12 = 4k-8   -> misma dimension: es UNION DE COMPONENTES.
    k >= 6: 2k+2 < 4k-8        -> cerrado propio: el punto generico no la preserva.

    python dimensiones.py
"""
from fractions import Fraction as F
from itertools import combinations

# ---------------------------------------------------------------- algebra exacta
def mul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def I3():
    return [[F(1) if i == j else F(0) for j in range(3)] for i in range(3)]

def add(A, B):
    return [[A[i][j] + B[i][j] for j in range(3)] for i in range(3)]

def smul(c, A):
    return [[c*A[i][j] for j in range(3)] for i in range(3)]

def outer(o, l):
    """columna x fila -> 3x3."""
    return [[o[i]*l[j] for j in range(3)] for i in range(3)]

def dot(l, o):
    return sum(l[i]*o[i] for i in range(3))

def rank(M):
    """Rango exacto por eliminacion de Gauss sobre Q."""
    M = [row[:] for row in M]
    filas, cols = len(M), len(M[0]) if M else 0
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, filas) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = F(1)/M[r][c]
        M[r] = [x*inv for x in M[r]]
        for i in range(filas):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j] - f*M[r][j] for j in range(cols)]
        r += 1
        if r == filas:
            break
    return r

def nullity(M, ncols):
    return ncols - rank(M)

# ---------------------------------------------------------------- la homologia
def H_of(l, o):
    d = dot(l, o)
    assert d != 0, "transversalidad degenerada"
    return add(I3(), smul(F(-2, 1)/d, outer(o, l)))

def dH_of(l, o, dl, do):
    """Diferencial de H respecto de (L,O), en direccion (dl,do). Exacto."""
    d = dot(l, o)
    t1 = smul(F(-2)/d, add(outer(do, l), outer(o, dl)))
    dd = dot(dl, o) + dot(l, do)
    t2 = smul(F(2)*dd/(d*d), outer(o, l))
    return add(t1, t2)

def producto(Ls, Os):
    P = I3()
    for i in range(len(Ls)):          # P = H_{k-1} ... H_0
        P = mul(H_of(Ls[i], Os[i]), P)
    return P

# ---------------------------------------------------------------- el jacobiano
def jacobiano(Ls, Os):
    """Matriz 8 x 6k de la condicion 'el producto es escalar', en el punto dado."""
    k = len(Ls)
    Hs = [H_of(Ls[i], Os[i]) for i in range(k)]
    # prefijos y sufijos del producto H_{k-1}...H_0
    pre = [I3()]                      # pre[i] = H_{i-1}...H_0
    for i in range(k):
        pre.append(mul(Hs[i], pre[-1]))
    # OJO con el indice: tras invertir, rev[i] = H_{k-1}...H_i, que INCLUYE H_i.
    # Lo que hace falta es H_{k-1}...H_{i+1}, o sea rev[i+1]. Un desplazamiento de
    # uno aqui pasa desapercibido (el rango sigue saliendo 8) y falsea la cuenta.
    rev = [I3()]
    for i in range(k-1, -1, -1):
        rev.append(mul(rev[-1], Hs[i]))
    rev = rev[::-1]                   # rev[i] = H_{k-1}...H_i
    suf = rev[1:] + [I3()]            # suf[i] = H_{k-1}...H_{i+1}
    P = pre[k]
    c = P[0][0]
    assert all(P[i][j] == (c if i == j else 0) for i in range(3) for j in range(3)), \
        "el punto no esta en la variedad"

    cols = []
    for i in range(k):
        for var in range(6):          # 3 de L_i, 3 de O_i
            dl = [F(1) if (var < 3 and t == var) else F(0) for t in range(3)]
            do = [F(1) if (var >= 3 and t == var-3) else F(0) for t in range(3)]
            dH = dH_of(Ls[i], Os[i], dl, do)
            dP = mul(suf[i], mul(dH, pre[i]))
            # parte sin traza de dP/c : 8 coordenadas independientes
            tr = dP[0][0] + dP[1][1] + dP[2][2]
            A = [[dP[a][b] - (tr/3 if a == b else 0) for b in range(3)] for a in range(3)]
            cols.append([A[0][1], A[0][2], A[1][0], A[1][2], A[2][0], A[2][1],
                         A[0][0], A[1][1]])
    # matriz 8 x 6k -> filas = las 8 ecuaciones
    return [[cols[j][e] for j in range(len(cols))] for e in range(8)]

# ---------------------------------------------------------------- estabilizador
def dim_estabilizador(Ls, Os):
    """dim del subgrupo de PGL3 que fija la configuracion ETIQUETADA.

    X en gl3 con  X O_i = a_i O_i  y  L_i X = b_i L_i  para todo i.
    Incognitas: 9 de X, k de a, k de b.  Se resta 1 por X = mult. de I.
    """
    k = len(Ls)
    nvar = 9 + 2*k
    filas = []
    for i in range(k):
        for r in range(3):            # (X O_i - a_i O_i)_r = 0
            fila = [F(0)]*nvar
            for cc in range(3):
                fila[3*r + cc] = Os[i][cc]
            fila[9 + i] = -Os[i][r]
            filas.append(fila)
        for cc in range(3):           # (L_i X - b_i L_i)_c = 0
            fila = [F(0)]*nvar
            for r in range(3):
                fila[3*r + cc] += Ls[i][r]
            fila[9 + k + i] = -Ls[i][cc]
            filas.append(fila)
    return nullity(filas, nvar) - 1   # -1: X proporcional a I actua trivialmente

# ---------------------------------------------------------------- estrato conica
def dim_conicas_invariantes(Ls, Os):
    """dim del espacio de conicas Q con L_i propto (Q O_i)^T para todo i."""
    filas = []
    idx = [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
    for i in range(len(Ls)):
        # Q O_i  x  L_i^T = 0  -> 2 condiciones independientes (producto vectorial)
        base = []
        for (a, b) in idx:
            Q = [[F(0)]*3 for _ in range(3)]
            Q[a][b] = F(1); Q[b][a] = F(1)
            QO = [sum(Q[r][c]*Os[i][c] for c in range(3)) for r in range(3)]
            L = Ls[i]
            cr = [QO[1]*L[2]-QO[2]*L[1], QO[2]*L[0]-QO[0]*L[2], QO[0]*L[1]-QO[1]*L[0]]
            base.append(cr)
        for r in range(3):
            filas.append([base[j][r] for j in range(6)])
    return nullity(filas, 6)

# ---------------------------------------------------------------- configuraciones
K5 = [(1,0,-3), (3,2,2), (-1,0,1), (-8,9,5), (0,-2,2)]
K6 = [(-2,3,3), (-3,-1,-1), (-2,-2,3), (1,-3,2), (2,1,-3), (3,4,3)]

def desde_normales(ns):
    Ls = [[F(x) for x in n] for n in ns]
    Os = [[F(x) for x in n] for n in ns]
    return Ls, Os

def analiza(nombre, ns):
    k = len(ns)
    Ls, Os = desde_normales(ns)
    P = producto(Ls, Os)
    esc = P[0][0]
    print(f"\n{'='*70}\n{nombre}   (k = {k})\n{'='*70}")
    print(f"  producto = {esc} * I   (escalar: {'SI' if all(P[i][j]==(esc if i==j else 0) for i in range(3) for j in range(3)) else 'NO'})")
    J = jacobiano(Ls, Os)
    r = rank(J)
    stab = dim_estabilizador(Ls, Os)
    dimV = 6*k - 2*k - r
    moduli = dimV - (8 - stab)
    conic = dim_conicas_invariantes(Ls, Os)
    print(f"  rango del jacobiano        = {r}      (max posible 8)")
    print(f"  dim del estabilizador      = {stab}")
    print(f"  dim de la variedad         = 6k - 2k - rango = {dimV}      (4k-8 = {4*k-8})")
    print(f"  dim del espacio de moduli  = {moduli}      (4k-16 = {4*k-16})")
    print(f"  dim de conicas invariantes = {conic}")
    print(f"  dim del estrato con conica = 2k+2 = {2*k+2}")
    if 2*k+2 == dimV:
        print("  -> MISMA DIMENSION: el estrato con conica es UNION DE COMPONENTES.")
    elif 2*k+2 < dimV:
        print("  -> ESTRICTAMENTE MENOR: el estrato con conica es un cerrado propio,")
        print("     luego el punto generico de esta componente NO preserva conica.")
    ok = (r == 8 and stab == 0 and dimV == 4*k-8 and moduli == 4*k-16)
    return ok

def main():
    print("DIMENSIONES EN ARITMETICA EXACTA SOBRE Q")
    print("Semicontinuidad: rango 8 y estabilizador 0 en UN punto => genericos en su componente.")
    r1 = analiza("Segundo ejemplo k=5 (datos enteros)", K5)
    r2 = analiza("Contraejemplo k=6 (datos enteros)", K6)
    print()
    print("="*70)
    print("RESULTADO:", "OK - dim moduli = 4k-16 certificado exactamente" if (r1 and r2)
          else "REVISAR")
    print("="*70)
    return r1 and r2


# Guardado bajo __main__ para que `import dimensiones` (lo hace estrato_conica.py
# para reutilizar el algebra exacta) no dispare toda la salida.
if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
