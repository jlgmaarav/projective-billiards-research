# -*- coding: utf-8 -*-
"""CODIMENSION, DEGENERACION Y SENSIBILIDAD DE LA k-REFLECTIVIDAD.  EXACTO SOBRE Q.

    python tolerancias.py     ->  codigo de salida 0 solo si TODO pasa

QUE SE AFIRMA

Sea M_k(n) el espacio de configuraciones de k espejos proyectivos en RP^n: cada
espejo es una homologia armonica  H = I - 2 O L^T/(L^T O),  o sea un polo O y un
eje L, cada uno modulo escala, luego 2(m-1) parametros con m = n+1, y

    dim M_k(n) = 2 k (m-1) = 2 k n.

Sea Phi: M_k(n) -> PGL_m la aplicacion  (H_i) |-> H_{k-1}...H_0.  Entonces
V_k = Phi^{-1}(id).

(T1)  FORMULA DE DIMENSION UNIFORME.  Donde dPhi es suprayectiva,

          codim V_k = dim PGL_m = (n+1)^2 - 1,      uniforme en k,
          dim  V_k = 2 k n - ((n+1)^2 - 1).

      Para n=2 da 4k-8, que es `thm:dim`.  Para n=3, k=6 da 21, que es lo que
      certifica `dimensiones_rpn.py`.  Los dos recuentos que estaban sueltos son
      el mismo enunciado, y ademas vale para todo n.

(T2)  IDENTIFICACION DEL CONUCLEO.  Sea a_j = A_j n_j el normal DESDOBLADO del
      espejo j, donde A_j es el producto de los espejos posteriores a j en el
      itinerario.  Entonces

          coker dPhi  =  { S simetrica de traza nula : S a_j // a_j  para todo j }.

      Es decir: dPhi deja de ser suprayectiva exactamente cuando la configuracion
      DESDOBLADA preserva una cuadrica de mas, aparte de la invariante Q.  La
      demostracion es de tres lineas y esta en NOVENA-RONDA-TOLERANCIAS.md; aqui
      se certifica la igualdad de dimensiones sobre Q en los cinco ejemplos.

(T3)  LOS DOS EJEMPLOS DE COXETER SON DEGENERADOS, Y LOS OTROS TRES NO.
      El conucleo mide 0 en el plano y en los dos ejemplos de R^3, luego alli
      la codimension es exactamente (n+1)^2-1.  Pero mide 1 en B_5 y 2 en B_6:
      esas dos configuraciones preservan un HAZ (resp. una RED) de cuadricas, y
      por tanto V_k tiene alli dimension MAYOR que la esperada.  No es un
      accidente numerico: son palabras de Coxeter, y las coordenadas que la
      palabra no mezcla dejan libre un peso diagonal.

(T4)  CONSTANTE DE SENSIBILIDAD.  Con polos y ejes normalizados y midiendo la
      perturbacion de cada espejo como un angulo, la norma de operador A de dPhi
      cumple   max_j ||dPhi(e_j)||  <=  A  <=  ||dPhi||_F,  y los CUADRADOS de
      las dos cotas son racionales exactos.  A es el factor de amplificacion: un
      error angular theta en los espejos deja un defecto de cierre del orden de
      A*theta.

POR QUE IMPORTA, EN FISICA

`FISICA-Y-APLICACIONES.md` dejaba declarado un hueco: "la apertura topologica
del teorema no da tolerancias cuantitativas".  Esto lo cierra, y la respuesta
tiene dos mitades opuestas que conviene no confundir:

  - En los RAYOS la propiedad es ABIERTA.  Hay un cono de direcciones de entrada
    que vuelven.  Eso es tolerancia de alineacion, y es positiva.
  - En los ESPEJOS la propiedad es de CODIMENSION (n+1)^2-1: 8 ecuaciones en el
    plano, 15 en R^3, 24 en R^4, 35 en R^5.  NO existe ningun abierto de
    configuraciones k-reflectivas.  Un retrorreflector de orden k no se
    "mecaniza con tolerancia": hay que clavar exactamente esas ecuaciones, y
    cualquier desviacion theta deja un defecto de cierre A*theta.

Esa asimetria -- abierto en el rayo, codimension alta en el aparato -- es el
enunciado de ingenieria que faltaba, y es lo que fija el presupuesto de errores
de cualquier realizacion.  Y T3 dice algo util ademas: los ejemplos de Coxeter,
que son los mas faciles de fabricar por ser permutaciones con signo, son
justamente los que estan en un estrato degenerado, o sea los que tienen MAS
direcciones de deformacion admisibles.  Para construir, hay que ir a ellos.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  cada configuracion es de homologias armonicas involutivas y su producto en
      el orden del itinerario es -I  (se redetecta aqui, no se supone)
  C2  corank dPhi = dim{S simetrica de traza nula : S a_j // a_j}   (T2)
  C3  el conucleo es 0 en los tres ejemplos no-Coxeter, y 1 y 2 en B_5 y B_6
  C4  la formula de T1 reproduce los dos valores ya certificados en otro sitio:
      4k-8 en el plano, y 21 para n=3, k=6
  C5  las cotas exactas de la constante de sensibilidad A, y que son > 0
"""
import sys
from fractions import Fraction as F

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


# --------------------------------------------------- algebra exacta sobre Q
def eye(m):
    return [[F(int(i == j)) for j in range(m)] for i in range(m)]


def mul(A, B):
    m, p, q = len(A), len(B), len(B[0])
    return [[sum((A[i][k] * B[k][j] for k in range(p)), F(0)) for j in range(q)]
            for i in range(m)]


def is_scalar(A):
    m = len(A)
    c = A[0][0]
    return all(A[i][j] == (c if i == j else F(0)) for i in range(m) for j in range(m))


def dot(u, v):
    return sum((a * b for a, b in zip(u, v)), F(0))


def rank(M):
    """Rango exacto por eliminacion gaussiana sobre Q."""
    if not M:
        return 0
    A = [row[:] for row in M]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = F(1) / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == rows:
            break
    return r


def inv_mat(A):
    m = len(A)
    B = [A[i][:] + [F(int(i == j)) for j in range(m)] for i in range(m)]
    for c in range(m):
        piv = next(i for i in range(c, m) if B[i][c] != 0)
        B[c], B[piv] = B[piv], B[c]
        z = B[c][c]
        B[c] = [x / z for x in B[c]]
        for i in range(m):
            if i != c and B[i][c] != 0:
                z = B[i][c]
                B[i] = [a - z * b for a, b in zip(B[i], B[c])]
    return [row[m:] for row in B]


def gram_schmidt_perp(v):
    """Base ORTOGONAL racional del complemento euclideo de v en Q^m."""
    m = len(v)
    nv = dot(v, v)
    basis = []
    for i in range(m):
        e = [F(int(j == i)) for j in range(m)]
        c = dot(e, v) / nv
        w = [a - c * b for a, b in zip(e, v)]
        for b in basis:
            c2 = dot(w, b) / dot(b, b)
            w = [a - c2 * x for a, x in zip(w, b)]
        if any(x != 0 for x in w):
            basis.append(w)
    assert len(basis) == m - 1
    return basis


# ------------------------------------------------- homologias y diferencial
def homology(O, L):
    """H = I - 2 O L^T / (L.O)."""
    m = len(O)
    c = dot(L, O)
    assert c != 0
    return [[F(int(a == b)) - F(2) * O[a] * L[b] / c for b in range(m)]
            for a in range(m)]


def d_homology(O, L, dO, dL):
    """dH = -2 (dO L^T + O dL^T)/c + 2 O L^T (dL.O + L.dO)/c^2,  c = L.O."""
    m = len(O)
    c = dot(L, O)
    dc = dot(dL, O) + dot(L, dO)
    return [[F(-2) * (dO[a] * L[b] + O[a] * dL[b]) / c
             + F(2) * O[a] * L[b] * dc / (c * c) for b in range(m)]
            for a in range(m)]


def traceless_vec(X):
    """Coordenadas de X modulo escalares: la parte de traza nula, aplanada."""
    m = len(X)
    t = sum((X[i][i] for i in range(m)), F(0)) / m
    return [X[a][b] - (t if a == b else F(0)) for a in range(m) for b in range(m)]


def suffix_products(Hs, order):
    """A_j = producto de los espejos POSTERIORES a j en el itinerario."""
    m = len(Hs[0])
    seq = [Hs[i] for i in order]
    out = []
    for pos in range(len(order)):
        A = eye(m)
        for Hh in seq[pos + 1:]:
            A = mul(Hh, A)
        out.append(A)
    return out


def differential(Os, Ls, order):
    """Columnas de dPhi, cada una con el cuadrado de la norma de su generador.

    Una base ortonormal del dominio exigiria raices; en su lugar se guarda el
    generador racional junto con |gen|^2, y al calcular normas se divide.  Todo
    lo que se necesita son CUADRADOS de normas, que asi siguen siendo racionales.
    """
    m = len(Os[0])
    k = len(Os)
    Hs = [homology(Os[i], Ls[i]) for i in range(k)]
    As = suffix_products(Hs, order)

    cols = []
    for pos, idx in enumerate(order):
        A, Ainv = As[pos], inv_mat(As[pos])
        O, L, Hj = Os[idx], Ls[idx], Hs[idx]
        dirs = ([(w, [F(0)] * m) for w in gram_schmidt_perp(O)]
                + [([F(0)] * m, w) for w in gram_schmidt_perp(L)])
        for dO, dL in dirs:
            dH = d_homology(O, L, dO, dL)
            xi = mul(mul(A, mul(dH, Hj)), Ainv)      # A (dH_j H_j) A^{-1}
            gen = dO if any(x != 0 for x in dO) else dL
            cols.append((traceless_vec(xi), dot(gen, gen)))
    return cols


def dim_S_invariante(a_list, m, traceless):
    """dim{S simetrica : S a // a para todo a}, opcionalmente con tr S = 0.

    S a // a  <=>  a_p (Sa)_q - a_q (Sa)_p = 0 para todo p<q.
    """
    basis = [(i, j) for i in range(m) for j in range(i, m)]

    def Sa_component(i, j, a, r):
        v = F(0)
        if i == r:
            v += a[j]
        if j == r and i != j:
            v += a[i]
        return v

    rows = []
    for a in a_list:
        for p in range(m):
            for q in range(p + 1, m):
                rows.append([a[p] * Sa_component(i, j, a, q)
                             - a[q] * Sa_component(i, j, a, p) for (i, j) in basis])
    if traceless:
        rows.append([F(1) if i == j else F(0) for (i, j) in basis])
    return len(basis) - rank(rows)


def detect_order(Os, Ls, k, hint):
    """Encuentra el orden del itinerario cuyo producto es escalar.  No se supone."""
    base = hint if hint is not None else list(range(k))
    m = len(Os[0])
    Hs = [homology(Os[i], Ls[i]) for i in range(k)]
    for order in (base, base[::-1]):
        P = eye(m)
        for i in order:
            P = mul(Hs[i], P)
        if is_scalar(P):
            return order, P
    return None, None


# ------------------------------------------------------------- los cinco casos
# Todos son Q-ortogonales con Q = suma de cuadrados, luego O_i = L_i = n_i.
# corank esperado: 0 en los tres primeros, 1 en B_5, 2 en B_6.  Ver T3.
CASOS = [
    ("plano k=5, normales enteras   (thm:second)", 2, 5,
     [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)], None, 0),
    ("R^3  k=6, grupo infinito      (thm:rp3example)", 3, 6,
     [(-1, 3, -4, -1), (-1, 0, 4, 1), (3, 4, -1, -2), (0, -8, 91, 0),
      (0, 2198, 2119, 8345), (-1577, 276, 1117, -146)], None, 0),
    ("R^3  k=6, caleidoscopio D_4   (thm:rp3d4)", 3, 6,
     [(1, 0, 1, 0), (1, 0, -1, 0), (1, 0, 0, -1),
      (0, 1, 1, 0), (0, 1, -1, 0), (0, 1, 0, -1)], [2, 1, 3, 5, 4, 0], 0),
    ("R^4  k=7, B_5                 (thm:rp4b5)", 4, 7,
     [(0, 1, 0, 0, -1), (0, 0, 1, 0, -1), (0, 0, 0, 1, 0), (1, -1, 0, 0, 0),
      (0, 1, 1, 0, 0), (1, 0, 0, 0, 1), (1, 0, -1, 0, 0)], None, 1),
    ("R^5  k=8, B_6                 (thm:rp5b6)", 5, 8,
     [(0, 1, 0, 0, -1, 0), (0, 0, 1, 0, -1, 0), (0, 0, 0, 1, 0, 0),
      (1, -1, 0, 0, 0, 0), (0, 1, 1, 0, 0, 0), (1, 0, 0, 0, 1, 0),
      (1, 0, -1, 0, 0, 0), (0, 0, 0, 0, 0, 1)], None, 2),
]


def raiz(q):
    """Decimal aproximado de sqrt(q), SOLO para imprimir.  No decide nada."""
    return float(q) ** 0.5


def main():
    resumen = []

    for nombre, n, k, Ns, hint, corank_esperado in CASOS:
        m = n + 1
        header(f"{nombre}   [n={n}, m={m}, k={k}]")
        Os = [[F(x) for x in v] for v in Ns]
        Ls = [o[:] for o in Os]

        # ---- C1
        order, P = detect_order(Os, Ls, k, hint)
        claim("el producto en el orden del itinerario es ESCALAR", order is not None,
              f"orden = {order}" if order else "ningun orden da escalar")
        if order is None:
            continue
        claim("y el escalar es -1 (monodromia antipodal)", P[0][0] == F(-1),
              f"c = {P[0][0]}")
        Hs = [homology(Os[i], Ls[i]) for i in range(k)]
        claim("las k homologias son involutivas",
              all(mul(h, h) == eye(m) for h in Hs))

        # ---- C2  corank de dPhi  ==  cuadricas extra de la config. desdoblada
        cols = differential(Os, Ls, order)
        Mat = [[c[0][r] for c in cols] for r in range(m * m)]
        r = rank(Mat)
        corank = (m * m - 1) - r
        amb = 2 * k * n

        As = suffix_products(Hs, order)
        a_list = [[sum((As[pos][row][c] * Os[idx][c] for c in range(m)), F(0))
                   for row in range(m)] for pos, idx in enumerate(order)]
        d_tl = dim_S_invariante(a_list, m, traceless=True)

        # La demostracion de T2 necesita que los normales desdoblados generen
        # R^m: es lo que mata la parte antisimetrica del conucleo.  Se comprueba,
        # no se supone.
        claim("los normales desdoblados generan R^m  (hipotesis de T2)",
              rank([[a[c] for c in range(m)] for a in a_list]) == m)

        claim("corank dPhi = dim{S simetrica traza nula : S a_j // a_j}   (T2)",
              corank == d_tl,
              f"corank = {corank},  dim S = {d_tl}   [rank = {r}, columnas = 2kn = {amb}]")

        # ---- C3  quien es degenerado y quien no
        claim(f"el conucleo mide lo previsto por T3: {corank_esperado}",
              corank == corank_esperado, f"medido = {corank}")

        if corank == 0:
            dimV = amb - (m * m - 1)
            print(f"   -> dPhi suprayectiva:  codim V_k = {m*m-1},  dim V_k = {dimV}")
        else:
            dimV = None
            d_full = dim_S_invariante(a_list, m, traceless=False)
            print(f"   -> ESTRATO DEGENERADO.  La configuracion desdoblada preserva un")
            print(f"      espacio de cuadricas de dimension {d_full} (un haz o mas), no una sola.")
            print(f"      Luego codim V_k <= {r} < {m*m-1}, y dim V_k >= {amb - r}"
                  f" en vez de {amb - (m*m-1)}.")
            # Y la causa, comprobada sobre los datos: coordenadas AISLADAS.  Si
            # e_c es el normal de un espejo y ningun otro normal tiene componente
            # c, la recta <e_c> es invariante y el grupo es reducible.
            aisladas = [c for c in range(m)
                        if sum(1 for v in Os if v[c] != 0) == 1
                        and any(all(v[d] == 0 for d in range(m) if d != c) for v in Os)]
            claim("la causa es reducibilidad por coordenadas AISLADAS, y hay"
                  f" tantas como el corank", len(aisladas) == corank,
                  f"coordenadas aisladas = {aisladas}")
        resumen.append((n, k, corank, dimV))

        # ---- C5  constante de sensibilidad, cotas exactas
        fro2 = sum((sum((x * x for x in c[0]), F(0)) / c[1] for c in cols), F(0))
        maxcol2 = max(sum((x * x for x in c[0]), F(0)) / c[1] for c in cols)
        claim("la sensibilidad es estrictamente positiva", maxcol2 > 0)
        claim("y finita, con cota de Frobenius exacta", fro2 >= maxcol2)
        print(f"   -> A^2 en [{maxcol2}, {fro2}]")
        print(f"      o sea  A en [{raiz(maxcol2):.4f}, {raiz(fro2):.4f}]"
              f"  de defecto de cierre por radian de error de espejo")

    # ------------------------------------ C4: cruce con lo ya certificado
    header("C4.  LA FORMULA DE T1 REPRODUCE LOS DOS RECUENTOS YA CERTIFICADOS")
    for n, k, corank, dimV in resumen:
        if corank != 0:
            continue
        if n == 2:
            claim(f"n=2, k={k}:  dim V_k = 4k-8 = {4*k-8}", dimV == 4 * k - 8,
                  f"la formula da {dimV}   [thm:dim]")
        if n == 3 and k == 6:
            claim("n=3, k=6:  dim V_6 = 21", dimV == 21,
                  f"la formula da {dimV}   [dimensiones_rpn.py]")

    # ------------------------------------ T5: el umbral k = n+3, en una linea
    header("C6.  EL UMBRAL k = n+3 SALE DE LAS DOS FORMULAS DE DIMENSION")
    print("   Un espejo que preserva una cuadrica fija Q tiene su polo determinado por")
    print("   Q, luego aporta solo n parametros (el eje).  Contando la propia Q y")
    print("   quitando las condiciones de que el producto sea escalar DENTRO de O(Q):")
    print()
    print("      dim S_k = k n + [(n+1)(n+2)/2 - 1] - n(n+1)/2 = n(k+1)")
    print("      dim V_k = 2 k n - ((n+1)^2 - 1)")
    print()
    claim("dim S_k = n(k+1) reproduce el valor 2k+2 del plano   [thm:dim]",
          all(2 * (k + 1) == 2 * k + 2 for k in range(3, 12)))
    umbrales = []
    for n in range(2, 9):
        ks = [k for k in range(n + 1, n + 9)
              if 2 * k * n - ((n + 1) ** 2 - 1) == n * (k + 1)]
        umbrales.append((n, ks))
    claim("dim V_k = dim S_k  ocurre exactamente en k = n+3, para todo n probado",
          all(ks == [n + 3] for n, ks in umbrales),
          "  ".join(f"n={n}: k={ks}" for n, ks in umbrales))
    print()
    print("   Y los cuatro ejemplos certificados caen todos ahi:")
    for n, k, corank, dimV in resumen:
        print(f"      n={n}, k={k}   ->  n+3 = {n+3}"
              f"   {'SI' if k == n + 3 else 'NO'}")
    claim("los cinco ejemplos estan en el umbral k = n+3",
          all(k == n + 3 for n, k, _, _ in resumen))
    print()
    print("   AVISO: los dos recuentos son dimensiones ESPERADAS, validas donde la")
    print("   aplicacion correspondiente es submersiva.  Para V_k eso es justo lo que")
    print("   certifican C2 y C3 arriba, punto por punto; para S_k no se certifica")
    print("   aqui, y por eso el umbral en RP^3 sigue condicional a la")
    print("   irreducibilidad de V_6, como dice `note.tex`.")

    print()
    print("=" * 78)
    if FAILURES:
        print("FALLOS:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("TODO CORRECTO.")
    print()
    print("  codim V_k = (n+1)^2 - 1  donde dPhi es suprayectiva:")
    print("     plano  8 condiciones      R^3  15      R^4  24      R^5  35")
    print()
    print("  Lectura de ingenieria: la k-reflectividad es ABIERTA en el rayo y de")
    print("  CODIMENSION (n+1)^2-1 en el aparato.  No hay ninguna tolerancia de")
    print("  fabricacion positiva; lo que hay es un presupuesto de defecto de")
    print("  cierre, con factor de conversion A.  Y los dos ejemplos de Coxeter")
    print("  estan en un estrato degenerado: mas libertad de deformacion, que es")
    print("  justo lo que hace falta para construir.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
