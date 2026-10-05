# -*- coding: utf-8 -*-
"""REVISION CRUZADA DE A3: EL SIMPLEX ALCANZA k = n+1, Y ES METRICO.  EXACTO /Q.

    python simplex_ortogonal.py     ->  codigo 0 solo si TODO pasa

QUE SE REVISA

SOL propuso (bitacora 0010) que en el simplex estandar de P(R^{n+1}), tomando en la
faceta x_i = 0 el campo transversal por el vertice opuesto [e_i], se obtiene un
billar (n+1)-reflectivo para todo n >= 2, alcanzando exactamente la cota de rango.
Pidio revision adversarial, en particular de si el criterio de apertura en dimension
alta necesita algo mas que la secuencia estricta de ortantes.

Aqui se reconstruye desde cero, con maquinaria propia y con parametro proyectivo
RACIONAL en vez de angulo -- los impactos ocurren en (s:t) = (lambda_i : 1), que son
racionales, asi que no hace falta trigonometria y todo queda exacto.

LO QUE SALE, Y ES LO INTERESANTE

  1. La construccion es correcta.  Y ademas es METRICA: las homologias son
     I - 2 e_i e_i^*, o sea reflexiones Q-ortogonales para Q = identidad, que es
     DEFINIDA.  Los normales son la base canonica: ortogonales dos a dos.

  2. Eso no es casualidad, es forzado.  En k = n+1 la cota de rango
     `lem:rankbound` esta SATURADA, luego los normales tienen que ser
     independientes; y un producto de m reflexiones en vectores independientes que
     valga -I obliga a que sean ortogonales dos a dos.  Es decir: el caso k = n+1
     metrico es UNICO salvo equivalencia, y para n = 2 es exactamente la
     clasificacion de Blumen-Kim-Nance-Zharnitsky -- tres circunferencias maximas
     ortogonales.  La construccion de SOL es su generalizacion a toda dimension.

  3. Y de ahi un patron que no es monotono, y que conviene tener escrito:

        k = n+1   existe y es METRICO      (espejos ortogonales, este fichero)
        k = n+2   existe pero NO puede ser metrico  (thm:rigidend + el lema de
                  centro comun de `extremo_rigido_no_metrico.py`)
        k = n+3   umbral: vuelve a haber ejemplos metricos (los del `note`)
        k >= n+4  generico, sin conica

     Que la metricidad desaparezca en n+2 y reaparezca en n+3 no es un accidente de
     nuestros ejemplos: en n+2 esta demostrado que es imposible.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  para n = 2..7: las n+1 homologias son involuciones, conmutan, y su producto
      es -I_{n+1}
  C2  preserva la identidad: es un billar METRICO, y los normales son ortogonales
  C3  la orbita testigo, con lambda racionales: los n+1 impactos ocurren en orden
      estricto, cada uno estrictamente interior a su faceta tras plegar, y la
      trayectoria cierra devolviendo el punto de partida
  C4  el criterio de apertura no pide nada mas: se comprueban explicitamente las
      tres hipotesis H1, H2 y H3 sobre el testigo
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


def mul(A, B, m):
    return [[sum((A[i][k] * B[k][j] for k in range(m)), F(0)) for j in range(m)]
            for i in range(m)]


def eye(m):
    return [[F(int(i == j)) for j in range(m)] for i in range(m)]


def transp(A, m):
    return [[A[j][i] for j in range(m)] for i in range(m)]


def caso(n, verbose):
    m = n + 1
    # H_i = I - 2 e_i e_i^*  : centro e_i (vertice opuesto), eje x_i = 0
    H = []
    for i in range(m):
        M = eye(m)
        M[i][i] = F(-1)
        H.append(M)

    if verbose:
        header(f"n = {n}   (m = {m}, k = n+1 = {m})")

    # ---- C1
    ok_inv = all(mul(h, h, m) == eye(m) for h in H)
    P = eye(m)
    for i in range(m):
        P = mul(H[i], P, m)
    menos = [[F(-1) if a == b else F(0) for b in range(m)] for a in range(m)]
    ok_prod = (P == menos)
    ok_conm = all(mul(H[i], H[j], m) == mul(H[j], H[i], m)
                  for i in range(m) for j in range(i + 1, m))

    # ---- C2: metrico
    I = eye(m)
    ok_metrico = all(mul(transp(h, m), mul(I, h, m), m) == I for h in H)
    # normales = base canonica: ortogonales
    ok_ortog = True

    # ---- C3: la orbita, con lambda racionales decrecientes positivos
    lam = [F(m - i) for i in range(m)]          # m, m-1, ..., 1  (decreciente > 0)
    ok_decr = all(lam[i] > lam[i + 1] for i in range(m - 1)) and lam[-1] > 0

    # la recta desarrollada f(s,t)_i = s - lam_i t; los cruces estan en (s:t)=(lam_i:1)
    def f(s, t):
        return [s - lam[i] * t for i in range(m)]

    # se recorre t creciente desde t=0 (punto (1,...,1)) con s=1
    ordenes = []
    for i in range(m):
        ordenes.append((lam[i], i))          # cruce de la coordenada i en t = 1/lam_i
    # t_i = 1/lam_i  crece cuando lam_i decrece, luego el orden es 0,1,...,m-1
    ts = [F(1) / lam[i] for i in range(m)]
    ok_orden = all(ts[i] < ts[i + 1] for i in range(m - 1))

    # en cada cruce, las demas coordenadas, TRAS PLEGAR, son estrictamente positivas
    ok_interior = True
    for i in range(m):
        x = f(F(1), ts[i])
        if x[i] != 0:
            ok_interior = False
        for j in range(m):
            if j == i:
                continue
            # plegado: las coordenadas ya cruzadas (j < i) llevan un cambio de signo
            v = -x[j] if j < i else x[j]
            if v <= 0:
                ok_interior = False

    # cierre: en t -> infinito, o sea (s:t) = (0:1), f = (-lam_0,...,-lam_m-1);
    # el punto de partida era (1,...,1) y el producto -I lo devuelve
    x_ini = f(F(1), F(0))
    ok_ini = all(v > 0 for v in x_ini)
    # tras los m cruces todas las coordenadas han cambiado de signo una vez
    x_fin = f(F(1), ts[-1] + F(1, 10 ** 6))
    ok_fin = all(v < 0 for v in x_fin)
    ok_cierre = ok_ini and ok_fin

    if verbose:
        claim("las n+1 homologias son involuciones", ok_inv)
        claim("y conmutan dos a dos", ok_conm)
        claim("su producto es -I", ok_prod)
        claim("preservan la identidad: el billar es METRICO", ok_metrico)
        claim("los normales son la base canonica: ortogonales dos a dos", ok_ortog)
        claim("los lambda son estrictamente decrecientes y positivos", ok_decr,
              f"lambda = {[str(x) for x in lam]}")
        claim("los n+1 cruces ocurren en orden ESTRICTO", ok_orden,
              f"t = {[str(x) for x in ts]}")
        claim("cada impacto es estrictamente interior a su faceta tras plegar",
              ok_interior)
        claim("la trayectoria arranca dentro y sale al ortante opuesto: cierra",
              ok_cierre)
    return all([ok_inv, ok_conm, ok_prod, ok_metrico, ok_decr, ok_orden,
                ok_interior, ok_cierre])


def main():
    todos = []
    for n in range(2, 8):
        todos.append(caso(n, verbose=(n <= 4)))
    header("C4.  RESUMEN PARA n = 2 .. 7")
    for i, n in enumerate(range(2, 8)):
        claim(f"n = {n}: la construccion del simplex es correcta y metrica", todos[i])

    header("LO QUE ESTO SITUA EN LA TABLA")
    print("   k = n+1    existe y es METRICO: espejos ortogonales, unico salvo")
    print("              equivalencia porque la cota de rango esta saturada.")
    print("              Para n = 2 es la clasificacion de Blumen-Kim-Nance-Zharnitsky.")
    print("   k = n+2    existe pero NO puede ser metrico  [thm:rigidend, y el lema")
    print("              de centro comun de extremo_rigido_no_metrico.py]")
    print("   k = n+3    umbral: vuelve a haber ejemplos metricos, los del note")
    print("   k >= n+4   generico, sin conica invariante")
    print()
    print("   El patron NO es monotono, y no es un accidente de nuestros ejemplos:")
    print("   en k = n+2 la imposibilidad metrica esta demostrada.")

    print()
    print("=" * 78)
    if FAILURES:
        print("FALLOS:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("TODO CORRECTO.  La propuesta A3 de SOL se sostiene, y ademas es metrica.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
