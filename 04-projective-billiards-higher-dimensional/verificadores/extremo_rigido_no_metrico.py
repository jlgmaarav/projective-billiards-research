# -*- coding: utf-8 -*-
"""EN EL EXTREMO RIGIDO k=n+2 EL BILLAR NO PUEDE SER METRICO.  EXACTO SOBRE Q.

    python extremo_rigido_no_metrico.py   ->  codigo 0 solo si TODO pasa

LA APARENTE CONTRADICCION QUE ESTO RESUELVE

Dos resultados del proyecto parecian chocar el 9-VIII-2026:

  `thm:rigidend`   Si Q es DEFINIDA POSITIVA y R_0...R_{k-1} son Q-reflexiones con
                   k = n+2 y producto escalar, entonces los normales no generan, los
                   espejos son concurrentes y la region es no acotada: NO HAY BILLAR.
  Tabla de SOL     Para n par existe un billar k-reflectivo para todo k >= n+1, y en
                   particular para k = n+2, con producto +I y grupo infinito;
                   certificado en `verificadores/rank5_k6.py` con (n,k) = (4,6).

No se contradicen, y la razon es la hipotesis: la construccion de SOL es de CENTRO
COMUN.  Pero «tiene centro comun» no es por si solo «no preserva ninguna cuadrica»,
que es lo que hace falta.  Eso es lo que se demuestra aqui.

EL LEMA, Y ES DE TRES LINEAS

  LEMA.  Sea una configuracion de homologias armonicas con CENTRO COMUN O y al menos
  DOS ejes distintos.  Entonces no preserva ninguna cuadrica NO DEGENERADA.

  Demostracion.  Si S es invariante y no degenerada, cada H_i es una involucion que
  preserva S, luego sus autoespacios +1 y -1 son S-ortogonales.  El autoespacio -1 es
  <O> y el +1 es ker L_i, asi que <SO, x> = 0 para todo x con L_i x = 0, es decir
  SO pertenece a <L_i>.  Como S es no degenerada y O no es cero, SO no es cero, luego
  L_i es proporcional a SO PARA TODO i: todos los ejes coinciden.  []

Notese que sale sin usar definitud: se excluye toda cuadrica no degenerada, no solo
las definidas.  Es mas fuerte de lo que hacia falta.

LA LECTURA ESTRUCTURAL, QUE ES LO QUE APORTA

  El extremo rigido k = n+2 NO esta vacio: esta vacio de billares METRICOS.  Existe
  un k-reflectivo con k = n+2, pero esta OBLIGADO a vivir en la rama 3 de la
  tricotomia de `rem:branches` -- sin conica invariante --, igual que los ejemplos
  euclideos de Fierobe.  `thm:rigidend` no dice que el extremo rigido este vacio;
  dice que alli se acaba la geometria metrica.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  la configuracion (4,6) tiene centro comun, ejes distintos y producto +I_5
  C2  el lema, sobre esos datos: se calcula el espacio de cuadricas invariantes y se
      comprueba que TODAS sus elementos son degenerados, exhibiendo el vector del
      nucleo comun
  C3  contraste con `thm:rp4b5`, que si es metrico: alli el espacio invariante
      contiene la identidad, que es definida
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


def nucleo(rows, nc):
    A = [r[:] for r in rows]
    nr = len(A)
    piv, r = [], 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = F(1) / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(nr):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    libres = [c for c in range(nc) if c not in piv]
    base = []
    for fc in libres:
        v = [F(0)] * nc
        v[fc] = F(1)
        for i, c in enumerate(piv):
            v[c] = -A[i][fc]
        base.append(v)
    return base


def mul(A, B, m):
    return [[sum((A[i][k] * B[k][j] for k in range(m)), F(0)) for j in range(m)]
            for i in range(m)]


def transp(A, m):
    return [[A[j][i] for j in range(m)] for i in range(m)]


def cuadricas_invariantes(Hs, m):
    """Base exacta de {S simetrica : H^T S H = S para toda H de Hs}."""
    idx = [(i, j) for i in range(m) for j in range(i, m)]

    def desde(v):
        S = [[F(0)] * m for _ in range(m)]
        for t, (i, j) in enumerate(idx):
            S[i][j] = v[t]
            S[j][i] = v[t]
        return S

    ecu = []
    for H in Hs:
        Ht = transp(H, m)
        cols = []
        for t in range(len(idx)):
            e = [F(0)] * len(idx)
            e[t] = F(1)
            S = desde(e)
            D = mul(Ht, mul(S, H, m), m)
            cols.append([D[a][b] - S[a][b] for a, b in idx])
        for a in range(len(idx)):
            ecu.append([cols[t][a] for t in range(len(idx))])
    return [desde(v) for v in nucleo(ecu, len(idx))]


def homology(L, O, m):
    c = sum((L[a] * O[a] for a in range(m)), F(0))
    assert c != 0
    return [[F(int(a == b)) - F(2) * O[a] * L[b] / c for b in range(m)]
            for a in range(m)]


# --------------------------------------------- datos de rank5_k6.py, (n,k)=(4,6)
N46 = [[F(1), F(0), F(0), F(0)], [F(0), F(1), F(0), F(0)],
       [F(0), F(0), F(1), F(0)], [F(0), F(0), F(0), F(1)],
       [F(-1), F(0), F(-1), F(0)], [F(0), F(-1), F(0), F(-1)]]
L46 = [[r[0], r[1], F(-1), r[2], r[3]] for r in N46]
O46 = [F(0), F(0), F(1), F(0), F(0)]
H46 = [homology(L, O46, 5) for L in L46]
WORD46 = [2, 1, 0, 5, 4, 3]

ROOTS_B5 = [(0, 1, 0, 0, -1), (0, 0, 1, 0, -1), (0, 0, 0, 1, 0), (1, -1, 0, 0, 0),
            (0, 1, 1, 0, 0), (1, 0, 0, 0, 1), (1, 0, -1, 0, 0)]


def refl(n, m):
    d = sum((F(x * x) for x in n), F(0))
    return [[F(int(a == b)) - F(2) * F(n[a]) * F(n[b]) / d for b in range(m)]
            for a in range(m)]


HB5 = [refl(n, 5) for n in ROOTS_B5]


def main():
    header("C1.  LA CONFIGURACION (n,k) = (4,6): CENTRO COMUN, EJES DISTINTOS, +I")
    print(f"   centro comun O = {[str(x) for x in O46]}")
    claim("los seis ejes son distintos dos a dos",
          len({tuple(L) for L in L46}) == 6)
    claim("todos son homologias armonicas involutivas",
          all(mul(H, H, 5) == [[F(int(i == j)) for j in range(5)] for i in range(5)]
              for H in H46))
    P = [[F(int(i == j)) for j in range(5)] for i in range(5)]
    for i in WORD46:
        P = mul(H46[i], P, 5)
    claim("el producto en el orden de la palabra es +I_5",
          P == [[F(int(i == j)) for j in range(5)] for i in range(5)])
    claim("y k = n+2: es el extremo rigido", 6 == 4 + 2, "n = 4, k = 6")

    header("C2.  EL LEMA SOBRE ESTOS DATOS: TODA CUADRICA INVARIANTE ES DEGENERADA")
    base = cuadricas_invariantes(H46, 5)
    print(f"   dimension del espacio de cuadricas invariantes: {len(base)}")
    claim("el espacio es no trivial", len(base) > 0)
    # el lema predice que todas son degeneradas; se comprueba exhibiendo un vector
    # que esta en el nucleo de TODAS a la vez
    filas = []
    for S in base:
        filas.extend([row[:] for row in S])
    ker = nucleo(filas, 5)
    claim("existe un vector no nulo en el nucleo de TODAS las cuadricas invariantes",
          len(ker) > 0, f"dim del nucleo comun = {len(ker)}")
    if ker:
        print(f"      un tal vector: {[str(x) for x in ker[0]]}")
    for i, S in enumerate(base):
        r = rank([row[:] for row in S])
        if r >= 5:
            claim(f"generador {i} degenerado", False, f"rango {r}")
    claim("luego NINGUNA combinacion es no degenerada: todas anulan ese vector",
          len(ker) > 0)
    print("   Es exactamente lo que dice el lema: con centro comun y dos ejes")
    print("   distintos, una cuadrica invariante no degenerada obligaria a que todos")
    print("   los ejes fueran proporcionales a S*O, o sea a que coincidieran.")

    header("C3.  CONTRASTE: thm:rp4b5 SI ES METRICO")
    baseB = cuadricas_invariantes(HB5, 5)
    print(f"   dimension del espacio invariante: {len(baseB)}")
    I5 = [[F(int(i == j)) for j in range(5)] for i in range(5)]
    filasB = []
    for S in baseB:
        filasB.extend([row[:] for row in S])
    kerB = nucleo(filasB, 5)
    claim("alli NO hay vector comun en el nucleo de todas", len(kerB) == 0,
          f"dim del nucleo comun = {len(kerB)}")
    claim("y la identidad es invariante: cuadrica DEFINIDA",
          all(mul(transp(H, 5), mul(I5, H, 5), 5) == I5 for H in HB5))
    print("   Luego el mismo calculo distingue las dos situaciones sin ambiguedad.")

    print()
    print("=" * 78)
    if FAILURES:
        print("FALLOS:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("TODO CORRECTO.")
    print()
    print("  El extremo rigido k = n+2 no esta vacio: esta vacio de billares")
    print("  METRICOS.  `thm:rigidend` y la tabla de SOL son las dos ciertas, y")
    print("  juntas dicen mas que cualquiera por separado: en k = n+2 la ausencia de")
    print("  conica invariante no es una eleccion, es OBLIGATORIA.  Es la rama 3 de")
    print("  rem:branches, la misma de los ejemplos euclideos de Fierobe.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
