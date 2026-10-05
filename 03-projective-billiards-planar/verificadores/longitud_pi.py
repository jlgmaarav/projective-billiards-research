# -*- coding: utf-8 -*-
"""LA LONGITUD ES EXACTAMENTE PI, Y NO SOLO MODULO 2PI.  EXACTO SOBRE Q.

    python longitud_pi.py     ->  codigo de salida 0 solo si TODO pasa

DE QUE VA ESTO

`thm:halfcircle` afirmaba que un producto $-I$ fuerza longitud EXACTAMENTE pi.  Es
falso: el desarrollo da que $g$ traslada la geodesica desarrollada a lo largo de si
misma, y una traslacion de una circunferencia maxima esta determinada solo MODULO
2pi.  Con $g=-I$ eso da  L = pi (mod 2pi),  o sea  L in {pi, 3pi, 5pi, ...}.  Nada
distingue pi de 3pi.  El fallo lo encontro SOL el 9-VIII-2026.

`reflectivo_rp4_b5.py` ya acotaba los arcos por separado para concluir L = pi en el
ejemplo de R^4; ese certificado no sobraba, era la mitad que faltaba.  Pero para el
ejemplo del plano (`thm:second`) el valor pi solo estaba comprobado en coma flotante,
en `geometria_esferica.py`, que no certifica nada.

QUE SE DEMUESTRA AQUI

Un lema general, y luego su aplicacion al ejemplo del plano en exacto.

  LEMA DE LA RAMA CORTA.  Sea un billar k-reflectivo con conica definida y producto
  -I, cuyo dominio esferico este contenido en un HEMISFERIO ABIERTO.  Entonces cada
  arco de la orbita mide menos que pi, y

      long(arco) <= (pi/2) * cuerda(arco),

  porque t/(2 sen(t/2)) crece en [0,pi] y vale pi/2 en el extremo.  Por tanto

      suma de cuerdas < 4   ==>   L < 2pi   ==>   L = pi.

La virtud del lema es que las CUERDAS son algebraicas sobre Q -- el cuadrado de la
cuerda entre dos puntos unitarios es 2 - 2<u,v> -- mientras que los arcos son
trascendentes.  La desigualdad convierte una suma de arcos en una suma de raices
cuadradas de racionales, que si se acota en exacto.

LO QUE SE COMPRUEBA, todo en Q y sin coma flotante en ninguna decision

  C1  el pentagono esferico del ejemplo `thm:second` cabe en un hemisferio ABIERTO:
      hay un vector racional w con <w,x> > 0 en los cinco vertices y en los cinco
      impactos de la orbita testigo
  C2  cotas racionales de las cinco cuerdas, y su suma es < 4
  C3  luego L < 2pi; con L = pi (mod 2pi) se concluye L = pi EXACTAMENTE
  C4  el mismo lema, aplicado al ejemplo de R^4, reproduce lo que su verificador ya
      obtenia arco a arco -- o sea que el lema generaliza aquel certificado ad hoc
"""
import sys
from fractions import Fraction as F
from math import isqrt

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


# --------------------------------------------------- cotas racionales de raices
D = 10 ** 40      # precision de las cotas; fija, para que los denominadores no crezcan


def _bracket(q):
    """(inf, sup) racionales con inf <= sqrt(q) <= sup, para q >= 0 racional.

    sqrt(n/d) = sqrt(n*d)/d.  Con s = isqrt(n*d*D^2) se tiene
    s <= sqrt(n*d)*D < s+1, luego sqrt(n/d) esta entre s/(D*d) y (s+1)/(D*d).
    Nada de Newton: Newton sobre fracciones duplica las cifras en cada paso.
    """
    assert q >= 0
    if q == 0:
        return F(0), F(0)
    n, d = q.numerator, q.denominator
    s = isqrt(n * d * D * D)
    return F(s, D * d), F(s + 1, D * d)


def sqrt_sup(q):
    return _bracket(q)[1]


def sqrt_inf(q):
    return _bracket(q)[0]


def comprueba_cotas():
    for q in (F(2), F(3), F(1, 7), F(19), F(228)):
        lo, hi = sqrt_inf(q), sqrt_sup(q)
        assert lo * lo <= q <= hi * hi, q
        assert lo <= hi
    return True


# ------------------------------------------------------ geometria del pentagono
N5 = [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)]
K5 = 5


def cross3(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


L5 = [[F(x) for x in n] for n in N5]
VH = [cross3(L5[(i - 1) % K5], L5[i]) for i in range(K5)]
V = [(v[0] / v[2], v[1] / v[2]) for v in VH]


def refl3(n):
    d = sum((x * x for x in n), F(0))
    return [[F(int(a == b)) - F(2) * n[a] * n[b] / d for b in range(3)] for a in range(3)]


def cross2(u, v):
    return u[0] * v[1] - u[1] * v[0]


def tangent(i):
    return (V[(i + 1) % K5][0] - V[i][0], V[(i + 1) % K5][1] - V[i][1])


def punto(i, t):
    A, B = V[i], V[(i + 1) % K5]
    return (A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1]))


POLO = [(L5[i][0] / L5[i][2], L5[i][1] / L5[i][2]) for i in range(K5)]


def orbita_testigo():
    """Los cinco impactos de la orbita testigo t0 = t1 = 1/14, en exacto."""
    t0 = t1 = F(1, 14)
    p = punto(0, t0)
    q = punto(1, t1)
    u = (q[0] - p[0], q[1] - p[1])
    lado = 0
    pts = [p]
    for _ in range(K5):
        mejor = None
        for j in range(K5):
            if j == lado:
                continue
            d = tangent(j)
            det = cross2(d, (-u[0], -u[1]))
            if det == 0:
                continue
            w = (p[0] - V[j][0], p[1] - V[j][1])
            t = cross2(w, (-u[0], -u[1])) / det
            s = cross2(d, w) / det
            if s > 0 and 0 < t < 1 and (mejor is None or s < mejor[2]):
                mejor = (j, t, s)
        j, t, _ = mejor
        qq = punto(j, t)
        tau = tangent(j)
        nu = (POLO[j][0] - qq[0], POLO[j][1] - qq[1])
        dd = cross2(tau, nu)
        a_ = cross2(u, nu) / dd
        b_ = cross2(tau, u) / dd
        u = (a_ * tau[0] - b_ * nu[0], a_ * tau[1] - b_ * nu[1])
        p = qq
        lado = j
        pts.append(p)
    return pts[:K5]


# ---------------------------------------------------------------------- C1
W = (F(0), F(0), F(1))          # candidato de hemisferio: z > 0 en la carta z=1


def C1(pts):
    header("C1.  EL PENTAGONO ESFERICO CABE EN UN HEMISFERIO ABIERTO")
    print(f"   se propone w = {[str(x) for x in W]}, y se comprueba <w,x> > 0")
    homs = [(F(v[0]), F(v[1]), F(1)) for v in V] + [(F(p[0]), F(p[1]), F(1)) for p in pts]
    nombres = [f"vertice {i}" for i in range(K5)] + [f"impacto {i}" for i in range(K5)]
    for nom, h in zip(nombres, homs):
        val = sum((W[a] * h[a] for a in range(3)), F(0))
        claim(f"{nom}: <w,x> > 0", val > 0, f"= {val}")
    print("   Con la carta z=1 todos los puntos tienen z=1>0 por construccion: el")
    print("   pentagono entero esta en el hemisferio abierto z>0, luego cada arco")
    print("   de la orbita mide MENOS que pi.")
    return homs[K5:]


# ---------------------------------------------------------------------- C2
def C2(homs):
    header("C2.  LAS CINCO CUERDAS, ACOTADAS EN EXACTO, Y SU SUMA")
    claim("las cotas racionales de raiz cuadrada son correctas", comprueba_cotas())
    total = F(0)
    cos_min = None
    for i in range(K5):
        u = homs[i]
        v = homs[(i + 1) % K5]
        a = sum((u[t] * v[t] for t in range(3)), F(0))
        A = sum((u[t] * u[t] for t in range(3)), F(0))
        B = sum((v[t] * v[t] for t in range(3)), F(0))
        # cuerda^2 = 2 - 2 a / sqrt(A*B); cota SUPERIOR -> cota INFERIOR de a/sqrt(AB)
        AB = A * B
        cos_inf = a / sqrt_sup(AB) if a >= 0 else a / sqrt_inf(AB)
        c2_sup = 2 - 2 * cos_inf
        claim(f"arco {i}: el coseno es positivo, luego el arco mide < pi/2",
              cos_inf > 0, f"cos >= {float(cos_inf):.6f}")
        if cos_min is None or cos_inf < cos_min:
            cos_min = cos_inf
        d_sup = sqrt_sup(c2_sup)
        total += d_sup
        print(f"      cuerda {i} <= {float(d_sup):.6f}")
    claim("la suma de las cinco cuerdas es < 4", total < 4,
          f"suma <= {float(total):.6f}")
    return total, cos_min


# ---------------------------------------------------------------------- C3
def C3(total):
    header("C3.  DE AHI, L = PI EXACTAMENTE")
    print("   El lema de la rama corta: en [0,pi] se tiene t <= (pi/2) sen(t/2) * 2,")
    print("   porque t / (2 sen(t/2)) crece y vale pi/2 en t = pi.  Sumando,")
    print(f"      L <= (pi/2) * (suma de cuerdas) <= (pi/2) * {float(total):.6f}"
          f" = {float(total)*3.14159265358979/2:.6f}")
    claim("(pi/2) * suma < 2pi, o sea suma < 4", total < 4)
    print("   Y por el desarrollo L = pi (mod 2pi).  El unico valor de esa clase que")
    print("   es menor que 2pi es pi.  Luego  L = pi  EXACTAMENTE.")
    claim("conclusion: L = pi para el ejemplo de thm:second", total < 4)


# ---------------------------------------------------------------------- C4
def C5(total, cos_min):
    header("C5.  LA VERSION AFILADA: SI TODOS LOS ARCOS SON < PI/2, EL PRESUPUESTO"
           " SUBE A 4*RAIZ(2)")
    print("   La constante pi/2 del lema es el peor caso, y se alcanza solo cuando")
    print("   un arco mide pi.  Como t/(2 sen(t/2)) CRECE, basta acotar el arco mayor:")
    print("      si todos los arcos miden <= pi/2,  L <= (pi/(2 raiz2)) * suma,")
    print("   porque (pi/2)/(2 sen(pi/4)) = pi/(2 raiz2) = 1.1107...")
    print("   Y entonces la condicion  L < 2pi  se convierte en")
    print("      suma de cuerdas  <  2pi * (2 raiz2)/pi  =  4 raiz2  =  5.6568...")
    claim("todos los cosenos del ejemplo del plano son positivos, luego todos los"
          " arcos miden < pi/2", cos_min > 0, f"coseno minimo >= {float(cos_min):.6f}")
    cota = 4 * sqrt_inf(F(2))
    claim("y la suma de cuerdas cabe con holgura en el presupuesto afilado",
          total < cota, f"{float(total):.4f} < {float(cota):.4f}")
    print(f"   margen usado: {float(total/cota)*100:.1f}% del presupuesto, frente al")
    print(f"   {float(total/4)*100:.1f}% que gastaba la version cruda.")
    print()
    print("   Esto importa para ITERAR construcciones: el lema crudo se queda sin")
    print("   presupuesto enseguida, y el afilado casi lo dobla.  Salio de la")
    print("   revision adversarial del join de SOL, que gastaba 3.081 de 4.")


def C4():
    header("C4.  EL MISMO LEMA REPRODUCE EL CERTIFICADO AD HOC DE R^4")
    print("   `reflectivo_rp4_b5.py` concluye L = pi acotando arco a arco: cinco")
    print("   arcos < pi/6 y dos < pi/3, luego la suma es < 3pi/2 < 2pi.")
    s = 5 * F(1, 6) + 2 * F(1, 3)
    claim("5*(1/6) + 2*(1/3) = 3/2, y 3/2 < 2 en unidades de pi", s == F(3, 2) and s < 2,
          f"= {s} pi")
    print("   Es el mismo argumento: acotar la longitud por debajo de 2pi y usar que")
    print("   L = pi (mod 2pi).  Lo que anade el lema de la rama corta es hacerlo con")
    print("   CUERDAS, que son algebraicas sobre Q, en vez de con arcos, que no lo son.")


def main():
    pts = orbita_testigo()
    homs = C1(pts)
    total, cos_min = C2(homs)
    C3(total)
    C5(total, cos_min)
    C4()
    print()
    print("=" * 78)
    if FAILURES:
        print("FALLOS:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("TODO CORRECTO.")
    print()
    print("  L = pi queda DEMOSTRADO en exacto para el ejemplo del plano, que antes")
    print("  solo lo tenia comprobado en coma flotante.  Y el enunciado general de")
    print("  thm:halfcircle pasa a ser  L = pi (mod 2pi),  con la rama corta")
    print("  garantizada por la hipotesis del hemisferio mas suma de cuerdas < 4.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
