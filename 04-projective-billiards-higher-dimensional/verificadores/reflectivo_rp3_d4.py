# -*- coding: utf-8 -*-
"""UN SEGUNDO BILLAR 6-REFLECTIVO EN R^3, ESTE CALEIDOSCOPICO Y CON DATOS MINIMOS.

    python reflectivo_rp3_d4.py     ->  codigo de salida 0 solo si TODO pasa

Seis espejos del grupo de reflexiones D_4, con normales

    n_0 = ( 1, 0,  1, 0)      n_3 = (0, 1,  1, 0)
    n_1 = ( 1, 0, -1, 0)      n_4 = (0, 1, -1, 0)
    n_2 = ( 1, 0,  0,-1)      n_5 = (0, 1,  0,-1)

y signos e = (-,-,+,-,-,+):  el dominio es
    { x : e_i (a_i . x + d_i) < 0 },   n_i = (a_i, d_i),
un politopo convexo ACOTADO de R^3.  El campo transversal en la cara i es el haz
por el polo O_i = n_i respecto de Q = x0^2+x1^2+x2^2+x3^2.

POR QUE IMPORTA, ademas de ser un segundo ejemplo

Es el analogo en dimension 3 de la familia lambda del plano: un CALEIDOSCOPIO
esferico visto en RP^3, con los seis espejos tomados de un grupo de reflexiones
FINITO de O(4), a saber D_4.  El otro ejemplo de R^3 (reflectivo_rp3.py) tiene
grupo de espejos infinito.  Asi que en dimension 3 se repite la dicotomia del
plano: un ejemplo caleidoscopico y otro que no lo es.

Y los datos son minimos: entradas en {0,1,-1}.

COMO SE ENCONTRO

Barriendo los subconjuntos de 6 espejos de D_4 (12), B_4 (16) y F_4 (24), con
TODOS los 2^6 signos, y para cada uno calculando TODOS los ordenes que dan
producto escalar y comparandolos con los itinerarios que la dinamica produce de
verdad.  La primera version del barrido fijaba e_0 = +1 y por eso no encontro
nada: ese bug se corrigio.

LO QUE SE COMPRUEBA, exacto sobre Q
  C1  las seis son reflexiones de O(Q) y el producto en el orden del itinerario
      es +-I                                                              (H1)
  C2  el politopo es ACOTADO (certificado mu > 0) y NO VACIO
  C3  la orbita testigo racional, itinerario 2,1,3,5,4,0, con desigualdades
      estrictas                                                           (H2)
  C4  el cierre: p_6 = p_0 y u_6 = kappa u_0 con kappa > 0                (H3)
  C5  el grupo de espejos es FINITO: es un caleidoscopio
"""
import sys
from fractions import Fraction as F
from itertools import product as iproduct

M, K = 4, 6
FAILURES = []

def claim(name, cond, detail=""):
    ok = bool(cond)
    print(f"   [{'OK  ' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok: FAILURES.append(name)
    return ok

def header(s):
    print(); print("=" * 78); print(s); print("=" * 78)

def eye(): return [[F(int(i == j)) for j in range(M)] for i in range(M)]
def sc(c): return [[F(c) if i == j else F(0) for j in range(M)] for i in range(M)]
def mul(A, B): return [[sum(A[i][t]*B[t][j] for t in range(M)) for j in range(M)] for i in range(M)]
def eqm(A, B): return all(A[i][j] == B[i][j] for i in range(M) for j in range(M))
def transp(A): return [list(c) for c in zip(*A)]
def mv(A, v): return [sum(A[i][j]*v[j] for j in range(M)) for i in range(M)]
def tr(A): return sum(A[i][i] for i in range(M))
def det4(A):
    from itertools import permutations
    s = F(0)
    for p in permutations(range(M)):
        sg = 1
        for i in range(M):
            for j in range(i+1, M):
                if p[i] > p[j]: sg = -sg
        t = F(sg)
        for i in range(M): t *= A[i][p[i]]
        s += t
    return s

def rank(rows):
    if not rows: return 0
    A = [r[:] for r in rows]; cols = len(A[0]); r = 0
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

N = [(1, 0, 1, 0), (1, 0, -1, 0), (1, 0, 0, -1),
     (0, 1, 1, 0), (0, 1, -1, 0), (0, 1, 0, -1)]
EPS = (-1, -1, 1, -1, -1, 1)
ITIN = [2, 1, 3, 5, 4, 0]        # las caras que se golpean, en orden
NS = [[F(x) for x in n] for n in N]

def refl(nv):
    d = sum(nv[a]*nv[a] for a in range(M))
    return [[F(int(a == b)) - F(2)*nv[a]*nv[b]/d for b in range(M)] for a in range(M)]

R = [refl(n) for n in NS]
A3 = [[F(EPS[i])*NS[i][j] for j in range(3)] for i in range(K)]
D3 = [F(EPS[i])*NS[i][3] for i in range(K)]
def lado(i, x): return sum(A3[i][j]*x[j] for j in range(3)) + D3[i]

def C1():
    header("C1.  SEIS REFLEXIONES DE O(Q), Y EL PRODUCTO EN EL ORDEN DEL ITINERARIO")
    for i in range(K):
        claim(f"R_{i}^2 = I", eqm(mul(R[i], R[i]), eye()))
        claim(f"tr R_{i} = 2", tr(R[i]) == F(2))
        claim(f"R_{i} en O(Q)", eqm(mul(transp(R[i]), R[i]), eye()))
    print(f"   itinerario: se golpean las caras {ITIN} en ese orden, partiendo de la 0.")
    print("   En la convencion de puntos el producto relevante es H_(a_6)...H_(a_1).")
    P = eye()
    for i in ITIN:                     # a_1 primero, aplicado a la derecha
        P = mul(R[i], P)
    esc = eqm(P, sc(-1)) or eqm(P, sc(1))
    claim("el producto en el orden del itinerario es +-I", esc,
          f"= {'-I' if eqm(P, sc(-1)) else ('+I' if eqm(P, sc(1)) else 'no escalar')}")
    claim("las seis normales tienen rango 4", rank([n[:] for n in NS]) == 4)
    caras = set()
    for n in NS:
        g = next(x for x in n if x != 0)
        caras.add(tuple(x/g for x in n))
    claim("las seis caras son distintas", len(caras) == K)

def C2():
    header("C2.  EL POLITOPO ES ACOTADO Y NO VACIO")
    claim("las a_i tienen rango 3", rank([r[:] for r in A3]) == 3)
    rows = [[A3[i][j] for i in range(K)] for j in range(3)]
    A = [r[:] for r in rows]; piv, r = [], 0
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
    mu = None
    for coef in iproduct(range(-6, 7), repeat=len(base)):
        if all(c == 0 for c in coef): continue
        v = [sum(F(coef[t])*base[t][i] for t in range(len(base))) for i in range(K)]
        if all(x > 0 for x in v): mu = v; break
    claim("existe mu > 0 con suma mu_i a_i = 0, luego el cono de recesion es {0}",
          mu is not None, f"mu = {[str(x) for x in mu]}" if mu else "")
    if mu:
        chk = [sum(mu[i]*A3[i][j] for i in range(K)) for j in range(3)]
        claim("y la identidad se cumple exactamente", all(x == 0 for x in chk))
    x0 = [F(1, 2), F(1, 2), F(-1, 4)]
    print(f"   punto interior x0 = ({x0[0]}, {x0[1]}, {x0[2]})")
    for i in range(K):
        v = lado(i, x0)
        claim(f"lado {i} < 0", v < 0, f"= {v}")

def C3_C4():
    header("C3.  LA ORBITA TESTIGO RACIONAL")
    print("   La cara 0 es  -(x1 + x3) = 0, luego x3 = -x1.")
    p = [F(4, 7), F(4, 5), F(-4, 7)]
    u = [F(15), F(5), F(3)]
    print(f"   p_0 = ({p[0]}, {p[1]}, {p[2]})")
    print(f"   u_0 = ({u[0]}, {u[1]}, {u[2]})")
    claim("p_0 esta en la cara 0", lado(0, p) == 0)
    for i in range(1, K):
        claim(f"p_0 estrictamente interior respecto de la cara {i}", lado(i, p) < 0,
              f"= {lado(i,p)}")
    claim("u_0 apunta hacia dentro", sum(A3[0][j]*u[j] for j in range(3)) < 0)
    cara, pa, ua, u0 = 0, p[:], u[:], u[:]
    for step in range(K):
        cands = []
        for i in range(K):
            if i == cara: continue
            den = sum(A3[i][j]*ua[j] for j in range(3))
            claim(f"paso {step+1}: la cara {i} no es paralela al rayo", den != 0) \
                if False else None
            if den == 0: continue
            s = -lado(i, pa)/den
            if s > 0:
                cands.append((s, i, [pa[j] + s*ua[j] for j in range(3)]))
        claim(f"paso {step+1}: hay impacto", len(cands) > 0)
        if not cands: return
        cands.sort(key=lambda t: t[0])
        s, i, q = cands[0]
        claim(f"paso {step+1}: el primer impacto es la cara {ITIN[step]}",
              i == ITIN[step], f"cara {i}, s = {s}")
        claim(f"paso {step+1}: s > 0 estricto", s > 0)
        claim(f"paso {step+1}: las perdedoras excluidas ESTRICTAMENTE",
              all(s2 > s for (s2, j2, q2) in cands[1:]))
        claim(f"paso {step+1}: el impacto es INTERIOR a la cara {i}",
              all(lado(j, q) < 0 for j in range(K) if j != i))
        claim(f"paso {step+1}: transversalidad no nula",
              sum(NS[i][t]*NS[i][t] for t in range(M)) != 0)
        wv = mv(R[i], [ua[0], ua[1], ua[2], F(0)])
        ua = [wv[j] - wv[3]*q[j] for j in range(3)]
        claim(f"paso {step+1}: direccion reflejada no nula", any(x != 0 for x in ua))
        pa, cara = q, i
    header("C4.  EL CIERRE")
    claim("p_6 = p_0 EXACTAMENTE", pa == p, f"p_6 = {[str(x) for x in pa]}")
    kap = None
    for j in range(3):
        if u0[j] != 0: kap = ua[j]/u0[j]; break
    claim("u_6 proporcional a u_0",
          all(ua[j]*u0[t] - ua[t]*u0[j] == 0 for j in range(3) for t in range(3)))
    claim("kappa > 0", kap is not None and kap > 0, f"kappa = {kap}")
    claim("seis caras distintas: periodo primitivo 6", len(set(ITIN)) == K)

def C5():
    header("C5.  ES UN CALEIDOSCOPIO: EL GRUPO DE ESPEJOS ES FINITO")
    print("   Los seis espejos son de D_4, que es un grupo de reflexiones FINITO de")
    print("   orden 192.  Se comprueba cerrando el grupo generado por busqueda en")
    print("   anchura, en PGL_4 y en aritmetica exacta.")
    def key(A):
        flat = [A[i][j] for i in range(M) for j in range(M)]
        for x in flat:
            if x != 0:
                if x < 0: flat = [-y for y in flat]
                break
        return tuple(flat)
    seen = {key(eye())}
    frontier = [eye()]
    CAP = 3000
    while frontier and len(seen) <= CAP:
        nxt = []
        for A in frontier:
            for g in R:
                B = mul(g, A); k = key(B)
                if k not in seen:
                    seen.add(k); nxt.append(B)
        frontier = nxt
    claim("el grupo es FINITO en PGL_4", len(seen) <= CAP, f"|G| = {len(seen)}")
    print(f"   |G| = {len(seen)} en PGL_4, es decir {2*len(seen)} en O(4).")
    print("   Frente al otro ejemplo de R^3, cuyo grupo es INFINITO (15 de 15 pares")
    print("   con orden infinito por Niven).  Se repite la dicotomia del plano.")

def main():
    print(__doc__)
    C1(); C2(); C3_C4(); C5()
    header("VEREDICTO")
    if FAILURES:
        print(f"   {len(FAILURES)} COMPROBACIONES FALLIDAS:")
        for f in FAILURES: print(f"     - {f}")
        sys.exit(1)
    print("   Todas las comprobaciones pasan, exacto sobre Q.")
    print()
    print("   SEGUNDO BILLAR 6-REFLECTIVO EN R^3, CALEIDOSCOPICO, CON DATOS EN")
    print("   {0, 1, -1}.  Es el analogo en dimension 3 de la familia lambda.")
    sys.exit(0)

if __name__ == "__main__":
    main()
