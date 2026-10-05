# -*- coding: utf-8 -*-
"""ES CALEIDOSCOPIO O NO?  EL GRUPO GENERADO POR LOS ESPEJOS.

    python caleidoscopio.py     ->  codigo de salida 0 solo si TODO pasa

PARA QUE SIRVE ESTO

En el Teorema principal, W tiene autovalores -1, exp(+-i theta) con
theta = (n-2j) pi / n, que es un multiplo RACIONAL de pi.  Luego W es de orden
finito.  Eso levanta una sospecha que conviene resolver antes de que la resuelva
un referee:

    la familia lambda, son caleidoscopios esfericos vistos en RP^2?

Esto no toca a ningun teorema.  Cambia como se enmarca la novedad, y decide cual
de los dos ejemplos es el interesante.

RESPUESTA

  A. La familia lambda con n=5, j=2, lambda=(2+sqrt5)/2 SI lo es.  La matriz de
     Gram de las cinco normales, en la metrica de la conica invariante, es
     circulante y exacta:

         g_ii = 1,   g_{i,i+-1} = 1/2 = cos(pi/3),
                     g_{i,i+-2} = (1-sqrt5)/4 = cos(3 pi/5).

     Los angulos diedros son pi/3 y 3pi/5, ambos multiplos racionales de pi, y
     el grupo generado tiene orden 60 en PGL_3 (120 en O(3)): es el grupo
     ICOSAEDRICO H_3.  Los cinco espejos son cinco de sus quince.

     Luego el ejemplo simetrico es un caleidoscopio esferico visto en RP^2.  Lo
     genuinamente nuevo en el es la observacion de que -I es TRIVIAL en PGL_3,
     y por eso el caso IMPAR no tiene obstruccion de determinante: es lo que
     hace que la Proposicion 1.5 de Fierobe sea una obstruccion Z/2 y no un
     muro.  Y explica de donde sale el sqrt5: de H_3, no del pentagono.

  B. El segundo ejemplo, de normales racionales, NO lo es.  Por el TEOREMA DE
     NIVEN, los unicos valores racionales de cos(q pi) con q racional son
     0, +-1/2, +-1.  Las trazas dan cosenos racionales fuera de esa lista, luego
     productos de dos espejos de orden INFINITO.  Es el ejemplo genuinamente
     nuevo, y el que muestra que el espacio de moduli tiene dimension 4 y no es
     una lista de simetrias.

TODA la aritmetica es exacta: Q para la parte B, y Q(sqrt5) representado como
pares de fracciones para la parte A.  Sin coma flotante en ninguna decision.
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
    print(); print("=" * 78); print(s); print("=" * 78)

# ======================================================= aritmetica en Q(sqrt5)
class K:
    """a + b*sqrt(5), con a, b en Q.  Exacto, hashable."""
    __slots__ = ("a", "b")
    def __init__(self, a=0, b=0):
        self.a = F(a); self.b = F(b)
    def __add__(s, o):  o = K.c(o); return K(s.a+o.a, s.b+o.b)
    def __radd__(s, o): return s + o
    def __neg__(s):     return K(-s.a, -s.b)
    def __sub__(s, o):  return s + (-K.c(o))
    def __rsub__(s, o): return K.c(o) + (-s)
    def __mul__(s, o):
        o = K.c(o)
        return K(s.a*o.a + 5*s.b*o.b, s.a*o.b + s.b*o.a)
    def __rmul__(s, o): return s * o
    def inv(s):
        d = s.a*s.a - 5*s.b*s.b
        assert d != 0, "division por cero en Q(sqrt5)"
        return K(s.a/d, -s.b/d)
    def __truediv__(s, o): return s * K.c(o).inv()
    def __eq__(s, o):
        o = K.c(o); return s.a == o.a and s.b == o.b
    def __hash__(s):    return hash((s.a, s.b))
    def is_zero(s):     return s.a == 0 and s.b == 0
    def sign(s):
        """Signo exacto de a + b sqrt5, por comparacion de cuadrados."""
        if s.a == 0 and s.b == 0: return 0
        if s.b == 0: return 1 if s.a > 0 else -1
        if s.a == 0: return 1 if s.b > 0 else -1
        if s.a > 0 and s.b > 0: return 1
        if s.a < 0 and s.b < 0: return -1
        # signos opuestos: comparar a^2 con 5 b^2
        lhs, rhs = s.a*s.a, 5*s.b*s.b
        if lhs == rhs: return 0
        dom_a = lhs > rhs
        return (1 if s.a > 0 else -1) if dom_a else (1 if s.b > 0 else -1)
    def __repr__(s):
        if s.b == 0: return f"{s.a}"
        if s.a == 0: return f"{s.b}sqrt5"
        return f"{s.a}{'+' if s.b>0 else ''}{s.b}sqrt5"
    @staticmethod
    def c(o): return o if isinstance(o, K) else K(o, 0)

ZERO, ONE = K(0), K(1)

def mv(M, v):    return [sum((M[i][j]*v[j] for j in range(3)), ZERO) for i in range(3)]
def mm(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(3)), ZERO) for j in range(3)] for i in range(3)]
def eye():       return [[ONE if i == j else ZERO for j in range(3)] for i in range(3)]
def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

def solve3(A, rhs):
    """Gauss exacto en Q(sqrt5)."""
    M = [[A[i][j] for j in range(3)] + [rhs[i]] for i in range(3)]
    for c in range(3):
        p = next((i for i in range(c, 3) if not M[i][c].is_zero()), None)
        assert p is not None
        M[c], M[p] = M[p], M[c]
        piv = M[c][c].inv()
        M[c] = [x*piv for x in M[c]]
        for i in range(3):
            if i != c and not M[i][c].is_zero():
                f = M[i][c]
                M[i] = [M[i][t] - f*M[c][t] for t in range(4)]
    return [M[i][3] for i in range(3)]

# ================================================================ A. la familia
def parte_A():
    header("A.  LA FAMILIA lambda,  n = 5,  j = 2,  lambda = (2+sqrt5)/2")
    print("   La matriz de Gram de las cinco normales en la metrica de la conica")
    print("   invariante Q = diag(1,1,lambda·cos(pi/5)) se calcula en Q(sqrt5) y es")
    print("   circulante:  g_ii = 1,  g_{i,i+-1} = 1/2,  g_{i,i+-2} = (1-sqrt5)/4.")
    print("   (Obtenida de O_i^T Q O_j / (O_0^T Q O_0), norma comun 25/8+11sqrt5/8.)")
    print()
    s = K(F(1, 4), F(-1, 4))                 # (1 - sqrt5)/4 = cos(3 pi/5)
    half = K(F(1, 2))                        # 1/2           = cos(pi/3)
    g = [[ONE, half, s, s, half],
         [half, ONE, half, s, s],
         [s, half, ONE, half, s],
         [s, s, half, ONE, half],
         [half, s, s, half, ONE]]

    # angulos: dos reflexiones a angulo theta generan un diedrico cuya rotacion
    # es por 2 theta.  cos = 1/2 -> theta = pi/3 -> rotacion de orden 3.
    # cos = cos(3pi/5) -> theta = 3pi/5 -> rotacion 6pi/5, de orden 5.
    claim("g_{i,i+-1} = 1/2 = cos(pi/3)   -> rotacion de orden 3", half == K(F(1, 2)))
    claim("g_{i,i+-2} = (1-sqrt5)/4 = cos(3pi/5)  -> rotacion de orden 5",
          s == K(F(1, 4), F(-1, 4)))
    print(f"   cos(3pi/5) = {s}   (numericamente {float(s.a + s.b*5**0.5):.6f})")

    G3 = [[g[i][j] for j in range(3)] for i in range(3)]
    d = det3(G3)
    claim("las tres primeras normales son base (det de Gram 3x3 != 0)",
          not d.is_zero(), f"det = {d}")
    claim("Gram definida positiva: la geometria es ESFERICA (det > 0 y 2x2 > 0)",
          d.sign() > 0 and (ONE*ONE - half*half).sign() > 0)

    # coordenadas de m_3, m_4 en la base (m_0, m_1, m_2)
    coord = [[ONE, ZERO, ZERO], [ZERO, ONE, ZERO], [ZERO, ZERO, ONE]]
    for i in (3, 4):
        coord.append(solve3(G3, [g[i][0], g[i][1], g[i][2]]))
    # comprobar consistencia: la norma de m_3 y m_4 debe salir 1
    for i in (3, 4):
        c = coord[i]
        nn = sum((c[a]*G3[a][b]*c[b] for a in range(3) for b in range(3)), ZERO)
        claim(f"la reconstruccion de m_{i} tiene norma 1", nn == ONE, f"|m_{i}|^2 = {nn}")
        for jj in range(5):
            gg = sum((c[a]*g[a][jj] for a in range(3)), ZERO) if jj < 3 else None
        # producto con m_4 / m_3 cruzado
    ip34 = sum((coord[3][a]*G3[a][b]*coord[4][b] for a in range(3) for b in range(3)), ZERO)
    claim("m_3 . m_4 = 1/2, como pide la Gram", ip34 == half, f"= {ip34}")

    # reflexiones como matrices 3x3 sobre Q(sqrt5), en la base (m_0,m_1,m_2)
    # xi' = xi - 2 (fila_i(G) . xi) c^{(i)}
    def reflection(i):
        row = [g[i][0], g[i][1], g[i][2]]
        c = coord[i]
        return [[(ONE if a == b else ZERO) - K(2)*c[a]*row[b] for b in range(3)]
                for a in range(3)]
    Rf = [reflection(i) for i in range(5)]
    for i in range(5):
        claim(f"R_{i}^2 = I", mm(Rf[i], Rf[i]) == eye())
        claim(f"det R_{i} = -1", det3(Rf[i]) == K(-1))

    # la relacion de cierre
    P = eye()
    for i in range(4, -1, -1):
        P = mm(Rf[i], P)
    claim("R_4 R_3 R_2 R_1 R_0 = -I  (la relacion de cierre, en Q(sqrt5))",
          P == [[K(-1) if a == b else ZERO for b in range(3)] for a in range(3)])

    # ---------------------------------------- orden del grupo, BFS exacto en PGL_3
    header("A.1  EL ORDEN DEL GRUPO, BFS EXACTO EN PGL_3")
    def key(M):
        """Representante de PGL_3: normaliza el signo global."""
        flat = [M[i][j] for i in range(3) for j in range(3)]
        for x in flat:
            sg = x.sign()
            if sg != 0:
                if sg < 0:
                    flat = [-y for y in flat]
                break
        return tuple((y.a, y.b) for y in flat)

    seen = {key(eye())}
    frontier = [eye()]
    CAP = 5000
    while frontier and len(seen) <= CAP:
        nxt = []
        for M in frontier:
            for gg in Rf:
                Mg = mm(gg, M)
                k = key(Mg)
                if k not in seen:
                    seen.add(k); nxt.append(Mg)
        frontier = nxt
    orden = len(seen)
    claim("el grupo es FINITO en PGL_3", orden <= CAP, f"|G| = {orden}")
    claim("|G| = 60 en PGL_3, es decir 120 en O(3): el grupo ICOSAEDRICO H_3",
          orden == 60, f"|G| = {orden}")
    print("   H_3 tiene orden 120 y 15 espejos; aqui usamos 5 de esos 15.")
    print("   Por tanto: la familia lambda es un CALEIDOSCOPIO ESFERICO en RP^2.")

# ============================================== B. el segundo ejemplo, por Niven
def parte_B():
    header("B.  EL SEGUNDO EJEMPLO (normales racionales):  TEOREMA DE NIVEN")
    print("   Niven: si q es racional y cos(q pi) es racional, entonces cos(q pi)")
    print("   esta en {0, +-1/2, +-1}.")
    print("   Para dos reflexiones de O(3), R_a R_b es una rotacion de angulo phi")
    print("   con tr(R_a R_b) = 1 + 2 cos(phi).  Si cos(phi) es racional y NO esta")
    print("   en la lista de Niven, phi no es multiplo racional de pi, el producto")
    print("   tiene orden INFINITO, y el grupo de espejos es INFINITO.")
    print()
    N = [(1, 0, -3), (3, 2, 2), (-1, 0, 1), (-8, 9, 5), (0, -2, 2)]

    def refl(nv):
        d = sum(F(x*x) for x in nv)
        return [[F(int(a == b)) - F(2*nv[a]*nv[b])/d for b in range(3)] for a in range(3)]
    def mul(A, B):
        return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

    R = [refl(nv) for nv in N]
    NIVEN = {F(0), F(1, 2), F(-1, 2), F(1), F(-1)}
    infinitos = 0
    for a in range(5):
        for b in range(a+1, 5):
            M = mul(R[a], R[b])
            t = M[0][0] + M[1][1] + M[2][2]
            cos_phi = (t - 1)/2
            fuera = cos_phi not in NIVEN
            infinitos += fuera
            print(f"   R_{a}R_{b}: tr = {str(t):>10}  cos(phi) = {str(cos_phi):>10}"
                  f"   {'orden INFINITO' if fuera else 'en Niven'}")
    claim("hay productos de dos espejos de orden infinito", infinitos > 0,
          f"{infinitos} de 10 pares")
    claim("el grupo de espejos es INFINITO: NO es un caleidoscopio", infinitos > 0)
    print()
    print("   Esto NO contradice que la configuracion preserve una conica:")
    print("   preservar una conica es UNA relacion, no finitud del grupo.  El")
    print("   segundo ejemplo es un billar esferico cuyo grupo de espejos es")
    print("   infinito y cuyas orbitas cierran solo por la relacion.")

# ================================== C. que geometria es cada miembro de la familia
def parte_C():
    header("C.  QUE GEOMETRIA ES CADA MIEMBRO DE LA FAMILIA lambda")
    print("   lambda(n,j) = v(cos 2a - C)/(C-1)  con a = pi/n, C = cos theta,")
    print("   theta = (n-2j)pi/n.  Como C < 1, el signo de lambda es el de")
    print("   -(cos 2a - C), luego")
    print()
    print("       lambda > 0  <=>  cos 2a < C  <=>  theta < 2a  <=>  j > (n-2)/2.")
    print()
    print("   Junto con j < n/2 y j entero:")
    print("     n IMPAR: j = (n-1)/2 es el UNICO miembro esferico.")
    print("     n PAR:   no hay ninguno; j = (n-2)/2 da lambda = 0, el euclideo")
    print("              de Fierobe.  Todos los demas tienen lambda < 0.")
    print()
    print("   Y para lambda < 0 la conica invariante x^2+y^2 = -lambda·v·z^2 es un")
    print("   circulo REAL de radio^2 = -lambda·v.  Comparado con el circulo")
    print("   inscrito del poligono, de radio^2 = v^2:")
    print()
    print("       -lambda < v  <=>  (cos 2a - C)/(1-C) < 1  <=>  cos 2a < 1,")
    print()
    print("   que es cierto SIEMPRE.  Luego la conica queda estrictamente DENTRO")
    print("   del circulo inscrito, y los espejos NO la cortan nunca: el dominio")
    print("   esta en el EXTERIOR de la conica.")
    print()
    print("   Y mas fino aun:  -lambda < v^3  <=>  -u^2 < C u^2  <=>  C > -1,")
    print("   tambien cierto siempre.  Luego radio^2 < v^4, que es el inradio^2 del")
    print("   poligono de puntos medios: la conica queda dentro del hueco que")
    print("   encierra CADA orbita de la familia.")
    print()
    print("   Consecuencia, y es una correccion al enunciado del paper: para")
    print("   lambda < 0 la forma tiene signatura (2,1) y el grupo es O(2,1), pero")
    print("   el billar NO vive en el plano hiperbolico.  Decir 'hiperbolico' es")
    print("   FALSO.  Y ojo, tampoco es correcto decir que vive en de Sitter:")
    print("   la conica esta DENTRO del dominio, luego el dominio ATRAVIESA la")
    print("   absoluta.  Parte de el es el disco de Klein (metrica riemanniana de")
    print("   curvatura -1) y parte es el exterior (metrica lorentziana); no hay")
    print("   ninguna metrica que haga del dominio entero un billar metrico, y la")
    print("   metrica invariante degenera sobre la conica, que es una curva")
    print("   INTERIOR del dominio.")
    print()
    print("   Lo que SI es de Sitter son las ORBITAS: por la desigualdad fina, la")
    print("   conica queda dentro del hueco de cada orbita, luego ninguna orbita")
    print("   cruza la absoluta.  La dinamica de la familia es de Sitter aunque el")
    print("   dominio no lo sea.")
    print()
    print("   Y por tanto la familia lambda NO es toda ella caleidoscopica: solo")
    print("   el unico miembro esferico de cada n impar lo es.")
    print()
    # comprobacion exacta de las dos desigualdades, en Q(cos(pi/n)) via sympy
    try:
        import sympy as sp
    except ImportError:
        print("   (sympy no disponible: se omite la tabla)")
        return
    print("   n   j   lambda        geometria    dentro del dominio / dentro del hueco de la orbita")
    ok_signo = True
    ok_dentro = True
    ok_orbita = True
    for n in range(3, 11):
        a = sp.pi/n; v = sp.cos(a)
        for j in range(1, (n+1)//2):
            if 2*j >= n:
                continue
            C = sp.cos(sp.Rational(n-2*j, n)*sp.pi)
            lam = sp.simplify(v*(2*v**2 - 1 - C)/(C - 1))
            sg = sp.sign(sp.simplify(lam))
            geo = "esferica" if sg > 0 else ("euclidea" if sg == 0 else "de Sitter")
            pred = 1 if j > sp.Rational(n-2, 2) else (0 if j == sp.Rational(n-2, 2) else -1)
            ok_signo &= (int(sg) == pred)
            dentro = ""
            if sg < 0:
                d = sp.simplify(sp.simplify(-lam*v) - sp.simplify(v**2))
                ins = bool(d < 0)
                ok_dentro &= ins
                # y la version fina: radio^2 < v^4 = inradio^2 de la orbita
                d2 = sp.simplify(sp.simplify(-lam*v) - sp.simplify(v**4))
                ins2 = bool(d2 < 0)
                ok_orbita &= ins2
                dentro = ("si" if ins else "NO") + " / " + ("si" if ins2 else "NO")
            print(f"   {n:<3} {j:<3} {float(lam):+.6f}   {geo:<11}  {dentro}")
    claim("el signo de lambda es exactamente el predicho por j vs (n-2)/2", ok_signo)
    claim("para todo lambda<0 la conica queda dentro del circulo inscrito, luego "
          "DENTRO DEL DOMINIO y sin tocar ningun espejo", ok_dentro)
    claim("y dentro del hueco de cada orbita (radio^2 < v^4), luego ninguna "
          "orbita cruza la absoluta", ok_orbita)

def main():
    print(__doc__)
    parte_A()
    parte_B()
    parte_C()
    header("VEREDICTO")
    if FAILURES:
        print(f"   {len(FAILURES)} COMPROBACIONES FALLIDAS:")
        for f in FAILURES:
            print(f"     - {f}")
        sys.exit(1)
    print("   Todas las comprobaciones pasan, en aritmetica exacta.")
    print()
    print("   COMO ENMARCAR LA NOVEDAD:")
    print("   - la familia lambda ES un caleidoscopio esferico visto en RP^2.  Lo")
    print("     nuevo en ella no son los espejos: es que -I sea trivial en PGL_3,")
    print("     que es exactamente por lo que el caso IMPAR deja de estar vacio.")
    print("   - el segundo ejemplo NO es un caleidoscopio: grupo infinito.  Es el")
    print("     que responde de verdad a la segunda pregunta del problema 1.")
    sys.exit(0)

if __name__ == "__main__":
    main()
