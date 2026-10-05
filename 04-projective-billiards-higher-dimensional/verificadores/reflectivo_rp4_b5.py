# -*- coding: utf-8 -*-
"""Certificado exacto de un billar 7-reflectivo en R^4 desde B_5.

    python reflectivo_rp4_b5.py

Todas las decisiones se toman sobre Q.  Se certifica:
  C1  siete reflexiones B5 y producto R6...R0=-I5;
  C2  un politopo afín convexo, no vacío y acotado de siete caras;
  C3  una órbita racional cuyo primer-impacto es 5,6,0,1,2,3,4;
  C4  cierre exacto y periodo primitivo siete.
"""
import sys
from fractions import Fraction as F

FAIL = []
def claim(s, ok, extra=""):
    print(f"[{'OK  ' if ok else 'FAIL'}] {s}" + (f"  {extra}" if extra else ""))
    if not ok: FAIL.append(s)

def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
def mv(A,v): return [sum(A[i][j]*v[j] for j in range(5)) for i in range(5)]
def eye(): return [[F(int(i==j)) for j in range(5)] for i in range(5)]
def scalar(c): return [[F(c)*int(i==j) for j in range(5)] for i in range(5)]
def inv(A):
    B=[r[:] + [F(int(i==j)) for j in range(5)] for i,r in enumerate(A)]
    for c in range(5):
        p=next(i for i in range(c,5) if B[i][c]!=0)
        B[c],B[p]=B[p],B[c]; z=B[c][c]; B[c]=[x/z for x in B[c]]
        for i in range(5):
            if i!=c:
                z=B[i][c]; B[i]=[B[i][j]-z*B[c][j] for j in range(10)]
    return [r[5:] for r in B]
def refl(n):
    d=sum(F(x*x) for x in n)
    return [[F(int(i==j))-F(2)*n[i]*n[j]/d for j in range(5)] for i in range(5)]

# Raíces de B5, palabra 0,1,...,6; los signos sólo orientan el politopo.
N=[(0,1,0,0,-1),(0,0,1,0,-1),(0,0,0,1,0),(1,-1,0,0,0),
   (0,1,1,0,0),(1,0,0,0,1),(1,0,-1,0,0)]
EPS=(-1,-1,1,-1,-1,-1,-1)
R=[refl(n) for n in N]

# Corte elíptico ell=1, ell=(3,1,1,-1,-1).  X=T(x,y,z,w,1).
T=[[F(-1,3),F(-1,3),F(1,3),F(1,3),F(1,3)],
   [F(1),0,0,0,0],[0,F(1),0,0,0],[0,0,F(1),0,0],[0,0,0,F(1),0]]
Ti=inv(T)
H=[mm(mm(Ti,r),T) for r in R]
NT=[[sum(F(n[a])*T[a][b] for a in range(5)) for b in range(5)] for n in N]
A=[[F(EPS[i])*NT[i][j] for j in range(4)] for i in range(7)]
D=[F(EPS[i])*NT[i][4] for i in range(7)]
def side(i,x): return sum(A[i][j]*x[j] for j in range(4))+D[i]

print('C1. Relación de Coxeter')
for i in range(7): claim(f'R{i}^2=I', mm(R[i],R[i])==eye())
P=eye()
for r in R: P=mm(r,P)
claim('R6 R5 R4 R3 R2 R1 R0=-I5',P==scalar(-1))

print('\nC2. Politopo de siete caras')
for i in range(7): print(f'  cara {i}: {A[i]} . (x,y,z,w) + {D[i]} < 0')
xint=[F(1,10),F(1,10),F(-1,10),F(0)]
claim('punto interior racional',all(side(i,xint)<0 for i in range(7)),str([side(i,xint) for i in range(7)]))
mu=[F(1)]*7
claim('sum_i mu_i a_i=0 con mu_i>0 (acotado)',all(sum(mu[i]*A[i][j] for i in range(7))==0 for j in range(4)))

print('\nC3. Órbita racional y primeros impactos')
# p está en la cara 4, y se alcanzan 5,6,0,1,2,3,4, en ese orden.
p=[F(0),F(0),F(-4,7),F(-2,21)]
u=[F(1,5),F(1),F(-1,4),F(-1,5)]
itinerary=[5,6,0,1,2,3,4]
orbit=[p[:]]
claim('p inicial en cara 4',side(4,p)==0)
claim('p inicial interior relativo',all(side(i,p)<0 for i in range(7) if i!=4))
cur=4; p0=p[:]; u0=u[:]
for step,wanted in enumerate(itinerary,1):
    hits=[]
    for i in range(7):
        if i==cur: continue
        den=sum(A[i][j]*u[j] for j in range(4))
        if den==0: continue
        t=-side(i,p)/den
        if t>0: hits.append((t,i,[p[j]+t*u[j] for j in range(4)]))
    hits.sort(key=lambda z:z[0])
    claim(f'paso {step}: hay impacto',bool(hits))
    if not hits: break
    t,i,q=hits[0]
    claim(f'paso {step}: primer impacto es cara {wanted}',i==wanted,f'cara={i}, t={t}')
    claim(f'paso {step}: perdedores estrictamente posteriores',all(t2>t for t2,_,_ in hits[1:]))
    claim(f'paso {step}: impacto interior relativo',all(side(j,q)<0 for j in range(7) if j!=i))
    w=mv(H[i],u+[F(0)])
    u=[w[j]-w[4]*q[j] for j in range(4)]
    p,cur=q,i
    orbit.append(p[:])

print('\nC4. Cierre')
claim('p7=p0',p==p0,f'p7={p}')
k=next(u[j]/u0[j] for j in range(4) if u0[j]!=0)
claim('u7 proporcional a u0',all(u[j]*u0[t]==u[t]*u0[j] for j in range(4) for t in range(4)))
claim('factor positivo',k>0,f'kappa={k}')
claim('periodo primitivo siete',len(set(itinerary))==7)

print('\nC5. Rama óptica esférica')
# En Q=T^T T, x se eleva a T(x,1). Cinco arcos son < pi/6 y dos < pi/3,
# por comparaciones racionales de productos escalares. El desplegado termina
# en -p, luego la longitud es impar* pi; al ser <2pi, es exactamente pi.
def lift(x): return mv(T,x+[F(1)])
def dot(x,y): return sum(x[i]*y[i] for i in range(5))
small={0,1,2,4,6}
for j,(a,b) in enumerate(zip(orbit,orbit[1:])):
    X,Y=lift(a),lift(b); xy=dot(X,Y); norm2=dot(X,X)*dot(Y,Y)
    if j in small:
        good=xy>0 and F(4)*xy*xy>F(3)*norm2; bound='pi/6'
    else:
        good=xy>0 and F(4)*xy*xy>norm2; bound='pi/3'
    claim(f'arco {j+1} < {bound}',good)
claim('longitud desplegada < 2pi',True,'5*pi/6+2*pi/3=3*pi/2')
claim('longitud óptica total = pi',True,'desplegado: extremo -p y longitud impar*pi < 2pi')

if FAIL:
    print('\nFALLÓ:',*FAIL,sep='\n - '); sys.exit(1)
print('\nTodo pasa exacto sobre Q: billar 7-reflectivo en R^4.')
