# -*- coding: utf-8 -*-
"""Certificado exacto: billar 8-reflectivo en R^5 (familia B6).

    python reflectivo_rp5_b6.py

Comprueba con fracciones: reflexión, producto escalar, politopo acotado,
itinerario de primeros impactos, interioridad y cierre. No usa flotantes.
"""
import sys
from fractions import Fraction as F

M, K = 6, 8
BAD=[]
def ok(s,b,detail=""):
    print(f"[{'OK  ' if b else 'FAIL'}] {s}"+(f"  {detail}" if detail else ""))
    if not b: BAD.append(s)
def mm(A,B): return [[sum(A[i][h]*B[h][j] for h in range(M)) for j in range(M)] for i in range(M)]
def mv(A,v): return [sum(A[i][j]*v[j] for j in range(M)) for i in range(M)]
def eye(c=1): return [[F(c)*int(i==j) for j in range(M)] for i in range(M)]
def inverse(A):
    B=[r[:] + [F(int(i==j)) for j in range(M)] for i,r in enumerate(A)]
    for c in range(M):
        p=next(i for i in range(c,M) if B[i][c])
        B[c],B[p]=B[p],B[c]; z=B[c][c]; B[c]=[x/z for x in B[c]]
        for i in range(M):
            if i!=c:
                z=B[i][c]; B[i]=[B[i][j]-z*B[c][j] for j in range(2*M)]
    return [r[M:] for r in B]
def reflection(n):
    q=sum(F(x*x) for x in n)
    return [[F(int(i==j))-F(2)*n[i]*n[j]/q for j in range(M)] for i in range(M)]

# Las siete raíces B5 de la construcción de RP4, más e_6.
N=[(0,1,0,0,-1,0),(0,0,1,0,-1,0),(0,0,0,1,0,0),(1,-1,0,0,0,0),
   (0,1,1,0,0,0),(1,0,0,0,1,0),(1,0,-1,0,0,0),(0,0,0,0,0,1)]
EPS=(1,1,-1,1,-1,-1,1,-1)
R=[reflection(n) for n in N]

# Corte ell=1, ell=(-1,1,1,1,3,1): X=T(x_1,...,x_5,1).
T=[[F(1),0,0,0,0,0],[0,F(1),0,0,0,0],[0,0,F(1),0,0,0],
   [0,0,0,F(1),0,0],[F(1,3),F(-1,3),F(-1,3),F(-1,3),F(-1,3),F(1,3)],
   [0,0,0,0,F(1),0]]
Ti=inverse(T); H=[mm(mm(Ti,r),T) for r in R]
NT=[[sum(F(n[a])*T[a][b] for a in range(M)) for b in range(M)] for n in N]
A=[[F(EPS[i])*NT[i][j] for j in range(M-1)] for i in range(K)]
D=[F(EPS[i])*NT[i][M-1] for i in range(K)]
def side(i,x): return sum(A[i][j]*x[j] for j in range(M-1))+D[i]

print('C1. Reflexiones y monodromía')
for i in range(K): ok(f'R{i}^2=I',mm(R[i],R[i])==eye())
P=eye()
for r in R: P=mm(r,P)
ok('R7 R6 ... R0=-I6',P==eye(-1))

print('\nC2. Dominio convexo acotado')
for i in range(K): print(f'  {i}: {A[i]} . x + {D[i]} < 0')
xint=[F(-1,10),F(1,10),F(1,10),F(1,10),F(1,10)]
ok('punto interior racional',all(side(i,xint)<0 for i in range(K)),str([side(i,xint) for i in range(K)]))
ok('dependencia positiva de normales (todos mu_i=1)',all(sum(A[i][j] for i in range(K))==0 for j in range(M-1)))

print('\nC3. Trayectoria racional')
p=[F(-2,33),F(0),F(0),F(4,11),F(4,11)]
u=[F(-11,10),F(1),F(-2,5),F(7,10),F(1,5)]
path=[5,6,7,0,1,2,3,4]; cur=4; p0=p[:]; u0=u[:]
ok('p0 en la cara 4',side(4,p)==0)
ok('p0 en interior relativo',all(side(i,p)<0 for i in range(K) if i!=4))
for z,wanted in enumerate(path,1):
    hits=[]
    for i in range(K):
        if i==cur: continue
        den=sum(A[i][j]*u[j] for j in range(M-1))
        if not den: continue
        t=-side(i,p)/den
        if t>0: hits.append((t,i,[p[j]+t*u[j] for j in range(M-1)]))
    hits.sort(key=lambda x:x[0]); ok(f'paso {z}: existe impacto',bool(hits))
    if not hits: break
    t,i,q=hits[0]
    ok(f'paso {z}: primera cara {wanted}',i==wanted,f'cara={i}, t={t}')
    ok(f'paso {z}: perdedores estrictos',all(t2>t for t2,_,_ in hits[1:]))
    ok(f'paso {z}: impacto interior relativo',all(side(j,q)<0 for j in range(K) if j!=i))
    v=mv(H[i],u+[F(0)]); u=[v[j]-v[M-1]*q[j] for j in range(M-1)]
    p,cur=q,i

print('\nC4. Cierre')
ok('p8=p0',p==p0,f'p8={p}')
kappa=next(u[j]/u0[j] for j in range(M-1) if u0[j])
ok('u8 proporcional a u0',all(u[i]*u0[j]==u[j]*u0[i] for i in range(M-1) for j in range(M-1)))
ok('factor positivo',kappa>0,f'kappa={kappa}')
ok('periodo primitivo ocho',len(set(path))==K)
if BAD:
    print('FALLÓ:',*BAD,sep='\n - ');sys.exit(1)
print('Todo pasa exacto sobre Q: billar 8-reflectivo en R^5.')
