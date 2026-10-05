"""Exact second-order obstruction at the paper 04 rank-five example.

Independent of the existing Jacobian checkers. A nonzero quadratic obstruction
on ker(dPhi) distinguishes Zariski tangent dimension from local dimension.
"""
import sympy as s
from pathlib import Path
import json, random

roots=list(map(s.Matrix,[(0,1,0,0,-1),(0,0,1,0,-1),(0,0,0,1,0),
    (1,-1,0,0,0),(0,1,1,0,0),(1,0,0,0,1),(1,0,-1,0,0)]))
m=5; E=s.eye(m); zero=s.zeros(m); parameters=[]; H=[]; derivative=[]
for i,n in enumerate(roots):
    q=(n.T*n)[0]; H.append(E-2*n*n.T/q)
    pivot=next(j for j in range(m) if n[j])
    ds=[]
    for kind in ('pole','axis'):
        for j in range(m):
            if j==pivot: continue
            e=E[:,j]
            dp=e if kind=='pole' else s.zeros(m,1)
            dl=e.T if kind=='axis' else s.zeros(1,m)
            q1=(n.T*dp+dl*n)[0]
            ds.append(-2*((dp*n.T+n*dl)/q-n*n.T*q1/q**2))
            parameters.append((i,kind,j))
    derivative.append(ds)
product=E
for a in H: product=a*product
assert product==-E
columns=[]
for i in range(7):
    for d in derivative[i]:
        p=E
        for j in range(7): p=(d if j==i else H[j])*p
        columns.append(s.Matrix(list(p)))
J=s.Matrix.hstack(*columns)
rank=J.rank(); basis=J.nullspace()
assert rank==23 and len(basis)==33
# Trace pairing with a traceless commutant of the unfolded mirror group.
Y=s.diag(1,1,1,-4,1)
assert all(Y*a==a*Y for a in H)
assert all(s.trace(Y*s.Matrix(m,m,v))==0 for v in columns)

def obstruction(v):
    dp=[s.zeros(m,1) for _ in roots];dl=[s.zeros(1,m) for _ in roots]
    for value,(i,kind,j) in zip(v,parameters):
        if kind=='pole':dp[i][j]=value
        else:dl[i][j]=value
    P0,P1,P2=E,zero,zero
    for n,p,l,h in zip(roots,dp,dl,H):
        q0=(n.T*n)[0];q1=(n.T*p+l*n)[0];q2=(l*p)[0]
        a0=n*n.T;a1=p*n.T+n*l;a2=p*l
        h1=-2*(a1/q0-a0*q1/q0**2)
        h2=-2*(a2/q0-a1*q1/q0**2+a0*(q1**2/q0**3-q2/q0**2))
        P0,P1,P2=h*P0,h*P1+h1*P0,h*P2+h1*P1+h2*P0
    assert P1==zero
    return s.trace(Y*P2)

trials=[]
for j in range(len(basis)):
    trials.append([int(i==j) for i in range(len(basis))])
for j in range(len(basis)):
    for k in range(j):
        trials.append([int(i in (j,k)) for i in range(len(basis))])
for coefficients in trials:
    v=sum((c*b for c,b in zip(coefficients,basis)),s.zeros(56,1))
    q=obstruction(v)
    if q: break
else: raise AssertionError('No quadratic obstruction found; no dimension verdict justified.')
result=dict(jacobian_rank=rank,tangent_dimension=56-rank,commutant=list(Y.diagonal()),
    nullspace_coefficients=coefficients,tangent_vector=list(v),quadratic_obstruction=q,
    conclusion='Nonzero quadratic equation after eliminating 23 regular equations: local dimension at most 32, not at least 33.')
t=s.Symbol('t')
pt=[n.copy() for n in roots]; lt=[n.T.copy() for n in roots]
for value,(i,kind,j) in zip(v,parameters):
    if kind=='pole':pt[i][j]+=t*value
    else:lt[i][j]+=t*value
exact=E
for p,l in zip(pt,lt):
    exact=((E-2*p*l/(l*p)[0])*exact).applyfunc(s.cancel)
scalar=s.factor(s.trace(Y*(exact+E)))
assert s.limit(scalar/t**2,t,0)==q
result['direct_rational_trace']=str(scalar)
result['nonzero_parameter_changes']=[(i,kind,j,str(value)) for value,(i,kind,j) in zip(v,parameters) if value]
Path(__file__).with_name('corank-obstruction.json').write_text(json.dumps(result,indent=2,default=str),encoding='utf8')
print(json.dumps(result,indent=2,default=str))
