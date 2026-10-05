# -*- coding: utf-8 -*-
"""Revisión adversarial numérica del levantamiento por paridad.

Construye independientemente varios (n,k), comprueba rango, acotación, caras
genuinas, sección regular, producto +I y unipotencia. La demostración principal
es analítica; este script busca errores de signo, índice o dimensión.
"""
import sys
import numpy as np


FAIL = []
TOL = 2e-10


def claim(label, ok, extra=""):
    print(f"[{'OK  ' if ok else 'FAIL'}] {label}" + (f"  {extra}" if extra else ""))
    if not ok:
        FAIL.append(label)


def parity_basis(k, n):
    angles = np.array([2*i*np.pi/k for i in range(k)])
    cols = [np.cos(angles), np.sin(angles)]
    if n == 2:
        return np.column_stack(cols)
    for parity in (0, 1):
        indices = list(range(parity, k, 2))
        anchor = indices[-1]
        for j in indices[:-1]:
            e = np.zeros(k)
            e[j] = 1
            e[anchor] = -1
            trial = np.column_stack(cols+[e])
            if np.linalg.matrix_rank(trial, tol=1e-10) > len(cols):
                cols.append(e)
            if len(cols) == n:
                return np.column_stack(cols)
    if len(cols) != n:
        raise RuntimeError("no se completó la base de paridad")
    return np.column_stack(cols)


def build(n, k):
    B = parity_basis(k, n)
    N = B.copy()
    O = np.r_[np.zeros(n), 1.]
    H = []
    for row in N:
        L = np.r_[row, -1.]
        H.append(np.eye(n+1)+2*np.outer(O, L))
    return B, N, H


cases = [(2,4), (2,8), (3,6), (3,8), (4,6), (4,12),
         (5,8), (5,12), (6,8), (7,10), (8,16), (14,20)]
for n, k in cases:
    print(f"\n(n,k)=({n},{k})")
    B,N,H = build(n,k)
    r = k//2
    claim("dimensión del subespacio elegido", B.shape == (k,n))
    claim("cada paridad suma cero",
          np.max(np.abs(B[0::2].sum(axis=0))) < TOL and
          np.max(np.abs(B[1::2].sum(axis=0))) < TOL)
    claim("las normales generan todo el dual", np.linalg.matrix_rank(N, TOL) == n)
    claim("relación positiva total", np.max(np.abs(N.sum(axis=0))) < TOL)
    claim("suma alternada nula",
          np.max(np.abs(sum(((-1)**i)*N[i] for i in range(k)))) < TOL)

    # Every inequality is exposed by its parity-block Fourier point.
    exposed = True
    for i in range(k):
        q = np.array([np.cos(i*np.pi/r), np.sin(i*np.pi/r)])
        point = np.zeros(n)
        point[:2] = q
        vals = N@point
        exposed &= abs(vals[i]-1) < TOL
        exposed &= np.max(np.delete(vals, i)) < 1-1e-8
    claim("todas las desigualdades son caras genuinas", exposed)

    # The diagonal Fourier plane cuts out the regular k-gon.
    section_ok = True
    for i in range(k):
        q = np.array([0.371, -0.219])
        point = np.zeros(n)
        point[:2] = q
        expected = np.array([np.cos(j*np.pi/r)*q[0]+np.sin(j*np.pi/r)*q[1]
                             for j in range(k)])
        section_ok &= np.max(np.abs(N@point-expected)) < TOL
        break
    claim("la sección tiene las k normales regulares", section_ok)

    P = np.eye(n+1)
    for A in H:
        P = A@P
    claim("producto canónico +I", np.max(np.abs(P-np.eye(n+1))) < 5e-9)
    U = H[1]@H[0]-np.eye(n+1)
    claim("par de espejos unipotente no trivial",
          np.max(np.abs(U)) > 1e-8 and np.max(np.abs(U@U)) < 5e-9)

if FAIL:
    print("\nFALLÓ:", *FAIL, sep="\n - ")
    sys.exit(1)
print("\nTodos los casos pasan: no se detectan errores de signo o dimensión.")
