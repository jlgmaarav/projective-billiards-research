# -*- coding: utf-8 -*-
"""Exact short-branch certificate for the RP3 D4 atom.

The developed lift has product ``-I``.  We compute its six consecutive
spherical chords over Q, use outward rational square-root brackets, and prove
that their sum is < 4.  Since every developed cosine is positive, every arc
is the minor arc; the standard chord bound then gives L < 2*pi, hence
L = pi (mod 2*pi) forces exactly one turn.
"""
from fractions import Fraction as F
from math import isqrt

import reflectivo_rp3_d4 as d4


PRECISION = 10**40


def sqrt_sup(q):
    if q == 0:
        return F(0)
    assert q > 0
    n, den = q.numerator, q.denominator
    z = isqrt(n * den * PRECISION * PRECISION)
    return F(z + 1, PRECISION * den)


def inverse(mx):
    b = [row[:] + [F(int(i == j)) for j in range(4)]
         for i, row in enumerate(mx)]
    for c in range(4):
        p = next(i for i in range(c, 4) if b[i][c] != 0)
        b[c], b[p] = b[p], b[c]
        z = b[c][c]
        b[c] = [x/z for x in b[c]]
        for i in range(4):
            if i != c and b[i][c] != 0:
                q = b[i][c]
                b[i] = [b[i][j] - q*b[c][j] for j in range(8)]
    return [row[4:] for row in b]


def simulate():
    p = [F(4, 7), F(4, 5), F(-4, 7)]
    u = [F(15), F(5), F(3)]
    face = 0
    points = [p[:]]
    for _ in range(d4.K):
        candidates = []
        for i in range(d4.K):
            if i == face:
                continue
            denominator = sum(d4.A3[i][j] * u[j] for j in range(3))
            if denominator == 0:
                continue
            step = -d4.lado(i, p) / denominator
            if step > 0:
                q = [p[j] + step*u[j] for j in range(3)]
                candidates.append((step, i, q))
        step, face, q = min(candidates)
        reflected = d4.mv(d4.R[face], [u[0], u[1], u[2], F(0)])
        u = [reflected[j] - reflected[3]*q[j] for j in range(3)]
        p = q
        points.append(p[:])
    assert points[-1] == points[0]
    return points


def developed(points):
    # The initial impact is on face 0; the next five faces are the first five
    # entries of the strict itinerary.  The sixth reflection gives X_6=-X_0.
    reflected_faces = [0] + d4.ITIN[:5]
    product = d4.eye()
    out = []
    for i, p in enumerate(points[:7]):
        if i:
            product = d4.mul(d4.R[reflected_faces[i - 1]], product)
        out.append(d4.mv(inverse(product), [p[0], p[1], p[2], F(1)]))
    assert out[-1] == [-x for x in out[0]]
    return out


def main():
    points = simulate()
    xs = developed(points)
    total = F(0)
    for i in range(6):
        a, b = xs[i], xs[i + 1]
        aa = sum(x*x for x in a)
        bb = sum(x*x for x in b)
        dot = sum(a[j]*b[j] for j in range(4))
        assert dot > 0
        cosine_lower = dot / sqrt_sup(aa * bb)
        assert cosine_lower > 0
        chord_sq_upper = 2 - 2*cosine_lower
        chord_upper = sqrt_sup(chord_sq_upper)
        total += chord_upper
        print(f"chord_{i}_upper={chord_upper}")
    assert total < 4
    print(f"sum_chords_upper={total}")
    print(f"sum_chords_upper_float={float(total):.12f}")
    print("UN GIRO CERTIFICADO: L < 2*pi y producto -I => L = pi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
