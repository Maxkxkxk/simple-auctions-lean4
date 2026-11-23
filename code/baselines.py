"""BASELINE reproductions (not new results).  Run:  python baselines.py
1. Rubinstein-type (two groups, limit 2):  exact finite instances + the limit value.
2. Huang (arXiv 2610.24973): exact two-item seed (2.02) and the 2.11 compound-Poisson seed
   (continuous; transcendental constants -> decimal evaluation of the closed forms).
3. Giannakopoulos-Hahn (arXiv 2609.35711, Thm D.1): the A-group with N_a = 1000, evaluated with
   rigorous directed rounding (integer fixed point, lower bounds of all pmf cells).
"""
from fractions import Fraction as Fr
from math import exp
import evaluator as E
try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


# ---------- 1. Rubinstein-type ----------
def rubinstein(L):
    """L^2 iid ER items with mean total 1, plus one rare item making SREV = 1 (value 1000)."""
    H = sum(Fr(1, j) for j in range(1, L + 1))
    d = 1 / (L * L * H)
    F = {j * d: (Fr(1, j) - Fr(1, j + 1) if j < L else Fr(1, L)) for j in range(1, L + 1)}
    E.check_dist(F)
    ZB = E.power(F, L * L)
    r = 1 - 1 / H
    M = Fr(1000)
    rare = {Fr(0): 1 - r / M, M: r / M}
    S = L * L * E.rev1(F) + E.rev1(rare)
    B = E.rev_sum(E.convolve(ZB, rare))
    R = E.rev_sum(ZB) + r          # price the concentrated group and the rare item separately
    return R / max(S, B)


# ---------- 2. Huang ----------
def huang_two_item():
    F = {Fr(0): Fr(3, 5), Fr(1): Fr(1, 5), Fr(2): Fr(1, 15), Fr(3): Fr(2, 15)}
    E.check_dist(F)
    s = E.srev([F, F])
    r = E.menu_revenue([F, F], [((1, 1), 3)])
    Z = E.convolve(E.total([F, F]), {Fr(2): Fr(1)})
    assert E.rev_sum(Z) == 2 and s == Fr(4, 5) and r == Fr(21, 25)
    return 2 + (r - s) / 2


def huang_211():
    Q = lambda t: 47 / 144 * exp(-t) + 11 / 18 * exp(-4 * t)
    revseed = .5 * Q(.5)
    n = 200000
    dx = 40 / n
    ell = lambda x: (exp(-x) - exp(-4 * x / 3) + exp(-4 * x) - exp(-48 * x)) / x
    tail = 0.
    best = 0.
    for k in range(n, 0, -1):
        x = (k - .5) * dx
        tail += ell(x) * dx
        best = max(best, (k - 1) * dx * tail)
    return revseed, best, 2 + (revseed - best) / .25, 2 + (Fr(7, 50) - Fr(9, 80)) / Fr(1, 4)


# ---------- 3. Giannakopoulos-Hahn ----------
BITS = 96
S = 2 * BITS + 24


def pack(co):
    if len(co) == 1:
        return mpz(co[0])
    h = len(co) // 2
    return pack(co[:h]) + (pack(co[h:]) << (S * h))


def unpack(F, n):
    m = (mpz(1) << S) - 1
    return [int((F >> (S * i)) & m) for i in range(n)]


def mulr(a, b, n, up):
    c = unpack(pack(a) * pack(b), n)
    return [(x + (1 << BITS) - 1 >> BITS) if up else x >> BITS for x in c]


def powr(a, k, n, up):
    r = [1 << BITS] + [0] * (n - 1)
    while k:
        if k & 1:
            r = mulr(r, a, n, up)
        k >>= 1
        if k:
            a = mulr(a, a, n, up)
    return r


def rnd(x, up):
    y = x * (1 << BITS)
    return -((-y.numerator) // y.denominator) if up else y.numerator // y.denominator


def gh():
    Na = 1000
    N = 185                      # S_A = 4.625 in units of 1/40
    e = 280
    groups = [(0, 5, 280, Fr(12)), (40, 40, 320, Fr(38)), (600, 600, 880, Fr(135))]
    # tail count r/k per group (r in lattice units): 3/10*40=12, 19/20*40=38, 27/8*40=135
    T = 2000
    laws = []
    for p, a, b, r in groups:
        t = lambda k: r / (Na * max(k, a)) if k <= b else Fr(0)
        f = [Fr(0)] * (b + 1)
        for k in range(a, b + 1):
            f[k] = t(k) - t(k + 1)
        f[0] = 1 - sum(f)
        assert min(f) >= 0
        laws.append((p, f))
    pad = lambda x, n: (x + [0] * n)[:n]
    Z = [1 << BITS] + [0] * (T - 1)
    for p, f in laws:
        Z = mulr(Z, powr(pad([rnd(x, False) for x in f], T), Na, T, False), T, False)
    C = 0
    worst = Fr(0)
    for t in range(T):
        if t:
            worst = max(worst, Fr((N + t) * ((1 << BITS) - C), N << BITS))
        C += Z[t]
    # t >= T: (N+t) Pr[Z>=t] <= E[(N+Z) 1{Z>=T}] = N + E Z - sum_{z<T} (N+z) Pr[Z=z]
    EZ = sum((Na * sum(k * x for k, x in enumerate(f)) for p, f in laws), Fr(0))
    tail = (N + EZ - Fr(sum((N + z) * x for z, x in enumerate(Z)), 1 << BITS)) / N
    # revenue, upper bounds of the loss terms
    g, h = [], []
    for p, f in laws:
        g.append(pad([rnd(sum(f[:p + 1]), True)] + [rnd(x, True) for x in f[p + 1:]], e))
        h.append(pad([rnd(x, True) for x in f[p:]], e))
        if p == 0:
            h[-1] = [0] * e
    pw = [powr(x, Na, e, True) for x in g]
    pw1 = [powr(x, Na - 1, e, True) for x in g]
    G = [1 << BITS] + [0] * (e - 1)
    for x in pw:
        G = mulr(G, x, e, True)
    loss = Fr(0)
    for i, (p, f) in enumerate(laws):
        q = mulr(h[i], pw1[i], e, True)
        for j in range(3):
            if j != i:
                q = mulr(q, pw[j], e, True)
        loss += Fr(p * Na * sum(q), 1 << BITS)
    pay = sum((Fr(p) * Na * sum(f[p:]) for p, f in laws if p), Fr(0))
    rev = e * (1 - Fr(sum(G), 1 << BITS)) + pay - loss
    return worst, tail, rev / N, 1 + rev / N


if __name__ == '__main__':
    print('1. Rubinstein-type, exact ratios:', {L: round(float(rubinstein(L)), 5) for L in (2, 3, 4, 5, 6)},
          '(-> 2 as L -> infinity)')
    print('2. Huang two-item seed, exact sup ratio =', huang_two_item())
    a, b, c, lo = huang_211()
    print('   Huang 2.11 seed: REVseed=%.6f SREVseed=%.6f ratio=%.6f, proved bound=%s' % (a, b, c, float(lo)))
    w, tl, r, ratio = gh()
    print('3. Giannakopoulos-Hahn: max_t (S+t)P(Z>=t)/S <= %.9f (t < 2000), tail bound (t >= 2000) %.9f, REV_A/S_A >= %.9f, ratio >= %.9f'
          % (float(w), float(tl), float(r), float(ratio)))
    print('   published: 2 + 0.88/4.625 = %.9f' % (2 + .88 / 4.625))
