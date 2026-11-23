"""Exact evaluator (fractions.Fraction only) for finite product distributions and explicit menus.

A distribution is a dict {value: probability}; an instance is a list of them.
  srev(D)                 exact global separate revenue (maximum over every support value)
  brev(D)                 exact global bundle revenue (exact convolution, maximum over every sum)
  menu_revenue(D, menu)   exact revenue of a finite menu [(allocation vector, payment), ...];
                          the buyer picks a utility-maximising entry (ties: larger payment)
  check_direct(D, alloc, pay)   exact check of EVERY IC and IR constraint of a direct mechanism on
                          the support; returns (ok, revenue)
Run the self tests:  python evaluator.py
"""
from fractions import Fraction as Fr
from itertools import product


def check_dist(F):
    assert all(v >= 0 and p >= 0 for v, p in F.items()) and sum(F.values()) == 1
    return True


def rev1(F):
    best = Fr(0)
    for p in F:
        best = max(best, p * sum(q for v, q in F.items() if v >= p))
    return best


def srev(D):
    return sum((rev1(F) for F in D), Fr(0))


def convolve(F, G):
    H = {}
    for a, p in F.items():
        for b, q in G.items():
            H[a + b] = H.get(a + b, Fr(0)) + p * q
    return H


def total(D):
    Z = {Fr(0): Fr(1)}
    for F in D:
        Z = convolve(Z, {v: p for v, p in F.items() if p})
    return Z


def power(F, n):
    R = {Fr(0): Fr(1)}
    while n:
        if n & 1:
            R = convolve(R, F)
        n >>= 1
        if n:
            F = convolve(F, F)
    return R


def rev_sum(Z):
    """exact monopoly revenue of a one dimensional law (sorted suffix sums)"""
    best = tail = Fr(0)
    for v in sorted(Z, reverse=True):
        tail += Z[v]
        best = max(best, v * tail)
    return best


def brev(D):
    return rev_sum(total(D))


def menu_revenue(D, menu):
    menu = [((Fr(0),) * len(D), Fr(0))] + [(tuple(map(Fr, a)), Fr(p)) for a, p in menu]
    R = Fr(0)
    for cell in product(*[list(F.items()) for F in D]):
        pr = Fr(1)
        for _, p in cell:
            pr *= p
        v = [x for x, _ in cell]
        R += pr * max((sum(a * x for a, x in zip(al, v)) - pay, pay) for al, pay in menu)[1]
    return R


def check_direct(D, alloc, pay):
    types = list(product(*[list(F) for F in D]))
    dot = lambda x, v: sum(a * b for a, b in zip(x, v))
    for t in types:
        if any(not 0 <= a <= 1 for a in alloc[t]) or dot(alloc[t], t) - pay[t] < 0:
            return False, None
        for s in types:
            if dot(alloc[t], t) - pay[t] < dot(alloc[s], t) - pay[s]:
                return False, None
    R = Fr(0)
    for t in types:
        pr = Fr(1)
        for F, x in zip(D, t):
            pr *= F[x]
        R += pr * pay[t]
    return True, R


if __name__ == '__main__':
    # Hart-Nisan style test: two iid items on {1,2}, prob 1/2 each
    F = {Fr(1): Fr(1, 2), Fr(2): Fr(1, 2)}
    check_dist(F)
    D = [F, F]
    assert srev(D) == 2 and brev(D) == Fr(9, 4)
    assert menu_revenue(D, [((1, 1), 3)]) == Fr(9, 4)
    al = {}
    py = {}
    for t in product(F, F):
        b = sum(t) >= 3
        al[t] = (int(b),) * 2
        py[t] = Fr(3 * b)
    assert check_direct(D, al, py) == (True, Fr(9, 4))
    py[(Fr(2), Fr(2))] = Fr(4)
    assert check_direct(D, al, py)[0] is False
    print('evaluator self tests ok')
