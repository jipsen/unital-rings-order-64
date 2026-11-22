#!/usr/bin/env python3
"""
Completeness check by the mass formula

    sum_{A in list} f(A) / |Aut(A)|  ==  T(n) / 2^(n(n-1)/2)

f(A) = number of complete flags A = V_1 > V_2 > ... > V_n > 0 with
A*V_i + V_i*A <= V_{i+1}; T(n) = number of strictly triangular associative
tables (from count_tables.c); |Aut(A)| from ModIsom (autsizes_n.txt).
A missing class makes the left side too small, a repeated class too large.
"""
import json
import subprocess
import sys
from fractions import Fraction

sys.path.insert(0, ".")
from enum_nilpotent import prod, rref, automorphisms, bits

data = json.load(open("nilpotent_F2_upto5.json"))


def span_elems(basis):
    out = [0]
    for b in basis:
        out += [x ^ b for x in out]
    return out


def adapted_flags(mult, V):
    """Number of adapted flag tails starting at subspace V (list of basis vecs)."""
    if not V:
        return 1
    n = len(mult)
    W = rref([prod(mult, 1 << i, v) for i in range(n) for v in V] +
             [prod(mult, v, 1 << i) for i in range(n) for v in V])
    # extend W to a basis of V: W + U
    U = []
    cur = list(W)
    for v in V:
        if len(rref(cur + [v])) > len(cur):
            cur = list(rref(cur + [v]))
            U.append(v)
    b = len(U)
    total = 0
    for c in range(1, 1 << b):  # nonzero functional on span(U)
        # kernel of c on span(U): vectors x in F_2^b with <c,x> = 0
        ker = [x for x in range(1 << b) if bin(x & c).count("1") % 2 == 0]
        kb = rref([_comb(U, x) for x in ker])
        H = rref(list(W) + list(kb))
        assert len(H) == len(V) - 1
        total += adapted_flags(mult, list(H))
    return total


def _comb(U, x):
    r = 0
    for t in bits(x):
        r ^= U[t]
    return r


ok = True
for n in sorted(data, key=int):
    n_int = int(n)
    algs = data[n]
    sizes = eval(open(f"autsizes_{n}.txt").read().replace("\n", ""))
    assert len(sizes) == len(algs)
    if n_int <= 4:  # independent brute-force |Aut| agrees with ModIsom?
        bf = [len(automorphisms(A)) for A in algs]
        assert bf == sizes, (n, bf, sizes)
    lhs = sum(Fraction(adapted_flags(A, [1 << i for i in range(n_int)]), s)
              for A, s in zip(algs, sizes))
    T = int(subprocess.check_output(["./count_tables", n]).decode())
    rhs = Fraction(T, 2 ** (n_int * (n_int - 1) // 2))
    print(f"dim {n}: classes {len(algs)}, sum f/|Aut| = {lhs}, T/|Borel| = {rhs}, "
          f"{'OK' if lhs == rhs else 'MISMATCH'}")
    ok &= lhs == rhs
print("ALL OK" if ok else "FAILED")
