#!/usr/bin/env python3
r"""
Local rings R with residue field F_q (q = 4 or 8) of order q * |J|.

Reduction.  Let zeta in R be a unit of order q-1 lifting a fixed generator of
F_q^*  (it exists, and all such lifts are conjugate under 1+J, by
Schur-Zassenhaus).  zeta is the Teichmueller element of the Galois subring
S = Z[zeta] = GR(char R, r) and satisfies h(zeta) = 0, where h is the 2-adic
(Hensel) lift of the minimal polynomial of the generator (q=4: h = x^2+x+1).
R = S + J with S /\ J = 2S, so
        R  =  (Z[x]/(h)  |x  J) / { (2s, -s*c) },      c = 2*1_R in J,
and R is determined by the marked ring (J, c) together with the additive maps
        L(y) = zeta*y,   Rt(y) = y*zeta     on J,
subject to
        h(L) = h(Rt) = 0,  L Rt = Rt L,  L(c) = Rt(c),
        L(xy) = L(x) y,  Rt(xy) = x Rt(y),  Rt(x) y = x L(y).
Isomorphism classes of R  <->  orbits of Aut(J,c) x <Frobenius> on such pairs
(conjugation psi L psi^-1, and Frobenius (L, Rt) -> (L^2, Rt^2), i.e. zeta -> zeta^2).
"""
import json
import sys
from itertools import product

sys.path.insert(0, ".")
from enum_local_F2 import Ring, bits, generators

DATA = json.load(open("local_F2_upto5.json"))


def marked_rings(m):
    if m == 0:
        return [Ring(0, [], [], 0)]
    return [Ring(r["m"], r["P"], r["Q"], r["c"]) for r in DATA[str(m)]]


def hensel_factor(q):
    """Monic h (coeff list a0..a_{r-1}, leading 1 implicit) mod 4 with
    h | x^(q-1) - 1 mod 4 and h = minimal polynomial of the chosen generator mod 2."""
    if q == 4:
        return [1, 1]          # x^2 + x + 1 (exact over Z)
    if q == 8:
        target = [1, 1, 0]     # x^3 + x + 1 mod 2
        sols = []
        for a in product(range(4), repeat=3):
            if [x % 2 for x in a] != target:
                continue
            # divide x^7 - 1 by x^3 + a2 x^2 + a1 x + a0 mod 4
            num = [-1, 0, 0, 0, 0, 0, 0, 1]
            for d in range(7, 2, -1):
                lead = num[d] % 4
                for k in range(3):
                    num[d - 3 + k] -= lead * a[k]
                num[d] -= lead
            if all(x % 4 == 0 for x in num[:3]):
                sols.append(list(a))
        assert len(sols) == 1
        return sols[0]
    raise ValueError


def additive_endos(J):
    """All additive endomorphisms of J (as value tables), via generator images."""
    m, n = J.m, 1 << J.m
    out = []
    imgs = [0] * m

    def rec(i):
        if i < 0:
            out.append(tuple(J.apply(imgs, x) for x in range(n)))
            return
        for hh in range(n):
            if J.add[hh][hh] == J.apply(imgs, J.P[i]):
                imgs[i] = hh
                rec(i - 1)
        imgs[i] = 0

    rec(m - 1)
    return out


def classify(J, q, auts):
    n = 1 << J.m
    h = hensel_factor(q)
    r = len(h)
    A, M = J.add, J.mul

    def times(k, y):
        k %= 64
        z = 0
        for _ in range(k):
            z = A[z][y]
        return z

    def hpoly_zero(L):
        for x in range(n):
            powers = [x]
            for _ in range(r):
                powers.append(L[powers[-1]])
            v = powers[r]
            for k in range(r):
                v = A[v][times(h[k], powers[k])]
            if v:
                return False
        return True

    if n == 1:
        pairs = [((0,), (0,))]
    else:
        endos = additive_endos(J)
        Ls = [L for L in endos if hpoly_zero(L) and
              all(L[M[x][y]] == M[L[x]][y] for x in range(n) for y in range(n))]
        Rs = [R for R in endos if hpoly_zero(R) and
              all(R[M[x][y]] == M[x][R[y]] for x in range(n) for y in range(n))]
        pairs = []
        for L in Ls:
            for R in Rs:
                if L[J.c] != R[J.c]:
                    continue
                if any(L[R[x]] != R[L[x]] for x in range(n)):
                    continue
                if any(M[R[x]][y] != M[x][L[y]] for x in range(n) for y in range(n)):
                    continue
                pairs.append((L, R))

    # group: Aut(J,c) (as value tables) and Frobenius
    def table(phi):
        return tuple(J.apply(phi, x) for x in range(n)) if J.m else (0,)

    G = [table(phi) for phi in auts]
    gens = [table(phi) for phi in generators(J, auts)] if J.m else []
    inv = lambda p: tuple(sorted(range(n), key=lambda x: p[x]))

    def conj(p, pinv, T):
        return tuple(p[T[pinv[x]]] for x in range(n))

    def frob(T):
        return tuple(T[T[x]] for x in range(n))

    def act(g, frobpow, LR):
        L, R = LR
        for _ in range(frobpow):
            L, R = frob(L), frob(R)
        p, pinv = g
        return (conj(p, pinv, L), conj(p, pinv, R))

    Ginv = [(p, inv(p)) for p in G]
    idx = {pr: i for i, pr in enumerate(pairs)}
    parent = list(range(len(pairs)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    # orbits under generators of Aut(J,c) and Frobenius
    ident = (tuple(range(n)), tuple(range(n)))
    gens_inv = [(p, inv(p)) for p in gens]
    for i, pr in enumerate(pairs):
        for g in gens_inv + [None]:
            img = act(ident, 1, pr) if g is None else act(g, 0, pr)
            j = idx[img]
            a, b = find(i), find(j)
            if a != b:
                parent[a] = b
    orbits = {}
    for i in range(len(pairs)):
        orbits.setdefault(find(i), []).append(i)
    reps = []
    order_frob = r  # Frobenius has order r on these data
    for o in orbits.values():
        rep = pairs[o[0]]
        stab = sum(1 for g in Ginv for f in range(order_frob) if act(g, f, rep) == rep)
        assert len(o) * stab == len(G) * order_frob, "orbit-stabilizer fails"
        reps.append(rep)
    return reps, len(pairs)


def build_R(J, q, L, Rt):
    """Explicit addition/multiplication tables of R (order q*|J|) for checking."""
    h = hensel_factor(q)
    r = len(h)
    n = 1 << J.m
    A, M = J.add, J.mul
    expJ = max(J.additive_order(x) for x in range(n))

    def times(k, y):
        k %= expJ
        z = 0
        for _ in range(k):
            z = A[z][y]
        return z

    def Lpow(j, y, T):
        for _ in range(j):
            y = T[y]
        return y

    def normalize(poly, x):
        poly = list(poly) + [0] * max(0, r - len(poly))
        for d in range(len(poly) - 1, r - 1, -1):      # reduce mod h
            lead = poly[d]
            if lead:
                for k in range(r):
                    poly[d - r + k] -= lead * h[k]
                poly[d] = 0
        a = 0
        for j in range(r):
            bit, t = poly[j] % 2, poly[j] // 2
            if bit:
                a |= 1 << j
            if t and n > 1:   # 2*t*zeta^j = t * (zeta^j c)
                x = A[x][times(t, Lpow(j, J.c, L))]
        return a, x

    elems = [(a, x) for a in range(q) for x in range(n)]
    code = {e: i for i, e in enumerate(elems)}
    N = len(elems)
    add = [[0] * N for _ in range(N)]
    mul = [[0] * N for _ in range(N)]
    for i, (a, x) in enumerate(elems):
        for j, (b, y) in enumerate(elems):
            pa = [((a >> k) & 1) + ((b >> k) & 1) for k in range(r)]
            add[i][j] = code[normalize(pa, A[x][y] if n > 1 else 0)]
            pm = [0] * (2 * r - 1)
            for k in range(r):
                for l in range(r):
                    pm[k + l] += ((a >> k) & 1) * ((b >> l) & 1)
            if n > 1:
                z = M[x][y]
                for k in bits(a):
                    z = A[z][Lpow(k, y, L)]
                for l in bits(b):
                    z = A[z][Lpow(l, x, Rt)]
            else:
                z = 0
            mul[i][j] = code[normalize(pm, z)]
    one = code[(1, 0)]
    return N, add, mul, one, n


def check_and_invariants(N, add, mul, one, nJ):
    for x in range(N):
        assert mul[one][x] == x and mul[x][one] == x
        for y in range(N):
            for z in range(N):
                assert mul[mul[x][y]][z] == mul[x][mul[y][z]]
                assert mul[x][add[y][z]] == add[mul[x][y]][mul[x][z]]
                assert mul[add[x][y]][z] == add[mul[x][z]][mul[y][z]]
    units = [x for x in range(N) if any(mul[x][y] == one for y in range(N))]
    nonunits = set(range(N)) - set(units)
    # local: non-units form an ideal of size nJ
    assert len(nonunits) == nJ
    for x in nonunits:
        for y in nonunits:
            assert add[x][y] in nonunits
        for y in range(N):
            assert mul[x][y] in nonunits and mul[y][x] in nonunits
    zero = 0
    ch, y = 1, one
    while y != zero:
        y = add[y][one]
        ch += 1
    Jl = sorted(nonunits)
    J2 = {mul[x][y] for x in Jl for y in Jl}
    center = sum(all(mul[x][y] == mul[y][x] for y in range(N)) for x in range(N))
    sq0 = sum(mul[x][x] == zero for x in Jl)
    zp = sum(mul[x][y] == zero for x in Jl for y in Jl)
    # elements of J killed by left/right multiplication with J
    lann = sum(all(mul[x][y] == zero for y in Jl) for x in Jl)
    rann = sum(all(mul[y][x] == zero for y in Jl) for x in Jl)
    comm = center == N
    return (ch, comm, center, len(J2), sq0, zp, lann, rann)


def run(q, m, verify=True):
    rings = marked_rings(m)
    results = []
    for J in rings:
        auts = J.automorphisms()
        reps, npairs = classify(J, q, auts)
        for (L, R) in reps:
            inv = None
            if verify:
                inv = check_and_invariants(*build_R(J, q, L, R))
            results.append((J, L, R, inv))
    return results


if __name__ == "__main__":
    for q, m in [(4, 0), (4, 2), (8, 0), (8, 3), (4, 4)]:
        res = run(q, m)
        invs = [x[3] for x in res]
        bychar = {}
        for inv in invs:
            key = (inv[0], inv[1])
            bychar[key] = bychar.get(key, 0) + 1
        print(f"residue F{q}, order {q * 2 ** m}: {len(res)} local rings; "
              f"{len(set(invs))} distinct invariant tuples")
        for (ch, cm), v in sorted(bychar.items()):
            print(f"    char {ch}, {'commutative' if cm else 'noncommutative'}: {v}")
        if q == 4 and m == 4:
            with open("residue_F4_order64.json", "w") as f:
                json.dump([{"J": J.to_json(), "L": list(L), "R": list(R), "inv": list(inv)}
                           for J, L, R, inv in res], f)
        if q == 8 and m == 3:
            with open("residue_F8_order64.json", "w") as f:
                json.dump([{"J": J.to_json(), "L": list(L), "R": list(R), "inv": list(inv)}
                           for J, L, R, inv in res], f)
