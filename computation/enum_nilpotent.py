#!/usr/bin/env python3
r"""
Enumerate nilpotent associative F_2-algebras of dimension <= NMAX up to isomorphism.

Why: a local ring R of characteristic 2 with R/J(R) = F_2 is F_2*1 + J with J a
nilpotent F_2-algebra, and R <-> J is a bijection on isomorphism classes.
So local F_2-algebras of order 2^(n+1) with residue field F_2  <->  nilpotent
F_2-algebras of dimension n.  For |R| = 64 we need n = 5.

Method (Skjelbred-Sund / de Graaf central extensions):
  * Every nilpotent A is uniquely A0 (+) F^k with A0 having no annihilator
    component.
  * A0 of dimension n with dim Ann(A0) = s >= 1 is B_theta, where B = A0/Ann(A0)
    (dim n-s), theta = (theta_1..theta_s) in Z^2(B,F) with the theta_i linearly
    independent in H^2(B,F), and Ann(B) /\ Ann(theta) = 0.
  * B_theta ~= B_eta  iff  <theta> and <eta> lie in the same Aut(B)-orbit of
    s-dimensional subspaces of H^2(B,F).

Vectors over F_2 are Python ints (bit k = coefficient of basis vector e_k).
An algebra of dimension n is a tuple-of-tuples mult[i][j] = e_i * e_j (bitmask).
"""
import json
import sys
from itertools import combinations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


# ---------- linear algebra over F_2 with bitmask vectors ----------

def rref(vecs):
    """Reduced echelon basis (pivot = highest bit); canonical for the span."""
    basis = []
    for v in vecs:
        for b in basis:
            if (v >> (b.bit_length() - 1)) & 1:
                v ^= b
        if v:
            p = v.bit_length() - 1
            basis = [b ^ v if (b >> p) & 1 else b for b in basis]
            basis.append(v)
    return tuple(sorted(basis, reverse=True))


def nullspace(rows, nvars):
    """Basis of {x : <row, x> = 0 for all rows} in F_2^nvars."""
    piv = {}  # pivot bit -> row
    for r in rows:
        for p, b in piv.items():
            if (r >> p) & 1:
                r ^= b
        if r:
            p = r.bit_length() - 1
            for q in list(piv):
                if (piv[q] >> p) & 1:
                    piv[q] ^= r
            piv[p] = r
    free = [i for i in range(nvars) if i not in piv]
    basis = []
    for f in free:
        x = 1 << f
        for p, r in piv.items():
            if (r >> f) & 1:
                x |= 1 << p
        basis.append(x)
    return basis


def bits(v):
    i = 0
    while v:
        if v & 1:
            yield i
        v >>= 1
        i += 1


def apply(M, v):
    """M = list of column images; returns M*v."""
    r = 0
    for i in bits(v):
        r ^= M[i]
    return r


class Coords:
    """Express vectors in a basis of a subspace (raises if not in the span)."""

    def __init__(self, basis):
        self.rows = []  # (vec, tag) echelon with tags recording combinations
        for idx, v in enumerate(basis):
            tag = 1 << idx
            for b, t in self.rows:
                if (v >> (b.bit_length() - 1)) & 1:
                    v ^= b
                    tag ^= t
            assert v, "basis not independent"
            self.rows.append((v, tag))
            self.rows.sort(key=lambda bt: -bt[0])
        # make fully reduced so order does not matter
        changed = True
        while changed:
            changed = False
            for i, (b, t) in enumerate(self.rows):
                p = b.bit_length() - 1
                for j, (c, u) in enumerate(self.rows):
                    if j != i and (c >> p) & 1:
                        self.rows[j] = (c ^ b, u ^ t)
                        changed = True

    def coords(self, v):
        tag = 0
        for b, t in self.rows:
            if (v >> (b.bit_length() - 1)) & 1:
                v ^= b
                tag ^= t
        assert v == 0, "vector not in span"
        return tag


# ---------- algebras ----------

def prod(mult, x, y):
    r = 0
    for i in bits(x):
        row = mult[i]
        for j in bits(y):
            r ^= row[j]
    return r


def is_associative(mult):
    n = len(mult)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if prod(mult, mult[i][j], 1 << k) != prod(mult, 1 << i, mult[j][k]):
                    return False
    return True


def annihilator(mult):
    n = len(mult)
    # x in Ann iff x*e_j = e_j*x = 0 for all j; linear conditions on x
    rows = []
    for j in range(n):
        for k in range(n):
            # coefficient of e_k in x*e_j and e_j*x as linear forms in x
            r1 = 0
            r2 = 0
            for i in range(n):
                if (mult[i][j] >> k) & 1:
                    r1 |= 1 << i
                if (mult[j][i] >> k) & 1:
                    r2 |= 1 << i
            rows += [r1, r2]
    return nullspace(rows, n)


def invertible_matrices(d):
    """All d x d invertible matrices over F_2, as tuples of column images."""
    out = []

    def rec(cols, span):
        if len(cols) == d:
            out.append(tuple(cols))
            return
        for v in range(1, 1 << d):
            if v not in span:
                rec(cols + [v], span | {s ^ v for s in span})

    rec([], {0})
    return out


_GL = {}


def automorphisms(mult):
    d = len(mult)
    if d not in _GL:
        _GL[d] = invertible_matrices(d)
    auts = []
    for phi in _GL[d]:
        ok = True
        for i in range(d):
            for j in range(d):
                if apply(phi, mult[i][j]) != prod(mult, phi[i], phi[j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(phi)
    return auts


def compose(phi, psi):
    return tuple(apply(phi, c) for c in psi)


def generators(group):
    """A small generating set of a finite group given as a list."""
    ident = tuple(1 << i for i in range(len(group[0])))
    gens = []
    H = {ident}
    for g in group:
        if g in H:
            continue
        gens.append(g)
        queue = list(H)  # re-close under the enlarged generating set
        while queue:
            x = queue.pop()
            for s in gens:
                y = compose(s, x)
                if y not in H:
                    H.add(y)
                    queue.append(y)
        if len(H) == len(group):
            break
    assert len(H) == len(group)
    return gens


def subspaces(h, s):
    """All s-dim subspaces of F_2^h in canonical RREF (tuple of ints)."""
    res = []
    for pivs in combinations(range(h - 1, -1, -1), s):  # descending pivots
        pivset = set(pivs)
        # row r: pivot pivs[r], free bits at positions < pivs[r] not pivots
        frees = [[q for q in range(p) if q not in pivset] for p in pivs]

        def rec(r, acc):
            if r == s:
                res.append(tuple(acc))
                return
            fr = frees[r]
            for m in range(1 << len(fr)):
                v = 1 << pivs[r]
                for t, q in enumerate(fr):
                    if (m >> t) & 1:
                        v |= 1 << q
                rec(r + 1, acc + [v])

        rec(0, [])
    return res


def central_extensions(B, s, auts):
    """Algebras without annihilator component with A/Ann(A) = B, dim Ann = s."""
    d = len(B)
    nv = d * d
    var = lambda i, j: 1 << (i * d + j)
    # cocycle equations theta(e_i e_j, e_k) = theta(e_i, e_j e_k)
    eqs = []
    for i in range(d):
        for j in range(d):
            for k in range(d):
                e = 0
                for m in bits(B[i][j]):
                    e ^= var(m, k)
                for m in bits(B[j][k]):
                    e ^= var(i, m)
                if e:
                    eqs.append(e)
    Z2 = nullspace(eqs, nv)
    # coboundaries delta f (x,y) = f(xy)
    B2 = []
    for k in range(d):
        v = 0
        for i in range(d):
            for j in range(d):
                if (B[i][j] >> k) & 1:
                    v |= var(i, j)
        B2.append(v)
    B2 = list(rref(B2))
    # complement of B2 in Z2
    comp = []
    span = list(B2)
    for z in Z2:
        if len(rref(span + [z])) > len(rref(span)):
            comp.append(z)
            span.append(z)
    h = len(comp)
    if h < s:
        return []
    co = Coords(B2 + comp)
    nb = len(B2)

    def h2(theta):
        return co.coords(theta) >> nb

    def act(phi, theta):
        # (phi.theta)(e_i, e_j) = theta(phi e_i, phi e_j)
        out = 0
        for i in range(d):
            for j in range(d):
                val = 0
                for a in bits(phi[i]):
                    for b in bits(phi[j]):
                        val ^= (theta >> (a * d + b)) & 1
                if val:
                    out |= var(i, j)
        return out

    gens = generators(auts)
    mats = []
    for phi in gens:
        for b in B2:  # sanity: coboundaries are preserved
            assert h2(act(phi, b)) == 0
        mats.append([h2(act(phi, c)) for c in comp])

    def lift(w):
        th = 0
        for t in bits(w):
            th ^= comp[t]
        return th

    annB = annihilator(B)
    annB_vecs = [v for v in rref(annB)]
    # all nonzero elements of Ann(B)
    annB_elems = set()
    for m in range(1, 1 << len(annB_vecs)):
        x = 0
        for t in bits(m):
            x ^= annB_vecs[t]
        annB_elems.add(x)

    def ann_ok(W):
        thetas = [lift(w) for w in W]
        for x in annB_elems:
            killed = True
            for th in thetas:
                for j in range(d):
                    if prod_form(th, x, 1 << j) or prod_form(th, 1 << j, x):
                        killed = False
                        break
                if not killed:
                    break
            if killed:
                return False
        return True

    def prod_form(th, x, y):
        val = 0
        for a in bits(x):
            for b in bits(y):
                val ^= (th >> (a * d + b)) & 1
        return val

    pts = [W for W in subspaces(h, s) if ann_ok(W)]
    index = {W: i for i, W in enumerate(pts)}
    parent = list(range(len(pts)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for W, i in index.items():
        for M in mats:
            img = rref([apply(M, w) for w in W])
            j = index[img]  # KeyError would mean the condition is not invariant
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    reps = {}
    for i, W in enumerate(pts):
        r = find(i)
        if r not in reps:
            reps[r] = W
    out = []
    for W in reps.values():
        thetas = [lift(w) for w in W]
        n = d + s
        mult = [[0] * n for _ in range(n)]
        for i in range(d):
            for j in range(d):
                v = B[i][j]
                for t, th in enumerate(thetas):
                    if (th >> (i * d + j)) & 1:
                        v |= 1 << (d + t)
                mult[i][j] = v
        out.append(tuple(tuple(r) for r in mult))
    return out


def pad(mult, n):
    m = len(mult)
    return tuple(tuple(list(mult[i]) + [0] * (n - m)) if i < m else tuple([0] * n)
                 for i in range(n))


def is_commutative(mult):
    n = len(mult)
    return all(mult[i][j] == mult[j][i] for i in range(n) for j in range(n))


def main():
    noann = {0: [()]}  # algebras without annihilator component, by dimension
    allalg = {0: [()]}
    for n in range(1, NMAX + 1):
        noann[n] = []
        for s in range(1, n + 1):
            for B in allalg[n - s]:
                if len(B) == 0:
                    continue  # H^2 of the zero-dim algebra is 0
                auts = automorphisms(B)
                noann[n] += central_extensions(B, s, auts)
        allalg[n] = [pad(A0, n) for m in range(n + 1) for A0 in noann[m]]
        for A in allalg[n]:
            assert is_associative(A)
        ncomm = sum(is_commutative(A) for A in allalg[n])
        print(f"dim {n}: no-ann-component {len(noann[n])}, total {len(allalg[n])}, "
              f"commutative {ncomm}", flush=True)
    with open(f"nilpotent_F2_upto{NMAX}.json", "w") as f:
        json.dump({n: [[list(r) for r in A] for A in allalg[n]] for n in allalg if n > 0}, f)


if __name__ == "__main__":
    main()
