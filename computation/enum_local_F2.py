#!/usr/bin/env python3
r"""
Local rings R of order 2^(m+1) with residue field F_2, all characteristics.

Reduction.  R is determined by the pair (J, c) with J = J(R) (a nilpotent ring of
order 2^m) and c = 2*1_R in J, which satisfies c*x = x*c = 2x for all x in J.
Conversely every such "marked" pair gives R = (Z + J)/<(2,-c)>, and
R ~= R'  iff  (J,c) ~= (J',c').  char R = 2 * (additive order of c).

Enumeration (central extensions by the full annihilator).
  For a marked pair, Z = Ann(J) is elementary abelian (2z = cz = 0), and
  (K, cbar) = (J/Z, c mod Z) is again a marked pair.  Fix a presentation of K
  and lifts g_i of its generators.  J is then determined by data in Z:
      2 g_i     = P_i  + tau_i        (additive extension)
      g_i g_j   = Q_ij + theta_ij     (multiplicative cocycle)
      c         = cbar + zc           (position of the marking)
  subject to linear conditions (well-defined, associative, c*x = x*c = 2x).
  Changing the lifts g_i -> g_i + h_i gives the coboundaries.  With
  H = cocycles / coboundaries and s = dim Z:
      iso classes of (J,c) with J/Ann(J) ~= (K,cbar), dim Ann(J) = s
        <->  Aut(K,cbar)-orbits of subspaces W of H with dim W <= s
             and Ann(K) /\ Ann(W) = 0.
  (A tuple of s classes up to GL(s) is the same as its span.)

Rings are stored as "pc presentations": generators g_0..g_{m-1}, every element
uniquely sum of a subset (bitmask), 2 g_i = P[i] (bits > i), g_i g_j = Q[i][j]
(bits > max(i,j) -- the construction keeps this triangular), marking c.
"""
import json
import sys
from itertools import combinations

sys.setrecursionlimit(10000)
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def bits(v):
    i = 0
    while v:
        if v & 1:
            yield i
        v >>= 1
        i += 1


def parity(v):
    return bin(v).count("1") & 1


# ---------------- F_2 linear algebra on bitmasks ----------------

def rref(vecs):
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
    piv = {}
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
    out = []
    for f in range(nvars):
        if f in piv:
            continue
        x = 1 << f
        for p, r in piv.items():
            if (r >> f) & 1:
                x |= 1 << p
        out.append(x)
    return out


class Coords:
    """Coordinates w.r.t. a basis of a subspace."""

    def __init__(self, basis):
        rows = []
        for idx, v in enumerate(basis):
            tag = 1 << idx
            for b, t in rows:
                if (v >> (b.bit_length() - 1)) & 1:
                    v ^= b
                    tag ^= t
            assert v
            p = v.bit_length() - 1
            rows = [(b ^ v, t ^ tag) if (b >> p) & 1 else (b, t) for b, t in rows]
            rows.append((v, tag))
        self.rows = rows

    def coords(self, v):
        tag = 0
        for b, t in self.rows:
            if (v >> (b.bit_length() - 1)) & 1:
                v ^= b
                tag ^= t
        assert v == 0, "not in span"
        return tag


def subspaces(h, s):
    res = []
    for pivs in combinations(range(h - 1, -1, -1), s):
        pivset = set(pivs)
        frees = [[q for q in range(p) if q not in pivset] for p in pivs]

        def rec(r, acc):
            if r == s:
                res.append(tuple(acc))
                return
            fr = frees[r]
            for msk in range(1 << len(fr)):
                v = 1 << pivs[r]
                for t, q in enumerate(fr):
                    if (msk >> t) & 1:
                        v |= 1 << q
                rec(r + 1, acc + [v])

        rec(0, [])
    return res


# ---------------- rings in pc presentation ----------------

class Ring:
    def __init__(self, m, P, Q, c):
        self.m, self.P, self.Q, self.c = m, list(P), [list(r) for r in Q], c
        n = 1 << m
        # addition with carry record: carry bit l flips each time g_l + g_l occurs
        self.add = [[0] * n for _ in range(n)]
        self.addcarry = [[0] * n for _ in range(n)]
        for a in range(n):
            for b in range(n):
                self.add[a][b], self.addcarry[a][b] = self._add(a, b)
        # multiplication: K-part, carries from summing the Q's, theta symbols used
        self.mul = [[0] * n for _ in range(n)]
        self.mulcarry = [[0] * n for _ in range(n)]
        self.multheta = [[0] * n for _ in range(n)]  # bitmask over pairs i*m+j
        for a in range(n):
            for b in range(n):
                acc, car, th = 0, 0, 0
                for i in bits(a):
                    for j in bits(b):
                        car ^= self.addcarry[acc][self.Q[i][j]]
                        acc = self.add[acc][self.Q[i][j]]
                        th ^= 1 << (i * m + j)
                self.mul[a][b], self.mulcarry[a][b], self.multheta[a][b] = acc, car, th

    def _add(self, x, y):
        carry = 0
        for l in bits(y):
            x, carry = self._add_gen(x, l, carry)
        return x, carry

    def _add_gen(self, x, l, carry):
        if not (x >> l) & 1:
            return x | (1 << l), carry
        x &= ~(1 << l)
        carry ^= 1 << l
        for q in bits(self.P[l]):
            x, carry = self._add_gen(x, q, carry)
        return x, carry

    def to_json(self):
        return {"m": self.m, "P": self.P, "Q": self.Q, "c": self.c}

    # ---- invariants / checks ----
    def is_ring(self):
        n = 1 << self.m
        A, M = self.add, self.mul
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if M[A[x][y]][z] != A[M[x][z]][M[y][z]]:
                        return False
                    if M[x][A[y][z]] != A[M[x][y]][M[x][z]]:
                        return False
                    if M[M[x][y]][z] != M[x][M[y][z]]:
                        return False
            if M[self.c][x] != A[x][x] or M[x][self.c] != A[x][x]:
                return False
        return True

    def is_commutative(self):
        n = 1 << self.m
        return all(self.mul[x][y] == self.mul[y][x] for x in range(n) for y in range(n))

    def additive_order(self, x):
        k, y = 1, x
        while y:
            y = self.add[y][x]
            k += 1
        return k if x else 1

    def characteristic(self):
        """char of the local ring R built from (J, c): 2 * order(c)."""
        return 2 * self.additive_order(self.c)

    def annihilator(self):
        n = 1 << self.m
        return [x for x in range(1, n)
                if all(self.mul[x][1 << j] == 0 and self.mul[1 << j][x] == 0
                       for j in range(self.m))]

    def apply(self, phi, x):
        r = 0
        for q in bits(x):
            r = self.add[r][phi[q]]
        return r

    def automorphisms(self):
        """Automorphisms fixing c, as tuples of generator images (backtracking
        from the deepest generator; every relation is checked once all the
        generators it involves have images)."""
        m = self.m
        if m == 0:
            return [()]
        A, M, P, Q = self.add, self.mul, self.P, self.Q
        n = 1 << m
        imgs = [0] * m
        out = []

        def rec(i, S):
            if i < 0:
                if self.apply(imgs, self.c) == self.c:
                    out.append(tuple(imgs))
                return
            Sset = set(S)
            for h in range(n):
                if h in Sset:
                    continue
                if A[h][h] != self.apply(imgs, P[i]):
                    continue
                imgs[i] = h
                ok = True
                for j in range(i, m):
                    hj = imgs[j]
                    if M[h][hj] != self.apply(imgs, Q[i][j]) or \
                       M[hj][h] != self.apply(imgs, Q[j][i]):
                        ok = False
                        break
                if ok:
                    rec(i - 1, S + [A[h][x] for x in S])
            imgs[i] = 0

        rec(m - 1, [0])
        return out


def compose(R, phi, psi):
    return tuple(R.apply(phi, x) for x in psi)


def generators(R, group):
    ident = tuple(1 << i for i in range(R.m))
    gens, H = [], {ident}
    for g in group:
        if g in H:
            continue
        gens.append(g)
        queue = list(H)
        while queue:
            x = queue.pop()
            for s in gens:
                y = compose(R, s, x)
                if y not in H:
                    H.add(y)
                    queue.append(y)
        if len(H) == len(group):
            break
    assert len(H) == len(group)
    return gens


# ---------------- extensions ----------------

def extensions(K, s, auts):
    """All marked (J,c) with Ann(J) of dim s and (J/Ann J, cbar) ~= K."""
    m = K.m
    TH = lambda i, j: 1 << (i * m + j)
    TAU = lambda i: 1 << (m * m + i)
    ZC = 1 << (m * m + m)
    N = m * m + m + 1

    def tau_of(carry):
        z = 0
        for l in bits(carry):
            z ^= TAU(l)
        return z

    # symbolic J-elements: (kpart, zform) ; zform = linear form in the data
    def sadd(x, y):
        return (K.add[x[0]][y[0]], x[1] ^ y[1] ^ tau_of(K.addcarry[x[0]][y[0]]))

    def smul(x, y):
        a, b = x[0], y[0]
        th = 0
        for pair in bits(K.multheta[a][b]):
            th ^= 1 << pair  # pair index i*m+j is exactly the theta symbol
        return (K.mul[a][b], th ^ tau_of(K.mulcarry[a][b]))

    g = [(1 << i, 0) for i in range(m)]
    two_g = [(K.P[i], TAU(i)) for i in range(m)]
    c = (K.c, ZC)

    eqs = []

    def eq(x, y):
        assert x[0] == y[0], "K is not a valid marked ring"
        if x[1] ^ y[1]:
            eqs.append(x[1] ^ y[1])

    for i in range(m):
        for j in range(m):
            gij = smul(g[i], g[j])
            eq(smul(two_g[i], g[j]), sadd(gij, gij))
            gji = smul(g[j], g[i])
            eq(smul(g[j], two_g[i]), sadd(gji, gji))
            for k in range(m):
                eq(smul(smul(g[i], g[j]), g[k]), smul(g[i], smul(g[j], g[k])))
        tw = sadd(g[i], g[i])
        eq(smul(c, g[i]), tw)
        eq(smul(g[i], c), tw)
    Zc = nullspace(eqs, N)

    # coboundaries: g_q -> g_q + h_q
    B = []
    for q in range(m):
        v = 0
        for i in range(m):
            if (K.P[i] >> q) & 1:
                v ^= TAU(i)
            for j in range(m):
                if (K.Q[i][j] >> q) & 1:
                    v ^= TH(i, j)
        if (K.c >> q) & 1:
            v ^= ZC
        B.append(v)
    B = list(rref(B))
    for b in B:
        assert all(parity(e & b) == 0 for e in eqs), "coboundary is not a cocycle"
    comp, span = [], list(B)
    for z in Zc:
        if len(rref(span + [z])) > len(rref(span)):
            comp.append(z)
            span.append(z)
    h = len(comp)
    co = Coords(B + comp)
    nb = len(B)
    to_H = lambda v: co.coords(v) >> nb

    # action of an automorphism phi of (K,cbar): new lifts u_i = phi(g_i)
    def transform_forms(phi):
        u = [(phi[i], 0) for i in range(m)]

        def ssum(mask):
            acc = (0, 0)
            for q in bits(mask):
                acc = sadd(acc, u[q])
            return acc

        forms = {}
        for i in range(m):
            two = sadd(u[i], u[i])
            rhs = ssum(K.P[i])
            assert two[0] == rhs[0]
            forms[TAU(i)] = two[1] ^ rhs[1]
            for j in range(m):
                pr = smul(u[i], u[j])
                rhs = ssum(K.Q[i][j])
                assert pr[0] == rhs[0]
                forms[TH(i, j)] = pr[1] ^ rhs[1]
        rhs = ssum(K.c)
        assert rhs[0] == K.c
        forms[ZC] = c[1] ^ rhs[1]
        return forms

    def transform(forms, d):
        out = 0
        for sym, f in forms.items():
            if parity(f & d):
                out |= sym
        return out

    mats = []
    for phi in generators(K, auts):
        forms = transform_forms(phi)
        for b in B:
            assert to_H(transform(forms, b)) == 0
        mats.append([to_H(transform(forms, cv)) for cv in comp])

    lift = lambda w: 0 if w == 0 else _lift(w, comp)

    # Ann condition: every nonzero x in Ann(K) must be "seen" by some w in W
    annforms = []
    for x in K.annihilator():
        fs = [smul((x, 0), g[j])[1] for j in range(m)] + \
             [smul(g[j], (x, 0))[1] for j in range(m)]
        annforms.append(fs)

    def ann_ok(W):
        lifts = [lift(w) for w in W]
        for fs in annforms:
            if not any(parity(f & L) for f in fs for L in lifts):
                return False
        return True

    results = []
    for t in range(0, min(s, h) + 1):
        pts = [W for W in subspaces(h, t) if ann_ok(W)]
        if not pts:
            continue
        index = {W: i for i, W in enumerate(pts)}
        parent = list(range(len(pts)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for W, i in index.items():
            for Mx in mats:
                img = rref([_apply_mat(Mx, w) for w in W])
                j = index[img]
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
        seen = set()
        for i, W in enumerate(pts):
            r = find(i)
            if r in seen:
                continue
            seen.add(r)
            results.append(build(K, s, [lift(w) for w in W]))
    return results


def _lift(w, comp):
    th = 0
    for t in bits(w):
        th ^= comp[t]
    return th


def _apply_mat(Mx, v):
    r = 0
    for i in bits(v):
        r ^= Mx[i]
    return r


def build(K, s, data):
    """J = K extended by Z = F_2^s with cocycle vectors data (len <= s)."""
    m = K.m
    n = m + s
    P = list(K.P) + [0] * s
    Q = [list(K.Q[i]) + [0] * s for i in range(m)] + [[0] * n for _ in range(s)]
    c = K.c
    for t, d in enumerate(data):
        zb = 1 << (m + t)
        for i in range(m):
            if (d >> (m * m + i)) & 1:
                P[i] |= zb
            for j in range(m):
                if (d >> (i * m + j)) & 1:
                    Q[i][j] |= zb
        if (d >> (m * m + m)) & 1:
            c |= zb
    return Ring(n, P, Q, c)


def main():
    levels = {0: [Ring(0, [], [], 0)]}
    summary = {}
    for mm in range(1, MMAX + 1):
        found = []
        for s in range(1, mm + 1):
            for K in levels[mm - s]:
                auts = K.automorphisms()
                found += extensions(K, s, auts)
        # sanity: every constructed ring is a valid marked nilpotent ring
        if mm <= 3:
            assert all(J.is_ring() for J in found)
        levels[mm] = found
        bychar = {}
        for J in found:
            ch = J.characteristic()
            key = (ch, J.is_commutative())
            bychar[key] = bychar.get(key, 0) + 1
        comm = sum(v for (ch, cm), v in bychar.items() if cm)
        print(f"|J| = 2^{mm} (rings of order {2 ** (mm + 1)}): {len(found)} local rings, "
              f"{comm} commutative", flush=True)
        for ch in sorted({k[0] for k in bychar}):
            a = bychar.get((ch, True), 0)
            b = bychar.get((ch, False), 0)
            print(f"    char {ch:3d}: {a + b:4d}  (commutative {a}, noncommutative {b})")
        summary[mm] = {str(k): v for k, v in bychar.items()}
    with open(f"local_F2_upto{MMAX}.json", "w") as f:
        json.dump({mm: [J.to_json() for J in levels[mm]] for mm in levels if mm > 0}, f)


if __name__ == "__main__":
    main()
