#!/usr/bin/env python3
"""
Independent checks of the list of marked rings (J, c)  <->  local rings with
residue field F_2:
  1. every listed (J,c) satisfies the ring axioms and the marking on all elements;
  2. no two listed rings are isomorphic (invariants, then an explicit
     isomorphism search for every pair with equal invariants);
  3. mass formula  sum f(J)/|Aut(J,c)| = T(m)/2^(m(m-1)/2), with T(m) from
     count_marked.c, so no class is missing.
"""
import json
import subprocess
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
data = json.load(open("local_F2_upto5.json"))
T = {}
for line in open("T_marked.txt"):
    m, t = line.split()
    T[int(m)] = int(t)

allok = True
for m in range(1, MMAX + 1):
    rings = data[str(m)]
    with open(f"rings_{m}.txt", "w") as f:
        f.write(f"{len(rings)}\n")
        for r in rings:
            f.write(f"{r['m']} {r['c']}\n{' '.join(map(str, r['P']))}\n")
            f.write(" ".join(str(x) for row in r["Q"] for x in row) + "\n")
    out = subprocess.check_output(["./ringcheck", "stats", f"rings_{m}.txt"]).decode().split("\n")
    stats = [list(map(int, l.split())) for l in out if l.strip()]
    assert len(stats) == len(rings)
    axioms_ok = all(s[0] == 1 for s in stats)
    mass = sum(Fraction(s[2], s[1]) for s in stats)
    target = Fraction(T[m], 2 ** (m * (m - 1) // 2))
    # duplicates: group by invariants (|Aut|, f, and structural counts)
    groups = defaultdict(list)
    for idx, s in enumerate(stats):
        groups[tuple(s[1:])].append(idx)
    pairs = [p for g in groups.values() for p in combinations(g, 2)]
    iso = 0
    if pairs:
        res = subprocess.run(["./ringcheck", "iso", f"rings_{m}.txt"],
                             input="\n".join(f"{a} {b}" for a, b in pairs),
                             capture_output=True, text=True).stdout.split("\n")
        iso = sum(int(l.split()[2]) for l in res if l.strip())
    ok = axioms_ok and mass == target and iso == 0
    allok &= ok
    print(f"order {2 ** (m + 1)}: {len(rings)} rings | axioms {'ok' if axioms_ok else 'FAIL'} | "
          f"{len(pairs)} same-invariant pairs, {iso} isomorphic | "
          f"mass {mass} vs T/2^{m*(m-1)//2} = {target} -> {'OK' if ok else 'MISMATCH'}")
print("ALL OK" if allok else "FAILED")
