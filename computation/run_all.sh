#!/usr/bin/env bash
# Reproduce every count and check.  Run from this directory.
# Needs python3, gcc; GAP with the ModIsom package for step 2
# (e.g. `pip install passagemath-gap`, which provides a `gap` executable).
set -euo pipefail
GAP=${GAP:-gap}

echo "== 1. nilpotent F2-algebras (char-2 local rings, residue F2)"
python3 enum_nilpotent.py 5

echo "== 2. ModIsom canonical forms and |Aut| (no duplicates)"
python3 make_algs_g.py
"$GAP" -q check_modisom.g

echo "== 3. mass-formula check for the F2-algebras"
gcc -O2 -o count_tables count_tables.c
python3 mass_check.py

echo "== 4. all local rings with residue field F2 (every characteristic)"
python3 enum_local_F2.py 5

echo "== 5. axioms, pairwise non-isomorphism, mass formula"
gcc -O2 -o count_marked count_marked.c
gcc -O2 -o ringcheck ringcheck.c
for m in 1 2 3 4 5; do echo "$m $(./count_marked $m)"; done > T_marked.txt
python3 mass_check_local.py 5

echo "== 6. local rings with residue field F4 and F8"
python3 enum_residue_Fq.py
python3 iso_residue.py
