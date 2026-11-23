# Unital rings of order 64

Goal: count the unital (associative) rings with 64 elements up to isomorphism.
The sequence is [OEIS A037291](https://oeis.org/A037291); its known terms
stop at n = 63 (a(32) = 208).

## Current status

| | count |
|---|---:|
| Local rings of order 64 | **1,789** |
| &nbsp;&nbsp;residue field F2 | 1,773 |
| &nbsp;&nbsp;residue field F4 | 11 |
| &nbsp;&nbsp;residue field F8 | 4 |
| &nbsp;&nbsp;residue field F64 | 1 |
| Commutative local rings of order 64 | 259 |
| **Commutative unital rings of order 64** | **462** |
| **All unital rings of order 64** | **2,134 + N** |

N is the number of indecomposable **non-local** rings of order 64 (R/J a
product of at least two fields). That is the one piece still missing, so
currently **a(64) ≥ 2,134**.

The commutative count is complete, because every commutative indecomposable
finite ring is local. (OEIS [A127707](https://oeis.org/A127707) also stops at
n = 63.)

## Reduction

1. **Indecomposables.** A finite unital ring splits uniquely (via central
   idempotents) into indecomposable rings, so the totals are the Euler
   transform of the indecomposable counts. From a(2^k) = 1, 4, 11, 50, 208
   (k = 1..5) the indecomposable counts are i_k = 1, 3, 7, 33, 137, and

   **a(64) = 345 + i_6.**

   Commutative: indecomposable counts 1, 3, 6, 21, 54 (from
   [A127707](https://oeis.org/A127707)), giving 203 + (commutative local rings
   of order 64).
2. **No matrix blocks.** If R/J(R) had a block M_2(F_2), lifting matrix units
   makes that corner M_2(S), of order |S|^4 (16 or 256, never 64), and gluing
   anything to it needs at least 128 elements. So an indecomposable ring of
   order 64 is local, or basic non-local (R/J a product of ≥ 2 fields).
3. **Local rings, residue field F2** ↔ *marked pairs* (J, c): J = J(R) is a
   nilpotent ring of order |R|/2 and c = 2·1 ∈ J satisfies cx = xc = 2x.
   Conversely R = (Z ⊕ J)/⟨(2, −c)⟩, and R ≅ R' iff (J, c) ≅ (J', c').
   char R = 2·ord(c); characteristic 2 is c = 0, i.e. J a nilpotent F2-algebra.
4. **Local rings, residue field F_q (q = 4, 8)** ↔ a marked pair (J, c) of
   order |R|/q together with the left/right action L, Rt of a Teichmüller unit
   ζ (order q − 1) on J: h(L) = h(Rt) = 0 for h the 2-adic lift of the minimal
   polynomial, L·Rt = Rt·L, L(c) = Rt(c), L(xy) = L(x)y, Rt(xy) = xRt(y),
   Rt(x)y = xL(y). Isomorphism classes are orbits of Aut(J, c) × ⟨ζ ↦ ζ²⟩.

## Results

### Local rings with residue field F2

| order | rings | commutative | char 2 | char 4 | char 8 | char 16 | char 32 | char 64 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 2 | 2 | 1 | 1 | | | | |
| 8 | 5 | 5 | 2 | 2 | 1 | | | |
| 16 | 22 | 18 | 7 | 12 | 2 | 1 | | |
| 32 | 105 | 53 | 27 | 61 | 14 | 2 | 1 | |
| **64** | **1773** | **251** | 404 | 1206 | 146 | 14 | 2 | 1 |

Order 64, commutative / noncommutative by characteristic: char 2: 38 / 366;
char 4: 137 / 1069; char 8: 61 / 85; char 16: 12 / 2; char 32: 2 / 0;
char 64: 1 / 0.

### Local rings with residue field F4, F8, F64

| residue | order | rings | char 2 (comm / noncomm) | char 4 (comm / noncomm) | char 8 |
|---|---:|---:|---|---|---:|
| F4 | 16 | 3 | 1 / 1 | 1 / 0 | |
| F4 | **64** | **11** | 2 / 3 | 2 / 3 | 1 (GR(8,2)) |
| F8 | **64** | **4** | 1 / 2 | 1 / 0 (GR(4,3)) | |
| F64 | 64 | 1 | 1 / 0 | | |

## Methods

* **Nilpotent F2-algebras** (`enum_nilpotent.py`): Skjelbred–Sund central
  extensions — algebras without annihilator component are B_θ for
  B = A/Ann(A), with θ running over Aut(B)-orbits of subspaces of H²(B, F2)
  satisfying Ann(B) ∩ Ann(θ) = 0.
* **Marked rings, all characteristics** (`enum_local_F2.py`): central
  extensions by the full annihilator Z = Ann(J) (elementary abelian, since
  2z = cz = 0). For (K, c̄) = (J/Z, c mod Z) in a fixed polycyclic
  presentation, J is given by additive data τ (2g_i = P_i + τ_i), product
  data θ (g_i g_j = Q_ij + θ_ij) and the offset z_c of c, subject to linear
  conditions; changing lifts gives the coboundaries. Classes ↔
  Aut(K, c̄)-orbits of subspaces W ≤ H with dim W ≤ dim Z and
  Ann(K) ∩ Ann(W) = 0.
* **Residue fields F4, F8** (`enum_residue_Fq.py`): enumerate all additive
  endomorphisms of J satisfying the conditions in step 4 above, then orbits.

## Verification

* **Known values reproduced.** Commutative local rings with residue F2:
  2, 5, 18, 53 for orders 4–32. Characteristic-2 subset: 1, 2, 7, 27, 404.
  Order 16 is consistent with the 13 noncommutative unital rings of order 16
  ([A127708](https://oeis.org/A127708)): 4 local with residue F2, 7 non-local
  (including M_2(F2)), F4[x;σ]/(x²), and T_2(F2) × F2. Residue F4 at order 16
  gives the expected 3; residue F8 at order 64 gives the expected 4.
* **No duplicates.** The 404 F2-algebras have pairwise distinct canonical
  forms under GAP/ModIsom. For the 1,773 marked rings, all 11,508 pairs with
  equal invariants were ruled out by exhaustive isomorphism search (which
  finds 60/60 randomly re-presented copies as a control). Same for the
  residue-F4/F8 rings.
* **Nothing missing (mass formula).**
  Σ f(J)/|Aut(J, c)| = T(m)/2^(m(m−1)/2) holds exactly for m = 1..5, where
  T(m) counts triangular marked presentations (2, 8, 112, 5824, 2292096;
  the C count matches brute force for m ≤ 3) and f(J) counts adapted chains
  of ideals. The analogous identity for F2-algebras holds for n = 1..5.
* **Axioms.** Every ring is checked on all elements (residue-F4/F8 rings are
  built as explicit 64-element tables and checked to be local with the right
  radical). Every orbit satisfies orbit–stabilizer against the full group.

## Reproducing

```
cd computation
./run_all.sh          # set GAP=/path/to/gap if needed
```

Needs Python 3 and gcc; step 2 needs GAP with ModIsom
(`pip install passagemath-gap` provides one). Total running time is a few
minutes.

## Data files (`computation/`)

* `nilpotent_F2_upto5.json` — `{n: [table, ...]}`; `table[i][j]` is the
  bitmask of e_i·e_j in the basis e_0..e_{n−1} of a nilpotent F2-algebra.
  Dimension 5 = the 404 characteristic-2 local rings of order 64.
* `local_F2_upto5.json` — `{m: [{m, P, Q, c}, ...]}`: marked rings in
  polycyclic presentation. Generators g_0..g_{m−1}; each element is a subset
  sum (bitmask); 2g_i = `P[i]`, g_i g_j = `Q[i][j]`, c = 2·1. The local ring is
  R = (Z ⊕ J)/⟨(2, −c)⟩ of order 2^(m+1). Key `"5"` = the 1,773 rings of order 64.
* `residue_F4_order64.json`, `residue_F8_order64.json` — each entry has the
  marked ring `J`, the value tables `L`, `R` of y ↦ ζy and y ↦ yζ on J, and
  invariants (char, commutative, |centre|, |J²|, #{x ∈ J: x² = 0},
  #{(x, y) ∈ J²: xy = 0}, |left ann of J in J|, |right ann|).
* `autsizes_n.txt`, `T_marked.txt` — automorphism group orders (ModIsom) and
  presentation counts used by the mass-formula checks.

## Remaining work

* Non-local indecomposable rings of order 64 (Peirce decomposition /
  quivers with Galois-ring coefficients). Validation targets: 7 at order 16
  and 31 at order 32 (from i = 33, 137).
* Optional cross-check of the local counts against Akers & Szabo,
  *A partial classification of local rings of order p⁶*, Involve 16 (2023)
  151–165.

## References

* OEIS A037291, A127707, A127708.
* B. Corbas, G. D. Williams, Rings of order p⁵, Parts I and II,
  J. Algebra 231 (2000).
* C. Akers, S. Szabo, A partial classification of local rings of order p⁶,
  Involve 16 (2023) 151–165.
* B. Eick, ModIsom (GAP package).
* W. A. de Graaf, Classification of nilpotent associative algebras of small
  dimension, IJAC 28 (2018) 133–161.

These are computational results and have not been peer reviewed.

## License

* **Code** (`computation/*.py`, `*.c`, `*.g`, `*.sh`): MIT License — see
  [LICENSE](LICENSE).
* **Data** (`computation/*.json`, `computation/*.txt`): Creative Commons
  Attribution 4.0 International (CC BY 4.0) — see [LICENSE-DATA](LICENSE-DATA).
