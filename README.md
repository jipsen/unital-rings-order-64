# Unital rings of order 64

Goal: count the unital (associative) rings with 64 elements up to isomorphism.
The sequence is [OEIS A037291](https://oeis.org/A037291); its known terms
stop at $`n = 63`$ ($`a(32) = 208`$). The code, data and descriptions in this repository 
were produced by Claude Opus 5.5 under a free Claude Team plan for scientists.

## Current status

| | count |
|---|---:|
| Local rings of order 64 | **1,789** |
| &nbsp;&nbsp;residue field $`\mathbb{F}_2`$ | 1,773 |
| &nbsp;&nbsp;residue field $`\mathbb{F}_4`$ | 11 |
| &nbsp;&nbsp;residue field $`\mathbb{F}_8`$ | 4 |
| &nbsp;&nbsp;residue field $`\mathbb{F}_{64}`$ | 1 |
| Commutative local rings of order 64 | 259 |
| **Commutative unital rings of order 64** | **462** |
| **All unital rings of order 64** | $`2134 + N`$ |

$`N`$ is the number of indecomposable **non-local** rings of order 64
($`R/J`$ a product of at least two fields). That is the one piece still
missing, so currently

```math
a(64) \;\ge\; 2134 .
```

The commutative count is complete, because every commutative indecomposable
finite ring is local. (OEIS [A127707](https://oeis.org/A127707) also stops at
$`n = 63`$.)

## Reduction

### 1. Indecomposables

A finite unital ring splits uniquely (via central idempotents) into
indecomposable rings, so the totals are the Euler transform of the
indecomposable counts $`i_k`$:

```math
\sum_{n \ge 0} a(2^n)\, x^n \;=\; \prod_{k \ge 1} \left(1 - x^k\right)^{-i_k}.
```

From $`a(2^k) = 1, 4, 11, 50, 208`$ ($`k = 1, \dots, 5`$) the indecomposable
counts are $`i_k = 1, 3, 7, 33, 137`$, and therefore

```math
a(64) \;=\; 345 + i_6 .
```

Commutative: the indecomposable counts are $`1, 3, 6, 21, 54`$ (from
[A127707](https://oeis.org/A127707)), giving $`203 + c_6`$, where $`c_6`$ is
the number of commutative local rings of order 64.

### 2. No matrix blocks

If $`R/J(R)`$ had a block $`M_2(\mathbb{F}_2)`$, lifting matrix units makes
that corner $`M_2(S)`$, of order $`\lvert S \rvert^4 \in \{16, 256\}`$ (never
64), and gluing anything to it needs at least 128 elements. So an
indecomposable ring of order 64 is local, or basic non-local ($`R/J`$ a
product of $`\ge 2`$ fields).

### 3. Local rings with residue field $`\mathbb{F}_2`$

These correspond to *marked pairs* $`(J, c)`$: $`J = J(R)`$ is a nilpotent ring
of order $`\lvert R \rvert / 2`$ and $`c = 2 \cdot 1_R \in J`$ satisfies
$`cx = xc = 2x`$ for all $`x \in J`$. Conversely

```math
R \;\cong\; (\mathbb{Z} \oplus J) \,/\, \langle (2, -c) \rangle,
\qquad
R \cong R' \iff (J, c) \cong (J', c'),
```

and $\text{char} R = 2 \text{ord}(c)$. Characteristic 2 is
$`c = 0`$, i.e. $`J`$ a nilpotent $`\mathbb{F}_2`$-algebra.

### 4. Local rings with residue field $`\mathbb{F}_q`$, $`q = 4, 8`$

These correspond to a marked pair $`(J, c)`$ of order $`\lvert R \rvert / q`$
together with the left and right multiplications $`\lambda(y) = \zeta y`$,
$`\rho(y) = y \zeta`$ by a Teichmüller unit $`\zeta`$ of order $`q - 1`$,
subject to

```math
h(\lambda) = h(\rho) = 0, \quad
\lambda\rho = \rho\lambda, \quad
\lambda(c) = \rho(c), \quad
\lambda(xy) = \lambda(x)\, y, \quad
\rho(xy) = x\, \rho(y), \quad
\rho(x)\, y = x\, \lambda(y),
```

where $`h`$ is the 2-adic lift of the minimal polynomial of
$`\bar\zeta \in \mathbb{F}_q`$ ($`h = x^2 + x + 1`$ for $`q = 4`$). Then

```math
R \;\cong\; \bigl(\mathbb{Z}[x]/(h) \ltimes J\bigr) \,/\, \{(2s, -s c)\},
```

and isomorphism classes are the orbits of
$`\text{Aut}(J, c) \times \langle \zeta \mapsto \zeta^2 \rangle`$.

## Results

### Local rings with residue field $`\mathbb{F}_2`$

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

### Local rings with residue field $`\mathbb{F}_4`$, $`\mathbb{F}_8`$, $`\mathbb{F}_{64}`$

| residue | order | rings | char 2 (comm / noncomm) | char 4 (comm / noncomm) | char 8 |
|---|---:|---:|---|---|---:|
| $`\mathbb{F}_4`$ | 16 | 3 | 1 / 1 | 1 / 0 | |
| $`\mathbb{F}_4`$ | **64** | **11** | 2 / 3 | 2 / 3 | 1 ($`\mathrm{GR}(8,2)`$) |
| $`\mathbb{F}_8`$ | **64** | **4** | 1 / 2 | 1 / 0 ($`\mathrm{GR}(4,3)`$) | |
| $`\mathbb{F}_{64}`$ | 64 | 1 | 1 / 0 | | |

## Methods

**Nilpotent $`\mathbb{F}_2`$-algebras** (`enum_nilpotent.py`).
Skjelbred–Sund central extensions. Algebras without an annihilator component
are $`B_\theta`$ with $`B = A / \text{Ann}(A)`$, where $`\theta`$ runs
over $`\text{Aut}(B)`$-orbits of subspaces of $`H^2(B, \mathbb{F}_2)`$
satisfying $`\text{Ann}(B) \cap \text{Ann}(\theta) = 0`$.

**Marked rings, all characteristics** (`enum_local_F2.py`).
Central extensions by the full annihilator $`Z = \text{Ann}(J)`$, which
is elementary abelian since $`2z = cz = 0`$. For
$`(K, \bar c) = (J/Z,\, c \bmod Z)`$ in a fixed polycyclic presentation, $`J`$
is given by additive data $`\tau`$, product data $`\theta`$ and the offset
$`z_c`$ of $`c`$:

```math
2 g_i = P_i + \tau_i, \qquad
g_i g_j = Q_{ij} + \theta_{ij}, \qquad
c = \bar c + z_c ,
```

subject to linear conditions; changing the lifts $`g_i`$ gives the
coboundaries. With $`H`$ = cocycles / coboundaries, isomorphism classes
correspond to $`\text{Aut}(K, \bar c)`$-orbits of subspaces
$`W \le H`$ with $`\dim W \le \dim Z`$ and
$`\text{Ann}(K) \cap \text{Ann}(W) = 0`$.

**Residue fields $`\mathbb{F}_4`$, $`\mathbb{F}_8`$** (`enum_residue_Fq.py`).
Enumerate all additive endomorphisms $`\lambda, \rho`$ of $`J`$ satisfying the
conditions in step 4 above, then take orbits.

## Verification

**Known values reproduced.** Commutative local rings with residue field
$`\mathbb{F}_2`$: $`2, 5, 18, 53`$ for orders 4–32. Characteristic-2 subset:
$`1, 2, 7, 27, 404`$. Order 16 is consistent with the 13 noncommutative unital
rings of order 16 ([A127708](https://oeis.org/A127708)): 4 local with residue
field $`\mathbb{F}_2`$, 7 non-local (including $`M_2(\mathbb{F}_2)`$),
$`\mathbb{F}_4[x;\sigma]/(x^2)`$, and $`T_2(\mathbb{F}_2) \times \mathbb{F}_2`$.
Residue field $`\mathbb{F}_4`$ at order 16 gives the expected 3; residue field
$`\mathbb{F}_8`$ at order 64 gives the expected 4.

**No duplicates.** The 404 $`\mathbb{F}_2`$-algebras have pairwise distinct
canonical forms under GAP/ModIsom. For the 1,773 marked rings, all 11,508 pairs
with equal invariants were ruled out by exhaustive isomorphism search (which
finds 60/60 randomly re-presented copies as a control). Same for the
residue-$`\mathbb{F}_4`$/$`\mathbb{F}_8`$ rings.

**Nothing missing (mass formula).** For $`m = 1, \dots, 5`$,

```math
\sum_{(J,c)} \frac{f(J)}{\lvert \text{Aut}(J, c) \rvert}
\;=\; \frac{T(m)}{2^{m(m-1)/2}}
```

holds exactly, where the sum runs over the listed marked rings of order
$`2^m`$, $`T(m)`$ counts triangular marked presentations
($`2, 8, 112, 5824, 2292096`$; the C count matches brute force for
$`m \le 3`$) and $`f(J)`$ counts adapted chains of ideals. The analogous
identity for $`\mathbb{F}_2`$-algebras holds for $`n = 1, \dots, 5`$.

**Axioms.** Every ring is checked on all elements (the residue-$`\mathbb{F}_4`$
and $`\mathbb{F}_8`$ rings are built as explicit 64-element tables and checked
to be local with the right radical). Every orbit satisfies orbit–stabilizer
against the full group.

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
  bitmask of $`e_i e_j`$ in the basis $`e_0, \dots, e_{n-1}`$ of a nilpotent
  $`\mathbb{F}_2`$-algebra. Dimension 5 gives the 404 characteristic-2 local
  rings of order 64.
* `local_F2_upto5.json` — `{m: [{m, P, Q, c}, ...]}`: marked rings in
  polycyclic presentation. Generators $`g_0, \dots, g_{m-1}`$; each element is
  a subset sum (bitmask); $`2 g_i`$ = `P[i]`, $`g_i g_j`$ = `Q[i][j]`,
  $`c = 2 \cdot 1`$. The local ring is
  $`R = (\mathbb{Z} \oplus J)/\langle (2, -c) \rangle`$ of order $`2^{m+1}`$.
  Key `"5"` holds the 1,773 rings of order 64.
* `residue_F4_order64.json`, `residue_F8_order64.json` — each entry has the
  marked ring `J`, the value tables `L` ($`\lambda : y \mapsto \zeta y`$) and
  `R` ($`\rho : y \mapsto y \zeta`$) on $`J`$, and the invariants
  ($`\text{char} R`$, commutative, $`\lvert Z(R) \rvert`$,
  $`\lvert J^2 \rvert`$, $`\#\{x \in J : x^2 = 0\}`$,
  $`\#\{(x, y) \in J \times J : xy = 0\}`$, and the sizes of the left and right
  annihilators of $`J`$ in $`J`$).
* `autsizes_n.txt`, `T_marked.txt` — automorphism group orders (ModIsom) and
  presentation counts used by the mass-formula checks.

## Remaining work

* Non-local indecomposable rings of order 64 (Peirce decomposition /
  quivers with Galois-ring coefficients). Validation targets: 7 at order 16
  and 31 at order 32 (from $`i_4 = 33`$, $`i_5 = 137`$).
* Optional cross-check of the local counts against Akers & Szabo,
  "A partial classification of local rings of order $`p^6`$", Involve 16
  (2023) 151–165.

## References

* OEIS A037291, A127707, A127708.
* B. Corbas, G. D. Williams, Rings of order $`p^5`$, Parts I and II,
  J. Algebra 231 (2000).
* C. Akers, S. Szabo, A partial classification of local rings of order
  $`p^6`$, Involve 16 (2023) 151–165.
* B. Eick, ModIsom (GAP package).
* W. A. de Graaf, Classification of nilpotent associative algebras of small
  dimension, IJAC 28 (2018) 133–161.

These are computational results and have not been peer reviewed.

## License

* **Code** (`computation/*.py`, `*.c`, `*.g`, `*.sh`): MIT License — see
  [LICENSE](LICENSE).
* **Data** (`computation/*.json`, `computation/*.txt`): Creative Commons
  Attribution 4.0 International (CC BY 4.0) — see [LICENSE-DATA](LICENSE-DATA).
