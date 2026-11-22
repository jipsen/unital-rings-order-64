/*
 * T(m) = number of "triangular marked presentations" of order 2^m:
 *   generators g_0..g_{m-1}, every element = sum of a subset (bitmask),
 *   2 g_i     = P[i]     with support in {i+1..m-1}
 *   g_i g_j   = Q[i][j]  with support in {max(i,j)+1..m-1}
 *   marking c = any bitmask
 * such that the bilinear product is well defined, associative, and
 * c*x = x*c = 2x for all x.
 *
 * Every marked nilpotent ring (J,c) of order 2^m arises from
 * f(J) * 2^(m(m-1)/2) / |Aut(J,c)| distinct presentations, where f(J) is the
 * number of chains J = I_0 > I_1 > ... > I_m = 0 of index-2 subgroups with
 * J I_k + I_k J + 2 I_k <= I_{k+1}.  Hence
 *     sum_{(J,c)} f(J)/|Aut(J,c)| = T(m) / 2^(m(m-1)/2).
 *
 * Search: c first, then levels L = m-1 .. 0 (level L holds P[L] and all
 * Q[i][j] with max(i,j) = L); a condition is tested as soon as every datum
 * it uses has been assigned.
 */
#include <stdio.h>
#include <stdlib.h>

static int m;
static unsigned P[8], Q[8][8], C;
static long long count;

static unsigned addgen(unsigned x, int l) {
    if (!((x >> l) & 1)) return x | (1u << l);
    x &= ~(1u << l);
    for (int q = l + 1; q < m; q++)
        if ((P[l] >> q) & 1) x = addgen(x, q);
    return x;
}
static unsigned add(unsigned x, unsigned y) {
    for (int l = 0; l < m; l++) if ((y >> l) & 1) x = addgen(x, l);
    return x;
}
static unsigned mul(unsigned x, unsigned y) {
    unsigned r = 0;
    for (int i = 0; i < m; i++) if ((x >> i) & 1)
        for (int j = 0; j < m; j++) if ((y >> j) & 1) r = add(r, Q[i][j]);
    return r;
}

static int check_level(int L) {
    unsigned gL = 1u << L;
    /* well-definedness: (2g_L) g_j = 2 (g_L g_j), g_j (2g_L) = 2 (g_j g_L) */
    for (int j = 0; j < m; j++) {
        unsigned gj = 1u << j;
        unsigned a = mul(P[L], gj), b = Q[L][j];
        if (a != add(b, b)) return 0;
        a = mul(gj, P[L]); b = Q[j][L];
        if (a != add(b, b)) return 0;
    }
    /* marking on g_L */
    unsigned two = P[L];
    if (mul(C, gL) != two || mul(gL, C) != two) return 0;
    /* associativity of triples whose shallowest datum is at level L */
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) for (int k = 0; k < m; k++) {
        int a = i > j ? i : j, b = j > k ? j : k;
        if ((a < b ? a : b) != L) continue;
        if (mul(Q[i][j], 1u << k) != mul(1u << i, Q[j][k])) return 0;
    }
    return 1;
}

static void level(int L) {
    if (L < 0) { count++; return; }
    int pairs[16][2], np = 0;
    for (int i = 0; i <= L; i++)
        for (int j = 0; j <= L; j++)
            if ((i > j ? i : j) == L) { pairs[np][0] = i; pairs[np][1] = j; np++; }
    int w = m - 1 - L;
    int total = (np + 1) * w;              /* products + power relation */
    for (long long s = 0; s < (1LL << total); s++) {
        unsigned msk = (1u << w) - 1;
        for (int t = 0; t < np; t++)
            Q[pairs[t][0]][pairs[t][1]] = (unsigned)((s >> (t * w)) & msk) << (L + 1);
        P[L] = (unsigned)((s >> (np * w)) & msk) << (L + 1);
        if (check_level(L)) level(L - 1);
    }
    for (int t = 0; t < np; t++) Q[pairs[t][0]][pairs[t][1]] = 0;
    P[L] = 0;
}

int main(int argc, char **argv) {
    m = atoi(argv[1]);
    long long total = 0;
    for (C = 0; C < (1u << m); C++) { count = 0; level(m - 1); total += count; }
    printf("%lld\n", total);
    return 0;
}
