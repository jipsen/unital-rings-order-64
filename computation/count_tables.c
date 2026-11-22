/*
 * T(n) = number of associative multiplication tables on F_2^n that are
 * strictly triangular: e_i * e_j lies in span{ e_k : k > max(i,j) }.
 * These are exactly the labelled nilpotent algebras for which the standard
 * flag span(e_i..e_{n-1}) is "adapted" (A*V_i + V_i*A <= V_{i+1}).
 *
 * Mass formula used for the completeness check:
 *   sum over iso classes A of f(A)/|Aut(A)| = T(n) / 2^(n(n-1)/2)
 * where f(A) = number of adapted complete flags of A.
 *
 * Products are assigned level by level (level L = max(i,j), from deep to
 * shallow); each associativity triple is checked as soon as every product
 * it involves is known.
 */
#include <stdio.h>
#include <stdlib.h>

static int n;
static unsigned p[8][8];           /* p[i][j] = e_i * e_j as bitmask */
static long long count;

static unsigned mulvec(unsigned x, int k, int left) {
    /* left: (x) * e_k ; else e_k * (x) */
    unsigned r = 0;
    for (int m = 0; m < n; m++)
        if ((x >> m) & 1) r ^= left ? p[m][k] : p[k][m];
    return r;
}

static int check_level(int L) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++) {
                int a = i > j ? i : j, b = j > k ? j : k;
                int lo = a < b ? a : b;
                if (lo != L) continue;
                unsigned lhs = mulvec(p[i][j], k, 1);   /* (e_i e_j) e_k */
                unsigned rhs = mulvec(p[j][k], i, 0);   /* e_i (e_j e_k) */
                if (lhs != rhs) return 0;
            }
    return 1;
}

static void level(int L) {
    if (L < 0) { count++; return; }
    int pairs[16][2], np = 0;
    for (int i = 0; i <= L; i++)
        for (int j = 0; j <= L; j++)
            if ((i > j ? i : j) == L) { pairs[np][0] = i; pairs[np][1] = j; np++; }
    int w = n - 1 - L;                     /* bits available: indices L+1..n-1 */
    int total = np * w;
    for (long long m = 0; m < (1LL << total); m++) {
        for (int t = 0; t < np; t++) {
            unsigned v = (unsigned)((m >> (t * w)) & ((1u << w) - 1));
            p[pairs[t][0]][pairs[t][1]] = v << (L + 1);
        }
        if (check_level(L)) level(L - 1);
    }
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) p[i][j] = 0;
    level(n - 2 >= 0 ? n - 2 : -1);
    printf("%lld\n", count);
    return 0;
}
