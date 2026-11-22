/*
 * For marked nilpotent rings (J,c) given by pc presentations:
 *   mode "stats": verify ring axioms + marking on all elements, and print
 *                 |Aut(J,c)|, f(J) (adapted chains of ideals) and invariants
 *   mode "iso"  : read pairs of indices, print 1 if isomorphic (as marked rings)
 * Input file: N, then per ring: m c, P[0..m-1], Q[0..m-1][0..m-1].
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXM 6
#define MAXN 64

typedef struct {
    int m, n;
    unsigned P[MAXM], Q[MAXM][MAXM], c;
    unsigned char add[MAXN][MAXN], mul[MAXN][MAXN];
} Ring;

static Ring *R;
static int NR;

static unsigned addgen(Ring *r, unsigned x, int l) {
    if (!((x >> l) & 1)) return x | (1u << l);
    x &= ~(1u << l);
    for (int q = 0; q < r->m; q++) if ((r->P[l] >> q) & 1) x = addgen(r, x, q);
    return x;
}
static void build(Ring *r) {
    r->n = 1 << r->m;
    for (unsigned a = 0; a < (unsigned)r->n; a++)
        for (unsigned b = 0; b < (unsigned)r->n; b++) {
            unsigned x = a;
            for (int l = 0; l < r->m; l++) if ((b >> l) & 1) x = addgen(r, x, l);
            r->add[a][b] = x;
        }
    for (unsigned a = 0; a < (unsigned)r->n; a++)
        for (unsigned b = 0; b < (unsigned)r->n; b++) {
            unsigned acc = 0;
            for (int i = 0; i < r->m; i++) if ((a >> i) & 1)
                for (int j = 0; j < r->m; j++) if ((b >> j) & 1) acc = r->add[acc][r->Q[i][j]];
            r->mul[a][b] = acc;
        }
}
static int axioms(Ring *r) {
    int n = r->n;
    for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) for (int z = 0; z < n; z++) {
        if (r->mul[r->add[x][y]][z] != r->add[r->mul[x][z]][r->mul[y][z]]) return 0;
        if (r->mul[x][r->add[y][z]] != r->add[r->mul[x][y]][r->mul[x][z]]) return 0;
        if (r->mul[r->mul[x][y]][z] != r->mul[x][r->mul[y][z]]) return 0;
    }
    for (int x = 0; x < n; x++)
        if (r->mul[r->c][x] != r->add[x][x] || r->mul[x][r->c] != r->add[x][x]) return 0;
    /* triangular presentation (needed for the mass formula bookkeeping) */
    for (int i = 0; i < r->m; i++) {
        if (r->P[i] & ((2u << i) - 1)) return 0;
        for (int j = 0; j < r->m; j++) {
            int mx = i > j ? i : j;
            if (r->Q[i][j] & ((2u << mx) - 1)) return 0;
        }
    }
    return 1;
}

/* ---------- homomorphism search A -> B (images of generators) ---------- */
static Ring *HA, *HB;
static unsigned img[MAXM];
static long long found;
static int stop_first;

static unsigned applyB(unsigned x) {
    unsigned r = 0;
    for (int q = 0; q < HA->m; q++) if ((x >> q) & 1) r = HB->add[r][img[q]];
    return r;
}
static void hsearch(int i, uint64_t S) {   /* S = image of <g_{i+1}..> as element set */
    if (stop_first && found) return;
    if (i < 0) { if (applyB(HA->c) == HB->c) found++; return; }
    for (unsigned h = 0; h < (unsigned)HB->n; h++) {
        if ((S >> h) & 1) continue;
        if (HB->add[h][h] != applyB(HA->P[i])) continue;
        img[i] = h;
        int ok = 1;
        for (int j = i; j < HA->m && ok; j++) {
            unsigned hj = img[j];
            if (HB->mul[h][hj] != applyB(HA->Q[i][j])) ok = 0;
            else if (HB->mul[hj][h] != applyB(HA->Q[j][i])) ok = 0;
        }
        if (ok) {
            uint64_t T = S;
            for (unsigned x = 0; x < (unsigned)HB->n; x++) if ((S >> x) & 1) T |= 1ULL << HB->add[h][x];
            hsearch(i - 1, T);
        }
        img[i] = 0;
    }
}
static long long count_homs(Ring *a, Ring *b, int first) {
    if (a->m != b->m) return 0;
    HA = a; HB = b; found = 0; stop_first = first;
    memset(img, 0, sizeof img);
    hsearch(a->m - 1, 1ULL);
    return found;
}

/* ---------- f(J): chains of ideals with index-2 steps ---------- */
static Ring *FR;
#define HS (1 << 16)
static uint64_t hkey[HS]; static long long hval[HS]; static char hused[HS];

static uint64_t closure(uint64_t gens) {
    uint64_t S = 1ULL;
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int x = 0; x < FR->n; x++) if ((S >> x) & 1)
            for (int g = 0; g < FR->n; g++) if ((gens >> g) & 1) {
                int y = FR->add[x][g];
                if (!((S >> y) & 1)) { S |= 1ULL << y; changed = 1; }
            }
    }
    return S;
}
static long long fcount(uint64_t I) {
    if (I == 1ULL) return 1;
    unsigned hh = (unsigned)((I * 0x9E3779B97F4A7C15ULL) >> 48) & (HS - 1);
    while (hused[hh]) { if (hkey[hh] == I) return hval[hh]; hh = (hh + 1) & (HS - 1); }
    uint64_t gens = 0;
    for (int y = 0; y < FR->n; y++) if ((I >> y) & 1) {
        gens |= 1ULL << FR->add[y][y];
        for (int x = 0; x < FR->n; x++) {
            gens |= 1ULL << FR->mul[x][y];
            gens |= 1ULL << FR->mul[y][x];
        }
    }
    uint64_t W = closure(gens);
    /* basis of I/W (elementary abelian since 2I <= W) */
    unsigned b[MAXM]; int d = 0;
    uint64_t S = W;
    while (S != I) {
        int x = 0;
        while (!((I >> x) & 1) || ((S >> x) & 1)) x++;
        b[d++] = x;
        uint64_t g = 0;
        for (int k = 0; k < d; k++) g |= 1ULL << b[k];
        for (int y = 0; y < FR->n; y++) if ((W >> y) & 1) g |= 1ULL << y;
        S = closure(g);
    }
    /* coordinates eps(x) of each element of I */
    int eps[MAXN];
    for (int e = 0; e < (1 << d); e++) {
        unsigned s = 0;
        for (int k = 0; k < d; k++) if ((e >> k) & 1) s = FR->add[s][b[k]];
        for (int w = 0; w < FR->n; w++) if ((W >> w) & 1) eps[FR->add[w][s]] = e;
    }
    long long tot = 0;
    for (int lam = 1; lam < (1 << d); lam++) {
        uint64_t H = 0;
        for (int x = 0; x < FR->n; x++)
            if (((I >> x) & 1) && !(__builtin_popcount(lam & eps[x]) & 1)) H |= 1ULL << x;
        tot += fcount(H);
    }
    hh = (unsigned)((I * 0x9E3779B97F4A7C15ULL) >> 48) & (HS - 1);
    while (hused[hh]) hh = (hh + 1) & (HS - 1);
    hused[hh] = 1; hkey[hh] = I; hval[hh] = tot;
    return tot;
}

/* ---------- invariants ---------- */
static void invariants(Ring *r, long long *inv) {
    int n = r->n;
    FR = r;
    uint64_t prods = 0;
    for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) prods |= 1ULL << r->mul[x][y];
    uint64_t J2 = closure(prods);
    uint64_t p3 = 0;
    for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) if ((J2 >> y) & 1) p3 |= 1ULL << r->mul[x][y];
    uint64_t J3 = closure(p3);
    int ann = 0, sq0 = 0, two0 = 0, zp = 0, comm = 1;
    for (int x = 0; x < n; x++) {
        int a = 1;
        for (int y = 0; y < n; y++) {
            if (r->mul[x][y] || r->mul[y][x]) a = 0;
            if (!r->mul[x][y]) zp++;
            if (r->mul[x][y] != r->mul[y][x]) comm = 0;
        }
        ann += a; sq0 += !r->mul[x][x]; two0 += !r->add[x][x];
    }
    int ordc = 1; { unsigned y = r->c; while (y) { y = r->add[y][r->c]; ordc++; } if (!r->c) ordc = 1; }
    inv[0] = __builtin_popcountll(J2); inv[1] = __builtin_popcountll(J3);
    inv[2] = ann; inv[3] = sq0; inv[4] = two0; inv[5] = zp; inv[6] = comm; inv[7] = ordc;
}

int main(int argc, char **argv) {
    FILE *fp = fopen(argv[2], "r");
    if (fscanf(fp, "%d", &NR) != 1) return 1;
    R = calloc(NR, sizeof(Ring));
    for (int k = 0; k < NR; k++) {
        Ring *r = &R[k];
        if (fscanf(fp, "%d %u", &r->m, &r->c) != 2) return 1;
        for (int i = 0; i < r->m; i++) if (fscanf(fp, "%u", &r->P[i]) != 1) return 1;
        for (int i = 0; i < r->m; i++) for (int j = 0; j < r->m; j++)
            if (fscanf(fp, "%u", &r->Q[i][j]) != 1) return 1;
        build(r);
    }
    fclose(fp);
    if (!strcmp(argv[1], "stats")) {
        for (int k = 0; k < NR; k++) {
            Ring *r = &R[k];
            int ok = axioms(r);
            long long aut = count_homs(r, r, 0);
            FR = r; memset(hused, 0, sizeof hused);
            long long f = fcount(r->n == 64 ? ~0ULL : ((1ULL << r->n) - 1));
            long long inv[8]; invariants(r, inv);
            printf("%d %lld %lld", ok, aut, f);
            for (int t = 0; t < 8; t++) printf(" %lld", inv[t]);
            printf("\n");
        }
    } else {  /* iso: pairs on stdin */
        int a, b;
        while (scanf("%d %d", &a, &b) == 2)
            printf("%d %d %d\n", a, b, count_homs(&R[a], &R[b], 1) > 0);
    }
    return 0;
}
