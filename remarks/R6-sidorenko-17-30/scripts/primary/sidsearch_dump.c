/* Exhaustive check of Prop. prop:complex (2),(3) for WL-discrete sphere triangulations.
   Input lines: V F a b c a b c ...  (faces). Output: FOUND lines + summary.
   Build: gcc -O2 -o sidsearch sidsearch.c ; Run: ./sidsearch < wl_12.txt          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int V, F, E, fa[64][3];
static u64 vmask[64], fadj[64];
static int eu[128], ev[128], ef1[128], ef2[128];
static u64 *cuts; static long ncuts, capc;
static long **elist; static long *elen;

static int connected(u64 S) {
  u64 r = S & (-S), prev = 0;
  while (r != prev) { prev = r; u64 x = r; while (x) { int f = __builtin_ctzll(x); x &= x - 1; r |= fadj[f] & S; } }
  return r == S;
}
static int good(u64 S) {
  u64 cov = 0, x = S; while (x) { int f = __builtin_ctzll(x); x &= x - 1; cov |= vmask[f]; }
  if (cov != ((1ULL << V) - 1)) return 0;
  return connected(S);
}
int main(int argc, char **argv) {
  long tri = 0, hit = 0, maxcuts = 0; capc = 1 << 22; cuts = malloc(capc * 8);
  while (scanf("%d %d", &V, &F) == 2) {
    for (int j = 0; j < F; j++) scanf("%d %d %d", &fa[j][0], &fa[j][1], &fa[j][2]);
    tri++;
    for (int j = 0; j < F; j++) { vmask[j] = (1ULL << fa[j][0]) | (1ULL << fa[j][1]) | (1ULL << fa[j][2]); fadj[j] = 0; }
    E = 0;
    for (int a = 0; a < F; a++) for (int b = a + 1; b < F; b++)
      if (__builtin_popcountll(vmask[a] & vmask[b]) == 2) { fadj[a] |= 1ULL << b; fadj[b] |= 1ULL << a; ef1[E] = a; ef2[E] = b; E++; }
    if (E != 3 * V - 6) { fprintf(stderr, "bad E\n"); return 1; }
    u64 allF = (F == 64) ? ~0ULL : ((1ULL << F) - 1);
    u64 allE = (1ULL << E) - 1;
    ncuts = 0;
    int k = V - 3, n = F - 1;           /* choose k of faces 1..F-1 */
    u64 c = (1ULL << k) - 1, lim = 1ULL << n;
    while (c < lim) {
      u64 S = 1ULL | (c << 1);
      if (good(S) && good(allF ^ S)) {
        u64 cut = 0;
        for (int e = 0; e < E; e++) if (((S >> ef1[e]) & 1) != ((S >> ef2[e]) & 1)) cut |= 1ULL << e;
        if (ncuts == capc) { capc *= 2; cuts = realloc(cuts, capc * 8); }
        cuts[ncuts++] = cut;
      }
      u64 t = c | (c - 1); c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(c) + 1));
    }
    if (ncuts > maxcuts) maxcuts = ncuts;
    if (getenv("DUMP")) { printf("T %ld %d %d", tri, V, E); for (long i = 0; i < ncuts; i++) printf(" %llx", cuts[i]); printf("\n"); for (int j = 0; j < F; j++) fprintf(stderr, ""); continue; }
    /* per-edge lists */
    elist = malloc(E * sizeof(long *)); elen = calloc(E, sizeof(long));
    for (int e = 0; e < E; e++) elist[e] = malloc((ncuts + 1) * sizeof(long));
    for (long i = 0; i < ncuts; i++) for (int e = 0; e < E; e++) if (cuts[i] >> e & 1) elist[e][elen[e]++] = i;
    int found = 0; long fi = -1, fj = -1, fk = -1;
    for (long i = 0; i < ncuts && !found; i++) for (long j = i; j < ncuts && !found; j++) {
      u64 u = cuts[i] | cuts[j], miss = allE & ~u;
      if (!miss) { /* any third cut works; need a doubly-cut edge */
        for (long kk = 0; kk < ncuts; kk++) if ((cuts[i] & cuts[j]) | (cuts[kk] & u)) { found = 1; fi = i; fj = j; fk = kk; break; }
        continue;
      }
      int e = __builtin_ctzll(miss);
      for (long t = 0; t < elen[e]; t++) { long kk = elist[e][t]; u64 ck = cuts[kk];
        if ((ck & miss) == miss && ((cuts[i] & cuts[j]) | (ck & u))) { found = 1; fi = i; fj = j; fk = kk; break; } }
    }
    if (found) {
      hit++;
      printf("FOUND V=%d cuts=%ld faces:", V, ncuts);
      for (int j = 0; j < F; j++) printf(" %d,%d,%d", fa[j][0], fa[j][1], fa[j][2]);
      printf(" | cuts %llx %llx %llx\n", cuts[fi], cuts[fj], cuts[fk]); fflush(stdout);
    }
    for (int e = 0; e < E; e++) free(elist[e]); free(elist); free(elen);
  }
  printf("SUMMARY triangulations=%ld satisfying(2)(3)=%ld max_splits=%ld\n", tri, hit, maxcuts);
  return 0;
}
