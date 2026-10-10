/* Independent exhaustive search for perfect mu_q-valued functions on a finite
 * abelian group G = Z_{n1} x Z_{n2} x Z_{n3}.  Written for this remark; it shares no
 * code with the upstream artifacts.
 *
 * f : G -> mu_q is perfect iff sum_g f(g) conj(f(g+s)) = 0 for all s != 0.
 * Encode f(g) = zeta_q^{a(g)}.  For q in {2,3,4} the value
 * sum_k c_k zeta_q^k (c_k = #{g : a(g)-a(g+s) = k mod q}) vanishes iff
 *   q=2: c0 = c1;   q=3: c0 = c1 = c2;   q=4: c0 = c2 and c1 = c3.
 * Normalisation a(0) = 0 (multiplying f by a constant keeps it perfect), so the
 * search visits q^{|G|-1} functions and reports the number of normalised perfect ones.
 *
 * usage: perfect_search q n1 [n2 [n3]]      (single thread)
 */
#include <stdio.h>
#include <stdlib.h>

static int q, n[3], N;
static int add_tab[64][64];

static int vanishes(const int *c) {
    if (q == 2) return c[0] == c[1];
    if (q == 3) return c[0] == c[1] && c[1] == c[2];
    if (q == 4) return c[0] == c[2] && c[1] == c[3];
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s q n1 [n2 [n3]]\n", argv[0]); return 2; }
    q = atoi(argv[1]);
    n[0] = atoi(argv[2]); n[1] = argc > 3 ? atoi(argv[3]) : 1; n[2] = argc > 4 ? atoi(argv[4]) : 1;
    N = n[0] * n[1] * n[2];
    if (N > 64 || (q != 2 && q != 3 && q != 4)) { fprintf(stderr, "unsupported\n"); return 2; }
    for (int g = 0; g < N; g++)
        for (int s = 0; s < N; s++) {
            int g0 = g % n[0], g1 = (g / n[0]) % n[1], g2 = g / (n[0] * n[1]);
            int s0 = s % n[0], s1 = (s / n[0]) % n[1], s2 = s / (n[0] * n[1]);
            add_tab[g][s] = (g0 + s0) % n[0] + n[0] * (((g1 + s1) % n[1]) + n[1] * ((g2 + s2) % n[2]));
        }
    int a[64] = {0};
    unsigned long long count = 0, visited = 0;
    for (;;) {
        visited++;
        int ok = 1;
        for (int s = 1; s < N && ok; s++) {
            int c[4] = {0, 0, 0, 0};
            for (int g = 0; g < N; g++) c[((a[g] - a[add_tab[g][s]]) % q + q) % q]++;
            ok = vanishes(c);
        }
        if (ok) count++;
        int i = 1;                      /* odometer over a[1..N-1] */
        while (i < N && ++a[i] == q) a[i++] = 0;
        if (i == N) break;
    }
    printf("q=%d G=Z%d x Z%d x Z%d |G|=%d visited=%llu perfect(normalised a(0)=0)=%llu\n",
           q, n[0], n[1], n[2], N, visited, count);
    return 0;
}
