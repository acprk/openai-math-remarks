// kron_stair.cpp -- direct computation of the Kronecker coefficients
//   g(rho_m, rho_m, lambda)  (mod a prime P)  for ALL lambda |- N = m(m+1)/2,
// together with the exact character values a(lambda) = chi^lambda(hhat_m),
// where hhat_m = (2m-1, 2m-5, 2m-9, ...) is the cycle type given by the
// diagonal hook lengths of the staircase rho_m = (m, m-1, ..., 1).
//
// This code is written from scratch for remark R15 and does NOT use the
// 205 recursion (Prop. 3.2 / 4.2 / 5.3) or any Littlewood-Richardson code.
// It uses only:
//   * g(rho,rho,lambda) = sum_{mu |- N} chi^rho(mu)^2 chi^lambda(mu) / z_mu;
//   * chi^rho(mu) = 0 unless every part of mu is odd (rho_m is a 2-core, so it
//     has no hook of even length; apply Murnaghan-Nakayama to an even part first);
//   * the Murnaghan-Nakayama rule  p_k s_nu = sum_{lambda/nu k-border strip}
//     (-1)^{ht} s_lambda, implemented with beta-sets.
//
// Algorithm.  Split every odd partition mu = head u tail, head = parts > T,
// tail = parts <= T (T odd).  Write V_head = p_head in the Schur basis
// (vector over partitions of |head|), Y_tail[nu] = <p_tail s_nu, s_rho>.
//   chi^rho(mu) = <V_head, Y_tail>,  z_mu = z_head z_tail,
//   F := sum_mu c_mu p_mu = sum_tail p_tail A_tail,
//   A_tail := sum_head c_{head u tail} V_head,
// and the outer sum is evaluated by Horner along the tail trie.
// Then F = sum_lambda g(rho,rho,lambda) s_lambda (mod P).
//
// Usage:  ./kron_stair m P T threads outfile
//   outfile: summary lines (#), one line per hook (a, g mod P), and COLLIDE lines for any
//   lambda with g == -eps*a (mod P), -eps*a >= 0 (none => chi_+ certified universal when eps=+1).
//   For N <= 45 also outfile.full: "parts g_mod_P a" for every lambda.
// Build:  g++ -O3 -march=native -fopenmp -o kron_stair kron_stair.cpp
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <omp.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
typedef __int128 i128;

static int N;
static ll P;
static vector<vector<ll>> CNT; // CNT[s][k] = #partitions of s with parts <= k
static vector<ll> PN;          // PN[s] = p(s)

static void init_counts(int n) {
    CNT.assign(n + 1, vector<ll>(n + 1, 0));
    for (int k = 0; k <= n; k++) CNT[0][k] = 1;
    for (int s = 1; s <= n; s++)
        for (int k = 1; k <= n; k++)
            CNT[s][k] = CNT[s][k - 1] + (k <= s ? CNT[s - k][k] : 0);
    PN.resize(n + 1);
    for (int s = 0; s <= n; s++) PN[s] = CNT[s][s > 0 ? s : 0];
}
// rank in reverse-lexicographic order: (s) has rank 0, (1^s) has rank p(s)-1
static inline ll rank_of(const int* a, int L, int s) {
    ll r = 0; int rem = s, prev = s;
    for (int i = 0; i < L; i++) {
        int x = a[i]; int top = prev < rem ? prev : rem;
        r += CNT[rem][top] - CNT[rem][x];
        rem -= x; prev = x;
    }
    return r;
}
static int unrank(int s, ll r, int* a) {
    int L = 0, rem = s, prev = s;
    while (rem > 0) {
        int top = prev < rem ? prev : rem;
        for (int y = top; y >= 1; y--) {
            ll c = CNT[rem - y][y];
            if (r < c) { a[L++] = y; rem -= y; prev = y; break; }
            r -= c;
        }
    }
    return L;
}
static inline int next_part(int* a, int L) { // next in reverse-lex order; returns new L (0 if none)
    int i = L - 1;
    while (i >= 0 && a[i] == 1) i--;
    if (i < 0) return 0;
    int R = (L - 1 - i) + 1; int x = a[i] - 1; a[i] = x; L = i + 1;
    while (R > x) { a[L++] = x; R -= x; }
    if (R > 0) a[L++] = R;
    return L;
}
// enumerate k-border strips of lambda (parts a[0..L-1]); callback(rank_of_result, sign)
// result partition of size s = |lambda| - k.
template <class F>
static inline void strips(const int* a, int L, int k, int s, F f) {
    int B[160]; ull bits[4] = {0, 0, 0, 0};
    for (int i = 0; i < L; i++) { B[i] = a[i] + (L - 1 - i); bits[B[i] >> 6] |= 1ULL << (B[i] & 63); }
    int C[160], b2[160];
    for (int i = 0; i < L; i++) {
        int b = B[i], nb = b - k;
        if (nb < 0) continue;
        if (bits[nb >> 6] >> (nb & 63) & 1ULL) continue;
        // sign: number of beta numbers strictly between nb and b (they have index > i)
        int cnt = 0;
        for (int j = i + 1; j < L && B[j] > nb; j++) cnt++;
        // new beta set: replace B[i] by nb and restore descending order
        for (int j = 0; j < L; j++) b2[j] = B[j];
        b2[i] = nb;
        for (int j = i; j + 1 < L && b2[j + 1] > b2[j]; j++) { int tmp = b2[j]; b2[j] = b2[j + 1]; b2[j + 1] = tmp; }
        int Lc = 0;
        for (int j = 0; j < L; j++) { int part = b2[j] - (L - 1 - j); if (part > 0) C[Lc++] = part; }
        f(rank_of(C, Lc, s), (cnt & 1) ? -1 : 1);
    }
}
static int NT = 1;
// W (size p(s+k)) = p_k * V (V size p(s)); mod P if md else exact
static void mn_pull(int s, int k, const ll* V, ll* W, bool md) {
    ll tot = PN[s + k];
    int nth = NT;
    #pragma omp parallel num_threads(nth)
    {
        int tid = omp_get_thread_num(), nt = omp_get_num_threads();
        ll lo = tot * tid / nt, hi = tot * (tid + 1) / nt;
        int a[200]; int L = unrank(s + k, lo, a);
        for (ll r = lo; r < hi; r++) {
            ll acc = 0;
            strips(a, L, k, s, [&](ll rr, int sg) { acc += sg > 0 ? V[rr] : -V[rr]; });
            if (md) { acc %= P; if (acc < 0) acc += P; }
            W[r] = acc;
            if (r + 1 < hi) L = next_part(a, L);
        }
    }
}
static ll pw(ll b, ll e) { ll r = 1; b %= P; if (b < 0) b += P; while (e) { if (e & 1) r = (i128)r * b % P; b = (i128)b * b % P; e >>= 1; } return r; }
static ll inv(ll x) { return pw(x, P - 2); }

// ---------- tails ----------
struct Tail { vector<int> parts; int size; ll zinv; int parent; int x; };
static vector<Tail> tails;               // tail trie nodes (parts nonincreasing, each <= T)
static vector<vector<int>> tails_by_size;
static int T;
static void build_tails(int idx) {
    Tail& t0 = tails[idx];
    int lim = t0.parts.empty() ? T : t0.parts.back();
    for (int x = (lim % 2 ? lim : lim - 1); x >= 1; x -= 2) {
        if (tails[idx].size + x > N) continue;
        Tail nt; nt.parts = tails[idx].parts; nt.parts.push_back(x); nt.size = tails[idx].size + x;
        nt.parent = idx; nt.x = x; nt.zinv = 0;
        tails.push_back(nt); build_tails((int)tails.size() - 1);
    }
}
static ll z_of(const vector<int>& parts) { // z mod P
    map<int, int> m; for (int x : parts) m[x]++;
    ll z = 1;
    for (auto& kv : m) { for (int i = 0; i < kv.second; i++) z = (i128)z * (kv.first % P) % P * ((i + 1) % P) % P; }
    return z;
}
// sparse Y_tail over partitions of N-|tail| (subshapes of rho), mod P
static vector<unordered_map<ll, ll>> Y;
// transposed MN: Yout[nu] += sign * Yin[lambda] for every k-strip lambda/nu ; lambda |- s+k
static void mn_push_sparse(int s, int k, const unordered_map<ll, ll>& Yin, unordered_map<ll, ll>& Yout) {
    int a[200];
    for (auto& kv : Yin) {
        if (kv.second == 0) continue;
        int L = unrank(s + k, kv.first, a);
        strips(a, L, k, s, [&](ll rr, int sg) {
            ll& v = Yout[rr]; v = (v + (sg > 0 ? kv.second : P - kv.second)) % P; });
    }
}
// ---------- heads ----------
static vector<vector<ll>> A; // A[tail idx] dense over partitions of N - |tail|
static long long nleaves = 0;
static void head_dfs(vector<int>& head, ll zh_inv, const vector<ll>& V, int s) {
    // leaves: combine with every tail of size N-s
    int j = N - s;
    if (j >= 0 && j < (int)tails_by_size.size()) {
        for (int ti : tails_by_size[j]) {
            const auto& Yt = Y[ti];
            if (Yt.empty()) continue;
            i128 x = 0;
            for (auto& kv : Yt) x += (i128)V[kv.first] * kv.second;
            ll xm = (ll)(x % P);
            ll c = (i128)xm * xm % P * zh_inv % P * tails[ti].zinv % P;
            nleaves++;
            if (c == 0) continue;
            vector<ll>& At = A[ti];
            if (At.empty()) At.assign(PN[s], 0);
            ll* ap = At.data(); const ll* vp = V.data(); ll len = PN[s];
            #pragma omp parallel for num_threads(NT) schedule(static)
            for (ll i = 0; i < len; i++) ap[i] = (ll)(((i128)c * vp[i] + ap[i]) % P);
        }
    }
    int lim = head.empty() ? N : head.back();
    for (int k = (lim % 2 ? lim : lim - 1); k > T; k -= 2) {
        if (s + k > N) continue;
        // a head of size s+k must leave room for some tail: N-(s+k) >= 0 always ok
        vector<ll> W(PN[s + k]);
        mn_pull(s, k, V.data(), W.data(), true);
        head.push_back(k);
        // update z: multiplicity of k in head
        int mult = 0; for (int x : head) if (x == k) mult++;
        ll zi = (i128)zh_inv * inv((ll)k * mult % P) % P;
        head_dfs(head, zi, W, s + k);
        head.pop_back();
    }
}
// R(sigma) = A_sigma + sum_x p_x R(sigma+x) ; returns vector over partitions of N-|sigma|
static vector<vector<int>> children;
static vector<ll> horner(int ti) {
    int deg = N - tails[ti].size;
    vector<ll> R(PN[deg], 0);
    if (!A[ti].empty()) R = A[ti];
    A[ti].clear(); A[ti].shrink_to_fit();
    for (int ci : children[ti]) {
        vector<ll> Rc = horner(ci);
        int x = tails[ci].x;
        vector<ll> W(PN[deg]);
        mn_pull(deg - x, x, Rc.data(), W.data(), true);
        for (ll i = 0; i < PN[deg]; i++) { R[i] += W[i]; if (R[i] >= P) R[i] -= P; }
    }
    return R;
}
int main(int argc, char** argv) {
    if (argc < 6) { fprintf(stderr, "usage: m P T threads outfile\n"); return 1; }
    int m = atoi(argv[1]); P = atoll(argv[2]); T = atoi(argv[3]); NT = atoi(argv[4]);
    N = m * (m + 1) / 2;
    init_counts(N);
    fprintf(stderr, "m=%d N=%d p(N)=%lld P=%lld T=%d threads=%d\n", m, N, PN[N], P, T, NT);
    // tails
    tails.push_back(Tail{{}, 0, 1, -1, 0});
    build_tails(0);
    children.assign(tails.size(), {});
    for (size_t i = 1; i < tails.size(); i++) children[tails[i].parent].push_back((int)i);
    tails_by_size.assign(N + 1, {});
    for (size_t i = 0; i < tails.size(); i++) { tails[i].zinv = inv(z_of(tails[i].parts)); tails_by_size[tails[i].size].push_back((int)i); }
    // Y
    Y.assign(tails.size(), {});
    {
        int rho[200]; for (int i = 0; i < m; i++) rho[i] = m - i;
        Y[0][rank_of(rho, m, N)] = 1;
        for (size_t i = 1; i < tails.size(); i++) { // parents precede children (DFS order)
            int s = N - tails[i].size; int x = tails[i].x;
            mn_push_sparse(s, x, Y[tails[i].parent], Y[i]);
            for (auto it = Y[i].begin(); it != Y[i].end();) { if (it->second == 0) it = Y[i].erase(it); else ++it; }
        }
    }
    fprintf(stderr, "tails: %zu nodes\n", tails.size());
    A.assign(tails.size(), {});
    vector<ll> V0(1, 1); vector<int> head;
    head_dfs(head, 1, V0, 0);
    fprintf(stderr, "heads done, leaves(mu) = %lld\n", nleaves);
    vector<ll> G = horner(0);
    // exact a(lambda) = chi^lambda(hhat)
    vector<ll> H(1, 1); int s = 0;
    vector<int> hh; for (int h = 2 * m - 1; h > 0; h -= 4) hh.push_back(h);
    for (int h : hh) { vector<ll> W(PN[s + h]); mn_pull(s, h, H.data(), W.data(), false); H.swap(W); s += h; }
    int d = hh.size(); int eps = (((N - d) / 2) % 2 == 0) ? 1 : -1;
    // checks
    FILE* fo = fopen(argv[5], "w");
    int a[200];
    ll nneg = 0, nzero = 0, npos = 0, ncollide = 0, ncollide0 = 0, nself = 0;
    ll rho_rank; { int rho[200]; for (int i = 0; i < m; i++) rho[i] = m - i; rho_rank = rank_of(rho, m, N); }
    fprintf(fo, "# m=%d N=%d eps=%d P=%lld hhat=", m, N, eps, P);
    for (int h : hh) fprintf(fo, "%d,", h);
    fprintf(fo, "\n# g(rho,rho,rho) mod P = %lld ; a(rho) = %lld\n", G[rho_rank], H[rho_rank]);
    // hooks: print g mod P and a
    for (ll r = 0; r < PN[N]; r++) {
        int L = unrank(N, r, a);
        bool hook = (L <= 1) || a[1] <= 1;
        ll av = H[r];
        if (av < 0) nneg++; else if (av == 0) nzero++; else npos++;
        ll need = -eps * av; // mult of chi^lambda-constituent in chi_+^2 is (g - need)/2 (or /4); g >= need is automatic
        if (need >= 0 && ((need % P) == G[r])) {
            if (need == 0) ncollide0++; else ncollide++;
            fprintf(fo, "COLLIDE rank=%lld a=%lld gmodP=%lld lam=", r, av, G[r]); for (int i = 0; i < L; i++) fprintf(fo, "%d ", a[i]); fprintf(fo, "\n"); }
        if (hook) { fprintf(fo, "hook %d^1^%d a=%lld gmodP=%lld\n", a[0], L - 1, av, G[r]); }
    }
    fprintf(fo, "# #lambda=%lld  a<0:%lld a=0:%lld a>0:%lld\n", PN[N], nneg, nzero, npos);
    fprintf(fo, "# lambdas with eps*a<0 and g == |a| (mod P): %lld   (0 => every such g > |a|)\n", ncollide);
    fprintf(fo, "# lambdas with a==0 and g == 0 (mod P): %lld\n", ncollide0);
    // checksum 1: sum_lambda g(rho,rho,lambda) * dim(lambda) == dim(rho)^2  (mod P)
    {
        vector<ll> fact(N + 1, 1); for (int i = 1; i <= N; i++) fact[i] = (i128)fact[i - 1] * i % P;
        auto dimmod = [&](const int* q, int L) -> ll {
            int lt[200] = {0}; for (int i = 0; i < L; i++) for (int j = 0; j < q[i]; j++) lt[j]++;
            ll h = 1; for (int i = 0; i < L; i++) for (int j = 0; j < q[i]; j++) h = (i128)h * (q[i] - j + lt[j] - i - 1) % P;
            return (i128)fact[N] * inv(h) % P; };
        int rho[200]; for (int i = 0; i < m; i++) rho[i] = m - i;
        ll dr = dimmod(rho, m); ll lhs = 0;
        for (ll r = 0; r < PN[N]; r++) { int L = unrank(N, r, a); lhs = (lhs + (i128)G[r] * dimmod(a, L)) % P; }
        fprintf(fo, "# check sum_lambda g*dim(lambda) - dim(rho)^2 == 0 mod P : %s\n", (lhs == (i128)dr * dr % P) ? "OK" : "FAIL");
        // checksum 1b: sum_lambda g * chi^lambda(hhat) == chi^rho(hhat)^2 = 1  (mod P)
        { ll s2 = 0; for (ll r = 0; r < PN[N]; r++) { ll hv = H[r] % P; if (hv < 0) hv += P; s2 = (s2 + (i128)G[r] * hv) % P; }
          fprintf(fo, "# check sum_lambda g*a(lambda) == 1 mod P : %s\n", s2 == 1 ? "OK" : "FAIL"); }
        // checksum 2: g(rho,rho,(N)) = 1 and g(rho,rho,(1^N)) = 1 (rho self-conjugate)
        fprintf(fo, "# g(rho,rho,(N)) mod P = %lld ; g(rho,rho,(1^N)) mod P = %lld\n", G[0], G[PN[N] - 1]);
    }
    // checksum of g vector
    ull hsh = 1469598103934665603ULL; for (ll r = 0; r < PN[N]; r++) { hsh ^= (ull)G[r]; hsh *= 1099511628211ULL; }
    fprintf(fo, "# fnv1a(g mod P) = %016llx\n", hsh);
    fclose(fo);
    // full dump for small N
    if (N <= 45) {
        string fn = string(argv[5]) + ".full";
        FILE* ff = fopen(fn.c_str(), "w");
        for (ll r = 0; r < PN[N]; r++) { int L = unrank(N, r, a); for (int i = 0; i < L; i++) fprintf(ff, "%d%s", a[i], i + 1 < L ? "," : ""); fprintf(ff, " %lld %lld\n", G[r], H[r]); }
        fclose(ff);
    }
    fprintf(stderr, "done\n");
    return 0;
}
