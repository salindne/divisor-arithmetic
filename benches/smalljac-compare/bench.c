/* Micro-benchmark for smalljac's genus-2 Jacobian group law, scalar and
 * batched, to compare against the Rust `divisor-arithmetic` cargo-bench numbers
 * for g2::ramified::not_char2. See ../../README.md ("smalljac comparison").
 *
 * Imaginary/ramified genus-2 model: y^2 = f(x), deg f = 5 (monic, depressed:
 * f[4] = 0 — REQUIRED, else hecurve reverts to the slow generic Cantor path and
 * the comparison is meaningless; the Rust not_char2 model also uses f4 = 0).
 *
 *  - scalar:  ctx = NULL  => one field inversion per op.
 *  - batched: ctx state machine + ff_parallel_invert => one inversion per batch
 *             (Montgomery's trick) — the throughput metric that matters for the
 *             generic-group algorithms smalljac is built for.
 *
 * Build: see README.md in this directory (needs an external smalljac + ff_poly
 * checkout; on arm64 apply the two portability patches listed there).
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "hecurve.h"
#include "cstd.h"

#define POOL  1100
#define BATCH 1024

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

static ff_t U[POOL][HECURVE_GENUS + 1], V[POOL][HECURVE_GENUS];
static ff_t f[HECURVE_DEGREE + 1];
static hecurve_ctx_t ctx[BATCH];
static ff_t RU[BATCH][HECURVE_GENUS + 1], RV[BATCH][HECURVE_GENUS];
static ff_t invs[BATCH];
static int idx[BATCH];

static int build_pool(unsigned long p) {
    ff_setup_ui(p);
    for (int i = 0; i < HECURVE_DEGREE + 1; i++) _ff_set_zero(f[i]);
    for (int i = 0; i < 4; i++) _ff_random(f[i]); /* f0..f3 random, f4 = 0 */
    _ff_set_one(f[5]);                            /* monic degree 5 */
    int n = 0;
    for (int tries = 0; tries < POOL * 50 && n < POOL; tries++) {
        ff_t u[HECURVE_GENUS + 1], v[HECURVE_GENUS];
        hecurve_random(u, v, f);
        if (_ff_zero(u[2]) || !_ff_one(u[2])) continue;
        if (!hecurve_verify(u, v, f)) continue;
        for (int i = 0; i < HECURVE_GENUS + 1; i++) _ff_set(U[n][i], u[i]);
        for (int i = 0; i < HECURVE_GENUS; i++) _ff_set(V[n][i], v[i]);
        n++;
    }
    return n;
}

static void run_fieldops(unsigned long p, const char *label) {
    ff_setup_ui(p);
    ff_t a, b, acc, t;
    _ff_set_ui(a, 12345 % p); _ff_set_ui(b, 67891 % p); _ff_set(acc, a);
    long iters = 50000000;
    double t0 = now_ns();
    for (long k = 0; k < iters; k++) ff_mult(acc, acc, b);
    double m_ns = (now_ns() - t0) / iters;
    unsigned long sink = _ff_get_ui(acc);
    _ff_set_ui(acc, 12345 % p);
    long it2 = 5000000;
    t0 = now_ns();
    for (long k = 0; k < it2; k++) { ff_invert(t, acc); _ff_mult(acc, t, b); }
    double i_ns = (now_ns() - t0) / it2;
    sink += _ff_get_ui(acc);
    static ff_t xs[BATCH], zs[BATCH];
    for (int i = 0; i < BATCH; i++) _ff_random(xs[i]);
    long it3 = 20000;
    t0 = now_ns();
    for (long k = 0; k < it3; k++) ff_parallel_invert(zs, xs, BATCH);
    double pi_ns = (now_ns() - t0) / ((double)it3 * BATCH);
    sink += _ff_get_ui(zs[0]);
    printf("fieldop %-12s p=%-20lu  M=%5.2f ns  I=%6.2f ns  batchedI=%5.2f ns/elem  (I/M=%.1f)\n",
           label, p, m_ns, i_ns, pi_ns, i_ns / m_ns);
    (void)sink;
}

static void run_scalar(unsigned long p, const char *label) {
    int n = build_pool(p);
    if (n < 8) { printf("scalar  %-12s p=%lu: only %d divisors\n", label, p, n); return; }
    ff_t ru[HECURVE_GENUS + 1], rv[HECURVE_GENUS];
    unsigned long sink = 0;
    long iters = 2000000;
    for (int i = 0; i < n - 1; i++) hecurve_g2_compose(ru, rv, U[i], V[i], U[i+1], V[i+1], f, 0);
    double t0 = now_ns();
    for (long k = 0; k < iters; k++) { int i = k % (n - 1);
        hecurve_g2_compose(ru, rv, U[i], V[i], U[i+1], V[i+1], f, 0); sink += _ff_get_ui(ru[0]); }
    double add_ns = (now_ns() - t0) / iters;
    for (int i = 0; i < n; i++) hecurve_g2_square(ru, rv, U[i], V[i], f, 0);
    t0 = now_ns();
    for (long k = 0; k < iters; k++) { int i = k % n;
        hecurve_g2_square(ru, rv, U[i], V[i], f, 0); sink += _ff_get_ui(ru[0]); }
    double dbl_ns = (now_ns() - t0) / iters;
    printf("scalar  %-12s p=%-20lu add=%7.1f ns  double=%7.1f ns\n", label, p, add_ns, dbl_ns);
    (void)sink;
}

static void run_batched(unsigned long p, const char *label) {
    int n = build_pool(p);
    if (n < BATCH + 1) { printf("batched %-12s p=%lu: only %d divisors\n", label, p, n); return; }
    unsigned long sink = 0;
    long rounds = 2000;

    double t0 = now_ns();
    for (long r = 0; r < rounds; r++) {
        int m = 0;
        for (int i = 0; i < BATCH; i++) {
            ctx[i].state = 0;
            int done = hecurve_g2_compose(RU[i], RV[i], U[i], V[i], U[i+1], V[i+1], f, &ctx[i]);
            if (!done) { invs[m] = ctx[i].invert; idx[m] = i; m++; }
        }
        ff_parallel_invert(invs, invs, m);
        for (int j = 0; j < m; j++) {
            int i = idx[j];
            _ff_set(ctx[i].invert, invs[j]);
            hecurve_g2_compose(RU[i], RV[i], U[i], V[i], U[i+1], V[i+1], f, &ctx[i]);
        }
        sink += _ff_get_ui(RU[0][0]);
    }
    double add_ns = (now_ns() - t0) / ((double)rounds * BATCH);

    t0 = now_ns();
    for (long r = 0; r < rounds; r++) {
        int m = 0;
        for (int i = 0; i < BATCH; i++) {
            ctx[i].state = 0;
            int done = hecurve_g2_square(RU[i], RV[i], U[i], V[i], f, &ctx[i]);
            if (!done) { invs[m] = ctx[i].invert; idx[m] = i; m++; }
        }
        ff_parallel_invert(invs, invs, m);
        for (int j = 0; j < m; j++) {
            int i = idx[j];
            _ff_set(ctx[i].invert, invs[j]);
            hecurve_g2_square(RU[i], RV[i], U[i], V[i], f, &ctx[i]);
        }
        sink += _ff_get_ui(RU[0][0]);
    }
    double dbl_ns = (now_ns() - t0) / ((double)rounds * BATCH);
    printf("batched %-12s p=%-20lu add=%7.1f ns  double=%7.1f ns  (N=%d)\n",
           label, p, add_ns, dbl_ns, BATCH);
    (void)sink;
}

int main(void) {
    run_fieldops(65521UL, "16-bit");
    run_fieldops(72057594037927931UL, "56-bit Fp56");
    run_scalar(65521UL, "16-bit");
    run_batched(65521UL, "16-bit");
    run_scalar(72057594037927931UL, "56-bit Fp56");
    run_batched(72057594037927931UL, "56-bit Fp56");
    return 0;
}
