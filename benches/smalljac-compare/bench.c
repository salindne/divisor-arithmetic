/* Micro-benchmark for smalljac's genus-2 Jacobian group law
 * (hecurve_g2_compose = ADD, hecurve_g2_square = DBL), to compare against the
 * Rust `divisor-arithmetic` cargo-bench numbers for g2::ramified::not_char2.
 *
 * Imaginary/ramified genus-2 model: y^2 = f(x), deg f = 5 (monic).
 * Times the affine path (ctx = NULL => one field inversion per op), which is
 * the same affine, one-inversion model the Rust explicit formulas use.
 *
 * Build: see README.md in this directory (requires an external smalljac +
 * ff_poly checkout; on arm64 the x86-64 inline asm must be replaced — the
 * README lists the two portable patches).
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "hecurve.h"
#include "cstd.h"

#define POOL 256

static double now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

static void run(unsigned long p, const char *label) {
    ff_setup_ui(p);

    /* monic degree-5 f with random coefficients f[0..4], f[5]=1, f[6]=0 */
    ff_t f[HECURVE_DEGREE + 1];
    for (int i = 0; i < HECURVE_DEGREE + 1; i++) _ff_set_zero(f[i]);
    for (int i = 0; i < 5; i++) _ff_random(f[i]);
    _ff_set_one(f[5]);

    /* pool of valid random weight-2 divisors */
    static ff_t U[POOL][HECURVE_GENUS + 1], V[POOL][HECURVE_GENUS];
    int n = 0;
    for (int tries = 0; tries < POOL * 50 && n < POOL; tries++) {
        ff_t u[HECURVE_GENUS + 1], v[HECURVE_GENUS];
        hecurve_random(u, v, f);
        if (_ff_zero(u[2]) || !_ff_one(u[2])) continue; /* keep monic deg-2 u */
        if (!hecurve_verify(u, v, f)) continue;
        for (int i = 0; i < HECURVE_GENUS + 1; i++) _ff_set(U[n][i], u[i]);
        for (int i = 0; i < HECURVE_GENUS; i++) _ff_set(V[n][i], v[i]);
        n++;
    }
    if (n < 8) { printf("%-12s p=%lu: only %d divisors, skipping\n", label, p, n); return; }

    ff_t ru[HECURVE_GENUS + 1], rv[HECURVE_GENUS];
    unsigned long sink = 0;
    long iters = 2000000;

    /* ---- ADD (compose distinct divisors) ---- */
    for (int i = 0; i < n - 1; i++) hecurve_g2_compose(ru, rv, U[i], V[i], U[i+1], V[i+1], f, 0);
    double t0 = now_ns();
    for (long k = 0; k < iters; k++) {
        int i = k % (n - 1);
        hecurve_g2_compose(ru, rv, U[i], V[i], U[i+1], V[i+1], f, 0);
        sink += _ff_get_ui(ru[0]);
    }
    double add_ns = (now_ns() - t0) / iters;

    /* ---- DBL (square a divisor) ---- */
    for (int i = 0; i < n; i++) hecurve_g2_square(ru, rv, U[i], V[i], f, 0);
    t0 = now_ns();
    for (long k = 0; k < iters; k++) {
        int i = k % n;
        hecurve_g2_square(ru, rv, U[i], V[i], f, 0);
        sink += _ff_get_ui(ru[0]);
    }
    double dbl_ns = (now_ns() - t0) / iters;

    printf("%-12s p=%-20lu add=%7.1f ns   double=%7.1f ns   (n=%d, sink=%lu)\n",
           label, p, add_ns, dbl_ns, n, sink);
}

int main(void) {
    run(8191UL, "13-bit");
    run(65521UL, "16-bit");
    run(1048583UL, "20-bit");
    run(1000000007UL, "30-bit");
    run((1UL << 31) - 1, "31-bit M31");
    /* 56-bit prime: matched-width comparison vs Rust PrimeField<72057594037927931> */
    run(72057594037927931UL, "56-bit Fp56");
    return 0;
}
