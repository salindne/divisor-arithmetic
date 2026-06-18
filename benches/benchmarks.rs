//! Benchmarks for genus 2 ramified divisor arithmetic operations.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use divisor_arithmetic::field::{BinaryExtField, Field, MontgomeryField, PrimeField};
use divisor_arithmetic::g2::ramified::{arbitrary, char2, not_char2};

// =============================================================================
// NCH2 Benchmarks (y² = f(x), not characteristic 2)
// =============================================================================

fn bench_not_char2_deg2_add<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();

    // Random curve constants
    let cc = not_char2::CurveConstants {
        f3: F::random(&mut rng),
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
    };

    // Random degree-2 divisors
    let d1 = not_char2::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );
    let d2 = not_char2::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );

    c.bench_with_input(
        BenchmarkId::new("not_char2_deg2_add", name),
        &(&d1, &d2, &cc),
        |b, (d1, d2, cc)| b.iter(|| not_char2::add(black_box(*d1), black_box(*d2), black_box(*cc))),
    );
}

fn bench_not_char2_deg2_dbl<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();

    let cc = not_char2::CurveConstants {
        f3: F::random(&mut rng),
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
    };

    let d = not_char2::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );

    c.bench_with_input(
        BenchmarkId::new("not_char2_deg2_dbl", name),
        &(&d, &cc),
        |b, (d, cc)| b.iter(|| not_char2::double(black_box(*d), black_box(*cc))),
    );
}

// =============================================================================
// ARB Benchmarks (y² + h(x)y = f(x), arbitrary characteristic)
// =============================================================================

fn bench_arbitrary_deg2_add<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();

    let cc = arbitrary::CurveConstants {
        f4: F::random(&mut rng),
        f3: F::random(&mut rng),
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
        h2: F::zero(),
        h1: F::random(&mut rng),
        h0: F::random(&mut rng),
    };

    let d1 = arbitrary::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );
    let d2 = arbitrary::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );

    c.bench_with_input(
        BenchmarkId::new("arbitrary_deg2_add", name),
        &(&d1, &d2, &cc),
        |b, (d1, d2, cc)| b.iter(|| arbitrary::add(black_box(*d1), black_box(*d2), black_box(*cc))),
    );
}

fn bench_arbitrary_deg2_dbl<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();

    let cc = arbitrary::CurveConstants {
        f4: F::random(&mut rng),
        f3: F::random(&mut rng),
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
        h2: F::zero(),
        h1: F::random(&mut rng),
        h0: F::random(&mut rng),
    };

    let d = arbitrary::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );

    c.bench_with_input(
        BenchmarkId::new("arbitrary_deg2_dbl", name),
        &(&d, &cc),
        |b, (d, cc)| b.iter(|| arbitrary::double(black_box(*d), black_box(*cc))),
    );
}

// =============================================================================
// CH2 Benchmarks (y² + h(x)y = f(x), characteristic 2)
// =============================================================================

fn bench_char2_deg2_add<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();

    let cc = char2::CurveConstants {
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
        h2: F::zero(),
        h1: F::random(&mut rng),
        h0: F::random(&mut rng),
    };

    let d1 = char2::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );
    let d2 = char2::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );

    c.bench_with_input(
        BenchmarkId::new("char2_deg2_add", name),
        &(&d1, &d2, &cc),
        |b, (d1, d2, cc)| b.iter(|| char2::add(black_box(*d1), black_box(*d2), black_box(*cc))),
    );
}

fn bench_char2_deg2_dbl<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();

    let cc = char2::CurveConstants {
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
        h2: F::zero(),
        h1: F::random(&mut rng),
        h0: F::random(&mut rng),
    };

    let d = char2::DivisorCoords::deg2(
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
        F::random(&mut rng),
    );

    c.bench_with_input(
        BenchmarkId::new("char2_deg2_dbl", name),
        &(&d, &cc),
        |b, (d, cc)| b.iter(|| char2::double(black_box(*d), black_box(*cc))),
    );
}

// =============================================================================
// SPLIT model benchmarks (two points at infinity, balanced divisors)
//
// Split divisors carry a balance weight and must be valid; we build a
// representative degree-2 pair via the generic Cantor reference, then time the
// explicit coordinate formulas on it.
// =============================================================================

use divisor_arithmetic::g2::split::{arbitrary as sarb, char2 as sch2, not_char2 as snch2};
use divisor_arithmetic::generic::split as gsplit;
use divisor_arithmetic::poly::Poly;

fn sqrt_trial<F: Field>(t: F, rng: &mut impl rand::Rng) -> Option<F> {
    if t.is_zero() {
        return Some(F::zero());
    }
    for _ in 0..300_000 {
        let c = F::random(rng);
        if c.square() == t {
            return Some(c);
        }
    }
    None
}

/// Build a valid degree-1 divisor at point `(a,b)` then compose two of them into
/// a degree-2 divisor via the generic oracle (negative basis).
fn build_deg2<F: Field>(
    f: &Poly<F>,
    h: &Poly<F>,
    vn: &Poly<F>,
    pts: &[(F, F)],
) -> gsplit::Divisor<F> {
    let deg1 = |a: F, b: F| {
        let u = Poly::from_coeffs(vec![-a, F::one()]);
        let r = (vn - &Poly::constant(b)).rem(&u);
        let vhat = vn - &r;
        let w = (f - &(&vhat * &(&vhat + h))).exact_div(&u);
        gsplit::Divisor::new(u, vhat, w, 0)
    };
    gsplit::add_neg(
        &deg1(pts[0].0, pts[0].1),
        &deg1(pts[1].0, pts[1].1),
        f,
        h,
        vn,
        2,
    )
}

// ---- nch2 split (h = 0) ----
fn bench_split_nch2<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();
    let (cc, d1, d2) = loop {
        let cc = snch2::precompute(
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
        );
        let (f, vn) = (cc.f_poly(), cc.vn());
        let h = Poly::zero();
        let mut pts = Vec::new();
        for _ in 0..400 {
            let a = F::random(&mut rng);
            if let Some(b) = sqrt_trial(f.eval(a), &mut rng) {
                pts.push((a, b));
                if pts.len() == 4 {
                    break;
                }
            }
        }
        if pts.len() < 4 {
            continue;
        }
        let g1 = build_deg2(&f, &h, &vn, &pts[0..2]);
        let g2 = build_deg2(&f, &h, &vn, &pts[2..4]);
        let to = |d: &gsplit::Divisor<F>| match d.u.deg() {
            2 => snch2::DivisorCoords::deg2(
                d.u.coeff(1),
                d.u.coeff(0),
                d.v.coeff(1),
                d.v.coeff(0),
                d.n,
            ),
            _ => snch2::DivisorCoords::deg1(d.u.coeff(0), d.v.coeff(1), d.v.coeff(0), d.n),
        };
        if g1.u.deg() == 2 && g2.u.deg() == 2 {
            break (cc, to(&g1), to(&g2));
        }
    };
    c.bench_with_input(
        BenchmarkId::new("split_nch2_add", name),
        &(d1, d2, cc),
        |b, (d1, d2, cc)| b.iter(|| snch2::add_neg(black_box(d1), black_box(d2), black_box(cc))),
    );
    c.bench_with_input(
        BenchmarkId::new("split_nch2_dbl", name),
        &(d1, cc),
        |b, (d1, cc)| b.iter(|| snch2::double_neg(black_box(d1), black_box(cc))),
    );
}

// ---- arbitrary characteristic split (h != 0) ----
fn bench_split_arb<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();
    let two_inv = (F::one() + F::one()).inv();
    let (cc, d1, d2) = loop {
        let h = [
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
        ];
        let y3 = F::random(&mut rng);
        let c3 = y3.double() + h[3];
        if c3.is_zero() {
            continue;
        }
        let f6 = y3 * y3 + h[3] * y3;
        if f6.is_zero() {
            continue;
        }
        let f = [
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            f6,
        ];
        let cc = sarb::precompute(f, h, y3);
        let (fp, hp, vn) = (cc.f_poly(), cc.h_poly(), cc.vn());
        let mut pts = Vec::new();
        for _ in 0..400 {
            let a = F::random(&mut rng);
            let (ha, fa) = (hp.eval(a), fp.eval(a));
            if let Some(s) = sqrt_trial(ha * ha + (fa + fa).double(), &mut rng) {
                pts.push((a, (s - ha) * two_inv));
                if pts.len() == 4 {
                    break;
                }
            }
        }
        if pts.len() < 4 {
            continue;
        }
        let g1 = build_deg2(&fp, &hp, &vn, &pts[0..2]);
        let g2 = build_deg2(&fp, &hp, &vn, &pts[2..4]);
        let to = |d: &gsplit::Divisor<F>| match d.u.deg() {
            2 => sarb::DivisorCoords::deg2(
                d.u.coeff(1),
                d.u.coeff(0),
                d.v.coeff(1),
                d.v.coeff(0),
                d.n,
            ),
            _ => sarb::DivisorCoords::deg1(d.u.coeff(0), d.v.coeff(1), d.v.coeff(0), d.n),
        };
        if g1.u.deg() == 2 && g2.u.deg() == 2 {
            break (cc, to(&g1), to(&g2));
        }
    };
    c.bench_with_input(
        BenchmarkId::new("split_arb_add", name),
        &(d1, d2, cc),
        |b, (d1, d2, cc)| b.iter(|| sarb::add_neg(black_box(d1), black_box(d2), black_box(cc))),
    );
    c.bench_with_input(
        BenchmarkId::new("split_arb_dbl", name),
        &(d1, cc),
        |b, (d1, cc)| b.iter(|| sarb::double_neg(black_box(d1), black_box(cc))),
    );
}

// ---- characteristic 2 split ----
fn bench_split_char2<const K: usize>(c: &mut Criterion, name: &str) {
    type GF<const K: usize> = BinaryExtField<K>;
    let mut rng = rand::thread_rng();
    let solve = |cc: &sch2::CurveConstants<GF<K>>, a: GF<K>| -> Option<GF<K>> {
        let (ha, fa) = (cc.h_poly().eval(a), cc.f_poly().eval(a));
        (0..(1u64 << K))
            .map(GF::<K>::new)
            .find(|&b| b * b + ha * b == fa)
    };
    let (cc, d1, d2) = loop {
        let beta = GF::<K>::random(&mut rng);
        let f6 = beta.square() + beta;
        if f6.is_zero() {
            continue;
        }
        let cc = sch2::precompute(
            GF::<K>::random(&mut rng),
            GF::<K>::random(&mut rng),
            GF::<K>::random(&mut rng),
            f6,
            GF::<K>::random(&mut rng),
            GF::<K>::random(&mut rng),
            beta,
        );
        let (fp, hp, vn) = (cc.f_poly(), cc.h_poly(), cc.vn());
        let mut pts = Vec::new();
        for ai in 0..(1u64 << K) {
            let a = GF::<K>::new(ai);
            if let Some(b) = solve(&cc, a) {
                pts.push((a, b));
                if pts.len() == 4 {
                    break;
                }
            }
        }
        if pts.len() < 4 {
            continue;
        }
        let g1 = build_deg2(&fp, &hp, &vn, &pts[0..2]);
        let g2 = build_deg2(&fp, &hp, &vn, &pts[2..4]);
        let to = |d: &gsplit::Divisor<GF<K>>| match d.u.deg() {
            2 => sch2::DivisorCoords::deg2(
                d.u.coeff(1),
                d.u.coeff(0),
                d.v.coeff(1),
                d.v.coeff(0),
                d.n,
            ),
            _ => sch2::DivisorCoords::deg1(d.u.coeff(0), d.v.coeff(1), d.v.coeff(0), d.n),
        };
        if g1.u.deg() == 2 && g2.u.deg() == 2 {
            break (cc, to(&g1), to(&g2));
        }
    };
    c.bench_with_input(
        BenchmarkId::new("split_char2_add", name),
        &(d1, d2, cc),
        |b, (d1, d2, cc)| b.iter(|| sch2::add_neg(black_box(d1), black_box(d2), black_box(cc))),
    );
    c.bench_with_input(
        BenchmarkId::new("split_char2_dbl", name),
        &(d1, cc),
        |b, (d1, cc)| b.iter(|| sch2::double_neg(black_box(d1), black_box(cc))),
    );
}

fn split_benchmarks(c: &mut Criterion) {
    bench_split_nch2::<PrimeField<8191>>(c, "F8191");
    bench_split_nch2::<PrimeField<65521>>(c, "F65521");
    bench_split_arb::<PrimeField<8191>>(c, "F8191");
    bench_split_arb::<PrimeField<65521>>(c, "F65521");
    bench_split_char2::<8>(c, "GF256");
    bench_split_char2::<16>(c, "GF65536");
}

// =============================================================================
// Field operation benchmarks
// =============================================================================

fn bench_field_ops<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();
    let a = F::random(&mut rng);
    let b = F::random(&mut rng);

    c.bench_with_input(
        BenchmarkId::new("field_mul", name),
        &(a, b),
        |bench, (a, b)| bench.iter(|| black_box(*a) * black_box(*b)),
    );

    c.bench_with_input(BenchmarkId::new("field_inv", name), &a, |bench, a| {
        bench.iter(|| black_box(*a).inv())
    });

    c.bench_with_input(BenchmarkId::new("field_square", name), &a, |bench, a| {
        bench.iter(|| black_box(*a).square())
    });
}

// =============================================================================
// Benchmark groups
// =============================================================================

fn not_char2_benchmarks(c: &mut Criterion) {
    bench_not_char2_deg2_add::<PrimeField<7>>(c, "F7");
    bench_not_char2_deg2_add::<PrimeField<8191>>(c, "F8191");
    bench_not_char2_deg2_add::<PrimeField<65521>>(c, "F65521");
    // ~56-bit prime: matched-width comparison vs smalljac (built for 57-bit primes).
    bench_not_char2_deg2_add::<PrimeField<72057594037927931>>(c, "Fp56");
    // Montgomery field at the same 56-bit prime (the "chase smalljac" field).
    bench_not_char2_deg2_add::<MontgomeryField<72057594037927931>>(c, "Mont56");

    bench_not_char2_deg2_dbl::<PrimeField<7>>(c, "F7");
    bench_not_char2_deg2_dbl::<PrimeField<8191>>(c, "F8191");
    bench_not_char2_deg2_dbl::<PrimeField<65521>>(c, "F65521");
    bench_not_char2_deg2_dbl::<PrimeField<72057594037927931>>(c, "Fp56");
    bench_not_char2_deg2_dbl::<MontgomeryField<72057594037927931>>(c, "Mont56");
}

fn arbitrary_benchmarks(c: &mut Criterion) {
    bench_arbitrary_deg2_add::<PrimeField<7>>(c, "F7");
    bench_arbitrary_deg2_add::<PrimeField<8191>>(c, "F8191");
    bench_arbitrary_deg2_add::<PrimeField<65521>>(c, "F65521");

    bench_arbitrary_deg2_dbl::<PrimeField<7>>(c, "F7");
    bench_arbitrary_deg2_dbl::<PrimeField<8191>>(c, "F8191");
    bench_arbitrary_deg2_dbl::<PrimeField<65521>>(c, "F65521");
}

fn char2_benchmarks(c: &mut Criterion) {
    bench_char2_deg2_add::<BinaryExtField<8>>(c, "GF256");
    bench_char2_deg2_add::<BinaryExtField<16>>(c, "GF65536");

    bench_char2_deg2_dbl::<BinaryExtField<8>>(c, "GF256");
    bench_char2_deg2_dbl::<BinaryExtField<16>>(c, "GF65536");
}

fn field_benchmarks(c: &mut Criterion) {
    bench_field_ops::<PrimeField<65521>>(c, "F65521");
    bench_field_ops::<PrimeField<72057594037927931>>(c, "Fp56");
    bench_field_ops::<MontgomeryField<72057594037927931>>(c, "Mont56");
    bench_field_ops::<BinaryExtField<8>>(c, "GF256");
    bench_field_ops::<BinaryExtField<16>>(c, "GF65536");
}

// =============================================================================
// Batched group law (Montgomery simultaneous inversion)
// =============================================================================

fn rand_deg2<F: Field>(rng: &mut impl rand::Rng) -> not_char2::DivisorCoords<F> {
    not_char2::DivisorCoords::deg2(
        F::random(rng),
        F::random(rng),
        F::random(rng),
        F::random(rng),
    )
}

/// Amortized throughput of the batched group law: one field inversion is shared
/// across a batch of `n` independent ops via `batch_invert`. This is the metric
/// that matters for smalljac-style generic-group algorithms (and matches its
/// `ctx` + `ff_parallel_invert` path). Criterion reports elements/sec; per-op
/// time = batch_time / n.
fn bench_not_char2_batched<F: Field>(c: &mut Criterion, name: &str, n: usize) {
    let mut rng = rand::thread_rng();
    let cc = not_char2::CurveConstants {
        f3: F::random(&mut rng),
        f2: F::random(&mut rng),
        f1: F::random(&mut rng),
        f0: F::random(&mut rng),
    };
    let pairs: Vec<_> = (0..n)
        .map(|_| (rand_deg2::<F>(&mut rng), rand_deg2::<F>(&mut rng)))
        .collect();
    let singles: Vec<_> = (0..n).map(|_| rand_deg2::<F>(&mut rng)).collect();

    let mut g = c.benchmark_group("not_char2_batched");
    g.throughput(Throughput::Elements(n as u64));
    g.bench_with_input(
        BenchmarkId::new("add_batch", format!("{name}/N={n}")),
        &pairs,
        |b, pairs| b.iter(|| not_char2::add_batch(black_box(pairs), black_box(&cc))),
    );
    g.bench_with_input(
        BenchmarkId::new("dbl_batch", format!("{name}/N={n}")),
        &singles,
        |b, singles| b.iter(|| not_char2::double_batch(black_box(singles), black_box(&cc))),
    );
    g.finish();
}

fn batched_benchmarks(c: &mut Criterion) {
    // Matched 56-bit width vs smalljac's batched (ctx + ff_parallel_invert) path,
    // plus 16-bit for the field-width contrast.
    bench_not_char2_batched::<PrimeField<72057594037927931>>(c, "Fp56", 1024);
    bench_not_char2_batched::<MontgomeryField<72057594037927931>>(c, "Mont56", 1024);
    bench_not_char2_batched::<PrimeField<65521>>(c, "F65521", 1024);
}

criterion_group!(
    benches,
    field_benchmarks,
    not_char2_benchmarks,
    arbitrary_benchmarks,
    char2_benchmarks,
    split_benchmarks,
    batched_benchmarks,
);
criterion_main!(benches);
