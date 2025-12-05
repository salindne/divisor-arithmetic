//! Benchmarks for genus 2 ramified divisor arithmetic operations.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use divisor_arithmetic::field::{BinaryExtField, Field, PrimeField};
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
    
    c.bench_with_input(BenchmarkId::new("not_char2_deg2_add", name), &(&d1, &d2, &cc), |b, (d1, d2, cc)| {
        b.iter(|| not_char2::add(black_box(*d1), black_box(*d2), black_box(*cc)))
    });
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
    
    c.bench_with_input(BenchmarkId::new("not_char2_deg2_dbl", name), &(&d, &cc), |b, (d, cc)| {
        b.iter(|| not_char2::double(black_box(*d), black_box(*cc)))
    });
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
    
    c.bench_with_input(BenchmarkId::new("arbitrary_deg2_add", name), &(&d1, &d2, &cc), |b, (d1, d2, cc)| {
        b.iter(|| arbitrary::add(black_box(*d1), black_box(*d2), black_box(*cc)))
    });
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
    
    c.bench_with_input(BenchmarkId::new("arbitrary_deg2_dbl", name), &(&d, &cc), |b, (d, cc)| {
        b.iter(|| arbitrary::double(black_box(*d), black_box(*cc)))
    });
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
    
    c.bench_with_input(BenchmarkId::new("char2_deg2_add", name), &(&d1, &d2, &cc), |b, (d1, d2, cc)| {
        b.iter(|| char2::add(black_box(*d1), black_box(*d2), black_box(*cc)))
    });
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
    
    c.bench_with_input(BenchmarkId::new("char2_deg2_dbl", name), &(&d, &cc), |b, (d, cc)| {
        b.iter(|| char2::double(black_box(*d), black_box(*cc)))
    });
}

// =============================================================================
// Field operation benchmarks
// =============================================================================

fn bench_field_ops<F: Field>(c: &mut Criterion, name: &str) {
    let mut rng = rand::thread_rng();
    let a = F::random(&mut rng);
    let b = F::random(&mut rng);
    
    c.bench_with_input(BenchmarkId::new("field_mul", name), &(a, b), |bench, (a, b)| {
        bench.iter(|| black_box(*a) * black_box(*b))
    });
    
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
    
    bench_not_char2_deg2_dbl::<PrimeField<7>>(c, "F7");
    bench_not_char2_deg2_dbl::<PrimeField<8191>>(c, "F8191");
    bench_not_char2_deg2_dbl::<PrimeField<65521>>(c, "F65521");
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
    bench_field_ops::<BinaryExtField<8>>(c, "GF256");
    bench_field_ops::<BinaryExtField<16>>(c, "GF65536");
}

criterion_group!(
    benches,
    field_benchmarks,
    not_char2_benchmarks,
    arbitrary_benchmarks,
    char2_benchmarks,
);
criterion_main!(benches);
