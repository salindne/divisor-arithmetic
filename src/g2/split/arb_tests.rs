//! Tests for the arbitrary-characteristic split formulas, over odd prime fields
//! with non-zero `h` (the new generality vs nch2). Explicit formulas are checked
//! against the now-h-aware generic Cantor oracle.

use crate::field::{Field, PrimeField};
use crate::generic::split::{self, Divisor};
use crate::poly::Poly;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use super::arbitrary::{
    add_neg, add_pos, double_neg, double_pos, precompute, CurveConstants, DivisorCoords, G,
};
use super::test_support::{from_generic, try_sqrt};

// ---------------------------------------------------------------------------
// Harness (arb: h ≠ 0)
// ---------------------------------------------------------------------------

fn derivative<F: Field>(p: &Poly<F>) -> Poly<F> {
    if p.deg() < 1 {
        return Poly::zero();
    }
    let mut c = Vec::with_capacity(p.coeffs.len() - 1);
    let mut i_fe = F::zero();
    for i in 1..p.coeffs.len() {
        i_fe += F::one();
        c.push(p.coeffs[i] * i_fe);
    }
    Poly::from_coeffs(c)
}

/// Random valid arbitrary-characteristic split curve over an odd field, with a
/// non-zero `h` and a known point-at-infinity root `y3`. Nonsingular ⇔ the
/// completed-square polynomial `4f + h²` is squarefree of degree 6.
fn random_curve_arb<F: Field, R: Rng>(rng: &mut R) -> CurveConstants<F> {
    let two = F::one() + F::one();
    loop {
        let h = [F::random(rng), F::random(rng), F::random(rng), F::random(rng)];
        let y3 = F::random(rng);
        let c3 = y3.double() + h[3]; // 2y3 + h3 must be nonzero (distinct ∞ points)
        if c3.is_zero() {
            continue;
        }
        let f6 = y3 * y3 + h[3] * y3; // makes y3 a root of x² + h3 x − f6
        if f6.is_zero() {
            continue;
        }
        let f = [
            F::random(rng), F::random(rng), F::random(rng), F::random(rng),
            F::random(rng), F::random(rng), f6,
        ];
        let cc = precompute(f, h, y3);

        // nonsingularity: D = 4f + h² squarefree, deg 6.
        let fp = cc.f_poly();
        let hp = cc.h_poly();
        let d_poly = &(&fp * &Poly::constant(two * two)) + &(&hp * &hp);
        if d_poly.deg() != 6 {
            continue;
        }
        let dp = derivative(&d_poly);
        if dp.is_zero() || d_poly.xgcd(&dp).0.deg() != 0 {
            continue;
        }
        return cc;
    }
}

/// Convert explicit coords → generic divisor with `w = (f − v(v+h))/u`.
fn to_generic_arb<F: Field>(d: &DivisorCoords<F>, cc: &CurveConstants<F>, basis: &Poly<F>) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let u = match d.degree() {
        0 => Poly::constant(F::one()),
        1 => Poly::from_coeffs(vec![d.u0, F::one()]),
        2 => Poly::from_coeffs(vec![d.u0, d.u1, F::one()]),
        _ => unreachable!(),
    };
    let v = Poly::from_coeffs(vec![d.v0, d.v1, basis.coeff(2), basis.coeff(3)]);
    let vh = &v + &h;
    let w = (&f - &(&v * &vh)).exact_div(&u);
    Divisor::new(u, v, w, d.n)
}

fn neutral_neg_arb<F: Field>(cc: &CurveConstants<F>) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vn = cc.vn();
    let w = &f - &(&vn * &(&vn + &h));
    Divisor::new(Poly::constant(F::one()), vn, w, 1)
}

fn deg0_neg_arb<F: Field>(cc: &CurveConstants<F>, n: i32) -> Divisor<F> {
    let mut d = neutral_neg_arb(cc);
    d.n = n;
    d
}

/// Solve `b² + h(a)·b = f(a)` for `b` (odd characteristic). `None` if no root.
fn find_b<F: Field, R: Rng>(cc: &CurveConstants<F>, a: F, rng: &mut R) -> Option<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let two_inv = (F::one() + F::one()).inv();
    let ha = h.eval(a);
    let fa = f.eval(a);
    let disc = ha * ha + (fa + fa).double(); // h(a)² + 4 f(a)
    let s = try_sqrt(disc, rng)?;
    Some((s - ha) * two_inv)
}

fn deg1_neg_arb<F: Field>(cc: &CurveConstants<F>, a: F, b: F, n: i32) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vn = cc.vn();
    let u = Poly::from_coeffs(vec![-a, F::one()]);
    let r = (&vn - &Poly::constant(b)).rem(&u);
    let vhat = &vn - &r;
    let w = (&f - &(&vhat * &(&vhat + &h))).exact_div(&u);
    Divisor::new(u, vhat, w, n)
}

fn random_point_arb<F: Field, R: Rng>(cc: &CurveConstants<F>, rng: &mut R) -> Option<(F, F)> {
    for _ in 0..200 {
        let a = F::random(rng);
        if let Some(b) = find_b(cc, a, rng) {
            return Some((a, b));
        }
    }
    None
}

fn random_valid_neg_arb<F: Field, R: Rng>(cc: &CurveConstants<F>, rng: &mut R) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vn = cc.vn();
    let mut d = neutral_neg_arb(cc);
    let mut produced = 0usize;
    let mut guard = 0usize;
    while produced < G {
        guard += 1;
        if guard > 100_000 {
            break;
        }
        let r1 = F::random(rng);
        let r0 = F::random(rng);
        if r1.is_zero() {
            let base = deg0_neg_arb(cc, rng.gen_range(0..=G as i32));
            d = split::add_neg(&d, &base, &f, &h, &vn, G);
            produced += 1;
        } else {
            let u0 = r0 * r1.inv();
            let a = -u0;
            if let Some(b) = find_b(cc, a, rng) {
                let base = deg1_neg_arb(cc, a, b, rng.gen_range(0..=(G as i32 - 1)));
                d = split::add_neg(&d, &base, &f, &h, &vn, G);
                produced += 1;
            }
        }
    }
    d
}

fn check_double_neg_arb<F: Field + std::fmt::Display>(d: &Divisor<F>, cc: &CurveConstants<F>) {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vn = cc.vn();
    let got = double_neg(&from_generic(d), cc);
    let expected = from_generic(&split::double_neg(d, &f, &h, &vn, G));
    assert_eq!(got, expected, "arb double_neg mismatch for {:?}", from_generic(d));
}

fn check_add_neg_arb<F: Field + std::fmt::Display>(
    d1: &Divisor<F>,
    d2: &Divisor<F>,
    cc: &CurveConstants<F>,
) {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vn = cc.vn();
    let got = add_neg(&from_generic(d1), &from_generic(d2), cc);
    let expected = from_generic(&split::add_neg(d1, d2, &f, &h, &vn, G));
    assert_eq!(
        got, expected,
        "arb add_neg mismatch for {:?} + {:?}",
        from_generic(d1),
        from_generic(d2)
    );
}

fn dbl_neg_arb<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve_arb::<PrimeField<P>, _>(&mut rng);
        let vn = cc.vn();
        for _ in 0..divisors {
            let d = random_valid_neg_arb(&cc, &mut rng);
            // round-trip incl. n
            assert_eq!(to_generic_arb(&from_generic(&d), &cc, &vn), d, "arb neg round-trip");
            check_double_neg_arb(&d, &cc);
        }
        for n in 0..=2 {
            check_double_neg_arb(&deg0_neg_arb(&cc, n), &cc);
        }
        for _ in 0..8 {
            if let Some((a, b)) = random_point_arb(&cc, &mut rng) {
                check_double_neg_arb(&deg1_neg_arb(&cc, a, b, 0), &cc);
                check_double_neg_arb(&deg1_neg_arb(&cc, a, b, 1), &cc);
            }
        }
    }
}

#[test]
fn arb_dbl_neg_f7() {
    dbl_neg_arb::<7>(7001, 12, 60);
}

#[test]
fn arb_dbl_neg_f31() {
    dbl_neg_arb::<31>(7002, 12, 60);
}

#[test]
fn arb_dbl_neg_f127() {
    dbl_neg_arb::<127>(7003, 8, 40);
}

/// All 18 arb-neg DBL branches (`ADBL00`..`ADBL17`) are exercised.
#[test]
fn arb_dbl_neg_branch_coverage() {
    dbl_neg_arb::<5>(7101, 60, 120);
    dbl_neg_arb::<7>(7102, 60, 120);
    dbl_neg_arb::<11>(7103, 40, 120);
    dbl_neg_arb::<31>(7104, 30, 120);
    let labels: Vec<String> = (0..=17).map(|i| format!("ADBL{i:02}")).collect();
    let refs: Vec<&str> = labels.iter().map(|s| s.as_str()).collect();
    let missing = super::coverage::missing(&refs);
    assert!(missing.is_empty(), "arb-neg DBL branches never exercised: {:?}", missing);
}

// ---------------------------------------------------------------------------
// Addition (arb, neg basis)
// ---------------------------------------------------------------------------

fn add_neg_arb<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve_arb::<PrimeField<P>, _>(&mut rng);
        for _ in 0..divisors {
            let d1 = random_valid_neg_arb(&cc, &mut rng);
            let d2 = random_valid_neg_arb(&cc, &mut rng);
            // doubling always valid; reuses double_neg, not the explicit ADD.
            check_double_neg_arb(&d1, &cc);
            // adding a divisor to itself is out of the explicit ADD's domain
            // (it divides by zero); the dispatcher routes those to double_neg.
            if from_generic(&d1) != from_generic(&d2) {
                check_add_neg_arb(&d1, &d2, &cc);
            }
        }
    }
}

#[test]
fn arb_add_neg_f7() {
    add_neg_arb::<7>(7201, 8, 400);
}

#[test]
fn arb_add_neg_f31() {
    add_neg_arb::<31>(7202, 8, 400);
}

#[test]
fn arb_add_neg_f127() {
    add_neg_arb::<127>(7203, 6, 200);
}

/// All 59 arb-neg ADD branches (`AADD00`..`AADD58`) are exercised by random
/// testing (the deg-0/deg-1 neutral short-circuits included, since
/// `random_valid_neg_arb` produces divisors of every degree and weight).
#[test]
fn arb_add_neg_branch_coverage() {
    add_neg_arb::<5>(7301, 60, 120);
    add_neg_arb::<7>(7302, 60, 120);
    add_neg_arb::<11>(7303, 40, 120);
    add_neg_arb::<31>(7304, 30, 120);
    let labels: Vec<String> = (0..=58).map(|i| format!("AADD{i:02}")).collect();
    let refs: Vec<&str> = labels.iter().map(|s| s.as_str()).collect();
    let missing = super::coverage::missing(&refs);
    assert!(missing.is_empty(), "arb-neg ADD branches never exercised: {:?}", missing);
}

// ---------------------------------------------------------------------------
// Harness (arb, POSITIVE basis: uses Vpl and the pos oracle)
// ---------------------------------------------------------------------------

fn neutral_pos_arb<F: Field>(cc: &CurveConstants<F>) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vpl = cc.vpl();
    let w = &f - &(&vpl * &(&vpl + &h));
    Divisor::new(Poly::constant(F::one()), vpl, w, 1)
}

fn deg0_pos_arb<F: Field>(cc: &CurveConstants<F>, n: i32) -> Divisor<F> {
    let mut d = neutral_pos_arb(cc);
    d.n = n;
    d
}

fn deg1_pos_arb<F: Field>(cc: &CurveConstants<F>, a: F, b: F, n: i32) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vpl = cc.vpl();
    let u = Poly::from_coeffs(vec![-a, F::one()]);
    let r = (&vpl - &Poly::constant(b)).rem(&u);
    let vhat = &vpl - &r;
    let w = (&f - &(&vhat * &(&vhat + &h))).exact_div(&u);
    Divisor::new(u, vhat, w, n)
}

fn random_valid_pos_arb<F: Field, R: Rng>(cc: &CurveConstants<F>, rng: &mut R) -> Divisor<F> {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vpl = cc.vpl();
    let mut d = neutral_pos_arb(cc);
    let mut produced = 0usize;
    let mut guard = 0usize;
    while produced < G {
        guard += 1;
        if guard > 100_000 {
            break;
        }
        let r1 = F::random(rng);
        let r0 = F::random(rng);
        if r1.is_zero() {
            let base = deg0_pos_arb(cc, rng.gen_range(0..=G as i32));
            d = split::add_pos(&d, &base, &f, &h, &vpl, G);
            produced += 1;
        } else {
            let u0 = r0 * r1.inv();
            let a = -u0;
            if let Some(b) = find_b(cc, a, rng) {
                let base = deg1_pos_arb(cc, a, b, rng.gen_range(0..=(G as i32 - 1)));
                d = split::add_pos(&d, &base, &f, &h, &vpl, G);
                produced += 1;
            }
        }
    }
    d
}

fn check_double_pos_arb<F: Field + std::fmt::Display>(d: &Divisor<F>, cc: &CurveConstants<F>) {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vpl = cc.vpl();
    let got = double_pos(&from_generic(d), cc);
    let expected = from_generic(&split::double_pos(d, &f, &h, &vpl, G));
    assert_eq!(got, expected, "arb double_pos mismatch for {:?}", from_generic(d));
}

fn check_add_pos_arb<F: Field + std::fmt::Display>(
    d1: &Divisor<F>,
    d2: &Divisor<F>,
    cc: &CurveConstants<F>,
) {
    let f = cc.f_poly();
    let h = cc.h_poly();
    let vpl = cc.vpl();
    let got = add_pos(&from_generic(d1), &from_generic(d2), cc);
    let expected = from_generic(&split::add_pos(d1, d2, &f, &h, &vpl, G));
    assert_eq!(
        got, expected,
        "arb add_pos mismatch for {:?} + {:?}",
        from_generic(d1),
        from_generic(d2)
    );
}

fn dbl_pos_arb<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve_arb::<PrimeField<P>, _>(&mut rng);
        let vpl = cc.vpl();
        for _ in 0..divisors {
            let d = random_valid_pos_arb(&cc, &mut rng);
            // round-trip incl. n
            assert_eq!(to_generic_arb(&from_generic(&d), &cc, &vpl), d, "arb pos round-trip");
            check_double_pos_arb(&d, &cc);
        }
        for n in 0..=2 {
            check_double_pos_arb(&deg0_pos_arb(&cc, n), &cc);
        }
        for _ in 0..8 {
            if let Some((a, b)) = random_point_arb(&cc, &mut rng) {
                check_double_pos_arb(&deg1_pos_arb(&cc, a, b, 0), &cc);
                check_double_pos_arb(&deg1_pos_arb(&cc, a, b, 1), &cc);
            }
        }
    }
}

#[test]
fn arb_dbl_pos_f7() {
    dbl_pos_arb::<7>(8001, 12, 60);
}

#[test]
fn arb_dbl_pos_f31() {
    dbl_pos_arb::<31>(8002, 12, 60);
}

#[test]
fn arb_dbl_pos_f127() {
    dbl_pos_arb::<127>(8003, 8, 40);
}

fn add_pos_arb<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve_arb::<PrimeField<P>, _>(&mut rng);
        for _ in 0..divisors {
            let d1 = random_valid_pos_arb(&cc, &mut rng);
            let d2 = random_valid_pos_arb(&cc, &mut rng);
            check_double_pos_arb(&d1, &cc);
            // adding a divisor to itself is out of the explicit ADD's domain.
            if from_generic(&d1) != from_generic(&d2) {
                check_add_pos_arb(&d1, &d2, &cc);
            }
        }
    }
}

#[test]
fn arb_add_pos_f7() {
    add_pos_arb::<7>(8201, 8, 400);
}

#[test]
fn arb_add_pos_f31() {
    add_pos_arb::<31>(8202, 8, 400);
}

#[test]
fn arb_add_pos_f127() {
    add_pos_arb::<127>(8203, 6, 200);
}

/// All arb-pos DBL (`APDBL00`..`APDBL17`) and ADD (`APADD00`..`APADD58`)
/// branches are exercised by random testing.
#[test]
fn arb_pos_branch_coverage() {
    dbl_pos_arb::<5>(8301, 60, 120);
    dbl_pos_arb::<7>(8302, 60, 120);
    dbl_pos_arb::<11>(8303, 40, 120);
    dbl_pos_arb::<31>(8304, 30, 120);
    add_pos_arb::<5>(8311, 60, 120);
    add_pos_arb::<7>(8312, 60, 120);
    add_pos_arb::<11>(8313, 40, 120);
    add_pos_arb::<31>(8314, 30, 120);

    let mut labels: Vec<String> = (0..=17).map(|i| format!("APDBL{i:02}")).collect();
    labels.extend((0..=58).map(|i| format!("APADD{i:02}")));
    let refs: Vec<&str> = labels.iter().map(|s| s.as_str()).collect();
    let missing = super::coverage::missing(&refs);
    assert!(missing.is_empty(), "arb-pos branches never exercised: {:?}", missing);
}
