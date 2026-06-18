//! Shared test helpers for the genus 2 split model: conversions between the
//! explicit [`DivisorCoords`] representation and the generic Cantor
//! [`Divisor`] oracle, plus generators for random valid balanced divisors.
//!
//! These mirror the Magma testers, which validate the explicit formulas against
//! `Add_SPLIT_NEG`/`Double_SPLIT_NEG` (ported to `crate::generic::split`).

#![allow(dead_code)]

use crate::field::Field;
use crate::generic::split::{self, Divisor};
use crate::poly::Poly;
use rand::Rng;

use super::not_char2::{precompute, CurveConstants, DivisorCoords, G};

// ---------------------------------------------------------------------------
// Conversions: explicit DivisorCoords  <->  generic Divisor
// ---------------------------------------------------------------------------

/// Convert explicit coords to a generic divisor, reconstructing the full
/// degree-(g+1) `v` from the `basis` polynomial's high coefficients.
///
/// `basis` is `Vn` for the negative basis, `Vpl` for the positive basis.
pub fn to_generic<F: Field>(
    d: &DivisorCoords<F>,
    cc: &CurveConstants<F>,
    basis: &Poly<F>,
) -> Divisor<F> {
    let f = cc.f_poly();
    let u = match d.degree() {
        0 => Poly::constant(F::one()),
        1 => Poly::from_coeffs(vec![d.u0, F::one()]),
        2 => Poly::from_coeffs(vec![d.u0, d.u1, F::one()]),
        _ => unreachable!(),
    };
    // v = basis_c3·x³ + basis_c2·x² + v1·x + v0
    let v = Poly::from_coeffs(vec![d.v0, d.v1, basis.coeff(2), basis.coeff(3)]);
    let v_sq = &v * &v;
    let w = (&f - &v_sq).exact_div(&u);
    Divisor::new(u, v, w, d.n)
}

/// Read explicit coords out of a generic divisor (basis-independent: the high
/// `v` coefficients are discarded, the low two and `n` are kept).
pub fn from_generic<F: Field>(d: &Divisor<F>) -> DivisorCoords<F> {
    match d.u.deg() {
        -1 | 0 => DivisorCoords::deg0(d.v.coeff(1), d.v.coeff(0), d.n),
        1 => DivisorCoords::deg1(d.u.coeff(0), d.v.coeff(1), d.v.coeff(0), d.n),
        2 => DivisorCoords::deg2(d.u.coeff(1), d.u.coeff(0), d.v.coeff(1), d.v.coeff(0), d.n),
        other => panic!("split g2: unexpected deg(u) = {other}"),
    }
}

// ---------------------------------------------------------------------------
// Random curve + divisor generation (mirrors RandomDivisor_SPLIT_NEG)
// ---------------------------------------------------------------------------

/// The formal derivative of a polynomial.
fn derivative<F: Field>(p: &Poly<F>) -> Poly<F> {
    if p.deg() < 1 {
        return Poly::zero();
    }
    let mut c = Vec::with_capacity(p.coeffs.len() - 1);
    let mut i_fe = F::zero();
    for i in 1..p.coeffs.len() {
        i_fe += F::one(); // = i as a field element
        c.push(p.coeffs[i] * i_fe);
    }
    Poly::from_coeffs(c)
}

/// Whether `f` defines a smooth genus-2 split curve: squarefree of degree 6
/// (for nch2 the leading coeff `f6 = 1` is a square, giving two points at ∞).
fn is_valid_curve<F: Field>(cc: &CurveConstants<F>) -> bool {
    let f = cc.f_poly();
    if f.deg() != 6 {
        return false;
    }
    let fp = derivative(&f);
    if fp.is_zero() {
        return false;
    }
    // squarefree ⇔ gcd(f, f') is a (nonzero) constant
    f.xgcd(&fp).0.deg() == 0
}

/// Random not-char-2 split curve constants (`f = x⁶ + f4x⁴ + f3x³ + f2x² + f1x + f0`),
/// rejecting singular (non-squarefree) curves like the Magma generator does.
pub fn random_curve<F: Field, R: Rng>(rng: &mut R) -> CurveConstants<F> {
    loop {
        let cc = precompute(
            F::random(rng),
            F::random(rng),
            F::random(rng),
            F::random(rng),
            F::random(rng),
        );
        if is_valid_curve(&cc) {
            return cc;
        }
    }
}

/// Try to find a square root of `target` by trial (works for the small test
/// fields). Returns `None` if `target` is a non-residue / not found.
pub fn try_sqrt<F: Field, R: Rng>(target: F, rng: &mut R) -> Option<F> {
    if target.is_zero() {
        return Some(F::zero());
    }
    for _ in 0..30000 {
        let c = F::random(rng);
        if c.square() == target {
            return Some(c);
        }
    }
    None
}

/// The balanced neutral divisor in the negative basis: `<1, Vn, f − Vn², 1>`.
pub fn neutral_neg<F: Field>(cc: &CurveConstants<F>) -> Divisor<F> {
    let f = cc.f_poly();
    let vn = cc.vn();
    let w = &f - &(&vn * &vn);
    // Neutral_SPLIT weight = Ceiling((deg(f) − 2)/2 / 2) = 1 for g = 2.
    Divisor::new(Poly::constant(F::one()), vn, w, 1)
}

/// Generate a random valid reduced divisor in the negative basis by composing
/// `g` random base divisors (degree-0 or degree-1) through the trusted oracle.
pub fn random_valid_neg<F: Field, R: Rng>(cc: &CurveConstants<F>, rng: &mut R) -> Divisor<F> {
    let f = cc.f_poly();
    let vn = cc.vn();
    let h = Poly::zero();
    let mut d = neutral_neg(cc);

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
            // Degree-0 base divisor with random weight in [0, g].
            let base = Divisor::new(
                Poly::constant(F::one()),
                vn.clone(),
                &f - &(&vn * &vn),
                rng.gen_range(0..=G as i32),
            );
            d = split::add_neg(&d, &base, &f, &h, &vn, G);
            produced += 1;
        } else {
            // Monic u = x + u0 (root a = −u0); need a point (a, b): b² = f(a).
            let u0 = r0 * r1.inv();
            let a = -u0;
            let fa = f.eval(a);
            if let Some(b) = try_sqrt(fa, rng) {
                let u = Poly::from_coeffs(vec![u0, F::one()]);
                // vhat = Vn − ((Vn − b) mod u)  (reduced basis)
                let r = (&vn - &Poly::constant(b)).rem(&u);
                let vhat = &vn - &r;
                let w = (&f - &(&vhat * &vhat)).exact_div(&u);
                let base = Divisor::new(u, vhat, w, rng.gen_range(0..=(G as i32 - 1)));
                d = split::add_neg(&d, &base, &f, &h, &vn, G);
                produced += 1;
            }
        }
    }
    d
}

/// A degree-0 divisor `<1, Vn, f − Vn², n>` in the negative basis.
pub fn deg0_neg<F: Field>(cc: &CurveConstants<F>, n: i32) -> Divisor<F> {
    let f = cc.f_poly();
    let vn = cc.vn();
    let w = &f - &(&vn * &vn);
    Divisor::new(Poly::constant(F::one()), vn, w, n)
}

/// Find a random curve point `(a, b)` with `b² = f(a)` (negative basis: no `h`).
pub fn random_point<F: Field, R: Rng>(cc: &CurveConstants<F>, rng: &mut R) -> Option<(F, F)> {
    let f = cc.f_poly();
    for _ in 0..200 {
        let a = F::random(rng);
        if let Some(b) = try_sqrt(f.eval(a), rng) {
            return Some((a, b));
        }
    }
    None
}

/// A degree-1 reduced divisor for the curve point `(a, b)` with weight `n`
/// (negative basis): `u = x − a`, `v = Vn − ((Vn − b) mod u)`.
pub fn deg1_neg<F: Field>(cc: &CurveConstants<F>, a: F, b: F, n: i32) -> Divisor<F> {
    let f = cc.f_poly();
    let vn = cc.vn();
    let u = Poly::from_coeffs(vec![-a, F::one()]);
    let r = (&vn - &Poly::constant(b)).rem(&u);
    let vhat = &vn - &r;
    let w = (&f - &(&vhat * &vhat)).exact_div(&u);
    Divisor::new(u, vhat, w, n)
}

// --- positive basis variants (use Vpl; oracle add_pos/double_pos) ---

/// The balanced neutral divisor in the positive basis: `<1, Vpl, f − Vpl², 1>`.
pub fn neutral_pos<F: Field>(cc: &CurveConstants<F>) -> Divisor<F> {
    let f = cc.f_poly();
    let vpl = cc.vpl();
    let w = &f - &(&vpl * &vpl);
    Divisor::new(Poly::constant(F::one()), vpl, w, 1)
}

/// A degree-0 divisor `<1, Vpl, f − Vpl², n>` in the positive basis.
pub fn deg0_pos<F: Field>(cc: &CurveConstants<F>, n: i32) -> Divisor<F> {
    let f = cc.f_poly();
    let vpl = cc.vpl();
    let w = &f - &(&vpl * &vpl);
    Divisor::new(Poly::constant(F::one()), vpl, w, n)
}

/// A degree-1 reduced divisor for point `(a, b)` with weight `n` (positive basis).
pub fn deg1_pos<F: Field>(cc: &CurveConstants<F>, a: F, b: F, n: i32) -> Divisor<F> {
    let f = cc.f_poly();
    let vpl = cc.vpl();
    let u = Poly::from_coeffs(vec![-a, F::one()]);
    let r = (&vpl - &Poly::constant(b)).rem(&u);
    let vhat = &vpl - &r;
    let w = (&f - &(&vhat * &vhat)).exact_div(&u);
    Divisor::new(u, vhat, w, n)
}

/// Random valid reduced divisor in the positive basis (composes via `add_pos`).
pub fn random_valid_pos<F: Field, R: Rng>(cc: &CurveConstants<F>, rng: &mut R) -> Divisor<F> {
    let f = cc.f_poly();
    let vpl = cc.vpl();
    let h = Poly::zero();
    let mut d = neutral_pos(cc);

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
            let base = Divisor::new(
                Poly::constant(F::one()),
                vpl.clone(),
                &f - &(&vpl * &vpl),
                rng.gen_range(0..=G as i32),
            );
            d = split::add_pos(&d, &base, &f, &h, &vpl, G);
            produced += 1;
        } else {
            let u0 = r0 * r1.inv();
            let a = -u0;
            if let Some(b) = try_sqrt(f.eval(a), rng) {
                let u = Poly::from_coeffs(vec![u0, F::one()]);
                let r = (&vpl - &Poly::constant(b)).rem(&u);
                let vhat = &vpl - &r;
                let w = (&f - &(&vhat * &vhat)).exact_div(&u);
                let base = Divisor::new(u, vhat, w, rng.gen_range(0..=(G as i32 - 1)));
                d = split::add_pos(&d, &base, &f, &h, &vpl, G);
                produced += 1;
            }
        }
    }
    d
}

/// Assert a generic divisor is a well-formed reduced balanced divisor:
/// `u` monic, `w = (f − v²)/u` exact, and `0 ≤ n ≤ g − deg(u)`.
pub fn assert_well_formed<F: Field + std::fmt::Display>(d: &Divisor<F>, cc: &CurveConstants<F>) {
    let f = cc.f_poly();
    let deg_u = d.u.deg().max(0);
    assert!(deg_u <= G as i32, "deg(u) = {} > g", deg_u);
    if d.u.deg() >= 1 {
        assert!(d.u.is_monic(), "u not monic: {:?}", d.u);
    }
    let v_sq = &d.v * &d.v;
    let (q, r) = (&f - &v_sq).div_rem(&d.u);
    assert!(r.is_zero(), "w = (f − v²)/u not exact");
    assert_eq!(q, d.w, "stored w disagrees with (f − v²)/u");
    assert!(
        d.n >= 0 && d.n <= G as i32 - deg_u,
        "balance weight n = {} out of range [0, {}]",
        d.n,
        G as i32 - deg_u
    );
}
