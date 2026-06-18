//! Field-operation counts (M = mul, S = square, I = inverse, A = add/sub/double)
//! for the generic-branch degree-2 `add` and `double` of each split variant.
//!
//! These are field-size-independent and are the metric the explicit-formula
//! literature uses for comparison (e.g. Lange; Erickson–Jacobson–Stein;
//! Costello–Lauter). Run with:
//!   `cargo test --release g2::split::op_counts -- --nocapture`
//!
//! Counts are obtained by running the *actual* explicit formulas over a
//! `CountingField` wrapper on valid degree-2 inputs (so they reflect the code,
//! not a hand count). Setup (curve precompute + building the divisors) is not
//! counted — only the measured operation.

// Builder return types are (CurveConstants, Divisor, Divisor) tuples.
#![allow(clippy::type_complexity)]

use std::cell::Cell;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::field::{BinaryExtField, Field, PrimeField};
use crate::generic::split::{self as gsplit, Divisor};
use crate::poly::Poly;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use super::test_support::{from_generic, try_sqrt};
use super::{arbitrary as sarb, char2 as sch2, not_char2 as snch2};

// --- op counters (M, S, I, A) ---
thread_local! {
    static OPS: Cell<[u64; 4]> = const { Cell::new([0; 4]) };
}
fn bump(i: usize) {
    OPS.with(|c| {
        let mut a = c.get();
        a[i] += 1;
        c.set(a);
    });
}
fn reset() {
    OPS.with(|c| c.set([0; 4]));
}
fn snap() -> [u64; 4] {
    OPS.with(|c| c.get())
}

/// A field wrapper that tallies field operations. Squaring, doubling, and
/// inversion are counted distinctly (the formulas call `.square()`/`.double()`/
/// `.inv()`), so the tallies match how the literature reports M/S/I.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Cf<F: Field>(F);

impl<F: Field> Add for Cf<F> {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        bump(3);
        Cf(self.0 + r.0)
    }
}
impl<F: Field> Sub for Cf<F> {
    type Output = Self;
    fn sub(self, r: Self) -> Self {
        bump(3);
        Cf(self.0 - r.0)
    }
}
impl<F: Field> Mul for Cf<F> {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        bump(0);
        Cf(self.0 * r.0)
    }
}
impl<F: Field> Div for Cf<F> {
    type Output = Self;
    fn div(self, r: Self) -> Self {
        bump(2); // an inversion
        bump(0); // and a multiply
        Cf(self.0 / r.0)
    }
}
impl<F: Field> Neg for Cf<F> {
    type Output = Self;
    fn neg(self) -> Self {
        Cf(-self.0) // negation is free (XOR / sign flip)
    }
}
impl<F: Field> AddAssign for Cf<F> {
    fn add_assign(&mut self, r: Self) {
        *self = *self + r;
    }
}
impl<F: Field> SubAssign for Cf<F> {
    fn sub_assign(&mut self, r: Self) {
        *self = *self - r;
    }
}
impl<F: Field> MulAssign for Cf<F> {
    fn mul_assign(&mut self, r: Self) {
        *self = *self * r;
    }
}
impl<F: Field> DivAssign for Cf<F> {
    fn div_assign(&mut self, r: Self) {
        *self = *self / r;
    }
}
impl<F: Field> Field for Cf<F> {
    fn zero() -> Self {
        Cf(F::zero())
    }
    fn one() -> Self {
        Cf(F::one())
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn inv(&self) -> Self {
        bump(2);
        Cf(self.0.inv())
    }
    fn square(&self) -> Self {
        bump(1);
        Cf(self.0.square())
    }
    fn double(&self) -> Self {
        bump(3);
        Cf(self.0.double())
    }
    fn random<R: Rng>(rng: &mut R) -> Self {
        Cf(F::random(rng))
    }
}

impl<F: Field + std::fmt::Display> std::fmt::Display for Cf<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

// --- valid degree-2 divisor builders over an arbitrary (counted) field ---

/// Lift point `(a,b)` to a reduced degree-1 divisor and compose two of them
/// into a degree-2 divisor via the generic oracle (negative basis).
fn deg2_from_points<F: Field>(f: &Poly<F>, h: &Poly<F>, vn: &Poly<F>, p: &[(F, F)]) -> Divisor<F> {
    let deg1 = |a: F, b: F| {
        let u = Poly::from_coeffs(vec![-a, F::one()]);
        let r = (vn - &Poly::constant(b)).rem(&u);
        let vhat = vn - &r;
        let w = (f - &(&vhat * &(&vhat + h))).exact_div(&u);
        Divisor::new(u, vhat, w, 0)
    };
    gsplit::add_neg(&deg1(p[0].0, p[0].1), &deg1(p[1].0, p[1].1), f, h, vn, 2)
}

// nch2 (odd characteristic, h = 0)
fn nch2_inputs() -> (
    snch2::CurveConstants<Cf<PrimeField<65521>>>,
    Divisor<Cf<PrimeField<65521>>>,
    Divisor<Cf<PrimeField<65521>>>,
) {
    type F = Cf<PrimeField<65521>>;
    let mut rng = StdRng::seed_from_u64(2024);
    loop {
        let cc = super::test_support::random_curve::<F, _>(&mut rng);
        let (f, vn, h) = (cc.f_poly(), cc.vn(), Poly::zero());
        let mut pts = Vec::new();
        for _ in 0..2000 {
            let a = F::random(&mut rng);
            if let Some(b) = try_sqrt(f.eval(a), &mut rng) {
                pts.push((a, b));
                if pts.len() == 4 {
                    break;
                }
            }
        }
        if pts.len() < 4 {
            continue;
        }
        let d1 = deg2_from_points(&f, &h, &vn, &pts[0..2]);
        let d2 = deg2_from_points(&f, &h, &vn, &pts[2..4]);
        if d1.u.deg() == 2 && d2.u.deg() == 2 {
            return (cc, d1, d2);
        }
    }
}

// arbitrary characteristic (h != 0), over an odd prime field
fn arb_inputs() -> (
    sarb::CurveConstants<Cf<PrimeField<65521>>>,
    Divisor<Cf<PrimeField<65521>>>,
    Divisor<Cf<PrimeField<65521>>>,
) {
    type F = Cf<PrimeField<65521>>;
    let mut rng = StdRng::seed_from_u64(2025);
    let two_inv = (F::one() + F::one()).inv();
    loop {
        let h = [
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
        ];
        let y3 = F::random(&mut rng);
        if (y3.double() + h[3]).is_zero() {
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
        for _ in 0..2000 {
            let a = F::random(&mut rng);
            let (ha, fa) = (hp.eval(a), fp.eval(a));
            if let Some(s) = try_sqrt(ha * ha + (fa + fa).double(), &mut rng) {
                pts.push((a, (s - ha) * two_inv));
                if pts.len() == 4 {
                    break;
                }
            }
        }
        if pts.len() < 4 {
            continue;
        }
        let d1 = deg2_from_points(&fp, &hp, &vn, &pts[0..2]);
        let d2 = deg2_from_points(&fp, &hp, &vn, &pts[2..4]);
        if d1.u.deg() == 2 && d2.u.deg() == 2 {
            return (cc, d1, d2);
        }
    }
}

// characteristic 2, over GF(2^8)
fn char2_inputs() -> (
    sch2::CurveConstants<Cf<BinaryExtField<8>>>,
    Divisor<Cf<BinaryExtField<8>>>,
    Divisor<Cf<BinaryExtField<8>>>,
) {
    type Inner = BinaryExtField<8>;
    type F = Cf<Inner>;
    let mut rng = StdRng::seed_from_u64(2026);
    let solve = |cc: &sch2::CurveConstants<F>, a: F| -> Option<F> {
        let (ha, fa) = (cc.h_poly().eval(a), cc.f_poly().eval(a));
        (0..256u64)
            .map(|i| Cf(Inner::new(i)))
            .find(|&b| b * b + ha * b == fa)
    };
    loop {
        let beta = F::random(&mut rng);
        let f6 = beta.square() + beta;
        if f6.is_zero() {
            continue;
        }
        let cc = sch2::precompute(
            F::random(&mut rng),
            F::random(&mut rng),
            F::random(&mut rng),
            f6,
            F::random(&mut rng),
            F::random(&mut rng),
            beta,
        );
        let (fp, hp, vn) = (cc.f_poly(), cc.h_poly(), cc.vn());
        let mut pts = Vec::new();
        for ai in 0..256u64 {
            let a = Cf(Inner::new(ai));
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
        let d1 = deg2_from_points(&fp, &hp, &vn, &pts[0..2]);
        let d2 = deg2_from_points(&fp, &hp, &vn, &pts[2..4]);
        if d1.u.deg() == 2 && d2.u.deg() == 2 {
            return (cc, d1, d2);
        }
    }
}

fn fmt_row(name: &str, [m, s, i, a]: [u64; 4]) -> String {
    format!("| {name:<22} | {m:>3} | {s:>3} | {i:>2} | {a:>3} |")
}

#[test]
fn op_count_table() {
    let mut rows = Vec::new();
    let mut measure = |name: &str, counts: [u64; 4]| {
        // The affine explicit formulas perform exactly one field inversion per
        // group operation; if this ever changes, a branch is mis-selected or the
        // formula regressed.
        assert_eq!(
            counts[2], 1,
            "{name}: expected exactly 1 inversion, got {counts:?}"
        );
        rows.push(fmt_row(name, counts));
    };

    let (cc, d1, d2) = nch2_inputs();
    let (c1, c2) = (from_generic(&d1), from_generic(&d2));
    reset();
    let _ = snch2::add_neg(&c1, &c2, &cc);
    measure("nch2  deg2 add", snap());
    reset();
    let _ = snch2::double_neg(&c1, &cc);
    measure("nch2  deg2 double", snap());

    let (cc, d1, d2) = arb_inputs();
    let (c1, c2) = (from_generic(&d1), from_generic(&d2));
    reset();
    let _ = sarb::add_neg(&c1, &c2, &cc);
    measure("arb   deg2 add", snap());
    reset();
    let _ = sarb::double_neg(&c1, &cc);
    measure("arb   deg2 double", snap());

    let (cc, d1, d2) = char2_inputs();
    let (c1, c2) = (from_generic(&d1), from_generic(&d2));
    reset();
    let _ = sch2::add_neg(&c1, &c2, &cc);
    measure("char2 deg2 add", snap());
    reset();
    let _ = sch2::double_neg(&c1, &cc);
    measure("char2 deg2 double", snap());

    println!("\nGeneric-branch degree-2 field-operation counts (negative basis):\n");
    println!("| operation              |   M |   S |  I |   A |");
    println!("|------------------------|-----|-----|----|-----|");
    for r in &rows {
        println!("{r}");
    }
    println!();
}
