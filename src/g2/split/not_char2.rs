//! Not-characteristic-2 genus 2 **split** model divisor arithmetic.
//!
//! Explicit formulas for divisor addition and doubling on split (two points at
//! infinity) genus 2 hyperelliptic curves over fields NOT of characteristic 2.
//!
//! ## Curve form
//! `y² = f(x)` where (after the standard split-model normalization)
//! - `f(x) = x⁶ + f4·x⁴ + f3·x³ + f2·x² + f1·x + f0` (monic, `f6 = 1`, `f5 = 0`)
//! - `h(x) = 0`
//!
//! ## Divisor representation
//! Divisors use the 4-coordinate balanced Mumford representation `(u, v, w, n)`
//! where `u` is monic of degree ≤ g = 2, `w = (f − v²)/u`, and `n` is the integer
//! **balance weight** (`0 ≤ n ≤ g − deg(u)`). The explicit formulas carry only the
//! low coordinates `(u2, u1, u0, v1, v0)` plus `n`; the high coefficients of `v`
//! (degree g+1 = 3 and g = 2) are fixed by the reduced-basis polynomial `Vpl`
//! (positive basis) or `Vn = −Vpl − h` (negative basis).
//!
//! Two reduced bases exist; functions are suffixed `_neg` / `_pos`.
//!
//! Based on: Sebastian Lindner, 2020 (Magma `*_splitG2_{UTL,ADD,DBL}.mag`).

use crate::field::Field;
use crate::poly::Poly;

/// Genus of this module's curves.
pub const G: usize = 2;

/// Curve constants and precomputations for a not-char-2 split genus 2 curve.
///
/// Built by [`precompute`] from the five free coefficients of
/// `f = x⁶ + f4·x⁴ + f3·x³ + f2·x² + f1·x + f0`.
#[derive(Clone, Copy, Debug)]
pub struct CurveConstants<F: Field> {
    // --- curve coefficients (f6 = 1, f5 = 0, h = 0) ---
    pub f0: F,
    pub f1: F,
    pub f2: F,
    pub f3: F,
    pub f4: F,

    // --- positive reduced basis Vpl = x³ + y1·x + y0 (y3 = 1, y2 = 0) ---
    pub y0: F,
    pub y1: F,
    // --- negative reduced basis Vn = −Vpl = −x³ + yn1·x + yn0 (yn3 = −1, yn2 = 0) ---
    pub yn0: F,
    pub yn1: F,

    // --- precomputations used by the explicit ADD/DBL formulas ---
    /// `1 / 2` (Magma `d5`); the not-char-2 formulas divide by 2.
    pub half: F,
    /// `d1 = f2 − y1²`  (leading coeff of `f − V²` reduced to degree ≤ 2)
    pub d1: F,
    /// `f1 / 2`  (Magma `d7` in neg / `d3` in pos)
    pub half_f1: F,
    /// `d1 / 2`  (Magma `d8` in neg / `d2` in pos)
    pub half_d1: F,

    // --- precomputed degree-0 "adjust" divisor `<adu, adv, ·>` ---
    // `adu = (f − V²)` made monic; shared by both bases.
    pub au2: F,
    pub au1: F,
    pub au0: F,
    pub audeg: i8,
    // `adv = adv1·x + adv0`, differs per basis.
    pub adv1_neg: F,
    pub adv0_neg: F,
    pub adv1_pos: F,
    pub adv0_pos: F,
}

impl<F: Field> CurveConstants<F> {
    /// The defining polynomial `f` as a [`Poly`] (degree 6, monic, `f5 = 0`).
    pub fn f_poly(&self) -> Poly<F> {
        Poly::from_coeffs(vec![
            self.f0,
            self.f1,
            self.f2,
            self.f3,
            self.f4,
            F::zero(),
            F::one(),
        ])
    }

    /// The positive reduced-basis polynomial `Vpl = x³ + y1·x + y0`.
    pub fn vpl(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.y0, self.y1, F::zero(), F::one()])
    }

    /// The negative reduced-basis polynomial `Vn = −Vpl = −x³ + yn1·x + yn0`.
    pub fn vn(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.yn0, self.yn1, F::zero(), -F::one()])
    }
}

/// Build [`CurveConstants`] from the free coefficients of
/// `f = x⁶ + f4·x⁴ + f3·x³ + f2·x² + f1·x + f0`.
///
/// Port of `Precompute` from `nch2_splitG2_UTL.mag` (both bases). Panics on a
/// characteristic-2 field (the not-char-2 formulas divide by 2).
pub fn precompute<F: Field>(f0: F, f1: F, f2: F, f3: F, f4: F) -> CurveConstants<F> {
    let two = F::one() + F::one();
    assert!(!two.is_zero(), "not_char2 split model requires characteristic ≠ 2");
    let half = two.inv();

    // Positive basis Vpl = x³ + y1 x + y0 with y1 = f4/2, y0 = f3/2.
    let y1 = f4 * half;
    let y0 = f3 * half;
    let yn1 = -y1;
    let yn0 = -y0;

    let d1 = f2 - y1.square();
    let half_f1 = f1 * half;
    let half_d1 = d1 * half;

    // Degree-0 adjust divisor: adu = (f − V²) made monic.
    // k2 = d1, k1 = f1 − 2·y0·y1, k0 = f0 − y0² (identical for both bases).
    let k2 = d1;
    let k1 = f1 - (y0 * y1).double();
    let k0 = f0 - y0.square();

    let (au2, au1, au0, audeg, adv1_neg, adv0_neg, adv1_pos, adv0_pos);

    if k2.is_zero() {
        if k1.is_zero() {
            // audeg 0
            au2 = F::zero();
            au1 = F::zero();
            au0 = F::one();
            audeg = 0;
            adv1_neg = yn1;
            adv0_neg = yn0;
            adv1_pos = y1;
            adv0_pos = y0;
        } else {
            // audeg 1
            let w1 = k1.inv();
            au2 = F::zero();
            au1 = F::one();
            au0 = k0 * w1;
            audeg = 1;
            let au0sq = au0.square();
            adv1_neg = yn1;
            adv0_neg = y0 - au0 * (f4 + au0sq.double());
            adv1_pos = y1;
            adv0_pos = -y0 + au0 * (y1.double() + au0sq.double());
        }
    } else {
        // audeg 2
        let w1 = k2.inv();
        au2 = F::one();
        au1 = k1 * w1;
        au0 = k0 * w1;
        audeg = 2;
        let w2 = -au1.double();
        adv1_neg = y1 - au0.double() - au1 * w2;
        adv0_neg = y0 - au0 * w2;
        adv1_pos = au0.double() + au1 * w2 - y1;
        adv0_pos = au0 * w2 - y0;
    }

    CurveConstants {
        f0,
        f1,
        f2,
        f3,
        f4,
        y0,
        y1,
        yn0,
        yn1,
        half,
        d1,
        half_f1,
        half_d1,
        au2,
        au1,
        au0,
        audeg,
        adv1_neg,
        adv0_neg,
        adv1_pos,
        adv0_pos,
    }
}

/// A divisor on a split genus 2 curve, in the 4-coordinate balanced Mumford
/// representation, carrying the integer balance weight `n`.
///
/// - `u2 == 1`: degree 2, `u = x² + u1·x + u0`
/// - `u2 == 0, u1 == 1`: degree 1, `u = x + u0`
/// - `u2 == 0, u1 == 0`: degree 0, `u = 1`
///
/// `v1, v0` are the low coefficients of the (degree 3) `v`; the high
/// coefficients come from the basis polynomial (`Vn` for neg, `Vpl` for pos).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DivisorCoords<F: Field> {
    pub u2: F,
    pub u1: F,
    pub u0: F,
    pub v1: F,
    pub v0: F,
    pub n: i32,
}

impl<F: Field> DivisorCoords<F> {
    /// Degree of `u` (0, 1, or 2).
    #[inline]
    pub fn degree(&self) -> usize {
        if self.u2.is_one() {
            2
        } else if self.u1.is_one() {
            1
        } else {
            0
        }
    }

    /// A degree-2 divisor `<x² + u1·x + u0, v1·x + v0, n>`.
    pub fn deg2(u1: F, u0: F, v1: F, v0: F, n: i32) -> Self {
        Self { u2: F::one(), u1, u0, v1, v0, n }
    }

    /// A degree-1 divisor `<x + u0, v1·x + v0, n>` (v1 holds the basis x-coeff).
    pub fn deg1(u0: F, v1: F, v0: F, n: i32) -> Self {
        Self { u2: F::zero(), u1: F::one(), u0, v1, v0, n }
    }

    /// A degree-0 divisor `<1, v1·x + v0, n>` (v1,v0 hold the basis coeffs).
    pub fn deg0(v1: F, v0: F, n: i32) -> Self {
        Self { u2: F::zero(), u1: F::zero(), u0: F::one(), v1, v0, n }
    }
}

// ===========================================================================
// Doubling — negative reduced basis
// Port of `nch2_splitG2_DBL.mag` (negReduced). Branch labels `DBLnn` from Magma
// are preserved in comments.
// ===========================================================================

/// Double a degree-1 divisor `<x + u0, v0, n=1>` (DWN case). (`Deg1DBLDWN`)
#[inline]
fn deg1_dbl_dwn_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (yn0, yn1, f1, f4, d8) = (cc.yn0, cc.yn1, cc.f1, cc.f4, cc.half_d1);

    // v := -V - h - ((-V-h - v) mod u)
    let u0sq = u0.square();
    let t0 = -u0 * (f4 + u0sq.double());
    let vp0 = v0 - t0;

    // d := (2v + h) mod u
    let d = v0 + vp0;
    if d.is_zero() {
        super::branch("DBL00");
        return DivisorCoords::deg0(yn1, yn0, 2);
    }

    let z0 = -vp0 - yn0;
    let t1 = z0 * u0;
    let z1 = d8 - t1;
    let t2 = u0 * (z1.double() - t1);
    let sp0 = f1 - f4 * vp0 - t2.double();

    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("DBL01");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        super::branch("DBL02");
        return DivisorCoords::deg1(upp0, yn1, -vp0, 0);
    }

    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = -z0 * w2 - u0;
    let upp0 = -yn1 + s0 * cc.half - z1 * w2 - u0 * upp1;
    let vpp1 = yn1 - s0;
    let vpp0 = -vp0 - s0 * u0;
    super::branch("DBL03");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-1 divisor `<x + u0, v0, n=0>` (UP case). (`Deg1DBLUP`)
#[inline]
fn deg1_dbl_up_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (yn0, yn1, f1, f4, d8) = (cc.yn0, cc.yn1, cc.f1, cc.f4, cc.half_d1);

    // d := (2v + h) mod u
    let u0sq = u0.square();
    let d = v0.double() + u0 * (f4 + u0sq.double());
    if d.is_zero() {
        super::branch("DBL04");
        return DivisorCoords::deg0(yn1, yn0, 0);
    }

    let z0 = v0 - yn0;
    let t1 = u0 * z0;
    let z1 = d8 - t1;
    let t2 = u0 * (z1.double() - t1);
    let sp0 = f1 + f4 * v0 - t2.double();

    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("DBL05");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        let upp2 = upp0.square();
        let vpp0 = -v0 - upp0 * (f4 + upp2.double());
        super::branch("DBL06");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }

    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = z0 * w2 - u0;
    let upp0 = z1 * w2 - yn1 - s0 * cc.half - u0 * upp1;
    let t0 = upp1.double();
    let vpp1 = upp1 * t0 - upp0.double() - yn1 - s0;
    let vpp0 = upp0 * t0 - v0 - s0 * u0;
    super::branch("DBL07");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-2 divisor `<x² + u1x + u0, v1x + v0, n=0>`. (`Deg2DBL`)
#[inline]
fn deg2_dbl_neg<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f2, yn0, yn1) = (cc.f2, cc.yn0, cc.yn1);

    // d := Resultant(u, 2v + h) via 2x2 system
    let t0 = u1.square();
    let u12 = u1.double();
    let t1 = t0 - u0;
    let t2 = t1 - v1;
    let m3 = t2.double();
    let m4 = v0.double() - u0 * u12;
    let m1 = m4 + m3 * u1;
    let m2 = -m3 * u0;
    let d = m4 * m1 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            super::branch("DBL08");
            return DivisorCoords::deg0(yn1, yn0, 1);
        }
        let b1 = -m3.inv();
        let t2b = v1 - yn1;
        let t3 = v0 - yn0;
        let k2 = t2b.double();
        let k1 = t3.double() - u1 * k2;
        let k0 = f2 - v1.square() - u0 * k2 - u1 * k1;
        let u0n = u1 - m4 * b1;
        let s0 = b1 * (k0 - u0n * (k1 - u0n * k2));
        let upp1 = u0n.double();
        let upp0 = u0n.square();
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0n * s0;
        super::branch("DBL09");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    // Compute s(d)/c3
    let z0 = v0 - yn0;
    let z1 = v1 - yn1;
    let r1 = z0 - z1 * u12;
    let r0 = (f2 - v1.square()) * cc.half - u1 * r1 - z1 * (t0 + u0.double());
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("DBL10");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let t1b = z1 * cc.half;
        let upp0 = -s0 - u1 + w4 * t1b;
        let t0b = upp0 * (s0 * u1 + t1b + yn1 - upp0 * (s0 + upp0)) - s0 * u0;
        let vpp0 = t0b.double() - v0;
        super::branch("DBL11");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }

    let w1 = sp1 * (sp1 - d); // sp1² + sp1·d
    if w1.is_zero() {
        let w0 = sp0 + sp1 * u1;
        if w0.is_zero() {
            super::branch("DBL12");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let w4 = s0 + u1;
        let t3 = z1 * cc.half;
        let upp0 = w3 * (w4 * s0 - t2 - t3);
        let t1b = upp0 * (t3 + yn1 + s0 * u1 + u0 - upp0 * w4) - s0 * u0;
        let vpp0 = t1b.double() - v0;
        super::branch("DBL13");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }

    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * d * dd;
    let s0 = w3 * sp0;
    let s1 = w3 * sp1;
    let t4 = s0 + u1;

    let upp1 = w4 * ((s0 + t4) * s1 - s0);
    let upp0 = w4 * (t4 * s0 - t2 * s1 - z1 * cc.half);

    let zz0 = upp0 - u0;
    let zz1 = upp1 - u1;
    let ww0 = zz0 * s0;
    let ww1 = zz1 * s1;
    let ww2 = upp1 - ww1;
    let t1b = (s0 + s1) * (zz0 + zz1) - ww0 - ww1 - upp0 + ww2 * upp1;
    let t0b = ww0 + ww2 * upp0;
    let vpp1 = t1b.double() - v1;
    let vpp0 = t0b.double() - v0;
    super::branch("DBL14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a divisor in the **negative** reduced basis. (`DBL` dispatcher)
#[inline]
pub fn double_neg<F: Field>(d: &DivisorCoords<F>, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    match d.degree() {
        2 => deg2_dbl_neg(d.u1, d.u0, d.v1, d.v0, cc),
        1 => {
            if d.n == 0 {
                deg1_dbl_up_neg(d.u0, d.v0, cc)
            } else {
                deg1_dbl_dwn_neg(d.u0, d.v0, cc)
            }
        }
        _ => {
            // degree 0: <1, V, n>; the adjust divisor `adu` has degree `audeg`.
            if d.n == 1 {
                super::branch("DBL15");
                DivisorCoords::deg0(cc.yn1, cc.yn0, 1)
            } else if d.n == 0 {
                super::branch("DBL16");
                // <adu, adv, 2 − audeg>
                DivisorCoords {
                    u2: cc.au2,
                    u1: cc.au1,
                    u0: cc.au0,
                    v1: cc.adv1_neg,
                    v0: cc.adv0_neg,
                    n: 2 - cc.audeg as i32,
                }
            } else {
                super::branch("DBL17");
                // <adu, V, 0>
                DivisorCoords {
                    u2: cc.au2,
                    u1: cc.au1,
                    u0: cc.au0,
                    v1: cc.yn1,
                    v0: cc.yn0,
                    n: 0,
                }
            }
        }
    }
}

// ===========================================================================
// Addition — negative reduced basis
// Port of `nch2_splitG2_ADD.mag` (negReduced). Branch labels `ADDnn` preserved.
// In every sub-function `(u0[,u1],v0[,v1])` is the first operand and
// `(up0[,up1],vp0[,vp1])` the second; `v1`/`vp1` are the low x-coefficients.
// ===========================================================================

/// `<1, V, 0>` + degree-1 divisor `<x + u0, v0, 1>`. (`Deg01ADDUP`)
#[inline]
fn deg01_add_up_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f4, yn0, yn1, d7, d8) = (cc.f4, cc.yn0, cc.yn1, cc.half_f1, cc.half_d1);
    let z3 = v0 - yn0;
    if z3.is_zero() {
        if d8.is_zero() {
            super::branch("ADD00");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = d8.inv();
        let up0 = w2 * (d7 - v0 * yn1) - u0;
        let t1 = up0.square();
        let vp0 = -up0 * (f4 + t1.double()) - v0;
        super::branch("ADD01");
        return DivisorCoords::deg1(up0, yn1, vp0, 1);
    }
    let w2 = z3.inv();
    let up1 = w2 * d8 - u0;
    let up0 = w2 * (d7 - v0 * yn1) - u0 * up1;
    let t0 = up1.square() - up0;
    let vp1 = t0.double() - yn1;
    let vp0 = up0 * up1.double() - v0;
    super::branch("ADD02");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 2>` + degree-1 divisor `<x + u0, v0, 0>`. (`Deg01ADDDWN`)
#[inline]
fn deg01_add_dwn_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f4, yn0, yn1, d7, d8) = (cc.f4, cc.yn0, cc.yn1, cc.half_f1, cc.half_d1);
    let t0 = u0.square();
    let vp0 = -v0 - u0 * (f4 + t0.double());
    let z3 = vp0 - yn0;
    if z3.is_zero() {
        if d8.is_zero() {
            super::branch("ADD03");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = d8.inv();
        let up0 = w2 * (d7 - vp0 * yn1) - u0;
        super::branch("ADD04");
        return DivisorCoords::deg1(up0, yn1, vp0, 0);
    }
    let w2 = z3.inv();
    let up1 = w2 * d8 - u0;
    let up0 = w2 * (d7 - vp0 * yn1) - u0 * up1;
    super::branch("ADD05");
    DivisorCoords::deg2(up1, up0, yn1, vp0, 0)
}

/// `<1, V, 0>` + degree-2 divisor. (`Deg02ADDUP`)
#[inline]
fn deg02_add_up_neg<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (yn0, yn1, f4, f2) = (cc.yn0, cc.yn1, cc.f4, cc.f2);
    let z1 = v1 - yn1;
    let z0 = v0 - yn0;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("ADD06");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z0.inv();
        let up0 = w2 * (f2 - v1.square()) * cc.half - u1;
        let t0 = up0.square();
        let vp0 = -up0 * (f4 + t0.double()) - v0;
        super::branch("ADD07");
        return DivisorCoords::deg1(up0, yn1, vp0, 1);
    }
    let w2 = z1.inv();
    let up1 = w2 * z0 - u1;
    let up0 = w2 * (f2 - v1.square()) * cc.half - u0 - u1 * up1;
    let t0 = up1.square() - up0;
    let t1 = up1.double();
    let vp1 = t0.double() - v1;
    let vp0 = t1 * up0 - v0;
    super::branch("ADD08");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 2>` + degree-2 divisor. (`Deg02ADDDWN`)
#[inline]
fn deg02_add_dwn_neg<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f2, yn0, yn1) = (cc.f2, cc.yn0, cc.yn1);
    let t1 = u1.double();
    let t2 = u0 - u1.square();
    let v1 = v1 + t2.double();
    let v0 = v0 - u0 * t1;
    let z1 = -v1 - yn1;
    let z0 = -v0 - yn0;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("ADD09");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z0.inv();
        let up0 = w2 * (f2 - v1.square()) * cc.half - u1;
        super::branch("ADD10");
        return DivisorCoords::deg1(up0, yn1, -v0, 0);
    }
    let w2 = z1.inv();
    let t2 = w2 * z0;
    let up1 = t2 - u1;
    let up0 = w2 * (f2 - v1.square()) * cc.half - u0 - up1 * u1;
    super::branch("ADD11");
    DivisorCoords::deg2(up1, up0, -v1, -v0, 0)
}

/// Two degree-1 divisors with `n + np = 1`. (`Deg1ADD`)
#[inline]
fn deg1_add_neg<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (yn0, yn1, f4, f1, d1) = (cc.yn0, cc.yn1, cc.f4, cc.f1, cc.d1);
    let d = u0 - up0;
    if d.is_zero() {
        let upp0 = u0.square();
        let upp1 = u0.double();
        let dw = v0 + vp0 + u0 * f4 + upp0 * upp1;
        if dw.is_zero() {
            super::branch("ADD12");
            return DivisorCoords::deg0(yn1, yn0, 1);
        }
        let t0 = v0 - yn0;
        let t2 = f1 + f4 * v0 - u0 * (d1.double() - upp1 * (t0 + t0 + t0));
        let s0 = t2 * dw.inv();
        let vpp1 = yn1 + s0;
        let vpp0 = v0 + s0 * u0;
        super::branch("ADD13");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let s0 = (vp0 - v0) * d.inv();
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;
    let vpp1 = yn1 + s0;
    let vpp0 = v0 + u0 * s0;
    super::branch("ADD14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 0` (UP). (`Deg1ADDUP`)
#[inline]
fn deg1_add_up_neg<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (yn0, yn1, d8, f4) = (cc.yn0, cc.yn1, cc.half_d1, cc.f4);
    let d = u0 - up0;
    if d.is_zero() {
        super::branch("ADD15");
        return DivisorCoords::deg0(yn1, yn0, 0);
    }
    let sp0 = vp0 - v0;
    let z0 = v0 - yn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("ADD16");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w1 = z0.inv();
        let upp0 = w1 * (d8 - z0 * u0) - up0;
        let t1 = upp0.square();
        let vpp0 = -upp0 * (f4 + t1.double()) - v0;
        super::branch("ADD17");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t1 = z0 * w1;
    let upp1 = t1 - up0;
    let upp0 = -yn1 + w1 * d8 - t1 * u0 - s0 * cc.half - up0 * upp1;
    let t0 = upp1.square() - upp0;
    let vpp1 = -yn1 - s0 + t0.double();
    let vpp0 = -v0 - s0 * u0 + upp0 * upp1.double();
    super::branch("ADD18");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 1` (DWN). (`Deg1ADDDWN`)
#[inline]
fn deg1_add_dwn_neg<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (yn0, yn1, d8, f4) = (cc.yn0, cc.yn1, cc.half_d1, cc.f4);
    let d = u0 - up0;
    if d.is_zero() {
        super::branch("ADD19");
        return DivisorCoords::deg0(yn1, yn0, 2);
    }
    let t0 = u0.square();
    let vh0 = -v0 - u0 * (f4 + t0.double());
    let t1 = up0.square();
    let sp0 = vp0 + vh0 + up0 * (f4 + t1.double());
    let z0 = vh0 - yn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("ADD20");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w1 = z0.inv();
        let upp0 = w1 * d8 - u0 - up0;
        super::branch("ADD21");
        return DivisorCoords::deg1(upp0, yn1, vh0, 0);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = z0 * w1;
    let upp1 = -t0 - up0;
    let upp0 = -yn1 + s0 * cc.half - w1 * d8 + t0 * u0 - up0 * upp1;
    let vpp1 = yn1 - s0;
    let vpp0 = vh0 - s0 * u0;
    super::branch("ADD22");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,1>` + degree-2 (UP). (`Deg12ADDUP`)
#[inline]
fn deg12_add_up_neg<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, yn0, yn1) = (cc.f2, cc.yn0, cc.yn1);
    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 - t0 + t1;
    if d.is_zero() {
        let dw = v0 + vp0 - u0 * (vp1 + yn1 - t1.double());
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let vpp0 = vp0 + upp0 * (yn1 - vp1);
            super::branch("ADD23");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let t2 = vp1 - yn1;
        let t3 = vp0 - yn0;
        let k2 = t2.double();
        let k1 = t3.double() - up1 * k2;
        let k0 = f2 - vp1.square() - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("ADD24");
                return DivisorCoords::deg0(yn1, yn0, 2);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 - up1 - u0;
            let t0 = upp0.square();
            let vpp0 = upp0 * (vp1 + yn1 - t0.double()) - vp0;
            super::branch("ADD25");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1;
        let upp1 = w4 * t2 - s0 * cc.half - u0;
        let upp0 = w4 * (vp0 - yn0 - t2 * up1) - (vh1 + vp1) * cc.half - u0 * upp1;
        let t1 = s0 + upp1.double();
        let vpp1 = t1 * upp1 - upp0.double() - vh1;
        let vpp0 = upp0 * t1 - s0 * up0 - vp0;
        super::branch("ADD26");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = vp1 - yn1;
    let sp0 = v0 - vp0 + z0 * u0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("ADD27");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (vp0 - yn0) - up1 - u0;
        let t0 = upp0.square();
        let vpp0 = upp0 * (vp1 + yn1 - t0.double()) - vp0;
        super::branch("ADD28");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1;
    let upp1 = w4 * z0 - s0 * cc.half - u0;
    let upp0 = w4 * (vp0 - yn0 - z0 * up1) - (vh1 + vp1) * cc.half - u0 * upp1;
    let t1 = s0 + upp1.double();
    let vpp1 = t1 * upp1 - upp0.double() - vh1;
    let vpp0 = upp0 * t1 - s0 * up0 - vp0;
    super::branch("ADD29");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,0>` + degree-2 (DWN). (`Deg12ADD`)
#[inline]
fn deg12_add_neg<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (yn0, yn1, f2) = (cc.yn0, cc.yn1, cc.f2);
    // vp := -V - h - ((-V-h - vp) mod up)
    let t2 = up1.double();
    let vp1 = vp1 + up0.double() - t2 * up1;
    let vp0 = vp0 - up0 * t2;
    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 - t0 + t1;
    if d.is_zero() {
        let dw = vp0 + v0 - u0 * (yn1 + vp1);
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let t0 = upp0.square();
            let vpp0 = vp0 + upp0 * (yn1 - vp1 - t0.double());
            super::branch("ADD30");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let t2 = -vp1 - yn1;
        let t3 = -vp0 - yn0;
        let k2 = t2.double();
        let k1 = t3.double() - up1 * k2;
        let k0 = f2 - vp1.square() - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("ADD31");
                return DivisorCoords::deg0(yn1, yn0, 0);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 - up1 - u0;
            let vpp0 = -upp0 * t2 - vp0;
            super::branch("ADD32");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1;
        let upp1 = s0 * cc.half - w4 * t2 - u0;
        let upp0 = (vp1 + vh1) * cc.half + w4 * (vp0 + yn0 + t2 * up1) - u0 * upp1;
        let vpp1 = upp1 * s0 - vh1;
        let vpp0 = s0 * (upp0 - up0) - vp0;
        super::branch("ADD33");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = -vp1 - yn1;
    let sp0 = v0 - vp0 - u0 * (yn1 - vp1 - t1.double());
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("ADD34");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z0.inv();
        let upp0 = -w2 * (vp0 + yn0) - up1 - u0;
        let vpp0 = -vp0 - upp0 * z0;
        super::branch("ADD35");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1;
    let upp1 = s0 * cc.half - w4 * z0 - u0;
    let upp0 = (vp1 + vh1) * cc.half + w4 * (vp0 + yn0 + z0 * up1) - u0 * upp1;
    let vpp1 = upp1 * s0 - vh1;
    let vpp0 = s0 * (upp0 - up0) - vp0;
    super::branch("ADD36");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-2 divisors. (`Deg2ADD`)
#[inline]
fn deg2_add_neg<F: Field>(
    u0: F, u1: F, v0: F, v1: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (yn0, yn1, f2) = (cc.yn0, cc.yn1, cc.f2);
    let m3 = up1 - u1;
    let m4 = u0 - up0;
    let m1 = m4 + up1 * m3;
    let m2 = -up0 * m3;
    let d = m1 * m4 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            // u = up
            let t0 = u1.double();
            let dw21 = vp1 + v1 + u0.double() - u1 * t0;
            let dw20 = vp0 + v0 - u0 * t0;
            if dw20.is_zero() && dw21.is_zero() {
                super::branch("ADD37");
                return DivisorCoords::deg0(yn1, yn0, 1);
            }
            let t2 = v1 - yn1;
            let t3 = v0 - yn0;
            let k2 = t2.double();
            let k1 = t3.double() - u1 * k2;
            let k0 = f2 - v1.square() - u0 * k2 - u1 * k1;
            let b2 = dw21.inv();
            let u0 = u1 - dw20 * b2;
            let s0 = b2 * (k0 - u0 * (k1 - u0 * k2));
            let upp1 = u0.double();
            let upp0 = u0.square();
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0 * s0;
            super::branch("ADD38");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let m2sq = m3.square();
        let m3cu = -m3 * m2sq;
        let t0 = m4.square();
        let t1 = vp1 + v1;
        let dw3 = m3cu * (vp0 + v0) - m4 * (m2sq * t1 - t0.double());
        if dw3.is_zero() {
            let a1 = -m3.inv();
            let s1cap = m4 * a1;
            let u0 = u1 - s1cap;
            let up0 = up1 - s1cap;
            let s0 = a1 * (vp0 - v0 - up0 * (vp1 - v1));
            let upp1 = u0 + up0;
            let upp0 = u0 * up0;
            let vpp1 = v1 + s0;
            let vpp0 = v0 + s0 * u0;
            super::branch("ADD39");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t2 = v1 - yn1;
        let t3 = v0 - yn0;
        let k2 = t2.double();
        let k1 = t3.double() - u1 * k2;
        let k0 = f2 - v1.square() - u0 * k2 - u1 * k1;
        let t2b = up0 - up1.square();
        let a12 = -m2sq * (t1 + t2b.double());
        let sp1 = a12 * (vp1 - v1) + m3cu * (k1 - up1 * k2);
        let sp0 = a12 * (vp0 - v0) + m3cu * (k0 - up0 * k2);
        let dd = dw3.square();
        if sp1.is_zero() {
            if sp0.is_zero() {
                super::branch("ADD40");
                return DivisorCoords::deg0(yn1, yn0, 2);
            }
            let w3 = (dw3 * sp0).inv();
            let s0 = sp0.square() * w3;
            let w4 = dd * w3;
            let t0 = s0 * u1;
            let upp0 = w4 * (v1 - yn1) - s0 * cc.half - up1;
            let vpp0 = upp0 * (t0 + yn1 + v1 - upp0 * (s0 + upp0.double())) - v0 - s0 * u0;
            super::branch("ADD41");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let w1 = dd - (sp1 - dw3).square();
        if w1.is_zero() {
            let t1 = u1.double();
            let w0 = sp0 + dw3 * t1;
            if w0.is_zero() {
                super::branch("ADD42");
                return DivisorCoords::deg0(yn1, yn0, 0);
            }
            let w2 = (dw3 * w0).inv();
            let s0 = sp0 * w0 * w2;
            let w3 = dd * w2;
            let t2 = s0 * u1 + u0.double() + v1 + yn1;
            let upp0 = w3 * t2 + s0 * cc.half - up1;
            let vpp0 = upp0 * (t2 - upp0 * (s0 + t1)) - v0 - s0 * u0;
            super::branch("ADD43");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let w2 = (dw3 * w1).inv();
        let w3 = w2 * w1;
        let w4 = w2 * dd * dw3;
        let s1 = sp1 * w3;
        let s0 = sp0 * w3;
        let l0 = s0 * u0;
        let t1 = s1 * u1;
        let l2 = s0 + t1;
        let l1 = (s0 + s1) * (u0 + u1) - l0 - t1;
        let vh1 = l1 + v1;
        let upp1 = w4 * (s0.double() - s1 * (l2 + s0)) - up1;
        let t0 = v1 - yn1;
        let upp0 = w4 * (t0.double() - s0 * l2 - s1 * (v1 + vh1)) - up0 - up1 * upp1;
        let t0 = s1 - (F::one() + F::one());
        let t1 = upp1 * t0;
        let t2 = l2 - t1;
        let t3 = upp0 * t2;
        let vpp1 = (upp0 + upp1) * (t0 + t2) - vh1 - t1 - t3;
        let vpp0 = t3 - v0 - l0;
        super::branch("ADD44");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    // d != 0 (frequent case)
    let r0 = vp0 - v0;
    let r1 = vp1 - v1;
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;
    let dd = d.square();
    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("ADD45");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let t0 = s0 * u1;
        let upp0 = w4 * (v1 - yn1) - s0 * cc.half - up1;
        let vpp0 = upp0 * (t0 + yn1 + v1 - upp0 * (s0 + upp0.double())) - v0 - s0 * u0;
        super::branch("ADD46");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w1 = dd - (sp1 - d).square();
    if w1.is_zero() {
        let t1 = u1.double();
        let w0 = sp0 + d * t1;
        if w0.is_zero() {
            super::branch("ADD47");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t2 = s0 * u1 + u0.double() + v1 + yn1;
        let upp0 = w3 * t2 + s0 * cc.half - up1;
        let vpp0 = upp0 * (t2 - upp0 * (s0 + t1)) - v0 - s0 * u0;
        super::branch("ADD48");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }
    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * dd * d;
    let s1 = sp1 * w3;
    let s0 = sp0 * w3;
    let l0 = s0 * u0;
    let t1 = s1 * u1;
    let l2 = s0 + t1;
    let t4 = (s0 + s1) * (u0 + u1) - l0 - t1 + v1;
    let upp1 = w4 * (s0.double() - s1 * (l2 + s0)) - up1;
    let t0 = v1 - yn1;
    let upp0 = w4 * (t0.double() - s0 * l2 - s1 * (v1 + t4)) - up0 - up1 * upp1;
    let t0 = s1 - (F::one() + F::one());
    let t1 = upp1 * t0;
    let t2 = l2 - t1;
    let t3 = upp0 * t2;
    let vpp1 = (upp0 + upp1) * (t0 + t2) - t4 - t1 - t3;
    let vpp0 = t3 - v0 - l0;
    super::branch("ADD49");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Add two divisors in the **negative** reduced basis. (`ADD` dispatcher)
#[inline]
pub fn add_neg<F: Field>(
    d1: &DivisorCoords<F>,
    d2: &DivisorCoords<F>,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let n = d1.n + d2.n - 1;
    match (d1.degree(), d2.degree()) {
        (2, 2) => deg2_add_neg(d1.u0, d1.u1, d1.v0, d1.v1, d2.u0, d2.u1, d2.v0, d2.v1, cc),
        (2, 1) => {
            if d2.n == 0 {
                deg12_add_up_neg(d2.u0, d2.v0, d1.u0, d1.u1, d1.v0, d1.v1, cc)
            } else {
                deg12_add_neg(d2.u0, d2.v0, d1.u0, d1.u1, d1.v0, d1.v1, cc)
            }
        }
        (2, 0) => match d2.n {
            0 => deg02_add_up_neg(d1.u0, d1.u1, d1.v0, d1.v1, cc),
            1 => {
                super::branch("ADD50");
                DivisorCoords { n: 0, ..*d1 }
            }
            _ => deg02_add_dwn_neg(d1.u0, d1.u1, d1.v0, d1.v1, cc),
        },
        (1, 2) => {
            if d1.n == 0 {
                deg12_add_up_neg(d1.u0, d1.v0, d2.u0, d2.u1, d2.v0, d2.v1, cc)
            } else {
                deg12_add_neg(d1.u0, d1.v0, d2.u0, d2.u1, d2.v0, d2.v1, cc)
            }
        }
        (1, 1) => match n {
            1 => deg1_add_dwn_neg(d1.u0, d1.v0, d2.u0, d2.v0, cc),
            -1 => deg1_add_up_neg(d1.u0, d1.v0, d2.u0, d2.v0, cc),
            _ => deg1_add_neg(d1.u0, d1.v0, d2.u0, d2.v0, cc),
        },
        (1, 0) => match d2.n {
            0 => {
                if d1.n == 0 {
                    deg01_add_up_neg(d1.u0, d1.v0, cc)
                } else {
                    super::branch("ADD51");
                    DivisorCoords { n: 0, ..*d1 }
                }
            }
            1 => {
                super::branch("ADD52");
                *d1
            }
            _ => {
                if d1.n == 0 {
                    super::branch("ADD53");
                    DivisorCoords { n: 1, ..*d1 }
                } else {
                    deg01_add_dwn_neg(d1.u0, d1.v0, cc)
                }
            }
        },
        (0, 2) => match d1.n {
            0 => deg02_add_up_neg(d2.u0, d2.u1, d2.v0, d2.v1, cc),
            1 => {
                super::branch("ADD54");
                DivisorCoords { n: 0, ..*d2 }
            }
            _ => deg02_add_dwn_neg(d2.u0, d2.u1, d2.v0, d2.v1, cc),
        },
        (0, 1) => match d1.n {
            0 => {
                if d2.n == 0 {
                    deg01_add_up_neg(d2.u0, d2.v0, cc)
                } else {
                    super::branch("ADD55");
                    DivisorCoords { n: 0, ..*d2 }
                }
            }
            1 => {
                super::branch("ADD56");
                *d2
            }
            _ => {
                if d2.n == 0 {
                    super::branch("ADD57");
                    DivisorCoords { n: 1, ..*d2 }
                } else {
                    deg01_add_dwn_neg(d2.u0, d2.v0, cc)
                }
            }
        },
        _ => {
            super::branch("ADD58");
            DivisorCoords::deg0(cc.yn1, cc.yn0, n)
        }
    }
}

// ===========================================================================
// Doubling — positive reduced basis
// Port of `nch2_splitG2_DBL.mag` (posReduced). Uses Vpl coeffs (y0,y1); the
// positive-basis divisor's v has leading coeff +1 (Vpl). Labels `PDBLnn`.
// ===========================================================================

/// Double a degree-1 divisor `<x + u0, v0, n=1>` (DWN). (`Deg1DBLDWN`, pos)
#[inline]
fn deg1_dbl_dwn_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f1, f4, y0, y1, d2) = (cc.f1, cc.f4, cc.y0, cc.y1, cc.half_d1);
    let u2 = u0.square();
    let d = v0.double() - u0 * (f4 + u2.double());
    if d.is_zero() {
        super::branch("PDBL00");
        return DivisorCoords::deg0(y1, y0, 2);
    }
    let z0 = y0 - v0;
    let t1 = u0 * z0;
    let z1 = d2 - t1;
    let t2 = u0 * (z1.double() - t1);
    let sp0 = f1 - f4 * v0 - t2.double();
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("PDBL01");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        let upp2 = upp0.square();
        let vpp0 = upp0 * (f4 + upp2.double()) - v0;
        super::branch("PDBL02");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = -z0 * w2 - u0;
    let upp0 = y1 + s0 * cc.half - z1 * w2 - u0 * upp1;
    let t0 = upp1.double();
    let vpp1 = upp0.double() - upp1 * t0 - y1 - s0;
    let vpp0 = -v0 - s0 * u0 - upp0 * t0;
    super::branch("PDBL03");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-1 divisor `<x + u0, v0, n=0>` (UP). (`Deg1DBLUP`, pos)
#[inline]
fn deg1_dbl_up_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f1, f4, y0, y1, d2) = (cc.f1, cc.f4, cc.y0, cc.y1, cc.half_d1);
    let u2 = u0.square();
    let t0 = u0 * (f4 + u2.double());
    let v0 = v0 - t0;
    let d = v0.double() + t0;
    if d.is_zero() {
        super::branch("PDBL04");
        return DivisorCoords::deg0(y1, y0, 0);
    }
    let z0 = y0 + v0;
    let t1 = u0 * z0;
    let z1 = d2 - t1;
    let t2 = u0 * (z1.double() - t1);
    let sp0 = f1 + f4 * v0 - t2.double();
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("PDBL05");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        super::branch("PDBL06");
        return DivisorCoords::deg1(upp0, y1, -v0, 1);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = z0 * w2 - u0;
    let upp0 = y1 + z1 * w2 - s0 * cc.half - u0 * upp1;
    let vpp1 = y1 - s0;
    let vpp0 = -v0 - s0 * u0;
    super::branch("PDBL07");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-2 divisor `<x² + u1x + u0, v1x + v0, n=0>`. (`Deg2DBL`, pos)
#[inline]
fn deg2_dbl_pos<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f2, y0, y1) = (cc.f2, cc.y0, cc.y1);
    let d4 = cc.f2 * cc.half; // f2/2
    let t0 = u1.square();
    let u12 = u1.double();
    let t1 = u0 - t0;
    let t2 = t1 - v1;
    let m3 = t2.double();
    let m4 = v0.double() + u0 * u12;
    let m1 = m4 + m3 * u1;
    let m2 = -m3 * u0;
    let d = m4 * m1 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            super::branch("PDBL08");
            return DivisorCoords::deg0(y1, y0, 1);
        }
        let b1 = -m3.inv();
        let t2b = y1 - v1;
        let t3 = y0 - v0;
        let k2 = t2b.double();
        let k1 = t3.double() - u1 * k2;
        let k0 = f2 - v1.square() - u0 * k2 - u1 * k1;
        let u0 = u1 - m4 * b1;
        let s0 = b1 * (k0 - u0 * (k1 - u0 * k2));
        let upp1 = u0.double();
        let upp0 = u0.square();
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0 * s0;
        super::branch("PDBL09");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let z0 = y0 - v0;
    let z1 = y1 - v1;
    let r1 = z0 - z1 * u12;
    let r0 = d4 - v1.square() * cc.half - u1 * z0 - z1 * (t1 + u0);
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("PDBL10");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = (d * sp0).inv();
        let s0 = sp0.square() * w2;
        let w3 = dd * w2;
        let t1b = z1 * cc.half;
        let upp0 = w3 * ((s0 - u1) * s0 - t1b);
        let t0 = upp0 * (s0 * u1 + t1b + v1 - upp0 * (s0 - upp0)) - s0 * u0;
        let vpp0 = t0.double() - v0;
        super::branch("PDBL11");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }

    let w1 = sp1 * (sp1 + d);
    if w1.is_zero() {
        let w0 = -sp0 - sp1 * u1;
        if w0.is_zero() {
            super::branch("PDBL12");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let w4 = s0 - u1;
        let t3 = z1 * cc.half;
        let upp0 = w3 * (w4 * s0 + t2 - t3);
        let t0 = upp0 * (t3 + v1 + s0 * u1 - u0 - upp0 * w4) - s0 * u0;
        let vpp0 = t0.double() - v0;
        super::branch("PDBL13");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }

    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * d * dd;
    let s0 = w3 * sp0;
    let s1 = w3 * sp1;
    let t3 = s0 - u1;
    let upp1 = w4 * ((s0 + t3) * s1 + s0);
    let upp0 = w4 * (t3 * s0 - t2 * s1 - z1 * cc.half);
    let zz0 = upp0 - u0;
    let zz1 = upp1 - u1;
    let ww0 = zz0 * s0;
    let ww1 = zz1 * s1;
    let ww2 = ww1 + upp1;
    let t1b = (s0 + s1) * (zz0 + zz1) - ww0 - ww1 + upp0 - ww2 * upp1;
    let t0 = ww0 - ww2 * upp0;
    let vpp1 = t1b.double() - v1;
    let vpp0 = t0.double() - v0;
    super::branch("PDBL14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a divisor in the **positive** reduced basis. (`DBL` dispatcher, pos)
#[inline]
pub fn double_pos<F: Field>(d: &DivisorCoords<F>, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    match d.degree() {
        2 => deg2_dbl_pos(d.u1, d.u0, d.v1, d.v0, cc),
        1 => {
            if d.n == 0 {
                deg1_dbl_up_pos(d.u0, d.v0, cc)
            } else {
                deg1_dbl_dwn_pos(d.u0, d.v0, cc)
            }
        }
        _ => {
            if d.n == 1 {
                super::branch("PDBL15");
                DivisorCoords::deg0(cc.y1, cc.y0, 1)
            } else if d.n == 0 {
                super::branch("PDBL16");
                // <adu, Vpl, 2 − audeg>
                DivisorCoords { u2: cc.au2, u1: cc.au1, u0: cc.au0, v1: cc.y1, v0: cc.y0, n: 2 - cc.audeg as i32 }
            } else {
                super::branch("PDBL17");
                // <adu, adv_pos, 0>
                DivisorCoords { u2: cc.au2, u1: cc.au1, u0: cc.au0, v1: cc.adv1_pos, v0: cc.adv0_pos, n: 0 }
            }
        }
    }
}

// ===========================================================================
// Addition — positive reduced basis
// Port of `nch2_splitG2_ADD.mag` (posReduced). Uses Vpl coeffs (y0,y1) and pos
// d-constants (d1 = f2−y1² = cc.d1, d2 = d1/2 = half_d1, d3 = f1/2 = half_f1,
// d4 = f2/2 = f2·half). Labels `PADDnn`.
// ===========================================================================

/// `<1, V, 2>` + degree-1 `<x+u0, v0, 0>`. (`Deg01ADDDWN`, pos)
#[inline]
fn deg01_add_dwn_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1, d1, d2, d3, f1, f4) = (cc.y0, cc.y1, cc.d1, cc.half_d1, cc.half_f1, cc.f1, cc.f4);
    let z0 = y0 - v0;
    if z0.is_zero() {
        if d1.is_zero() {
            super::branch("PADD00");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = d1.inv();
        let up0 = w2 * (f1 - v0 * f4) - u0;
        let t1 = up0.square();
        let vp0 = up0 * (f4 + t1.double()) - v0;
        super::branch("PADD01");
        return DivisorCoords::deg1(up0, y1, vp0, 0);
    }
    let w2 = z0.inv();
    let up1 = w2 * d2 - u0;
    let up0 = w2 * (d3 - v0 * y1) - u0 * up1;
    let t0 = up0 - up1.square();
    let t1 = up1.double();
    let vp1 = t0.double() - y1;
    let vp0 = -up0 * t1 - v0;
    super::branch("PADD02");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 0>` + degree-1 `<x+u0, v0, 0>`. (`Deg01ADDUP`, pos)
#[inline]
fn deg01_add_up_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1, d1, d2, d3, f1, f4) = (cc.y0, cc.y1, cc.d1, cc.half_d1, cc.half_f1, cc.f1, cc.f4);
    let t0 = u0.square();
    let vp0 = u0 * (f4 + t0.double()) - v0;
    let z0 = y0 - vp0;
    if z0.is_zero() {
        if d1.is_zero() {
            super::branch("PADD03");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = d1.inv();
        let up0 = w2 * (f1 - vp0 * f4) - u0;
        super::branch("PADD04");
        return DivisorCoords::deg1(up0, y1, vp0, 1);
    }
    let w2 = z0.inv();
    let up1 = w2 * d2 - u0;
    let up0 = w2 * (d3 - vp0 * y1) - u0 * up1;
    super::branch("PADD05");
    DivisorCoords::deg2(up1, up0, y1, vp0, 0)
}

/// `<1, V, 2>` + degree-2 divisor. (`Deg02ADDDWN`, pos)
#[inline]
fn deg02_add_dwn_pos<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1, f4) = (cc.y0, cc.y1, cc.f4);
    let d4 = cc.f2 * cc.half;
    let z0 = y0 - v0;
    let z1 = y1 - v1;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("PADD06");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = z0.inv();
        let up0 = w2 * (d4 - v1.square() * cc.half) - u1;
        let t0 = up0.square();
        let vp0 = up0 * (f4 + t0.double()) - v0;
        super::branch("PADD07");
        return DivisorCoords::deg1(up0, y1, vp0, 0);
    }
    let w2 = z1.inv();
    let up1 = w2 * z0 - u1;
    let up0 = w2 * (d4 - v1.square() * cc.half) - u0 - u1 * up1;
    let t0 = up0 - up1.square();
    let t1 = up1.double();
    let vp1 = t0.double() - v1;
    let vp0 = -v0 - t1 * up0;
    super::branch("PADD08");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 0>` + degree-2 divisor. (`Deg02ADDUP`, pos)
#[inline]
fn deg02_add_up_pos<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1) = (cc.y0, cc.y1);
    let d4 = cc.f2 * cc.half;
    let t1 = u1.double();
    let t2 = u1.square() - u0;
    let v1 = v1 + t2.double();
    let v0 = v0 + u0 * t1;
    let z0 = v0 + y0;
    let z1 = v1 + y1;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("PADD09");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = z0.inv();
        let up0 = w2 * (d4 - v1.square() * cc.half) - u1;
        super::branch("PADD10");
        return DivisorCoords::deg1(up0, y1, -v0, 1);
    }
    let w2 = z1.inv();
    let up1 = w2 * z0 - u1;
    let up0 = w2 * (d4 - v1.square() * cc.half) - u0 - u1 * up1;
    super::branch("PADD11");
    DivisorCoords::deg2(up1, up0, -v1, -v0, 0)
}

/// Two degree-1 divisors with `n + np = 1`. (`Deg1ADD`, pos)
#[inline]
fn deg1_add_pos<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1, d1, f1, f4) = (cc.y0, cc.y1, cc.d1, cc.f1, cc.f4);
    let d = u0 - up0;
    if d.is_zero() {
        let upp1 = u0.double();
        let upp0 = u0.square();
        let dw = v0 + vp0 - u0 * f4 - upp0 * upp1;
        if dw.is_zero() {
            super::branch("PADD12");
            return DivisorCoords::deg0(y1, y0, 1);
        }
        let t0 = y0 - v0;
        let t2 = f1 - f4 * v0 - u0 * (d1.double() - upp1 * (t0 + t0 + t0));
        let s0 = t2 * dw.inv();
        let vpp1 = y1 + s0;
        let vpp0 = v0 + s0 * u0;
        super::branch("PADD13");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let s0 = (vp0 - v0) * d.inv();
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;
    let vpp1 = y1 + s0;
    let vpp0 = v0 + u0 * s0;
    super::branch("PADD14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 1` (DWN). (`Deg1ADDDWN`, pos)
#[inline]
fn deg1_add_dwn_pos<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1, d2, f4) = (cc.y0, cc.y1, cc.half_d1, cc.f4);
    let d = u0 - up0;
    if d.is_zero() {
        super::branch("PADD15");
        return DivisorCoords::deg0(y1, y0, 2);
    }
    let sp0 = vp0 - v0;
    let z0 = y0 - v0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("PADD16");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w1 = z0.inv();
        let upp0 = w1 * (d2 - u0 * z0) - up0;
        let t1 = upp0.square();
        let vpp0 = upp0 * (f4 + t1.double()) - v0;
        super::branch("PADD17");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t1 = z0 * w1;
    let upp1 = -t1 - up0;
    let upp0 = y1 + s0 * cc.half - w1 * d2 + u0 * t1 - up0 * upp1;
    let t0 = -upp1.double();
    let vpp1 = -y1 - s0 + upp0.double() + upp1 * t0;
    let vpp0 = -v0 - s0 * u0 + upp0 * t0;
    super::branch("PADD18");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 0` (UP). (`Deg1ADDUP`, pos)
#[inline]
fn deg1_add_up_pos<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (y0, y1, d2, f4) = (cc.y0, cc.y1, cc.half_d1, cc.f4);
    let d = u0 - up0;
    if d.is_zero() {
        super::branch("PADD19");
        return DivisorCoords::deg0(y1, y0, 0);
    }
    let t0 = u0.square();
    let vt0 = v0 - u0 * (t0.double() + f4);
    let t1 = up0.square();
    let sp0 = vp0 - vt0 - up0 * (t1.double() + f4);
    let z0 = y0 + vt0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("PADD20");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w1 = z0.inv();
        let upp0 = w1 * (d2 - u0 * z0) - up0;
        super::branch("PADD21");
        return DivisorCoords::deg1(upp0, y1, -vt0, 1);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = z0 * w1;
    let upp1 = t0 - up0;
    let upp0 = y1 - s0 * cc.half + w1 * d2 - u0 * t0 - up0 * upp1;
    let vpp1 = y1 - s0;
    let vpp0 = -vt0 - s0 * u0;
    super::branch("PADD22");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,1>` + degree-2 (DWN). (`Deg12ADD`, pos)
#[inline]
fn deg12_add_pos<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, y0, y1) = (cc.f2, cc.y0, cc.y1);
    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 - t0 + t1;
    if d.is_zero() {
        let dw = v0 + vp0 - u0 * (vp1 + y1 + t1.double());
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let vpp0 = vp0 + upp0 * (y1 - vp1);
            super::branch("PADD23");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let t2 = y1 - vp1;
        let t3 = y0 - vp0;
        let k2 = t2.double();
        let k1 = t3.double() - up1 * k2;
        let k0 = f2 - vp1.square() - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("PADD24");
                return DivisorCoords::deg0(y1, y0, 0);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 - up1 - u0;
            let vpp0 = upp0 * (y1 + vp1 + upp0 * upp0.double()) - vp0;
            super::branch("PADD25");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let t0 = s0 * up1 + vp1;
        let upp1 = s0 * cc.half - w4 * t2 - u0;
        let upp0 = (t0 + vp1) * cc.half - w4 * (y0 - vp0 - t2 * up1) - u0 * upp1;
        let t1 = s0 - upp1.double();
        let vpp1 = upp1 * t1 + upp0.double() - t0;
        let vpp0 = upp0 * t1 - s0 * up0 - vp0;
        super::branch("PADD26");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = y1 - vp1;
    let sp0 = v0 - vp0 - z0 * u0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("PADD27");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (y0 - vp0) - up1 - u0;
        let t1 = upp0.square();
        let vpp0 = upp0 * (y1 + vp1 + t1.double()) - vp0;
        super::branch("PADD28");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1;
    let upp1 = s0 * cc.half - w4 * z0 - u0;
    let upp0 = (vh1 + vp1) * cc.half - w4 * (y0 - vp0 - z0 * up1) - u0 * upp1;
    let t1 = s0 - upp1.double();
    let vpp1 = upp0.double() + t1 * upp1 - vh1;
    let vpp0 = upp0 * t1 - s0 * up0 - vp0;
    super::branch("PADD29");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,0>` + degree-2 (UP). (`Deg12ADDUP`, pos)
#[inline]
fn deg12_add_up_pos<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, y0, y1) = (cc.f2, cc.y0, cc.y1);
    // vp := -V - h - ((-V-h - vp) mod up)
    let t2 = up1.double();
    let vp1 = vp1 - up0.double() + t2 * up1;
    let vp0 = vp0 + up0 * t2;
    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 - t0 + t1;
    if d.is_zero() {
        let dw = vp0 + v0 - u0 * (y1 + vp1);
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let t0 = upp0.square();
            let vpp0 = vp0 + upp0 * (y1 - vp1 + t0.double());
            super::branch("PADD30");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let t2 = vp1 + y1;
        let t3 = vp0 + y0;
        let k2 = t2.double();
        let k1 = t3.double() - up1 * k2;
        let k0 = f2 - vp1.square() - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("PADD31");
                return DivisorCoords::deg0(y1, y0, 2);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 - up1 - u0;
            let vpp0 = upp0 * t2 - vp0;
            super::branch("PADD32");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let t0 = s0 * up1;
        let upp1 = w4 * t2 - s0 * cc.half - u0;
        let upp0 = w4 * (t3 - t2 * up1) - vp1 - t0 * cc.half - u0 * upp1;
        let vpp1 = upp1 * s0 - t0 - vp1;
        let vpp0 = s0 * (upp0 - up0) - vp0;
        super::branch("PADD33");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let k2 = y1 + vp1;
    let sp0 = v0 - vp0 - u0 * (y1 - vp1 + t1.double());
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("PADD34");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = k2.inv();
        let upp0 = w2 * (y0 + vp0) - up1 - u0;
        let vpp0 = upp0 * k2 - vp0;
        super::branch("PADD35");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1;
    let upp1 = w4 * k2 - s0 * cc.half - u0;
    let upp0 = w4 * (y0 + vp0 - k2 * up1) - (vp1 + vh1) * cc.half - u0 * upp1;
    let vpp1 = upp1 * s0 - vh1;
    let vpp0 = s0 * (upp0 - up0) - vp0;
    super::branch("PADD36");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-2 divisors. (`Deg2ADD`, pos)
#[inline]
fn deg2_add_pos<F: Field>(
    u0: F, u1: F, v0: F, v1: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, y0, y1) = (cc.f2, cc.y0, cc.y1);
    let m3 = up1 - u1;
    let m4 = u0 - up0;
    let m1 = m4 + up1 * m3;
    let m2 = -up0 * m3;
    let d = m1 * m4 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            let t0 = u1.double();
            let dw21 = vp1 + v1 - u0.double() + u1 * t0;
            let dw20 = vp0 + v0 + u0 * t0;
            if dw20.is_zero() && dw21.is_zero() {
                super::branch("PADD37");
                return DivisorCoords::deg0(y1, y0, 1);
            }
            let t2 = y1 - v1;
            let t3 = y0 - v0;
            let k2 = t2.double();
            let k1 = t3.double() - u1 * k2;
            let k0 = f2 - v1.square() - u0 * k2 - u1 * k1;
            let b2 = dw21.inv();
            let u0 = u1 - dw20 * b2;
            let s0 = b2 * (k0 - u0 * (k1 - u0 * k2));
            let upp1 = u0.double();
            let upp0 = u0.square();
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0 * s0;
            super::branch("PADD38");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let m2sq = m3.square();
        let m3cu = -m3 * m2sq;
        let t0 = m4.square();
        let t1 = vp1 + v1;
        let dw3 = m3cu * (vp0 + v0) - m4 * (m2sq * t1 + t0.double());
        if dw3.is_zero() {
            let a1 = -m3.inv();
            let s1cap = m4 * a1;
            let u0 = u1 - s1cap;
            let up0 = up1 - s1cap;
            let s0 = a1 * (vp0 - v0 - up0 * (vp1 - v1));
            let upp1 = u0 + up0;
            let upp0 = u0 * up0;
            let vpp1 = v1 + s0;
            let vpp0 = v0 + s0 * u0;
            super::branch("PADD39");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t2 = y1 - v1;
        let t3 = y0 - v0;
        let k2 = t2.double();
        let k1 = t3.double() - u1 * k2;
        let k0 = f2 - v1.square() - u0 * k2 - u1 * k1;
        let t2b = up1.square() - up0;
        let a12 = -m2sq * (t1 + t2b.double());
        let sp1 = a12 * (vp1 - v1) + m3cu * (k1 - up1 * k2);
        let sp0 = a12 * (vp0 - v0) + m3cu * (k0 - up0 * k2);
        let dd = dw3.square();
        if sp1.is_zero() {
            if sp0.is_zero() {
                super::branch("PADD40");
                return DivisorCoords::deg0(y1, y0, 0);
            }
            let w2 = (dw3 * sp0).inv();
            let s0 = sp0.square() * w2;
            let w3 = dd * w2;
            let t0 = s0 * u1;
            let upp0 = w3 * (s0.square() * cc.half - y1 + v1) - up1;
            let vpp0 = -v0 - s0 * u0 + upp0 * (t0 + y1 + v1 - upp0 * (s0 - upp0.double()));
            super::branch("PADD41");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let w1 = dd - (sp1 + dw3).square();
        if w1.is_zero() {
            let t1 = -u1.double();
            let t2a = sp0 + dw3 * t1;
            let w0 = t2a.double();
            if w0.is_zero() {
                super::branch("PADD42");
                return DivisorCoords::deg0(y1, y0, 2);
            }
            let w2 = (dw3 * w0).inv();
            let s0 = sp0 * w0 * w2;
            let w3 = dd * w2;
            let t2 = s0 * u1 - u0.double() + v1 + y1;
            let t3 = s0 + t1;
            let upp0 = w3 * (t2.double() - s0 * t3) - up1;
            let vpp0 = upp0 * (t2 - upp0 * t3) - v0 - s0 * u0;
            super::branch("PADD43");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let w2 = (dw3 * w1).inv();
        let w3 = w2 * w1;
        let w4 = w2 * dd * dw3;
        let s1 = sp1 * w3;
        let s0 = sp0 * w3;
        let l0 = s0 * u0;
        let t1 = s1 * u1;
        let l2 = s0 + t1;
        let l1 = (s0 + s1) * (u0 + u1) - l0 - t1;
        let upp1 = -w4 * (s1 * (l2 + s0) + s0.double()) - up1;
        let t0 = y1 - v1;
        let t4 = l1 + v1;
        let upp0 = w4 * (t0.double() - s0 * l2 - s1 * (t4 + v1)) - up0 - up1 * upp1;
        let t0 = (F::one() + F::one()) + s1;
        let t1 = upp1 * t0;
        let t2 = l2 - t1;
        let t3 = upp0 * t2;
        let vpp1 = (upp0 + upp1) * (t0 + t2) - t4 - t1 - t3;
        let vpp0 = t3 - v0 - l0;
        super::branch("PADD44");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let r0 = vp0 - v0;
    let r1 = vp1 - v1;
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;
    let dd = d.square();
    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("PADD45");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = (d * sp0).inv();
        let s0 = sp0.square() * w2;
        let w3 = dd * w2;
        let t0 = s0 * u1;
        let t1 = y1 - v1;
        let upp0 = w3 * (s0.square() * cc.half - t1) - up1;
        let vpp0 = -v0 - s0 * u0 + upp0 * (t0 + y1 + v1 - upp0 * (s0 - upp0.double()));
        super::branch("PADD46");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w1 = dd - (sp1 + d).square();
    if w1.is_zero() {
        let t1 = -u1.double();
        let t2a = sp0 + d * t1;
        let w0 = t2a.double();
        if w0.is_zero() {
            super::branch("PADD47");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t0 = s0 * u1;
        let t2 = t0 - u0.double() + v1 + y1;
        let t3 = s0 + t1;
        let upp0 = w3 * (t2.double() - s0 * t3) - up1;
        let vpp0 = upp0 * (t2 - upp0 * t3) - v0 - s0 * u0;
        super::branch("PADD48");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }
    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * dd * d;
    let s1 = sp1 * w3;
    let s0 = sp0 * w3;
    let l0 = s0 * u0;
    let t1 = s1 * u1;
    let l2 = s0 + t1;
    let l1 = (s0 + s1) * (u0 + u1) - l0 - t1;
    let upp1 = -w4 * (s1 * (l2 + s0) + s0.double()) - up1;
    let t0 = y1 - v1;
    let t4 = l1 + v1;
    let upp0 = w4 * (t0.double() - s0 * l2 - s1 * (t4 + v1)) - up0 - up1 * upp1;
    let t0 = (F::one() + F::one()) + s1;
    let t1 = upp1 * t0;
    let t2 = l2 - t1;
    let t3 = upp0 * t2;
    let vpp1 = (upp0 + upp1) * (t0 + t2) - t4 - t1 - t3;
    let vpp0 = t3 - v0 - l0;
    super::branch("PADD49");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Add two divisors in the **positive** reduced basis. (`ADD` dispatcher, pos)
#[inline]
pub fn add_pos<F: Field>(
    d1: &DivisorCoords<F>,
    d2: &DivisorCoords<F>,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    match (d1.degree(), d2.degree()) {
        (2, 2) => deg2_add_pos(d1.u0, d1.u1, d1.v0, d1.v1, d2.u0, d2.u1, d2.v0, d2.v1, cc),
        (2, 1) => {
            if d2.n == 0 {
                deg12_add_up_pos(d2.u0, d2.v0, d1.u0, d1.u1, d1.v0, d1.v1, cc)
            } else {
                deg12_add_pos(d2.u0, d2.v0, d1.u0, d1.u1, d1.v0, d1.v1, cc)
            }
        }
        (2, 0) => match d2.n {
            0 => deg02_add_up_pos(d1.u0, d1.u1, d1.v0, d1.v1, cc),
            1 => {
                super::branch("PADD50");
                DivisorCoords { n: 0, ..*d1 }
            }
            _ => deg02_add_dwn_pos(d1.u0, d1.u1, d1.v0, d1.v1, cc),
        },
        (1, 2) => {
            if d1.n == 0 {
                deg12_add_up_pos(d1.u0, d1.v0, d2.u0, d2.u1, d2.v0, d2.v1, cc)
            } else {
                deg12_add_pos(d1.u0, d1.v0, d2.u0, d2.u1, d2.v0, d2.v1, cc)
            }
        }
        (1, 1) => match d1.n + d2.n {
            0 => deg1_add_up_pos(d1.u0, d1.v0, d2.u0, d2.v0, cc),
            2 => deg1_add_dwn_pos(d1.u0, d1.v0, d2.u0, d2.v0, cc),
            _ => deg1_add_pos(d1.u0, d1.v0, d2.u0, d2.v0, cc),
        },
        (1, 0) => match d2.n {
            0 => {
                if d1.n == 0 {
                    deg01_add_up_pos(d1.u0, d1.v0, cc)
                } else {
                    super::branch("PADD51");
                    DivisorCoords { n: 0, ..*d1 }
                }
            }
            1 => {
                super::branch("PADD52");
                *d1
            }
            _ => {
                if d1.n == 0 {
                    super::branch("PADD53");
                    DivisorCoords { n: 1, ..*d1 }
                } else {
                    deg01_add_dwn_pos(d1.u0, d1.v0, cc)
                }
            }
        },
        (0, 2) => match d1.n {
            0 => deg02_add_up_pos(d2.u0, d2.u1, d2.v0, d2.v1, cc),
            1 => {
                super::branch("PADD54");
                DivisorCoords { n: 0, ..*d2 }
            }
            _ => deg02_add_dwn_pos(d2.u0, d2.u1, d2.v0, d2.v1, cc),
        },
        (0, 1) => match d1.n {
            0 => {
                if d2.n == 0 {
                    deg01_add_up_pos(d2.u0, d2.v0, cc)
                } else {
                    super::branch("PADD55");
                    DivisorCoords { n: 0, ..*d2 }
                }
            }
            1 => {
                super::branch("PADD56");
                *d2
            }
            _ => {
                if d2.n == 0 {
                    super::branch("PADD57");
                    DivisorCoords { n: 1, ..*d2 }
                } else {
                    deg01_add_dwn_pos(d2.u0, d2.v0, cc)
                }
            }
        },
        _ => {
            super::branch("PADD58");
            DivisorCoords::deg0(cc.y1, cc.y0, d1.n + d2.n - 1)
        }
    }
}
