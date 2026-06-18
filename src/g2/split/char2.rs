//! Characteristic-2 genus 2 **split** model divisor arithmetic.
//!
//! Explicit formulas for split genus 2 hyperelliptic curves over `GF(2^k)`,
//! in the normalized form `h(x) = x³ + h1·x + h0`, `f(x) = f6·x⁶ + f2·x² + f1·x + f0`.
//! These use XOR-based char-2 arithmetic (no division by 2).
//!
//! Divisor representation and balance weight `n` are as in [`super::not_char2`].
//! Reduced-basis polynomials: `Vpl = x³ + y1·x + y0` (`y3=1`-ish; actually `y3`
//! is a root of `x²+x+f6`), `Vn = Vpl + h` (`= −Vpl − h` in char 2).
//!
//! Ported from the Magma `ch2_splitG2_{UTL,DBL,ADD}.mag` (neg + pos reduced).

// Explicit divisor formulas take the divisor's coordinates as separate scalars.
#![allow(clippy::too_many_arguments)]

use crate::field::Field;
use crate::poly::Poly;

pub use super::not_char2::{DivisorCoords, G};

/// Curve constants and precomputations for a char-2 split genus 2 curve.
/// For the **positive** basis, `c0 = yn0` and `c1 = yn1` (`= y_i + h_i`).
#[derive(Clone, Copy, Debug)]
pub struct CurveConstants<F: Field> {
    pub f0: F,
    pub f1: F,
    pub f2: F,
    pub f6: F,
    pub h0: F,
    pub h1: F,
    // Vpl = y3 x³ + y1 x + y0 (y2 = 0); Vn = yn3 x³ + yn1 x + yn0.
    pub y0: F,
    pub y1: F,
    pub y3: F,
    pub yn0: F,
    pub yn1: F,
    pub yn3: F,
    pub d0: F,
    pub d1: F,
    pub au0: F,
    pub au1: F,
    pub au2: F,
    pub audeg: i8,
    pub adv0n: F,
    pub adv1n: F,
    pub adv0p: F,
    pub adv1p: F,
}

impl<F: Field> CurveConstants<F> {
    pub fn f_poly(&self) -> Poly<F> {
        Poly::from_coeffs(vec![
            self.f0,
            self.f1,
            self.f2,
            F::zero(),
            F::zero(),
            F::zero(),
            self.f6,
        ])
    }
    /// `h = x³ + h1·x + h0` (h3 = 1, h2 = 0).
    pub fn h_poly(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.h0, self.h1, F::zero(), F::one()])
    }
    /// Positive reduced-basis polynomial `Vpl = y3 x³ + y1 x + y0`.
    pub fn vpl(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.y0, self.y1, F::zero(), self.y3])
    }
    /// Negative reduced-basis polynomial `Vn = −Vpl − h = yn3 x³ + yn1 x + yn0`.
    pub fn vn(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.yn0, self.yn1, F::zero(), self.yn3])
    }
}

/// Build [`CurveConstants`] from the char-2 curve coefficients and a chosen
/// root `y3` of `x² + x + f6` (`y3² + y3 + f6 = 0`; the two roots `y3, y3+1`
/// are the points at infinity). Port of the Magma char-2 `Precompute`.
pub fn precompute<F: Field>(f0: F, f1: F, f2: F, f6: F, h0: F, h1: F, y3: F) -> CurveConstants<F> {
    let _ = f6; // f6 only constrains y3 (already supplied)
    let one = F::one();
    let yn3 = y3 + one; // −y3 − h3 with h3 = 1
    let y1 = h1 * y3;
    let yn1 = y1 + h1;
    let y0 = h0 * y3;
    let yn0 = y0 + h0;
    let d1 = f2 + y1 * yn1;
    let d0 = f1 + h0 * y1;

    let k1 = f1 + yn1 * y0 + yn0 * y1;
    let k0 = f0 + y0 * yn0;

    let (au0, au1, au2, audeg, adv0n, adv1n, adv0p, adv1p);
    if d1.is_zero() {
        if k1.is_zero() {
            au2 = F::zero();
            au1 = F::zero();
            au0 = one;
            audeg = 0;
            adv1n = yn1;
            adv0n = yn0;
            adv1p = y1;
            adv0p = y0;
        } else {
            let w1 = k1.inv();
            au2 = F::zero();
            au1 = one;
            au0 = k0 * w1;
            audeg = 1;
            adv1n = yn1;
            adv0n = y0 + au0 * (h1 + au0.square());
            adv1p = y1;
            adv0p = yn0 + au0 * (h1 + au0.square()); // c0 = y0 + h0 = yn0
        }
    } else {
        let w1 = d1.inv();
        au2 = one;
        au1 = k1 * w1;
        au0 = k0 * w1;
        audeg = 2;
        adv1n = y1 + au0 + au1.square();
        adv0n = y0 + au0 * au1;
        adv1p = yn1 + au0 + au1.square(); // c1 = y1 + h1 = yn1
        adv0p = yn0 + au0 * au1; // c0 = yn0
    }

    CurveConstants {
        f0,
        f1,
        f2,
        f6,
        h0,
        h1,
        y0,
        y1,
        y3,
        yn0,
        yn1,
        yn3,
        d0,
        d1,
        au0,
        au1,
        au2,
        audeg,
        adv0n,
        adv1n,
        adv0p,
        adv1p,
    }
}

// ===========================================================================
// Doubling — negative reduced basis. Port of `ch2_splitG2_DBL.mag` (neg).
// Labels `CDBLnn`. In char 2: `-x == x`, `x - y == x + y`, `x + x == 0`.
// Neg c-constants are inlined: c3 = 1, c2 = 0, c1 = h1, c0 = h0.
// ===========================================================================

/// Double a degree-1 divisor `<x + u0, v0, n=1>` (DWN). (`Deg1DBLDWN`)
#[inline]
fn deg1_dbl_dwn_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1) = (cc.h0, cc.h1, cc.yn0, cc.yn1);
    let (d0, d1) = (cc.d0, cc.d1);

    // v := -V-h - ((-V-h - v) mod u);
    let t0 = u0.square();
    let t1 = u0 * (h1 + t0);
    let vp0 = v0 + t1;
    // d := (2*v + h) mod u;
    let vh0 = vp0 + h0;
    let d = t1 + h0;
    if d.is_zero() {
        super::branch("CDBL00");
        return DivisorCoords::deg0(yn1, yn0, 2);
    }
    // k := ExactQuotient(f - v*(v + h), u);
    let z0 = vh0 + yn0;
    let z1 = d1 + z0 * u0;
    let sp0 = d0 + h1 * vp0 + t0 * z0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("CDBL01");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 + u0;
        super::branch("CDBL02");
        return DivisorCoords::deg1(upp0, yn1, vh0, 0);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = z0 * w2 + u0;
    let upp0 = h1 + s0 + z1 * w2 + u0 * upp1;
    let vpp1 = yn1 + s0;
    let vpp0 = vh0 + s0 * u0;
    super::branch("CDBL03");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-1 divisor `<x + u0, v0, n=0>` (UP). (`Deg1DBLUP`)
#[inline]
fn deg1_dbl_up_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1, y1) = (cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);

    // d := (v + v + h) mod u;
    let t0 = u0.square();
    let d = h0 + u0 * (h1 + t0);
    if d.is_zero() {
        super::branch("CDBL04");
        return DivisorCoords::deg0(yn1, yn0, 0);
    }
    let vh0 = v0 + h0;
    let z0 = v0 + yn0;
    let z1 = d1 + z0 * u0;
    let sp0 = d0 + h1 * vh0 + t0 * z0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("CDBL05");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 + u0;
        let vpp0 = vh0 + upp0 * (h1 + upp0.square());
        super::branch("CDBL06");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = z0 * w2 + u0;
    let upp0 = h1 + s0 + z1 * w2 + u0 * upp1;
    let vpp1 = y1 + s0 + upp0 + upp1.square();
    let vpp0 = vh0 + s0 * u0 + upp0 * upp1;
    super::branch("CDBL07");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-2 divisor `<x² + u1x + u0, v1x + v0, n=0>`. (`Deg2DBL`)
#[inline]
fn deg2_dbl_neg<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f2, h0, h1, yn0, yn1, y1) = (cc.f2, cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);

    let t0 = v0 + h0;
    let t1 = v1 + h1;
    let t2 = u1.square();
    let m3 = t2 + u0 + h1;
    let m4 = h0 + u0 * u1;
    let m1 = m4 + m3 * u1;
    let m2 = m3 * u0;
    let d = m4 * m1 + m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            super::branch("CDBL08");
            return DivisorCoords::deg0(yn1, yn0, 1);
        }
        let b1 = m3.inv();
        let k2 = v1 + yn1;
        let t3 = v0 + yn0;
        let k1 = t3 + u1 * k2;
        let k0 = f2 + v1 * t1 + u0 * k2 + u1 * k1;
        let u0 = u1 + m4 * b1;
        let s0 = b1 * (k0 + u0 * (k1 + u0 * k2));
        let upp1 = u0.double(); // u² : upp1 = u0+u0 = 0 in char 2
        let upp0 = u0.square();
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0 * s0;
        super::branch("CDBL09");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    // Compute s(d)/c3
    let k2p = v1 + yn1;
    let r1 = v0 + yn0;
    let r0 = f2 + v1 * t1 + u1 * r1 + k2p * t2;
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("CDBL10");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let upp0 = s0 + u1 + w4 * k2p;
        let vpp0 = upp0 * (s0 * (upp0 + u1) + v1 + y1 + upp0.square()) + s0 * u0 + t0;
        super::branch("CDBL11");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }

    let w1 = sp1 * (sp1 + d);
    if w1.is_zero() {
        let w0 = sp0 + sp1 * u1;
        if w0.is_zero() {
            super::branch("CDBL12");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t1 = s0 * u1;
        let upp0 = w3 * (s0.square() + t1 + m3 + k2p);
        let vpp0 = t0 + s0 * u0 + upp0 * (v1 + y1 + t1 + u0 + upp0 * (s0 + u1));
        super::branch("CDBL13");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }

    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * d * dd;
    let s0 = w3 * sp0;
    let s1 = w3 * sp1;
    let upp1 = w4 * (u1 * s1 + s0);
    let upp0 = w4 * (s0 * (s0 + u1) + m3 * s1 + k2p);
    let z0 = upp0 + u0;
    let z1 = upp1 + u1;
    let w0 = z0 * s0;
    let w1 = z1 * s1;
    let w2 = upp1 + w1;
    // `t1` here is the original `t1 = v1 + h1` (Magma never reassigns it on this
    // path; the inner `t1 := s0*u1` lives in returned branches).
    let vpp1 = t1 + (s0 + s1) * (z0 + z1) + w0 + w1 + upp0 + w2 * upp1;
    let vpp0 = t0 + w0 + w2 * upp0;
    super::branch("CDBL14");
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
            if d.n == 1 {
                super::branch("CDBL15");
                DivisorCoords::deg0(cc.yn1, cc.yn0, 1)
            } else if d.n == 0 {
                super::branch("CDBL16");
                DivisorCoords {
                    u2: cc.au2,
                    u1: cc.au1,
                    u0: cc.au0,
                    v1: cc.adv1n,
                    v0: cc.adv0n,
                    n: 2 - cc.audeg as i32,
                }
            } else {
                super::branch("CDBL17");
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
// Addition — negative reduced basis. Port of `ch2_splitG2_ADD.mag` (neg).
// Labels `CADDnn`. `(u0[,u1],v0[,v1])` is the first operand,
// `(up0[,up1],vp0[,vp1])` the second. Neg c-constants inlined (c1 = h1, etc.).
// ===========================================================================

/// `<1, V, 0>` + degree-1 divisor `<x + u0, v0, 1>`. (`Deg01ADDUP`)
#[inline]
fn deg01_add_up_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1, y1) = (cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);

    let z3 = v0 + yn0;
    if z3.is_zero() {
        let z2 = d1;
        if z2.is_zero() {
            super::branch("CADD00");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z2.inv();
        let vh0 = v0 + h0;
        let up0 = w2 * (d0 + vh0 * h1) + u0;
        let vp0 = up0 * (up0.square() + h1) + vh0;
        super::branch("CADD01");
        return DivisorCoords::deg1(up0, yn1, vp0, 1);
    }
    let w2 = z3.inv();
    let vh0 = v0 + h0;
    let up1 = w2 * d1 + u0;
    let up0 = w2 * (d0 + vh0 * h1) + u0 * up1;
    let vp1 = y1 + up1.square() + up0;
    let vp0 = vh0 + up1 * up0;
    super::branch("CADD02");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 2>` + degree-1 divisor `<x + u0, v0, 0>`. (`Deg01ADDDWN`)
#[inline]
fn deg01_add_dwn_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1) = (cc.h0, cc.h1, cc.yn0, cc.yn1);
    let (d0, d1) = (cc.d0, cc.d1);

    let v0 = v0 + u0 * (h1 + u0.square());
    let vp0 = v0 + h0;
    let z3 = vp0 + yn0;
    if z3.is_zero() {
        let z2 = d1;
        if z2.is_zero() {
            super::branch("CADD03");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z2.inv();
        let up0 = w2 * (d0 + v0 * h1) + u0;
        super::branch("CADD04");
        return DivisorCoords::deg1(up0, yn1, vp0, 0);
    }
    let w2 = z3.inv();
    let up1 = w2 * d1 + u0;
    let up0 = w2 * (d0 + v0 * h1) + u0 * up1;
    super::branch("CADD05");
    DivisorCoords::deg2(up1, up0, yn1, vp0, 0)
}

/// `<1, V, 0>` + degree-2 divisor. (`Deg02ADDUP`)
#[inline]
fn deg02_add_up_neg<F: Field>(
    u0: F,
    u1: F,
    v0: F,
    v1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, yn0, yn1) = (cc.f2, cc.h0, cc.h1, cc.yn0, cc.yn1);

    let z1 = v1 + yn1;
    let z0 = v0 + yn0;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("CADD06");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z0.inv();
        let up0 = w2 * (f2 + v1 * (h1 + v1)) + u1;
        let vp0 = up0 * (up0.square() + h1) + v0 + h0;
        super::branch("CADD07");
        return DivisorCoords::deg1(up0, yn1, vp0, 1);
    }
    let w2 = z1.inv();
    let t1 = h1 + v1;
    let t2 = w2 * z0;
    let up1 = t2 + u1;
    let up0 = w2 * (f2 + v1 * t1) + u0 + u1 * up1;
    let vp1 = up1.square() + up0 + t1;
    let vp0 = up1 * up0 + h0 + v0;
    super::branch("CADD08");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 2>` + degree-2 divisor. (`Deg02ADDDWN`)
#[inline]
fn deg02_add_dwn_neg<F: Field>(
    u0: F,
    u1: F,
    v0: F,
    v1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, yn0, yn1) = (cc.f2, cc.h0, cc.h1, cc.yn0, cc.yn1);

    let v1 = v1 + u1.square() + u0;
    let v0 = v0 + u1 * u0;
    let vp1 = v1 + h1;
    let vp0 = v0 + h0;
    let z1 = vp1 + yn1;
    let z0 = vp0 + yn0;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("CADD09");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z0.inv();
        let up0 = w2 * (f2 + v1 * vp1) + u1;
        super::branch("CADD10");
        return DivisorCoords::deg1(up0, yn1, vp0, 0);
    }
    let w2 = z1.inv();
    let t2 = w2 * z0;
    let up1 = t2 + u1;
    let up0 = w2 * (f2 + v1 * vp1) + u0 + up1 * u1;
    super::branch("CADD11");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// Two degree-1 divisors with `n + np = 1`. (`Deg1ADD`)
#[inline]
fn deg1_add_neg<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1) = (cc.h0, cc.h1, cc.yn0, cc.yn1);
    let d0 = cc.d0;

    let d = u0 + up0;
    if d.is_zero() {
        let upp0 = u0.square();
        let vh0 = v0 + h0;
        let dw = vh0 + vp0 + u0 * (h1 + upp0);
        if dw.is_zero() {
            super::branch("CADD12");
            return DivisorCoords::deg0(yn1, yn0, 1);
        }
        let t2 = d0 + h1 * vh0 + upp0 * (v0 + yn0);
        let s0 = t2 * dw.inv();
        let vpp1 = yn1 + s0;
        let vpp0 = v0 + s0 * u0;
        super::branch("CADD13");
        return DivisorCoords::deg2(u0.double(), upp0, vpp1, vpp0, 0);
    }
    let s0 = (vp0 + v0) * d.inv();
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;
    let vpp1 = yn1 + s0;
    let vpp0 = v0 + u0 * s0;
    super::branch("CADD14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 0` (UP). (`Deg1ADDUP`)
#[inline]
fn deg1_add_up_neg<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1, y1) = (cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);
    // Magma `Deg1ADDUP` declares d0 but doesn't use it; only d1 is needed.
    let d1 = cc.d1;

    let d = u0 + up0;
    if d.is_zero() {
        super::branch("CADD15");
        return DivisorCoords::deg0(yn1, yn0, 0);
    }
    let sp0 = vp0 + v0;
    let z0 = v0 + yn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("CADD16");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w1 = z0.inv();
        let upp0 = w1 * (d1 + z0 * u0) + up0;
        let vpp0 = upp0 * (upp0.square() + h1) + v0 + h0;
        super::branch("CADD17");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t1 = z0 * w1;
    let upp1 = t1 + up0;
    let upp0 = h1 + w1 * d1 + t1 * u0 + s0 + up0 * upp1;
    let vpp1 = y1 + s0 + upp1.square() + upp0;
    let vpp0 = v0 + h0 + s0 * u0 + upp0 * upp1;
    super::branch("CADD18");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 1` (DWN). (`Deg1ADDDWN`)
#[inline]
fn deg1_add_dwn_neg<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1) = (cc.h0, cc.h1, cc.yn0, cc.yn1);
    let d1 = cc.d1;

    let d = u0 + up0;
    if d.is_zero() {
        super::branch("CADD19");
        return DivisorCoords::deg0(yn1, yn0, 2);
    }
    let vt0 = v0 + u0 * (h1 + u0.square());
    let sp0 = vp0 + vt0 + up0 * (h1 + up0.square());
    let vh0 = vt0 + h0;
    let z0 = vh0 + yn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("CADD20");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w1 = z0.inv();
        let upp0 = w1 * d1 + u0 + up0;
        super::branch("CADD21");
        return DivisorCoords::deg1(upp0, yn1, vh0, 0);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = z0 * w1;
    let upp1 = t0 + up0;
    let upp0 = h1 + s0 + w1 * d1 + t0 * u0 + up0 * upp1;
    let vpp1 = yn1 + s0;
    let vpp0 = vh0 + s0 * u0;
    super::branch("CADD22");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,1>` + degree-2 (UP). (`Deg12ADDUP`)
#[inline]
fn deg12_add_up_neg<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    up1: F,
    vp0: F,
    vp1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, yn0, yn1, y1) = (cc.f2, cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);

    let t0 = u0 * up1;
    let d = up0 + t0 + u0.square();
    if d.is_zero() {
        let dw = v0 + vp0 + h0 + u0 * (vp1 + y1 + u0.square());
        if dw.is_zero() {
            let upp0 = up1 + u0;
            let vpp0 = vp0 + upp0 * (yn1 + vp1);
            super::branch("CADD23");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let k2 = vp1 + yn1;
        let t3 = vp0 + yn0;
        let k1 = t3 + up1 * k2;
        let k0 = f2 + vp1 * (vp1 + h1) + up0 * k2 + up1 * k1;
        let sp0 = k0 + u0 * (k1 + u0 * k2);
        if sp0.is_zero() {
            if k2.is_zero() {
                super::branch("CADD24");
                return DivisorCoords::deg0(yn1, yn0, 2);
            }
            let w2 = k2.inv();
            let upp0 = w2 * t3 + up1 + u0;
            let vpp0 = upp0 * (vp1 + y1 + upp0.square()) + vp0 + h0;
            super::branch("CADD25");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1 + h1;
        let upp1 = w4 * k2 + s0 + u0;
        let upp0 = w4 * (vp0 + yn0 + k2 * up1) + vh1 + vp1 + u0 * upp1;
        let t1 = s0 + upp1;
        let vpp1 = upp1 * t1 + upp0 + vh1;
        let vpp0 = upp0 * t1 + s0 * up0 + vp0 + h0;
        super::branch("CADD26");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = vp1 + yn1;
    let sp0 = v0 + vp0 + z0 * u0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("CADD27");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (vp0 + yn0) + up1 + u0;
        let vpp0 = upp0 * (vp1 + y1 + upp0.square()) + vp0 + h0;
        super::branch("CADD28");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1 + h1;
    let upp1 = w4 * z0 + s0 + u0;
    let upp0 = w4 * (vp0 + yn0 + z0 * up1) + vh1 + vp1 + u0 * upp1;
    let t1 = s0 + upp1;
    let vpp1 = upp1 * t1 + upp0 + vh1;
    let vpp0 = upp0 * t1 + s0 * up0 + vp0 + h0;
    super::branch("CADD29");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,0>` + degree-2 (DWN). (`Deg12ADD`)
#[inline]
fn deg12_add_neg<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    up1: F,
    vp0: F,
    vp1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, yn0, yn1, y1) = (cc.f2, cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);

    // vp := -V-h - ((-V-h - vp) mod up)
    let vp1 = vp1 + up0 + up1.square();
    let vp0 = vp0 + up0 * up1;

    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 + t0 + t1;
    if d.is_zero() {
        let dw = vp0 + v0 + h0 + u0 * (y1 + vp1);
        if dw.is_zero() {
            let upp0 = up1 + u0;
            let vpp0 = vp0 + upp0 * (yn1 + vp1 + upp0.square());
            super::branch("CADD30");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let vh1 = vp1 + h1;
        let vh0 = vp0 + h0;
        let k2 = vh1 + yn1;
        let t3 = vh0 + yn0;
        let k1 = t3 + up1 * k2;
        let k0 = f2 + vp1 * vh1 + up0 * k2 + up1 * k1;
        let sp0 = k0 + u0 * (k1 + u0 * k2);
        if sp0.is_zero() {
            if k2.is_zero() {
                super::branch("CADD31");
                return DivisorCoords::deg0(yn1, yn0, 0);
            }
            let w2 = k2.inv();
            let upp0 = w2 * t3 + up1 + u0;
            let vpp0 = upp0 * k2 + vh0;
            super::branch("CADD32");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vh1;
        let upp1 = s0 + w4 * k2 + u0;
        let upp0 = vp1 + vh1 + w4 * (vh0 + yn0 + k2 * up1) + u0 * upp1;
        let vpp1 = upp1 * s0 + vh1;
        let vpp0 = s0 * (upp0 + up0) + vh0;
        super::branch("CADD33");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let vh1 = vp1 + h1;
    let vh0 = vp0 + h0;
    let z0 = vh1 + yn1;
    let sp0 = v0 + vp0 + u0 * (yn1 + vp1 + t1);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("CADD34");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (vh0 + yn0) + up1 + u0;
        let vpp0 = vh0 + upp0 * z0;
        super::branch("CADD35");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vh1;
    let upp1 = s0 + w4 * z0 + u0;
    let upp0 = vp1 + vh1 + w4 * (vh0 + yn0 + z0 * up1) + u0 * upp1;
    let vpp1 = upp1 * s0 + vh1;
    let vpp0 = s0 * (upp0 + up0) + h0 + vp0;
    super::branch("CADD36");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-2 divisors. (`Deg2ADD`)
#[inline]
fn deg2_add_neg<F: Field>(
    u0: F,
    u1: F,
    v0: F,
    v1: F,
    up0: F,
    up1: F,
    vp0: F,
    vp1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, yn0, yn1, y1) = (cc.f2, cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);

    let m3 = up1 + u1;
    let m4 = u0 + up0;
    let m1 = m4 + up1 * m3;
    let m2 = up0 * m3;
    let d = m1 * m4 + m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            // u = up
            let t1 = v1 + h1;
            let dw21 = vp1 + t1 + u1.square() + u0;
            let dw20 = vp0 + v0 + h0 + u0 * u1;
            if dw20.is_zero() && dw21.is_zero() {
                super::branch("CADD37");
                return DivisorCoords::deg0(yn1, yn0, 1);
            }
            let k2 = v1 + yn1;
            let t3 = v0 + yn0;
            let k1 = t3 + u1 * k2;
            let k0 = f2 + v1 * t1 + u0 * k2 + u1 * k1;
            let b2 = dw21.inv();
            let u0 = u1 + dw20 * b2;
            let s0 = b2 * (k0 + u0 * (k1 + u0 * k2));
            let upp1 = u0.double(); // u² : upp1 = u0+u0 = 0 in char 2
            let upp0 = u0.square();
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0 * s0;
            super::branch("CADD38");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t0 = v0 + h0;
        let vh1 = v1 + h1;
        let m2sq = m3.square();
        let m3cu = m3 * m2sq;
        let dw3 = m3cu * (vp0 + t0) + m4 * (m2sq * (vp1 + vh1) + m4.square());
        if dw3.is_zero() {
            let a1 = m3.inv();
            let s1cap = m4 * a1;
            let u0 = u1 + s1cap;
            let up0 = up1 + s1cap;
            let s0 = a1 * (vp0 + v0 + up0 * (vp1 + v1));
            let upp1 = u0 + up0;
            let upp0 = u0 * up0;
            let vpp1 = v1 + s0;
            let vpp0 = v0 + s0 * u0;
            super::branch("CADD39");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let k2 = v1 + yn1;
        let t3 = v0 + yn0;
        let k1 = t3 + u1 * k2;
        let k0 = f2 + v1 * vh1 + u0 * k2 + u1 * k1;
        let a12 = m2sq * (vp1 + vh1 + up0 + up1.square());
        let sp1 = a12 * (vp1 + v1) + m3cu * (k1 + up1 * k2);
        let sp0 = a12 * (vp0 + v0) + m3cu * (k0 + up0 * k2);
        let dd = dw3.square();
        if sp1.is_zero() {
            if sp0.is_zero() {
                super::branch("CADD40");
                return DivisorCoords::deg0(yn1, yn0, 2);
            }
            let w3 = (dw3 * sp0).inv();
            let s0 = sp0.square() * w3;
            let w4 = dd * w3;
            let upp0 = s0 + w4 * (v1 + yn1) + up1;
            let vpp0 = upp0 * (s0 * (u1 + upp0) + y1 + v1 + upp0.square()) + v0 + h0 + s0 * u0;
            super::branch("CADD41");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let w1 = sp1 * (sp1 + dw3);
        if w1.is_zero() {
            let w0 = sp0 + dw3 * u1;
            if w0.is_zero() {
                super::branch("CADD42");
                return DivisorCoords::deg0(yn1, yn0, 0);
            }
            let w2 = (dw3 * w0).inv();
            let s0 = sp0 * w0 * w2;
            let w3 = dd * w2;
            let t2 = s0 * u1 + u0 + v1 + y1;
            let upp0 = w3 * t2 + s0 + up1;
            let vpp0 = upp0 * (t2 + upp0 * (s0 + u1)) + v0 + h0 + s0 * u0;
            super::branch("CADD43");
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
        let l1 = (s0 + s1) * (u0 + u1) + l0 + t1;
        let t4 = l1 + h1;
        let upp1 = w4 * (s1 * t1 + s0) + up1;
        let upp0 = w4 * (s0 * l2 + s1 * t4 + v1 + yn1) + up0 + up1 * upp1;
        let t1 = upp1 + s1 * upp1;
        let t2 = l2 + t1;
        let t3 = upp0 * t2;
        let vpp1 = (upp0 + upp1) * (s1 + F::one() + t2) + v1 + t4 + t1 + t3;
        let vpp0 = t3 + v0 + h0 + l0;
        super::branch("CADD44");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let r0 = vp0 + v0;
    let r1 = vp1 + v1;
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("CADD45");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let upp0 = s0 + w4 * (v1 + yn1) + up1;
        let vpp0 = upp0 * (s0 * u1 + y1 + v1 + upp0 * (s0 + upp0)) + v0 + h0 + s0 * u0;
        super::branch("CADD46");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w1 = sp1 * (sp1 + d);
    if w1.is_zero() {
        let w0 = sp0 + d * u1;
        if w0.is_zero() {
            super::branch("CADD47");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t2 = s0 * u1 + u0 + v1 + y1;
        let upp0 = w3 * t2 + s0 + up1;
        let vpp0 = upp0 * (t2 + upp0 * (s0 + u1)) + v0 + h0 + s0 * u0;
        super::branch("CADD48");
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
    let l1 = (s0 + s1) * (u0 + u1) + l0 + t1;
    let t4 = l1 + h1;
    let upp1 = w4 * (s1 * t1 + s0) + up1;
    let upp0 = w4 * (s0 * l2 + s1 * t4 + v1 + yn1) + up0 + up1 * upp1;
    let t1 = upp1 + s1 * upp1;
    let t2 = l2 + t1;
    let t3 = upp0 * t2;
    let vpp1 = (upp0 + upp1) * (s1 + F::one() + t2) + v1 + t4 + t1 + t3;
    let vpp0 = t3 + v0 + h0 + l0;
    super::branch("CADD49");
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
                super::branch("CADD50");
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
                    super::branch("CADD51");
                    DivisorCoords { n: 0, ..*d1 }
                }
            }
            1 => {
                super::branch("CADD52");
                *d1
            }
            _ => {
                if d1.n == 0 {
                    super::branch("CADD53");
                    DivisorCoords { n: 1, ..*d1 }
                } else {
                    deg01_add_dwn_neg(d1.u0, d1.v0, cc)
                }
            }
        },
        (0, 2) => match d1.n {
            0 => deg02_add_up_neg(d2.u0, d2.u1, d2.v0, d2.v1, cc),
            1 => {
                super::branch("CADD54");
                DivisorCoords { n: 0, ..*d2 }
            }
            _ => deg02_add_dwn_neg(d2.u0, d2.u1, d2.v0, d2.v1, cc),
        },
        (0, 1) => match d1.n {
            0 => {
                if d2.n == 0 {
                    deg01_add_up_neg(d2.u0, d2.v0, cc)
                } else {
                    super::branch("CADD55");
                    DivisorCoords { n: 0, ..*d2 }
                }
            }
            1 => {
                super::branch("CADD56");
                *d2
            }
            _ => {
                if d2.n == 0 {
                    super::branch("CADD57");
                    DivisorCoords { n: 1, ..*d2 }
                } else {
                    deg01_add_dwn_neg(d2.u0, d2.v0, cc)
                }
            }
        },
        _ => {
            super::branch("CADD58");
            DivisorCoords::deg0(cc.yn1, cc.yn0, n)
        }
    }
}

// ===========================================================================
// Doubling — positive reduced basis. Port of `ch2_splitG2_DBL.mag` (pos).
// Labels `CPDBLnn`. Flat-packed `ccs`: c0 = yn0, c1 = yn1, c5 = h1.
// ===========================================================================

/// Double a degree-1 divisor `<x + u0, v0, n=1>` (DWN). (`Deg1DBLDWN`, pos)
#[inline]
fn deg1_dbl_dwn_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);
    let c1 = cc.yn1; // pos c1 = y1 + h1 = yn1

    // d := h mod u;
    let d = h0 + u0 * (h1 + u0.square());
    if d.is_zero() {
        super::branch("CPDBL00");
        return DivisorCoords::deg0(y1, y0, 2);
    }
    // k := ExactQuotient(f - v*h - v^2, u);
    let k2 = y0 + v0;
    let t0 = u0 * k2;
    let k1 = d1 + t0;
    let sp0 = d0 + h1 * v0 + u0 * t0;
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("CPDBL01");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w1 = k2.inv();
        let upp0 = k1 * w1 + u0;
        let vpp0 = v0 + h0 + upp0 * (h1 + upp0.square());
        super::branch("CPDBL02");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = k2 * w2 + u0;
    let upp0 = h1 + s0 + k1 * w2 + u0 * upp1;
    let vpp1 = c1 + s0 + upp0 + upp1.square();
    let vpp0 = v0 + h0 + s0 * u0 + upp0 * upp1;
    super::branch("CPDBL03");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-1 divisor `<x + u0, v0, n=0>` (UP). (`Deg1DBLUP`, pos)
#[inline]
fn deg1_dbl_up_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);

    // vp := -V - h + (V + h + v) mod u;
    let t0 = u0 * (h1 + u0.square());
    let v0 = v0 + t0;
    let vp0 = v0 + h0;
    let d = h0 + t0;
    if d.is_zero() {
        super::branch("CPDBL04");
        return DivisorCoords::deg0(y1, y0, 0);
    }
    let k2 = y0 + vp0;
    let t0 = u0 * k2;
    let k1 = d1 + t0;
    let sp0 = d0 + h1 * vp0 + u0 * t0;
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("CPDBL05");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w1 = k2.inv();
        let upp0 = k1 * w1 + u0;
        super::branch("CPDBL06");
        return DivisorCoords::deg1(upp0, y1, vp0, 1);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = k2 * w2 + u0;
    let upp0 = h1 + s0 + k1 * w2 + u0 * upp1;
    let vpp1 = y1 + s0;
    let vpp0 = vp0 + s0 * u0;
    super::branch("CPDBL07");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-2 divisor `<x² + u1x + u0, v1x + v0, n=0>`. (`Deg2DBL`, pos)
#[inline]
fn deg2_dbl_pos<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (f2, h0, h1, y0, y1) = (cc.f2, cc.h0, cc.h1, cc.y0, cc.y1);
    let c1 = cc.yn1; // pos c1 = yn1

    let t0 = v0 + h0;
    let t1 = v1 + h1;
    let t2 = u1.square();
    let m3 = h1 + t2 + u0;
    let m4 = h0 + u0 * u1;
    let m1 = m4 + m3 * u1;
    let m2 = m3 * u0;
    let d = m4 * m1 + m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            super::branch("CPDBL08");
            return DivisorCoords::deg0(y1, y0, 1);
        }
        let b1 = m3.inv();
        let k2 = y1 + v1;
        let t3 = y0 + v0;
        let k1 = t3 + u1 * k2;
        let k0 = f2 + v1 * t1 + u0 * k2 + u1 * k1;
        let u0 = u1 + m4 * b1;
        let s0 = b1 * (k0 + u0 * (k1 + u0 * k2));
        let upp1 = u0.double(); // u² : x² + (u0+u0)x + u0²  (upp1 = 0 in char 2)
        let upp0 = u0.square();
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0 * s0;
        super::branch("CPDBL09");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let r1 = y0 + v0;
    let k2p = y1 + v1;
    let r0 = f2 + v1 * t1 + u1 * r1 + k2p * t2;
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("CPDBL10");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = (d * sp0).inv();
        let s0 = sp0.square() * w2;
        let w3 = dd * w2;
        let t1 = s0 * u1;
        let upp0 = w3 * (s0.square() + t1 + k2p);
        let vpp0 = t0 + s0 * u0 + upp0 * (t1 + c1 + v1 + upp0 * (s0 + upp0));
        super::branch("CPDBL11");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }

    let w1 = sp1 * (sp1 + d);
    if w1.is_zero() {
        let w0 = sp0 + sp1 * u1;
        if w0.is_zero() {
            super::branch("CPDBL12");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t1 = s0 * u1;
        let upp0 = w3 * (s0.square() + t1 + m3 + k2p);
        let vpp0 = t0 + s0 * u0 + upp0 * (c1 + v1 + t1 + u0 + upp0 * (s0 + u1));
        super::branch("CPDBL13");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }

    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * d * dd;
    let s0 = w3 * sp0;
    let s1 = w3 * sp1;
    let upp1 = w4 * (s1 * u1 + s0);
    let upp0 = w4 * (s0 * (s0 + u1) + m3 * s1 + k2p);
    let z0 = upp0 + u0;
    let z1 = upp1 + u1;
    let w0 = z0 * s0;
    let w1 = z1 * s1;
    let w2 = w1 + upp1;
    // `t1` here is the original `t1 = v1 + h1` (unchanged on this path).
    let vpp1 = t1 + (s0 + s1) * (z0 + z1) + w0 + w1 + upp0 + w2 * upp1;
    let vpp0 = t0 + w0 + w2 * upp0;
    super::branch("CPDBL14");
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
                super::branch("CPDBL15");
                DivisorCoords::deg0(cc.y1, cc.y0, 1)
            } else if d.n == 0 {
                super::branch("CPDBL16");
                DivisorCoords {
                    u2: cc.au2,
                    u1: cc.au1,
                    u0: cc.au0,
                    v1: cc.y1,
                    v0: cc.y0,
                    n: 2 - cc.audeg as i32,
                }
            } else {
                super::branch("CPDBL17");
                DivisorCoords {
                    u2: cc.au2,
                    u1: cc.au1,
                    u0: cc.au0,
                    v1: cc.adv1p,
                    v0: cc.adv0p,
                    n: 0,
                }
            }
        }
    }
}

// ===========================================================================
// Addition — positive reduced basis. Port of `ch2_splitG2_ADD.mag` (pos).
// Labels `CPADDnn`. Flat-packed `ccs`: c0 = yn0, c1 = yn1, c5 = h1.
// ===========================================================================

/// `<1, V, 2>` + degree-1 divisor `<x + u0, v0, 1>`. (`Deg01ADDDWN`, pos)
#[inline]
fn deg01_add_dwn_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);
    let c1 = cc.yn1;

    // z = (f - h*v - v^2)/c3; z1 := d1; z0 := y0 + v0;
    let z0 = y0 + v0;
    if z0.is_zero() {
        if d1.is_zero() {
            super::branch("CPADD00");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let up0 = d1.inv() * (d0 + v0 * h1) + u0;
        let vp0 = h0 + v0 + up0 * (h1 + up0.square());
        super::branch("CPADD01");
        return DivisorCoords::deg1(up0, y1, vp0, 0);
    }
    let w2 = z0.inv();
    let up1 = w2 * d1 + u0;
    let up0 = w2 * (d0 + v0 * h1) + u0 * up1;
    let vp1 = up1.square() + up0 + c1;
    let vp0 = up0 * up1 + v0 + h0;
    super::branch("CPADD02");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 0>` + degree-1 divisor `<x + u0, v0, 0>`. (`Deg01ADDUP`, pos)
#[inline]
fn deg01_add_up_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);
    let c0 = cc.yn0;

    // vp:= -V-h -((-V-h - v) mod u);
    let v0 = v0 + u0 * (h1 + u0.square());
    let vp0 = v0 + h0;
    let z0 = c0 + v0;
    if z0.is_zero() {
        if d1.is_zero() {
            super::branch("CPADD03");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = d1.inv();
        let up0 = w2 * (d0 + vp0 * h1) + u0;
        super::branch("CPADD04");
        return DivisorCoords::deg1(up0, y1, vp0, 1);
    }
    let w2 = z0.inv();
    let up1 = w2 * d1 + u0;
    let up0 = w2 * (d0 + vp0 * h1) + u0 * up1;
    super::branch("CPADD05");
    DivisorCoords::deg2(up1, up0, y1, vp0, 0)
}

/// `<1, V, 2>` + degree-2 divisor. (`Deg02ADDDWN`, pos)
#[inline]
fn deg02_add_dwn_pos<F: Field>(
    u0: F,
    u1: F,
    v0: F,
    v1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, y0, y1) = (cc.f2, cc.h0, cc.h1, cc.y0, cc.y1);

    let z0 = y0 + v0;
    let z1 = y1 + v1;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("CPADD06");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = z0.inv();
        let up0 = w2 * (f2 + v1 * (h1 + v1)) + u1;
        let vp0 = v0 + h0 + up0 * (h1 + up0.square());
        super::branch("CPADD07");
        return DivisorCoords::deg1(up0, y1, vp0, 0);
    }
    let t1 = h1 + v1;
    let w2 = z1.inv();
    let up1 = w2 * z0 + u1;
    let up0 = w2 * (f2 + v1 * t1) + u0 + u1 * up1;
    let vp1 = t1 + up0 + up1.square();
    let vp0 = h0 + v0 + up1 * up0;
    super::branch("CPADD08");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 0>` + degree-2 divisor. (`Deg02ADDUP`, pos)
#[inline]
fn deg02_add_up_pos<F: Field>(
    u0: F,
    u1: F,
    v0: F,
    v1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, y0, y1) = (cc.f2, cc.h0, cc.h1, cc.y0, cc.y1);
    let (c0, c1) = (cc.yn0, cc.yn1);

    // v1 := -V-h -(-V-h - v1) mod u1;
    let v1 = v1 + u0 + u1.square();
    let v0 = v0 + u1 * u0;
    let z0 = c0 + v0;
    let z1 = c1 + v1;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("CPADD09");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = z0.inv();
        let up0 = w2 * (f2 + v1 * (v1 + h1)) + u1;
        let vp0 = v0 + h0;
        super::branch("CPADD10");
        return DivisorCoords::deg1(up0, y1, vp0, 1);
    }
    let w2 = z1.inv();
    let vp1 = v1 + h1;
    let up1 = w2 * z0 + u1;
    let up0 = w2 * (f2 + v1 * vp1) + u0 + u1 * up1;
    let vp0 = v0 + h0;
    super::branch("CPADD11");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// Two degree-1 divisors with `n + np = 1`. (`Deg1ADD`, pos)
#[inline]
fn deg1_add_pos<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d0, d1) = (cc.d0, cc.d1);

    // d = upp1 = u mod up
    let upp1 = u0 + up0;
    let sp0 = vp0 + v0;
    if upp1.is_zero() {
        let upp0 = u0.square();
        let dw = v0 + vp0 + h0 + u0 * h1 + upp0 * u0;
        if dw.is_zero() {
            super::branch("CPADD12");
            return DivisorCoords::deg0(y1, y0, 1);
        }
        let upp1 = u0 + u0;
        let t0 = y0 + v0;
        let t1 = d0 + h1 * v0 + u0 * (d1 + d1 + u0 * (t0 + t0 + t0));
        let s0 = t1 * dw.inv();
        let vpp1 = y1 + s0;
        let vpp0 = v0 + s0 * u0;
        super::branch("CPADD13");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let s0 = sp0 * upp1.inv();
    let upp0 = u0 * up0;
    let vpp1 = y1 + s0;
    let vpp0 = v0 + u0 * s0;
    super::branch("CPADD14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 1` (DWN). (`Deg1ADDDWN`, pos)
#[inline]
fn deg1_add_dwn_pos<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let d1 = cc.d1;
    let c1 = cc.yn1;

    let d = u0 + up0;
    if d.is_zero() {
        super::branch("CPADD15");
        return DivisorCoords::deg0(y1, y0, 2);
    }
    let sp0 = vp0 + v0;
    let k2 = y0 + v0;
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("CPADD16");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w1 = k2.inv();
        let upp0 = w1 * (d1 + u0 * k2) + up0;
        let vpp0 = v0 + h0 + upp0 * (h1 + upp0.square());
        super::branch("CPADD17");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = k2 * w1;
    let upp1 = t0 + up0;
    let upp0 = h1 + s0 + w1 * d1 + t0 * u0 + up0 * upp1;
    let vpp1 = c1 + s0 + upp0 + upp1.square();
    let vpp0 = v0 + h0 + s0 * u0 + upp0 * upp1;
    super::branch("CPADD18");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 0` (UP). (`Deg1ADDUP`, pos)
#[inline]
fn deg1_add_up_pos<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let d1 = cc.d1;

    let d = u0 + up0;
    if d.is_zero() {
        super::branch("CPADD19");
        return DivisorCoords::deg0(y1, y0, 0);
    }
    let vt0 = v0 + u0 * (u0.square() + h1);
    let sp0 = vp0 + vt0 + up0 * (up0.square() + h1);
    let vn0 = vt0 + h0;
    let k2 = y0 + vn0;
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("CPADD20");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w1 = k2.inv();
        let upp0 = w1 * (d1 + u0 * k2) + up0;
        super::branch("CPADD21");
        return DivisorCoords::deg1(upp0, y1, vn0, 1);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = k2 * w1;
    let upp1 = t0 + up0;
    let upp0 = h1 + s0 + w1 * d1 + t0 * u0 + up0 * upp1;
    let vpp1 = y1 + s0;
    let vpp0 = vn0 + s0 * u0;
    super::branch("CPADD22");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,1>` + degree-2 (DWN). (`Deg12ADD`, pos)
#[inline]
fn deg12_add_pos<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    up1: F,
    vp0: F,
    vp1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, y0, y1) = (cc.f2, cc.h0, cc.h1, cc.y0, cc.y1);
    // Magma pos `Deg12ADD` declares c0 but doesn't use it; only c1 is needed.
    let c1 = cc.yn1;

    let t1 = u0.square();
    let d = up0 + u0 * up1 + t1;
    if d.is_zero() {
        let vh0 = vp0 + h0;
        let dw = v0 + vp0 + h0 + u0 * (vp1 + c1 + t1);
        if dw.is_zero() {
            let upp0 = up1 + u0;
            let vpp0 = vp0 + upp0 * (y1 + vp1);
            super::branch("CPADD23");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let vh1 = vp1 + h1;
        let k2 = y1 + vp1;
        let t2 = y0 + vp0;
        let k1 = t2 + up1 * k2;
        let k0 = f2 + vp1 * vh1 + up0 * k2 + up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if k2.is_zero() {
                super::branch("CPADD24");
                return DivisorCoords::deg0(y1, y0, 0);
            }
            let w2 = k2.inv();
            let upp0 = w2 * t2 + up1 + u0;
            let vpp0 = upp0 * (c1 + vp1 + upp0.square()) + vh0;
            super::branch("CPADD25");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vh1;
        let upp1 = s0 + w4 * k2 + u0;
        let upp0 = vh1 + vp1 + w4 * (y0 + vp0 + k2 * up1) + u0 * upp1;
        let t1 = s0 + upp1;
        let vpp1 = upp0 + t1 * upp1 + vh1;
        let vpp0 = upp0 * t1 + s0 * up0 + vh0;
        super::branch("CPADD26");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let k2 = y1 + vp1;
    let sp0 = v0 + vp0 + k2 * u0;
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("CPADD27");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = k2.inv();
        let upp0 = w2 * (y0 + vp0) + up1 + u0;
        let vpp0 = upp0 * (c1 + vp1 + upp0.square()) + vp0 + h0;
        super::branch("CPADD28");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let t0 = s0 * up0 + h0;
    let t1 = s0 * up1 + h1;
    let t2 = w4 * k2;
    let upp1 = s0 + t2 + u0;
    let upp0 = t1 + w4 * (y0 + vp0) + t2 * up1 + u0 * upp1;
    let t2 = s0 + upp1;
    let vpp1 = upp0 + upp1 * t2 + vp1 + t1;
    let vpp0 = upp0 * t2 + vp0 + t0;
    super::branch("CPADD29");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,0>` + degree-2 (UP). (`Deg12ADDUP`, pos)
#[inline]
fn deg12_add_up_pos<F: Field>(
    u0: F,
    v0: F,
    up0: F,
    up1: F,
    vp0: F,
    vp1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, y0, y1) = (cc.f2, cc.h0, cc.h1, cc.y0, cc.y1);
    let (c0, c1) = (cc.yn0, cc.yn1);

    // vp := -V-h - ((-V-h -vp) mod up);
    let vp1 = vp1 + up0 + up1.square();
    let vp0 = vp0 + up0 * up1;

    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 + t0 + t1;
    if d.is_zero() {
        let vh0 = vp0 + h0;
        let dw = vh0 + v0 + u0 * (c1 + vp1);
        if dw.is_zero() {
            let upp0 = up1 + u0;
            let vpp0 = vp0 + upp0 * (y1 + vp1 + upp0.square());
            super::branch("CPADD30");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let k2 = vp1 + c1;
        let t3 = vp0 + c0;
        let k1 = t3 + up1 * k2;
        let k0 = f2 + vp1 * y1 + k2 * (up0 + vp1) + up1 * k1;
        let sp0 = k0 + u0 * (k1 + u0 * k2);
        if sp0.is_zero() {
            if k2.is_zero() {
                super::branch("CPADD31");
                return DivisorCoords::deg0(y1, y0, 2);
            }
            let w2 = k2.inv();
            let upp0 = w2 * (vp0 + c0) + up1 + u0;
            let vpp0 = upp0 * k2 + vp0 + h0;
            super::branch("CPADD32");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let t0 = s0 * up1 + h1;
        let upp1 = w4 * k2 + s0 + u0;
        let upp0 = w4 * (c0 + vp0 + k2 * up1) + t0 + u0 * upp1;
        let vpp1 = upp1 * s0 + t0 + vp1;
        let vpp0 = s0 * (upp0 + up0) + vh0;
        super::branch("CPADD33");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let k2 = c1 + vp1;
    let sp0 = v0 + vp0 + u0 * (k2 + h1 + t1);
    if sp0.is_zero() {
        if k2.is_zero() {
            super::branch("CPADD34");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = k2.inv();
        let upp0 = w2 * (vp0 + c0) + up1 + u0;
        let vpp0 = upp0 * k2 + vp0 + h0;
        super::branch("CPADD35");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let t1 = s0 * up1 + h1;
    let upp1 = s0 + w4 * k2 + u0;
    let upp0 = t1 + w4 * (c0 + vp0 + k2 * up1) + u0 * upp1;
    let vpp1 = upp1 * s0 + t1 + vp1;
    let vpp0 = s0 * (upp0 + up0) + h0 + vp0;
    super::branch("CPADD36");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-2 divisors. (`Deg2ADD`, pos)
#[inline]
fn deg2_add_pos<F: Field>(
    u0: F,
    u1: F,
    v0: F,
    v1: F,
    up0: F,
    up1: F,
    vp0: F,
    vp1: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (f2, h0, h1, y0, y1) = (cc.f2, cc.h0, cc.h1, cc.y0, cc.y1);
    let c1 = cc.yn1;

    let m3 = up1 + u1;
    let m4 = u0 + up0;
    let m1 = m4 + up1 * m3;
    let m2 = up0 * m3;
    let d = m1 * m4 + m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            // u = up
            let t1 = v1 + h1;
            let dw21 = vp1 + t1 + u0 + u1.square();
            let dw20 = vp0 + v0 + h0 + u0 * u1;
            if dw20.is_zero() && dw21.is_zero() {
                super::branch("CPADD37");
                return DivisorCoords::deg0(y1, y0, 1);
            }
            let k2 = y1 + v1;
            let t3 = y0 + v0;
            let k1 = t3 + u1 * k2;
            let k0 = f2 + v1 * t1 + u0 * k2 + u1 * k1;
            let b2 = dw21.inv();
            let u0 = u1 + dw20 * b2;
            let s0 = b2 * (k0 + u0 * (k1 + u0 * k2));
            let upp1 = u0 + u0;
            let upp0 = u0.square();
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0 * s0;
            super::branch("CPADD38");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t0 = v0 + h0;
        let vh1 = v1 + h1;
        let m2cap = m3.square();
        let m3cap = m3 * m2cap;
        let dw3 = m3cap * (vp0 + t0) + m4 * (m2cap * (vp1 + vh1) + m4.square());
        if dw3.is_zero() {
            let a1 = m3.inv();
            let s1cap = m4 * a1;
            let u0 = u1 + s1cap;
            let up0 = up1 + s1cap;
            let s0 = a1 * (vp0 + v0 + up0 * (vp1 + v1));
            let upp1 = u0 + up0;
            let upp0 = u0 * up0;
            let vpp1 = v1 + s0;
            let vpp0 = v0 + s0 * u0;
            super::branch("CPADD39");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let k2 = y1 + v1;
        let t3 = y0 + v0;
        let k1 = t3 + u1 * k2;
        let k0 = f2 + v1 * vh1 + u0 * k2 + u1 * k1;
        let a12 = m2cap * (vp1 + vh1 + up0 + up1.square());
        let sp1 = a12 * (vp1 + v1) + m3cap * (k1 + up1 * k2);
        let sp0 = a12 * (vp0 + v0) + m3cap * (k0 + up0 * k2);
        let dd = dw3.square();
        if sp1.is_zero() {
            if sp0.is_zero() {
                super::branch("CPADD40");
                return DivisorCoords::deg0(y1, y0, 0);
            }
            let w2 = (dw3 * sp0).inv();
            let s0 = sp0 * sp0 * w2;
            let w3 = dd * w2;
            let upp0 = w3 * (s0.square() + y1 + v1) + up1;
            let vpp0 = v0 + h0 + s0 * u0 + upp0 * (s0 * u1 + c1 + v1 + upp0 * (s0 + upp0));
            super::branch("CPADD41");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let w1 = sp1 * (sp1 + dw3);
        if w1.is_zero() {
            let w0 = sp0 + dw3 * u1;
            if w0.is_zero() {
                super::branch("CPADD42");
                return DivisorCoords::deg0(y1, y0, 2);
            }
            let w2 = (dw3 * w0).inv();
            let s0 = sp0 * w0 * w2;
            let w3 = dd * w2;
            let t0 = s0 * u1;
            let t2 = t0 + u0 + v1 + c1;
            let t3 = s0 + u1;
            let upp0 = w3 * (s0 * t3 + t2) + up1;
            let vpp0 = v0 + h0 + s0 * u0 + upp0 * (t2 + upp0 * t3);
            super::branch("CPADD43");
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
        let l1 = (s0 + s1) * (u0 + u1) + l0 + t1;
        let t4 = l1 + h1;
        let upp1 = w4 * (s1 * t1 + s0) + up1;
        let upp0 = w4 * (s0 * l2 + s1 * t4 + y1 + v1) + up0 + up1 * upp1;
        let t0 = F::one() + s1;
        let t1 = upp1 * t0;
        let t2 = l2 + t1;
        let t3 = upp0 * t2;
        let vpp1 = v1 + t4 + (upp0 + upp1) * (t0 + t2) + t1 + t3;
        let vpp0 = v0 + h0 + l0 + t3;
        super::branch("CPADD44");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let r0 = vp0 + v0;
    let r1 = vp1 + v1;
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("CPADD45");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = (d * sp0).inv();
        let s0 = sp0 * sp0 * w2;
        let w3 = dd * w2;
        let upp0 = w3 * (s0.square() + y1 + v1) + up1;
        let vpp0 = v0 + h0 + s0 * u0 + upp0 * (s0 * u1 + c1 + v1 + upp0 * (s0 + upp0));
        super::branch("CPADD46");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w1 = sp1 * (sp1 + d);
    if w1.is_zero() {
        let w0 = sp0 + d * u1;
        if w0.is_zero() {
            super::branch("CPADD47");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t0 = s0 * u1;
        let t2 = t0 + u0 + v1 + c1;
        let t3 = s0 + u1;
        let upp0 = w3 * (s0 * t3 + t2) + up1;
        let vpp0 = v0 + h0 + s0 * u0 + upp0 * (t2 + upp0 * t3);
        super::branch("CPADD48");
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
    let l1 = (s0 + s1) * (u0 + u1) + l0 + t1;
    let upp1 = w4 * (s1 * t1 + s0) + up1;
    let upp0 = w4 * (s0 * l2 + s1 * (l1 + h1) + y1 + v1) + up0 + up1 * upp1;
    let t0 = F::one() + s1;
    let t1 = upp1 * t0;
    let t2 = l2 + t1;
    let t3 = upp0 * t2;
    let vpp1 = v1 + h1 + l1 + (upp0 + upp1) * (t0 + t2) + t1 + t3;
    let vpp0 = v0 + h0 + l0 + t3;
    super::branch("CPADD49");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Add two divisors in the **positive** reduced basis. (`ADD` dispatcher, pos)
#[inline]
pub fn add_pos<F: Field>(
    d1: &DivisorCoords<F>,
    d2: &DivisorCoords<F>,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let n = d1.n + d2.n - 1;
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
                super::branch("CPADD50");
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
        (1, 1) => match n {
            1 => deg1_add_dwn_pos(d1.u0, d1.v0, d2.u0, d2.v0, cc),
            -1 => deg1_add_up_pos(d1.u0, d1.v0, d2.u0, d2.v0, cc),
            _ => deg1_add_pos(d1.u0, d1.v0, d2.u0, d2.v0, cc),
        },
        (1, 0) => match d2.n {
            0 => {
                if d1.n == 0 {
                    deg01_add_up_pos(d1.u0, d1.v0, cc)
                } else {
                    super::branch("CPADD51");
                    DivisorCoords { n: 0, ..*d1 }
                }
            }
            1 => {
                super::branch("CPADD52");
                *d1
            }
            _ => {
                if d1.n == 0 {
                    super::branch("CPADD53");
                    DivisorCoords { n: 1, ..*d1 }
                } else {
                    deg01_add_dwn_pos(d1.u0, d1.v0, cc)
                }
            }
        },
        (0, 2) => match d1.n {
            0 => deg02_add_up_pos(d2.u0, d2.u1, d2.v0, d2.v1, cc),
            1 => {
                super::branch("CPADD54");
                DivisorCoords { n: 0, ..*d2 }
            }
            _ => deg02_add_dwn_pos(d2.u0, d2.u1, d2.v0, d2.v1, cc),
        },
        (0, 1) => match d1.n {
            0 => {
                if d2.n == 0 {
                    deg01_add_up_pos(d2.u0, d2.v0, cc)
                } else {
                    super::branch("CPADD55");
                    DivisorCoords { n: 0, ..*d2 }
                }
            }
            1 => {
                super::branch("CPADD56");
                *d2
            }
            _ => {
                if d2.n == 0 {
                    super::branch("CPADD57");
                    DivisorCoords { n: 1, ..*d2 }
                } else {
                    deg01_add_dwn_pos(d2.u0, d2.v0, cc)
                }
            }
        },
        _ => {
            super::branch("CPADD58");
            DivisorCoords::deg0(cc.y1, cc.y0, n)
        }
    }
}
