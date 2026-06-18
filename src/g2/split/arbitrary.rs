//! Arbitrary-characteristic genus 2 **split** model divisor arithmetic.
//!
//! Explicit formulas for split (two points at infinity) genus 2 hyperelliptic
//! curves over **any** field, including characteristic 2 and 3:
//! `y² + h(x)·y = f(x)` with `deg f = 6`, `deg h ≤ 3`.
//!
//! Divisor representation and balance weight `n` are as in [`super::not_char2`];
//! the high coefficients of `v` come from the reduced-basis polynomial
//! `Vpl = y3·x³ + y2·x² + y1·x + y0` (positive) or `Vn = −Vpl − h` (negative).
//!
//! Ported from the Magma `arb_splitG2_{UTL,DBL,ADD}.mag` (neg + pos reduced).

use crate::field::Field;
use crate::poly::Poly;

pub use super::not_char2::{DivisorCoords, G};

/// Curve constants and precomputations for an arbitrary-characteristic split
/// genus 2 curve. Built by [`precompute`].
#[derive(Clone, Copy, Debug)]
pub struct CurveConstants<F: Field> {
    // f = f6·x⁶ + … + f0, h = h3·x³ + … + h0
    pub f0: F, pub f1: F, pub f2: F, pub f3: F, pub f4: F, pub f5: F, pub f6: F,
    pub h0: F, pub h1: F, pub h2: F, pub h3: F,
    // Vpl = y3 x³ + y2 x² + y1 x + y0; Vn = −Vpl − h = yn3 x³ + …
    pub y0: F, pub y1: F, pub y2: F, pub y3: F,
    pub yn0: F, pub yn1: F, pub yn2: F, pub yn3: F,
    // c-constants: shared c2,c3,c4; basis-specific c0,c1 and pos-only c5.
    pub c0n: F, pub c1n: F, // neg: c_i = 2·y_i + h_i
    pub c0p: F, pub c1p: F, pub c5p: F, // pos: c0,c1 = y_i + h_i, c5 = 2y1 + h1
    pub c2: F, pub c3: F, pub c4: F,
    // d-precomputations d0..d9 (shared between bases).
    pub d0: F, pub d1: F, pub d2: F, pub d3: F, pub d4: F,
    pub d5: F, pub d6: F, pub d7: F, pub d8: F, pub d9: F,
    // degree-0 adjust divisor: shared adu (au), per-basis adv.
    pub au0: F, pub au1: F, pub au2: F, pub audeg: i8,
    pub adv0n: F, pub adv1n: F,
    pub adv0p: F, pub adv1p: F,
}

impl<F: Field> CurveConstants<F> {
    pub fn f_poly(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.f0, self.f1, self.f2, self.f3, self.f4, self.f5, self.f6])
    }
    pub fn h_poly(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.h0, self.h1, self.h2, self.h3])
    }
    /// Positive reduced-basis polynomial `Vpl`.
    pub fn vpl(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.y0, self.y1, self.y2, self.y3])
    }
    /// Negative reduced-basis polynomial `Vn = −Vpl − h`.
    pub fn vn(&self) -> Poly<F> {
        Poly::from_coeffs(vec![self.yn0, self.yn1, self.yn2, self.yn3])
    }
}

/// Build [`CurveConstants`] from the curve coefficients and a chosen root
/// `y3` of `x² + h3·x − f6` (the leading coefficient of `Vpl`; the two roots
/// are the curve's points at infinity). Port of the Magma `Precompute`.
#[allow(clippy::too_many_arguments)]
pub fn precompute<F: Field>(f: [F; 7], h: [F; 4], y3: F) -> CurveConstants<F> {
    let [f0, f1, f2, f3, f4, f5, f6] = f;
    let [h0, h1, h2, h3] = h;

    let yn3 = -y3 - h3;
    let c3 = y3 - yn3; // 2·y3 + h3
    let d5 = c3.inv();
    let d4 = f4 - h1 * y3;
    let d3 = f3 - h0 * y3;

    let y2 = (f5 - y3 * h2) * d5;
    let yn2 = -y2 - h2;
    let c2 = y2 - yn2; // 2·y2 + h2

    let y1 = (d4 + y2 * yn2) * d5;
    let yn1 = -y1 - h1;
    let c1n = y1 - yn1; // 2·y1 + h1

    let y0 = (d3 - y2 * c1n - y1 * h2) * d5;
    let yn0 = -y0 - h0;
    let c0n = y0 - yn0; // 2·y0 + h0
    let c4 = c2 + c3;

    // positive-basis c-constants
    let c0p = y0 + h0;
    let c1p = y1 + h1;
    let c5p = c1p + y1; // 2·y1 + h1 = c1n

    let d2 = f2 + yn2 * y0 + yn0 * y2;
    let d1 = d2 + y1 * yn1;
    let d0 = f1 - h0 * y1;
    let d6 = d5 * c2;
    let d7 = d5 * d0;
    let d8 = d5 * d1;
    let d9 = d5 * c1n;

    // degree-0 adjust divisor (adu shared; adv per basis).
    let k2 = d1;
    let k1 = f1 + yn1 * y0 + yn0 * y1;
    let k0 = f0 + y0 * yn0;

    let (au0, au1, au2, audeg, adv0n, adv1n, adv0p, adv1p);
    if k2.is_zero() {
        if k1.is_zero() {
            au2 = F::zero(); au1 = F::zero(); au0 = F::one(); audeg = 0;
            adv1n = yn1; adv0n = yn0;
            adv1p = y1; adv0p = y0;
        } else {
            let w1 = k1.inv();
            au2 = F::zero(); au1 = F::one(); au0 = k0 * w1; audeg = 1;
            adv1n = yn1;
            adv0n = y0 - au0 * (c1n - au0 * (c2 - au0 * c3));
            adv1p = y1;
            adv0p = au0 * (y1 + c1p - au0 * (c2 - au0 * c3)) - c0p;
        }
    } else {
        let w1 = k2.inv();
        au2 = F::one(); au1 = k1 * w1; au0 = k0 * w1; audeg = 2;
        let w2 = c2 - au1 * c3;
        adv1n = y1 - au0 * c3 - au1 * w2;
        adv0n = y0 - au0 * w2;
        adv1p = au0 * c3 + au1 * w2 - c1p;
        adv0p = au0 * w2 - c0p;
    }

    CurveConstants {
        f0, f1, f2, f3, f4, f5, f6, h0, h1, h2, h3,
        y0, y1, y2, y3, yn0, yn1, yn2, yn3,
        c0n, c1n, c0p, c1p, c5p, c2, c3, c4,
        d0, d1, d2, d3, d4, d5, d6, d7, d8, d9,
        au0, au1, au2, audeg, adv0n, adv1n, adv0p, adv1p,
    }
}

// ===========================================================================
// Doubling — negative reduced basis. Port of `arb_splitG2_DBL.mag` (neg).
// Labels `ADBLnn`.
// ===========================================================================

/// Double a degree-1 divisor `<x + u0, v0, n=1>` (DWN). (`Deg1DBLDWN`)
#[inline]
fn deg1_dbl_dwn_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, yn0, yn1) = (cc.h0, cc.yn0, cc.yn1);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);
    let (d0, d5, d6, d8, d9) = (cc.d0, cc.d5, cc.d6, cc.d8, cc.d9);

    let t0 = u0 * c3;
    let t1 = -u0 * (c1 - u0 * (c2 - t0));
    let vp0 = v0 - t1;
    let vh0 = -vp0 - h0;
    let d = v0 - vh0;
    if d.is_zero() {
        super::branch("ADBL00");
        return DivisorCoords::deg0(yn1, yn0, 2);
    }
    let z0 = vh0 - yn0;
    let z1 = d8 + z0 * (d6 - u0);
    let sp0 = d0 - c1 * vp0 - t0 * (z1.double() - u0 * z0);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("ADBL01");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        super::branch("ADBL02");
        return DivisorCoords::deg1(upp0, yn1, vh0, 0);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = d6 - z0 * w2 - u0;
    let upp0 = d9 + s0 * d5 - z1 * w2 - u0 * upp1;
    let vpp1 = yn1 - s0;
    let vpp0 = vh0 - s0 * u0;
    super::branch("ADBL03");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-1 divisor `<x + u0, v0, n=0>` (UP). (`Deg1DBLUP`)
#[inline]
fn deg1_dbl_up_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y1, yn0, yn1) = (cc.h0, cc.y1, cc.yn0, cc.yn1);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);
    let (d0, d5, d6, d8, d9) = (cc.d0, cc.d5, cc.d6, cc.d8, cc.d9);

    let t0 = u0 * c3;
    let vh0 = -v0 - h0;
    let d = v0 - vh0 + u0 * (c1 - u0 * (c2 - t0));
    if d.is_zero() {
        super::branch("ADBL04");
        return DivisorCoords::deg0(yn1, yn0, 0);
    }
    let z0 = v0 - yn0;
    let z1 = d8 + z0 * (d6 - u0);
    let sp0 = d0 - c1 * vh0 - t0 * (z1.double() - u0 * z0);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("ADBL05");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        let vpp0 = -upp0 * (c1 - upp0 * (c2 - upp0 * c3)) + vh0;
        super::branch("ADBL06");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = d6 + z0 * w2 - u0;
    let upp0 = d9 - s0 * d5 + z1 * w2 - u0 * upp1;
    let t0 = upp1 * c3;
    let t1 = c2 - t0;
    let t2 = upp0 * t1;
    let vpp1 = y1 - s0 - (upp0 + upp1) * (c3 + t1) + t0 + t2;
    let vpp0 = vh0 - s0 * u0 - t2;
    super::branch("ADBL07");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-2 divisor `<x² + u1x + u0, v1x + v0, n=0>`. (`Deg2DBL`)
#[inline]
fn deg2_dbl_neg<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y1) = (cc.h0, cc.h1, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (yn0, yn1) = (cc.yn0, cc.yn1);
    let (c2, c3, c4) = (cc.c2, cc.c3, cc.c4);

    let t0 = v0 + h0;
    let t1 = v1 + h1;
    let t2 = u1.square();
    let t3 = c3 * u1 - c2;
    let m3 = c3 * (t2 - u0) - c2 * u1 - t1 - v1;
    let m4 = t0 + v0 - u0 * t3;
    let m1 = m4 + m3 * u1;
    let m2 = -m3 * u0;
    let d = m4 * m1 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            super::branch("ADBL08");
            return DivisorCoords::deg0(yn1, yn0, 1);
        }
        let b1 = -m3.inv();
        let t2 = v1 - yn1;
        let t3 = v0 - yn0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - u1 * k2;
        let k0 = d2 + t4 - v1 * t1 - u0 * k2 - u1 * k1;
        let u0 = u1 - m4 * b1;
        let s0 = b1 * (k0 - u0 * (k1 - u0 * k2));
        let upp1 = u0.double();
        let upp0 = u0.square();
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0 * s0;
        super::branch("ADBL09");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let z0 = v0 - yn0;
    let z1 = v1 - yn1;
    let r1 = z0 - z1 * (u1.double() - d6);
    let r0 = d6 * z0 + d5 * (d2 - v1 * t1) - u1 * r1 - z1 * (t2 + u0.double());
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("ADBL10");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let upp0 = -s0 + d6 - u1 + w4 * d5 * z1;
        let vpp0 = c3 * (upp0 * (s0 * u1 + d5 * (v1 - y1) - upp0 * (s0 - d6 + upp0)) - s0 * u0) - t0;
        super::branch("ADBL11");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }

    let w1 = sp1 * (sp1 - d);
    if w1.is_zero() {
        let t1b = d6 - u1;
        let w0 = sp0 - sp1 * t1b;
        if w0.is_zero() {
            super::branch("ADBL12");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t2 = s0 - t1b;
        let upp0 = w3 * (t2 * s0 - d5 * (m3 + z1));
        let t1c = v1 - y1 + c3 * (s0 * u1 + u0 - upp0 * t2);
        let vpp0 = upp0 * t1c - t0 - s0 * u0 * c3;
        super::branch("ADBL13");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }

    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * d * dd;
    let s0 = w3 * sp0;
    let s1 = w3 * sp1;
    let t4 = s0 - d6 + u1;
    let upp1 = w4 * ((s0 + t4) * s1 - s0);
    let upp0 = w4 * (t4 * s0 - d5 * (m3 * s1 + z1));
    let zz0 = upp0 - u0;
    let zz1 = upp1 - u1;
    let ww0 = zz0 * s0;
    let ww1 = zz1 * s1;
    let ww2 = upp1 - d6 - ww1;
    let vpp1 = c3 * ((s0 + s1) * (zz0 + zz1) - ww0 - ww1 - upp0 + ww2 * upp1) - t1;
    let vpp0 = c3 * (ww0 + ww2 * upp0) - t0;
    super::branch("ADBL14");
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
                super::branch("ADBL15");
                DivisorCoords::deg0(cc.yn1, cc.yn0, 1)
            } else if d.n == 0 {
                super::branch("ADBL16");
                DivisorCoords { u2: cc.au2, u1: cc.au1, u0: cc.au0, v1: cc.adv1n, v0: cc.adv0n, n: 2 - cc.audeg as i32 }
            } else {
                super::branch("ADBL17");
                DivisorCoords { u2: cc.au2, u1: cc.au1, u0: cc.au0, v1: cc.yn1, v0: cc.yn0, n: 0 }
            }
        }
    }
}

// ===========================================================================
// Addition — negative reduced basis. Port of `arb_splitG2_ADD.mag` (neg).
// Labels `AADDnn`. In every sub-function `(u0[,u1],v0[,v1])` is the first
// operand and `(up0[,up1],vp0[,vp1])` the second; `v1`/`vp1` are the low
// x-coefficients (the high ones come from the `Vn` basis).
// ===========================================================================

/// `<1, V, 0>` + degree-1 divisor `<x + u0, v0, 1>`. (`Deg01ADDUP`)
#[inline]
fn deg01_add_up_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y1, yn0, yn1) = (cc.h0, cc.y1, cc.yn0, cc.yn1);
    let (d6, d7, d8, d9) = (cc.d6, cc.d7, cc.d8, cc.d9);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);

    let z3 = v0 - yn0;
    if z3.is_zero() {
        let z2 = d8 + z3 * d6;
        if z2.is_zero() {
            super::branch("AADD00");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z2.inv();
        let vh0 = v0 + h0;
        let up0 = w2 * (d7 + vh0 * d9) - u0;
        let vp0 = up0 * (up0 * (c2 - c3 * up0) - c1) - vh0;
        super::branch("AADD01");
        return DivisorCoords::deg1(up0, yn1, vp0, 1);
    }
    let w2 = z3.inv();
    let vh0 = v0 + h0;
    let up1 = w2 * d8 + d6 - u0;
    let up0 = w2 * (d7 + vh0 * d9) - u0 * up1;
    let t1 = up1 * c3;
    let t2 = c2 - t1;
    let t3 = up0 * t2;
    let vp1 = y1 - (up1 + up0) * (t2 + c3) + t1 + t3;
    let vp0 = -vh0 - t3;
    super::branch("AADD02");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 2>` + degree-1 divisor `<x + u0, v0, 0>`. (`Deg01ADDDWN`)
#[inline]
fn deg01_add_dwn_neg<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, yn0, yn1) = (cc.h0, cc.yn0, cc.yn1);
    let (d6, d7, d8, d9) = (cc.d6, cc.d7, cc.d8, cc.d9);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);

    let v0 = v0 + u0 * (c1 - u0 * (c2 - c3 * u0));
    let vp0 = -v0 - h0;
    let z3 = vp0 - yn0;
    if z3.is_zero() {
        let z2 = d8 + z3 * d6;
        if z2.is_zero() {
            super::branch("AADD03");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z2.inv();
        let up0 = w2 * (d7 - v0 * d9) - u0;
        super::branch("AADD04");
        return DivisorCoords::deg1(up0, yn1, vp0, 0);
    }
    let w2 = z3.inv();
    let up1 = w2 * d8 + d6 - u0;
    let up0 = w2 * (d7 - v0 * d9) - u0 * up1;
    super::branch("AADD05");
    DivisorCoords::deg2(up1, up0, yn1, vp0, 0)
}

/// `<1, V, 0>` + degree-2 divisor. (`Deg02ADDUP`)
#[inline]
fn deg02_add_up_neg<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1) = (cc.h0, cc.h1, cc.yn0, cc.yn1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c1, c2, c3, c4) = (cc.c1n, cc.c2, cc.c3, cc.c4);

    let z1 = v1 - yn1;
    let z0 = v0 - yn0;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("AADD06");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z0.inv();
        let up0 = w2 * d5 * (d2 - v1 * (h1 + v1)) + d6 - u1;
        let vp0 = up0 * (up0 * (c2 - c3 * up0) - c1) - v0 - h0;
        super::branch("AADD07");
        return DivisorCoords::deg1(up0, yn1, vp0, 1);
    }
    let w2 = z1.inv();
    let t1 = h1 + v1;
    let t2 = w2 * z0;
    let up1 = t2 + d6 - u1;
    let up0 = w2 * d5 * (d2 - v1 * t1) + d6 * t2 - u0 - u1 * up1;
    let t2 = c3 * up1;
    let t3 = (t2 - c2) * up0;
    let vp1 = t2 - (c4 - t2) * (up0 + up1) - t3 - t1;
    let vp0 = t3 - h0 - v0;
    super::branch("AADD08");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 2>` + degree-2 divisor. (`Deg02ADDDWN`)
#[inline]
fn deg02_add_dwn_neg<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1) = (cc.h0, cc.h1, cc.yn0, cc.yn1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c2, c3) = (cc.c2, cc.c3);

    let t0 = u1 * c3;
    let t1 = t0 - c2;
    let t2 = u0 * t1;
    let v1 = v1 - (u0 + u1) * (-c3 + t1) + t2 - t0;
    let v0 = v0 - t2;
    let vp1 = -v1 - h1;
    let vp0 = -v0 - h0;
    let z1 = vp1 - yn1;
    let z0 = vp0 - yn0;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("AADD09");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z0.inv();
        let up0 = w2 * d5 * (d2 + v1 * vp1) + d6 - u1;
        super::branch("AADD10");
        return DivisorCoords::deg1(up0, yn1, vp0, 0);
    }
    let w2 = z1.inv();
    let t2 = w2 * z0;
    let up1 = t2 + d6 - u1;
    let up0 = w2 * d5 * (d2 + v1 * vp1) + d6 * t2 - u0 - up1 * u1;
    super::branch("AADD11");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// Two degree-1 divisors with `n + np = 1`. (`Deg1ADD`)
#[inline]
fn deg1_add_neg<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, yn0, yn1) = (cc.h0, cc.yn0, cc.yn1);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);
    let (d0, d1) = (cc.d0, cc.d1);

    let d = u0 - up0;
    if d.is_zero() {
        let upp0 = u0.square();
        let t3 = c3 * u0;
        let vh0 = v0 + h0;
        let dw = vh0 + vp0 + u0 * c1 - upp0 * (c2 - t3);
        if dw.is_zero() {
            super::branch("AADD12");
            return DivisorCoords::deg0(yn1, yn0, 1);
        }
        let upp1 = u0.double();
        let t0 = v0 - yn0;
        let t1 = d1 + c2 * t0;
        let t2 = d0 + c1 * vh0 - u0 * (t1.double() - t3 * (t0 + t0 + t0));
        let s0 = t2 * dw.inv();
        let vpp1 = yn1 + s0;
        let vpp0 = v0 + s0 * u0;
        super::branch("AADD13");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let s0 = (vp0 - v0) * d.inv();
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;
    let vpp1 = yn1 + s0;
    let vpp0 = v0 + u0 * s0;
    super::branch("AADD14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 0` (UP). (`Deg1ADDUP`)
#[inline]
fn deg1_add_up_neg<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y1, yn0, yn1) = (cc.h0, cc.y1, cc.yn0, cc.yn1);
    let (d5, d6, d8, d9) = (cc.d5, cc.d6, cc.d8, cc.d9);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);

    let d = u0 - up0;
    if d.is_zero() {
        super::branch("AADD15");
        return DivisorCoords::deg0(yn1, yn0, 0);
    }
    let sp0 = vp0 - v0;
    let z0 = v0 - yn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("AADD16");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w1 = z0.inv();
        let upp0 = w1 * (d8 + z0 * (d6 - u0)) - up0;
        let vpp0 = upp0 * (upp0 * (c2 - upp0 * c3) - c1) - v0 - h0;
        super::branch("AADD17");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t1 = z0 * w1;
    let upp1 = d6 + t1 - up0;
    let upp0 = d9 + w1 * d8 + t1 * (d6 - u0) - s0 * d5 - up0 * upp1;
    let t0 = upp1 * c3;
    let t1 = c2 - t0;
    let t2 = upp0 * t1;
    let vpp1 = y1 - s0 - (upp0 + upp1) * (c3 + t1) + t0 + t2;
    let vpp0 = -v0 - h0 - s0 * u0 - t2;
    super::branch("AADD18");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 1` (DWN). (`Deg1ADDDWN`)
#[inline]
fn deg1_add_dwn_neg<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, yn0, yn1) = (cc.h0, cc.yn0, cc.yn1);
    let (d5, d6, d8, d9) = (cc.d5, cc.d6, cc.d8, cc.d9);
    let (c1, c2, c3) = (cc.c1n, cc.c2, cc.c3);

    let d = u0 - up0;
    if d.is_zero() {
        super::branch("AADD19");
        return DivisorCoords::deg0(yn1, yn0, 2);
    }
    let vt0 = v0 + u0 * (c1 - u0 * (c2 - u0 * c3));
    let sp0 = vp0 - vt0 + up0 * (c1 - up0 * (c2 - up0 * c3));
    let vh0 = -vt0 - h0;
    let z0 = vh0 - yn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("AADD20");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w1 = z0.inv();
        let upp0 = w1 * d8 + d6 - u0 - up0;
        super::branch("AADD21");
        return DivisorCoords::deg1(upp0, yn1, vh0, 0);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = z0 * w1;
    let upp1 = d6 - t0 - up0;
    let upp0 = d9 + s0 * d5 - w1 * d8 - t0 * (d6 - u0) - up0 * upp1;
    let vpp1 = yn1 - s0;
    let vpp0 = vh0 - s0 * u0;
    super::branch("AADD22");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,1>` + degree-2 (UP). (`Deg12ADDUP`)
#[inline]
fn deg12_add_up_neg<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1, y1) = (cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c2, c3, c4) = (cc.c2, cc.c3, cc.c4);

    let t0 = u0 * up1;
    let d = up0 - t0 + u0.square();
    if d.is_zero() {
        let dw = v0 + vp0 + h0 - u0 * (vp1 - y1 + u0 * (c2 - c3 * u0));
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let vpp0 = vp0 + upp0 * (yn1 - vp1);
            super::branch("AADD23");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let t2 = vp1 - yn1;
        let t3 = vp0 - yn0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - up1 * k2;
        let k0 = d2 + t4 - vp1 * (vp1 + h1) - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("AADD24");
                return DivisorCoords::deg0(yn1, yn0, 2);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 + d6 - up1 - u0;
            let vpp0 = upp0 * (vp1 - y1 + upp0 * (c2 - upp0 * c3)) - vp0 - h0;
            super::branch("AADD25");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1 + h1;
        let upp1 = d6 + w4 * t2 - d5 * s0 - u0;
        let upp0 = w4 * (vp0 - yn0 + t2 * (d6 - up1)) - d5 * (vh1 + vp1) - u0 * upp1;
        let t0 = upp1 * c3;
        let t1 = c2 - s0 - t0;
        let t2 = upp0 * t1;
        let vpp1 = t0 + t2 - (upp1 + upp0) * (c3 + t1) - vh1;
        let vpp0 = -t2 - s0 * up0 - vp0 - h0;
        super::branch("AADD26");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = vp1 - yn1;
    let sp0 = v0 - vp0 + z0 * u0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("AADD27");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (vp0 - yn0) + d6 - up1 - u0;
        let vpp0 = upp0 * (vp1 - y1 + upp0 * (c2 - upp0 * c3)) - vp0 - h0;
        super::branch("AADD28");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1 + h1;
    let upp1 = d6 + w4 * z0 - d5 * s0 - u0;
    let upp0 = w4 * (vp0 - yn0 + z0 * (d6 - up1)) - d5 * (vh1 + vp1) - u0 * upp1;
    let t0 = upp1 * c3;
    let t1 = c2 - s0 - t0;
    let t2 = upp0 * t1;
    let vpp1 = t0 + t2 - (upp1 + upp0) * (c3 + t1) - vh1;
    let vpp0 = -t2 - s0 * up0 - vp0 - h0;
    super::branch("AADD29");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,0>` + degree-2 (DWN). (`Deg12ADD`)
#[inline]
fn deg12_add_neg<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, yn0, yn1, y1) = (cc.h0, cc.h1, cc.yn0, cc.yn1, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c2, c3, c4) = (cc.c2, cc.c3, cc.c4);

    // vp := -V - h - ((-V-h - vp) mod up)
    let t2 = up1 * c3;
    let t3 = c2 - t2;
    let t4 = up0 * t3;
    let vp1 = vp1 + (up0 + up1) * (c3 + t3) - t2 - t4;
    let vp0 = vp0 + t4;

    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 - t0 + t1;
    if d.is_zero() {
        let dw = vp0 + v0 + h0 + u0 * (y1 - vp1);
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let vpp0 = vp0 + upp0 * (yn1 - vp1 + upp0 * (c2 - c3 * upp0));
            super::branch("AADD30");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let vh1 = -vp1 - h1;
        let vh0 = -vp0 - h0;
        let t2 = vh1 - yn1;
        let t3 = vh0 - yn0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - up1 * k2;
        let k0 = d2 + t4 + vp1 * vh1 - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("AADD31");
                return DivisorCoords::deg0(yn1, yn0, 0);
            }
            let w2 = t2.inv();
            let upp0 = d6 + w2 * t3 - up1 - u0;
            let vpp0 = -upp0 * t2 + vh0;
            super::branch("AADD32");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1 + h1;
        let upp1 = d5 * s0 + d6 - w4 * t2 - u0;
        let upp0 = d5 * (vp1 + vh1) - w4 * (vh0 - yn0 + t2 * (d6 - up1)) - u0 * upp1;
        let vpp1 = upp1 * s0 - vh1;
        let vpp0 = s0 * (upp0 - up0) + vh0;
        super::branch("AADD33");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let vh1 = -vp1 - h1;
    let vh0 = -vp0 - h0;
    let z0 = vh1 - yn1;
    let sp0 = v0 - vp0 - u0 * (yn1 - vp1 + u0 * c2 - t1 * c3);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("AADD34");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (vh0 - yn0) + d6 - up1 - u0;
        let vpp0 = vh0 - upp0 * z0;
        super::branch("AADD35");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 0);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 - vh1;
    let upp1 = d5 * s0 + d6 - w4 * z0 - u0;
    let upp0 = d5 * (vp1 + vh1) - w4 * (vh0 - yn0 + z0 * (d6 - up1)) - u0 * upp1;
    let vpp1 = upp1 * s0 - vh1;
    let vpp0 = s0 * (upp0 - up0) - h0 - vp0;
    super::branch("AADD36");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-2 divisors. (`Deg2ADD`)
#[inline]
fn deg2_add_neg<F: Field>(
    u0: F, u1: F, v0: F, v1: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y1, yn0, yn1) = (cc.h0, cc.h1, cc.y1, cc.yn0, cc.yn1);
    let (c2, c3, c4) = (cc.c2, cc.c3, cc.c4);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);

    let m3 = up1 - u1;
    let m4 = u0 - up0;
    let m1 = m4 + up1 * m3;
    let m2 = -up0 * m3;
    let d = m1 * m4 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            // u = up
            let t1 = v1 + h1;
            let t2 = u1 * c3;
            let t3 = c2 - t2;
            let t4 = u0 * t3;
            let dw21 = vp1 + t1 + (u0 + u1) * (c3 + t3) - t2 - t4;
            let dw20 = vp0 + v0 + h0 + t4;
            if dw20.is_zero() && dw21.is_zero() {
                super::branch("AADD37");
                return DivisorCoords::deg0(yn1, yn0, 1);
            }
            let t2 = v1 - yn1;
            let t3 = v0 - yn0;
            let t4 = c2 * t3;
            let k2 = c3 * t2;
            let k1 = c4 * (t3 + t2) - k2 - t4 - u1 * k2;
            let k0 = d2 + t4 - v1 * t1 - u0 * k2 - u1 * k1;
            let b2 = dw21.inv();
            let u0 = u1 - dw20 * b2;
            let s0 = b2 * (k0 - u0 * (k1 - u0 * k2));
            let upp1 = u0.double();
            let upp0 = u0.square();
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0 * s0;
            super::branch("AADD38");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t0 = v0 + h0;
        let vh1 = v1 + h1;
        let m2sq = m3.square();
        let m3cu = -m3 * m2sq;
        let dw3 = m3cu * (vp0 + t0) - m4 * (m2sq * (vp1 + vh1) - m4 * (m3 * c2 + m4 * c3));
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
            super::branch("AADD39");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t2 = v1 - yn1;
        let t3 = v0 - yn0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - u1 * k2;
        let k0 = d2 + t4 - v1 * vh1 - u0 * k2 - u1 * k1;
        let t2 = c2 - up1 * c3;
        let a12 = -m2sq * (vp1 + vh1 + up0 * c3 + up1 * t2);
        let sp1 = a12 * (vp1 - v1) + m3cu * (k1 - up1 * k2);
        let sp0 = a12 * (vp0 - v0) + m3cu * (k0 - up0 * k2);
        let dd = dw3.square();
        if sp1.is_zero() {
            if sp0.is_zero() {
                super::branch("AADD40");
                return DivisorCoords::deg0(yn1, yn0, 2);
            }
            let w3 = (dw3 * sp0).inv();
            let s0 = sp0.square() * w3;
            let w4 = dd * w3;
            let t0 = s0 * u1;
            let upp0 = d6 - s0 * d5 + w4 * (v1 - yn1) - up1;
            let vpp0 = upp0 * (t0 - y1 + v1 + upp0 * (c2 - s0 - c3 * upp0)) - v0 - h0 - s0 * u0;
            super::branch("AADD41");
            return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
        }
        let w1 = -sp1 * (sp1 - c3 * dw3);
        if w1.is_zero() {
            let t1 = c2 - c3 * u1;
            let w0 = sp0 - dw3 * t1;
            if w0.is_zero() {
                super::branch("AADD42");
                return DivisorCoords::deg0(yn1, yn0, 0);
            }
            let w2 = (dw3 * w0).inv();
            let s0 = sp0 * w0 * w2;
            let w3 = dd * w2;
            let t2 = s0 * u1 + c3 * u0 + v1 - y1;
            let upp0 = w3 * t2 + s0 * d5 - up1;
            let vpp0 = upp0 * (t2 - upp0 * (s0 - t1)) - v0 - h0 - s0 * u0;
            super::branch("AADD43");
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
        let vh1 = l1 + v1 + h1;
        let t2 = c2 - l2;
        let upp1 = w4 * (s1 * (t2 - s0) + s0 * c3) - up1;
        let upp0 = w4 * (s0 * t2 - s1 * (v1 + vh1) + c3 * (v1 - yn1)) - up0 - up1 * upp1;
        let t0 = s1 - c3;
        let t1 = upp1 * t0;
        let t2 = -t2 - t1;
        let t3 = upp0 * t2;
        let vpp1 = (upp0 + upp1) * (t0 + t2) - vh1 - t1 - t3;
        let vpp0 = t3 - v0 - h0 - l0;
        super::branch("AADD44");
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
            super::branch("AADD45");
            return DivisorCoords::deg0(yn1, yn0, 2);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let t0 = s0 * u1;
        let upp0 = d6 - s0 * d5 + w4 * (v1 - yn1) - up1;
        let vpp0 = upp0 * (t0 - y1 + v1 + upp0 * (c2 - s0 - c3 * upp0)) - v0 - h0 - s0 * u0;
        super::branch("AADD46");
        return DivisorCoords::deg1(upp0, yn1, vpp0, 1);
    }
    let w1 = -sp1 * (sp1 - c3 * d);
    if w1.is_zero() {
        let t1 = c2 - c3 * u1;
        let w0 = sp0 - d * t1;
        if w0.is_zero() {
            super::branch("AADD47");
            return DivisorCoords::deg0(yn1, yn0, 0);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t2 = s0 * u1 + c3 * u0 + v1 - y1;
        let upp0 = w3 * t2 + s0 * d5 - up1;
        let vpp0 = upp0 * (t2 - upp0 * (s0 - t1)) - v0 - h0 - s0 * u0;
        super::branch("AADD48");
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
    let l1 = (s0 + s1) * (u0 + u1) - l0 - t1;
    let vh1 = l1 + v1 + h1;
    let t2 = c2 - l2;
    let upp1 = w4 * (s1 * (t2 - s0) + s0 * c3) - up1;
    let upp0 = w4 * (s0 * t2 - s1 * (v1 + vh1) + c3 * (v1 - yn1)) - up0 - up1 * upp1;
    let t0 = s1 - c3;
    let t1 = upp1 * t0;
    let t2 = -t2 - t1;
    let t3 = upp0 * t2;
    let vpp1 = (upp0 + upp1) * (t0 + t2) - vh1 - t1 - t3;
    let vpp0 = t3 - v0 - h0 - l0;
    super::branch("AADD49");
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
                super::branch("AADD50");
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
                    super::branch("AADD51");
                    DivisorCoords { n: 0, ..*d1 }
                }
            }
            1 => {
                super::branch("AADD52");
                *d1
            }
            _ => {
                if d1.n == 0 {
                    super::branch("AADD53");
                    DivisorCoords { n: 1, ..*d1 }
                } else {
                    deg01_add_dwn_neg(d1.u0, d1.v0, cc)
                }
            }
        },
        (0, 2) => match d1.n {
            0 => deg02_add_up_neg(d2.u0, d2.u1, d2.v0, d2.v1, cc),
            1 => {
                super::branch("AADD54");
                DivisorCoords { n: 0, ..*d2 }
            }
            _ => deg02_add_dwn_neg(d2.u0, d2.u1, d2.v0, d2.v1, cc),
        },
        (0, 1) => match d1.n {
            0 => {
                if d2.n == 0 {
                    deg01_add_up_neg(d2.u0, d2.v0, cc)
                } else {
                    super::branch("AADD55");
                    DivisorCoords { n: 0, ..*d2 }
                }
            }
            1 => {
                super::branch("AADD56");
                *d2
            }
            _ => {
                if d2.n == 0 {
                    super::branch("AADD57");
                    DivisorCoords { n: 1, ..*d2 }
                } else {
                    deg01_add_dwn_neg(d2.u0, d2.v0, cc)
                }
            }
        },
        _ => {
            super::branch("AADD58");
            DivisorCoords::deg0(cc.yn1, cc.yn0, n)
        }
    }
}

// ===========================================================================
// Doubling — positive reduced basis. Port of `arb_splitG2_DBL.mag` (pos).
// Labels `APDBLnn`. Uses Vpl coeffs (y0,y1) and pos c-constants (c0p,c1p,c5p).
// ===========================================================================

/// Double a degree-1 divisor `<x + u0, v0, n=1>` (DWN). (`Deg1DBLDWN`, pos)
#[inline]
fn deg1_dbl_dwn_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    let (c1, c2, c3, c5) = (cc.c1p, cc.c2, cc.c3, cc.c5p);
    let (d0, d5, d6, d8, d9) = (cc.d0, cc.d5, cc.d6, cc.d8, cc.d9);

    let t0 = u0 * c3;
    let vh0 = v0 + h0;
    let d = v0 + vh0 - u0 * (c5 - u0 * (c2 - t0));
    if d.is_zero() {
        super::branch("APDBL00");
        return DivisorCoords::deg0(y1, y0, 2);
    }
    let z0 = y0 - v0;
    let z1 = d8 + z0 * (d6 - u0);
    let sp0 = d0 - c5 * v0 - t0 * (z1.double() - u0 * z0);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("APDBL01");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        let vpp0 = upp0 * (c5 - upp0 * (c2 - upp0 * c3)) - vh0;
        super::branch("APDBL02");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = d6 - z0 * w2 - u0;
    let upp0 = d9 + s0 * d5 - z1 * w2 - u0 * upp1;
    let t0 = upp1 * c3;
    let t1 = c2 - t0;
    let t2 = upp0 * t1;
    let vpp1 = (upp0 + upp1) * (c3 + t1) - t0 - t2 - c1 - s0;
    let vpp0 = t2 - vh0 - s0 * u0;
    super::branch("APDBL03");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-1 divisor `<x + u0, v0, n=0>` (UP). (`Deg1DBLUP`, pos)
#[inline]
fn deg1_dbl_up_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    let (c2, c3, c5) = (cc.c2, cc.c3, cc.c5p);
    let (d0, d5, d6, d8, d9) = (cc.d0, cc.d5, cc.d6, cc.d8, cc.d9);

    let t0 = u0 * c3;
    let t1 = u0 * (c5 - u0 * (c2 - t0));
    let v0 = v0 - t1;
    let vp0 = -v0 - h0;
    let d = v0 - vp0 + t1;
    if d.is_zero() {
        super::branch("APDBL04");
        return DivisorCoords::deg0(y1, y0, 0);
    }
    let z0 = y0 - vp0;
    let z1 = d8 + z0 * (d6 - u0);
    let sp0 = d0 - c5 * vp0 - t0 * (z1.double() - u0 * z0);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("APDBL05");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w1 = z0.inv();
        let upp0 = z1 * w1 - u0;
        super::branch("APDBL06");
        return DivisorCoords::deg1(upp0, y1, vp0, 1);
    }
    let w1 = (d * sp0).inv();
    let w2 = d.square() * w1;
    let s0 = w1 * sp0.square();
    let upp1 = d6 + z0 * w2 - u0;
    let upp0 = d9 - s0 * d5 + z1 * w2 - u0 * upp1;
    let vpp1 = y1 - s0;
    let vpp0 = vp0 - s0 * u0;
    super::branch("APDBL07");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Double a degree-2 divisor `<x² + u1x + u0, v1x + v0, n=0>`. (`Deg2DBL`, pos)
#[inline]
fn deg2_dbl_pos<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c1, c2, c3, c4) = (cc.c1p, cc.c2, cc.c3, cc.c4);

    let t0 = v0 + h0;
    let t1 = v1 + h1;
    let t2 = u1.square();
    let t3 = c2 - c3 * u1;
    let m3 = c2 * u1 - c3 * (t2 - u0) - t1 - v1;
    let m4 = t0 + v0 - u0 * t3;
    let m1 = m4 + m3 * u1;
    let m2 = -m3 * u0;
    let d = m4 * m1 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            super::branch("APDBL08");
            return DivisorCoords::deg0(y1, y0, 1);
        }
        let b1 = -m3.inv();
        let t2 = y1 - v1;
        let t3 = y0 - v0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - u1 * k2;
        let k0 = d2 + t4 - v1 * t1 - u0 * k2 - u1 * k1;
        let u0 = u1 - m4 * b1;
        let s0 = b1 * (k0 - u0 * (k1 - u0 * k2));
        let upp1 = u0.double();
        let upp0 = u0.square();
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0 * s0;
        super::branch("APDBL09");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let z0 = y0 - v0;
    let z1 = y1 - v1;
    let r1 = z0 - z1 * (u1.double() - d6);
    let r0 = d6 * z0 + d5 * (d2 - v1 * t1) - u1 * r1 - z1 * (t2 + u0.double());
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("APDBL10");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let upp0 = s0 + d6 - u1 - w4 * d5 * z1;
        let vpp0 = c3 * (upp0 * (s0 * u1 + d5 * (c1 + v1) - upp0 * (d6 + s0 - upp0)) - s0 * u0) - t0;
        super::branch("APDBL11");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }

    let w1 = sp1 * (sp1 + d);
    if w1.is_zero() {
        let t1 = d6 - u1;
        let w0 = sp1 * t1 - sp0;
        if w0.is_zero() {
            super::branch("APDBL12");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t2 = s0 + t1;
        let upp0 = w3 * (t2 * s0 + d5 * (m3 - z1));
        let t1 = c1 + v1 + c3 * (s0 * u1 - u0 - upp0 * t2);
        let vpp0 = upp0 * t1 - t0 - s0 * u0 * c3;
        super::branch("APDBL13");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }

    let w2 = (d * w1).inv();
    let w3 = w2 * w1;
    let w4 = w2 * d * dd;
    let s0 = w3 * sp0;
    let s1 = w3 * sp1;
    let t4 = s0 + d6 - u1;
    let upp1 = w4 * ((s0 + t4) * s1 + s0);
    let upp0 = w4 * (t4 * s0 - d5 * (m3 * s1 + z1));
    let z0 = upp0 - u0;
    let z1 = upp1 - u1;
    let w0 = z0 * s0;
    let w1 = z1 * s1;
    let w2 = d6 - w1 - upp1;
    let vpp1 = c3 * ((s0 + s1) * (z0 + z1) - w0 - w1 + upp0 + w2 * upp1) - t1;
    let vpp0 = c3 * (w0 + w2 * upp0) - t0;
    super::branch("APDBL14");
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
                super::branch("APDBL15");
                DivisorCoords::deg0(cc.y1, cc.y0, 1)
            } else if d.n == 0 {
                super::branch("APDBL16");
                // <adu, Vpl, 2 − audeg>
                DivisorCoords { u2: cc.au2, u1: cc.au1, u0: cc.au0, v1: cc.y1, v0: cc.y0, n: 2 - cc.audeg as i32 }
            } else {
                super::branch("APDBL17");
                // <adu, adv_pos, 0>
                DivisorCoords { u2: cc.au2, u1: cc.au1, u0: cc.au0, v1: cc.adv1p, v0: cc.adv0p, n: 0 }
            }
        }
    }
}

// ===========================================================================
// Addition — positive reduced basis. Port of `arb_splitG2_ADD.mag` (pos).
// Labels `APADDnn`. In every sub-function `(u0[,u1],v0[,v1])` is the first
// operand and `(up0[,up1],vp0[,vp1])` the second; `v1`/`vp1` are the low
// x-coefficients (the high ones come from the `Vpl` basis).
// ===========================================================================

/// `<1, V, 2>` + degree-1 divisor `<x + u0, v0, 0>`. (`Deg01ADDDWN`, pos)
#[inline]
fn deg01_add_dwn_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    let (d6, d7, d8, d9) = (cc.d6, cc.d7, cc.d8, cc.d9);
    let (c1, c2, c3, c5) = (cc.c1p, cc.c2, cc.c3, cc.c5p);

    let z0 = y0 - v0;
    if z0.is_zero() {
        let z1 = d8 + z0 * d6;
        if z1.is_zero() {
            super::branch("APADD00");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = z1.inv();
        let up0 = w2 * (d7 - v0 * d9) - u0;
        let vp0 = up0 * (c5 + up0 * (c3 * up0 - c2)) - h0 - v0;
        super::branch("APADD01");
        return DivisorCoords::deg1(up0, y1, vp0, 0);
    }
    let w2 = z0.inv();
    let up1 = w2 * d8 + d6 - u0;
    let up0 = w2 * (d7 - v0 * d9) - u0 * up1;
    let t1 = up1 * c3;
    let t2 = c2 - t1;
    let t3 = up0 * t2;
    let vp1 = (up1 + up0) * (c3 + t2) - t1 - t3 - c1;
    let vp0 = t3 - v0 - h0;
    super::branch("APADD02");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 0>` + degree-1 divisor `<x + u0, v0, 0>`. (`Deg01ADDUP`, pos)
#[inline]
fn deg01_add_up_pos<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    let (d6, d7, d8, d9) = (cc.d6, cc.d7, cc.d8, cc.d9);
    let (c0, c2, c3, c5) = (cc.c0p, cc.c2, cc.c3, cc.c5p);

    let v0 = v0 - u0 * (c5 + u0 * (c3 * u0 - c2));
    let vp0 = -v0 - h0;
    let z0 = c0 + v0;
    if z0.is_zero() {
        let z1 = d8 + z0 * d6;
        if z1.is_zero() {
            super::branch("APADD03");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = z1.inv();
        let up0 = w2 * (d7 - vp0 * d9) - u0;
        super::branch("APADD04");
        return DivisorCoords::deg1(up0, y1, vp0, 1);
    }
    let w2 = z0.inv();
    let up1 = w2 * d8 + d6 - u0;
    let up0 = w2 * (d7 - vp0 * d9) - u0 * up1;
    super::branch("APADD05");
    DivisorCoords::deg2(up1, up0, y1, vp0, 0)
}

/// `<1, V, 2>` + degree-2 divisor. (`Deg02ADDDWN`, pos)
#[inline]
fn deg02_add_dwn_pos<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c2, c3, c4, c5) = (cc.c2, cc.c3, cc.c4, cc.c5p);

    let z0 = y0 - v0;
    let z1 = y1 - v1;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("APADD06");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = z0.inv();
        let up0 = w2 * d5 * (d2 - v1 * (h1 + v1)) + d6 - u1;
        let vp0 = up0 * (c5 - up0 * (c2 - c3 * up0)) - v0 - h0;
        super::branch("APADD07");
        return DivisorCoords::deg1(up0, y1, vp0, 0);
    }
    let w2 = z1.inv();
    let t1 = h1 + v1;
    let t2 = w2 * z0;
    let up1 = t2 + d6 - u1;
    let up0 = w2 * d5 * (d2 - v1 * t1) + d6 * t2 - u0 - u1 * up1;
    let t2 = c3 * up1;
    let t3 = (c2 - t2) * up0;
    let vp1 = (c4 - t2) * (up0 + up1) - t1 - t2 - t3;
    let vp0 = t3 - h0 - v0;
    super::branch("APADD08");
    DivisorCoords::deg2(up1, up0, vp1, vp0, 0)
}

/// `<1, V, 0>` + degree-2 divisor. (`Deg02ADDUP`, pos)
#[inline]
fn deg02_add_up_pos<F: Field>(u0: F, u1: F, v0: F, v1: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c0, c1, c2, c3) = (cc.c0p, cc.c1p, cc.c2, cc.c3);

    let t0 = u1 * c3;
    let t1 = c2 - t0;
    let t2 = u0 * t1;
    let v1 = v1 - (u0 + u1) * (c3 + t1) + t2 + t0;
    let v0 = v0 - t2;
    let z0 = v0 + c0;
    let z1 = v1 + c1;
    if z1.is_zero() {
        if z0.is_zero() {
            super::branch("APADD09");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = z0.inv();
        let up0 = w2 * d5 * (d2 - v1 * (v1 + h1)) + d6 - u1;
        let vp0 = -v0 - h0;
        super::branch("APADD10");
        return DivisorCoords::deg1(up0, y1, vp0, 1);
    }
    let w2 = z1.inv();
    let t1 = -v1 - h1;
    let t2 = w2 * z0;
    let up1 = t2 + d6 - u1;
    let up0 = w2 * d5 * (d2 + v1 * t1) + d6 * t2 - u0 - up1 * u1;
    let vp0 = -v0 - h0;
    super::branch("APADD11");
    DivisorCoords::deg2(up1, up0, t1, vp0, 0)
}

/// Two degree-1 divisors with `n + np = 1`. (`Deg1ADD`, pos)
#[inline]
fn deg1_add_pos<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    // Magma `Deg1ADD` declares c1 but doesn't use it; kept for line-by-line parity.
    let (_c1, c2, c3, c5) = (cc.c1p, cc.c2, cc.c3, cc.c5p);
    let (d0, d1) = (cc.d0, cc.d1);

    let d = u0 - up0;
    if d.is_zero() {
        let upp0 = u0.square();
        let t3 = c3 * u0;
        let dw = v0 + vp0 + h0 - u0 * c5 + upp0 * (c2 - t3);
        if dw.is_zero() {
            super::branch("APADD12");
            return DivisorCoords::deg0(y1, y0, 1);
        }
        let upp1 = u0.double();
        let t0 = y0 - v0;
        let t1 = d1 + c2 * t0;
        let t2 = d0 - c5 * v0 - u0 * (t1.double() - t3 * (t0 + t0 + t0));
        let s0 = t2 * dw.inv();
        let vpp1 = y1 + s0;
        let vpp0 = v0 + s0 * u0;
        super::branch("APADD13");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let s0 = (vp0 - v0) * d.inv();
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;
    let vpp1 = y1 + s0;
    let vpp0 = v0 + u0 * s0;
    super::branch("APADD14");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 1` (DWN). (`Deg1ADDDWN`, pos)
#[inline]
fn deg1_add_dwn_pos<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    let (d5, d6, d8, d9) = (cc.d5, cc.d6, cc.d8, cc.d9);
    let (c1, c2, c3, c5) = (cc.c1p, cc.c2, cc.c3, cc.c5p);

    let d = u0 - up0;
    if d.is_zero() {
        super::branch("APADD15");
        return DivisorCoords::deg0(y1, y0, 2);
    }
    let sp0 = vp0 - v0;
    let z0 = y0 - v0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("APADD16");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w1 = z0.inv();
        let upp0 = w1 * d8 + d6 - u0 - up0;
        let vpp0 = -v0 - h0 + upp0 * (c5 - upp0 * (c2 - upp0 * c3));
        super::branch("APADD17");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t1 = z0 * w1;
    let upp1 = d6 - t1 - up0;
    let upp0 = d9 - w1 * d8 - t1 * (d6 - u0) + s0 * d5 - up0 * upp1;
    let t0 = upp1 * c3;
    let t1 = c2 - t0;
    let t2 = upp0 * t1;
    let vpp1 = -c1 - s0 + (upp0 + upp1) * (c3 + t1) - t0 - t2;
    let vpp0 = -v0 - h0 - s0 * u0 + t2;
    super::branch("APADD18");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-1 divisors, both `n = 0` (UP). (`Deg1ADDUP`, pos)
#[inline]
fn deg1_add_up_pos<F: Field>(u0: F, v0: F, up0: F, vp0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let (h0, y0, y1) = (cc.h0, cc.y0, cc.y1);
    let (d5, d6, d8, d9) = (cc.d5, cc.d6, cc.d8, cc.d9);
    let (c2, c3, c5) = (cc.c2, cc.c3, cc.c5p);

    let d = u0 - up0;
    if d.is_zero() {
        super::branch("APADD19");
        return DivisorCoords::deg0(y1, y0, 0);
    }
    let vt0 = v0 + u0 * (u0 * (c2 - u0 * c3) - c5);
    let sp0 = vp0 - vt0 + up0 * (up0 * (c2 - up0 * c3) - c5);
    let vn0 = -vt0 - h0;
    let z0 = y0 - vn0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("APADD20");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w1 = z0.inv();
        let upp0 = w1 * (d8 + z0 * (d6 - u0)) - up0;
        super::branch("APADD21");
        return DivisorCoords::deg1(upp0, y1, vn0, 1);
    }
    let w0 = (sp0 * d).inv();
    let w1 = w0 * d.square();
    let s0 = w0 * sp0.square();
    let t0 = z0 * w1;
    let upp1 = d6 + t0 - up0;
    let upp0 = d9 - s0 * d5 + w1 * d8 + t0 * (d6 - u0) - up0 * upp1;
    let vpp1 = y1 - s0;
    let vpp0 = vn0 - s0 * u0;
    super::branch("APADD22");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,1>` + degree-2 (DWN). (`Deg12ADD`, pos)
#[inline]
fn deg12_add_pos<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c1, c2, c3, c4) = (cc.c1p, cc.c2, cc.c3, cc.c4);

    let t0 = u0 * up1;
    let d = up0 - t0 + u0.square();
    if d.is_zero() {
        let dw = v0 + vp0 + h0 - u0 * (vp1 + c1 - u0 * (c2 - c3 * u0));
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let vpp0 = vp0 + upp0 * (y1 - vp1);
            super::branch("APADD23");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let t2 = y1 - vp1;
        let t3 = y0 - vp0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - up1 * k2;
        let k0 = d2 + t4 - vp1 * (vp1 + h1) - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("APADD24");
                return DivisorCoords::deg0(y1, y0, 0);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 + d6 - up1 - u0;
            let vpp0 = upp0 * (c1 + vp1 - upp0 * (c2 - upp0 * c3)) - vp0 - h0;
            super::branch("APADD25");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1 + h1;
        let upp1 = d5 * s0 + d6 - w4 * t2 - u0;
        let upp0 = d5 * (vh1 + vp1) - w4 * (y0 - vp0 + t2 * (d6 - up1)) - u0 * upp1;
        let t0 = upp1 * c3;
        let t1 = c2 + s0 - t0;
        let t2 = upp0 * t1;
        let vpp1 = (upp1 + upp0) * (c3 + t1) - t0 - t2 - vh1;
        let vpp0 = t2 - s0 * up0 - vp0 - h0;
        super::branch("APADD26");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = y1 - vp1;
    let sp0 = v0 - vp0 - z0 * u0;
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("APADD27");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (y0 - vp0) + d6 - up1 - u0;
        let vpp0 = upp0 * (c1 + vp1 - upp0 * (c2 - upp0 * c3)) - vp0 - h0;
        super::branch("APADD28");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1 + h1;
    let upp1 = d5 * s0 + d6 - w4 * z0 - u0;
    let upp0 = d5 * (vh1 + vp1) - w4 * (y0 - vp0 + z0 * (d6 - up1)) - u0 * upp1;
    let t0 = upp1 * c3;
    let t1 = c2 + s0 - t0;
    let t2 = upp0 * t1;
    let vpp1 = (upp1 + upp0) * (c3 + t1) - t0 - t2 - vh1;
    let vpp0 = t2 - s0 * up0 - vp0 - h0;
    super::branch("APADD29");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// degree-1 `<x+u0,v0,0>` + degree-2 (UP). (`Deg12ADDUP`, pos)
#[inline]
fn deg12_add_up_pos<F: Field>(
    u0: F, v0: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);
    let (c0, c1, c2, c3, c4) = (cc.c0p, cc.c1p, cc.c2, cc.c3, cc.c4);

    let t2 = up1 * c3;
    let t3 = c2 - t2;
    let t4 = up0 * t3;
    let vp1 = vp1 - (up0 + up1) * (c3 + t3) + t2 + t4;
    let vp0 = vp0 - t4;

    let t0 = u0 * up1;
    let t1 = u0.square();
    let d = up0 - t0 + t1;
    if d.is_zero() {
        let dw = vp0 + v0 + h0 - u0 * (c1 + vp1);
        if dw.is_zero() {
            let upp0 = up1 - u0;
            let vpp0 = vp0 + upp0 * (y1 - vp1 - upp0 * (c2 - c3 * upp0));
            super::branch("APADD30");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let t2 = vp1 + c1;
        let t3 = vp0 + c0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - up1 * k2;
        let k0 = d2 + t4 - vp1 * (vp1 + h1) - up0 * k2 - up1 * k1;
        let sp0 = k0 - u0 * (k1 - u0 * k2);
        if sp0.is_zero() {
            if t2.is_zero() {
                super::branch("APADD31");
                return DivisorCoords::deg0(y1, y0, 2);
            }
            let w2 = t2.inv();
            let upp0 = w2 * t3 + d6 - up1 - u0;
            let vpp0 = upp0 * t2 - vp0 - h0;
            super::branch("APADD32");
            return DivisorCoords::deg1(upp0, y1, vpp0, 1);
        }
        let w2 = (sp0 * dw).inv();
        let w3 = w2 * sp0;
        let w4 = dw.square() * w2;
        let s0 = w3 * sp0;
        let vh1 = s0 * up1 + vp1 + h1;
        let upp1 = w4 * t2 - d5 * s0 + d6 - u0;
        let upp0 = w4 * (c0 + vp0 + t2 * (d6 - up1)) - d5 * (vp1 + vh1) - u0 * upp1;
        let vpp1 = upp1 * s0 - vh1;
        let vpp0 = s0 * (upp0 - up0) - h0 - vp0;
        super::branch("APADD33");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }
    let z0 = c1 + vp1;
    let sp0 = v0 - vp0 - u0 * (y1 - vp1 - u0 * c2 + t1 * c3);
    if sp0.is_zero() {
        if z0.is_zero() {
            super::branch("APADD34");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = z0.inv();
        let upp0 = w2 * (vp0 + c0) + d6 - up1 - u0;
        let vpp0 = upp0 * z0 - vp0 - h0;
        super::branch("APADD35");
        return DivisorCoords::deg1(upp0, y1, vpp0, 1);
    }
    let w2 = (sp0 * d).inv();
    let w3 = w2 * sp0;
    let w4 = d.square() * w2;
    let s0 = w3 * sp0;
    let vh1 = s0 * up1 + vp1 + h1;
    let upp1 = w4 * z0 - d5 * s0 + d6 - u0;
    let upp0 = w4 * (c0 + vp0 + z0 * (d6 - up1)) - d5 * (vp1 + vh1) - u0 * upp1;
    let vpp1 = upp1 * s0 - vh1;
    let vpp0 = s0 * (upp0 - up0) - h0 - vp0;
    super::branch("APADD36");
    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0)
}

/// Two degree-2 divisors. (`Deg2ADD`, pos)
#[inline]
fn deg2_add_pos<F: Field>(
    u0: F, u1: F, v0: F, v1: F, up0: F, up1: F, vp0: F, vp1: F, cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let (h0, h1, y0, y1) = (cc.h0, cc.h1, cc.y0, cc.y1);
    let (c1, c2, c3, c4) = (cc.c1p, cc.c2, cc.c3, cc.c4);
    let (d2, d5, d6) = (cc.d2, cc.d5, cc.d6);

    let m3 = up1 - u1;
    let m4 = u0 - up0;
    let m1 = m4 + up1 * m3;
    let m2 = -up0 * m3;
    let d = m1 * m4 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            // u = up
            let t1 = v1 + h1;
            let t2 = u1 * c3;
            let t3 = c2 - t2;
            let t4 = u0 * t3;
            let dw21 = vp1 + t1 - (u0 + u1) * (c3 + t3) + t2 + t4;
            let dw20 = vp0 + v0 + h0 - t4;
            if dw20.is_zero() && dw21.is_zero() {
                super::branch("APADD37");
                return DivisorCoords::deg0(y1, y0, 1);
            }
            let t2 = y1 - v1;
            let t3 = y0 - v0;
            let t4 = c2 * t3;
            let k2 = c3 * t2;
            let k1 = c4 * (t3 + t2) - k2 - t4 - u1 * k2;
            let k0 = d2 + t4 - v1 * t1 - u0 * k2 - u1 * k1;
            let b2 = dw21.inv();
            let u0 = u1 - dw20 * b2;
            let s0 = b2 * (k0 - u0 * (k1 - u0 * k2));
            let upp1 = u0.double();
            let upp0 = u0.square();
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0 * s0;
            super::branch("APADD38");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t0 = v0 + h0;
        let vh1 = v1 + h1;
        let m2sq = m3.square();
        let m3cu = -m3 * m2sq;
        let dw3 = m3cu * (vp0 + t0) - m4 * (m2sq * (vp1 + vh1) + m4 * (m3 * c2 + m4 * c3));
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
            super::branch("APADD39");
            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
        }
        let t2 = y1 - v1;
        let t3 = y0 - v0;
        let t4 = c2 * t3;
        let k2 = c3 * t2;
        let k1 = c4 * (t3 + t2) - k2 - t4 - u1 * k2;
        let k0 = d2 + t4 - v1 * vh1 - u0 * k2 - u1 * k1;
        let t2 = c2 - up1 * c3;
        let a12 = -m2sq * (vp1 + vh1 - up0 * c3 - up1 * t2);
        let sp1 = a12 * (vp1 - v1) + m3cu * (k1 - up1 * k2);
        let sp0 = a12 * (vp0 - v0) + m3cu * (k0 - up0 * k2);
        let dd = dw3.square();
        if sp1.is_zero() {
            if sp0.is_zero() {
                super::branch("APADD40");
                return DivisorCoords::deg0(y1, y0, 0);
            }
            let w3 = (dw3 * sp0).inv();
            let s0 = sp0.square() * w3;
            let w4 = dd * w3;
            let t0 = s0 * u1;
            let upp0 = w4 * (s0 * (s0 * d5 + d6) - y1 + v1) - up1;
            let vpp0 = upp0 * (t0 + c1 + v1 - upp0 * (c2 + s0 - c3 * upp0)) - v0 - h0 - s0 * u0;
            super::branch("APADD41");
            return DivisorCoords::deg1(upp0, y1, vpp0, 0);
        }
        let w1 = -sp1 * (sp1 + c3 * dw3);
        if w1.is_zero() {
            let t1 = c2 - c3 * u1;
            let w0 = sp0 + dw3 * t1;
            if w0.is_zero() {
                super::branch("APADD42");
                return DivisorCoords::deg0(y1, y0, 2);
            }
            let w2 = (dw3 * w0).inv();
            let s0 = sp0 * w0 * w2;
            let w3 = dd * w2;
            let t0 = s0 * u1;
            let t2 = t0 - c3 * u0 + v1 + c1;
            let t3 = s0 + t1;
            let upp0 = w3 * (t2 - s0 * d5 * t3) - up1;
            let vpp0 = upp0 * (t2 - upp0 * t3) - v0 - h0 - s0 * u0;
            super::branch("APADD43");
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
        let t2 = l2 + c2;
        let l1 = (s0 + s1) * (u0 + u1) - l0 - t1;
        let upp1 = -w4 * (s1 * (t2 + s0) + s0 * c3) - up1;
        let upp0 = -w4 * (s0 * t2 + s1 * (l1 + v1 + vh1) - c3 * (y1 - v1)) - up0 - up1 * upp1;
        let t0 = c3 + s1;
        let t1 = upp1 * t0;
        let t2 = t2 - t1;
        let t3 = upp0 * t2;
        let vpp1 = (upp0 + upp1) * (t0 + t2) - vh1 - l1 - t1 - t3;
        let vpp0 = t3 - v0 - h0 - l0;
        super::branch("APADD44");
        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0, 0);
    }

    let r0 = vp0 - v0;
    let r1 = vp1 - v1;
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;
    let dd = d.square();

    if sp1.is_zero() {
        if sp0.is_zero() {
            super::branch("APADD45");
            return DivisorCoords::deg0(y1, y0, 0);
        }
        let w3 = (d * sp0).inv();
        let s0 = sp0.square() * w3;
        let w4 = dd * w3;
        let t0 = s0 * u1;
        let upp0 = w4 * (s0 * (s0 * d5 + d6) - y1 + v1) - up1;
        let vpp0 = upp0 * (t0 + c1 + v1 - upp0 * (c2 + s0 - c3 * upp0)) - v0 - h0 - s0 * u0;
        super::branch("APADD46");
        return DivisorCoords::deg1(upp0, y1, vpp0, 0);
    }
    let w1 = -sp1 * (sp1 + c3 * d);
    if w1.is_zero() {
        let t1 = c2 - c3 * u1;
        let w0 = sp0 + d * t1;
        if w0.is_zero() {
            super::branch("APADD47");
            return DivisorCoords::deg0(y1, y0, 2);
        }
        let w2 = (d * w0).inv();
        let s0 = sp0 * w0 * w2;
        let w3 = dd * w2;
        let t0 = s0 * u1;
        let t2 = t0 - c3 * u0 + v1 + c1;
        let t3 = s0 + t1;
        let upp0 = w3 * (t2 - s0 * d5 * t3) - up1;
        let vpp0 = upp0 * (t2 - upp0 * t3) - v0 - h0 - s0 * u0;
        super::branch("APADD48");
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
    let vh1 = v1 + h1;
    let t2 = l2 + c2;
    let upp1 = -w4 * (s1 * (t2 + s0) + s0 * c3) - up1;
    let upp0 = -w4 * (s0 * t2 + s1 * (l1 + v1 + vh1) - c3 * (y1 - v1)) - up0 - up1 * upp1;
    let t0 = c3 + s1;
    let t1 = upp1 * t0;
    let t2 = t2 - t1;
    let t3 = upp0 * t2;
    let vpp1 = (upp0 + upp1) * (t0 + t2) - vh1 - l1 - t1 - t3;
    let vpp0 = t3 - v0 - h0 - l0;
    super::branch("APADD49");
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
                super::branch("APADD50");
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
                    super::branch("APADD51");
                    DivisorCoords { n: 0, ..*d1 }
                }
            }
            1 => {
                super::branch("APADD52");
                *d1
            }
            _ => {
                if d1.n == 0 {
                    super::branch("APADD53");
                    DivisorCoords { n: 1, ..*d1 }
                } else {
                    deg01_add_dwn_pos(d1.u0, d1.v0, cc)
                }
            }
        },
        (0, 2) => match d1.n {
            0 => deg02_add_up_pos(d2.u0, d2.u1, d2.v0, d2.v1, cc),
            1 => {
                super::branch("APADD54");
                DivisorCoords { n: 0, ..*d2 }
            }
            _ => deg02_add_dwn_pos(d2.u0, d2.u1, d2.v0, d2.v1, cc),
        },
        (0, 1) => match d1.n {
            0 => {
                if d2.n == 0 {
                    deg01_add_up_pos(d2.u0, d2.v0, cc)
                } else {
                    super::branch("APADD55");
                    DivisorCoords { n: 0, ..*d2 }
                }
            }
            1 => {
                super::branch("APADD56");
                *d2
            }
            _ => {
                if d2.n == 0 {
                    super::branch("APADD57");
                    DivisorCoords { n: 1, ..*d2 }
                } else {
                    deg01_add_dwn_pos(d2.u0, d2.v0, cc)
                }
            }
        },
        _ => {
            super::branch("APADD58");
            DivisorCoords::deg0(cc.y1, cc.y0, n)
        }
    }
}
