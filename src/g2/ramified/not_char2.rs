//! Not characteristic 2 genus 2 ramified model divisor arithmetic.
//!
//! Explicit formulas for divisor addition and doubling on ramified genus 2
//! hyperelliptic curves over fields NOT of characteristic 2.
//!
//! ## Curve Form
//! `y² = f(x)` where:
//! - `f(x) = x⁵ + f₃x³ + f₂x² + f₁x + f₀` (monic, degree 5, f₄ = 0)
//! - `h(x) = 0` (no h polynomial since char ≠ 2)
//!
//! These are simplified formulas compared to the arbitrary characteristic case.
//!
//! Based on: Sebastian Lindner, 2019

use crate::field::Field;

/// Curve constants for a not-char-2 ramified genus 2 curve.
///
/// Represents `y² = f(x)` where:
/// - `f(x) = x⁵ + f3*x³ + f2*x² + f1*x + f0`
/// - `h(x) = 0`
#[derive(Clone, Copy, Debug)]
pub struct CurveConstants<F: Field> {
    pub f3: F,
    pub f2: F,
    pub f1: F,
    pub f0: F,
}

/// Result of divisor operations.
///
/// Represents a divisor `(u(x), v(x))` where:
/// - If `u2 == 1`: degree 2 divisor with `u = x² + u1*x + u0`, `v = v1*x + v0`
/// - If `u2 == 0` and `u1 == 1`: degree 1 divisor with `u = x + u0`, `v = v0`
/// - If `u2 == 0` and `u1 == 0`: identity (u = 1)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DivisorCoords<F: Field> {
    pub u2: F,
    pub u1: F,
    pub u0: F,
    pub v1: F,
    pub v0: F,
}

impl<F: Field> DivisorCoords<F> {
    /// Create the identity divisor
    pub fn identity() -> Self {
        Self {
            u2: F::zero(),
            u1: F::zero(),
            u0: F::one(),
            v1: F::zero(),
            v0: F::zero(),
        }
    }

    /// Create a degree 1 divisor
    pub fn deg1(u0: F, v0: F) -> Self {
        Self {
            u2: F::zero(),
            u1: F::one(),
            u0,
            v1: F::zero(),
            v0,
        }
    }

    /// Create a degree 2 divisor  
    pub fn deg2(u1: F, u0: F, v1: F, v0: F) -> Self {
        Self {
            u2: F::one(),
            u1,
            u0,
            v1,
            v0,
        }
    }

    /// Get the degree of this divisor
    pub fn degree(&self) -> usize {
        if self.u2.is_one() {
            2
        } else if self.u1.is_one() {
            1
        } else {
            0
        }
    }

    /// Check if this is the identity
    pub fn is_identity(&self) -> bool {
        self.degree() == 0
    }
}

/// Add two degree 1 divisors (not char 2).
///
/// Input: `D1 = (x + u0, v0)`, `D2 = (x + up0, vp0)`
/// Output: `D3 = D1 + D2`
#[inline]
pub fn deg1_add<F: Field>(u0: F, v0: F, up0: F, vp0: F) -> DivisorCoords<F> {
    // d := u mod up = u0 - up0
    let d = u0 - up0;

    if d.is_zero() {
        return DivisorCoords::identity();
    }

    // s := (vp - v) / d
    let w1 = d.inv();
    let s0 = w1 * (vp0 - v0);

    // upp = u * up
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;

    // vpp := (v + u*s) mod upp
    let vpp0 = s0 * u0 + v0;

    DivisorCoords::deg2(upp1, upp0, s0, vpp0)
}

/// Add a degree 1 and degree 2 divisor (not char 2).
///
/// Input: `D1 = (x + u0, v0)`, `D2 = (x² + up1*x + up0, vp1*x + vp0)`
/// Output: `D3 = D1 + D2`
#[inline]
pub fn deg12_add<F: Field>(
    u0: F,
    v0: F,
    up1: F,
    up0: F,
    vp1: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let CurveConstants {
        f3,
        f2,
        f1: _,
        f0: _,
    } = *cc;

    // d := up mod u = up0 - u0*(up1 - u0)
    let d = up0 - u0 * (up1 - u0);

    if d.is_zero() {
        // dw := (v + vp) mod u  (h = 0)
        let dw = v0 + vp0 - u0 * vp1;

        if dw.is_zero() {
            // Result is degree 1
            let upp0 = up1 - u0;
            let vpp0 = vp0 - vp1 * upp0;
            return DivisorCoords::deg1(upp0, vpp0);
        }

        // k := ExactQuotient(f - vp², up)
        // k2 = -up1
        let k1 = f3 - up0 + up1.square();
        let k0 = f2 - vp1.square() + up1 * (up0 - k1);

        // s := dw^-1 * k mod u
        let w1 = dw.inv();
        let s0 = w1 * (k0 - u0 * (k1 + u0 * (up1 + u0)));

        // upp := ExactQuotient(-s*(s*up + 2*vp) + k, u)
        let t0 = s0 * up1 + vp1;
        let upp1 = -up1 - u0 - s0.square();
        let upp0 = k1 - s0 * (t0 + vp1) - u0 * upp1;

        // vpp := (-s*up - vp) mod upp
        let vpp1 = s0 * upp1 - t0;
        let vpp0 = s0 * (upp0 - up0) - vp0;

        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
    }

    // General case
    let w1 = d.inv();
    let s0 = w1 * (v0 - vp0 + vp1 * u0);

    // k := ExactQuotient(f - vp², up)
    let t0 = s0 * up1 + vp1;
    let upp1 = -s0.square() - u0 - up1;
    let upp0 = f3 + up1.square() - up0 - s0 * (t0 + vp1) - upp1 * u0;

    // vpp := (-s*up - vp) mod upp
    let vpp1 = upp1 * s0 - t0;
    let vpp0 = s0 * (upp0 - up0) - vp0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Add two degree 2 divisors (not char 2).
///
/// Input: `D1 = (x² + u1*x + u0, v1*x + v0)`, `D2 = (x² + up1*x + up0, vp1*x + vp0)`
/// Output: `D3 = D1 + D2`
#[inline]
pub fn deg2_add<F: Field>(
    u1: F,
    u0: F,
    v1: F,
    v0: F,
    up1: F,
    up0: F,
    vp1: F,
    vp0: F,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let CurveConstants {
        f3,
        f2,
        f1: _,
        f0: _,
    } = *cc;

    // d := Resultant(u, up) computed with 2x2 system
    let m3 = up1 - u1;
    let m4 = u0 - up0;
    let m1 = m4 + up1 * m3;
    let m2 = -up0 * m3;
    let d = m1 * m4 - m2 * m3;

    // Special case: d = 0
    if d.is_zero() {
        if m3.is_zero() {
            // u = up (same polynomial)
            let dw21 = vp1 + v1;
            let dw20 = vp0 + v0;

            if dw20.is_zero() && dw21.is_zero() {
                return DivisorCoords::identity();
            }

            // b2 := dw21^-1
            let b2 = dw21.inv();

            // k := ExactQuotient(f - v², u)
            // k2 = -u1
            let k1 = f3 - u0 + u1.square();
            let k0 = f2 - v1.square() + u1 * (u0 - k1);

            // u := ExactQuotient(u, dw2 * b2)
            let u0_new = u1 - dw20 * b2;

            // s := b2 * k mod u
            let s0 = b2 * (k0 - u0_new * (k1 + u0_new * (u1 + u0_new)));

            // upp := u²
            let upp1 = u0_new + u0_new;
            let upp0 = u0_new.square();

            // vpp := (v + u*s) mod upp
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0_new * s0;

            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
        }

        // m3 != 0: GCD is degree 1
        let t1 = v1 + vp1;
        let m3_sq = m3.square();
        let dw3 = m3_sq * (vp0 + v0) + m4 * m3 * t1;

        if dw3.is_zero() {
            // GCD divides (v + vp)
            let a1 = -m3.inv();
            let s1 = m4 * a1;

            let u0_new = u1 - s1;
            let up0_new = up1 - s1;

            let s0 = a1 * (vp0 - v0 - up0_new * (vp1 - v1));

            let upp1 = u0_new + up0_new;
            let upp0 = u0_new * up0_new;

            let vpp1 = v1 + s0;
            let vpp0 = v0 + s0 * u0_new;

            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
        }

        // General GCD = 1 case with d = 0
        let k2 = -u1;
        let k1 = f3 - u0 + u1.square();
        let k0 = f2 - v1.square() + u1 * (u0 - k1);

        // a12 = -a1*a2 with weight m3*dw3
        let a12 = -m3_sq * t1;

        // s := (a2*a1*(vp - v) + b2*k) mod up; with weight m3*dw3
        let t0 = u1 + up1;
        let t2 = m3 * m3_sq;
        let sp1 = t2 * (k1 - up0 + up1 * t0) - a12 * (vp1 - v1);
        let sp0 = t2 * (k0 + up0 * t0) - a12 * (vp0 - v0);
        let d_new = m3 * dw3;

        return deg2_add_common(u1, u0, v1, v0, up1, up0, m1, m3, sp1, sp0, d_new, cc);
    }

    // General case: d != 0
    let r0 = vp0 - v0;
    let r1 = vp1 - v1;
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;

    deg2_add_common(u1, u0, v1, v0, up1, up0, m1, m3, sp1, sp0, d, cc)
}

/// Common computation for deg2_add after computing sp0, sp1, d
#[inline]
fn deg2_add_common<F: Field>(
    u1: F,
    u0: F,
    v1: F,
    v0: F,
    up1: F,
    _up0: F,
    m1: F,
    m3: F,
    sp1: F,
    sp0: F,
    d: F,
    _cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    if sp1.is_zero() {
        let w1 = d.inv();
        let s0 = sp0 * w1;
        let upp0 = -u1 - up1 - s0.square();

        let t1 = s0 * (u1 - upp0) + v1;
        let vpp0 = upp0 * t1 - s0 * u0 - v0;

        return DivisorCoords::deg1(upp0, vpp0);
    }

    let w1 = (d * sp1).inv();
    let w2 = w1 * d;
    let w3 = w2 * d;
    let w4 = w3.square();
    let s1 = w1 * sp1.square();
    let spp0 = sp0 * w2;

    let t1 = spp0 - m3;
    let t2 = t1 - w4;
    let t3 = w3 * v1;
    let upp1 = spp0 + t2;
    let upp0 = spp0 * (t1 - m3) + m1 + t3 + t3 + w4 * (u1 + up1);

    let t0 = upp0 - u0;
    let t1 = u1 - upp1;
    let vpp1 = s1 * (t1 * t2 + t0) - v1;
    let vpp0 = s1 * (spp0 * t0 + upp0 * t1) - v0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Double a degree 1 divisor (not char 2).
///
/// Input: `D = (x + u0, v0)`
/// Output: `2D`
#[inline]
pub fn deg1_dbl<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let CurveConstants { f3, f2, f1, f0: _ } = *cc;

    // d := 2v mod u
    let d = v0 + v0;

    if d.is_zero() {
        return DivisorCoords::identity();
    }

    // b1 := d^-1
    let w1 = d.inv();

    // upp := u²
    let upp1 = u0 + u0;
    let upp0 = u0.square();

    // k := ExactQuotient(f - v², u)
    // s := b1*k mod u
    let t0 = upp0 + f3;
    let t1 = t0 + upp0;
    let t2 = upp0 * (t0 + t1 + t1) - f2 * upp1 + f1;
    let vpp1 = t2 * w1;
    let vpp0 = vpp1 * u0 + v0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Double a degree 2 divisor (not char 2).
///
/// Input: `D = (x² + u1*x + u0, v1*x + v0)`
/// Output: `2D`
#[inline]
pub fn deg2_dbl<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let CurveConstants {
        f3,
        f2,
        f1: _,
        f0: _,
    } = *cc;

    // Resultant(u, 2v) computed with 2x2 system
    let m3 = -v1 - v1;
    let m4 = v0 + v0;
    let m1 = m4 + m3 * u1;
    let m2 = -m3 * u0;
    let d = m4 * m1 - m2 * m3;

    // Special case: d = 0
    if d.is_zero() {
        if m3.is_zero() {
            return DivisorCoords::identity();
        }

        // b1 := LeadingCoefficient(dw1)^-1 = -m3^-1
        let b1 = -m3.inv();

        // k := ExactQuotient(f - v², u)
        // k2 = -u1
        let k1 = f3 - u0 + u1.square();
        let k0 = f2 - v1.square() + u0 * u1 - u1 * k1;

        // u := ExactQuotient(u, dw1*b1)
        // s := b1*k mod u
        let u0_new = u1 - m4 * b1;
        let s0 = b1 * (k0 - u0_new * (k1 + u0_new * (u1 + u0_new)));

        // upp := u²
        let upp1 = u0_new + u0_new;
        let upp0 = u0_new.square();

        // vpp := (v + u*s) mod upp
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0_new * s0;

        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
    }

    // General case
    // k := ExactQuotient(f - v², u)
    let t0 = u1.square();
    let t1 = f3 + t0;
    let t2 = u0 + u0;
    let t3 = t1 - t2;
    let r1 = t0 + t0 + t3;
    let r0 = u1 * (t2 - t3) + f2 - v1.square();

    // s := k/(2v) mod u; 2x2 system
    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;

    if sp1.is_zero() {
        let w1 = d.inv();
        let s0 = sp0 * w1;

        let upp0 = -s0.square() - u1 - u1;
        let t1 = s0 * (u1 - upp0) + v1;
        let vpp0 = upp0 * t1 - s0 * u0 - v0;

        return DivisorCoords::deg1(upp0, vpp0);
    }

    let w1 = (d * sp1).inv();
    let w2 = w1 * d;
    let w3 = w2 * d;
    let w4 = w3.square();
    let s1 = w1 * sp1.square();
    let spp0 = sp0 * w2;

    let t2 = spp0 - w4;
    let t3 = w3 * v1 + w4 * u1;
    let upp1 = spp0 + t2;
    let upp0 = spp0.square() + t3 + t3;

    let t0 = upp0 - u0;
    let t1 = u1 - upp1;
    let vpp1 = s1 * (t1 * t2 + t0) - v1;
    let vpp0 = s1 * (spp0 * t0 + upp0 * t1) - v0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Add two divisors of arbitrary degree.
#[inline]
pub fn add<F: Field>(
    d1: &DivisorCoords<F>,
    d2: &DivisorCoords<F>,
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    match (d1.degree(), d2.degree()) {
        (0, _) => *d2,
        (_, 0) => *d1,
        (1, 1) => deg1_add(d1.u0, d1.v0, d2.u0, d2.v0),
        (1, 2) => deg12_add(d1.u0, d1.v0, d2.u1, d2.u0, d2.v1, d2.v0, cc),
        (2, 1) => deg12_add(d2.u0, d2.v0, d1.u1, d1.u0, d1.v1, d1.v0, cc),
        (2, 2) => deg2_add(d1.u1, d1.u0, d1.v1, d1.v0, d2.u1, d2.u0, d2.v1, d2.v0, cc),
        _ => unreachable!(),
    }
}

/// Double a divisor of arbitrary degree.
#[inline]
pub fn double<F: Field>(d: &DivisorCoords<F>, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    match d.degree() {
        0 => DivisorCoords::identity(),
        1 => deg1_dbl(d.u0, d.v0, cc),
        2 => deg2_dbl(d.u1, d.u0, d.v1, d.v0, cc),
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;

    type F7 = PrimeField<7>;

    fn make_curve() -> CurveConstants<F7> {
        // f(x) = x^5 + 2x^3 + x^2 + 4x + 1 (h = 0 for not char 2)
        CurveConstants {
            f3: F7::new(2),
            f2: F7::new(1),
            f1: F7::new(4),
            f0: F7::new(1),
        }
    }

    #[test]
    fn test_identity() {
        let id = DivisorCoords::<F7>::identity();
        assert_eq!(id.degree(), 0);
        assert!(id.is_identity());
    }

    #[test]
    fn test_deg1_creation() {
        let d = DivisorCoords::deg1(F7::new(3), F7::new(5));
        assert_eq!(d.degree(), 1);
    }

    #[test]
    fn test_deg2_creation() {
        let d = DivisorCoords::deg2(F7::new(1), F7::new(2), F7::new(3), F7::new(4));
        assert_eq!(d.degree(), 2);
    }

    #[test]
    fn test_add_identity() {
        let cc = make_curve();
        let id = DivisorCoords::identity();
        let d = DivisorCoords::deg1(F7::new(3), F7::new(5));

        assert_eq!(add(&id, &d, &cc), d);
        assert_eq!(add(&d, &id, &cc), d);
    }

    #[test]
    fn test_deg1_add_different() {
        let d1 = DivisorCoords::deg1(F7::new(1), F7::new(2));
        let d2 = DivisorCoords::deg1(F7::new(3), F7::new(4));
        let cc = make_curve();

        let result = add(&d1, &d2, &cc);
        assert_eq!(result.degree(), 2);
    }

    #[test]
    fn test_deg1_add_same() {
        let d1 = DivisorCoords::deg1(F7::new(2), F7::new(3));
        let d2 = DivisorCoords::deg1(F7::new(2), F7::new(4));

        let result = deg1_add(d1.u0, d1.v0, d2.u0, d2.v0);
        assert!(result.is_identity());
    }

    #[test]
    fn test_double_identity() {
        let cc = make_curve();
        let id = DivisorCoords::identity();
        assert_eq!(double(&id, &cc), id);
    }
}
