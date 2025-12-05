//! Arbitrary characteristic genus 2 ramified model divisor arithmetic.
//!
//! Explicit formulas for divisor addition and doubling on ramified genus 2
//! hyperelliptic curves (one point at infinity) over arbitrary fields.
//!
//! ## Curve Form
//! `y² + h(x)y = f(x)` where:
//! - `f(x) = x⁵ + f₄x⁴ + f₃x³ + f₂x² + f₁x + f₀` (monic, degree 5)
//! - `h(x) = h₂x² + h₁x + h₀` where `h₂ ∈ {0,1}`
//!
//! ## Mumford Representation
//! Divisors are represented as `D = (u(x), v(x))` where:
//! - `u(x)` is monic of degree ≤ 2
//! - `deg(v) < deg(u)`
//! - `u | (v² + vh - f)`
//!
//! Based on: Sebastian Lindner, 2019

use crate::field::Field;

/// Curve constants for a ramified genus 2 curve.
///
/// Represents `y² + h(x)y = f(x)` where:
/// - `f(x) = x⁵ + f4*x⁴ + f3*x³ + f2*x² + f1*x + f0`
/// - `h(x) = h2*x² + h1*x + h0`
#[derive(Clone, Copy, Debug)]
pub struct CurveConstants<F: Field> {
    pub f4: F,
    pub f3: F,
    pub f2: F,
    pub f1: F,
    pub f0: F,
    pub h2: F,
    pub h1: F,
    pub h0: F,
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

/// Add two degree 1 divisors.
///
/// Input: `D1 = (x + u0, v0)`, `D2 = (x + up0, vp0)`
/// Output: `D3 = D1 + D2`
#[inline]
pub fn deg1_add<F: Field>(u0: F, v0: F, up0: F, vp0: F) -> DivisorCoords<F> {
    // d := u mod up = u0 - up0
    let d = u0 - up0;

    if d.is_zero() {
        // Same roots - result is identity
        return DivisorCoords::identity();
    }

    // s := (vp - v) / d
    let w1 = d.inv();
    let s0 = w1 * (vp0 - v0);

    // upp = u * up = (x + u0)(x + up0) = x² + (u0+up0)x + u0*up0
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;

    // vpp := (v + u*s) mod upp
    // vpp1 = s0
    let vpp0 = s0 * u0 + v0;

    DivisorCoords::deg2(upp1, upp0, s0, vpp0)
}

/// Add a degree 1 and degree 2 divisor.
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
        f4,
        f3,
        f2,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // d := up mod u = up0 - u0*(up1 - u0)
    let d = up0 - u0 * (up1 - u0);

    if d.is_zero() {
        // dw := (v + vp + h) mod u
        let vh1 = h1 + vp1;
        let dw = h0 + v0 + h2 * u0.square() - u0 * vh1 + vp0;

        if dw.is_zero() {
            // Result is degree 1: u = up / (x + u0)
            let upp0 = up1 - u0;
            let vpp0 = vp0 - vp1 * upp0;
            return DivisorCoords::deg1(upp0, vpp0);
        }

        // k := ExactQuotient(f - vp*(vp + h), up)
        let k2 = f4 - up1;
        let k1 = f3 - vp1 * h2 - up0 - up1 * k2;
        let k0 = f2 - vp1 * vh1 - vp0 * h2 - up0 * k2 - up1 * k1;

        // s := dw^-1 * k mod u
        let w1 = dw.inv();
        let s0 = w1 * (k0 - u0 * (k1 - u0 * (k2 - u0)));

        // upp := ExactQuotient(-s*(s*up + 2*vp + h) + k, u)
        let t0 = s0 * up1 + vh1;
        let upp1 = k2 - s0.square() - s0 * h2 - u0;
        let upp0 = k1 - s0 * (t0 + vp1) - u0 * upp1;

        // vpp := (-s*up - vp - h) mod upp
        let vpp1 = upp1 * h2 + s0 * upp1 - t0;
        let vpp0 = upp0 * h2 + s0 * (upp0 - up0) - h0 - vp0;

        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
    }

    // General case
    // s := d^-1 * (v - vp) mod u
    let w1 = d.inv();
    let s0 = w1 * (v0 - vp0 + vp1 * u0);

    // k := ExactQuotient(f - vp*(vp + h), up)
    let k2 = f4 - up1;
    let k1 = f3 - vp1 * h2 - up0 - up1 * k2;

    // vt := -s*up - vp - h
    let vh1 = h1 + vp1;
    let t0 = s0 * up1 + vh1;
    let upp1 = k2 - s0.square() - s0 * h2 - u0;
    let upp0 = k1 - s0 * (t0 + vp1) - u0 * upp1;

    // vpp := vt mod upp
    let vpp1 = upp1 * h2 + s0 * upp1 - t0;
    let vpp0 = upp0 * h2 + s0 * (upp0 - up0) - h0 - vp0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Add two degree 2 divisors.
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
        f4,
        f3,
        f2,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
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
            let t1 = v1 + h1;
            let dw21 = vp1 + t1 - h2 * u1;
            let dw20 = vp0 + h0 + v0 - h2 * u0;

            if dw20.is_zero() && dw21.is_zero() {
                // vp = -v - h, result is identity
                return DivisorCoords::identity();
            }

            // b2 := dw21^-1
            let b2 = dw21.inv();

            // k := ExactQuotient(f - v*(v + h), u)
            let k2 = f4 - u1;
            let k1 = f3 - v1 * h2 - u0 - u1 * k2;
            let k0 = f2 - v1 * t1 - v0 * h2 - u0 * k2 - u1 * k1;

            // u := ExactQuotient(u, dw2 * b2)
            let u0_new = u1 - dw20 * b2;

            // s := b2 * k mod u
            let s0 = b2 * (k0 - u0_new * (k1 - u0_new * (k2 - u0_new)));

            // upp := u²
            let upp1 = u0_new + u0_new;
            let upp0 = u0_new.square();

            // vpp := (v + u*s) mod upp
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0_new * s0;

            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
        }

        // m3 != 0: GCD is degree 1
        let t1 = v1 + h1;
        let m3_sq = m3.square();
        let dw3 = m3_sq * (vp0 + v0 + h0) + m4 * (m3 * (vp1 + t1) + m4 * h2);

        if dw3.is_zero() {
            // GCD divides (v + vp + h)
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
        let k2 = f4 - u1;
        let k1 = f3 - v1 * h2 - u0 - u1 * k2;
        let k0 = f2 - v1 * t1 - v0 * h2 - u0 * k2 - u1 * k1;

        let a12 = m3_sq * (h2 * up1 - vp1 - t1);

        let t2 = k2 - up1;
        let t0 = m3 * m3_sq;
        let sp1 = t0 * (k1 - up0 - up1 * t2) - a12 * (vp1 - v1);
        let sp0 = t0 * (k0 - up0 * t2) - a12 * (vp0 - v0);
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
    cc: &CurveConstants<F>,
) -> DivisorCoords<F> {
    let CurveConstants {
        f4,
        f3: _,
        f2: _,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    if sp1.is_zero() {
        let t0 = d.inv();
        let s0 = sp0 * t0;
        let upp0 = f4 - u1 - up1 - s0.square() - h2 * s0;

        let t1 = s0 * (u1 - upp0) + h1 + v1 - h2 * upp0;
        let vpp0 = upp0 * t1 - s0 * u0 - v0 - h0;

        return DivisorCoords::deg1(upp0, vpp0);
    }

    let w1 = (d * sp1).inv();
    let w2 = w1 * d;
    let w3 = w2 * d;
    let w4 = w3.square();
    let s1 = w1 * sp1.square();
    let spp0 = sp0 * w2;

    let t3 = spp0 - m3;
    let t2 = t3 - w4 + h2 * w3;
    let t4 = h1 + v1;
    let upp1 = spp0 + t2;
    let upp0 = spp0 * (t3 - m3) + m1 + w3 * (h2 * (spp0 - up1) + t4 + v1) + w4 * (u1 + up1 - f4);

    let t0 = upp0 - u0;
    let t1 = u1 - upp1;
    let vpp1 = s1 * (t1 * t2 + t0) - t4 + h2 * upp1;
    let vpp0 = s1 * (spp0 * t0 + upp0 * t1) - v0 - h0 + h2 * upp0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Double a degree 1 divisor.
///
/// Input: `D = (x + u0, v0)`
/// Output: `2D`
#[inline]
pub fn deg1_dbl<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let CurveConstants {
        f4,
        f3,
        f2,
        f1,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // d := (2v + h) mod u
    let upp0 = u0.square();
    let d = v0 + v0 + h2 * upp0 - h1 * u0 + h0;

    if d.is_zero() {
        return DivisorCoords::identity();
    }

    let w1 = d.inv();
    let upp1 = u0 + u0;

    // k := ExactQuotient(f - v*(v + h), u)
    // s := k / d mod u
    let t1 = upp0 + f3;
    let t2 = t1 - f4 * upp1 + upp0;
    let vpp1 = w1 * (upp0 * (t2 + t2 + t1) - f2 * upp1 + f1 - v0 * (h1 - h2 * upp1));
    let vpp0 = vpp1 * u0 + v0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Double a degree 2 divisor.
///
/// Input: `D = (x² + u1*x + u0, v1*x + v0)`
/// Output: `2D`
#[inline]
pub fn deg2_dbl<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let CurveConstants {
        f4,
        f3,
        f2,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // d := Resultant(u, 2v + h) computed with 2x2 system
    let vh1 = v1 + h1;
    let vh0 = v0 + h0;
    let m3 = h2 * u1 - v1 - vh1;
    let m4 = v0 + vh0 - h2 * u0;
    let m1 = m4 + m3 * u1;
    let m2 = -m3 * u0;
    let d = m4 * m1 - m2 * m3;

    if d.is_zero() {
        if m3.is_zero() {
            return DivisorCoords::identity();
        }

        // b1 := -m3^-1
        let b1 = -m3.inv();

        // k := ExactQuotient(f - v*(v + h), u)
        let k2 = f4 - u1;
        let k1 = f3 - v1 * h2 - u0 - u1 * k2;
        let k0 = f2 - v1 * vh1 - v0 * h2 - u0 * k2 - u1 * k1;

        // u := ExactQuotient(u, dw1 * b1)
        let u0_new = u1 - m4 * b1;
        let s0 = b1 * (k0 - u0_new * (k1 - u0_new * (k2 - u0_new)));

        let upp1 = u0_new + u0_new;
        let upp0 = u0_new.square();

        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0_new * s0;

        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
    }

    // General case
    let t0 = u1.square();
    let t1 = f3 + t0 - h2 * v1;
    let t2 = u0 + u0;
    let t3 = f4 * u1;
    let t4 = f4 * u0;
    let t5 = t0 - t3;
    let t6 = t1 - t2;
    let r1 = t5 + t5 + t6;
    let r0 = u1 * (t2 - t6 + t3) + f2 - v1 * vh1 - t4 - t4 - h2 * v0;

    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;

    if sp1.is_zero() {
        let w1 = d.inv();
        let s0 = sp0 * w1;
        let upp0 = f4 - s0.square() - s0 * h2 - u1 - u1;
        let t1 = s0 * (u1 - upp0) - h2.square() * upp0 + vh1;
        let vpp0 = upp0 * t1 - vh0 - s0 * u0;

        return DivisorCoords::deg1(upp0, vpp0);
    }

    let w1 = (d * sp1).inv();
    let w2 = w1 * d;
    let w3 = w2 * d;
    let w4 = w3.square();
    let s1 = w1 * sp1.square();
    let spp0 = sp0 * w2;

    let t2 = spp0 - w4 + h2 * w3;
    let upp0 = spp0.square() + w3 * (h2 * (spp0 - u1) + v1 + vh1) + w4 * (u1 + u1 - f4);
    let upp1 = spp0 + t2;

    let t0 = upp0 - u0;
    let t1 = u1 - upp1;
    let vpp1 = s1 * (t1 * t2 + t0) - vh1 + h2 * upp1;
    let vpp0 = s1 * (spp0 * t0 + upp0 * t1) - vh0 + h2 * upp0;

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
        // f(x) = x^5 + 3x^4 + 2x^3 + x^2 + 4x + 1
        // h(x) = x + 2
        CurveConstants {
            f4: F7::new(3),
            f3: F7::new(2),
            f2: F7::new(1),
            f1: F7::new(4),
            f0: F7::new(1),
            h2: F7::zero(),
            h1: F7::new(1),
            h0: F7::new(2),
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
        assert!(!d.is_identity());
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
        // Add two degree 1 divisors with different u0 values
        let d1 = DivisorCoords::deg1(F7::new(1), F7::new(2));
        let d2 = DivisorCoords::deg1(F7::new(3), F7::new(4));
        let cc = make_curve();

        let result = add(&d1, &d2, &cc);
        assert_eq!(result.degree(), 2);
    }

    #[test]
    fn test_deg1_add_same() {
        // Add two degree 1 divisors with same u0 (should give identity)
        let d1 = DivisorCoords::deg1(F7::new(2), F7::new(3));
        let d2 = DivisorCoords::deg1(F7::new(2), F7::new(4)); // same u0, different v0

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
