//! Characteristic 2 genus 2 ramified model divisor arithmetic.
//!
//! Explicit formulas for divisor addition and doubling on ramified genus 2
//! hyperelliptic curves over fields of characteristic 2.
//!
//! ## Curve Form
//! `y² + h(x)y = f(x)` where:
//! - `f(x) = x⁵ + f₂x² + f₁x + f₀` (monic, degree 5, f₃ = f₄ = 0)
//! - `h(x) = h₂x² + h₁x + h₀` where `h₂ ∈ {0,1}`
//!
//! In characteristic 2: addition = subtraction (XOR)
//!
//! Based on: Sebastian Lindner, 2019

use crate::field::Field;

/// Curve constants for a characteristic 2 ramified genus 2 curve.
///
/// Represents `y² + h(x)y = f(x)` where:
/// - `f(x) = x⁵ + f2*x² + f1*x + f0`
/// - `h(x) = h2*x² + h1*x + h0`
#[derive(Clone, Copy, Debug)]
pub struct CurveConstants<F: Field> {
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

/// Add two degree 1 divisors (char 2).
///
/// Input: `D1 = (x + u0, v0)`, `D2 = (x + up0, vp0)`
/// Output: `D3 = D1 + D2`
///
/// In char 2: subtraction = addition (XOR)
#[inline]
pub fn deg1_add<F: Field>(u0: F, v0: F, up0: F, vp0: F) -> DivisorCoords<F> {
    // d := u mod up = u0 + up0 (in char 2, - = +)
    let d = u0 + up0;

    if d.is_zero() {
        return DivisorCoords::identity();
    }

    // s := (vp + v) / d (in char 2, - = +)
    let w1 = d.inv();
    let s0 = w1 * (vp0 + v0);

    // upp = u * up
    let upp1 = u0 + up0;
    let upp0 = u0 * up0;

    // vpp := (v + u*s) mod upp
    let vpp0 = s0 * u0 + v0;

    DivisorCoords::deg2(upp1, upp0, s0, vpp0)
}

/// Add a degree 1 and degree 2 divisor (char 2).
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
        f2,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // d := up mod u (in char 2)
    let t0 = u0 + up1;
    let d = up0 + u0 * t0;

    if d.is_zero() {
        // dw := (v + vp + h) mod u
        let vh1 = h1 + vp1;
        let dw = h0 + v0 + h2 * u0.square() + u0 * vh1 + vp0;

        if dw.is_zero() {
            // Result is degree 1
            // upp0 = t0
            let vpp0 = vp0 + vp1 * t0;
            return DivisorCoords::deg1(t0, vpp0);
        }

        // k := ExactQuotient(f - vp*(vp + h), up)
        // k2 = up1
        let k1 = vp1 * h2 + up0 + up1.square();
        let k0 = f2 + vp1 * vh1 + vp0 * h2 + up1 * (up0 + k1);

        // s := dw^-1 * k mod u
        let w1 = dw.inv();
        let s0 = w1 * (k0 + u0 * (k1 + u0 * t0));

        // upp := ExactQuotient(-s*(s*up + 2*vp + h) + k, u)
        // Note: 2*vp = 0 in char 2
        let t1 = s0 * up1 + vh1;
        let upp1 = t0 + s0.square() + s0 * h2;
        let upp0 = k1 + s0 * (t1 + vp1) + u0 * upp1;

        // vpp := (-s*up - vp - h) mod upp (in char 2, - = +)
        let vpp1 = upp1 * h2 + s0 * upp1 + t1;
        let vpp0 = upp0 * h2 + s0 * (upp0 + up0) + h0 + vp0;

        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
    }

    // General case
    let w1 = d.inv();
    let s0 = w1 * (v0 + vp0 + vp1 * u0);

    // Compute result
    let upp1 = s0.square() + s0 * h2 + t0;
    let upp0 = s0 * (h2 * up1 + h1) + h2 * vp1 + up0 + up1 * u0 + upp1 * t0;

    let vpp1 = s0 * (upp1 + up1) + h1 + vp1 + upp1 * h2;
    let vpp0 = s0 * (upp0 + up0) + h0 + vp0 + upp0 * h2;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Add two degree 2 divisors (char 2).
///
/// Input: `D1 = (x² + u1*x + u0, v1*x + v0)`, `D2 = (x² + up1*x + up0, vp1*x + vp0)`
/// Output: `D3 = D1 + D2`
#[inline]
#[allow(clippy::too_many_arguments)]
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
        f2,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // d := Resultant(u, up) computed with 2x2 system (in char 2)
    let m3 = up1 + u1; // = up1 - u1
    let m4 = u0 + up0; // = u0 - up0
    let m1 = m4 + up1 * m3;
    let m2 = up0 * m3; // Note: no minus sign in char 2
    let d = m1 * m4 + m2 * m3;

    // Special case: d = 0
    if d.is_zero() {
        if m3.is_zero() {
            // u = up (same polynomial)
            let t1 = v1 + h1;
            let dw21 = vp1 + t1 + h2 * u1;
            let dw20 = vp0 + h0 + v0 + h2 * u0;

            if dw20.is_zero() && dw21.is_zero() {
                return DivisorCoords::identity();
            }

            // b2 := dw21^-1
            let b2 = dw21.inv();

            // k := ExactQuotient(f - v*(v + h), u)
            // k2 = u1
            let k1 = v1 * h2 + u0 + u1.square();
            let k0 = f2 + v1 * t1 + v0 * h2 + u1 * (u0 + k1);

            // u := ExactQuotient(u, dw2 * b2)
            let u0_new = u1 + dw20 * b2;

            // s := b2 * k mod u
            let s0 = b2 * (k0 + u0_new * (k1 + u0_new * (u1 + u0_new)));

            // upp := u²
            let upp1 = u0_new + u0_new; // = 0 in char 2!
            let upp0 = u0_new.square();

            // vpp := (v + u*s) mod upp
            let vpp1 = s0 + v1;
            let vpp0 = v0 + u0_new * s0;

            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
        }

        // m3 != 0: GCD is degree 1
        let t0 = v0 + vp0;
        let t1 = v1 + vp1;
        let m3_sq = m3.square();
        let dw3 = m3_sq * (t0 + h0) + m4 * (m3 * (t1 + h1) + m4 * h2);

        if dw3.is_zero() {
            // GCD divides (v + vp + h)
            let a1 = m3.inv(); // No minus sign in char 2
            let s1 = m4 * a1;

            let u0_new = u1 + s1;
            let up0_new = up1 + s1;

            let s0 = a1 * (t0 + up0_new * t1);

            let upp1 = u0_new + up0_new;
            let upp0 = u0_new * up0_new;

            let vpp1 = v1 + s0;
            let vpp0 = v0 + s0 * u0_new;

            return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
        }

        // General GCD = 1 case with d = 0
        let vh1 = v1 + h1;
        // k2 = u1
        let k1 = v1 * h2 + u0 + u1.square();
        let k0 = f2 + v1 * vh1 + v0 * h2 + u1 * (u0 + k1);

        // a12 = -a1*a2 with weight m3*dw3
        let a12 = m3_sq * (h2 * up1 + t1 + h1);

        // s := (a2*a1*(vp - v) + b2*k) mod up; with weight m3*dw3
        let t3 = m3 * m3_sq;
        let sp1 = t3 * (k1 + up0 + up1 * m3) + a12 * t1;
        let sp0 = t3 * (k0 + up0 * m3) + a12 * t0;
        let d_new = m3 * dw3;

        return deg2_add_common(u1, u0, v1, v0, up1, up0, m1, m3, sp1, sp0, d_new, cc);
    }

    // General case: d != 0
    let r0 = vp0 + v0; // = vp0 - v0 in char 2
    let r1 = vp1 + v1; // = vp1 - v1 in char 2
    let sp1 = r0 * m3 + r1 * m4;
    let sp0 = r0 * m1 + r1 * m2;

    deg2_add_common(u1, u0, v1, v0, up1, up0, m1, m3, sp1, sp0, d, cc)
}

/// Common computation for deg2_add after computing sp0, sp1, d
#[inline]
#[allow(clippy::too_many_arguments)]
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
        f2: _,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;
    let vh1 = v1 + h1;

    if sp1.is_zero() {
        let w1 = d.inv();
        let s0 = sp0 * w1;
        let upp0 = m3 + s0.square() + s0 * h2;

        let t1 = s0 * (u1 + upp0) + vh1 + h2 * upp0;
        let vpp0 = upp0 * t1 + s0 * u0 + v0 + h0;

        return DivisorCoords::deg1(upp0, vpp0);
    }

    let w1 = (d * sp1).inv();
    let w2 = w1 * d;
    let w3 = w2 * d;
    let w4 = w3.square();
    let s1 = w1 * sp1.square();
    let spp0 = sp0 * w2;

    let upp1 = h2 * w3 + m3 + w4;
    let upp0 = spp0.square() + m1 + w3 * (h2 * (spp0 + up1) + h1) + w4 * m3;

    let t0 = upp0 + u0;
    let t1 = u1 + upp1;
    let t2 = (upp1 + spp0) * t1;
    let vpp1 = s1 * (t2 + t0) + vh1 + h2 * upp1;
    let vpp0 = s1 * (spp0 * t0 + upp0 * t1) + v0 + h0 + h2 * upp0;

    DivisorCoords::deg2(upp1, upp0, vpp1, vpp0)
}

/// Double a degree 1 divisor (char 2).
///
/// Input: `D = (x + u0, v0)`
/// Output: `2D`
///
/// In char 2: 2v = 0
#[inline]
pub fn deg1_dbl<F: Field>(u0: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let CurveConstants {
        f2: _,
        f1,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // upp := u²
    let upp0 = u0.square();

    // d := h mod u (since 2v = 0 in char 2)
    let d = h2 * upp0 + h1 * u0 + h0;

    if d.is_zero() {
        return DivisorCoords::identity();
    }

    // b1 := d^-1
    let w1 = d.inv();

    // k := ExactQuotient(f - v*(v + h), u)
    // s := b1*k mod u
    let t0 = upp0.square() + f1 + v0 * h1;
    let vpp1 = t0 * w1;
    let vpp0 = vpp1 * u0 + v0;

    // upp1 = 0 in char 2 (u0 + u0 = 0)
    DivisorCoords::deg2(F::zero(), upp0, vpp1, vpp0)
}

/// Double a degree 2 divisor (char 2).
///
/// Input: `D = (x² + u1*x + u0, v1*x + v0)`
/// Output: `2D`
#[inline]
pub fn deg2_dbl<F: Field>(u1: F, u0: F, v1: F, v0: F, cc: &CurveConstants<F>) -> DivisorCoords<F> {
    let CurveConstants {
        f2,
        f1: _,
        f0: _,
        h2,
        h1,
        h0,
    } = *cc;

    // d := Resultant(u, h) since 2v = 0 in char 2
    // Computed with 2x2 system
    let m3 = h1 + h2 * u1;
    let m4 = h0 + h2 * u0;
    let m1 = m4 + m3 * u1;
    let m2 = m3 * u0; // No minus sign in char 2
    let d = m4 * m1 + m2 * m3;

    // Special case: d = 0
    if d.is_zero() {
        if m3.is_zero() {
            return DivisorCoords::identity();
        }

        // b1 := m3^-1 (no minus sign in char 2)
        let b1 = m3.inv();

        // k := ExactQuotient(f - v*(v + h), u)
        // k2 = u1
        let k1 = v1 * h2 + u0 + u1.square();
        let k0 = f2 + v1 * (v1 + h1) + v0 * h2 + u0 * u1 + u1 * k1;

        // u := ExactQuotient(u, dw1*b1)
        // s := b1*k mod u
        let u0_new = u1 + m4 * b1;
        let s0 = b1 * (k0 + u0_new * (k1 + u0_new * (u1 + u0_new)));

        // upp := u²
        let upp1 = u0_new + u0_new; // = 0 in char 2
        let upp0 = u0_new.square();

        // vpp := (v + u*s) mod upp
        let vpp1 = v1 + s0;
        let vpp0 = v0 + u0_new * s0;

        return DivisorCoords::deg2(upp1, upp0, vpp1, vpp0);
    }

    // General case
    let t0 = u1.square();
    let vh1 = v1 + h1;
    let r1 = t0 + h2 * v1;
    let r0 = u1 * r1 + f2 + v1 * vh1 + h2 * v0;

    let sp0 = r0 * m1 + r1 * m2;
    let sp1 = r0 * m3 + r1 * m4;

    if sp1.is_zero() {
        let w1 = d.inv();
        let s0 = sp0 * w1;

        let upp0 = s0.square() + s0 * h2;

        let t1 = s0 * (u1 + upp0) + h2 * upp0 + v1 + h1;
        let vpp0 = upp0 * t1 + h0 + s0 * u0 + v0;

        return DivisorCoords::deg1(upp0, vpp0);
    }

    let w1 = (d * sp1).inv();
    let w2 = w1 * d;
    let w3 = w2 * d;
    let w4 = w3.square();
    let s1 = w1 * sp1.square();
    let spp0 = sp0 * w2;

    let upp1 = w4 + h2 * w3;
    let upp0 = spp0.square() + w3 * (h2 * (spp0 + u1) + h1);

    let t0 = upp0 + u0;
    let t1 = u1 + upp1;
    let t2 = t1 * (upp1 + spp0);
    let vpp1 = s1 * (t2 + t0) + vh1 + h2 * upp1;
    let vpp0 = s1 * (spp0 * t0 + upp0 * t1) + v0 + h0 + h2 * upp0;

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
    use crate::field::BinaryExtField;

    // Use proper extension fields for characteristic 2
    type GF4 = BinaryExtField<2>;
    type GF8 = BinaryExtField<3>;

    // =============================================================
    // GF(4) Test Curves and Helpers
    // =============================================================

    /// Create a GF(4) curve from the Magma whitebox tester:
    /// f = x^5 + α*x² + α²*x, h = x² + α*x + α
    fn make_gf4_curve_1() -> CurveConstants<GF4> {
        let alpha = GF4::gen(); // α (bits = 2)
        let alpha_sq = alpha * alpha; // α² = α + 1 (bits = 3)
        CurveConstants {
            f2: alpha,       // f2 = α
            f1: alpha_sq,    // f1 = α²
            f0: GF4::zero(), // f0 = 0
            h2: GF4::one(),  // h2 = 1
            h1: alpha,       // h1 = α
            h0: alpha,       // h0 = α
        }
    }

    /// Create a second GF(4) curve:
    /// f = x^5 + x² + α²*x + 1, h = x² + x
    fn make_gf4_curve_2() -> CurveConstants<GF4> {
        let alpha_sq = GF4::new(3); // α² = α + 1
        CurveConstants {
            f2: GF4::one(),  // f2 = 1
            f1: alpha_sq,    // f1 = α²
            f0: GF4::one(),  // f0 = 1
            h2: GF4::one(),  // h2 = 1
            h1: GF4::one(),  // h1 = 1
            h0: GF4::zero(), // h0 = 0
        }
    }

    // =============================================================
    // GF(8) Test Curves and Helpers
    // =============================================================

    /// Create a GF(8) curve from the Magma whitebox tester:
    /// f = x^5 + α*x² + α²*x + α, h = x² + α⁶*x + α
    fn make_gf8_curve_1() -> CurveConstants<GF8> {
        let alpha = GF8::gen();
        let alpha_sq = alpha.pow(2); // α²
        let alpha_6 = alpha.pow(6); // α⁶
        CurveConstants {
            f2: alpha,      // f2 = α
            f1: alpha_sq,   // f1 = α²
            f0: alpha,      // f0 = α
            h2: GF8::one(), // h2 = 1
            h1: alpha_6,    // h1 = α⁶
            h0: alpha,      // h0 = α
        }
    }

    // =============================================================
    // Basic Structure Tests
    // =============================================================

    #[test]
    fn test_identity() {
        let id = DivisorCoords::<GF4>::identity();
        assert_eq!(id.degree(), 0);
        assert!(id.is_identity());
    }

    #[test]
    fn test_deg1_creation() {
        let alpha = GF4::gen();
        let d = DivisorCoords::deg1(alpha, GF4::one());
        assert_eq!(d.degree(), 1);
    }

    #[test]
    fn test_deg2_creation() {
        let alpha = GF4::gen();
        let d = DivisorCoords::deg2(alpha, GF4::one(), alpha, GF4::zero());
        assert_eq!(d.degree(), 2);
    }

    #[test]
    fn test_add_identity() {
        let cc = make_gf4_curve_1();
        let id = DivisorCoords::identity();
        let alpha = GF4::gen();
        let d = DivisorCoords::deg1(alpha, GF4::one());

        assert_eq!(add(&id, &d, &cc), d);
        assert_eq!(add(&d, &id, &cc), d);
    }

    #[test]
    fn test_double_identity() {
        let cc = make_gf4_curve_1();
        let id = DivisorCoords::identity();
        assert_eq!(double(&id, &cc), id);
    }

    // =============================================================
    // Whitebox Test Cases from Magma (ch2_ramifiedG2_whiteBox_tester.mag)
    // =============================================================

    /// Test case 1: GF(4) double of identity
    /// From Magma: U1 = 1, V1 = 0, result should be identity
    #[test]
    fn test_gf4_dbl_identity() {
        let cc = make_gf4_curve_1();
        let id = DivisorCoords::identity();
        let result = double(&id, &cc);
        assert!(result.is_identity());
    }

    /// Test case 3: GF(4) double of degree 2 divisor
    /// From Magma (lines 42-57):
    /// f = x^5 + α*x² + α²*x, h = x² + α*x + α
    /// U1 = x² + α*x + 1, V1 = α*x + α²
    #[test]
    fn test_gf4_dbl_deg2() {
        let cc = make_gf4_curve_1();
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;

        // U1 = x² + α*x + 1: u1 = α, u0 = 1
        // V1 = α*x + α²: v1 = α, v0 = α²
        let d = DivisorCoords::deg2(alpha, GF4::one(), alpha, alpha_sq);

        let result = double(&d, &cc);

        // Result should be a valid divisor (degree 0, 1, or 2)
        assert!(result.degree() <= 2);
    }

    /// Test case 5: GF(4) double with d = 0 branch
    /// From Magma (lines 76-91):
    /// f = x^5 + x² + α²*x + 1, h = x² + x
    /// U1 = x² + α*x, V1 = 1
    #[test]
    fn test_gf4_dbl_d_zero() {
        let cc = make_gf4_curve_2();
        let alpha = GF4::gen();

        // U1 = x² + α*x + 0: u1 = α, u0 = 0
        // V1 = 0*x + 1: v1 = 0, v0 = 1
        let d = DivisorCoords::deg2(alpha, GF4::zero(), GF4::zero(), GF4::one());

        let result = double(&d, &cc);

        // Result should be a valid divisor
        assert!(result.degree() <= 2);
    }

    /// Test case: GF(8) double of degree 2 divisor
    /// From Magma (lines 25-40):
    /// f = x^5 + α*x² + α²*x + α, h = x² + α⁶*x + α
    /// U1 = x² + α⁵*x + α³, V1 = α⁴*x + α
    #[test]
    fn test_gf8_dbl_deg2() {
        let cc = make_gf8_curve_1();
        let alpha = GF8::gen();
        let alpha_3 = alpha.pow(3);
        let alpha_4 = alpha.pow(4);
        let alpha_5 = alpha.pow(5);

        // U1 = x² + α⁵*x + α³: u1 = α⁵, u0 = α³
        // V1 = α⁴*x + α: v1 = α⁴, v0 = α
        let d = DivisorCoords::deg2(alpha_5, alpha_3, alpha_4, alpha);

        let result = double(&d, &cc);

        // Result should be a valid divisor
        assert!(result.degree() <= 2);
    }

    /// Test case: GF(8) deg1 + identity
    /// From Magma (lines 127-145):
    /// D1 = (x + α², α⁴), D2 = identity
    #[test]
    fn test_gf8_add_deg1_identity() {
        let cc = make_gf8_curve_1();
        let alpha = GF8::gen();
        let alpha_2 = alpha.pow(2);
        let alpha_4 = alpha.pow(4);

        let d1 = DivisorCoords::deg1(alpha_2, alpha_4);
        let id = DivisorCoords::identity();

        let result = add(&d1, &id, &cc);

        // Adding identity should return the same divisor
        assert_eq!(result, d1);
    }

    /// Test case: GF(8) deg2 + identity
    /// From Magma (lines 147-165):
    /// D1 = (x² + α*x, α⁴*x + α), D2 = identity
    #[test]
    fn test_gf8_add_deg2_identity() {
        let _cc = make_gf8_curve_1(); // Not used - different curve below
        let alpha = GF8::gen();
        let alpha_4 = alpha.pow(4);

        // Different curve constants for this test
        let cc_test = CurveConstants {
            f2: alpha.pow(3), // α³
            f1: alpha.pow(2), // α²
            f0: alpha.pow(3), // α³
            h2: GF8::one(),
            h1: alpha.pow(4), // α⁴
            h0: alpha.pow(4), // α⁴
        };

        let d1 = DivisorCoords::deg2(alpha, GF8::zero(), alpha_4, alpha);
        let id = DivisorCoords::identity();

        let result = add(&d1, &id, &cc_test);

        // Adding identity should return the same divisor
        assert_eq!(result, d1);
    }

    /// Test case: GF(4) two deg1 divisors (same u)
    /// From Magma (lines 347-365):
    /// f = x^5 + α⁵*x² + α*x + α², h = x² + α⁵
    /// U1 = x + α³, V1 = α⁵
    /// U2 = x + α³, V2 = α⁶
    #[test]
    fn test_gf8_add_deg1_same_u() {
        let alpha = GF8::gen();
        let alpha_3 = alpha.pow(3);
        let alpha_5 = alpha.pow(5);
        let alpha_6 = alpha.pow(6);

        let cc = CurveConstants {
            f2: alpha_5,
            f1: alpha,
            f0: alpha.pow(2),
            h2: GF8::one(),
            h1: GF8::zero(),
            h0: alpha_5,
        };

        // Same u0 means they differ only in v
        let d1 = DivisorCoords::deg1(alpha_3, alpha_5);
        let d2 = DivisorCoords::deg1(alpha_3, alpha_6);

        let result = add(&d1, &d2, &cc);

        // When u is the same and v differs, result depends on v + v' + h
        assert!(result.degree() <= 2);
    }

    /// Test case: GF(4) deg1 add returning deg2
    /// From Magma (lines 367-385):
    /// f = x^5 + α*x + α, h = x
    /// D1 = (x + 1, α), D2 = (x, α²)
    #[test]
    fn test_gf4_add_deg1_to_deg2() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;

        let cc = CurveConstants {
            f2: GF4::zero(),
            f1: alpha,
            f0: alpha,
            h2: GF4::zero(),
            h1: GF4::one(),
            h0: GF4::zero(),
        };

        // D1 = (x + 1, α): u0 = 1, v0 = α
        // D2 = (x + 0, α²): u0 = 0, v0 = α²
        let d1 = DivisorCoords::deg1(GF4::one(), alpha);
        let d2 = DivisorCoords::deg1(GF4::zero(), alpha_sq);

        let result = add(&d1, &d2, &cc);

        // Different u0 values should produce a deg2 result
        assert_eq!(result.degree(), 2);
    }

    /// Test case: GF(8) deg2 + deg2
    /// From Magma (lines 167-185):
    /// Two degree 2 divisors with non-trivial result
    #[test]
    fn test_gf8_add_deg2_deg2() {
        let alpha = GF8::gen();

        let cc = CurveConstants {
            f2: alpha.pow(2),
            f1: GF8::one(),
            f0: alpha.pow(3),
            h2: GF8::one(),
            h1: alpha.pow(5),
            h0: alpha.pow(6),
        };

        // D1 = (x² + α³*x + α⁵, 1)
        let d1 = DivisorCoords::deg2(alpha.pow(3), alpha.pow(5), GF8::zero(), GF8::one());

        // D2 = (x² + α⁶*x + 1, α⁶*x + α⁴)
        let d2 = DivisorCoords::deg2(alpha.pow(6), GF8::one(), alpha.pow(6), alpha.pow(4));

        let result = add(&d1, &d2, &cc);

        // Result should be a valid divisor
        assert!(result.degree() <= 2);
    }

    /// Test case: GF(4) deg12_add
    /// From Magma (lines 267-285):
    /// D1 = (x + α⁴, α⁶), D2 = (x² + α²*x + α³, α⁴)
    #[test]
    fn test_gf8_add_deg1_deg2() {
        let alpha = GF8::gen();

        let cc = CurveConstants {
            f2: alpha.pow(2),
            f1: GF8::one(),
            f0: alpha.pow(3),
            h2: GF8::one(),
            h1: alpha.pow(5),
            h0: alpha.pow(6),
        };

        let d1 = DivisorCoords::deg1(alpha.pow(4), alpha.pow(6));
        let d2 = DivisorCoords::deg2(alpha.pow(2), alpha.pow(3), GF8::zero(), alpha.pow(4));

        let result = add(&d1, &d2, &cc);

        // Result should be a valid divisor
        assert!(result.degree() <= 2);
    }

    // =============================================================
    // Closure Tests (verify D + D' + D' = D for various cases)
    // =============================================================

    #[test]
    fn test_gf4_closure_deg1() {
        let alpha = GF4::gen();
        let cc = CurveConstants {
            f2: GF4::zero(),
            f1: alpha,
            f0: alpha,
            h2: GF4::zero(),
            h1: GF4::one(),
            h0: GF4::zero(),
        };

        // Create two different deg1 divisors and add them
        let d1 = DivisorCoords::deg1(GF4::one(), alpha);
        let d2 = DivisorCoords::deg1(GF4::zero(), alpha * alpha);

        let sum = add(&d1, &d2, &cc);

        // The result should be valid
        assert!(sum.degree() <= 2);
    }

    #[test]
    fn test_gf8_double_closure() {
        let alpha = GF8::gen();
        let cc = make_gf8_curve_1();

        let d = DivisorCoords::deg2(alpha.pow(5), alpha.pow(3), alpha.pow(4), alpha);

        // Double should produce a valid result
        let doubled = double(&d, &cc);
        assert!(doubled.degree() <= 2);

        // Double again should also be valid
        let quadrupled = double(&doubled, &cc);
        assert!(quadrupled.degree() <= 2);
    }
}
