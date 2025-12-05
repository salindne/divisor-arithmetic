//! Generic ramified model divisor arithmetic for arbitrary genus.
//!
//! Ramified curves have one point at infinity:
//! - `y² = f(x)` where `deg(f) = 2g+1`
//!
//! Divisors are represented as `(u, v, w)` where:
//! - `u` is monic with `deg(u) ≤ g`
//! - `deg(v) < deg(u)`
//! - `w = (f - v²) / u`
//!
//! Based on: Sebastian Lindner, 2019

use crate::field::Field;
use crate::poly::Poly;

/// A divisor on a ramified hyperelliptic curve.
///
/// Represented as `(u, v, w)` where `w = (f - v²) / u`.
#[derive(Clone, PartialEq, Eq)]
pub struct Divisor<F: Field> {
    pub u: Poly<F>,
    pub v: Poly<F>,
    pub w: Poly<F>,
}

impl<F: Field + std::fmt::Display> std::fmt::Debug for Divisor<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Divisor")
            .field("u", &self.u)
            .field("v", &self.v)
            .field("w", &self.w)
            .finish()
    }
}

impl<F: Field> Divisor<F> {
    /// Create a new divisor
    pub fn new(u: Poly<F>, v: Poly<F>, w: Poly<F>) -> Self {
        Self { u, v, w }
    }

    /// Create the neutral (identity) divisor for curve f
    pub fn neutral(f: &Poly<F>) -> Self {
        Self {
            u: Poly::constant(F::one()),
            v: Poly::zero(),
            w: f.clone(),
        }
    }

    /// Get the degree of this divisor (degree of u)
    pub fn degree(&self) -> i32 {
        self.u.deg()
    }
}

/// Add two divisors on a ramified curve.
///
/// Implements the specialized Cantor algorithm from `Add_RAM` in Magma.
pub fn add<F: Field>(d1: &Divisor<F>, d2: &Divisor<F>, f: &Poly<F>, g: usize) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;
    let u2 = &d2.u;
    let v2 = &d2.v;

    // S, a1, b1 := XGCD(u1, u2)
    let (s, a1, _b1) = u1.xgcd(u2);

    // K := a1*(v2 - v1) mod u2
    let v_diff = v2 - v1;
    let mut k = (a1 * v_diff.clone()).rem(u2);

    let mut u1 = u1.clone();
    let mut u2 = u2.clone();
    let mut w1 = w1.clone();

    if !s.is_one() {
        // S, a2, b2 := XGCD(S, v2 + v1)
        let v_sum = v2 + v1;
        let (s2, a2, b2) = s.xgcd(&v_sum);

        // K := (a2*K + b2*w1) mod u2
        k = (a2 * k + b2 * w1.clone()).rem(&u2);

        if !s2.is_one() {
            // u1 := u1 / S
            u1 = u1.exact_div(&s2);
            // u2 := u2 / S
            u2 = u2.exact_div(&s2);
            // w1 := w1 * S
            w1 = w1 * s2;
        }
    }

    // T := u1 * K
    let t = &u1 * &k;
    // u := u1 * u2
    let mut u = &u1 * &u2;
    // v := v1 + T
    let mut v = v1 + &t;
    // w := (w1 - K*(v1 + v)) / u2
    let v1_plus_v = v1 + &v;
    let mut w = (w1 - k * v1_plus_v).exact_div(&u2);

    // Normalize or Reduce
    if u.deg() <= g as i32 {
        if v.deg() >= u.deg() {
            let (q, r) = v.div_rem(&u);
            // w := w + q*(v + r)
            let v_plus_r = &v + &r;
            w = w + q * v_plus_r;
            v = r;
        }
    } else {
        while u.deg() > g as i32 {
            let tu = w.clone();
            let neg_v = -v.clone();
            let (q, tv) = neg_v.div_rem(&w);
            // w := u + q*(tv - v)
            let tv_minus_v = &tv - &v;
            w = u + q * tv_minus_v;
            v = tv;
            u = tu;
        }
        // Make u monic
        let lc = u.leading_coeff();
        w = w * lc;
        u = u.make_monic();
    }

    Divisor::new(u, v, w)
}

/// Double a divisor on a ramified curve.
///
/// Implements `Double_RAM` from Magma.
pub fn double<F: Field>(d1: &Divisor<F>, f: &Poly<F>, g: usize) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;

    // t1 := v1 + v1
    let t1 = v1 + v1;

    // S, a1, b1 := XGCD(u1, t1)
    let (s, _a1, b1) = u1.xgcd(&t1);

    // K := b1*w1 mod u1
    let mut k = (b1 * w1.clone()).rem(u1);

    let mut u1 = u1.clone();
    let mut w1 = w1.clone();

    if !s.is_one() {
        // u1 := u1 / S
        u1 = u1.exact_div(&s);
        // w1 := w1 * S
        w1 = w1 * s;
    }

    // T := u1 * K
    let t = &u1 * &k;
    // u := u1^2
    let mut u = &u1 * &u1;
    // v := v1 + T
    let mut v = v1 + &t;
    // w := (w1 - K*(t1 + T)) / u1
    let t1_plus_t = &t1 + &t;
    let mut w = (w1 - k * t1_plus_t).exact_div(&u1);

    // Normalize or Reduce
    if u.deg() <= g as i32 {
        if v.deg() >= u.deg() {
            let (q, r) = v.div_rem(&u);
            let v_plus_r = &v + &r;
            w = w + q * v_plus_r;
            v = r;
        }
    } else {
        while u.deg() > g as i32 {
            let tu = w.clone();
            let neg_v = -v.clone();
            let (q, tv) = neg_v.div_rem(&w);
            let tv_minus_v = &tv - &v;
            w = u + q * tv_minus_v;
            v = tv;
            u = tu;
        }
        let lc = u.leading_coeff();
        w = w * lc;
        u = u.make_monic();
    }

    Divisor::new(u, v, w)
}

/// NUCOMP algorithm for adding two divisors.
///
/// Implements `Nucomp_RAM` from Magma - uses continued fractions for reduction.
pub fn nucomp<F: Field>(d1: &Divisor<F>, d2: &Divisor<F>, f: &Poly<F>, g: usize) -> Divisor<F> {
    // Ensure d1 has larger or equal degree
    let (d1, d2) = if d1.u.deg() < d2.u.deg() {
        (d2, d1)
    } else {
        (d1, d2)
    };

    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;
    let u2 = &d2.u;
    let v2 = &d2.v;

    let t2 = v2 - v1;

    // S, a1, b1 := XGCD(u1, u2)
    let (s, a1, _b1) = u1.xgcd(u2);

    // K := a1*t2 mod u2
    let mut k = (a1 * t2.clone()).rem(u2);

    let mut u1 = u1.clone();
    let mut u2 = u2.clone();
    let mut w1 = w1.clone();

    if !s.is_one() {
        let v_sum = v2 + v1;
        let (s2, a2, b2) = s.xgcd(&v_sum);

        k = (a2 * k + b2 * w1.clone()).rem(&u2);

        if !s2.is_one() {
            u1 = u1.exact_div(&s2);
            u2 = u2.exact_div(&s2);
            w1 = w1 * s2;
        }
    }

    let combined_deg = u2.deg() + u1.deg();

    if combined_deg <= g as i32 {
        // No NUCOMP needed, use regular addition
        let t = &u1 * &k;
        let mut u = &u1 * &u2;
        let mut v = v1 + &t;
        let v1_plus_v = v1 + &v;
        let mut w = (w1 - k.clone() * v1_plus_v).exact_div(&u2);

        if v.deg() >= u.deg() {
            let (q, r) = v.div_rem(&u);
            let v_plus_r = &v + &r;
            w = w + q * v_plus_r;
            v = r;
        }

        return Divisor::new(u, v, w);
    }

    // NUCOMP with continued fractions
    let mut rp = u2.clone();
    let mut r = k.clone();
    let mut cp = Poly::zero();
    let mut c = -Poly::constant(F::one());
    let mut l: i32 = -1;

    let bound = (u2.deg() - u1.deg() + g as i32) / 2;
    while r.deg() > bound {
        let (q, rn) = rp.div_rem(&r);
        rp = r;
        r = rn;
        let cn = cp.clone() - q * c.clone();
        cp = c;
        c = cn;
        l = -l;
    }

    // Reconstruct divisor from continued fraction
    let t3 = &u1 * &r;
    let m1 = (t3.clone() + t2 * c.clone()).exact_div(&u2);
    let v2_plus_v1 = v2 + v1;
    let m2 = (&r * &v2_plus_v1 + &w1 * &c).exact_div(&u2);

    let mut u = if l > 0 {
        &r * &m1 - &c * &m2
    } else {
        &c * &m2 - &r * &m1
    };

    let z = (t3 + cp * u.clone()).exact_div(&c);
    let v = (&z - v1).rem(&u);

    // Make u monic
    u = u.make_monic();

    // Compute w
    let v_sq = &v * &v;
    let w = (f - &v_sq).exact_div(&u);

    Divisor::new(u, v, w)
}

/// NUDUPLE algorithm for doubling a divisor.
///
/// Implements `Nuduple_RAM` from Magma.
pub fn nuduple<F: Field>(d1: &Divisor<F>, f: &Poly<F>, g: usize) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;

    let t2 = v1 + v1;

    // S, a1, b1 := XGCD(u1, t2)
    let (s, _a1, b1) = u1.xgcd(&t2);

    // K := b1*w1 mod u1
    let mut k = (b1 * w1.clone()).rem(u1);

    let mut u1 = u1.clone();
    let mut w1 = w1.clone();

    if !s.is_one() {
        u1 = u1.exact_div(&s);
        w1 = w1 * s;
    }

    if 2 * u1.deg() <= g as i32 {
        // No NUCOMP needed
        let t = &u1 * &k;
        let mut u = &u1 * &u1;
        let mut v = v1 + &t;
        let t2_plus_t = &t2 + &t;
        let mut w = (w1 - k.clone() * t2_plus_t).exact_div(&u1);

        if v.deg() >= u.deg() {
            let (q, r) = v.div_rem(&u);
            let v_plus_r = &v + &r;
            w = w + q * v_plus_r;
            v = r;
        }

        return Divisor::new(u, v, w);
    }

    // NUCOMP with continued fractions
    let mut rp = u1.clone();
    let mut r = k.clone();
    let mut cp = Poly::zero();
    let mut c = -Poly::constant(F::one());
    let mut l: i32 = -1;

    let bound = g as i32 / 2;
    while r.deg() > bound {
        let (q, rn) = rp.div_rem(&r);
        rp = r;
        r = rn;
        let cn = cp.clone() - q * c.clone();
        cp = c;
        c = cn;
        l = -l;
    }

    // Reconstruct divisor
    let m2 = (&r * &t2 + &w1 * &c).exact_div(&u1);

    let r_sq = &r * &r;
    let mut u = if l > 0 {
        r_sq.clone() - &c * &m2
    } else {
        &c * &m2 - r_sq
    };

    let z = (&u1 * &r + cp * u.clone()).exact_div(&c);
    let v = (&z - v1).rem(&u);

    u = u.make_monic();

    let v_sq = &v * &v;
    let w = (f - &v_sq).exact_div(&u);

    Divisor::new(u, v, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;

    type F7 = PrimeField<7>;

    fn poly(coeffs: &[u64]) -> Poly<F7> {
        Poly::from_coeffs(coeffs.iter().map(|&c| F7::new(c)).collect())
    }

    /// Create a genus 2 ramified curve over F7
    /// y² = f(x) where deg(f) = 5
    fn test_curve_g2() -> Poly<F7> {
        // f = x^5 + 3x^4 + 2x^3 + x^2 + 4x + 1
        poly(&[1, 4, 1, 2, 3, 1])
    }

    /// Find a valid divisor by finding a point on the curve.
    /// For x = a, we need f(a) to be a square in F7.
    /// If f(a) = b², then the divisor (x - a, b, w) is valid.
    fn find_valid_divisor(f: &Poly<F7>) -> Option<Divisor<F7>> {
        // Try to find x where f(x) is a quadratic residue
        for a_val in 0..7u64 {
            let a = F7::new(a_val);
            let fa = f.eval(a);

            // Check if fa is a square in F7
            // a is a square mod 7 if a^3 = 1 (or a = 0)
            if fa.is_zero() || fa.pow(3).is_one() {
                // Find square root: in F7, for a ≠ 0, sqrt(a) = a^2 (since a^4 = a for squares)
                // Actually in F7, squares are {0, 1, 2, 4} and sqrt can be computed
                let b = if fa.is_zero() {
                    F7::zero()
                } else {
                    // For p = 7 ≡ 3 (mod 4), sqrt(a) = a^((p+1)/4) = a^2
                    fa.pow(2)
                };

                // Verify it's actually a square root
                if b * b == fa {
                    // u = x - a
                    let u = poly(&[(7 - a_val) % 7, 1]); // x - a = x + (-a)
                    let v = Poly::constant(b);
                    let v_sq = &v * &v;
                    let diff = f - &v_sq;

                    // Check divisibility
                    let (_, rem) = diff.div_rem(&u);
                    if rem.is_zero() {
                        let w = diff.exact_div(&u);
                        return Some(Divisor::new(u, v, w));
                    }
                }
            }
        }
        None
    }

    #[test]
    fn test_neutral() {
        let f = test_curve_g2();
        let neutral = Divisor::neutral(&f);

        assert!(neutral.u.is_one());
        assert!(neutral.v.is_zero());
        assert_eq!(neutral.w, f);
    }

    #[test]
    fn test_add_neutral() {
        let f = test_curve_g2();
        let g = 2;
        let neutral = Divisor::neutral(&f);

        // Create a valid divisor
        let d = find_valid_divisor(&f).expect("Should find valid divisor");

        // D + 0 should give a valid reduced divisor
        let result = add(&d, &neutral, &f, g);
        assert!(result.u.deg() <= g as i32);

        // Verify it's still valid: w = (f - v²) / u
        let v_sq = &result.v * &result.v;
        let expected_w = (f - v_sq).exact_div(&result.u);
        assert_eq!(result.w, expected_w);
    }

    #[test]
    fn test_double_vs_add() {
        let f = test_curve_g2();
        let g = 2;

        // Create a valid divisor
        let d = find_valid_divisor(&f).expect("Should find valid divisor");

        // 2D via double should equal D + D via add
        let doubled = double(&d, &f, g);
        let added = add(&d, &d, &f, g);

        assert_eq!(doubled.u, added.u);
        assert_eq!(doubled.v, added.v);
    }

    #[test]
    fn test_nucomp_vs_add() {
        let f = test_curve_g2();
        let g = 2;

        // Use neutral and double to create divisors
        let d1 = find_valid_divisor(&f).expect("Should find valid divisor");
        let d2 = double(&d1, &f, g);

        // NUCOMP and regular add should give same result
        let via_add = add(&d1, &d2, &f, g);
        let via_nucomp = nucomp(&d1, &d2, &f, g);

        assert_eq!(via_add.u, via_nucomp.u);
        assert_eq!(via_add.v, via_nucomp.v);
    }

    #[test]
    fn test_nuduple_vs_double() {
        let f = test_curve_g2();
        let g = 2;

        let d = find_valid_divisor(&f).expect("Should find valid divisor");

        let via_double = double(&d, &f, g);
        let via_nuduple = nuduple(&d, &f, g);

        assert_eq!(via_double.u, via_nuduple.u);
        assert_eq!(via_double.v, via_nuduple.v);
    }
}
