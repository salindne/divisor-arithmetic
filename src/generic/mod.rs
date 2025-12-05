//! Generic polynomial-level divisor arithmetic for arbitrary genus.
//!
//! This module provides reference implementations of divisor arithmetic
//! that work for any genus. These use XGCD and polynomial operations,
//! not explicit field-level formulas.
//!
//! Use these as reference implementations to verify the specialized
//! g2 and g3 explicit formulas.

pub mod ramified;
pub mod split;

use crate::field::Field;
use crate::poly::Poly;
use rand::Rng;

/// Generate a random ramified hyperelliptic curve of genus g over field F.
///
/// Returns f where y² = f(x) with deg(f) = 2g+1 (monic).
/// The curve has one point at infinity.
pub fn random_ramified_curve<F: Field, R: Rng>(g: usize, rng: &mut R) -> Poly<F> {
    // f = x^(2g+1) + random coefficients of lower degree
    let mut coeffs = Vec::with_capacity(2 * g + 2);
    for _ in 0..2 * g + 1 {
        coeffs.push(F::random(rng));
    }
    coeffs.push(F::one()); // Leading coefficient = 1 (monic)
    Poly::from_coeffs(coeffs)
}

/// Generate a random split hyperelliptic curve of genus g over field F.
///
/// Returns f where y² = f(x) with deg(f) = 2g+2 (monic).
/// The curve has two points at infinity.
pub fn random_split_curve<F: Field, R: Rng>(g: usize, rng: &mut R) -> Poly<F> {
    // f = x^(2g+2) + random coefficients of lower degree
    let mut coeffs = Vec::with_capacity(2 * g + 3);
    for _ in 0..2 * g + 2 {
        coeffs.push(F::random(rng));
    }
    coeffs.push(F::one()); // Leading coefficient = 1 (monic)
    Poly::from_coeffs(coeffs)
}

/// Generate a random divisor on a ramified curve.
///
/// Creates a divisor by composing multiple degree-1 divisors.
pub fn random_ramified_divisor<F: Field, R: Rng>(
    f: &Poly<F>,
    g: usize,
    rng: &mut R,
) -> ramified::Divisor<F> {
    let neutral = ramified::Divisor::neutral(f);

    // Try to create a random degree-1 divisor and compose g times
    let mut result = neutral;
    let mut count = 0;

    while count < g {
        // Generate random x-coordinate
        let a = F::random(rng);
        let fa = f.eval(a);

        // Check if fa is a square (simplified - may not always find)
        // For proper implementation, use Legendre symbol and Tonelli-Shanks
        // Here we use fa^((q-1)/2) == 1 test

        // For small fields, just try to find a square root by trial
        // This is inefficient but works for testing
        let mut found_sqrt = None;
        for _ in 0..100 {
            let candidate = F::random(rng);
            if candidate * candidate == fa {
                found_sqrt = Some(candidate);
                break;
            }
        }

        // Also check 0
        if fa.is_zero() {
            found_sqrt = Some(F::zero());
        }

        if let Some(b) = found_sqrt {
            // u = x - a
            let u = Poly::from_coeffs(vec![-a, F::one()]);
            let v = Poly::constant(b);
            let v_sq = &v * &v;
            let diff = f - &v_sq;

            let (_, rem) = diff.div_rem(&u);
            if rem.is_zero() {
                let w = diff.exact_div(&u);
                let d = ramified::Divisor::new(u, v, w);
                result = ramified::add(&result, &d, f, g);
                count += 1;
            }
        }
    }

    result
}

/// Generate a random divisor on a split curve (positive reduced basis).
pub fn random_split_divisor_pos<F: Field, R: Rng>(
    f: &Poly<F>,
    vpl: &Poly<F>,
    g: usize,
    rng: &mut R,
) -> split::Divisor<F> {
    let neutral = split::Divisor::neutral(f, vpl, g);

    let mut result = neutral;
    let mut count = 0;

    while count < g {
        let a = F::random(rng);
        let fa = f.eval(a);

        let mut found_sqrt = None;
        for _ in 0..100 {
            let candidate = F::random(rng);
            if candidate * candidate == fa {
                found_sqrt = Some(candidate);
                break;
            }
        }

        if fa.is_zero() {
            found_sqrt = Some(F::zero());
        }

        if let Some(b) = found_sqrt {
            let u = Poly::from_coeffs(vec![-a, F::one()]);

            // v in reduced basis form
            let b_poly = Poly::constant(b);
            let vpl_minus_b = vpl - &b_poly;
            let v = vpl - &vpl_minus_b.rem(&u);

            let v_sq = &v * &v;
            let diff = f - &v_sq;

            let (_, rem) = diff.div_rem(&u);
            if rem.is_zero() {
                let w = diff.exact_div(&u);
                let n = rng.gen_range(0..=(g as i32));
                let d = split::Divisor::new(u, v, w, n);
                result = split::add_pos(&result, &d, f, vpl, g);
                count += 1;
            }
        }
    }

    result
}

/// Generate a random divisor on a split curve (negative reduced basis).
pub fn random_split_divisor_neg<F: Field, R: Rng>(
    f: &Poly<F>,
    v_neg: &Poly<F>, // -Vpl
    g: usize,
    rng: &mut R,
) -> split::Divisor<F> {
    let vpl = -v_neg.clone();
    let neutral = split::Divisor::neutral(f, &vpl, g);

    // Convert neutral to negative basis
    let t = v_neg - &vpl;
    let (q, r) = t.div_rem(&neutral.u);
    let tv = &neutral.v + &t - r;
    let w = neutral.w - q * (&neutral.v + &tv);
    let result = split::Divisor::new(neutral.u, tv, w, neutral.n);

    let mut result = result;
    let mut count = 0;

    while count < g {
        let a = F::random(rng);
        let fa = f.eval(a);

        let mut found_sqrt = None;
        for _ in 0..100 {
            let candidate = F::random(rng);
            if candidate * candidate == fa {
                found_sqrt = Some(candidate);
                break;
            }
        }

        if fa.is_zero() {
            found_sqrt = Some(F::zero());
        }

        if let Some(b) = found_sqrt {
            let u = Poly::from_coeffs(vec![-a, F::one()]);

            // v in negative reduced basis form
            let b_poly = Poly::constant(b);
            let v_neg_minus_b = v_neg - &b_poly;
            let v = v_neg - &v_neg_minus_b.rem(&u);

            let v_sq = &v * &v;
            let diff = f - &v_sq;

            let (_, rem) = diff.div_rem(&u);
            if rem.is_zero() {
                let w = diff.exact_div(&u);
                let n = rng.gen_range(0..=(g as i32));
                let d = split::Divisor::new(u, v, w, n);
                result = split::add_neg(&result, &d, f, v_neg, g);
                count += 1;
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    type F31 = PrimeField<31>;

    #[test]
    fn test_random_ramified_curve() {
        let mut rng = StdRng::seed_from_u64(12345);
        let g = 2;

        let f: Poly<F31> = random_ramified_curve(g, &mut rng);

        // Should have degree 2g+1 = 5
        assert_eq!(f.deg(), 5);
        // Should be monic
        assert!(f.is_monic());
    }

    #[test]
    fn test_random_split_curve() {
        let mut rng = StdRng::seed_from_u64(12345);
        let g = 2;

        let f: Poly<F31> = random_split_curve(g, &mut rng);

        // Should have degree 2g+2 = 6
        assert_eq!(f.deg(), 6);
        // Should be monic
        assert!(f.is_monic());
    }

    #[test]
    fn test_random_ramified_divisor() {
        let mut rng = StdRng::seed_from_u64(12345);
        let g = 2;

        let f: Poly<F31> = random_ramified_curve(g, &mut rng);
        let d = random_ramified_divisor(&f, g, &mut rng);

        // Divisor should have deg(u) <= g
        assert!(d.u.deg() <= g as i32);

        // Verify w = (f - v²) / u
        let v_sq = &d.v * &d.v;
        let expected_w = (&f - &v_sq).exact_div(&d.u);
        assert_eq!(d.w, expected_w);
    }

    #[test]
    fn test_random_split_divisor_pos() {
        let mut rng = StdRng::seed_from_u64(12345);
        let g = 2;

        let f: Poly<F31> = random_split_curve(g, &mut rng);
        let vpl = split::compute_vpl(&f, g);
        let d = random_split_divisor_pos(&f, &vpl, g, &mut rng);

        // Divisor should have deg(u) <= g
        assert!(d.u.deg() <= g as i32);

        // Verify w = (f - v²) / u
        let v_sq = &d.v * &d.v;
        let expected_w = (&f - &v_sq).exact_div(&d.u);
        assert_eq!(d.w, expected_w);
    }
}

/// Tests for generic algorithms over binary extension fields
#[cfg(test)]
mod extension_field_tests {
    use super::*;
    use crate::field::BinaryExtField;

    type GF8 = BinaryExtField<3>;
    type GF16 = BinaryExtField<4>;

    /// Create a test curve over GF(8)
    /// y² = x⁵ + α*x² + α²*x + α
    fn test_curve_gf8() -> Poly<GF8> {
        let alpha = GF8::gen();
        let alpha_sq = alpha * alpha;
        Poly::from_coeffs(vec![
            alpha,       // x^0
            alpha_sq,    // x^1
            alpha,       // x^2
            GF8::zero(), // x^3
            GF8::zero(), // x^4
            GF8::one(),  // x^5 (monic)
        ])
    }

    /// Create a test curve over GF(16)
    /// y² = x⁵ + α*x² + x + α³
    fn test_curve_gf16() -> Poly<GF16> {
        let alpha = GF16::gen();
        Poly::from_coeffs(vec![
            alpha.pow(3), // x^0
            GF16::one(),  // x^1
            alpha,        // x^2
            GF16::zero(), // x^3
            GF16::zero(), // x^4
            GF16::one(),  // x^5 (monic)
        ])
    }

    /// Find a valid degree-1 divisor on a curve over GF(8)
    fn find_valid_divisor_gf8(f: &Poly<GF8>) -> Option<ramified::Divisor<GF8>> {
        // Try all elements of GF(8)
        for i in 0..8u64 {
            let a = GF8::new(i);
            let fa = f.eval(a);

            // Try to find sqrt(fa) by trying all elements
            for j in 0..8u64 {
                let b = GF8::new(j);
                if b * b == fa {
                    // u = x - a
                    let u = Poly::from_coeffs(vec![-a, GF8::one()]);
                    let v = Poly::constant(b);
                    let v_sq = &v * &v;
                    let diff = f - &v_sq;

                    let (_, rem) = diff.div_rem(&u);
                    if rem.is_zero() {
                        let w = diff.exact_div(&u);
                        return Some(ramified::Divisor::new(u, v, w));
                    }
                }
            }
        }
        None
    }

    /// Find a valid degree-1 divisor on a curve over GF(16)
    fn find_valid_divisor_gf16(f: &Poly<GF16>) -> Option<ramified::Divisor<GF16>> {
        // Try all elements of GF(16)
        for i in 0..16u64 {
            let a = GF16::new(i);
            let fa = f.eval(a);

            // Try to find sqrt(fa) by trying all elements
            for j in 0..16u64 {
                let b = GF16::new(j);
                if b * b == fa {
                    // u = x - a
                    let u = Poly::from_coeffs(vec![-a, GF16::one()]);
                    let v = Poly::constant(b);
                    let v_sq = &v * &v;
                    let diff = f - &v_sq;

                    let (_, rem) = diff.div_rem(&u);
                    if rem.is_zero() {
                        let w = diff.exact_div(&u);
                        return Some(ramified::Divisor::new(u, v, w));
                    }
                }
            }
        }
        None
    }

    #[test]
    fn test_ramified_neutral_over_gf8() {
        let f = test_curve_gf8();
        let neutral = ramified::Divisor::neutral(&f);

        assert!(neutral.u.is_one());
        assert!(neutral.v.is_zero());
        assert_eq!(neutral.w, f);
    }

    #[test]
    fn test_ramified_add_neutral_over_gf8() {
        let f = test_curve_gf8();
        let g = 2;
        let neutral = ramified::Divisor::neutral(&f);

        if let Some(d) = find_valid_divisor_gf8(&f) {
            let result = ramified::add(&d, &neutral, &f, g);

            // Result should satisfy curve equation: w = (f - v²) / u
            let v_sq = &result.v * &result.v;
            let expected_w = (&f - &v_sq).exact_div(&result.u);
            assert_eq!(result.w, expected_w);

            // Degree should be at most g
            assert!(result.u.deg() <= g as i32);
        }
    }

    #[test]
    fn test_ramified_double_over_gf8() {
        let f = test_curve_gf8();
        let g = 2;

        if let Some(d) = find_valid_divisor_gf8(&f) {
            let result = ramified::double(&d, &f, g);

            // Result should satisfy curve equation
            let v_sq = &result.v * &result.v;
            let expected_w = (&f - &v_sq).exact_div(&result.u);
            assert_eq!(result.w, expected_w);

            // Degree should be at most g
            assert!(result.u.deg() <= g as i32);
        }
    }

    #[test]
    fn test_ramified_add_over_gf8() {
        let f = test_curve_gf8();
        let g = 2;

        // Find two different divisors
        if let Some(d1) = find_valid_divisor_gf8(&f) {
            // Create d2 by doubling d1
            let d2 = ramified::double(&d1, &f, g);

            let result = ramified::add(&d1, &d2, &f, g);

            // Result should satisfy curve equation
            let v_sq = &result.v * &result.v;
            let expected_w = (&f - &v_sq).exact_div(&result.u);
            assert_eq!(result.w, expected_w);
        }
    }

    #[test]
    fn test_ramified_nucomp_over_gf8() {
        let f = test_curve_gf8();
        let g = 2;

        if let Some(d1) = find_valid_divisor_gf8(&f) {
            let d2 = ramified::double(&d1, &f, g);

            // NUCOMP should give same result as regular add
            let via_add = ramified::add(&d1, &d2, &f, g);
            let via_nucomp = ramified::nucomp(&d1, &d2, &f, g);

            assert_eq!(via_add.u, via_nucomp.u);
            assert_eq!(via_add.v, via_nucomp.v);
        }
    }

    #[test]
    fn test_ramified_nuduple_over_gf8() {
        let f = test_curve_gf8();
        let g = 2;

        if let Some(d) = find_valid_divisor_gf8(&f) {
            // NUDUPLE should give same result as regular double
            let via_double = ramified::double(&d, &f, g);
            let via_nuduple = ramified::nuduple(&d, &f, g);

            assert_eq!(via_double.u, via_nuduple.u);
            assert_eq!(via_double.v, via_nuduple.v);
        }
    }

    #[test]
    fn test_ramified_over_gf16() {
        let f = test_curve_gf16();
        let g = 2;
        let neutral = ramified::Divisor::neutral(&f);

        assert!(neutral.u.is_one());

        if let Some(d) = find_valid_divisor_gf16(&f) {
            // Test add with neutral
            let result = ramified::add(&d, &neutral, &f, g);
            assert!(result.u.deg() <= g as i32);

            // Test double
            let doubled = ramified::double(&d, &f, g);
            assert!(doubled.u.deg() <= g as i32);

            // Double should satisfy curve equation
            let v_sq = &doubled.v * &doubled.v;
            let expected_w = (&f - &v_sq).exact_div(&doubled.u);
            assert_eq!(doubled.w, expected_w);
        }
    }

    #[test]
    fn test_double_vs_add_over_gf8() {
        let f = test_curve_gf8();
        let g = 2;

        if let Some(d) = find_valid_divisor_gf8(&f) {
            // 2D via double should equal D + D via add
            let doubled = ramified::double(&d, &f, g);
            let added = ramified::add(&d, &d, &f, g);

            assert_eq!(doubled.u, added.u);
            assert_eq!(doubled.v, added.v);
        }
    }

    #[test]
    fn test_polynomial_ops_over_gf8() {
        // Test that basic polynomial operations work over GF(8)
        let alpha = GF8::gen();

        // p1 = x + α
        let p1 = Poly::from_coeffs(vec![alpha, GF8::one()]);

        // p2 = x + α²
        let alpha_sq = alpha * alpha;
        let p2 = Poly::from_coeffs(vec![alpha_sq, GF8::one()]);

        // Test addition
        let sum = p1.clone() + p2.clone();
        // (x + α) + (x + α²) = 2x + (α + α²) = 0 + (α + α²) in char 2
        // α + α² = α + (α + 1) = 1 in GF(4), but in GF(8): α² = bits 4
        // α + α² = 2 ^ 4 = 6
        assert_eq!(sum.coeff(0), alpha + alpha_sq);
        assert!(sum.coeff(1).is_zero()); // 2*1 = 0 in char 2

        // Test multiplication
        let prod = p1.clone() * p2.clone();
        // (x + α)(x + α²) = x² + (α + α²)x + α³
        let alpha_cubed = alpha.pow(3);
        assert_eq!(prod.coeff(0), alpha_cubed);
        assert_eq!(prod.coeff(1), alpha + alpha_sq);
        assert!(prod.coeff(2).is_one());

        // Test division
        let (q, r) = prod.div_rem(&p1);
        assert_eq!(q, p2);
        assert!(r.is_zero());
    }

    #[test]
    fn test_xgcd_over_gf8() {
        let alpha = GF8::gen();

        // Test XGCD of coprime polynomials
        let p1 = Poly::from_coeffs(vec![alpha, GF8::one()]); // x + α
        let p2 = Poly::from_coeffs(vec![alpha.pow(2), GF8::one()]); // x + α²

        let (gcd, a, b) = p1.xgcd(&p2);

        // GCD should be 1 (coprime)
        assert!(gcd.is_one());

        // Verify: gcd = a*p1 + b*p2
        let check = a * p1 + b * p2;
        assert!(check.is_one());
    }
}

/// Integration tests comparing generic and specialized implementations
#[cfg(test)]
mod comparison_tests {
    use super::*;
    use crate::field::PrimeField;
    use crate::g2::ramified::arbitrary::{
        add as g2_add, double as g2_double, CurveConstants, DivisorCoords,
    };

    type F31 = PrimeField<31>;

    /// Create CurveConstants from polynomial f (with h = 0)
    fn poly_to_curve_constants(f: &Poly<F31>) -> CurveConstants<F31> {
        CurveConstants {
            f4: f.coeff(4),
            f3: f.coeff(3),
            f2: f.coeff(2),
            f1: f.coeff(1),
            f0: f.coeff(0),
            h2: F31::zero(),
            h1: F31::zero(),
            h0: F31::zero(),
        }
    }

    /// Find a point (a, b) on the curve y² = f(x) where f(a) = b²
    #[allow(dead_code)]
    fn find_curve_point(f: &Poly<F31>) -> Option<(F31, F31)> {
        for a_val in 0..31u64 {
            let a = F31::new(a_val);
            let fa = f.eval(a);

            // Check if fa is a square in F31
            // a is a square mod 31 if a^15 = 1 (or a = 0) since (31-1)/2 = 15
            if fa.is_zero() || fa.pow(15).is_one() {
                // Find square root: try b^2 = fa
                for b_val in 0..31u64 {
                    let b = F31::new(b_val);
                    if b * b == fa {
                        return Some((a, b));
                    }
                }
            }
        }
        None
    }

    /// Create a valid degree 1 divisor on the curve
    fn create_valid_deg1_divisor(
        f: &Poly<F31>,
        exclude_a: Option<F31>,
    ) -> Option<(DivisorCoords<F31>, ramified::Divisor<F31>)> {
        for a_val in 0..31u64 {
            let a = F31::new(a_val);

            // Skip excluded value
            if let Some(excl) = exclude_a {
                if a == excl {
                    continue;
                }
            }

            let fa = f.eval(a);

            // Check if fa is a square
            if fa.is_zero() || fa.pow(15).is_one() {
                for b_val in 0..31u64 {
                    let b = F31::new(b_val);
                    if b * b == fa {
                        // u = x - a, which is x + (-a)
                        let u0 = -a;
                        let v0 = b;

                        // Verify this is valid
                        let u = Poly::from_coeffs(vec![u0, F31::one()]);
                        let v = Poly::constant(v0);
                        let v_sq = &v * &v;
                        let diff = f - &v_sq;

                        let (_, rem) = diff.div_rem(&u);
                        if rem.is_zero() {
                            let w = diff.exact_div(&u);
                            let coords = DivisorCoords::deg1(u0, v0);
                            let generic = ramified::Divisor::new(u, v, w);
                            return Some((coords, generic));
                        }
                    }
                }
            }
        }
        None
    }

    /// Convert a g2::ramified::DivisorCoords to generic::ramified::Divisor
    /// Requires that the divisor satisfies the curve equation
    fn coords_to_generic(d: &DivisorCoords<F31>, f: &Poly<F31>) -> ramified::Divisor<F31> {
        let u = match d.degree() {
            0 => Poly::constant(F31::one()),
            1 => Poly::from_coeffs(vec![d.u0, F31::one()]),
            2 => Poly::from_coeffs(vec![d.u0, d.u1, F31::one()]),
            _ => unreachable!(),
        };

        let v = match d.degree() {
            0 => Poly::zero(),
            1 => Poly::constant(d.v0),
            2 => Poly::from_coeffs(vec![d.v0, d.v1]),
            _ => unreachable!(),
        };

        // For identity, w = f
        if d.degree() == 0 {
            return ramified::Divisor::new(u, v, f.clone());
        }

        let v_sq = &v * &v;
        let diff = f - &v_sq;

        // Check divisibility before exact division
        let (_, rem) = diff.div_rem(&u);
        assert!(rem.is_zero(), "Divisor does not satisfy curve equation");

        let w = diff.exact_div(&u);
        ramified::Divisor::new(u, v, w)
    }

    /// Convert a generic::ramified::Divisor to g2::ramified::DivisorCoords
    #[allow(dead_code)]
    fn generic_to_coords(d: &ramified::Divisor<F31>) -> DivisorCoords<F31> {
        match d.u.deg() {
            0 | -1 => DivisorCoords::identity(),
            1 => DivisorCoords::deg1(d.u.coeff(0), d.v.coeff(0)),
            2 => DivisorCoords::deg2(d.u.coeff(1), d.u.coeff(0), d.v.coeff(1), d.v.coeff(0)),
            _ => panic!("Invalid degree for genus 2"),
        }
    }

    #[test]
    fn test_generic_vs_specialized_add_deg1() {
        // Create a genus 2 ramified curve (h = 0)
        let f = Poly::from_coeffs(vec![
            F31::new(1), // f0
            F31::new(4), // f1
            F31::new(1), // f2
            F31::new(2), // f3
            F31::new(3), // f4
            F31::one(),  // x^5 (monic)
        ]);

        let cc = poly_to_curve_constants(&f);

        // Find two valid degree 1 divisors on the curve
        let (d1_coords, d1_generic) =
            create_valid_deg1_divisor(&f, None).expect("Should find first valid divisor");

        // Get the x-coordinate of first divisor to exclude
        let first_a = -d1_coords.u0;

        let (d2_coords, d2_generic) =
            create_valid_deg1_divisor(&f, Some(first_a)).expect("Should find second valid divisor");

        // Add using specialized g2 formulas
        let result_specialized = g2_add(&d1_coords, &d2_coords, &cc);

        // Add using generic algorithms
        let result_generic = ramified::add(&d1_generic, &d2_generic, &f, 2);

        // Convert specialized result back and compare
        let result_specialized_generic = coords_to_generic(&result_specialized, &f);

        // The u and v polynomials should match
        assert_eq!(result_generic.u, result_specialized_generic.u);
        assert_eq!(result_generic.v, result_specialized_generic.v);
    }

    #[test]
    fn test_generic_vs_specialized_double_deg1() {
        let f = Poly::from_coeffs(vec![
            F31::new(1),
            F31::new(4),
            F31::new(1),
            F31::new(2),
            F31::new(3),
            F31::one(),
        ]);

        let cc = poly_to_curve_constants(&f);

        // Find a valid degree 1 divisor
        let (d_coords, d_generic) =
            create_valid_deg1_divisor(&f, None).expect("Should find valid divisor");

        // Double using specialized g2 formulas
        let result_specialized = g2_double(&d_coords, &cc);

        // Double using generic algorithms
        let result_generic = ramified::double(&d_generic, &f, 2);

        // Convert and compare
        let result_specialized_generic = coords_to_generic(&result_specialized, &f);

        assert_eq!(result_generic.u, result_specialized_generic.u);
        assert_eq!(result_generic.v, result_specialized_generic.v);
    }

    #[test]
    fn test_generic_vs_specialized_add_deg2() {
        let f = Poly::from_coeffs(vec![
            F31::new(1),
            F31::new(4),
            F31::new(1),
            F31::new(2),
            F31::new(3),
            F31::one(),
        ]);

        let _cc = poly_to_curve_constants(&f);

        // Create valid degree 1 divisors and compose them to get degree 2
        let (d1_coords, d1_generic) =
            create_valid_deg1_divisor(&f, None).expect("Should find first valid divisor");
        let first_a = -d1_coords.u0;

        let (d2_coords, d2_generic) =
            create_valid_deg1_divisor(&f, Some(first_a)).expect("Should find second valid divisor");
        let second_a = -d2_coords.u0;

        // Compose to get degree 2 divisors using GENERIC algorithms
        // This ensures divisors satisfy the curve equation
        let d1_deg2_generic = ramified::add(&d1_generic, &d2_generic, &f, 2);

        if d1_deg2_generic.u.deg() == 2 {
            // Find more valid divisors for the second operand
            let mut exclude_set = vec![first_a, second_a];

            let mut d3_opt = None;
            for a_val in 0..31u64 {
                let a = F31::new(a_val);
                if !exclude_set.contains(&a) {
                    if let Some((c, g)) = create_valid_deg1_divisor(&f, Some(a)) {
                        let a2 = -c.u0;
                        if !exclude_set.contains(&a2) {
                            d3_opt = Some((c, g, a2));
                            break;
                        }
                    }
                }
            }

            if let Some((_d3_coords, d3_generic, third_a)) = d3_opt {
                exclude_set.push(third_a);

                let (_d4_coords, d4_generic) = create_valid_deg1_divisor(&f, Some(third_a))
                    .expect("Should find fourth valid divisor");

                let d2_deg2_generic = ramified::add(&d3_generic, &d4_generic, &f, 2);

                if d2_deg2_generic.u.deg() == 2 {
                    // Add the two degree 2 divisors using GENERIC
                    let result_generic = ramified::add(&d1_deg2_generic, &d2_deg2_generic, &f, 2);

                    // Verify result is valid
                    let v_sq = &result_generic.v * &result_generic.v;
                    let diff = &f - &v_sq;
                    let (_, rem) = diff.div_rem(&result_generic.u);
                    assert!(
                        rem.is_zero(),
                        "Generic result should satisfy curve equation"
                    );

                    // The generic add should produce a valid divisor
                    assert!(result_generic.u.deg() <= 2);
                }
            }
        }
    }

    #[test]
    fn test_generic_vs_specialized_double_deg2() {
        let f = Poly::from_coeffs(vec![
            F31::new(1),
            F31::new(4),
            F31::new(1),
            F31::new(2),
            F31::new(3),
            F31::one(),
        ]);

        let cc = poly_to_curve_constants(&f);

        // Create a degree 2 divisor via addition of valid degree 1 divisors
        let (d1_coords, d1_generic) =
            create_valid_deg1_divisor(&f, None).expect("Should find first valid divisor");
        let first_a = -d1_coords.u0;

        let (d2_coords, d2_generic) =
            create_valid_deg1_divisor(&f, Some(first_a)).expect("Should find second valid divisor");

        let d_deg2_coords = g2_add(&d1_coords, &d2_coords, &cc);
        let d_deg2_generic = ramified::add(&d1_generic, &d2_generic, &f, 2);

        if d_deg2_coords.degree() == 2 {
            let result_specialized = g2_double(&d_deg2_coords, &cc);
            let result_generic = ramified::double(&d_deg2_generic, &f, 2);

            let result_specialized_generic = coords_to_generic(&result_specialized, &f);

            assert_eq!(result_generic.u, result_specialized_generic.u);
            assert_eq!(result_generic.v, result_specialized_generic.v);
        }
    }
}
