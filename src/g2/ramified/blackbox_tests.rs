//! Blackbox tests for genus 2 ramified divisor arithmetic.
//!
//! These tests randomly generate valid divisors on random curves over progressively
//! larger fields, similar to the Magma random testers. All divisors satisfy the
//! curve equation and results are verified against each other.

use crate::field::{BinaryExtField, Field, PrimeField};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

// =============================================================================
// NCH2 Tests - y² = f(x) curves
// =============================================================================

mod not_char2_blackbox {
    use super::*;
    use crate::g2::ramified::not_char2::{add, double, CurveConstants, DivisorCoords};

    /// Try to find b such that b² = target in the field
    fn sqrt<F: Field>(target: F) -> Option<F> {
        if target.is_zero() {
            return Some(F::zero());
        }

        // For small fields, try all elements
        // This works for fields up to ~65536 elements
        let mut rng = rand::thread_rng();
        for _ in 0..10000 {
            let b = F::random(&mut rng);
            if b * b == target {
                return Some(b);
            }
        }
        None
    }

    /// Find a valid point (a, b) on curve y² = f(x)
    fn find_curve_point<F: Field, R: Rng>(f: &[F; 6], rng: &mut R) -> Option<(F, F)> {
        for _ in 0..100 {
            let a = F::random(rng);
            // Evaluate f(a) = a^5 + f4*a^4 + f3*a^3 + f2*a^2 + f1*a + f0
            let a2 = a * a;
            let a3 = a2 * a;
            let a4 = a3 * a;
            let a5 = a4 * a;
            let fa = a5 + f[4] * a4 + f[3] * a3 + f[2] * a2 + f[1] * a + f[0];

            if let Some(b) = sqrt(fa) {
                return Some((a, b));
            }
        }
        None
    }

    /// Create a degree-1 divisor from point (a, b): u = x - a, v = b
    fn deg1_from_point<F: Field>(a: F, b: F) -> DivisorCoords<F> {
        // u = x + (-a), so u0 = -a
        DivisorCoords::deg1(-a, b)
    }

    /// Generate a random valid divisor by composing curve points
    fn random_valid_divisor<F: Field, R: Rng>(
        f: &[F; 6],
        cc: &CurveConstants<F>,
        rng: &mut R,
    ) -> DivisorCoords<F> {
        // Randomly choose degree 0, 1, or 2
        match rng.gen_range(0..10) {
            0 => DivisorCoords::identity(),
            1..=3 => {
                // Degree 1: single point
                if let Some((a, b)) = find_curve_point(f, rng) {
                    deg1_from_point(a, b)
                } else {
                    DivisorCoords::identity()
                }
            }
            _ => {
                // Degree 2: add two points
                if let Some((a1, b1)) = find_curve_point(f, rng) {
                    if let Some((a2, b2)) = find_curve_point(f, rng) {
                        let d1 = deg1_from_point(a1, b1);
                        let d2 = deg1_from_point(a2, b2);
                        add(&d1, &d2, cc)
                    } else {
                        deg1_from_point(a1, b1)
                    }
                } else {
                    DivisorCoords::identity()
                }
            }
        }
    }

    fn random_nch2_curve<F: Field, R: Rng>(rng: &mut R) -> ([F; 6], CurveConstants<F>) {
        let f = [
            F::random(rng), // f0
            F::random(rng), // f1
            F::random(rng), // f2
            F::random(rng), // f3
            F::zero(),      // f4 = 0 for nch2
            F::one(),       // leading coeff (monic)
        ];
        let cc = CurveConstants {
            f3: f[3],
            f2: f[2],
            f1: f[1],
            f0: f[0],
        };
        (f, cc)
    }

    fn run_test<F: Field>(seed: u64, curves: usize, divisors_per_curve: usize) {
        let mut rng = StdRng::seed_from_u64(seed);

        for _ in 0..curves {
            let (f, cc): ([F; 6], CurveConstants<F>) = random_nch2_curve(&mut rng);

            for _ in 0..divisors_per_curve {
                let d1 = random_valid_divisor(&f, &cc, &mut rng);
                let d2 = random_valid_divisor(&f, &cc, &mut rng);

                // Test identity
                let id = DivisorCoords::identity();
                assert_eq!(add(&d1, &id, &cc), d1, "D + 0 != D");
                assert_eq!(add(&id, &d1, &cc), d1, "0 + D != D");

                // Test double of identity
                assert!(double(&id, &cc).is_identity(), "2*0 != 0");

                // Test add with distinct divisors (like Magma tests)
                if d1 != d2 {
                    let sum = add(&d1, &d2, &cc);
                    assert!(sum.degree() <= 2, "ADD result degree > 2");

                    // Test commutativity
                    let sum_rev = add(&d2, &d1, &cc);
                    assert_eq!(sum, sum_rev, "D1 + D2 != D2 + D1");
                }

                // Test doubling
                let dbl = double(&d1, &cc);
                assert!(dbl.degree() <= 2, "DBL result degree > 2");
            }
        }
    }

    #[test]
    fn test_nch2_f7() {
        run_test::<PrimeField<7>>(1001, 10, 2500);
    }

    #[test]
    fn test_nch2_f31() {
        run_test::<PrimeField<31>>(1002, 10, 2500);
    }

    #[test]
    fn test_nch2_f127() {
        run_test::<PrimeField<127>>(1003, 5, 2500);
    }

    #[test]
    fn test_nch2_f8191() {
        run_test::<PrimeField<8191>>(1004, 3, 2500);
    }
}

// =============================================================================
// ARB Tests - y² + h(x)y = f(x) curves
// =============================================================================

mod arbitrary_blackbox {
    use super::*;
    use crate::g2::ramified::arbitrary::{add, double, CurveConstants, DivisorCoords};

    /// Try to find b such that b² + h*b = target (quadratic in char != 2)
    /// Completing the square: (b + h/2)² = target + (h/2)²
    fn solve_quadratic<F: Field>(h_val: F, target: F) -> Option<F> {
        if target.is_zero() && h_val.is_zero() {
            return Some(F::zero());
        }

        // Try random values
        let mut rng = rand::thread_rng();
        for _ in 0..10000 {
            let b = F::random(&mut rng);
            if b * b + h_val * b == target {
                return Some(b);
            }
        }
        None
    }

    /// Find a valid point (a, b) on curve y² + h(a)y = f(a)
    fn find_curve_point<F: Field, R: Rng>(f: &[F; 6], h: &[F; 3], rng: &mut R) -> Option<(F, F)> {
        for _ in 0..100 {
            let a = F::random(rng);

            // Evaluate f(a)
            let a2 = a * a;
            let a3 = a2 * a;
            let a4 = a3 * a;
            let a5 = a4 * a;
            let fa = a5 + f[4] * a4 + f[3] * a3 + f[2] * a2 + f[1] * a + f[0];

            // Evaluate h(a)
            let ha = h[2] * a2 + h[1] * a + h[0];

            // Solve b² + ha*b = fa
            if let Some(b) = solve_quadratic(ha, fa) {
                return Some((a, b));
            }
        }
        None
    }

    fn deg1_from_point<F: Field>(a: F, b: F) -> DivisorCoords<F> {
        DivisorCoords::deg1(-a, b)
    }

    fn random_valid_divisor<F: Field, R: Rng>(
        f: &[F; 6],
        h: &[F; 3],
        cc: &CurveConstants<F>,
        rng: &mut R,
    ) -> DivisorCoords<F> {
        match rng.gen_range(0..10) {
            0 => DivisorCoords::identity(),
            1..=3 => {
                if let Some((a, b)) = find_curve_point(f, h, rng) {
                    deg1_from_point(a, b)
                } else {
                    DivisorCoords::identity()
                }
            }
            _ => {
                if let Some((a1, b1)) = find_curve_point(f, h, rng) {
                    if let Some((a2, b2)) = find_curve_point(f, h, rng) {
                        let d1 = deg1_from_point(a1, b1);
                        let d2 = deg1_from_point(a2, b2);
                        add(&d1, &d2, cc)
                    } else {
                        deg1_from_point(a1, b1)
                    }
                } else {
                    DivisorCoords::identity()
                }
            }
        }
    }

    fn random_arb_curve<F: Field, R: Rng>(rng: &mut R) -> ([F; 6], [F; 3], CurveConstants<F>) {
        let f = [
            F::random(rng), // f0
            F::random(rng), // f1
            F::random(rng), // f2
            F::random(rng), // f3
            F::random(rng), // f4
            F::one(),       // leading coeff
        ];
        let h = [
            F::random(rng), // h0
            F::random(rng), // h1
            if rng.gen_bool(0.3) {
                F::one()
            } else {
                F::zero()
            }, // h2
        ];
        let cc = CurveConstants {
            f4: f[4],
            f3: f[3],
            f2: f[2],
            f1: f[1],
            f0: f[0],
            h2: h[2],
            h1: h[1],
            h0: h[0],
        };
        (f, h, cc)
    }

    fn run_test<F: Field>(seed: u64, curves: usize, divisors_per_curve: usize) {
        let mut rng = StdRng::seed_from_u64(seed);

        for _ in 0..curves {
            let (f, h, cc): ([F; 6], [F; 3], CurveConstants<F>) = random_arb_curve(&mut rng);

            for _ in 0..divisors_per_curve {
                let d1 = random_valid_divisor(&f, &h, &cc, &mut rng);
                let d2 = random_valid_divisor(&f, &h, &cc, &mut rng);

                // Test identity
                let id = DivisorCoords::identity();
                assert_eq!(add(&d1, &id, &cc), d1, "D + 0 != D");
                assert_eq!(add(&id, &d1, &cc), d1, "0 + D != D");

                // Test double of identity
                assert!(double(&id, &cc).is_identity(), "2*0 != 0");

                // Test add with distinct divisors (like Magma tests)
                if d1 != d2 {
                    let sum = add(&d1, &d2, &cc);
                    assert!(sum.degree() <= 2, "ADD result degree > 2");

                    // Test commutativity
                    let sum_rev = add(&d2, &d1, &cc);
                    assert_eq!(sum, sum_rev, "D1 + D2 != D2 + D1");
                }

                // Test doubling
                let dbl = double(&d1, &cc);
                assert!(dbl.degree() <= 2, "DBL result degree > 2");
            }
        }
    }

    #[test]
    fn test_arb_f7() {
        run_test::<PrimeField<7>>(2001, 10, 2500);
    }

    #[test]
    fn test_arb_f31() {
        run_test::<PrimeField<31>>(2002, 10, 2500);
    }

    #[test]
    fn test_arb_f127() {
        run_test::<PrimeField<127>>(2003, 5, 2500);
    }

    #[test]
    fn test_arb_f8191() {
        run_test::<PrimeField<8191>>(2004, 3, 2500);
    }
}

// =============================================================================
// CH2 Tests - y² + h(x)y = f(x) curves over characteristic 2 fields
// =============================================================================

mod char2_blackbox {
    use super::*;
    use crate::g2::ramified::char2::{add, double, CurveConstants, DivisorCoords};

    /// Solve b² + h*b = target in characteristic 2
    /// This is the Artin-Schreier equation: b² + b = target/h² (after substitution)
    fn solve_char2_quadratic<F: Field>(h_val: F, target: F) -> Option<F> {
        if h_val.is_zero() {
            // b² = target, need square root
            // In char 2, squaring is a bijection, so we can try values
            let mut rng = rand::thread_rng();
            for _ in 0..1000 {
                let b = F::random(&mut rng);
                if b * b == target {
                    return Some(b);
                }
            }
            if target.is_zero() {
                return Some(F::zero());
            }
            return None;
        }

        // h != 0: try random values
        let mut rng = rand::thread_rng();
        for _ in 0..10000 {
            let b = F::random(&mut rng);
            // In char 2: b² + h*b = b*(b + h)
            if b * b + h_val * b == target {
                return Some(b);
            }
        }
        None
    }

    fn find_curve_point<F: Field, R: Rng>(f: &[F; 6], h: &[F; 3], rng: &mut R) -> Option<(F, F)> {
        for _ in 0..100 {
            let a = F::random(rng);

            // f(a) = a^5 + f2*a^2 + f1*a + f0 (f4=f3=0 for ch2)
            let a2 = a * a;
            let a5 = a2 * a2 * a;
            let fa = a5 + f[2] * a2 + f[1] * a + f[0];

            // h(a)
            let ha = h[2] * a2 + h[1] * a + h[0];

            if let Some(b) = solve_char2_quadratic(ha, fa) {
                return Some((a, b));
            }
        }
        None
    }

    fn deg1_from_point<F: Field>(a: F, b: F) -> DivisorCoords<F> {
        DivisorCoords::deg1(-a, b)
    }

    fn random_valid_divisor<F: Field, R: Rng>(
        f: &[F; 6],
        h: &[F; 3],
        cc: &CurveConstants<F>,
        rng: &mut R,
    ) -> DivisorCoords<F> {
        match rng.gen_range(0..10) {
            0 => DivisorCoords::identity(),
            1..=3 => {
                if let Some((a, b)) = find_curve_point(f, h, rng) {
                    deg1_from_point(a, b)
                } else {
                    DivisorCoords::identity()
                }
            }
            _ => {
                if let Some((a1, b1)) = find_curve_point(f, h, rng) {
                    if let Some((a2, b2)) = find_curve_point(f, h, rng) {
                        let d1 = deg1_from_point(a1, b1);
                        let d2 = deg1_from_point(a2, b2);
                        add(&d1, &d2, cc)
                    } else {
                        deg1_from_point(a1, b1)
                    }
                } else {
                    DivisorCoords::identity()
                }
            }
        }
    }

    fn random_ch2_curve<F: Field, R: Rng>(rng: &mut R) -> ([F; 6], [F; 3], CurveConstants<F>) {
        let f = [
            F::random(rng), // f0
            F::random(rng), // f1
            F::random(rng), // f2
            F::zero(),      // f3 = 0 for ch2
            F::zero(),      // f4 = 0 for ch2
            F::one(),       // leading coeff
        ];
        let h = [
            F::random(rng), // h0
            F::random(rng), // h1
            if rng.gen_bool(0.5) {
                F::one()
            } else {
                F::zero()
            }, // h2
        ];
        let cc = CurveConstants {
            f2: f[2],
            f1: f[1],
            f0: f[0],
            h2: h[2],
            h1: h[1],
            h0: h[0],
        };
        (f, h, cc)
    }

    fn run_test<F: Field>(seed: u64, curves: usize, divisors_per_curve: usize) {
        let mut rng = StdRng::seed_from_u64(seed);

        for _ in 0..curves {
            let (f, h, cc): ([F; 6], [F; 3], CurveConstants<F>) = random_ch2_curve(&mut rng);

            for _ in 0..divisors_per_curve {
                let d1 = random_valid_divisor(&f, &h, &cc, &mut rng);
                let d2 = random_valid_divisor(&f, &h, &cc, &mut rng);

                // Test identity
                let id = DivisorCoords::identity();
                assert_eq!(add(&d1, &id, &cc), d1, "D + 0 != D");
                assert_eq!(add(&id, &d1, &cc), d1, "0 + D != D");

                // Test double of identity
                assert!(double(&id, &cc).is_identity(), "2*0 != 0");

                // Test add with distinct divisors (like Magma tests)
                if d1 != d2 {
                    let sum = add(&d1, &d2, &cc);
                    assert!(sum.degree() <= 2, "ADD result degree > 2");

                    // Test commutativity
                    let sum_rev = add(&d2, &d1, &cc);
                    assert_eq!(sum, sum_rev, "D1 + D2 != D2 + D1");
                }

                // Test doubling
                let dbl = double(&d1, &cc);
                assert!(dbl.degree() <= 2, "DBL result degree > 2");
            }
        }
    }

    #[test]
    fn test_ch2_gf8() {
        run_test::<BinaryExtField<3>>(3001, 10, 2500);
    }

    #[test]
    fn test_ch2_gf256() {
        run_test::<BinaryExtField<8>>(3002, 10, 2500);
    }

    #[test]
    fn test_ch2_gf65536() {
        run_test::<BinaryExtField<16>>(3003, 5, 2500);
    }
}
