//! Integration tests comparing g2 specialized formulas against generic polynomial arithmetic.
//!
//! These tests verify that the explicit coordinate formulas in the ramified module
//! produce identical results to the generic polynomial-level algorithms.

use crate::field::{Field, PrimeField};
use crate::g2::ramified::{arbitrary, not_char2};
use crate::generic::ramified as generic_ram;
use crate::poly::Poly;

type F7 = PrimeField<7>;
type F11 = PrimeField<11>;
type F23 = PrimeField<23>;

/// Helper to create a polynomial from coefficients
fn poly<const P: u64>(coeffs: &[u64]) -> Poly<PrimeField<P>> {
    Poly::from_coeffs(coeffs.iter().map(|&c| PrimeField::<P>::new(c)).collect())
}

/// Convert specialized g2 ramified divisor to generic divisor format
fn arbitrary_to_generic<F: crate::field::Field>(
    d: &arbitrary::DivisorCoords<F>,
    f: &Poly<F>,
) -> generic_ram::Divisor<F> {
    let u = match d.degree() {
        0 => Poly::constant(F::one()),
        1 => Poly::from_coeffs(vec![d.u0, F::one()]),
        2 => Poly::from_coeffs(vec![d.u0, d.u1, F::one()]),
        _ => unreachable!(),
    };

    let v = match d.degree() {
        0 => Poly::zero(),
        1 => Poly::constant(d.v0),
        2 => Poly::from_coeffs(vec![d.v0, d.v1]),
        _ => unreachable!(),
    };

    // For generic ramified (y² = f(x), h = 0): w = (f - v²) / u
    let v_sq = &v * &v;
    let w = (f - &v_sq).exact_div(&u);

    generic_ram::Divisor::new(u, v, w)
}

/// Convert generic divisor to specialized g2 ramified format
fn generic_to_arbitrary<F: crate::field::Field>(
    d: &generic_ram::Divisor<F>,
) -> arbitrary::DivisorCoords<F> {
    let deg = d.u.deg();
    match deg {
        0 | -1 => arbitrary::DivisorCoords::identity(),
        1 => {
            let u0 = d.u.coeff(0);
            let v0 = d.v.coeff(0);
            arbitrary::DivisorCoords::deg1(u0, v0)
        }
        2 => {
            let u1 = d.u.coeff(1);
            let u0 = d.u.coeff(0);
            let v1 = d.v.coeff(1);
            let v0 = d.v.coeff(0);
            arbitrary::DivisorCoords::deg2(u1, u0, v1, v0)
        }
        _ => panic!("Expected degree <= 2 for g2"),
    }
}

/// Convert nch2 divisor to generic format (h = 0 case)
fn not_char2_to_generic<F: crate::field::Field>(
    d: &not_char2::DivisorCoords<F>,
    f: &Poly<F>,
) -> generic_ram::Divisor<F> {
    let u = match d.degree() {
        0 => Poly::constant(F::one()),
        1 => Poly::from_coeffs(vec![d.u0, F::one()]),
        2 => Poly::from_coeffs(vec![d.u0, d.u1, F::one()]),
        _ => unreachable!(),
    };

    let v = match d.degree() {
        0 => Poly::zero(),
        1 => Poly::constant(d.v0),
        2 => Poly::from_coeffs(vec![d.v0, d.v1]),
        _ => unreachable!(),
    };

    // w = (f - v²) / u
    let v_sq = &v * &v;
    let w = (f - &v_sq).exact_div(&u);

    generic_ram::Divisor::new(u, v, w)
}

/// Convert generic divisor to nch2 format
fn generic_to_not_char2<F: crate::field::Field>(
    d: &generic_ram::Divisor<F>,
) -> not_char2::DivisorCoords<F> {
    let deg = d.u.deg();
    match deg {
        0 | -1 => not_char2::DivisorCoords::identity(),
        1 => {
            let u0 = d.u.coeff(0);
            let v0 = d.v.coeff(0);
            not_char2::DivisorCoords::deg1(u0, v0)
        }
        2 => {
            let u1 = d.u.coeff(1);
            let u0 = d.u.coeff(0);
            let v1 = d.v.coeff(1);
            let v0 = d.v.coeff(0);
            not_char2::DivisorCoords::deg2(u1, u0, v1, v0)
        }
        _ => panic!("Expected degree <= 2 for g2"),
    }
}

/// Find a valid degree-1 divisor on a ramified g2 curve y² = f(x)
fn find_deg1_divisor<const P: u64>(
    f: &Poly<PrimeField<P>>,
) -> Option<generic_ram::Divisor<PrimeField<P>>> {
    // Try to find x where f(x) is a square
    for a_val in 0..P {
        let a = PrimeField::<P>::new(a_val);
        let fa = f.eval(a);

        // Check if fa is a quadratic residue using Euler's criterion
        // a^((p-1)/2) = 1 means a is a square (or a = 0)
        if fa.is_zero() {
            // v = 0, u = x - a
            let neg_a = -a;
            let u = Poly::from_coeffs(vec![neg_a, PrimeField::<P>::one()]);
            let v = Poly::zero();
            let w = f.exact_div(&u);
            return Some(generic_ram::Divisor::new(u, v, w));
        }

        let exp = (P - 1) / 2;
        if fa.pow(exp).is_one() {
            // fa is a square, find its root
            // For p ≡ 3 (mod 4): sqrt(a) = a^((p+1)/4)
            // For p ≡ 1 (mod 4): use Tonelli-Shanks (simplified brute force here)
            for b_val in 0..P {
                let b = PrimeField::<P>::new(b_val);
                if b * b == fa {
                    let neg_a = -a;
                    let u = Poly::from_coeffs(vec![neg_a, PrimeField::<P>::one()]);
                    let v = Poly::constant(b);
                    let v_sq = &v * &v;
                    let diff = f - &v_sq;

                    // Verify divisibility
                    let (_, rem) = diff.div_rem(&u);
                    if rem.is_zero() {
                        let w = diff.exact_div(&u);
                        return Some(generic_ram::Divisor::new(u, v, w));
                    }
                }
            }
        }
    }
    None
}

// ============================================================================
// Tests for arbitrary characteristic (arbitrary) module
// ============================================================================

mod arbitrary_tests {
    use super::*;
    use crate::field::Field;

    /// Create a test curve for arb: y² + h(x)y = f(x)
    /// For testing with generic (which uses h=0), we set h=0 here too
    fn make_arbitrary_curve_f7() -> (arbitrary::CurveConstants<F7>, Poly<F7>) {
        // f(x) = x^5 + 3x^4 + 2x^3 + x^2 + 4x + 1
        // h(x) = 0 (to match generic ramified)
        let cc = arbitrary::CurveConstants {
            f4: F7::new(3),
            f3: F7::new(2),
            f2: F7::new(1),
            f1: F7::new(4),
            f0: F7::new(1),
            h2: F7::zero(),
            h1: F7::zero(),
            h0: F7::zero(),
        };
        let f = poly::<7>(&[1, 4, 1, 2, 3, 1]);
        (cc, f)
    }

    #[test]
    fn test_arbitrary_add_vs_generic_add() {
        let (cc, f) = make_arbitrary_curve_f7();
        let g = 2;

        // Find valid divisors
        let d1_gen = find_deg1_divisor::<7>(&f).expect("Should find divisor");
        let d2_gen = generic_ram::double(&d1_gen, &f, g);

        // Convert to specialized format
        let d1_arb = generic_to_arbitrary(&d1_gen);
        let d2_arb = generic_to_arbitrary(&d2_gen);

        // Compute via both methods
        let result_gen = generic_ram::add(&d1_gen, &d2_gen, &f, g);
        let result_arb = arbitrary::add(&d1_arb, &d2_arb, &cc);

        // Convert arb result to generic for comparison
        let result_arb_gen = arbitrary_to_generic(&result_arb, &f);

        // Compare u and v polynomials
        assert_eq!(result_gen.u, result_arb_gen.u, "u polynomials differ");
        assert_eq!(result_gen.v, result_arb_gen.v, "v polynomials differ");
    }

    #[test]
    fn test_arbitrary_double_vs_generic_double() {
        let (cc, f) = make_arbitrary_curve_f7();
        let g = 2;

        let d_gen = find_deg1_divisor::<7>(&f).expect("Should find divisor");
        let d_arb = generic_to_arbitrary(&d_gen);

        // Double via both methods
        let result_gen = generic_ram::double(&d_gen, &f, g);
        let result_arb = arbitrary::double(&d_arb, &cc);

        let result_arb_gen = arbitrary_to_generic(&result_arb, &f);

        assert_eq!(result_gen.u, result_arb_gen.u, "u polynomials differ");
        assert_eq!(result_gen.v, result_arb_gen.v, "v polynomials differ");
    }

    #[test]
    fn test_arbitrary_add_associativity() {
        let (cc, f) = make_arbitrary_curve_f7();
        let _g = 2;

        let d1 = find_deg1_divisor::<7>(&f).expect("Should find divisor");
        let d1_arb = generic_to_arbitrary(&d1);
        let d2_arb = arbitrary::double(&d1_arb, &cc);
        let d3_arb = arbitrary::double(&d2_arb, &cc);

        // (D1 + D2) + D3
        let sum12 = arbitrary::add(&d1_arb, &d2_arb, &cc);
        let left = arbitrary::add(&sum12, &d3_arb, &cc);

        // D1 + (D2 + D3)
        let sum23 = arbitrary::add(&d2_arb, &d3_arb, &cc);
        let right = arbitrary::add(&d1_arb, &sum23, &cc);

        let left_gen = arbitrary_to_generic(&left, &f);
        let right_gen = arbitrary_to_generic(&right, &f);

        assert_eq!(left_gen.u, right_gen.u, "Associativity failed for u");
        assert_eq!(left_gen.v, right_gen.v, "Associativity failed for v");
    }

    #[test]
    fn test_arbitrary_nucomp_vs_specialized_add() {
        let (cc, f) = make_arbitrary_curve_f7();
        let g = 2;

        let d1_gen = find_deg1_divisor::<7>(&f).expect("Should find divisor");
        let d2_gen = generic_ram::double(&d1_gen, &f, g);

        let d1_arb = generic_to_arbitrary(&d1_gen);
        let d2_arb = generic_to_arbitrary(&d2_gen);

        // Use nucomp (generic) vs specialized add (arb)
        let result_nucomp = generic_ram::nucomp(&d1_gen, &d2_gen, &f, g);
        let result_arb = arbitrary::add(&d1_arb, &d2_arb, &cc);
        let result_arb_gen = arbitrary_to_generic(&result_arb, &f);

        assert_eq!(result_nucomp.u, result_arb_gen.u, "u polynomials differ");
        assert_eq!(result_nucomp.v, result_arb_gen.v, "v polynomials differ");
    }
}

// ============================================================================
// Tests for not characteristic 2 (not_char2) module
// ============================================================================

mod not_char2_tests {
    use super::*;

    /// Create a test curve for nch2: y² = f(x) (h = 0, f4 = 0)
    fn make_not_char2_curve_f11() -> (not_char2::CurveConstants<F11>, Poly<F11>) {
        // f(x) = x^5 + 2x^3 + 3x^2 + 4x + 5 (f4 = 0 for nch2)
        let cc = not_char2::CurveConstants {
            f3: F11::new(2),
            f2: F11::new(3),
            f1: F11::new(4),
            f0: F11::new(5),
        };
        // For generic, f must have degree 5 with leading coeff 1
        let f = poly::<11>(&[5, 4, 3, 2, 0, 1]);
        (cc, f)
    }

    #[test]
    fn test_not_char2_add_vs_generic_add() {
        let (cc, f) = make_not_char2_curve_f11();
        let g = 2;

        let d1_gen = find_deg1_divisor::<11>(&f).expect("Should find divisor");
        let d2_gen = generic_ram::double(&d1_gen, &f, g);

        let d1_nch2 = generic_to_not_char2(&d1_gen);
        let d2_nch2 = generic_to_not_char2(&d2_gen);

        let result_gen = generic_ram::add(&d1_gen, &d2_gen, &f, g);
        let result_nch2 = not_char2::add(&d1_nch2, &d2_nch2, &cc);
        let result_nch2_gen = not_char2_to_generic(&result_nch2, &f);

        assert_eq!(result_gen.u, result_nch2_gen.u, "u polynomials differ");
        assert_eq!(result_gen.v, result_nch2_gen.v, "v polynomials differ");
    }

    #[test]
    fn test_not_char2_double_vs_generic_double() {
        let (cc, f) = make_not_char2_curve_f11();
        let g = 2;

        let d_gen = find_deg1_divisor::<11>(&f).expect("Should find divisor");
        let d_nch2 = generic_to_not_char2(&d_gen);

        let result_gen = generic_ram::double(&d_gen, &f, g);
        let result_nch2 = not_char2::double(&d_nch2, &cc);
        let result_nch2_gen = not_char2_to_generic(&result_nch2, &f);

        assert_eq!(result_gen.u, result_nch2_gen.u, "u polynomials differ");
        assert_eq!(result_gen.v, result_nch2_gen.v, "v polynomials differ");
    }

    #[test]
    fn test_not_char2_add_commutativity() {
        let (cc, f) = make_not_char2_curve_f11();
        let g = 2;

        let d1_gen = find_deg1_divisor::<11>(&f).expect("Should find divisor");
        let d2_gen = generic_ram::double(&d1_gen, &f, g);

        let d1_nch2 = generic_to_not_char2(&d1_gen);
        let d2_nch2 = generic_to_not_char2(&d2_gen);

        let sum1 = not_char2::add(&d1_nch2, &d2_nch2, &cc);
        let sum2 = not_char2::add(&d2_nch2, &d1_nch2, &cc);

        let sum1_gen = not_char2_to_generic(&sum1, &f);
        let sum2_gen = not_char2_to_generic(&sum2, &f);

        assert_eq!(sum1_gen.u, sum2_gen.u, "Commutativity failed for u");
        assert_eq!(sum1_gen.v, sum2_gen.v, "Commutativity failed for v");
    }

    #[test]
    fn test_not_char2_nuduple_vs_specialized_double() {
        let (cc, f) = make_not_char2_curve_f11();
        let g = 2;

        let d_gen = find_deg1_divisor::<11>(&f).expect("Should find divisor");
        let d_nch2 = generic_to_not_char2(&d_gen);

        // Use nuduple (generic) vs specialized double (nch2)
        let result_nuduple = generic_ram::nuduple(&d_gen, &f, g);
        let result_nch2 = not_char2::double(&d_nch2, &cc);
        let result_nch2_gen = not_char2_to_generic(&result_nch2, &f);

        assert_eq!(result_nuduple.u, result_nch2_gen.u, "u polynomials differ");
        assert_eq!(result_nuduple.v, result_nch2_gen.v, "v polynomials differ");
    }
}

// ============================================================================
// Tests across larger primes
// ============================================================================

mod larger_prime_tests {
    use super::*;

    fn make_curve_f23() -> (not_char2::CurveConstants<F23>, Poly<F23>) {
        let cc = not_char2::CurveConstants {
            f3: F23::new(5),
            f2: F23::new(7),
            f1: F23::new(11),
            f0: F23::new(13),
        };
        let f = poly::<23>(&[13, 11, 7, 5, 0, 1]);
        (cc, f)
    }

    #[test]
    fn test_not_char2_chain_of_additions_f23() {
        let (cc, f) = make_curve_f23();
        let g = 2;

        let d1_gen = find_deg1_divisor::<23>(&f).expect("Should find divisor");
        let d1_nch2 = generic_to_not_char2(&d1_gen);

        // Build a chain: D, 2D, 3D, 4D, 5D using specialized formulas
        let d2 = not_char2::double(&d1_nch2, &cc);
        let d3 = not_char2::add(&d2, &d1_nch2, &cc);
        let d4 = not_char2::double(&d2, &cc);
        let d5 = not_char2::add(&d4, &d1_nch2, &cc);

        // Build the same chain using generic formulas
        let d2_gen = generic_ram::double(&d1_gen, &f, g);
        let d3_gen = generic_ram::add(&d2_gen, &d1_gen, &f, g);
        let d4_gen = generic_ram::double(&d2_gen, &f, g);
        let d5_gen = generic_ram::add(&d4_gen, &d1_gen, &f, g);

        // Compare final results
        let d5_from_nch2 = not_char2_to_generic(&d5, &f);
        assert_eq!(d5_gen.u, d5_from_nch2.u, "5D u polynomials differ");
        assert_eq!(d5_gen.v, d5_from_nch2.v, "5D v polynomials differ");

        // Also compare intermediate results
        let d3_from_nch2 = not_char2_to_generic(&d3, &f);
        assert_eq!(d3_gen.u, d3_from_nch2.u, "3D u polynomials differ");
        assert_eq!(d3_gen.v, d3_from_nch2.v, "3D v polynomials differ");
    }
}
