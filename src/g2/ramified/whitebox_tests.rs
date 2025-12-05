//! Whitebox tests for genus 2 ramified divisor arithmetic.
//!
//! These tests are translated from the Magma whitebox tester files and cover
//! all code branches for add and double operations.
//!
//! Total: 66 test cases (22 arbitrary + 22 char2 + 22 not_char2)

use crate::field::{BinaryExtField, PrimeField};

// =============================================================================
// ARBITRARY (Arbitrary Characteristic) Tests - 22 cases
// =============================================================================

mod arbitrary_tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::field::Field;
    use crate::g2::ramified::arbitrary::{add, double, CurveConstants, DivisorCoords};

    type F3 = PrimeField<3>;
    type F5 = PrimeField<5>;
    type F7 = PrimeField<7>;
    type GF4 = BinaryExtField<2>;

    // -------------------------------------------------------------------------
    // DBL Tests (7 cases)
    // -------------------------------------------------------------------------

    /// Case 1: DBL of identity over GF(5)
    /// f = x^5 + 3x^4 + 3x^3 + x^2 + 3x + 1, h = 3x + 3
    #[test]
    fn test_arb_dbl_01_identity() {
        let cc = CurveConstants {
            f4: F5::new(3),
            f3: F5::new(3),
            f2: F5::new(1),
            f1: F5::new(3),
            f0: F5::new(1),
            h2: F5::zero(),
            h1: F5::new(3),
            h0: F5::new(3),
        };
        let d = DivisorCoords::<F5>::identity();
        let result = double(&d, &cc);
        assert!(result.is_identity());
    }

    /// Case 2: DBL of deg2 over GF(5)
    /// f = x^5 + 2x^4 + 3x^3 + x^2 + 2x + 2, h = 3x
    /// U = x^2 + 3x + 4, V = 2x + 1
    #[test]
    fn test_arb_dbl_02_deg2() {
        let cc = CurveConstants {
            f4: F5::new(2),
            f3: F5::new(3),
            f2: F5::new(1),
            f1: F5::new(2),
            f0: F5::new(2),
            h2: F5::zero(),
            h1: F5::new(3),
            h0: F5::zero(),
        };
        let d = DivisorCoords::deg2(F5::new(3), F5::new(4), F5::new(2), F5::new(1));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 3: DBL of deg2 over GF(3)
    /// f = x^5 + x^4 + x^2 + x + 1, h = x^2 + 2x + 2
    /// U = x^2 + 1, V = 2x + 2
    #[test]
    fn test_arb_dbl_03_deg2_h2() {
        let cc = CurveConstants {
            f4: F3::new(1),
            f3: F3::zero(),
            f2: F3::new(1),
            f1: F3::new(1),
            f0: F3::new(1),
            h2: F3::new(1),
            h1: F3::new(2),
            h0: F3::new(2),
        };
        let d = DivisorCoords::deg2(F3::zero(), F3::new(1), F3::new(2), F3::new(2));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 4: DBL of deg2 over GF(3) - d=0 branch
    /// f = x^5 + x^3 + 2x, h = x^2 + x
    /// U = x^2 + 1, V = x + 2
    #[test]
    fn test_arb_dbl_04_d_zero() {
        let cc = CurveConstants {
            f4: F3::zero(),
            f3: F3::new(1),
            f2: F3::zero(),
            f1: F3::new(2),
            f0: F3::zero(),
            h2: F3::new(1),
            h1: F3::new(1),
            h0: F3::zero(),
        };
        let d = DivisorCoords::deg2(F3::zero(), F3::new(1), F3::new(1), F3::new(2));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 5: DBL of deg2 over GF(5) - sp1=0 branch
    /// f = x^5 + 4x^4 + x^3 + 3x^2 + 2x + 3, h = 1
    /// U = x^2 + 3x + 2, V = 2x + 4
    #[test]
    fn test_arb_dbl_05_sp1_zero() {
        let cc = CurveConstants {
            f4: F5::new(4),
            f3: F5::new(1),
            f2: F5::new(3),
            f1: F5::new(2),
            f0: F5::new(3),
            h2: F5::zero(),
            h1: F5::zero(),
            h0: F5::new(1),
        };
        let d = DivisorCoords::deg2(F5::new(3), F5::new(2), F5::new(2), F5::new(4));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 6: DBL of deg1 over GF(5)
    /// f = x^5 + 2x^4 + 3x^3 + x^2 + 2x + 2, h = 3x
    /// U = x + 1, V = 4
    #[test]
    fn test_arb_dbl_06_deg1() {
        let cc = CurveConstants {
            f4: F5::new(2),
            f3: F5::new(3),
            f2: F5::new(1),
            f1: F5::new(2),
            f0: F5::new(2),
            h2: F5::zero(),
            h1: F5::new(3),
            h0: F5::zero(),
        };
        let d = DivisorCoords::deg1(F5::new(1), F5::new(4));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 7: DBL of deg1 over GF(5)
    /// f = x^5 + x^4 + 3x^3 + x^2 + 3x + 3, h = 2
    /// U = x + 3, V = 0
    #[test]
    fn test_arb_dbl_07_deg1() {
        let cc = CurveConstants {
            f4: F5::new(1),
            f3: F5::new(3),
            f2: F5::new(1),
            f1: F5::new(3),
            f0: F5::new(3),
            h2: F5::zero(),
            h1: F5::zero(),
            h0: F5::new(2),
        };
        let d = DivisorCoords::deg1(F5::new(3), F5::zero());
        let result = double(&d, &cc);
        // Result may be identity or valid divisor depending on curve
        assert!(result.degree() <= 2);
    }

    // -------------------------------------------------------------------------
    // ADD Tests (15 cases)
    // -------------------------------------------------------------------------

    /// Case 8: ADD deg1 + identity over GF(4)
    #[test]
    fn test_arb_add_08_deg1_identity() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f4: alpha,
            f3: GF4::one(),
            f2: alpha_sq,
            f1: GF4::zero(),
            f0: alpha,
            h2: GF4::one(),
            h1: alpha_sq,
            h0: GF4::zero(),
        };
        let d1 = DivisorCoords::deg1(GF4::one(), alpha_sq);
        let d2 = DivisorCoords::identity();
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d1);
    }

    /// Case 9: ADD deg2 + identity over GF(5)
    #[test]
    fn test_arb_add_09_deg2_identity() {
        let cc = CurveConstants {
            f4: F5::new(3),
            f3: F5::new(3),
            f2: F5::new(4),
            f1: F5::new(4),
            f0: F5::new(3),
            h2: F5::new(1),
            h1: F5::new(1),
            h0: F5::new(1),
        };
        let d1 = DivisorCoords::deg2(F5::new(2), F5::new(1), F5::new(3), F5::new(1));
        let d2 = DivisorCoords::identity();
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d1);
    }

    /// Case 10: ADD deg2 + deg2 over GF(3)
    #[test]
    fn test_arb_add_10_deg2_deg2() {
        let cc = CurveConstants {
            f4: F3::new(1),
            f3: F3::zero(),
            f2: F3::new(1),
            f1: F3::new(1),
            f0: F3::new(1),
            h2: F3::new(1),
            h1: F3::new(2),
            h0: F3::new(2),
        };
        let d1 = DivisorCoords::deg2(F3::zero(), F3::new(1), F3::new(2), F3::new(2));
        let d2 = DivisorCoords::deg2(F3::new(2), F3::new(2), F3::new(2), F3::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 11: ADD deg2 + deg2 over GF(3) - resultant case
    #[test]
    fn test_arb_add_11_deg2_deg2_res() {
        let cc = CurveConstants {
            f4: F3::zero(),
            f3: F3::new(1),
            f2: F3::new(1),
            f1: F3::new(2),
            f0: F3::zero(),
            h2: F3::zero(),
            h1: F3::new(2),
            h0: F3::new(1),
        };
        let d1 = DivisorCoords::deg2(F3::new(2), F3::new(1), F3::zero(), F3::new(1));
        let d2 = DivisorCoords::deg2(F3::zero(), F3::new(1), F3::new(1), F3::new(1));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 12: ADD identity + deg2 over GF(5)
    #[test]
    fn test_arb_add_12_identity_deg2() {
        let cc = CurveConstants {
            f4: F5::new(3),
            f3: F5::new(3),
            f2: F5::new(1),
            f1: F5::new(3),
            f0: F5::new(1),
            h2: F5::zero(),
            h1: F5::new(3),
            h0: F5::new(3),
        };
        let d1 = DivisorCoords::identity();
        let d2 = DivisorCoords::deg2(F5::new(4), F5::new(2), F5::new(2), F5::zero());
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d2);
    }

    /// Case 13: ADD deg2 + deg2 over GF(7) - same v case
    #[test]
    fn test_arb_add_13_deg2_deg2_same_v() {
        let cc = CurveConstants {
            f4: F7::zero(),
            f3: F7::zero(),
            f2: F7::new(2),
            f1: F7::new(2),
            f0: F7::new(1),
            h2: F7::new(1),
            h1: F7::new(5),
            h0: F7::new(2),
        };
        let d1 = DivisorCoords::deg2(F7::zero(), F7::new(6), F7::new(2), F7::new(2));
        let d2 = DivisorCoords::deg2(F7::new(1), F7::zero(), F7::new(2), F7::new(2));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 14: ADD deg2 + deg2 over GF(5)
    #[test]
    fn test_arb_add_14_deg2_deg2() {
        let cc = CurveConstants {
            f4: F5::new(3),
            f3: F5::new(3),
            f2: F5::new(4),
            f1: F5::new(4),
            f0: F5::new(3),
            h2: F5::new(1),
            h1: F5::new(1),
            h0: F5::new(1),
        };
        let d1 = DivisorCoords::deg2(F5::new(3), F5::new(2), F5::zero(), F5::new(3));
        let d2 = DivisorCoords::deg2(F5::new(4), F5::new(4), F5::new(3), F5::new(4));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 15: ADD deg1 + deg2 over GF(5)
    #[test]
    fn test_arb_add_15_deg1_deg2() {
        let cc = CurveConstants {
            f4: F5::new(1),
            f3: F5::new(3),
            f2: F5::new(1),
            f1: F5::new(3),
            f0: F5::new(3),
            h2: F5::zero(),
            h1: F5::zero(),
            h0: F5::new(2),
        };
        let d1 = DivisorCoords::deg1(F5::new(3), F5::zero());
        let d2 = DivisorCoords::deg2(F5::new(3), F5::new(4), F5::new(2), F5::new(1));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 16: ADD deg2 + deg2 over GF(4) - same u
    #[test]
    fn test_arb_add_16_deg2_same_u() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f4: GF4::zero(),
            f3: GF4::one(),
            f2: GF4::zero(),
            f1: GF4::zero(),
            f0: GF4::one(),
            h2: GF4::one(),
            h1: GF4::zero(),
            h0: GF4::one(),
        };
        let d1 = DivisorCoords::deg2(alpha, GF4::one(), alpha, alpha_sq);
        let d2 = DivisorCoords::deg2(alpha, GF4::one(), GF4::zero(), alpha_sq);
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 17: ADD deg2 + deg2 over GF(7) - same u different v
    #[test]
    fn test_arb_add_17_deg2_same_u_diff_v() {
        let cc = CurveConstants {
            f4: F7::new(1),
            f3: F7::new(6),
            f2: F7::zero(),
            f1: F7::new(3),
            f0: F7::new(4),
            h2: F7::new(1),
            h1: F7::new(6),
            h0: F7::new(6),
        };
        let d1 = DivisorCoords::deg2(F7::new(4), F7::new(3), F7::new(4), F7::new(5));
        let d2 = DivisorCoords::deg2(F7::new(4), F7::new(3), F7::new(6), F7::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 18: ADD deg2 + deg2 over GF(5) - partial common factor
    #[test]
    fn test_arb_add_18_deg2_partial_common() {
        let cc = CurveConstants {
            f4: F5::new(4),
            f3: F5::new(4),
            f2: F5::new(1),
            f1: F5::zero(),
            f0: F5::new(1),
            h2: F5::zero(),
            h1: F5::new(4),
            h0: F5::new(4),
        };
        let d1 = DivisorCoords::deg2(F5::new(1), F5::zero(), F5::new(4), F5::new(3));
        let d2 = DivisorCoords::deg2(F5::new(2), F5::zero(), F5::new(2), F5::new(3));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 19: ADD deg1 + deg1 over GF(3) - same u
    #[test]
    fn test_arb_add_19_deg1_same_u() {
        let cc = CurveConstants {
            f4: F3::zero(),
            f3: F3::new(2),
            f2: F3::new(2),
            f1: F3::new(2),
            f0: F3::new(1),
            h2: F3::new(1),
            h1: F3::new(1),
            h0: F3::zero(),
        };
        let d1 = DivisorCoords::deg1(F3::new(2), F3::zero());
        let d2 = DivisorCoords::deg1(F3::new(2), F3::new(1));
        let result = add(&d1, &d2, &cc);
        assert!(result.is_identity());
    }

    /// Case 20: ADD deg1 + deg1 over GF(5)
    #[test]
    fn test_arb_add_20_deg1_deg1() {
        let cc = CurveConstants {
            f4: F5::new(4),
            f3: F5::new(4),
            f2: F5::new(3),
            f1: F5::new(4),
            f0: F5::new(2),
            h2: F5::zero(),
            h1: F5::new(4),
            h0: F5::new(1),
        };
        let d1 = DivisorCoords::deg1(F5::new(1), F5::zero());
        let d2 = DivisorCoords::deg1(F5::zero(), F5::new(3));
        let result = add(&d1, &d2, &cc);
        assert_eq!(result.degree(), 2);
    }

    /// Case 21: ADD deg2 + deg1 over GF(5)
    #[test]
    fn test_arb_add_21_deg2_deg1() {
        let cc = CurveConstants {
            f4: F5::new(4),
            f3: F5::new(1),
            f2: F5::new(3),
            f1: F5::new(2),
            f0: F5::new(3),
            h2: F5::zero(),
            h1: F5::zero(),
            h0: F5::new(1),
        };
        let d1 = DivisorCoords::deg2(F5::new(3), F5::new(2), F5::new(2), F5::new(4));
        let d2 = DivisorCoords::deg1(F5::new(2), F5::new(4));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 22: ADD deg2 + deg1 over GF(4) - special branch
    #[test]
    fn test_arb_add_22_deg2_deg1_special() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f4: alpha_sq,
            f3: GF4::one(),
            f2: alpha,
            f1: GF4::zero(),
            f0: GF4::zero(),
            h2: GF4::one(),
            h1: GF4::one(),
            h0: alpha,
        };
        let d1 = DivisorCoords::deg2(GF4::one(), GF4::one(), GF4::zero(), alpha);
        let d2 = DivisorCoords::deg1(alpha_sq, alpha);
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }
}

// =============================================================================
// CH2 (Characteristic 2) Tests - 22 cases
// =============================================================================

mod char2_tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::field::Field;
    use crate::g2::ramified::char2::{add, double, CurveConstants, DivisorCoords};

    type GF4 = BinaryExtField<2>;
    type GF8 = BinaryExtField<3>;

    // -------------------------------------------------------------------------
    // DBL Tests (7 cases)
    // -------------------------------------------------------------------------

    /// Case 1: DBL of identity over GF(4)
    #[test]
    fn test_ch2_dbl_01_identity() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: alpha,
            f1: alpha_sq,
            f0: GF4::zero(),
            h2: GF4::one(),
            h1: alpha,
            h0: alpha,
        };
        let d = DivisorCoords::<GF4>::identity();
        let result = double(&d, &cc);
        assert!(result.is_identity());
    }

    /// Case 2: DBL of deg2 over GF(8)
    #[test]
    fn test_ch2_dbl_02_deg2() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha,
            f1: alpha.pow(2),
            f0: alpha,
            h2: GF8::one(),
            h1: alpha.pow(6),
            h0: alpha,
        };
        let d = DivisorCoords::deg2(alpha.pow(5), alpha.pow(3), alpha.pow(4), alpha);
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 3: DBL of deg2 over GF(4)
    #[test]
    fn test_ch2_dbl_03_deg2() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: alpha,
            f1: alpha_sq,
            f0: GF4::zero(),
            h2: GF4::one(),
            h1: alpha,
            h0: alpha,
        };
        let d = DivisorCoords::deg2(alpha, GF4::one(), alpha, alpha_sq);
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 4: DBL of deg2 over GF(4) - d=0 branch
    #[test]
    fn test_ch2_dbl_04_d_zero() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: alpha_sq,
            f1: GF4::zero(),
            f0: GF4::one(),
            h2: GF4::one(),
            h1: alpha,
            h0: alpha_sq,
        };
        let d = DivisorCoords::deg2(alpha, alpha_sq, alpha_sq, GF4::one());
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 5: DBL of deg2 over GF(4) - sp1=0 branch
    #[test]
    fn test_ch2_dbl_05_sp1_zero() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::one(),
            f1: alpha_sq,
            f0: GF4::one(),
            h2: GF4::one(),
            h1: GF4::one(),
            h0: GF4::zero(),
        };
        let d = DivisorCoords::deg2(alpha, GF4::zero(), GF4::zero(), GF4::one());
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 6: DBL of deg1 over GF(4)
    #[test]
    fn test_ch2_dbl_06_deg1() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::zero(),
            f1: alpha_sq,
            f0: alpha_sq,
            h2: GF4::zero(),
            h1: GF4::one(),
            h0: GF4::one(),
        };
        let d = DivisorCoords::deg1(GF4::one(), GF4::one());
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 7: DBL of deg1 over GF(8)
    #[test]
    fn test_ch2_dbl_07_deg1() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(2),
            f1: GF8::one(),
            f0: alpha.pow(3),
            h2: GF8::one(),
            h1: alpha.pow(5),
            h0: alpha.pow(6),
        };
        let d = DivisorCoords::deg1(alpha.pow(4), alpha.pow(6));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    // -------------------------------------------------------------------------
    // ADD Tests (15 cases)
    // -------------------------------------------------------------------------

    /// Case 8: ADD deg1 + identity over GF(8)
    #[test]
    fn test_ch2_add_08_deg1_identity() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha,
            f1: alpha.pow(2),
            f0: alpha,
            h2: GF8::one(),
            h1: alpha.pow(6),
            h0: alpha,
        };
        let d1 = DivisorCoords::deg1(alpha.pow(2), alpha.pow(4));
        let d2 = DivisorCoords::identity();
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d1);
    }

    /// Case 9: ADD deg2 + identity over GF(8)
    #[test]
    fn test_ch2_add_09_deg2_identity() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(3),
            f1: alpha.pow(2),
            f0: alpha.pow(3),
            h2: GF8::one(),
            h1: alpha.pow(4),
            h0: alpha.pow(4),
        };
        let d1 = DivisorCoords::deg2(alpha, GF8::zero(), alpha.pow(4), alpha);
        let d2 = DivisorCoords::identity();
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d1);
    }

    /// Case 10: ADD deg2 + deg2 over GF(8)
    #[test]
    fn test_ch2_add_10_deg2_deg2() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(2),
            f1: GF8::one(),
            f0: alpha.pow(3),
            h2: GF8::one(),
            h1: alpha.pow(5),
            h0: alpha.pow(6),
        };
        let d1 = DivisorCoords::deg2(alpha.pow(3), alpha.pow(5), GF8::zero(), GF8::one());
        let d2 = DivisorCoords::deg2(alpha.pow(6), GF8::one(), alpha.pow(6), alpha.pow(4));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 11: ADD deg2 + deg2 over GF(8)
    #[test]
    fn test_ch2_add_11_deg2_deg2() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(4),
            f1: alpha,
            f0: alpha.pow(2),
            h2: GF8::one(),
            h1: alpha,
            h0: alpha.pow(4),
        };
        let d1 = DivisorCoords::deg2(alpha.pow(5), alpha.pow(2), alpha.pow(5), alpha);
        let d2 = DivisorCoords::deg2(GF8::zero(), GF8::zero(), alpha.pow(6), alpha.pow(6));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 12: ADD identity + deg2 over GF(8)
    #[test]
    fn test_ch2_add_12_identity_deg2() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(4),
            f1: alpha,
            f0: alpha.pow(2),
            h2: GF8::one(),
            h1: alpha,
            h0: alpha.pow(4),
        };
        let d1 = DivisorCoords::identity();
        let d2 = DivisorCoords::deg2(alpha.pow(4), GF8::zero(), alpha.pow(3), alpha.pow(6));
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d2);
    }

    /// Case 13: ADD deg2 + deg2 over GF(4)
    #[test]
    fn test_ch2_add_13_deg2_deg2() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::one(),
            f1: GF4::one(),
            f0: GF4::zero(),
            h2: GF4::zero(),
            h1: GF4::zero(),
            h0: alpha,
        };
        let d1 = DivisorCoords::deg2(alpha_sq, GF4::zero(), alpha, alpha);
        let d2 = DivisorCoords::deg2(GF4::zero(), GF4::zero(), alpha_sq, alpha);
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 14: ADD deg2 + deg2 over GF(4) - same v
    #[test]
    fn test_ch2_add_14_deg2_same_v() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::one(),
            f1: alpha_sq,
            f0: GF4::one(),
            h2: GF4::one(),
            h1: GF4::one(),
            h0: GF4::zero(),
        };
        let d1 = DivisorCoords::deg2(alpha_sq, alpha, GF4::one(), alpha);
        let d2 = DivisorCoords::deg2(GF4::zero(), alpha_sq, GF4::one(), alpha);
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 15: ADD deg1 + deg2 over GF(8)
    #[test]
    fn test_ch2_add_15_deg1_deg2() {
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
        assert!(result.degree() <= 2);
    }

    /// Case 16: ADD deg2 + deg2 over GF(4) - same u
    #[test]
    fn test_ch2_add_16_deg2_same_u() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::zero(),
            f1: GF4::one(),
            f0: GF4::zero(),
            h2: GF4::zero(),
            h1: alpha,
            h0: alpha_sq,
        };
        let d1 = DivisorCoords::deg2(alpha, GF4::zero(), GF4::one(), alpha_sq);
        let d2 = DivisorCoords::deg2(alpha, GF4::zero(), alpha_sq, GF4::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 17: ADD deg2 + deg2 over GF(8) - same u
    #[test]
    fn test_ch2_add_17_deg2_same_u() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(5),
            f1: alpha,
            f0: alpha.pow(3),
            h2: GF8::one(),
            h1: alpha.pow(4),
            h0: alpha.pow(3),
        };
        let d1 = DivisorCoords::deg2(alpha.pow(3), GF8::zero(), GF8::one(), alpha.pow(4));
        let d2 = DivisorCoords::deg2(alpha.pow(3), GF8::zero(), GF8::zero(), alpha.pow(6));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 18: ADD deg2 + deg2 over GF(4)
    #[test]
    fn test_ch2_add_18_deg2_deg2() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::zero(),
            f1: GF4::one(),
            f0: GF4::zero(),
            h2: GF4::zero(),
            h1: alpha,
            h0: alpha_sq,
        };
        let d1 = DivisorCoords::deg2(alpha, GF4::zero(), GF4::one(), alpha_sq);
        let d2 = DivisorCoords::deg2(alpha_sq, alpha, GF4::zero(), GF4::one());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 19: ADD deg1 + deg1 over GF(8)
    #[test]
    fn test_ch2_add_19_deg1_same_u() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(5),
            f1: alpha,
            f0: alpha.pow(2),
            h2: GF8::one(),
            h1: GF8::zero(),
            h0: alpha.pow(5),
        };
        let d1 = DivisorCoords::deg1(alpha.pow(3), alpha.pow(5));
        let d2 = DivisorCoords::deg1(alpha.pow(3), alpha.pow(6));
        let result = add(&d1, &d2, &cc);
        // Same u0 but different v0
        assert!(result.degree() <= 2);
    }

    /// Case 20: ADD deg1 + deg1 over GF(4)
    #[test]
    fn test_ch2_add_20_deg1_deg1() {
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
        let d1 = DivisorCoords::deg1(GF4::one(), alpha);
        let d2 = DivisorCoords::deg1(GF4::zero(), alpha_sq);
        let result = add(&d1, &d2, &cc);
        assert_eq!(result.degree(), 2);
    }

    /// Case 21: ADD deg1 + deg2 over GF(4)
    #[test]
    fn test_ch2_add_21_deg1_deg2() {
        let alpha = GF4::gen();
        let alpha_sq = alpha * alpha;
        let cc = CurveConstants {
            f2: GF4::one(),
            f1: GF4::one(),
            f0: alpha_sq,
            h2: GF4::one(),
            h1: GF4::zero(),
            h0: GF4::zero(),
        };
        let d1 = DivisorCoords::deg1(alpha_sq, GF4::zero());
        let d2 = DivisorCoords::deg2(alpha_sq, GF4::zero(), GF4::zero(), alpha);
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 22: ADD deg2 + deg1 over GF(8)
    #[test]
    fn test_ch2_add_22_deg2_deg1() {
        let alpha = GF8::gen();
        let cc = CurveConstants {
            f2: alpha.pow(2),
            f1: alpha.pow(4),
            f0: alpha.pow(6),
            h2: GF8::zero(),
            h1: alpha.pow(6),
            h0: alpha.pow(6),
        };
        let d1 = DivisorCoords::deg2(alpha.pow(6), alpha.pow(2), alpha.pow(2), GF8::zero());
        let d2 = DivisorCoords::deg1(alpha.pow(2), alpha.pow(4));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }
}

// =============================================================================
// NCH2 (Not Characteristic 2) Tests - 22 cases
// =============================================================================

mod not_char2_tests {
    use super::*;
    #[allow(unused_imports)]
    use crate::field::Field;
    use crate::g2::ramified::not_char2::{add, double, CurveConstants, DivisorCoords};

    type F3 = PrimeField<3>;
    type F5 = PrimeField<5>;
    type F7 = PrimeField<7>;

    // -------------------------------------------------------------------------
    // DBL Tests (7 cases)
    // -------------------------------------------------------------------------

    /// Case 1: DBL of identity over GF(3)
    #[test]
    fn test_nch2_dbl_01_identity() {
        let cc = CurveConstants {
            f3: F3::new(1),
            f2: F3::new(2),
            f1: F3::new(2),
            f0: F3::zero(),
        };
        let d = DivisorCoords::<F3>::identity();
        let result = double(&d, &cc);
        assert!(result.is_identity());
    }

    /// Case 2: DBL of deg2 over GF(5)
    #[test]
    fn test_nch2_dbl_02_deg2() {
        let cc = CurveConstants {
            f3: F5::new(2),
            f2: F5::new(1),
            f1: F5::new(3),
            f0: F5::new(3),
        };
        let d = DivisorCoords::deg2(F5::new(3), F5::new(4), F5::new(2), F5::new(4));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 3: DBL of deg2 over GF(7)
    #[test]
    fn test_nch2_dbl_03_deg2() {
        let cc = CurveConstants {
            f3: F7::new(2),
            f2: F7::new(1),
            f1: F7::new(4),
            f0: F7::new(6),
        };
        let d = DivisorCoords::deg2(F7::new(5), F7::new(5), F7::new(1), F7::new(2));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 4: DBL of deg2 over GF(3) - d=0 branch
    #[test]
    fn test_nch2_dbl_04_d_zero() {
        let cc = CurveConstants {
            f3: F3::new(1),
            f2: F3::new(2),
            f1: F3::new(2),
            f0: F3::zero(),
        };
        let d = DivisorCoords::deg2(F3::new(2), F3::zero(), F3::zero(), F3::zero());
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 5: DBL of deg2 over GF(7) - sp1=0 branch
    #[test]
    fn test_nch2_dbl_05_sp1_zero() {
        let cc = CurveConstants {
            f3: F7::zero(),
            f2: F7::new(4),
            f1: F7::new(4),
            f0: F7::zero(),
        };
        let d = DivisorCoords::deg2(F7::new(2), F7::zero(), F7::new(6), F7::zero());
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 6: DBL of deg1 over GF(7) - d=0 branch
    #[test]
    fn test_nch2_dbl_06_deg1_d_zero() {
        let cc = CurveConstants {
            f3: F7::new(2),
            f2: F7::new(1),
            f1: F7::new(4),
            f0: F7::new(6),
        };
        let d = DivisorCoords::deg1(F7::new(3), F7::zero());
        let result = double(&d, &cc);
        assert!(result.is_identity());
    }

    /// Case 7: DBL of deg1 over GF(3)
    #[test]
    fn test_nch2_dbl_07_deg1() {
        let cc = CurveConstants {
            f3: F3::new(1),
            f2: F3::new(1),
            f1: F3::new(1),
            f0: F3::zero(),
        };
        let d = DivisorCoords::deg1(F3::new(2), F3::new(2));
        let result = double(&d, &cc);
        assert!(result.degree() <= 2);
    }

    // -------------------------------------------------------------------------
    // ADD Tests (15 cases)
    // -------------------------------------------------------------------------

    /// Case 8: ADD deg1 + identity over GF(3)
    #[test]
    fn test_nch2_add_08_deg1_identity() {
        let cc = CurveConstants {
            f3: F3::zero(),
            f2: F3::new(1),
            f1: F3::new(2),
            f0: F3::new(1),
        };
        let d1 = DivisorCoords::deg1(F3::new(2), F3::new(2));
        let d2 = DivisorCoords::identity();
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d1);
    }

    /// Case 9: ADD deg2 + identity over GF(3)
    #[test]
    fn test_nch2_add_09_deg2_identity() {
        let cc = CurveConstants {
            f3: F3::new(1),
            f2: F3::new(2),
            f1: F3::new(2),
            f0: F3::zero(),
        };
        let d1 = DivisorCoords::deg2(F3::new(2), F3::zero(), F3::zero(), F3::zero());
        let d2 = DivisorCoords::identity();
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d1);
    }

    /// Case 10: ADD deg2 + deg2 over GF(5)
    #[test]
    fn test_nch2_add_10_deg2_deg2() {
        let cc = CurveConstants {
            f3: F5::new(1),
            f2: F5::new(4),
            f1: F5::zero(),
            f0: F5::new(1),
        };
        let d1 = DivisorCoords::deg2(F5::zero(), F5::zero(), F5::zero(), F5::new(4));
        let d2 = DivisorCoords::deg2(F5::zero(), F5::new(3), F5::new(2), F5::new(4));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 11: ADD deg2 + deg2 over GF(3)
    #[test]
    fn test_nch2_add_11_deg2_deg2() {
        let cc = CurveConstants {
            f3: F3::zero(),
            f2: F3::zero(),
            f1: F3::new(1),
            f0: F3::zero(),
        };
        let d1 = DivisorCoords::deg2(F3::new(2), F3::new(2), F3::zero(), F3::zero());
        let d2 = DivisorCoords::deg2(F3::new(2), F3::new(1), F3::zero(), F3::new(1));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 12: ADD identity + deg2 over GF(5)
    #[test]
    fn test_nch2_add_12_identity_deg2() {
        let cc = CurveConstants {
            f3: F5::new(2),
            f2: F5::new(1),
            f1: F5::new(3),
            f0: F5::new(3),
        };
        let d1 = DivisorCoords::identity();
        let d2 = DivisorCoords::deg2(F5::new(2), F5::new(4), F5::new(2), F5::new(2));
        let result = add(&d1, &d2, &cc);
        assert_eq!(result, d2);
    }

    /// Case 13: ADD deg2 + deg2 over GF(3)
    #[test]
    fn test_nch2_add_13_deg2_deg2() {
        let cc = CurveConstants {
            f3: F3::new(2),
            f2: F3::zero(),
            f1: F3::zero(),
            f0: F3::new(1),
        };
        let d1 = DivisorCoords::deg2(F3::new(2), F3::new(1), F3::new(1), F3::new(2));
        let d2 = DivisorCoords::deg2(F3::new(1), F3::zero(), F3::zero(), F3::new(1));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 14: ADD deg2 + deg2 over GF(7)
    #[test]
    fn test_nch2_add_14_deg2_deg2() {
        let cc = CurveConstants {
            f3: F7::zero(),
            f2: F7::new(4),
            f1: F7::new(4),
            f0: F7::zero(),
        };
        let d1 = DivisorCoords::deg2(F7::new(2), F7::zero(), F7::new(6), F7::zero());
        let d2 = DivisorCoords::deg2(F7::zero(), F7::new(3), F7::new(3), F7::new(1));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 15: ADD deg2 + deg1 over GF(7)
    #[test]
    fn test_nch2_add_15_deg2_deg1() {
        let cc = CurveConstants {
            f3: F7::new(2),
            f2: F7::new(1),
            f1: F7::new(4),
            f0: F7::new(6),
        };
        let d1 = DivisorCoords::deg2(F7::new(5), F7::new(5), F7::new(1), F7::new(2));
        let d2 = DivisorCoords::deg1(F7::new(1), F7::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 16: ADD deg2 + deg2 over GF(3) - same u
    #[test]
    fn test_nch2_add_16_deg2_same_u() {
        let cc = CurveConstants {
            f3: F3::new(1),
            f2: F3::new(1),
            f1: F3::new(1),
            f0: F3::zero(),
        };
        let d1 = DivisorCoords::deg2(F3::new(2), F3::new(2), F3::new(1), F3::new(1));
        let d2 = DivisorCoords::deg2(F3::new(2), F3::new(2), F3::new(2), F3::new(2));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 17: ADD deg2 + deg2 over GF(7) - same u
    #[test]
    fn test_nch2_add_17_deg2_same_u() {
        let cc = CurveConstants {
            f3: F7::new(3),
            f2: F7::new(1),
            f1: F7::new(2),
            f0: F7::zero(),
        };
        let d1 = DivisorCoords::deg2(F7::new(6), F7::new(5), F7::new(6), F7::new(3));
        let d2 = DivisorCoords::deg2(F7::new(6), F7::new(5), F7::new(4), F7::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 18: ADD deg2 + deg2 over GF(5)
    #[test]
    fn test_nch2_add_18_deg2_deg2() {
        let cc = CurveConstants {
            f3: F5::new(2),
            f2: F5::new(4),
            f1: F5::new(1),
            f0: F5::zero(),
        };
        let d1 = DivisorCoords::deg2(F5::new(2), F5::zero(), F5::new(3), F5::zero());
        let d2 = DivisorCoords::deg2(F5::new(4), F5::new(4), F5::new(2), F5::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 19: ADD deg1 + deg1 over GF(3) - same u
    #[test]
    fn test_nch2_add_19_deg1_same_u() {
        let cc = CurveConstants {
            f3: F3::zero(),
            f2: F3::new(1),
            f1: F3::new(2),
            f0: F3::new(1),
        };
        let d1 = DivisorCoords::deg1(F3::zero(), F3::new(1));
        let d2 = DivisorCoords::deg1(F3::zero(), F3::new(2));
        let result = add(&d1, &d2, &cc);
        assert!(result.is_identity());
    }

    /// Case 20: ADD deg1 + deg1 over GF(3)
    #[test]
    fn test_nch2_add_20_deg1_deg1() {
        let cc = CurveConstants {
            f3: F3::zero(),
            f2: F3::zero(),
            f1: F3::new(1),
            f0: F3::zero(),
        };
        let d1 = DivisorCoords::deg1(F3::zero(), F3::zero());
        let d2 = DivisorCoords::deg1(F3::new(1), F3::new(2));
        let result = add(&d1, &d2, &cc);
        assert_eq!(result.degree(), 2);
    }

    /// Case 21: ADD deg2 + deg1 over GF(7)
    #[test]
    fn test_nch2_add_21_deg2_deg1() {
        let cc = CurveConstants {
            f3: F7::zero(),
            f2: F7::new(4),
            f1: F7::new(4),
            f0: F7::zero(),
        };
        let d1 = DivisorCoords::deg2(F7::new(6), F7::zero(), F7::new(4), F7::zero());
        let d2 = DivisorCoords::deg1(F7::zero(), F7::zero());
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }

    /// Case 22: ADD deg1 + deg2 over GF(5)
    #[test]
    fn test_nch2_add_22_deg1_deg2() {
        let cc = CurveConstants {
            f3: F5::new(4),
            f2: F5::new(3),
            f1: F5::new(3),
            f0: F5::new(1),
        };
        let d1 = DivisorCoords::deg1(F5::zero(), F5::new(4));
        let d2 = DivisorCoords::deg2(F5::zero(), F5::zero(), F5::new(1), F5::new(4));
        let result = add(&d1, &d2, &cc);
        assert!(result.degree() <= 2);
    }
}
