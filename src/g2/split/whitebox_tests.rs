//! Whitebox tests for genus 2 split (nch2, negative basis).
//!
//! The input vectors in [`super::wb_vectors`] are extracted from the Magma
//! `nch2_splitG2_whiteBox_tester.mag` (GF(5) and GF(7) blocks; GF(9) blocks are
//! skipped — the crate's prime fields can't represent them). Each case feeds an
//! in-scope divisor (pair) chosen to cover a specific formula branch, and the
//! explicit result is checked against the generic Cantor oracle.

use crate::field::PrimeField;
use crate::generic::split::{self, Divisor};
use crate::poly::Poly;

use super::not_char2::{add_neg, add_pos, double_neg, double_pos, precompute, CurveConstants, G};
use super::test_support::from_generic;
use super::wb_vectors::{Case, CASES_GF5, CASES_GF7, CASES_POS_GF5, CASES_POS_GF7};

fn build_divisor<const P: u64>(
    cc: &CurveConstants<PrimeField<P>>,
    u: &[i64],
    v: &[i64],
    n: i32,
) -> Divisor<PrimeField<P>> {
    let fe = |x: i64| PrimeField::<P>::from_i64(x);
    let up = Poly::from_coeffs(u.iter().map(|&c| fe(c)).collect());
    let vp = Poly::from_coeffs(v.iter().map(|&c| fe(c)).collect());
    let f = cc.f_poly();
    let w = (&f - &(&vp * &vp)).exact_div(&up);
    Divisor::new(up, vp, w, n)
}

fn run_cases<const P: u64>(cases: &[Case]) {
    for (i, c) in cases.iter().enumerate() {
        let cc = precompute::<PrimeField<P>>(
            PrimeField::from_i64(c.f[0]),
            PrimeField::from_i64(c.f[1]),
            PrimeField::from_i64(c.f[2]),
            PrimeField::from_i64(c.f[3]),
            PrimeField::from_i64(c.f[4]),
        );
        let f = cc.f_poly();
        let vn = cc.vn();
        let d1 = build_divisor::<P>(&cc, c.u1, c.v1, c.n1);

        if c.op == b'D' {
            let got = double_neg(&from_generic(&d1), &cc);
            let expected = from_generic(&split::double_neg(
                &d1,
                &f,
                &crate::poly::Poly::zero(),
                &vn,
                G,
            ));
            assert_eq!(got, expected, "GF({P}) whitebox DBL case #{i}");
        } else {
            let d2 = build_divisor::<P>(&cc, c.u2, c.v2, c.n2);
            let got = add_neg(&from_generic(&d1), &from_generic(&d2), &cc);
            let expected = from_generic(&split::add_neg(
                &d1,
                &d2,
                &f,
                &crate::poly::Poly::zero(),
                &vn,
                G,
            ));
            assert_eq!(got, expected, "GF({P}) whitebox ADD case #{i}");
        }
    }
}

/// Positive-basis whitebox runner (uses Vpl, add_pos/double_pos).
fn run_cases_pos<const P: u64>(cases: &[Case]) {
    for (i, c) in cases.iter().enumerate() {
        let cc = precompute::<PrimeField<P>>(
            PrimeField::from_i64(c.f[0]),
            PrimeField::from_i64(c.f[1]),
            PrimeField::from_i64(c.f[2]),
            PrimeField::from_i64(c.f[3]),
            PrimeField::from_i64(c.f[4]),
        );
        let f = cc.f_poly();
        let vpl = cc.vpl();
        let d1 = build_divisor::<P>(&cc, c.u1, c.v1, c.n1);
        if c.op == b'D' {
            let got = double_pos(&from_generic(&d1), &cc);
            let expected = from_generic(&split::double_pos(
                &d1,
                &f,
                &crate::poly::Poly::zero(),
                &vpl,
                G,
            ));
            assert_eq!(got, expected, "GF({P}) pos whitebox DBL case #{i}");
        } else {
            let d2 = build_divisor::<P>(&cc, c.u2, c.v2, c.n2);
            let got = add_pos(&from_generic(&d1), &from_generic(&d2), &cc);
            let expected = from_generic(&split::add_pos(
                &d1,
                &d2,
                &f,
                &crate::poly::Poly::zero(),
                &vpl,
                G,
            ));
            assert_eq!(got, expected, "GF({P}) pos whitebox ADD case #{i}");
        }
    }
}

#[test]
fn whitebox_gf5() {
    run_cases::<5>(CASES_GF5);
}

#[test]
fn whitebox_gf7() {
    run_cases::<7>(CASES_GF7);
}

#[test]
fn whitebox_pos_gf5() {
    run_cases_pos::<5>(CASES_POS_GF5);
}

#[test]
fn whitebox_pos_gf7() {
    run_cases_pos::<7>(CASES_POS_GF7);
}

/// Self-contained branch coverage: the whitebox vectors (in-scope special cases)
/// plus a random sweep must together exercise every `ADD`/`DBL` formula branch.
#[test]
fn branch_coverage() {
    use super::test_support::{random_curve, random_valid_neg};
    use rand::{rngs::StdRng, SeedableRng};

    use super::test_support::random_valid_pos;

    // Deterministic in-scope special-case vectors (both bases).
    run_cases::<5>(CASES_GF5);
    run_cases::<7>(CASES_GF7);
    run_cases_pos::<5>(CASES_POS_GF5);
    run_cases_pos::<7>(CASES_POS_GF7);

    // Random sweep over small fields (hits the common/generic branches).
    for (seed, base) in [(909u64, 5u64), (910, 7)] {
        let mut rng = StdRng::seed_from_u64(seed);
        for _ in 0..200 {
            macro_rules! sweep {
                ($P:literal) => {{
                    let cc = random_curve::<PrimeField<$P>, _>(&mut rng);
                    for _ in 0..60 {
                        let (a, b) = (
                            random_valid_neg(&cc, &mut rng),
                            random_valid_neg(&cc, &mut rng),
                        );
                        double_neg(&from_generic(&a), &cc);
                        let (ca, cb) = (from_generic(&a), from_generic(&b));
                        if ca != cb {
                            add_neg(&ca, &cb, &cc);
                        }
                        let (a, b) = (
                            random_valid_pos(&cc, &mut rng),
                            random_valid_pos(&cc, &mut rng),
                        );
                        double_pos(&from_generic(&a), &cc);
                        let (ca, cb) = (from_generic(&a), from_generic(&b));
                        if ca != cb {
                            add_pos(&ca, &cb, &cc);
                        }
                    }
                }};
            }
            if base == 5 {
                sweep!(5);
            } else {
                sweep!(7);
            }
        }
    }

    let mut all: Vec<String> = Vec::new();
    for p in ["", "P"] {
        all.extend((0..=17).map(|i| format!("{p}DBL{i:02}")));
        all.extend((0..=58).map(|i| format!("{p}ADD{i:02}")));
    }
    let refs: Vec<&str> = all.iter().map(|s| s.as_str()).collect();
    let missing = super::coverage::missing(&refs);
    assert!(
        missing.is_empty(),
        "formula branches never exercised: {:?}",
        missing
    );
}
