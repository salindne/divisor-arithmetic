//! Blackbox tests for genus 2 split divisor arithmetic.
//!
//! Random valid balanced divisors are generated on random curves over growing
//! fields and the explicit formulas are checked against the generic Cantor
//! oracle (`crate::generic::split`), mirroring the Magma random testers.
//!
//! Phase A0 (this file, currently): proves the conversion + oracle + random
//! divisor harness is sound *before* any explicit formula is ported — the
//! decisive correctness anchor for everything that follows.

use crate::field::{Field, PrimeField};
use crate::generic::split;
use rand::rngs::StdRng;
use rand::SeedableRng;

use super::not_char2::{add_neg, add_pos, double_neg, double_pos, G};
use super::test_support::{
    assert_well_formed, deg0_neg, deg0_pos, deg1_neg, deg1_pos, from_generic, neutral_neg,
    random_curve, random_point, random_valid_neg, random_valid_pos, to_generic,
};

// ---------------------------------------------------------------------------
// Phase A0: harness self-validation (negative basis)
// ---------------------------------------------------------------------------

/// Generated random divisors are well-formed, and the explicit-coords
/// round-trip is the identity (including the balance weight `n`).
fn harness_roundtrip<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve::<PrimeField<P>, _>(&mut rng);
        let vn = cc.vn();
        for _ in 0..divisors {
            let d = random_valid_neg(&cc, &mut rng);
            assert_well_formed(&d, &cc);

            // Divisor -> coords -> Divisor must be the identity.
            let coords = from_generic(&d);
            let back = to_generic(&coords, &cc, &vn);
            assert_eq!(back, d, "neg round-trip changed the divisor");
            assert_eq!(coords.n, d.n, "round-trip lost the balance weight n");
        }
    }
}

/// The oracle itself (add_neg / double_neg) produces well-formed divisors, and
/// their coords round-trip — so the oracle is usable as the correctness anchor.
fn harness_oracle_valid<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve::<PrimeField<P>, _>(&mut rng);
        let f = cc.f_poly();
        let vn = cc.vn();
        for _ in 0..divisors {
            let d1 = random_valid_neg(&cc, &mut rng);
            let d2 = random_valid_neg(&cc, &mut rng);

            // D + 0 = D and 0 + D = D.
            let neutral = neutral_neg(&cc);
            let sum_id = split::add_neg(&d1, &neutral, &f, &crate::poly::Poly::zero(), &vn, G);
            assert_eq!(from_generic(&sum_id), from_generic(&d1), "D + 0 != D");

            // add and double produce well-formed reduced divisors.
            let sum = split::add_neg(&d1, &d2, &f, &crate::poly::Poly::zero(), &vn, G);
            assert_well_formed(&sum, &cc);
            let dbl = split::double_neg(&d1, &f, &crate::poly::Poly::zero(), &vn, G);
            assert_well_formed(&dbl, &cc);

            // commutativity of the oracle (sanity).
            let sum_rev = split::add_neg(&d2, &d1, &f, &crate::poly::Poly::zero(), &vn, G);
            assert_eq!(from_generic(&sum), from_generic(&sum_rev), "D1+D2 != D2+D1");

            // doubling equals adding to self.
            let dd = split::add_neg(&d1, &d1, &f, &crate::poly::Poly::zero(), &vn, G);
            assert_eq!(from_generic(&dbl), from_generic(&dd), "2D != D+D");
        }
    }
}

// ---------------------------------------------------------------------------
// Phase B: explicit double_neg == generic oracle, on (u, v, n)
// ---------------------------------------------------------------------------

/// Assert the explicit doubling matches the generic Cantor oracle for `d`.
fn check_double_neg<F: Field + std::fmt::Display>(
    d: &split::Divisor<F>,
    cc: &super::not_char2::CurveConstants<F>,
) {
    let f = cc.f_poly();
    let vn = cc.vn();
    let coords = from_generic(d);
    let got = double_neg(&coords, cc);
    let expected = from_generic(&split::double_neg(
        d,
        &f,
        &crate::poly::Poly::zero(),
        &vn,
        G,
    ));
    assert_eq!(got, expected, "double_neg mismatch for input {:?}", coords);
}

fn dbl_crosscheck<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve::<PrimeField<P>, _>(&mut rng);

        // Random (mostly degree-2) divisors.
        for _ in 0..divisors {
            check_double_neg(&random_valid_neg(&cc, &mut rng), &cc);
        }

        // Targeted degree-0 inputs with every valid weight n ∈ {0,1,2}.
        for n in 0..=2 {
            check_double_neg(&deg0_neg(&cc, n), &cc);
        }

        // Targeted degree-1 inputs (UP: n=0, DWN: n=1) from random points.
        for _ in 0..8 {
            if let Some((a, b)) = random_point(&cc, &mut rng) {
                check_double_neg(&deg1_neg(&cc, a, b, 0), &cc);
                check_double_neg(&deg1_neg(&cc, a, b, 1), &cc);
            }
        }
    }
}

#[test]
fn dbl_crosscheck_f7() {
    dbl_crosscheck::<7>(101, 12, 60);
}

#[test]
fn dbl_crosscheck_f31() {
    dbl_crosscheck::<31>(102, 12, 60);
}

#[test]
fn dbl_crosscheck_f127() {
    dbl_crosscheck::<127>(103, 8, 40);
}

#[test]
fn dbl_crosscheck_f8191() {
    dbl_crosscheck::<8191>(104, 3, 20);
}

/// All 18 DBL branches (`DBL00`..`DBL17`) are exercised by the cross-checks.
#[test]
fn dbl_branch_coverage() {
    // Run a broad sweep over small fields (more special-case collisions).
    dbl_crosscheck::<5>(201, 60, 120);
    dbl_crosscheck::<7>(202, 60, 120);
    dbl_crosscheck::<11>(203, 40, 120);
    dbl_crosscheck::<31>(204, 20, 120);
    let expected: Vec<&'static str> = (0..=17).map(|i| LABELS[i]).collect();
    let missing = super::coverage::missing(&expected);
    assert!(
        missing.is_empty(),
        "DBL branches never exercised: {:?}",
        missing
    );
}

const LABELS: [&str; 18] = [
    "DBL00", "DBL01", "DBL02", "DBL03", "DBL04", "DBL05", "DBL06", "DBL07", "DBL08", "DBL09",
    "DBL10", "DBL11", "DBL12", "DBL13", "DBL14", "DBL15", "DBL16", "DBL17",
];

// ---------------------------------------------------------------------------
// Phase D: positive basis — double_pos == generic oracle
// ---------------------------------------------------------------------------

fn check_double_pos<F: Field + std::fmt::Display>(
    d: &split::Divisor<F>,
    cc: &super::not_char2::CurveConstants<F>,
) {
    let f = cc.f_poly();
    let vpl = cc.vpl();
    let got = double_pos(&from_generic(d), cc);
    let expected = from_generic(&split::double_pos(
        d,
        &f,
        &crate::poly::Poly::zero(),
        &vpl,
        G,
    ));
    assert_eq!(
        got,
        expected,
        "double_pos mismatch for {:?}",
        from_generic(d)
    );
}

fn dbl_crosscheck_pos<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve::<PrimeField<P>, _>(&mut rng);
        let vpl = cc.vpl();
        for _ in 0..divisors {
            let d = random_valid_pos(&cc, &mut rng);
            assert_well_formed(&d, &cc);
            // round-trip incl. n in the positive basis
            assert_eq!(
                to_generic(&from_generic(&d), &cc, &vpl),
                d,
                "pos round-trip changed divisor"
            );
            check_double_pos(&d, &cc);
        }
        for n in 0..=2 {
            check_double_pos(&deg0_pos(&cc, n), &cc);
        }
        for _ in 0..8 {
            if let Some((a, b)) = random_point(&cc, &mut rng) {
                check_double_pos(&deg1_pos(&cc, a, b, 0), &cc);
                check_double_pos(&deg1_pos(&cc, a, b, 1), &cc);
            }
        }
    }
}

#[test]
fn dbl_pos_crosscheck_f7() {
    dbl_crosscheck_pos::<7>(501, 12, 60);
}

#[test]
fn dbl_pos_crosscheck_f31() {
    dbl_crosscheck_pos::<31>(502, 12, 60);
}

#[test]
fn dbl_pos_crosscheck_f127() {
    dbl_crosscheck_pos::<127>(503, 8, 40);
}

fn check_add_pos<F: Field + std::fmt::Display>(
    d1: &split::Divisor<F>,
    d2: &split::Divisor<F>,
    cc: &super::not_char2::CurveConstants<F>,
) {
    let f = cc.f_poly();
    let vpl = cc.vpl();
    let got = add_pos(&from_generic(d1), &from_generic(d2), cc);
    let expected = from_generic(&split::add_pos(
        d1,
        d2,
        &f,
        &crate::poly::Poly::zero(),
        &vpl,
        G,
    ));
    assert_eq!(
        got,
        expected,
        "add_pos mismatch:\n  d1={:?}\n  d2={:?}",
        from_generic(d1),
        from_generic(d2)
    );
}

/// Faithful port of the Magma pos random tester: random pos divisor pairs.
fn add_crosscheck_pos<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve::<PrimeField<P>, _>(&mut rng);
        for _ in 0..divisors {
            let d1 = random_valid_pos(&cc, &mut rng);
            let d2 = random_valid_pos(&cc, &mut rng);
            check_double_pos(&d1, &cc);
            if from_generic(&d1) != from_generic(&d2) {
                check_add_pos(&d1, &d2, &cc);
            }
        }
    }
}

#[test]
fn add_pos_crosscheck_f7() {
    add_crosscheck_pos::<7>(601, 8, 400);
}

#[test]
fn add_pos_crosscheck_f31() {
    add_crosscheck_pos::<31>(602, 8, 400);
}

#[test]
fn add_pos_crosscheck_f127() {
    add_crosscheck_pos::<127>(603, 5, 300);
}

// ---------------------------------------------------------------------------
// Phase C: explicit add_neg == generic oracle, on (u, v, n)
// ---------------------------------------------------------------------------

/// Assert the explicit addition matches the generic Cantor oracle for `d1 + d2`.
fn check_add_neg<F: Field + std::fmt::Display>(
    d1: &split::Divisor<F>,
    d2: &split::Divisor<F>,
    cc: &super::not_char2::CurveConstants<F>,
) {
    let f = cc.f_poly();
    let vn = cc.vn();
    let got = add_neg(&from_generic(d1), &from_generic(d2), cc);
    let expected = from_generic(&split::add_neg(
        d1,
        d2,
        &f,
        &crate::poly::Poly::zero(),
        &vn,
        G,
    ));
    assert_eq!(
        got,
        expected,
        "add_neg mismatch:\n  d1={:?}\n  d2={:?}",
        from_generic(d1),
        from_generic(d2)
    );
}

/// Faithful port of the Magma `nch2_splitG2_random.mag` tester: generate two
/// random canonical reduced divisors and check ADD (when distinct) + DBL.
fn add_crosscheck<const P: u64>(seed: u64, curves: usize, divisors: usize) {
    let mut rng = StdRng::seed_from_u64(seed);
    for _ in 0..curves {
        let cc = random_curve::<PrimeField<P>, _>(&mut rng);
        for _ in 0..divisors {
            let d1 = random_valid_neg(&cc, &mut rng);
            let d2 = random_valid_neg(&cc, &mut rng);
            check_double_neg(&d1, &cc);
            if from_generic(&d1) != from_generic(&d2) {
                check_add_neg(&d1, &d2, &cc);
            }
        }
    }
}

#[test]
fn add_crosscheck_f7() {
    add_crosscheck::<7>(301, 8, 400);
}

#[test]
fn add_crosscheck_f31() {
    add_crosscheck::<31>(302, 8, 400);
}

#[test]
fn add_crosscheck_f127() {
    add_crosscheck::<127>(303, 5, 300);
}

#[test]
fn add_crosscheck_f8191() {
    add_crosscheck::<8191>(304, 2, 150);
}

#[test]
fn a0_roundtrip_f7() {
    harness_roundtrip::<7>(1, 8, 40);
}

#[test]
fn a0_roundtrip_f31() {
    harness_roundtrip::<31>(2, 8, 40);
}

#[test]
fn a0_roundtrip_f127() {
    harness_roundtrip::<127>(3, 6, 30);
}

#[test]
fn a0_oracle_valid_f7() {
    harness_oracle_valid::<7>(11, 8, 40);
}

#[test]
fn a0_oracle_valid_f31() {
    harness_oracle_valid::<31>(12, 8, 40);
}

#[test]
fn a0_oracle_valid_f127() {
    harness_oracle_valid::<127>(13, 6, 30);
}
