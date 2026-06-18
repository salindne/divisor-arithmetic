//! Generic split model divisor arithmetic for arbitrary genus.
//!
//! Split curves have two points at infinity:
//! - `y² + h(x)·y = f(x)` where `deg(f) = 2g+2`
//!
//! Divisors are represented as `(u, v, w, n)` where:
//! - `u` is monic with `deg(u) ≤ g`
//! - `v` has a specific form related to the Vpl polynomial
//! - `w = (f - v·(v + h)) / u`
//! - `n` is the balancing weight
//!
//! These are faithful ports of `Add/Double/Adjust_SPLIT_{NEG,POS}` from
//! `reduced_basis_arithmetic.mag` and serve as the reference oracle for the
//! explicit g2 formulas. The basis polynomial passed in is `Vn = −Vpl − h` for
//! the negative basis and `Vpl` for the positive basis (matching the Magma
//! testers). `h` is threaded throughout, so these are correct for any `h`.
//!
//! Based on: Sebastian Lindner, 2019

use crate::field::Field;
use crate::poly::Poly;

/// A divisor on a split hyperelliptic curve.
///
/// Represented as `(u, v, w, n)` with balancing weight `n`.
#[derive(Clone, PartialEq, Eq)]
pub struct Divisor<F: Field> {
    pub u: Poly<F>,
    pub v: Poly<F>,
    pub w: Poly<F>,
    pub n: i32,
}

impl<F: Field + std::fmt::Display> std::fmt::Debug for Divisor<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Divisor")
            .field("u", &self.u)
            .field("v", &self.v)
            .field("w", &self.w)
            .field("n", &self.n)
            .finish()
    }
}

impl<F: Field> Divisor<F> {
    /// Create a new divisor
    pub fn new(u: Poly<F>, v: Poly<F>, w: Poly<F>, n: i32) -> Self {
        Self { u, v, w, n }
    }

    /// Create the neutral (identity) divisor `<1, V, f − V·(V + h), ⌈g/2⌉>`.
    pub fn neutral(f: &Poly<F>, v_basis: &Poly<F>, h: &Poly<F>, g: usize) -> Self {
        let vh = v_basis + h;
        let w = f - &(v_basis * &vh);
        Self {
            u: Poly::constant(F::one()),
            v: v_basis.clone(),
            w,
            n: ((g + 1) / 2) as i32,
        }
    }

    /// Get the degree of this divisor (degree of u)
    pub fn degree(&self) -> i32 {
        self.u.deg()
    }
}

/// Compute the Vpl polynomial for a split curve (negative basis caller should
/// negate and subtract `h`). Correct for monic `f` over odd characteristic;
/// for arbitrary/char-2 curves callers pass `Vpl` in from their own precompute.
pub fn compute_vpl<F: Field>(f: &Poly<F>, g: usize) -> Poly<F> {
    let fl = f.coeff(2 * g + 2);
    // Leading term of Vpl solves fl − x² = 0. For monic f this is 1.
    let vl = if fl.is_one() { F::one() } else { fl };
    let mut vpl = Poly::monomial(vl, g + 1);
    let two = F::one() + F::one();
    let dinv = (two * vl).inv();
    for i in (0..=g).rev() {
        let vpl_sq = &vpl * &vpl;
        let diff = f - &vpl_sq;
        let coeff = diff.coeff(g + 1 + i);
        vpl = vpl + Poly::monomial(dinv * coeff, i);
    }
    vpl
}

// ---------------------------------------------------------------------------
// Shared compose / finish helpers (identical between Add and Double except the
// composition step, and between neg and pos except the reduce sign test).
// ---------------------------------------------------------------------------

/// Composition step for `Add` (semi-reduced `(u, v, w, n)`), per Magma.
fn compose_add<F: Field>(d1: &Divisor<F>, d2: &Divisor<F>, h: &Poly<F>, g: usize) -> Divisor<F> {
    let v1 = &d1.v;
    let t1 = v1 + h; // v1 + h
    let (s, a1, _b1) = d1.u.xgcd(&d2.u);
    let mut k = (a1 * (&d2.v - v1)).rem(&d2.u);

    let mut u1 = d1.u.clone();
    let mut u2 = d2.u.clone();
    let mut w1 = d1.w.clone();
    let mut s_deg = s.deg();

    if !s.is_one() {
        let v2_plus_t1 = &d2.v + &t1; // v2 + v1 + h
        let (s2, a2, b2) = s.xgcd(&v2_plus_t1);
        if !s2.is_one() {
            u1 = u1.exact_div(&s2);
            u2 = u2.exact_div(&s2);
            k = (a2 * k + b2 * w1.clone()).rem(&u2);
            w1 = w1 * s2.clone();
            s_deg = s2.deg();
        } else {
            k = (a2 * k + b2 * w1.clone()).rem(&u2);
            s_deg = 0;
        }
    } else {
        s_deg = 0;
    }

    let t = &u1 * &k;
    let u = &u1 * &u2;
    let v = v1 + &t;
    let t1_plus_v = &t1 + &v; // (v1 + h) + v
    let w = (w1 - k * t1_plus_v).exact_div(&u2);
    let n = d1.n + d2.n + s_deg - ((g as i32 + 1) / 2);
    Divisor::new(u, v, w, n)
}

/// Composition step for `Double`, per Magma.
fn compose_double<F: Field>(d1: &Divisor<F>, h: &Poly<F>, g: usize) -> Divisor<F> {
    let v1 = &d1.v;
    let t1 = &(v1 + v1) + h; // 2·v1 + h
    let (s, _a1, b1) = d1.u.xgcd(&t1);

    let mut u1 = d1.u.clone();
    let mut w1 = d1.w.clone();
    let k;
    let s_deg;
    if !s.is_one() {
        u1 = u1.exact_div(&s);
        k = (b1 * w1.clone()).rem(&u1);
        w1 = w1 * s.clone();
        s_deg = s.deg();
    } else {
        k = (b1 * w1.clone()).rem(&u1);
        s_deg = 0;
    }

    let t = &u1 * &k;
    let u = &u1 * &u1;
    let v = v1 + &t;
    let t1_plus_t = &t1 + &t; // (2v1 + h) + T
    let w = (w1 - k * t1_plus_t).exact_div(&u1);
    let n = 2 * d1.n + s_deg - ((g as i32 + 1) / 2);
    Divisor::new(u, v, w, n)
}

/// Normalize + reduce a semi-reduced divisor to `deg(u) ≤ g+1`, monic. `pos`
/// selects the positive-basis reduce sign test. `vb` is the basis polynomial.
fn finish<F: Field>(mut d: Divisor<F>, h: &Poly<F>, vb: &Poly<F>, g: usize, pos: bool) -> Divisor<F> {
    let g1 = g as i32 + 1;
    let lc_v = vb.leading_coeff(); // lc(V)
    let lc_neg = (-(vb.clone()) - h.clone()).leading_coeff(); // lc(−V − h)

    // Normalize
    if d.v.deg() >= d.u.deg() {
        let (q, r) = (vb - &d.v).div_rem(&d.u);
        let tv = vb - &r;
        let vh = &d.v + h;
        let vht = &vh + &tv;
        d.w = d.w - q * vht;
        d.v = tv;
    }

    // Reduce
    while d.u.deg() > g1 {
        let vdeg = d.v.deg();
        let lcv = d.v.leading_coeff();
        // neg: first test lc(−V−h), second lc(V); pos: swapped.
        let (first, second) = if pos { (lc_v, lc_neg) } else { (lc_neg, lc_v) };
        if vdeg == g1 && lcv == first {
            d.n += d.u.deg() - g1;
        } else if vdeg == g1 && lcv == second {
            d.n += g1 - d.w.deg();
        } else {
            d.n += (d.u.deg() - d.w.deg()) / 2;
        }

        let ou = d.u.clone();
        d.u = d.w.clone();
        let vv = vb + &d.v;
        let vvh = &vv + h; // V + v + h
        let (q, r) = vvh.div_rem(&d.u);
        let tv = vb - &r;
        d.w = ou - q * (&tv - &d.v);
        d.v = tv;
    }

    let lc = d.u.leading_coeff();
    d.w = d.w * lc;
    d.u = d.u.make_monic();
    d
}

// ---------------------------------------------------------------------------
// Negative reduced basis
// ---------------------------------------------------------------------------

/// `Adjust_SPLIT_NEG`: reduce the balance weight into `0 ≤ n ≤ g − deg(u)`.
pub fn adjust_neg<F: Field>(mut d: Divisor<F>, _f: &Poly<F>, h: &Poly<F>, v_neg: &Poly<F>, g: usize) -> Divisor<F> {
    let g_i32 = g as i32;

    if d.n < 0 {
        // UP adjust
        while d.n < 0 {
            let ou = d.u.clone();
            d.u = d.w.clone();
            let vv = v_neg + &d.v;
            let vvh = &vv + h;
            let (q, r) = vvh.div_rem(&d.u);
            let tv = v_neg - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
            d.n = d.n + g_i32 + 1 - d.u.deg();
        }
        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    } else if d.n > g_i32 - d.u.deg() {
        // DWN adjust: convert to positive basis Vp = −V − h, step, convert back.
        let v_pos = -(v_neg.clone()) - h.clone();
        let t = &v_pos - v_neg;
        let (q, r) = t.div_rem(&d.u);
        let tv = &(&d.v + &t) - &r;
        let vh = &d.v + h;
        let vht = &vh + &tv;
        d.w = d.w - q * vht;
        d.v = tv;

        while d.n > g_i32 - d.u.deg() + 1 {
            d.n = d.n + d.u.deg() - (g_i32 + 1);
            let ou = d.u.clone();
            d.u = d.w.clone();
            let vv = &v_pos + &d.v;
            let vvh = &vv + h;
            let (q, r) = vvh.div_rem(&d.u);
            let tv = &v_pos - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
        }

        // Final step back into negative reduced.
        d.n = d.n + d.u.deg() - (g_i32 + 1);
        let ou = d.u.clone();
        d.u = d.w.clone();
        let vv = v_neg + &d.v;
        let vvh = &vv + h;
        let (q, r) = vvh.div_rem(&d.u);
        let tv = v_neg - &r;
        d.w = ou - q * (&tv - &d.v);
        d.v = tv;

        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    }
    d
}

/// `Add_SPLIT_NEG`.
pub fn add_neg<F: Field>(d1: &Divisor<F>, d2: &Divisor<F>, f: &Poly<F>, h: &Poly<F>, v_neg: &Poly<F>, g: usize) -> Divisor<F> {
    let d = compose_add(d1, d2, h, g);
    let d = finish(d, h, v_neg, g, false);
    adjust_neg(d, f, h, v_neg, g)
}

/// `Double_SPLIT_NEG`.
pub fn double_neg<F: Field>(d1: &Divisor<F>, f: &Poly<F>, h: &Poly<F>, v_neg: &Poly<F>, g: usize) -> Divisor<F> {
    let d = compose_double(d1, h, g);
    let d = finish(d, h, v_neg, g, false);
    adjust_neg(d, f, h, v_neg, g)
}

// ---------------------------------------------------------------------------
// Positive reduced basis
// ---------------------------------------------------------------------------

/// `Adjust_SPLIT_POS`: reduce the balance weight into `0 ≤ n ≤ g − deg(u)`.
pub fn adjust_pos<F: Field>(mut d: Divisor<F>, _f: &Poly<F>, h: &Poly<F>, v_pos: &Poly<F>, g: usize) -> Divisor<F> {
    let g_i32 = g as i32;

    if d.n > g_i32 - d.u.deg() {
        // DWN adjust
        while d.n > g_i32 - d.u.deg() {
            d.n = d.n + d.u.deg() - (g_i32 + 1);
            let ou = d.u.clone();
            d.u = d.w.clone();
            let vv = v_pos + &d.v;
            let vvh = &vv + h;
            let (q, r) = vvh.div_rem(&d.u);
            let tv = v_pos - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
        }
        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    } else if d.n < 0 {
        // Convert to negative basis Vp = −V − h, step, convert back up.
        let v_neg = -(v_pos.clone()) - h.clone();
        let t = &v_neg - v_pos;
        let (q, r) = t.div_rem(&d.u);
        let tv = &(&d.v + &t) - &r;
        let vh = &d.v + h;
        let vht = &vh + &tv;
        d.w = d.w - q * vht;
        d.v = tv;

        while d.n < -1 {
            let ou = d.u.clone();
            d.u = d.w.clone();
            let vv = &v_neg + &d.v;
            let vvh = &vv + h;
            let (q, r) = vvh.div_rem(&d.u);
            let tv = &v_neg - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
            d.n = d.n + g_i32 + 1 - d.u.deg();
        }

        // UP adjust back into positive reduced.
        let ou = d.u.clone();
        d.u = d.w.clone();
        let vv = v_pos + &d.v;
        let vvh = &vv + h;
        let (q, r) = vvh.div_rem(&d.u);
        let tv = v_pos - &r;
        d.w = ou - q * (&tv - &d.v);
        d.v = tv;
        d.n = d.n + g_i32 + 1 - d.u.deg();

        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    }
    d
}

/// `Add_SPLIT_POS`.
pub fn add_pos<F: Field>(d1: &Divisor<F>, d2: &Divisor<F>, f: &Poly<F>, h: &Poly<F>, v_pos: &Poly<F>, g: usize) -> Divisor<F> {
    let d = compose_add(d1, d2, h, g);
    let d = finish(d, h, v_pos, g, true);
    adjust_pos(d, f, h, v_pos, g)
}

/// `Double_SPLIT_POS`.
pub fn double_pos<F: Field>(d1: &Divisor<F>, f: &Poly<F>, h: &Poly<F>, v_pos: &Poly<F>, g: usize) -> Divisor<F> {
    let d = compose_double(d1, h, g);
    let d = finish(d, h, v_pos, g, true);
    adjust_pos(d, f, h, v_pos, g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;

    type F7 = PrimeField<7>;

    fn poly(coeffs: &[u64]) -> Poly<F7> {
        Poly::from_coeffs(coeffs.iter().map(|&c| F7::new(c)).collect())
    }

    /// Create a genus 2 split curve over F7 (h = 0).
    fn test_curve_g2_split() -> Poly<F7> {
        poly(&[1, 2, 4, 1, 3, 2, 1])
    }

    fn find_valid_split_divisor(f: &Poly<F7>, vpl: &Poly<F7>, _g: usize) -> Option<Divisor<F7>> {
        for a_val in 0..7u64 {
            let a = F7::new(a_val);
            let fa = f.eval(a);
            if fa.is_zero() || fa.pow(3).is_one() {
                let b = if fa.is_zero() { F7::zero() } else { fa.pow(2) };
                if b * b == fa {
                    let u = poly(&[(7 - a_val) % 7, 1]);
                    let b_poly = Poly::constant(b);
                    let vpl_minus_b = vpl - &b_poly;
                    let v = vpl - &vpl_minus_b.rem(&u);
                    let v_sq = &v * &v;
                    let diff = f - &v_sq;
                    let (_, rem) = diff.div_rem(&u);
                    if rem.is_zero() {
                        let w = diff.exact_div(&u);
                        return Some(Divisor::new(u, v, w, 0));
                    }
                }
            }
        }
        None
    }

    #[test]
    fn test_compute_vpl() {
        let f = test_curve_g2_split();
        let g = 2;
        let vpl = compute_vpl(&f, g);
        assert_eq!(vpl.deg(), 3);
        let vpl_sq = &vpl * &vpl;
        let diff = f - vpl_sq;
        assert!(diff.deg() <= g as i32);
    }

    #[test]
    fn test_neutral_split() {
        let f = test_curve_g2_split();
        let g = 2;
        let h = Poly::zero();
        let vpl = compute_vpl(&f, g);
        let neutral = Divisor::neutral(&f, &vpl, &h, g);
        assert!(neutral.u.is_one());
        assert_eq!(neutral.v, vpl);
    }

    #[test]
    fn test_double_pos_consistency() {
        let f = test_curve_g2_split();
        let g = 2;
        let h = Poly::zero();
        let vpl = compute_vpl(&f, g);
        if let Some(d) = find_valid_split_divisor(&f, &vpl, g) {
            let doubled = double_pos(&d, &f, &h, &vpl, g);
            assert!(doubled.u.deg() <= g as i32);
            let v_sq = &doubled.v * &doubled.v;
            let expected_w = (f - v_sq).exact_div(&doubled.u);
            assert_eq!(doubled.w, expected_w);
        }
    }

    #[test]
    fn test_double_neg_consistency() {
        let f = test_curve_g2_split();
        let g = 2;
        let h = Poly::zero();
        let vpl = compute_vpl(&f, g);
        let v_neg = -vpl.clone();
        if let Some(mut d) = find_valid_split_divisor(&f, &vpl, g) {
            let t = &v_neg - &vpl;
            let (q, r) = t.div_rem(&d.u);
            let tv = &d.v + &t - r;
            d.w = d.w - q * (&d.v + &tv);
            d.v = tv;
            let doubled = double_neg(&d, &f, &h, &v_neg, g);
            assert!(doubled.u.deg() <= g as i32);
            let v_sq = &doubled.v * &doubled.v;
            let expected_w = (f - v_sq).exact_div(&doubled.u);
            assert_eq!(doubled.w, expected_w);
        }
    }

    #[test]
    fn test_add_pos_neutral() {
        let f = test_curve_g2_split();
        let g = 2;
        let h = Poly::zero();
        let vpl = compute_vpl(&f, g);
        let neutral = Divisor::neutral(&f, &vpl, &h, g);
        if let Some(d) = find_valid_split_divisor(&f, &vpl, g) {
            let result = add_pos(&d, &neutral, &f, &h, &vpl, g);
            assert!(result.u.deg() <= g as i32);
            let v_sq = &result.v * &result.v;
            let expected_w = (f - v_sq).exact_div(&result.u);
            assert_eq!(result.w, expected_w);
        }
    }

    #[test]
    fn test_double_pos_vs_add_pos() {
        let f = test_curve_g2_split();
        let g = 2;
        let h = Poly::zero();
        let vpl = compute_vpl(&f, g);
        if let Some(d) = find_valid_split_divisor(&f, &vpl, g) {
            let doubled = double_pos(&d, &f, &h, &vpl, g);
            let added = add_pos(&d, &d, &f, &h, &vpl, g);
            assert_eq!(doubled.u, added.u);
        }
    }
}
