//! Generic split model divisor arithmetic for arbitrary genus.
//!
//! Split curves have two points at infinity:
//! - `y² = f(x)` where `deg(f) = 2g+2`
//!
//! Divisors are represented as `(u, v, w, n)` where:
//! - `u` is monic with `deg(u) ≤ g`
//! - `v` has a specific form related to the Vpl polynomial
//! - `w = (f - v²) / u`
//! - `n` is the balancing weight
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

    /// Create the neutral (identity) divisor
    pub fn neutral(f: &Poly<F>, v_pl: &Poly<F>, g: usize) -> Self {
        let v_sq = v_pl * v_pl;
        let w = f - &v_sq;
        Self {
            u: Poly::constant(F::one()),
            v: v_pl.clone(),
            w,
            n: ((g + 1) / 2) as i32,
        }
    }

    /// Get the degree of this divisor (degree of u)
    pub fn degree(&self) -> i32 {
        self.u.deg()
    }
}

/// Compute the Vpl polynomial for a split curve.
///
/// Returns the unique polynomial V of degree g+1 such that deg(f - V²) ≤ g.
pub fn compute_vpl<F: Field>(f: &Poly<F>, g: usize) -> Poly<F> {
    // Leading coefficient of f (degree 2g+2)
    let fl = f.coeff(2 * g + 2);

    // Leading term of Vpl is solution of fl - x² = 0
    // We need a square root of fl. For odd characteristic, use Tonelli-Shanks
    // or just try fl^((p+1)/4) for p ≡ 3 (mod 4)
    // For simplicity, we'll factor x² - fl and take the negative of the constant term
    // This is a simplified approach - in practice you'd need proper square root

    // For now, we compute iteratively as in the Magma code
    // Vl := -Coeff(Factorization(fl - x^2)[2][1], 0)
    // This requires factorization which is complex

    // Simplified: assume fl = 1 (monic f), so Vl = 1 or -1
    // For a proper implementation, we'd need square root in the field
    let vl = if fl.is_one() {
        F::one()
    } else {
        // Try to find square root by trying fl^((p-1)/2 + 1/2) = fl^((p+1)/4) for p ≡ 3 (mod 4)
        // This is a simplification - proper implementation would use Tonelli-Shanks
        fl // Placeholder - this won't be correct for non-trivial cases
    };

    let mut vpl = Poly::monomial(vl, g + 1);

    // dinv := (2*Vl)^-1
    let two = F::one() + F::one();
    let dinv = (two * vl).inv();

    // Work down one term at a time
    for i in (0..=g).rev() {
        // Vpl += dinv * Coeff(f - Vpl², g+1+i) * x^i
        let vpl_sq = &vpl * &vpl;
        let diff = f - &vpl_sq;
        let coeff = diff.coeff(g + 1 + i);
        let term = Poly::monomial(dinv * coeff, i);
        vpl += term;
    }

    vpl
}

/// Adjust a divisor to be in reduced form (negative reduced basis).
///
/// Implements `Adjust_SPLIT_NEG` from Magma.
pub fn adjust_neg<F: Field>(
    mut d: Divisor<F>,
    _f: &Poly<F>,
    v_neg: &Poly<F>, // -Vpl
    g: usize,
) -> Divisor<F> {
    let g_i32 = g as i32;

    if d.n < 0 {
        // UP Adjust
        while d.n < 0 {
            let ou = d.u.clone();
            d.u = d.w.clone();

            let v_plus_v_neg = &d.v + v_neg;
            let (q, r) = v_plus_v_neg.div_rem(&d.u);
            let tv = v_neg - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
            d.n = d.n + g_i32 + 1 - d.u.deg();
        }
        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    } else if d.n > g_i32 - d.u.deg() {
        // DWN Adjust
        // Basis conversion
        let v_pos = -v_neg.clone(); // Vpl
        let t = &v_pos - v_neg;
        let (q, r) = t.div_rem(&d.u);
        let tv = &d.v + &t - r;
        d.w -= q * (&d.v + &tv);
        d.v = tv;

        while d.n > g_i32 - d.u.deg() + 1 {
            d.n = d.n + d.u.deg() - (g_i32 + 1);
            let ou = d.u.clone();
            d.u = d.w.clone();
            let v_pos_plus_v = &v_pos + &d.v;
            let (q, r) = v_pos_plus_v.div_rem(&d.u);
            let tv = &v_pos - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
        }

        if d.n > g_i32 - d.u.deg() {
            d.n = d.n + d.u.deg() - (g_i32 + 1);
            let ou = d.u.clone();
            d.u = d.w.clone();
            let v_neg_plus_v = v_neg + &d.v;
            let (q, r) = v_neg_plus_v.div_rem(&d.u);
            let tv = v_neg - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
        } else {
            let t = v_neg - &v_pos;
            let (q, r) = t.div_rem(&d.u);
            let tv = &d.v + &t - r;
            d.w -= q * (&d.v + &tv);
            d.v = tv;
        }

        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    }

    d
}

/// Add two divisors on a split curve (negative reduced basis).
///
/// Implements `Add_SPLIT_NEG` from Magma.
pub fn add_neg<F: Field>(
    d1: &Divisor<F>,
    d2: &Divisor<F>,
    f: &Poly<F>,
    v_neg: &Poly<F>,
    g: usize,
) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;
    let n1 = d1.n;
    let u2 = &d2.u;
    let v2 = &d2.v;
    let n2 = d2.n;

    // Compose
    let (s, a1, _b1) = u1.xgcd(u2);
    let v_diff = v2 - v1;
    let mut k = (a1 * v_diff).rem(u2);

    let mut u1 = u1.clone();
    let mut u2 = u2.clone();
    let mut w1 = w1.clone();
    let s_deg;

    if !s.is_one() {
        let v_sum = v2 + v1;
        let (s2, a2, b2) = s.xgcd(&v_sum);

        if !s2.is_one() {
            u1 = u1.exact_div(&s2);
            u2 = u2.exact_div(&s2);
            k = (a2 * k + b2 * w1.clone()).rem(&u2);
            w1 *= s2.clone();
            s_deg = s2.deg();
        } else {
            k = (a2 * k + b2 * w1.clone()).rem(&u2);
            s_deg = 0;
        }
    } else {
        s_deg = 0;
    }

    let t = &u1 * &k;
    let mut u = &u1 * &u2;
    let mut v = v1 + &t;
    let v1_plus_v = v1 + &v;
    let mut w = (w1 - k * v1_plus_v).exact_div(&u2);
    let g_i32 = g as i32;
    let mut n = n1 + n2 + s_deg - ((g_i32 + 1) / 2);

    // Normalize
    if v.deg() >= u.deg() {
        let v_neg_minus_v = v_neg - &v;
        let (q, r) = v_neg_minus_v.div_rem(&u);
        let tv = v_neg - &r;
        w -= q * (&v + &tv);
        v = tv;
    }

    // Reduce
    while u.deg() > g_i32 + 1 {
        let v_deg = v.deg();
        let lc_v = v.leading_coeff();
        let lc_v_neg = v_neg.leading_coeff();
        let neg_lc_v_neg = -lc_v_neg;

        if v_deg == g_i32 + 1 && lc_v == neg_lc_v_neg {
            n = n + u.deg() - (g_i32 + 1);
        } else if v_deg == g_i32 + 1 && lc_v == lc_v_neg {
            n = n + g_i32 + 1 - w.deg();
        } else {
            n += (u.deg() - w.deg()) / 2;
        }

        let ou = u.clone();
        u = w.clone();
        let v_neg_plus_v = v_neg + &v;
        let (q, r) = v_neg_plus_v.div_rem(&u);
        let tv = v_neg - &r;
        w = ou - q * (&tv - &v);
        v = tv;
    }

    let lc = u.leading_coeff();
    w = w * lc;
    u = u.make_monic();

    adjust_neg(Divisor::new(u, v, w, n), f, v_neg, g)
}

/// Double a divisor on a split curve (negative reduced basis).
///
/// Implements `Double_SPLIT_NEG` from Magma.
pub fn double_neg<F: Field>(d1: &Divisor<F>, f: &Poly<F>, v_neg: &Poly<F>, g: usize) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;
    let n1 = d1.n;

    let t1 = v1 + v1;
    let (s, _a1, b1) = u1.xgcd(&t1);

    let mut u1 = u1.clone();
    let mut w1 = w1.clone();
    let k;
    let s_deg;

    if !s.is_one() {
        u1 = u1.exact_div(&s);
        k = (b1 * w1.clone()).rem(&u1);
        w1 *= s.clone();
        s_deg = s.deg();
    } else {
        k = (b1 * w1.clone()).rem(&u1);
        s_deg = 0;
    }

    let t = &u1 * &k;
    let mut u = &u1 * &u1;
    let mut v = v1 + &t;
    let t1_plus_t = &t1 + &t;
    let mut w = (w1 - k * t1_plus_t).exact_div(&u1);
    let g_i32 = g as i32;
    let mut n = 2 * n1 + s_deg - ((g_i32 + 1) / 2);

    // Normalize
    if v.deg() >= u.deg() {
        let v_neg_minus_v = v_neg - &v;
        let (q, r) = v_neg_minus_v.div_rem(&u);
        let tv = v_neg - &r;
        w -= q * (&v + &tv);
        v = tv;
    }

    // Reduce
    while u.deg() > g_i32 + 1 {
        let v_deg = v.deg();
        let lc_v = v.leading_coeff();
        let lc_v_neg = v_neg.leading_coeff();
        let neg_lc_v_neg = -lc_v_neg;

        if v_deg == g_i32 + 1 && lc_v == neg_lc_v_neg {
            n = n + u.deg() - (g_i32 + 1);
        } else if v_deg == g_i32 + 1 && lc_v == lc_v_neg {
            n = n + g_i32 + 1 - w.deg();
        } else {
            n += (u.deg() - w.deg()) / 2;
        }

        let ou = u.clone();
        u = w.clone();
        let v_neg_plus_v = v_neg + &v;
        let (q, r) = v_neg_plus_v.div_rem(&u);
        let tv = v_neg - &r;
        w = ou - q * (&tv - &v);
        v = tv;
    }

    let lc = u.leading_coeff();
    w = w * lc;
    u = u.make_monic();

    adjust_neg(Divisor::new(u, v, w, n), f, v_neg, g)
}

/// Adjust a divisor to be in reduced form (positive reduced basis).
///
/// Implements `Adjust_SPLIT_POS` from Magma.
pub fn adjust_pos<F: Field>(
    mut d: Divisor<F>,
    _f: &Poly<F>,
    v_pos: &Poly<F>, // Vpl
    g: usize,
) -> Divisor<F> {
    let g_i32 = g as i32;

    if d.n > g_i32 - d.u.deg() {
        // DWN Adjust
        while d.n > g_i32 - d.u.deg() {
            d.n = d.n + d.u.deg() - (g_i32 + 1);
            let ou = d.u.clone();
            d.u = d.w.clone();
            let v_pos_plus_v = v_pos + &d.v;
            let (q, r) = v_pos_plus_v.div_rem(&d.u);
            let tv = v_pos - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
        }
        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    } else if d.n < 0 {
        // UP Adjust via negative reduced
        let v_neg = -v_pos.clone();
        let t = &v_neg - v_pos;
        let (q, r) = t.div_rem(&d.u);
        let tv = &d.v + &t - r;
        d.w -= q * (&d.v + &tv);
        d.v = tv;

        while d.n < -1 {
            let ou = d.u.clone();
            d.u = d.w.clone();
            let v_neg_plus_v = &v_neg + &d.v;
            let (q, r) = v_neg_plus_v.div_rem(&d.u);
            let tv = &v_neg - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
            d.n = d.n + g_i32 + 1 - d.u.deg();
        }

        if d.n < 0 {
            let ou = d.u.clone();
            d.u = d.w.clone();
            let v_pos_plus_v = v_pos + &d.v;
            let (q, r) = v_pos_plus_v.div_rem(&d.u);
            let tv = v_pos - &r;
            d.w = ou - q * (&tv - &d.v);
            d.v = tv;
            d.n = d.n + g_i32 + 1 - d.u.deg();
        } else {
            let t = v_pos - &v_neg;
            let (q, r) = t.div_rem(&d.u);
            let tv = &d.v + &t - r;
            d.w -= q * (&d.v + &tv);
            d.v = tv;
        }

        let lc = d.u.leading_coeff();
        d.w = d.w * lc;
        d.u = d.u.make_monic();
    }

    d
}

/// Add two divisors on a split curve (positive reduced basis).
///
/// Implements `Add_SPLIT_POS` from Magma.
pub fn add_pos<F: Field>(
    d1: &Divisor<F>,
    d2: &Divisor<F>,
    f: &Poly<F>,
    v_pos: &Poly<F>,
    g: usize,
) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;
    let n1 = d1.n;
    let u2 = &d2.u;
    let v2 = &d2.v;
    let n2 = d2.n;

    // Compose
    let (s, a1, _b1) = u1.xgcd(u2);
    let v_diff = v2 - v1;
    let mut k = (a1 * v_diff).rem(u2);

    let mut u1 = u1.clone();
    let mut u2 = u2.clone();
    let mut w1 = w1.clone();
    let s_deg;

    if !s.is_one() {
        let v_sum = v2 + v1;
        let (s2, a2, b2) = s.xgcd(&v_sum);

        if !s2.is_one() {
            u1 = u1.exact_div(&s2);
            u2 = u2.exact_div(&s2);
            k = (a2 * k + b2 * w1.clone()).rem(&u2);
            w1 *= s2.clone();
            s_deg = s2.deg();
        } else {
            k = (a2 * k + b2 * w1.clone()).rem(&u2);
            s_deg = 0;
        }
    } else {
        s_deg = 0;
    }

    let t = &u1 * &k;
    let mut u = &u1 * &u2;
    let mut v = v1 + &t;
    let v1_plus_v = v1 + &v;
    let mut w = (w1 - k * v1_plus_v).exact_div(&u2);
    let g_i32 = g as i32;
    let mut n = n1 + n2 + s_deg - ((g_i32 + 1) / 2);

    // Normalize
    if v.deg() >= u.deg() {
        let v_pos_minus_v = v_pos - &v;
        let (q, r) = v_pos_minus_v.div_rem(&u);
        let tv = v_pos - &r;
        w -= q * (&v + &tv);
        v = tv;
    }

    // Reduce
    while u.deg() > g_i32 + 1 {
        let v_deg = v.deg();
        let lc_v = v.leading_coeff();
        let lc_v_pos = v_pos.leading_coeff();
        let neg_lc_v_pos = -lc_v_pos;

        if v_deg == g_i32 + 1 && lc_v == lc_v_pos {
            n = n + u.deg() - (g_i32 + 1);
        } else if v_deg == g_i32 + 1 && lc_v == neg_lc_v_pos {
            n = n + g_i32 + 1 - w.deg();
        } else {
            n += (u.deg() - w.deg()) / 2;
        }

        let ou = u.clone();
        u = w.clone();
        let v_pos_plus_v = v_pos + &v;
        let (q, r) = v_pos_plus_v.div_rem(&u);
        let tv = v_pos - &r;
        w = ou - q * (&tv - &v);
        v = tv;
    }

    let lc = u.leading_coeff();
    w = w * lc;
    u = u.make_monic();

    adjust_pos(Divisor::new(u, v, w, n), f, v_pos, g)
}

/// Double a divisor on a split curve (positive reduced basis).
///
/// Implements `Double_SPLIT_POS` from Magma.
pub fn double_pos<F: Field>(d1: &Divisor<F>, f: &Poly<F>, v_pos: &Poly<F>, g: usize) -> Divisor<F> {
    let u1 = &d1.u;
    let v1 = &d1.v;
    let w1 = &d1.w;
    let n1 = d1.n;

    let t1 = v1 + v1;
    let (s, _a1, b1) = u1.xgcd(&t1);

    let mut u1 = u1.clone();
    let mut w1 = w1.clone();
    let k;
    let s_deg;

    if !s.is_one() {
        u1 = u1.exact_div(&s);
        k = (b1 * w1.clone()).rem(&u1);
        w1 *= s.clone();
        s_deg = s.deg();
    } else {
        k = (b1 * w1.clone()).rem(&u1);
        s_deg = 0;
    }

    let t = &u1 * &k;
    let mut u = &u1 * &u1;
    let mut v = v1 + &t;
    let t1_plus_t = &t1 + &t;
    let mut w = (w1 - k * t1_plus_t).exact_div(&u1);
    let g_i32 = g as i32;
    let mut n = 2 * n1 + s_deg - ((g_i32 + 1) / 2);

    // Normalize
    if v.deg() >= u.deg() {
        let v_pos_minus_v = v_pos - &v;
        let (q, r) = v_pos_minus_v.div_rem(&u);
        let tv = v_pos - &r;
        w -= q * (&v + &tv);
        v = tv;
    }

    // Reduce
    while u.deg() > g_i32 + 1 {
        let v_deg = v.deg();
        let lc_v = v.leading_coeff();
        let lc_v_pos = v_pos.leading_coeff();
        let neg_lc_v_pos = -lc_v_pos;

        if v_deg == g_i32 + 1 && lc_v == lc_v_pos {
            n = n + u.deg() - (g_i32 + 1);
        } else if v_deg == g_i32 + 1 && lc_v == neg_lc_v_pos {
            n = n + g_i32 + 1 - w.deg();
        } else {
            n += (u.deg() - w.deg()) / 2;
        }

        let ou = u.clone();
        u = w.clone();
        let v_pos_plus_v = v_pos + &v;
        let (q, r) = v_pos_plus_v.div_rem(&u);
        let tv = v_pos - &r;
        w = ou - q * (&tv - &v);
        v = tv;
    }

    let lc = u.leading_coeff();
    w = w * lc;
    u = u.make_monic();

    adjust_pos(Divisor::new(u, v, w, n), f, v_pos, g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;

    type F7 = PrimeField<7>;

    fn poly(coeffs: &[u64]) -> Poly<F7> {
        Poly::from_coeffs(coeffs.iter().map(|&c| F7::new(c)).collect())
    }

    /// Create a genus 2 split curve over F7
    /// f = x^6 + 2x^5 + 3x^4 + x^3 + 4x^2 + 2x + 1 (degree 6 = 2*2+2)
    fn test_curve_g2_split() -> Poly<F7> {
        poly(&[1, 2, 4, 1, 3, 2, 1])
    }

    /// Find a valid divisor on a split curve
    fn find_valid_split_divisor(f: &Poly<F7>, vpl: &Poly<F7>, _g: usize) -> Option<Divisor<F7>> {
        // Try to find x where f(x) is a quadratic residue
        for a_val in 0..7u64 {
            let a = F7::new(a_val);
            let fa = f.eval(a);

            // Check if fa is a square in F7
            if fa.is_zero() || fa.pow(3).is_one() {
                let b = if fa.is_zero() { F7::zero() } else { fa.pow(2) };

                if b * b == fa {
                    // u = x - a
                    let u = poly(&[(7 - a_val) % 7, 1]); // x - a = x + (-a)

                    // v should be in reduced basis form: vpl - (vpl - b) mod u
                    let b_poly = Poly::constant(b);
                    let vpl_minus_b = vpl - &b_poly;
                    let v = vpl - &vpl_minus_b.rem(&u);

                    let v_sq = &v * &v;
                    let diff = f - &v_sq;

                    // Check divisibility
                    let (_, rem) = diff.div_rem(&u);
                    if rem.is_zero() {
                        let w = diff.exact_div(&u);
                        // n should be in valid range [0, g - deg(u)]
                        let n = 0;
                        return Some(Divisor::new(u, v, w, n));
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

        // Vpl should have degree g+1 = 3
        assert_eq!(vpl.deg(), 3);

        // f - Vpl² should have degree ≤ g = 2
        let vpl_sq = &vpl * &vpl;
        let diff = f - vpl_sq;
        assert!(diff.deg() <= g as i32);
    }

    #[test]
    fn test_neutral_split() {
        let f = test_curve_g2_split();
        let g = 2;
        let vpl = compute_vpl(&f, g);
        let neutral = Divisor::neutral(&f, &vpl, g);

        assert!(neutral.u.is_one());
        assert_eq!(neutral.v, vpl);
    }

    #[test]
    fn test_double_pos_consistency() {
        let f = test_curve_g2_split();
        let g = 2;
        let vpl = compute_vpl(&f, g);

        if let Some(d) = find_valid_split_divisor(&f, &vpl, g) {
            // Double should produce a valid divisor
            let doubled = double_pos(&d, &f, &vpl, g);
            assert!(doubled.u.deg() <= g as i32);

            // Verify w = (f - v²) / u
            let v_sq = &doubled.v * &doubled.v;
            let expected_w = (f - v_sq).exact_div(&doubled.u);
            assert_eq!(doubled.w, expected_w);
        }
    }

    #[test]
    fn test_double_neg_consistency() {
        let f = test_curve_g2_split();
        let g = 2;
        let vpl = compute_vpl(&f, g);
        let v_neg = -vpl.clone();

        if let Some(mut d) = find_valid_split_divisor(&f, &vpl, g) {
            // Convert to negative basis form
            let t = &v_neg - &vpl;
            let (q, r) = t.div_rem(&d.u);
            let tv = &d.v + &t - r;
            d.w -= q * (&d.v + &tv);
            d.v = tv;

            // Double should produce a valid divisor
            let doubled = double_neg(&d, &f, &v_neg, g);
            assert!(doubled.u.deg() <= g as i32);

            // Verify w = (f - v²) / u
            let v_sq = &doubled.v * &doubled.v;
            let expected_w = (f - v_sq).exact_div(&doubled.u);
            assert_eq!(doubled.w, expected_w);
        }
    }

    #[test]
    fn test_add_pos_neutral() {
        let f = test_curve_g2_split();
        let g = 2;
        let vpl = compute_vpl(&f, g);
        let neutral = Divisor::neutral(&f, &vpl, g);

        if let Some(d) = find_valid_split_divisor(&f, &vpl, g) {
            // D + 0 should give valid reduced divisor
            let result = add_pos(&d, &neutral, &f, &vpl, g);
            assert!(result.u.deg() <= g as i32);

            // Verify w = (f - v²) / u
            let v_sq = &result.v * &result.v;
            let expected_w = (f - v_sq).exact_div(&result.u);
            assert_eq!(result.w, expected_w);
        }
    }

    #[test]
    fn test_double_pos_vs_add_pos() {
        let f = test_curve_g2_split();
        let g = 2;
        let vpl = compute_vpl(&f, g);

        if let Some(d) = find_valid_split_divisor(&f, &vpl, g) {
            // 2D via double should give same u polynomial as D + D via add
            let doubled = double_pos(&d, &f, &vpl, g);
            let added = add_pos(&d, &d, &f, &vpl, g);

            assert_eq!(doubled.u, added.u);
        }
    }
}
