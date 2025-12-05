//! Polynomial arithmetic over finite fields.
//!
//! Provides a `Poly` type for univariate polynomials and operations needed
//! for divisor arithmetic.

use crate::field::Field;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A univariate polynomial over a field F.
///
/// Coefficients are stored in ascending order: `coeffs[i]` is the coefficient of x^i.
/// The polynomial is kept normalized (no trailing zeros except for the zero polynomial).
#[derive(Clone, PartialEq, Eq)]
pub struct Poly<F: Field> {
    /// Coefficients in ascending order of degree
    pub coeffs: Vec<F>,
}

impl<F: Field> Poly<F> {
    /// Create the zero polynomial
    pub fn zero() -> Self {
        Self { coeffs: vec![] }
    }

    /// Create a constant polynomial
    pub fn constant(c: F) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self { coeffs: vec![c] }
        }
    }

    /// Create a monomial c * x^n
    pub fn monomial(c: F, n: usize) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            let mut coeffs = vec![F::zero(); n + 1];
            coeffs[n] = c;
            Self { coeffs }
        }
    }

    /// Create x^n (monic monomial)
    pub fn x_pow(n: usize) -> Self {
        Self::monomial(F::one(), n)
    }

    /// Create from coefficients (ascending order)
    pub fn from_coeffs(coeffs: Vec<F>) -> Self {
        let mut p = Self { coeffs };
        p.normalize();
        p
    }

    /// Normalize by removing trailing zeros
    fn normalize(&mut self) {
        while let Some(c) = self.coeffs.last() {
            if c.is_zero() {
                self.coeffs.pop();
            } else {
                break;
            }
        }
    }

    /// Check if this is the zero polynomial
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Get the degree (returns None for zero polynomial)
    pub fn degree(&self) -> Option<usize> {
        if self.is_zero() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Get the degree as i32 (returns -1 for zero polynomial, matching Magma convention)
    pub fn deg(&self) -> i32 {
        self.degree().map_or(-1, |d| d as i32)
    }

    /// Get coefficient of x^n (returns 0 if n > degree)
    pub fn coeff(&self, n: usize) -> F {
        if n < self.coeffs.len() {
            self.coeffs[n]
        } else {
            F::zero()
        }
    }

    /// Get the leading coefficient (returns 0 for zero polynomial)
    pub fn leading_coeff(&self) -> F {
        self.coeffs.last().copied().unwrap_or(F::zero())
    }

    /// Check if this polynomial is monic
    pub fn is_monic(&self) -> bool {
        self.leading_coeff().is_one()
    }

    /// Make this polynomial monic (divide by leading coefficient)
    pub fn make_monic(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let lc = self.leading_coeff();
        let lc_inv = lc.inv();
        Self::from_coeffs(self.coeffs.iter().map(|c| *c * lc_inv).collect())
    }

    /// Evaluate the polynomial at a point
    pub fn eval(&self, x: F) -> F {
        // Horner's method
        let mut result = F::zero();
        for c in self.coeffs.iter().rev() {
            result = result * x + *c;
        }
        result
    }

    /// Polynomial division with remainder: self = q * divisor + r
    /// Returns (quotient, remainder)
    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        assert!(!divisor.is_zero(), "Division by zero polynomial");

        if self.is_zero() {
            return (Self::zero(), Self::zero());
        }

        let divisor_deg = divisor.degree().unwrap();
        let self_deg = match self.degree() {
            Some(d) => d,
            None => return (Self::zero(), Self::zero()),
        };

        if self_deg < divisor_deg {
            return (Self::zero(), self.clone());
        }

        let lc_inv = divisor.leading_coeff().inv();
        let mut remainder = self.coeffs.clone();
        let mut quotient = vec![F::zero(); self_deg - divisor_deg + 1];

        for i in (divisor_deg..=self_deg).rev() {
            let q_coeff = remainder[i] * lc_inv;
            quotient[i - divisor_deg] = q_coeff;
            for j in 0..=divisor_deg {
                remainder[i - divisor_deg + j] -= q_coeff * divisor.coeffs[j];
            }
        }

        (Self::from_coeffs(quotient), Self::from_coeffs(remainder))
    }

    /// Polynomial modulo: self mod divisor
    pub fn rem(&self, divisor: &Self) -> Self {
        self.div_rem(divisor).1
    }

    /// Exact division (panics if there's a remainder)
    pub fn exact_div(&self, divisor: &Self) -> Self {
        let (q, r) = self.div_rem(divisor);
        assert!(r.is_zero(), "Exact division failed: non-zero remainder");
        q
    }

    /// Extended GCD: returns (gcd, a, b) such that gcd = a*self + b*other
    ///
    /// The GCD is made monic (leading coefficient = 1).
    /// This matches Magma's XGCD function.
    pub fn xgcd(&self, other: &Self) -> (Self, Self, Self) {
        if self.is_zero() {
            if other.is_zero() {
                return (Self::zero(), Self::zero(), Self::zero());
            }
            let lc = other.leading_coeff();
            let lc_inv = lc.inv();
            return (other.clone() * lc_inv, Self::zero(), Self::constant(lc_inv));
        }
        if other.is_zero() {
            let lc = self.leading_coeff();
            let lc_inv = lc.inv();
            return (self.clone() * lc_inv, Self::constant(lc_inv), Self::zero());
        }

        // Extended Euclidean algorithm
        let mut r0 = self.clone();
        let mut r1 = other.clone();
        let mut a0 = Self::constant(F::one());
        let mut a1 = Self::zero();
        let mut b0 = Self::zero();
        let mut b1 = Self::constant(F::one());

        while !r1.is_zero() {
            let (q, r) = r0.div_rem(&r1);

            let new_a = a0.clone() - q.clone() * a1.clone();
            let new_b = b0.clone() - q * b1.clone();

            r0 = r1;
            r1 = r;
            a0 = a1;
            a1 = new_a;
            b0 = b1;
            b1 = new_b;
        }

        // Make GCD monic
        if !r0.is_zero() {
            let lc = r0.leading_coeff();
            let lc_inv = lc.inv();
            r0 = r0 * lc_inv;
            a0 = a0 * lc_inv;
            b0 = b0 * lc_inv;
        }

        (r0, a0, b0)
    }

    /// Check if this polynomial is one (the constant 1)
    pub fn is_one(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0].is_one()
    }

    /// Scalar division (divide all coefficients by a scalar)
    pub fn scalar_div(&self, scalar: F) -> Self {
        let inv = scalar.inv();
        Self::from_coeffs(self.coeffs.iter().map(|c| *c * inv).collect())
    }
}

impl<F: Field> Default for Poly<F> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<F: Field> Add for Poly<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = Vec::with_capacity(max_len);
        for i in 0..max_len {
            coeffs.push(self.coeff(i) + rhs.coeff(i));
        }
        Self::from_coeffs(coeffs)
    }
}

impl<F: Field> Add for &Poly<F> {
    type Output = Poly<F>;

    fn add(self, rhs: Self) -> Poly<F> {
        self.clone() + rhs.clone()
    }
}

impl<F: Field> AddAssign for Poly<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl<F: Field> Sub for Poly<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let max_len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = Vec::with_capacity(max_len);
        for i in 0..max_len {
            coeffs.push(self.coeff(i) - rhs.coeff(i));
        }
        Self::from_coeffs(coeffs)
    }
}

impl<F: Field> Sub for &Poly<F> {
    type Output = Poly<F>;

    fn sub(self, rhs: Self) -> Poly<F> {
        self.clone() - rhs.clone()
    }
}

impl<F: Field> SubAssign for Poly<F> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl<F: Field> Neg for Poly<F> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::from_coeffs(self.coeffs.iter().map(|c| -*c).collect())
    }
}

impl<F: Field> Mul for Poly<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        let new_len = self.coeffs.len() + rhs.coeffs.len() - 1;
        let mut coeffs = vec![F::zero(); new_len];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                coeffs[i + j] += *a * *b;
            }
        }
        Self::from_coeffs(coeffs)
    }
}

impl<F: Field> Mul for &Poly<F> {
    type Output = Poly<F>;

    fn mul(self, rhs: Self) -> Poly<F> {
        self.clone() * rhs.clone()
    }
}

impl<F: Field> MulAssign for Poly<F> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

/// Scalar multiplication
impl<F: Field> Mul<F> for Poly<F> {
    type Output = Self;

    fn mul(self, scalar: F) -> Self {
        if scalar.is_zero() {
            Self::zero()
        } else {
            Self::from_coeffs(self.coeffs.iter().map(|c| *c * scalar).collect())
        }
    }
}

impl<F: Field + fmt::Display> fmt::Debug for Poly<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut first = true;
        for (i, c) in self.coeffs.iter().enumerate().rev() {
            if c.is_zero() {
                continue;
            }
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            match i {
                0 => write!(f, "{:?}", c)?,
                1 => {
                    if c.is_one() {
                        write!(f, "x")?
                    } else {
                        write!(f, "{:?}*x", c)?
                    }
                }
                _ => {
                    if c.is_one() {
                        write!(f, "x^{}", i)?
                    } else {
                        write!(f, "{:?}*x^{}", c, i)?
                    }
                }
            }
        }
        Ok(())
    }
}

impl<F: Field + fmt::Display> fmt::Display for Poly<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::PrimeField;

    type F7 = PrimeField<7>;

    fn poly(coeffs: &[u64]) -> Poly<F7> {
        Poly::from_coeffs(coeffs.iter().map(|&c| F7::new(c)).collect())
    }

    #[test]
    fn test_basic() {
        let p = poly(&[1, 2, 3]); // 3x² + 2x + 1
        assert_eq!(p.degree(), Some(2));
        assert_eq!(p.coeff(0).value(), 1);
        assert_eq!(p.coeff(1).value(), 2);
        assert_eq!(p.coeff(2).value(), 3);
        assert_eq!(p.coeff(3).value(), 0);
    }

    #[test]
    fn test_add_sub() {
        let p1 = poly(&[1, 2, 3]); // 3x² + 2x + 1
        let p2 = poly(&[4, 5]); // 5x + 4

        let sum = p1.clone() + p2.clone();
        assert_eq!(sum.coeff(0).value(), 5); // 1 + 4
        assert_eq!(sum.coeff(1).value(), 0); // 2 + 5 = 7 ≡ 0
        assert_eq!(sum.coeff(2).value(), 3);

        let diff = p1 - p2;
        assert_eq!(diff.coeff(0).value(), 4); // 1 - 4 = -3 ≡ 4
        assert_eq!(diff.coeff(1).value(), 4); // 2 - 5 = -3 ≡ 4
    }

    #[test]
    fn test_mul() {
        let p1 = poly(&[1, 1]); // x + 1
        let p2 = poly(&[6, 1]); // x + 6 = x - 1 (mod 7)

        let prod = p1 * p2; // (x+1)(x-1) = x² - 1
        assert_eq!(prod.coeff(0).value(), 6); // -1 ≡ 6
        assert_eq!(prod.coeff(1).value(), 0);
        assert_eq!(prod.coeff(2).value(), 1);
    }

    #[test]
    fn test_div_rem() {
        // x² + 2x + 1 divided by x + 1
        let p1 = poly(&[1, 2, 1]); // x² + 2x + 1 = (x+1)²
        let p2 = poly(&[1, 1]); // x + 1

        let (q, r) = p1.div_rem(&p2);
        assert_eq!(q.coeff(0).value(), 1);
        assert_eq!(q.coeff(1).value(), 1); // q = x + 1
        assert!(r.is_zero()); // remainder is 0
    }

    #[test]
    fn test_eval() {
        let p = poly(&[1, 2, 3]); // 3x² + 2x + 1
        let x = F7::new(2);
        // 3*4 + 2*2 + 1 = 12 + 4 + 1 = 17 ≡ 3 (mod 7)
        assert_eq!(p.eval(x).value(), 3);
    }

    #[test]
    fn test_xgcd_coprime() {
        // GCD of coprime polynomials should be 1
        let p1 = poly(&[1, 1]); // x + 1
        let p2 = poly(&[2, 1]); // x + 2

        let (gcd, a, b) = p1.xgcd(&p2);

        // GCD should be 1 (monic constant)
        assert!(gcd.is_one());

        // Verify: gcd = a*p1 + b*p2
        let check = a * p1 + b * p2;
        assert!(check.is_one());
    }

    #[test]
    fn test_xgcd_common_factor() {
        // (x+1)² and (x+1)(x+2) have GCD = x+1
        let p1 = poly(&[1, 2, 1]); // x² + 2x + 1 = (x+1)²
        let p2 = poly(&[2, 3, 1]); // x² + 3x + 2 = (x+1)(x+2)

        let (gcd, a, b) = p1.xgcd(&p2);

        // GCD should be x + 1 (monic)
        assert_eq!(gcd.degree(), Some(1));
        assert!(gcd.is_monic());
        assert_eq!(gcd.coeff(0).value(), 1); // x + 1

        // Verify: gcd = a*p1 + b*p2
        let check = a * p1 + b * p2;
        assert_eq!(check, gcd);
    }

    #[test]
    fn test_xgcd_zero() {
        let p1 = poly(&[1, 1]); // x + 1
        let zero = Poly::<F7>::zero();

        let (gcd, a, _b) = p1.xgcd(&zero);

        // GCD(p, 0) = p (made monic)
        assert!(gcd.is_monic());
        assert_eq!(gcd.degree(), Some(1));

        // a * p1 = gcd
        let check = a * p1;
        assert_eq!(check, gcd);
    }

    #[test]
    fn test_is_one() {
        let one = Poly::constant(F7::one());
        let two = Poly::constant(F7::new(2));
        let x = poly(&[0, 1]);

        assert!(one.is_one());
        assert!(!two.is_one());
        assert!(!x.is_one());
    }
}
