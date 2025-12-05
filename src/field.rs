//! Finite field arithmetic traits and implementations.
//!
//! Provides a generic `Field` trait and a simple `PrimeField` implementation
//! for fields of prime order.

use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Trait for finite field elements.
///
/// Implementors must provide field arithmetic operations.
pub trait Field:
    Sized
    + Clone
    + Copy
    + Debug
    + Default
    + PartialEq
    + Eq
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
    + Neg<Output = Self>
{
    /// The zero element (additive identity)
    fn zero() -> Self;

    /// The one element (multiplicative identity)
    fn one() -> Self;

    /// Check if this element is zero
    fn is_zero(&self) -> bool;

    /// Check if this element is one
    #[inline]
    fn is_one(&self) -> bool {
        *self == Self::one()
    }

    /// Compute the multiplicative inverse (panics if zero)
    fn inv(&self) -> Self;

    /// Compute self^exp using square-and-multiply
    #[inline]
    fn pow(&self, exp: u64) -> Self {
        if exp == 0 {
            return Self::one();
        }
        let mut base = *self;
        let mut result = Self::one();
        let mut e = exp;
        while e > 0 {
            if e & 1 == 1 {
                result *= base;
            }
            base *= base;
            e >>= 1;
        }
        result
    }

    /// Square this element
    #[inline]
    fn square(&self) -> Self {
        *self * *self
    }

    /// Double this element (add to itself)
    #[inline]
    fn double(&self) -> Self {
        *self + *self
    }

    /// Generate a random field element
    fn random<R: rand::Rng>(rng: &mut R) -> Self;
}

/// A simple prime field implementation using u64 arithmetic.
/// Suitable for small primes where p² fits in u128.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PrimeField<const P: u64> {
    value: u64,
}

impl<const P: u64> PrimeField<P> {
    /// Create a new field element, reducing modulo P
    #[inline]
    pub fn new(value: u64) -> Self {
        Self { value: value % P }
    }

    /// Create from a signed integer
    #[inline]
    pub fn from_i64(value: i64) -> Self {
        if value >= 0 {
            Self::new(value as u64)
        } else {
            Self::new((P as i64 + (value % P as i64)) as u64)
        }
    }

    /// Get the raw value
    #[inline]
    pub fn value(&self) -> u64 {
        self.value
    }

    /// Extended Euclidean algorithm for modular inverse
    fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
        if a == 0 {
            (b, 0, 1)
        } else {
            let (gcd, x, y) = Self::extended_gcd(b % a, a);
            (gcd, y - (b / a) * x, x)
        }
    }
}

impl<const P: u64> Field for PrimeField<P> {
    #[inline]
    fn zero() -> Self {
        Self { value: 0 }
    }

    #[inline]
    fn one() -> Self {
        Self { value: 1 }
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.value == 0
    }

    #[inline]
    fn inv(&self) -> Self {
        assert!(!self.is_zero(), "Cannot invert zero");
        let (_, x, _) = Self::extended_gcd(self.value as i128, P as i128);
        let inv = ((x % P as i128) + P as i128) % P as i128;
        Self { value: inv as u64 }
    }

    #[inline]
    fn random<R: rand::Rng>(rng: &mut R) -> Self {
        Self::new(rng.gen::<u64>() % P)
    }
}

impl<const P: u64> Add for PrimeField<P> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let sum = (self.value as u128 + rhs.value as u128) % P as u128;
        Self { value: sum as u64 }
    }
}

impl<const P: u64> AddAssign for PrimeField<P> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const P: u64> Sub for PrimeField<P> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let diff = if self.value >= rhs.value {
            self.value - rhs.value
        } else {
            P - (rhs.value - self.value)
        };
        Self { value: diff }
    }
}

impl<const P: u64> SubAssign for PrimeField<P> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const P: u64> Mul for PrimeField<P> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let prod = (self.value as u128 * rhs.value as u128) % P as u128;
        Self { value: prod as u64 }
    }
}

impl<const P: u64> MulAssign for PrimeField<P> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const P: u64> Div for PrimeField<P> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl<const P: u64> DivAssign for PrimeField<P> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<const P: u64> Neg for PrimeField<P> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        if self.value == 0 {
            self
        } else {
            Self {
                value: P - self.value,
            }
        }
    }
}

impl<const P: u64> From<u64> for PrimeField<P> {
    #[inline]
    fn from(value: u64) -> Self {
        Self::new(value)
    }
}

impl<const P: u64> From<i64> for PrimeField<P> {
    #[inline]
    fn from(value: i64) -> Self {
        Self::from_i64(value)
    }
}

impl<const P: u64> std::fmt::Display for PrimeField<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

/// Binary extension field GF(2^k).
///
/// Elements are represented as polynomials over GF(2) modulo an irreducible polynomial.
/// The representation uses a u64 bitfield: `a_0 + a_1*x + ... + a_{k-1}*x^{k-1}`
/// where bit i represents coefficient a_i.
///
/// Standard irreducible polynomials:
/// - GF(4):  x² + x + 1
/// - GF(8):  x³ + x + 1
/// - GF(16): x⁴ + x + 1
///
/// The generator element is α = x (bits = 2).
#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub struct BinaryExtField<const K: usize> {
    bits: u64,
}

impl<const K: usize> BinaryExtField<K> {
    /// Get the irreducible polynomial for this field.
    /// Returns the polynomial as a bit pattern (including the leading x^k term).
    ///
    /// These are minimal-weight irreducible polynomials over GF(2).
    /// Preference is given to trinomials (x^k + x^m + 1) when they exist,
    /// otherwise pentanomials are used.
    ///
    /// Sources:
    /// - NIST FIPS 186-4 (for cryptographic fields)
    /// - "Table of Low-Weight Binary Irreducible Polynomials" (HP Labs)
    /// - Lidl & Niederreiter, "Finite Fields"
    const fn irreducible() -> u64 {
        // Bit pattern represents: x^K + (lower terms)
        // e.g., 0b111 for K=2 means x² + x + 1
        match K {
            1 => 0b11,                         // x + 1 (trivial: GF(2))
            2 => 0b111,                        // x² + x + 1
            3 => 0b1011,                       // x³ + x + 1
            4 => 0b10011,                      // x⁴ + x + 1
            5 => 0b100101,                     // x⁵ + x² + 1
            6 => 0b1000011,                    // x⁶ + x + 1
            7 => 0b10000011,                   // x⁷ + x + 1
            8 => 0b100011011,                  // x⁸ + x⁴ + x³ + x + 1
            9 => 0b1000010001,                 // x⁹ + x⁴ + 1
            10 => 0b10000001001,               // x¹⁰ + x³ + 1
            11 => 0b100000000101,              // x¹¹ + x² + 1
            12 => 0b1000001010011,             // x¹² + x⁶ + x⁴ + x + 1
            13 => 0b10000000011011,            // x¹³ + x⁴ + x³ + x + 1
            14 => 0b100000000101011,           // x¹⁴ + x⁵ + x³ + x + 1
            15 => 0b1000000000000011,          // x¹⁵ + x + 1
            16 => 0b10001000000001011,         // x¹⁶ + x¹² + x³ + x + 1
            17 => 0b100000000000001001,        // x¹⁷ + x³ + 1
            18 => 0b1000000000010000001,       // x¹⁸ + x⁷ + 1
            19 => 0b10000000000000100111,      // x¹⁹ + x⁵ + x² + x + 1
            20 => 0b100000000000000001001,     // x²⁰ + x³ + 1
            21 => 0b1000000000000000000101,    // x²¹ + x² + 1
            22 => 0b10000000000000000000011,   // x²² + x + 1
            23 => 0b100000000000000000100001,  // x²³ + x⁵ + 1
            24 => 0b1000000000000000010000111, // x²⁴ + x⁷ + x² + x + 1
            // Beyond K=24, the u64 bit patterns get unwieldy but still work up to K=63
            // For larger fields, consider using u128 or a polynomial representation
            _ => panic!("Unsupported field extension degree (supported: K=1-24)"),
        }
    }

    /// Create a new field element from a bit pattern.
    #[inline]
    pub const fn new(bits: u64) -> Self {
        // Mask to K bits
        let mask = (1u64 << K) - 1;
        Self { bits: bits & mask }
    }

    /// Create the generator element α (represented as x, i.e., bits = 2).
    /// This corresponds to `FF.1` in Magma.
    #[inline]
    pub const fn gen() -> Self {
        Self { bits: 2 }
    }

    /// Get the raw bit representation.
    #[inline]
    pub const fn bits(&self) -> u64 {
        self.bits
    }

    /// Carry-less multiplication of two polynomials over GF(2).
    /// Returns the unreduced product (may have degree up to 2k-2).
    #[inline]
    fn clmul(a: u64, b: u64) -> u64 {
        let mut result = 0u64;
        let mut shifted_a = a;
        let mut b_temp = b;

        while b_temp != 0 {
            if b_temp & 1 != 0 {
                result ^= shifted_a;
            }
            shifted_a <<= 1;
            b_temp >>= 1;
        }
        result
    }

    /// Reduce a polynomial modulo the irreducible polynomial.
    #[inline]
    fn reduce(mut value: u64) -> u64 {
        let irred = Self::irreducible();
        // Reduce starting from the highest possible degree (2k-2) down to k
        for i in (K..2 * K - 1).rev() {
            if value & (1u64 << i) != 0 {
                // XOR with irreducible shifted to align with bit i
                value ^= irred << (i - K);
            }
        }
        value & ((1u64 << K) - 1)
    }

    /// Compute self^exp using square-and-multiply.
    #[inline]
    pub fn pow(&self, exp: u64) -> Self {
        if exp == 0 {
            return Self::new(1);
        }
        let mut base = *self;
        let mut result = Self::new(1);
        let mut e = exp;
        while e > 0 {
            if e & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            e >>= 1;
        }
        result
    }
}

impl<const K: usize> Field for BinaryExtField<K> {
    #[inline]
    fn zero() -> Self {
        Self { bits: 0 }
    }

    #[inline]
    fn one() -> Self {
        Self { bits: 1 }
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.bits == 0
    }

    #[inline]
    fn inv(&self) -> Self {
        assert!(!self.is_zero(), "Cannot invert zero");
        // Use Fermat's little theorem: a^(-1) = a^(2^k - 2)
        let exp = (1u64 << K) - 2;
        self.pow(exp)
    }

    #[inline]
    fn random<R: rand::Rng>(rng: &mut R) -> Self {
        let mask = (1u64 << K) - 1;
        Self::new(rng.gen::<u64>() & mask)
    }
}

impl<const K: usize> Add for BinaryExtField<K> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        // In GF(2^k), addition is XOR
        Self {
            bits: self.bits ^ rhs.bits,
        }
    }
}

impl<const K: usize> AddAssign for BinaryExtField<K> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.bits ^= rhs.bits;
    }
}

impl<const K: usize> Sub for BinaryExtField<K> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        // In GF(2^k), subtraction is the same as addition (XOR)
        Self {
            bits: self.bits ^ rhs.bits,
        }
    }
}

impl<const K: usize> SubAssign for BinaryExtField<K> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.bits ^= rhs.bits;
    }
}

impl<const K: usize> Mul for BinaryExtField<K> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let product = Self::clmul(self.bits, rhs.bits);
        Self {
            bits: Self::reduce(product),
        }
    }
}

impl<const K: usize> MulAssign for BinaryExtField<K> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const K: usize> Div for BinaryExtField<K> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl<const K: usize> DivAssign for BinaryExtField<K> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<const K: usize> Neg for BinaryExtField<K> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        // In characteristic 2, -a = a
        self
    }
}

impl<const K: usize> From<u64> for BinaryExtField<K> {
    #[inline]
    fn from(value: u64) -> Self {
        Self::new(value)
    }
}

impl<const K: usize> std::fmt::Display for BinaryExtField<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Display as α^n or polynomial form
        if self.bits == 0 {
            write!(f, "0")
        } else if self.bits == 1 {
            write!(f, "1")
        } else {
            // Display as sum of powers of α
            let mut first = true;
            for i in (0..K).rev() {
                if self.bits & (1u64 << i) != 0 {
                    if !first {
                        write!(f, "+")?;
                    }
                    first = false;
                    match i {
                        0 => write!(f, "1")?,
                        1 => write!(f, "α")?,
                        _ => write!(f, "α^{}", i)?,
                    }
                }
            }
            Ok(())
        }
    }
}

impl<const K: usize> std::fmt::Debug for BinaryExtField<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GF(2^{})[{:#x}]", K, self.bits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type F7 = PrimeField<7>;

    #[test]
    fn test_basic_arithmetic() {
        let a = F7::new(3);
        let b = F7::new(5);

        assert_eq!((a + b).value(), 1); // 3 + 5 = 8 ≡ 1 (mod 7)
        assert_eq!((a - b).value(), 5); // 3 - 5 = -2 ≡ 5 (mod 7)
        assert_eq!((a * b).value(), 1); // 3 * 5 = 15 ≡ 1 (mod 7)
    }

    #[test]
    fn test_inverse() {
        let a = F7::new(3);
        let inv_a = a.inv();
        assert_eq!((a * inv_a).value(), 1);

        let b = F7::new(5);
        let inv_b = b.inv();
        assert_eq!((b * inv_b).value(), 1);
    }

    #[test]
    fn test_division() {
        let a = F7::new(6);
        let b = F7::new(2);
        assert_eq!((a / b).value(), 3); // 6 / 2 = 3
    }

    #[test]
    fn test_power() {
        let a = F7::new(3);
        assert_eq!(a.pow(0).value(), 1);
        assert_eq!(a.pow(1).value(), 3);
        assert_eq!(a.pow(2).value(), 2); // 9 ≡ 2 (mod 7)
        assert_eq!(a.pow(6).value(), 1); // Fermat's little theorem
    }

    // =============================================================
    // BinaryExtField tests
    // =============================================================

    type GF4 = BinaryExtField<2>;
    type GF8 = BinaryExtField<3>;

    #[test]
    fn test_gf4_basic() {
        // GF(4) = {0, 1, α, α+1} where α² + α + 1 = 0
        let zero = GF4::zero();
        let one = GF4::one();
        let alpha = GF4::gen(); // α = x, bits = 2
        let alpha_plus_1 = GF4::new(3); // α + 1, bits = 3

        assert!(zero.is_zero());
        assert!(one.is_one());
        assert_eq!(alpha.bits(), 2);
        assert_eq!(alpha_plus_1.bits(), 3);
    }

    #[test]
    fn test_gf4_addition() {
        // In GF(4), addition is XOR
        let alpha = GF4::gen();
        let one = GF4::one();

        // α + 1 = α + 1 (bits: 10 ^ 01 = 11)
        let result = alpha + one;
        assert_eq!(result.bits(), 3);

        // α + α = 0
        let result = alpha + alpha;
        assert!(result.is_zero());

        // (α + 1) + 1 = α
        let alpha_plus_1 = GF4::new(3);
        let result = alpha_plus_1 + one;
        assert_eq!(result, alpha);
    }

    #[test]
    fn test_gf4_multiplication() {
        // In GF(4), α² + α + 1 = 0, so α² = α + 1
        let alpha = GF4::gen();

        // α * α = α² = α + 1
        let alpha_sq = alpha * alpha;
        assert_eq!(alpha_sq.bits(), 3); // α + 1 = bits 11 = 3

        // α * (α + 1) = α² + α = (α + 1) + α = 1
        let alpha_plus_1 = GF4::new(3);
        let result = alpha * alpha_plus_1;
        assert!(result.is_one());

        // 1 * α = α
        let one = GF4::one();
        assert_eq!(one * alpha, alpha);
    }

    #[test]
    fn test_gf4_inverse() {
        // In GF(4):
        // 1^(-1) = 1
        // α^(-1) = α + 1 (since α * (α+1) = 1)
        // (α+1)^(-1) = α

        let one = GF4::one();
        let alpha = GF4::gen();
        let alpha_plus_1 = GF4::new(3);

        assert_eq!(one.inv(), one);
        assert_eq!(alpha.inv(), alpha_plus_1);
        assert_eq!(alpha_plus_1.inv(), alpha);

        // Verify: a * a^(-1) = 1
        assert!((alpha * alpha.inv()).is_one());
        assert!((alpha_plus_1 * alpha_plus_1.inv()).is_one());
    }

    #[test]
    fn test_gf8_basic() {
        // GF(8) = GF(2^3), α³ + α + 1 = 0
        let alpha = GF8::gen();

        // α³ = α + 1
        let alpha_cubed = alpha * alpha * alpha;
        assert_eq!(alpha_cubed.bits(), 3); // α + 1
    }

    #[test]
    fn test_gf8_powers_of_alpha() {
        // In GF(8) with α³ + α + 1 = 0:
        // α^0 = 1
        // α^1 = α
        // α^2 = α²
        // α^3 = α + 1
        // α^4 = α² + α
        // α^5 = α² + α + 1
        // α^6 = α² + 1
        // α^7 = 1 (generator of multiplicative group)

        let alpha = GF8::gen();

        assert_eq!(alpha.pow(0).bits(), 1); // 1
        assert_eq!(alpha.pow(1).bits(), 2); // α
        assert_eq!(alpha.pow(2).bits(), 4); // α²
        assert_eq!(alpha.pow(3).bits(), 3); // α + 1
        assert_eq!(alpha.pow(4).bits(), 6); // α² + α
        assert_eq!(alpha.pow(5).bits(), 7); // α² + α + 1
        assert_eq!(alpha.pow(6).bits(), 5); // α² + 1
        assert_eq!(alpha.pow(7).bits(), 1); // 1 (order divides 2^3-1=7)
    }

    #[test]
    fn test_gf8_all_inverses() {
        // Verify that every non-zero element has an inverse
        for i in 1..8u64 {
            let a = GF8::new(i);
            let a_inv = a.inv();
            assert!((a * a_inv).is_one(), "Inverse failed for bits={}", i);
        }
    }

    #[test]
    fn test_gf4_all_inverses() {
        // Verify that every non-zero element has an inverse
        for i in 1..4u64 {
            let a = GF4::new(i);
            let a_inv = a.inv();
            assert!((a * a_inv).is_one(), "Inverse failed for bits={}", i);
        }
    }

    #[test]
    fn test_gf8_division() {
        let alpha = GF8::gen();
        let alpha_sq = alpha * alpha;

        // α² / α = α
        let result = alpha_sq / alpha;
        assert_eq!(result, alpha);

        // α / α = 1
        let result = alpha / alpha;
        assert!(result.is_one());
    }

    #[test]
    fn test_gf4_negation() {
        // In characteristic 2, -a = a
        let alpha = GF4::gen();
        assert_eq!(-alpha, alpha);

        let alpha_plus_1 = GF4::new(3);
        assert_eq!(-alpha_plus_1, alpha_plus_1);
    }

    #[test]
    fn test_gf8_subtraction() {
        // In characteristic 2, subtraction = addition
        let alpha = GF8::gen();
        let one = GF8::one();

        let sum = alpha + one;
        let diff = alpha - one;
        assert_eq!(sum, diff);
    }

    // =============================================================
    // Larger extension field tests
    // =============================================================

    type GF16 = BinaryExtField<4>;
    type GF256 = BinaryExtField<8>;
    type GF65536 = BinaryExtField<16>;

    #[test]
    fn test_gf16_basic() {
        // GF(16) = GF(2^4), irreducible: x⁴ + x + 1
        let alpha = GF16::gen();

        // α⁴ = α + 1 (from x⁴ + x + 1 = 0)
        let alpha_4 = alpha.pow(4);
        assert_eq!(alpha_4.bits(), 3); // α + 1 = bits 0011

        // Verify multiplicative group order: α^15 = 1
        assert!(alpha.pow(15).is_one());
    }

    #[test]
    fn test_gf16_all_inverses() {
        for i in 1..16u64 {
            let a = GF16::new(i);
            let a_inv = a.inv();
            assert!((a * a_inv).is_one(), "Inverse failed for bits={}", i);
        }
    }

    #[test]
    fn test_gf256_basic() {
        // GF(256) = GF(2^8), irreducible: x⁸ + x⁴ + x³ + x + 1
        let alpha = GF256::gen();

        // Verify multiplicative group order: α^255 = 1
        assert!(alpha.pow(255).is_one());

        // Verify α^256 = α
        assert_eq!(alpha.pow(256), alpha);
    }

    #[test]
    fn test_gf256_sample_inverses() {
        // Test a sampling of inverses in GF(256)
        for i in [1, 2, 17, 42, 128, 200, 254, 255].iter() {
            let a = GF256::new(*i);
            let a_inv = a.inv();
            assert!((a * a_inv).is_one(), "Inverse failed for bits={}", i);
        }
    }

    #[test]
    fn test_gf65536_basic() {
        // GF(65536) = GF(2^16)
        let alpha = GF65536::gen();

        // Verify multiplicative group order: α^(2^16-1) = α^65535 = 1
        assert!(alpha.pow(65535).is_one());
    }

    #[test]
    fn test_gf65536_sample_inverses() {
        // Test a sampling of inverses in GF(2^16)
        for i in [1, 2, 1000, 12345, 32768, 65534, 65535].iter() {
            let a = GF65536::new(*i);
            let a_inv = a.inv();
            assert!((a * a_inv).is_one(), "Inverse failed for bits={}", i);
        }
    }

    #[test]
    fn test_gf65536_arithmetic() {
        let alpha = GF65536::gen();
        let a = alpha.pow(1000);
        let b = alpha.pow(2000);

        // Test associativity: (a * b) * b^-1 = a
        let ab = a * b;
        let result = ab * b.inv();
        assert_eq!(result, a);

        // Test distribution: a * (b + 1) = a*b + a
        let one = GF65536::one();
        let left = a * (b + one);
        let right = a * b + a;
        assert_eq!(left, right);
    }
}
