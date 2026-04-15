//! Algebraic signatures used by commutation and simplification analyses.

use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PauliBasis {
    X, Y, Z,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SymbolicFraction {
    pub num: i64,
    pub den: i64,
}

impl SymbolicFraction {
    pub fn new(num: i64, den: i64) -> Option<Self> {
        if den == 0 {
            None
        } else {
            Some(SymbolicFraction { num, den }.simplify())
        }
    }

    #[inline]
    pub fn to_radians(self) -> f64 {
        (self.num as f64 / self.den as f64) * PI
    }

    pub fn simplify(mut self) -> Self {
        if self.den == 0 {
            return self;
        }
        if self.den < 0 {
            self.num = self.num.wrapping_neg();
            self.den = self.den.wrapping_neg();
        }
        let a = self.num.unsigned_abs();
        let b = self.den.unsigned_abs();
        let g = gcd_u64(a, b);
        if g > 1 {
            self.num /= g as i64;
            self.den /= g as i64;
        }
        self
    }

    pub fn reduce(self) -> Self {
        let simplified = self.simplify();
        if simplified.den == 0 {
            return simplified;
        }
        let period = 2i64.saturating_mul(simplified.den);
        let mut n = simplified.num.rem_euclid(period);
        if n > simplified.den {
            n -= period;
        }
        SymbolicFraction {
            num: n,
            den: simplified.den,
        }
        .simplify()
    }
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.max(1)
}

#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum SymbolicAngle {
    Rational(SymbolicFraction),
    Float(f64),
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CommutationSignature {
    Diagonal(PauliBasis),
    CompositeDiagonal(Vec<(usize, PauliBasis)>),
    Clifford,
    Generic,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fraction_simplify() {
        let f = SymbolicFraction::new(4, 8).unwrap();
        assert_eq!(f, SymbolicFraction { num: 1, den: 2 });
    }

    #[test]
    fn test_fraction_simplify_gcd() {
        let f = SymbolicFraction::new(6, 9).unwrap();
        assert_eq!(f, SymbolicFraction { num: 2, den: 3 });
    }

    #[test]
    fn test_fraction_negative_denominator() {
        let f = SymbolicFraction::new(3, -6).unwrap();
        assert_eq!(f, SymbolicFraction { num: -1, den: 2 });
    }

    #[test]
    fn test_fraction_zero_den_sentinel() {
        assert!(SymbolicFraction::new(1, 0).is_none());
    }
}