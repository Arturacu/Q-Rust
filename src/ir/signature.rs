//! Algebraic signatures used by commutation and simplification analyses.

use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PauliBasis {
    X,
    Y,
    Z,
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

    pub fn simplify(self) -> Self {
        if self.den == 0 {
            return self;
        }
        Self::from_wide(i128::from(self.num), i128::from(self.den)).unwrap_or(self)
    }

    pub fn reduce(self) -> Self {
        if self.den == 0 {
            return self;
        }
        // A denominator's magnitude and its doubled period may both exceed
        // i64::MAX. Saturation changes the period; wrapping negation changes
        // the sign of fractions involving i64::MIN.
        let mut numerator = i128::from(self.num);
        let mut denominator = i128::from(self.den);
        if denominator < 0 {
            numerator = -numerator;
            denominator = -denominator;
        }
        let period = 2 * denominator;
        numerator = numerator.rem_euclid(period);
        if numerator > denominator {
            numerator -= period;
        }
        Self::from_wide(numerator, denominator).unwrap_or(self)
    }

    /// Normalize in wide arithmetic, preserving an exact negative-denominator
    /// representation when the positive-denominator orientation cannot fit.
    fn from_wide(mut numerator: i128, mut denominator: i128) -> Option<Self> {
        let divisor = gcd_u128(numerator.unsigned_abs(), denominator.unsigned_abs()) as i128;
        numerator /= divisor;
        denominator /= divisor;
        if denominator < 0 {
            numerator = -numerator;
            denominator = -denominator;
        }
        for (num, den) in [(numerator, denominator), (-numerator, -denominator)] {
            if let (Ok(num), Ok(den)) = (i64::try_from(num), i64::try_from(den)) {
                return Some(Self { num, den });
            }
        }
        // Callers retain their exact original fraction if no representation
        // fits, rather than silently narrowing or changing its value.
        None
    }
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
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

    #[test]
    fn reduction_preserves_period_when_twice_denominator_exceeds_i64() {
        let fraction = SymbolicFraction::new(i64::MAX - 1, i64::MAX - 2).unwrap();
        assert_eq!(
            fraction.reduce(),
            SymbolicFraction {
                num: -i64::MAX + 3,
                den: i64::MAX - 2,
            }
        );
        let negative = SymbolicFraction::new(-i64::MAX + 1, i64::MAX - 2).unwrap();
        assert_eq!(
            negative.reduce(),
            SymbolicFraction {
                num: i64::MAX - 3,
                den: i64::MAX - 2,
            }
        );
    }

    #[test]
    fn minimum_integer_sign_normalization_preserves_exact_fraction() {
        for (num, den, expected_num, expected_den) in [
            (i64::MIN, -1, i64::MIN, -1),
            (1, i64::MIN, 1, i64::MIN),
            (-1, i64::MIN, -1, i64::MIN),
            (i64::MIN, i64::MIN, 1, 1),
            (i64::MIN, 2, i64::MIN / 2, 1),
            (2, i64::MIN, -1, 1_i64 << 62),
            (0, i64::MIN, 0, 1),
        ] {
            let simplified = SymbolicFraction::new(num, den).unwrap();
            assert_eq!(
                simplified,
                SymbolicFraction {
                    num: expected_num,
                    den: expected_den
                }
            );
            // Exact rational equality, independent of floating-point radians.
            assert_eq!(
                i128::from(num) * i128::from(simplified.den),
                i128::from(simplified.num) * i128::from(den)
            );
            assert_eq!(simplified.simplify(), simplified);
        }
    }

    #[test]
    fn modular_reduction_handles_minimum_numerators_and_denominators() {
        for (num, den, expected_num, expected_den) in [
            (i64::MIN, -1, 0, 1),
            (i64::MIN, 1, 0, 1),
            (i64::MIN, -3, 2, 3),
            (1, i64::MIN, 1, i64::MIN),
            (-1, i64::MIN, -1, i64::MIN),
            (i64::MIN, i64::MIN, 1, 1),
            (i64::MAX, i64::MIN, i64::MAX, i64::MIN),
        ] {
            let fraction = SymbolicFraction { num, den };
            let reduced = fraction.reduce();
            assert_eq!(
                reduced,
                SymbolicFraction {
                    num: expected_num,
                    den: expected_den
                }
            );
            assert_eq!(reduced.reduce(), reduced);
        }
        let invalid = SymbolicFraction {
            num: i64::MIN,
            den: 0,
        };
        assert_eq!(invalid.simplify(), invalid);
        assert_eq!(invalid.reduce(), invalid);
    }
}
