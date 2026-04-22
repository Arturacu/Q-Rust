// The foundational axes of the Pauli matrices.

/// The foundational axes of the Pauli matrices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PauliBasis {
    X,
    Y,
    Z,
}

/// Represents a fractional angle expressed as `numerator / denominator * π`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SymbolicFraction {
    pub num: i64,
    pub den: i64,
}

impl SymbolicFraction {
    /// Evaluates the fraction to a floating-point value in radians.
    pub fn to_radians(self) -> f64 {
        (self.num as f64 / self.den as f64) * std::f64::consts::PI
    }

    /// Reduces the fraction modulo 2π, mapping it into [-π, π].
    pub fn reduce(mut self) -> Self {
        if self.den == 0 {
            return self;
        }
        // Simplified canonicalization:
        let period = 2 * self.den;
        self.num %= period;
        if self.num > self.den {
            self.num -= period;
        } else if self.num < -self.den {
            self.num += period;
        }
        // Try to simplify based on powers of 2 (typical in quantum circuits)
        while self.num % 2 == 0 && self.den % 2 == 0 && self.den > 1 {
            self.num /= 2;
            self.den /= 2;
        }
        self
    }
}

/// Strongly-typed canonical rotation parameters for use in optimizations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SymbolicAngle {
    /// Exact rational multiple of π
    Rational(SymbolicFraction),
    /// Inexact f64 parameterized angle, requiring epsilon bounded matching
    Float(f64),
}

/// Structural hint for an operator's algebraic commutativity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommutationSignature {
    /// Operates diagonally across its entire support in a shared Pauli Tensor Basis.
    /// Example: RZ, Z, S, T are Z-diagonal. RX, X are X-diagonal.
    Diagonal(PauliBasis),

    /// Entangling gate exhibiting independent local diagonal symmetries per-wire.
    /// E.g., CX is Z-diagonal on the Control (local idx 0), and X-diagonal on Target (local idx 1).
    CompositeDiagonal(Vec<(usize, PauliBasis)>),

    /// Pure Clifford basis transformers (e.g. H) acting natively on Pauli frames.
    Clifford,

    /// Unresolved or non-Clifford parametric gates structurally blocking trivial topological commutativity.
    Generic,
}
