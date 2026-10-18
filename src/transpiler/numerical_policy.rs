//! Shared numerical acceptance policy for optimization replacements.
//!
//! A local phase-aligned Frobenius residual bounds the local operator norm.
//! For unitary replacements the sum bounds composition error in exact
//! arithmetic, including when the local operation is embedded on more wires.
//! We additionally charge a floating-point guard per evaluated operation.
//! This is a conservative engineering check, not an interval-arithmetic proof
//! of transcendental functions, every analytic decomposition, or the verifier.
//! Never interpret this ledger as a certificate of circuit equivalence.

use super::property_set::PropertySet;
use crate::ir::{GateDefinition, Operation};
use nalgebra::DMatrix;
use num_complex::Complex64;

pub const CIRCUIT_ERROR_ALLOWANCE: f64 = 1e-9;
pub(crate) const LOCAL_REWRITE_LIMIT: f64 = 1e-12;
const KEY: &str = "numerical_budget";

/// Cumulative numerical replacements across passes and both pipeline stages.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NumericalBudget {
    pub accepted: usize,
    pub rejected: usize,
    pub residual_sum: f64,
    /// Residuals plus a guard for floating-point matrix work.
    pub charged_error: f64,
}

pub(crate) fn transfer(source: &mut PropertySet, destination: &mut PropertySet) {
    if let Some(budget) = source.remove::<NumericalBudget>(KEY) {
        destination.insert(KEY, *budget);
    }
}

pub(crate) fn snapshot(properties: &PropertySet) -> NumericalBudget {
    properties
        .get::<NumericalBudget>(KEY)
        .cloned()
        .unwrap_or_default()
}

pub(crate) fn reserve(
    properties: &mut PropertySet,
    residual: f64,
    local_limit: f64,
    evaluated_operations: usize,
) -> bool {
    let mut budget = snapshot(properties);
    let guard = 32.0 * f64::EPSILON * evaluated_operations.max(1) as f64;
    let charge = residual + guard;
    let next = budget.charged_error + charge;
    let accept = residual.is_finite()
        && residual >= 0.0
        && residual <= local_limit
        && next.is_finite()
        && next <= CIRCUIT_ERROR_ALLOWANCE;
    if accept {
        budget.accepted += 1;
        budget.residual_sum += residual;
        budget.charged_error = next;
    } else {
        budget.rejected += 1;
    }
    properties.insert(KEY, budget);
    accept
}

pub(crate) fn phase_residual(
    original: &DMatrix<Complex64>,
    candidate: &DMatrix<Complex64>,
) -> Option<f64> {
    if original.shape() != candidate.shape()
        || original
            .iter()
            .chain(candidate.iter())
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return None;
    }
    let overlap: Complex64 = original
        .iter()
        .zip(candidate.iter())
        .map(|(a, b)| a.conj() * b)
        .sum();
    let magnitude = overlap.norm();
    if !magnitude.is_finite() || magnitude <= f64::EPSILON {
        return None;
    }
    let phase = overlap / magnitude;
    let residual = original
        .iter()
        .zip(candidate.iter())
        .map(|(a, b)| (b - a * phase).norm_sqr())
        .sum::<f64>()
        .sqrt();
    residual.is_finite().then_some(residual)
}

/// Only called for supported named gates on the same ordered one/two wires.
/// The independent regression oracle does not reuse these gate definitions.
pub(crate) fn accept_rewrite(
    properties: &mut PropertySet,
    originals: &[&Operation],
    replacement: Option<&Operation>,
) -> bool {
    let Some(first) = originals.first() else {
        return false;
    };
    let wires = first.qubits();
    if !(1..=2).contains(&wires.len()) {
        return false;
    }
    let dimension = 1 << wires.len();
    let matrix = |op: &Operation| {
        if let Operation::Gate {
            name,
            params,
            qubits,
        } = op
        {
            if qubits == wires && params.iter().all(|p| p.is_finite()) {
                let value = name.unitary(params);
                return (value.shape() == (dimension, dimension)).then_some(value);
            }
        }
        None
    };
    let mut original = DMatrix::identity(dimension, dimension);
    for op in originals {
        let Some(value) = matrix(op) else {
            return false;
        };
        original = value * original;
    }
    let candidate = match replacement {
        Some(op) => {
            let Some(value) = matrix(op) else {
                return false;
            };
            value
        }
        None => DMatrix::identity(dimension, dimension),
    };
    let residual = phase_residual(&original, &candidate).unwrap_or(f64::INFINITY);
    reserve(
        properties,
        residual,
        LOCAL_REWRITE_LIMIT,
        originals.len() + 1,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allowance_is_shared_and_nonfinite_or_exhausted_requests_fail_closed() {
        let mut properties = PropertySet::new();
        for _ in 0..9 {
            assert!(reserve(&mut properties, 1e-10, 1e-10, 1));
        }
        assert!(!reserve(&mut properties, 1e-10, 1e-10, 1));
        for bad in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(!reserve(&mut properties, bad, 1e-10, 1));
        }
        let mut lower = PropertySet::new();
        let before = snapshot(&properties);
        transfer(&mut properties, &mut lower);
        assert_eq!(snapshot(&lower), before);
        assert!(!reserve(&mut lower, 1e-10, 1e-10, 1));
    }
}
