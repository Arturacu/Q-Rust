//! Quantum Shannon Decomposition dispatcher + numerical search for small
//! unitaries (QSearch-style stub).

use crate::ir::{Circuit, GateType, Operation};
use crate::transpiler::synthesis::kak::KakSynthesizer;
use crate::transpiler::synthesis::zyz::ZyzSynthesizer;
use crate::transpiler::synthesis::Synthesizer;
use nalgebra::DMatrix;
use num_complex::Complex;

/// Dispatcher to ZYZ (N=1) or KAK (N=2); returns `None` for N≥3.
#[derive(Debug, Clone, Copy)]
pub struct QsdSynthesizer;

impl Synthesizer for QsdSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, basis: &[GateType]) -> Option<Circuit> {
        let rows = unitary.nrows();
        if rows != unitary.ncols() || rows == 0 || !rows.is_power_of_two() {
            return None;
        }
        let n = (rows as f64).log2() as usize;
        match n {
            1 => ZyzSynthesizer.synthesize(unitary, basis),
            2 => KakSynthesizer.synthesize(unitary, basis),
            _ => None,
        }
    }
}

/// QSearch-style numerical synthesis via a gradient-free Nelder–Mead
/// optimizer over a fixed single-qubit ansatz (for 1-qubit) or two-qubit
/// ansatz (for 2-qubit). This is a *stub* implementation — correctness is
/// only guaranteed for 1-qubit inputs where it essentially reproduces ZYZ.
#[derive(Debug, Clone, Copy)]
pub struct QSearchSynthesizer {
    pub max_iter: usize,
    pub tol: f64,
}

impl Default for QSearchSynthesizer {
    fn default() -> Self {
        Self {
            max_iter: 500,
            tol: 1e-8,
        }
    }
}

impl Synthesizer for QSearchSynthesizer {
    fn synthesize(&self, unitary: &DMatrix<Complex<f64>>, basis: &[GateType]) -> Option<Circuit> {
        let rows = unitary.nrows();
        if rows != unitary.ncols() || !rows.is_power_of_two() {
            return None;
        }
        let n = (rows as f64).log2() as usize;
        match n {
            1 => qsearch_1q(unitary, self.max_iter, self.tol),
            2 => {
                // Fall back to KAK — QSearch over 2q is expensive; prefer
                // analytical decomposition when available.
                KakSynthesizer.synthesize(unitary, basis)
            }
            _ => None,
        }
    }
}

fn qsearch_1q(target: &DMatrix<Complex<f64>>, max_iter: usize, tol: f64) -> Option<Circuit> {
    // Ansatz: single U(theta, phi, lambda). Nelder-Mead over 3-dim.
    let cost = |params: &[f64]| -> f64 {
        let u = crate::transpiler::synthesis::zyz::u_to_matrix(params[0], params[1], params[2]);
        // Frobenius distance up to a global phase:
        // Use 1 - |tr(U† V)|/d
        let t = target;
        let mut acc = Complex::new(0.0, 0.0);
        for i in 0..2 {
            for j in 0..2 {
                acc += u[i][j].conj() * t[(i, j)];
            }
        }
        1.0 - (acc.norm() / 2.0)
    };

    // Initial simplex
    let mut simplex: Vec<Vec<f64>> = vec![
        vec![0.1, 0.0, 0.0],
        vec![1.2, 0.0, 0.0],
        vec![0.1, 1.3, 0.0],
        vec![0.1, 0.0, 1.4],
    ];
    let mut fs: Vec<f64> = simplex.iter().map(|p| cost(p)).collect();

    for _ in 0..max_iter {
        // Sort by cost ascending
        let mut idx: Vec<usize> = (0..simplex.len()).collect();
        idx.sort_by(|&a, &b| fs[a].partial_cmp(&fs[b]).unwrap());
        let simplex_sorted: Vec<_> = idx.iter().map(|&i| simplex[i].clone()).collect();
        let fs_sorted: Vec<_> = idx.iter().map(|&i| fs[i]).collect();
        simplex = simplex_sorted;
        fs = fs_sorted;

        if fs[0] < tol {
            break;
        }

        // Centroid of all but worst
        let n = simplex[0].len();
        let mut centroid = vec![0.0; n];
        for p in &simplex[..simplex.len() - 1] {
            for i in 0..n {
                centroid[i] += p[i];
            }
        }
        for c in &mut centroid {
            *c /= (simplex.len() - 1) as f64;
        }
        let worst = simplex.last().unwrap().clone();
        let f_worst = *fs.last().unwrap();
        let f_best = fs[0];
        let f_second_worst = fs[fs.len() - 2];

        // Reflection
        let alpha = 1.0;
        let xr: Vec<f64> = centroid
            .iter()
            .zip(&worst)
            .map(|(c, w)| c + alpha * (c - w))
            .collect();
        let fr = cost(&xr);

        if f_best <= fr && fr < f_second_worst {
            *simplex.last_mut().unwrap() = xr;
            *fs.last_mut().unwrap() = fr;
            continue;
        }

        // Expansion
        if fr < f_best {
            let gamma = 2.0;
            let xe: Vec<f64> = centroid
                .iter()
                .zip(&worst)
                .map(|(c, w)| c + gamma * (c - w))
                .collect();
            let fe = cost(&xe);
            if fe < fr {
                *simplex.last_mut().unwrap() = xe;
                *fs.last_mut().unwrap() = fe;
            } else {
                *simplex.last_mut().unwrap() = xr;
                *fs.last_mut().unwrap() = fr;
            }
            continue;
        }

        // Contraction
        let rho = 0.5;
        let xc: Vec<f64> = centroid
            .iter()
            .zip(&worst)
            .map(|(c, w)| c + rho * (w - c))
            .collect();
        let fc = cost(&xc);
        if fc < f_worst {
            *simplex.last_mut().unwrap() = xc;
            *fs.last_mut().unwrap() = fc;
            continue;
        }

        // Shrink
        let sigma = 0.5;
        let best = simplex[0].clone();
        for i in 1..simplex.len() {
            for j in 0..n {
                simplex[i][j] = best[j] + sigma * (simplex[i][j] - best[j]);
            }
            fs[i] = cost(&simplex[i]);
        }
    }

    let best_idx = (0..simplex.len()).min_by(|&a, &b| fs[a].partial_cmp(&fs[b]).unwrap())?;
    let params = simplex[best_idx].clone();
    if fs[best_idx] > 1e-3 {
        return None;
    }
    let mut c = Circuit::new(1, 0);
    c.add_op(Operation::Gate {
        name: GateType::U,
        qubits: vec![0],
        params,
    });
    Some(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qsd_returns_none_for_n3() {
        let id8 = DMatrix::<Complex<f64>>::identity(8, 8);
        assert!(QsdSynthesizer.synthesize(&id8, &[]).is_none());
    }

    #[test]
    fn test_qsearch_identity_1q() {
        let id = DMatrix::<Complex<f64>>::identity(2, 2);
        let res = QSearchSynthesizer::default().synthesize(&id, &[]);
        assert!(res.is_some());
    }
}
