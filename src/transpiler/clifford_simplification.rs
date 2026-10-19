//! Bounded, exact cross-block Clifford simplification.
//!
//! A CX is the controlled-Pauli involution
//! `C(P,Q) = (I + P + Q - PQ)/2`, initially `P = Z_control`,
//! `Q = X_target`. Conjugating through a local Clifford changes only these
//! signed axes. A third-wire interaction can be crossed only if it commutes
//! with *both* axes; thus its action on the transported interaction is exactly
//! the identity. Both endpoint interactions may move toward an interior cut:
//! for a middle segment `M = R L`, the original `B R L A` becomes
//! `R (R† B R) (L A L†) L`. Each crossed gate is checked independently.
//! Trying the reachable cuts, rather than only the second endpoint, preserves
//! parallelism that a legal but poorly placed rewrite could destroy.
//!
//! Synthesis uses integer signed-Pauli arithmetic and the 24 one-qubit
//! Cliffords. Exhausting the 24 x 24 right local frames determines whether
//! `U = L CX R` with local `L,R`; no floating-point gate recognition occurs.
//! Optional SWAP factorizations preserve the actual unitary here and are
//! permitted only for callers that subsequently absorb output permutations.
//! Every accepted rewrite strictly reduces CX cost and does not increase
//! total or CX depth, after that same permutation absorption when enabled.
//! Two bounded scans and at most 96 whole-circuit candidate evaluations avoid
//! quadratic folding/depth work on large inputs. Exhaustion retains the valid
//! improvements already accepted; this is deliberately not a fixed-point or
//! globally optimal Clifford optimizer.

use super::pass::Pass;
use super::property_set::PropertySet;
use crate::ir::{Circuit, GateType, Operation};
use std::collections::{HashMap, VecDeque};
use std::sync::{Mutex, OnceLock};

const LOOKAHEAD: usize = 96;
const SWEEPS: usize = 2;
const MAX_CANDIDATE_EVALUATIONS: usize = 96;

/// Only enable SWAP replacements when a subsequent pass absorbs their output
/// permutations. The default retains the exact fixed output-wire ordering.
#[derive(Clone, Copy, Debug, Default)]
pub struct CliffordSimplificationPass {
    pub allow_swap_replacements: bool,
}

/// Cumulative diagnostics under `clifford_simplification_stats`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CliffordSimplificationStats {
    pub candidates: usize,
    pub rewrites: usize,
    pub cx_removed: usize,
    pub depth_rejections: usize,
    /// Reachable insertion cuts considered, before synthesis and cost checks.
    pub insertion_cuts: usize,
    /// Accepted rewrites whose insertion lies strictly between both endpoints.
    pub interior_rewrites: usize,
    pub candidate_evaluations: usize,
    pub budget_exhausted: bool,
}

/// `i^phase X^x Z^z`; multiplication includes the exact ZX sign.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Pauli {
    x: u8,
    z: u8,
    phase: u8,
}

impl Pauli {
    const I: Self = Self {
        x: 0,
        z: 0,
        phase: 0,
    };
    fn hermitian(x: u8, z: u8) -> Self {
        Self {
            x,
            z,
            phase: (x & z).count_ones() as u8 % 4,
        }
    }
    fn mul(self, other: Self) -> Self {
        Self {
            x: self.x ^ other.x,
            z: self.z ^ other.z,
            phase: (self.phase + other.phase + 2 * (self.z & other.x).count_ones() as u8) % 4,
        }
    }
    fn neg(mut self) -> Self {
        self.phase = (self.phase + 2) % 4;
        self
    }
    fn anticommutes(self, other: Self) -> bool {
        ((self.x & other.z).count_ones() + (self.z & other.x).count_ones()) % 2 != 0
    }
}

/// Images of X0, Z0, X1, Z1. Equal tables mean equal unitaries up to a
/// *single* global phase, because these generators span the matrix algebra.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Clifford([Pauli; 4]);

impl Clifford {
    fn identity() -> Self {
        Self([
            Pauli::hermitian(1, 0),
            Pauli::hermitian(0, 1),
            Pauli::hermitian(2, 0),
            Pauli::hermitian(0, 2),
        ])
    }
    fn apply(self, p: Pauli) -> Pauli {
        let mut out = Pauli {
            phase: p.phase,
            ..Pauli::I
        };
        for q in 0..2 {
            if p.x & (1 << q) != 0 {
                out = out.mul(self.0[2 * q]);
            }
        }
        for q in 0..2 {
            if p.z & (1 << q) != 0 {
                out = out.mul(self.0[2 * q + 1]);
            }
        }
        out
    }
    /// Matrix composition `self * right`.
    fn compose(self, right: Self) -> Self {
        Self(right.0.map(|p| self.apply(p)))
    }
    fn inverse(self) -> Self {
        let mut out = Self::identity();
        for (index, target) in Self::identity().0.into_iter().enumerate() {
            for x in 0..4 {
                for z in 0..4 {
                    let p = Pauli::hermitian(x, z);
                    let image = self.apply(p);
                    if image.x == target.x && image.z == target.z {
                        out.0[index] = if image.phase == target.phase {
                            p
                        } else {
                            p.neg()
                        };
                    }
                }
            }
        }
        out
    }
    fn controlled(p: Pauli, q: Pauli) -> Self {
        debug_assert!(!p.anticommutes(q));
        Self(Self::identity().0.map(|r| {
            let a = r.anticommutes(p);
            let b = r.anticommutes(q);
            let mut image = r;
            if b {
                image = image.mul(p);
            }
            if a {
                image = image.mul(q);
            }
            if a && b {
                image = image.neg();
            }
            image
        }))
    }
    fn swap() -> Self {
        let id = Self::identity();
        Self([id.0[2], id.0[3], id.0[0], id.0[1]])
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SmallGate {
    H(usize),
    S(usize),
    Sdg(usize),
    X(usize),
    Y(usize),
    Z(usize),
    SX(usize),
    CX(usize),
    Swap,
}

impl SmallGate {
    fn table(self) -> Clifford {
        let mut c = Clifford::identity();
        let (q, x, z) = match self {
            Self::H(q) => (q, Pauli::hermitian(0, 1 << q), Pauli::hermitian(1 << q, 0)),
            Self::S(q) => (
                q,
                Pauli::hermitian(1 << q, 1 << q),
                Pauli::hermitian(0, 1 << q),
            ),
            Self::Sdg(q) => (
                q,
                Pauli::hermitian(1 << q, 1 << q).neg(),
                Pauli::hermitian(0, 1 << q),
            ),
            Self::X(q) => (
                q,
                Pauli::hermitian(1 << q, 0),
                Pauli::hermitian(0, 1 << q).neg(),
            ),
            Self::Y(q) => (
                q,
                Pauli::hermitian(1 << q, 0).neg(),
                Pauli::hermitian(0, 1 << q).neg(),
            ),
            Self::Z(q) => (
                q,
                Pauli::hermitian(1 << q, 0).neg(),
                Pauli::hermitian(0, 1 << q),
            ),
            Self::SX(q) => (
                q,
                Pauli::hermitian(1 << q, 0),
                Pauli::hermitian(1 << q, 1 << q).neg(),
            ),
            Self::CX(q) => {
                return Clifford::controlled(
                    Pauli::hermitian(0, 1 << q),
                    Pauli::hermitian(1 << (1 - q), 0),
                )
            }
            Self::Swap => return Clifford::swap(),
        };
        c.0[2 * q] = x;
        c.0[2 * q + 1] = z;
        c
    }
    fn inverse_table(self) -> Clifford {
        match self {
            Self::S(q) => Self::Sdg(q).table(),
            Self::Sdg(q) => Self::S(q).table(),
            // SX² = X, hence SX† = X SX up to the irrelevant shared phase.
            Self::SX(q) => Self::X(q).table().compose(self.table()),
            _ => self.table(),
        }
    }
    fn on_wire(self, q: usize) -> Self {
        match self {
            Self::H(_) => Self::H(q),
            Self::S(_) => Self::S(q),
            Self::Sdg(_) => Self::Sdg(q),
            Self::X(_) => Self::X(q),
            Self::Y(_) => Self::Y(q),
            Self::Z(_) => Self::Z(q),
            Self::SX(_) => Self::SX(q),
            _ => unreachable!("only local words are relocated"),
        }
    }
    fn operation(self, wires: [usize; 2]) -> Operation {
        let (name, qubits) = match self {
            Self::H(q) => (GateType::H, vec![wires[q]]),
            Self::S(q) => (GateType::S, vec![wires[q]]),
            Self::Sdg(q) => (GateType::Sdg, vec![wires[q]]),
            Self::X(q) => (GateType::X, vec![wires[q]]),
            Self::Y(q) => (GateType::Y, vec![wires[q]]),
            Self::Z(q) => (GateType::Z, vec![wires[q]]),
            Self::SX(q) => (GateType::SX, vec![wires[q]]),
            Self::CX(q) => (GateType::CX, vec![wires[q], wires[1 - q]]),
            Self::Swap => (GateType::SWAP, wires.to_vec()),
        };
        Operation::Gate {
            name,
            qubits,
            params: vec![],
        }
    }
}

fn local_gate(name: &GateType, q: usize) -> Option<SmallGate> {
    Some(match name {
        GateType::H => SmallGate::H(q),
        GateType::S => SmallGate::S(q),
        GateType::Sdg => SmallGate::Sdg(q),
        GateType::X => SmallGate::X(q),
        GateType::Y => SmallGate::Y(q),
        GateType::Z => SmallGate::Z(q),
        GateType::SX => SmallGate::SX(q),
        _ => return None,
    })
}

struct LocalFrames {
    words: HashMap<Clifford, Vec<SmallGate>>,
    pairs: Vec<(Clifford, Vec<SmallGate>)>,
}

fn frames() -> &'static LocalFrames {
    static FRAMES: OnceLock<LocalFrames> = OnceLock::new();
    FRAMES.get_or_init(|| {
        let generators = [
            SmallGate::H(0),
            SmallGate::S(0),
            SmallGate::Sdg(0),
            SmallGate::X(0),
            SmallGate::Y(0),
            SmallGate::Z(0),
            SmallGate::SX(0),
        ];
        let mut queue = VecDeque::from([(Clifford::identity(), vec![])]);
        let mut words = HashMap::from([(Clifford::identity(), vec![])]);
        let mut ordered = Vec::new();
        while let Some((table, word)) = queue.pop_front() {
            ordered.push((table, word.clone()));
            for gate in generators {
                let next = gate.table().compose(table);
                if let std::collections::hash_map::Entry::Vacant(entry) = words.entry(next) {
                    let mut sequence = word.clone();
                    sequence.push(gate);
                    entry.insert(sequence.clone());
                    queue.push_back((next, sequence));
                }
            }
        }
        assert_eq!(ordered.len(), 24);
        let mut pairs = Vec::with_capacity(576);
        for (left, a) in &ordered {
            for (_, b) in &ordered {
                let mut table = *left;
                let mut word = a.clone();
                for gate in b {
                    let relocated = gate.on_wire(1);
                    table = relocated.table().compose(table);
                    word.push(relocated);
                }
                pairs.push((table.inverse(), word));
            }
        }
        LocalFrames { words, pairs }
    })
}

fn local_word(table: Clifford) -> Option<Vec<SmallGate>> {
    let mut output = Vec::new();
    for q in 0..2 {
        let mut local = Clifford::identity();
        for k in 0..2 {
            let p = table.0[2 * q + k];
            if (p.x | p.z) & !(1 << q) != 0 {
                return None;
            }
            local.0[k] = Pauli {
                x: p.x >> q,
                z: p.z >> q,
                phase: p.phase,
            };
        }
        output.extend(frames().words.get(&local)?.iter().map(|g| g.on_wire(q)));
    }
    Some(output)
}

fn word_cost(word: &[SmallGate]) -> (usize, usize) {
    let mut depth = [0; 2];
    for gate in word {
        match *gate {
            SmallGate::CX(_) | SmallGate::Swap => {
                let d = depth[0].max(depth[1]) + 1;
                depth = [d, d];
            }
            SmallGate::H(q)
            | SmallGate::S(q)
            | SmallGate::Sdg(q)
            | SmallGate::X(q)
            | SmallGate::Y(q)
            | SmallGate::Z(q)
            | SmallGate::SX(q) => depth[q] += 1,
        }
    }
    (depth[0].max(depth[1]), word.len())
}

/// Exact bounded search, cached by the encountered controlled-Pauli product.
/// Each endpoint has 36 signed local-axis choices. At most 36² products and
/// their SWAP factorizations are encountered; the full Clifford group is never
/// generated. The lookahead and candidate budget bound calls per invocation.
fn synthesize(table: Clifford) -> Option<Vec<SmallGate>> {
    type Cache = Mutex<HashMap<Clifford, Option<Vec<SmallGate>>>>;
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(result) = cache
        .lock()
        .expect("Clifford cache poisoned")
        .get(&table)
        .cloned()
    {
        return result;
    }
    let result = if let Some(word) = local_word(table) {
        Some(word)
    } else {
        let mut best: Option<Vec<SmallGate>> = None;
        for (right_inverse, right_word) in &frames().pairs {
            for direction in 0..2 {
                let cx = SmallGate::CX(direction);
                let left = table.compose(*right_inverse).compose(cx.table());
                if let Some(left_word) = local_word(left) {
                    let mut candidate = right_word.clone();
                    candidate.push(cx);
                    candidate.extend(left_word);
                    if best
                        .as_ref()
                        .map_or(true, |old| word_cost(&candidate) < word_cost(old))
                    {
                        best = Some(candidate);
                    }
                }
            }
        }
        best
    };
    cache
        .lock()
        .expect("Clifford cache poisoned")
        .insert(table, result.clone());
    result
}

/// Remove redundant local frames exposed by deleting an interaction. Only
/// disjoint operations are reordered; all unsupported and dynamic operations
/// flush the frames on their wires (barriers/dynamic operations flush all).
fn fold_local<'a>(
    circuit: &Circuit,
    operations: impl IntoIterator<Item = &'a Operation>,
) -> Circuit {
    let mut out = Circuit::new(circuit.num_qubits, circuit.num_cbits);
    out.custom_gates = circuit.custom_gates.clone();
    out.operations.reserve(circuit.operations.len());
    let mut pending = vec![Clifford::identity(); circuit.num_qubits];
    let flush = |q: usize, pending: &mut Vec<Clifford>, out: &mut Circuit| {
        if pending[q] == Clifford::identity() {
            return;
        }
        for gate in &frames().words[&pending[q]] {
            out.operations.push(gate.operation([q, q]));
        }
        pending[q] = Clifford::identity();
    };
    for op in operations {
        if let Operation::Gate {
            name,
            qubits,
            params,
        } = op
        {
            if params.is_empty() && qubits.len() == 1 {
                if *name == GateType::ID {
                    continue;
                }
                if let Some(gate) = local_gate(name, 0) {
                    pending[qubits[0]] = gate.table().compose(pending[qubits[0]]);
                    continue;
                }
            }
            for &q in qubits {
                flush(q, &mut pending, &mut out);
            }
        } else {
            for q in 0..circuit.num_qubits {
                flush(q, &mut pending, &mut out);
            }
        }
        out.operations.push(op.clone());
    }
    for q in 0..circuit.num_qubits {
        flush(q, &mut pending, &mut out);
    }
    out
}

/// (CX-equivalent cost, depth, two-qubit depth). Local gates preserve causal
/// dependencies in the latter even though they have zero entangling weight.
fn cost(circuit: &Circuit, elide: bool) -> (usize, usize, usize) {
    // The caller validates once and enables elision only for static circuits.
    // Project wire labels while scoring, without cloning/revalidating a whole
    // circuit. This is the same permutation update as the extraction pass.
    let mut permutation: Vec<_> = if elide {
        (0..circuit.num_qubits).collect()
    } else {
        Vec::new()
    };
    let mut count = 0;
    let mut layers = vec![[0, 0]; circuit.num_qubits];
    for op in &circuit.operations {
        if elide {
            if let Operation::Gate {
                name: GateType::SWAP,
                qubits,
                ..
            } = op
            {
                permutation.swap(qubits[0], qubits[1]);
                continue;
            }
        }
        let wires = op.qubits();
        let weight = match op {
            Operation::Gate {
                name: GateType::CX | GateType::CY | GateType::CZ,
                ..
            } => 1,
            Operation::Gate {
                name: GateType::SWAP,
                ..
            } => 3,
            _ => 0,
        };
        count += weight;
        let map = |q: usize| if elide { permutation[q] } else { q };
        let mut next = [0, 0];
        for &q in wires {
            let current = layers[map(q)];
            next[0] = next[0].max(current[0]);
            next[1] = next[1].max(current[1]);
        }
        next[0] += usize::from(!matches!(op, Operation::Barrier { .. }));
        next[1] += weight;
        for &q in wires {
            layers[map(q)] = next;
        }
    }
    (
        count,
        layers.iter().map(|v| v[0]).max().unwrap_or(0),
        layers.iter().map(|v| v[1]).max().unwrap_or(0),
    )
}

/// Fold the exact edited stream without constructing and then discarding an
/// intermediate deep clone. The slices preserve remove-end/splice-cut/
/// remove-start ordering, including either boundary cut.
fn folded_replacement(
    circuit: &Circuit,
    start: usize,
    end: usize,
    cut: usize,
    replacement: &[Operation],
) -> Circuit {
    fold_local(
        circuit,
        circuit.operations[..start]
            .iter()
            .chain(&circuit.operations[start + 1..cut])
            .chain(replacement)
            .chain(&circuit.operations[cut..end])
            .chain(&circuit.operations[end + 1..]),
    )
}

fn axis_commutes_with_interaction(
    axis: Pauli,
    wire: usize,
    name: &GateType,
    qubits: &[usize],
    pair: [usize; 2],
) -> bool {
    let Some(position) = qubits.iter().position(|&q| q == pair[wire]) else {
        return true;
    };
    let axis_on_wire = Pauli {
        x: axis.x >> wire,
        z: axis.z >> wire,
        phase: axis.phase,
    };
    let interaction_axis = match (name, position) {
        (GateType::CX | GateType::CY | GateType::CZ, 0) | (GateType::CZ, 1) => {
            Pauli::hermitian(0, 1)
        }
        (GateType::CX, 1) => Pauli::hermitian(1, 0),
        (GateType::CY, 1) => Pauli::hermitian(1, 1),
        _ => return false,
    };
    !axis_on_wire.anticommutes(interaction_axis)
}

/// All unsupported gates are fences, including operations disjoint from the
/// selected pair. Endpoint search may pass a noncommuting interaction, but the
/// independent forward/backward transports below cannot cross it.
fn supported_transport_operation(op: &Operation) -> bool {
    match op {
        Operation::Gate {
            name,
            qubits,
            params,
        } if params.is_empty() => {
            (qubits.len() == 1 && (*name == GateType::ID || local_gate(name, 0).is_some()))
                || (qubits.len() == 2 && matches!(name, GateType::CX | GateType::CY | GateType::CZ))
        }
        _ => false,
    }
}

fn transport_axes(
    (p, q): (Pauli, Pauli),
    op: &Operation,
    pair: [usize; 2],
    backwards: bool,
) -> Option<(Pauli, Pauli)> {
    let Operation::Gate {
        name,
        qubits,
        params,
    } = op
    else {
        return None;
    };
    if !params.is_empty() {
        return None;
    }
    if qubits.len() == 1 {
        if *name == GateType::ID {
            return Some((p, q));
        }
        let gate = local_gate(name, 0)?;
        if let Some(wire) = pair.iter().position(|w| *w == qubits[0]) {
            let transform = if backwards {
                gate.on_wire(wire).inverse_table()
            } else {
                gate.on_wire(wire).table()
            };
            return Some((transform.apply(p), transform.apply(q)));
        }
        return Some((p, q));
    }
    if qubits.len() == 2
        && matches!(name, GateType::CX | GateType::CY | GateType::CZ)
        && axis_commutes_with_interaction(p, 0, name, qubits, pair)
        && axis_commutes_with_interaction(q, 1, name, qubits, pair)
    {
        Some((p, q))
    } else {
        None
    }
}

impl Pass for CliffordSimplificationPass {
    fn name(&self) -> &str {
        "CliffordSimplificationPass"
    }

    fn run(&self, circuit: &Circuit, properties: &mut PropertySet) -> Circuit {
        self.run_internal(circuit, properties, cost, folded_replacement)
    }
}

impl CliffordSimplificationPass {
    // Static function arguments keep the production path monomorphized while
    // tests can compare the same search against the frozen allocating helpers.
    fn run_internal(
        &self,
        circuit: &Circuit,
        properties: &mut PropertySet,
        score: impl Fn(&Circuit, bool) -> (usize, usize, usize),
        replace: impl Fn(&Circuit, usize, usize, usize, &[Operation]) -> Circuit,
    ) -> Circuit {
        if super::validate_circuit_structure(circuit).is_err() {
            return circuit.clone();
        }
        let mut output = circuit.clone();
        let elide = self.allow_swap_replacements
            && circuit
                .operations
                .iter()
                .all(|op| matches!(op, Operation::Gate { .. } | Operation::Barrier { .. }));
        let mut quality = score(&output, elide);
        let mut stats = properties
            .get::<CliffordSimplificationStats>("clifford_simplification_stats")
            .cloned()
            .unwrap_or_default();
        let mut evaluations = 0;
        let mut left = Vec::with_capacity(LOOKAHEAD);
        'sweeps: for _ in 0..SWEEPS {
            let before = stats.rewrites;
            let mut start = 0;
            while start < output.operations.len() {
                let pair = match &output.operations[start] {
                    Operation::Gate {
                        name: GateType::CX,
                        qubits,
                        params,
                    } if params.is_empty() => [qubits[0], qubits[1]],
                    _ => {
                        start += 1;
                        continue;
                    }
                };
                // Only the first later CX on this pair is considered. An
                // unsupported operation anywhere in the bounded window fences
                // the search, even if neither endpoint uses its wires.
                let mut endpoint = None;
                for end in start + 1..output.operations.len().min(start + LOOKAHEAD + 1) {
                    let op = &output.operations[end];
                    if !supported_transport_operation(op) {
                        break;
                    }
                    if let Operation::Gate {
                        name: GateType::CX,
                        qubits,
                        ..
                    } = op
                    {
                        if qubits.iter().all(|wire| pair.contains(wire)) {
                            endpoint = Some((end, qubits[0] != pair[0]));
                            break;
                        }
                    }
                }
                let Some((end, reversed)) = endpoint else {
                    start += 1;
                    continue;
                };
                stats.candidates += 1;
                // left[cut-start-1] is L A L† immediately before `cut`.
                left.clear();
                left.push((Pauli::hermitian(0, 1), Pauli::hermitian(2, 0)));
                for cut in start + 1..end {
                    let Some(next) =
                        transport_axes(*left.last().unwrap(), &output.operations[cut], pair, false)
                    else {
                        break;
                    };
                    left.push(next);
                }
                // Right transport uses inverse local conjugation: R† B R.
                // Store one axis on each original wire even for reversed CX.
                let mut right = if reversed {
                    (Pauli::hermitian(1, 0), Pauli::hermitian(0, 2))
                } else {
                    (Pauli::hermitian(0, 1), Pauli::hermitian(2, 0))
                };
                'cuts: for cut in (start + 1..=end).rev() {
                    if cut < end {
                        let Some(next) = transport_axes(right, &output.operations[cut], pair, true)
                        else {
                            break;
                        };
                        right = next;
                    }
                    let Some(&(p, q)) = left.get(cut - start - 1) else {
                        continue;
                    };
                    stats.insertion_cuts += 1;
                    let product =
                        Clifford::controlled(right.0, right.1).compose(Clifford::controlled(p, q));
                    for swap in [false, true] {
                        if swap && !elide {
                            continue;
                        }
                        let target = if swap {
                            product.compose(Clifford::swap())
                        } else {
                            product
                        };
                        let Some(mut word) = synthesize(target) else {
                            continue;
                        };
                        if evaluations == MAX_CANDIDATE_EVALUATIONS {
                            stats.budget_exhausted = true;
                            break 'sweeps;
                        }
                        evaluations += 1;
                        stats.candidate_evaluations += 1;
                        if swap {
                            word.insert(0, SmallGate::Swap);
                        }
                        let replacement: Vec<_> =
                            word.into_iter().map(|g| g.operation(pair)).collect();
                        let candidate = replace(&output, start, end, cut, &replacement);
                        let next = score(&candidate, elide);
                        if next.0 < quality.0 && next.1 <= quality.1 && next.2 <= quality.2 {
                            stats.rewrites += 1;
                            stats.interior_rewrites += usize::from(cut > start + 1 && cut < end);
                            stats.cx_removed += quality.0 - next.0;
                            quality = next;
                            output = candidate;
                            break 'cuts;
                        }
                        stats.depth_rejections += 1;
                    }
                }
                start += 1;
            }
            if stats.rewrites == before {
                break;
            }
        }
        properties.insert("clifford_simplification_stats", stats);
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    use num_complex::Complex64 as C;

    // Frozen PR17 helpers intentionally retain their old allocations and
    // independent Circuit::depth/permutation-extraction calls for parity.
    fn legacy_fold_local(circuit: &Circuit) -> Circuit {
        let mut out = circuit.clone();
        out.operations.clear();
        let mut pending = vec![Clifford::identity(); circuit.num_qubits];
        let flush = |q: usize, pending: &mut Vec<Clifford>, out: &mut Circuit| {
            for gate in &frames().words[&pending[q]] {
                out.operations.push(gate.operation([q, q]));
            }
            pending[q] = Clifford::identity();
        };
        for op in &circuit.operations {
            if let Operation::Gate {
                name,
                qubits,
                params,
            } = op
            {
                if params.is_empty() && qubits.len() == 1 {
                    if *name == GateType::ID {
                        continue;
                    }
                    if let Some(gate) = local_gate(name, 0) {
                        pending[qubits[0]] = gate.table().compose(pending[qubits[0]]);
                        continue;
                    }
                }
                for &q in qubits {
                    flush(q, &mut pending, &mut out);
                }
            } else {
                for q in 0..circuit.num_qubits {
                    flush(q, &mut pending, &mut out);
                }
            }
            out.operations.push(op.clone());
        }
        for q in 0..circuit.num_qubits {
            flush(q, &mut pending, &mut out);
        }
        out
    }

    fn legacy_cost(circuit: &Circuit, elide: bool) -> (usize, usize, usize) {
        let projected;
        let c = if elide {
            projected = super::super::optimization::ExtractOutputPermutationPass
                .run(circuit, &mut PropertySet::new());
            &projected
        } else {
            circuit
        };
        let mut count = 0;
        let mut layers = vec![0; c.num_qubits];
        for op in &c.operations {
            let wires = op.qubits();
            let weight = match op {
                Operation::Gate {
                    name: GateType::CX | GateType::CY | GateType::CZ,
                    ..
                } => 1,
                Operation::Gate {
                    name: GateType::SWAP,
                    ..
                } => 3,
                _ => 0,
            };
            count += weight;
            let d = wires.iter().map(|&q| layers[q]).max().unwrap_or(0) + weight;
            for &q in wires {
                layers[q] = d;
            }
        }
        (count, c.depth(), layers.into_iter().max().unwrap_or(0))
    }

    fn legacy_replacement(
        circuit: &Circuit,
        start: usize,
        end: usize,
        cut: usize,
        replacement: &[Operation],
    ) -> Circuit {
        let mut candidate = circuit.clone();
        candidate.operations.remove(end);
        candidate
            .operations
            .splice(cut..cut, replacement.iter().cloned());
        candidate.operations.remove(start);
        legacy_fold_local(&candidate)
    }

    // Literal matrices: these do not call gate_def, the simulator, numerical
    // synthesis or the production signed-Pauli action/compose functions.
    fn dense_pauli(p: Pauli) -> DMatrix<C> {
        let mut matrix = DMatrix::zeros(4, 4);
        let phase = [
            C::new(1.0, 0.0),
            C::new(0.0, 1.0),
            C::new(-1.0, 0.0),
            C::new(0.0, -1.0),
        ][p.phase as usize];
        for column in 0..4 {
            let sign = if ((column as u8 & p.z).count_ones() % 2) == 0 {
                1.0
            } else {
                -1.0
            };
            matrix[(column ^ p.x as usize, column)] = phase * sign;
        }
        matrix
    }

    fn dense_gate(gate: SmallGate) -> DMatrix<C> {
        let one = C::new(1.0, 0.0);
        let zero = C::new(0.0, 0.0);
        let i = C::new(0.0, 1.0);
        let mut matrix = DMatrix::zeros(4, 4);
        match gate {
            SmallGate::CX(control) => {
                for column in 0..4 {
                    let row = if column & (1 << control) != 0 {
                        column ^ (1 << (1 - control))
                    } else {
                        column
                    };
                    matrix[(row, column)] = one;
                }
            }
            SmallGate::Swap => {
                for column in 0..4 {
                    matrix[(((column & 1) << 1) | (column >> 1), column)] = one;
                }
            }
            _ => {
                let (q, u) = match gate {
                    SmallGate::H(q) => (
                        q,
                        [
                            [one / 2.0_f64.sqrt(), one / 2.0_f64.sqrt()],
                            [one / 2.0_f64.sqrt(), -one / 2.0_f64.sqrt()],
                        ],
                    ),
                    SmallGate::S(q) => (q, [[one, zero], [zero, i]]),
                    SmallGate::Sdg(q) => (q, [[one, zero], [zero, -i]]),
                    SmallGate::X(q) => (q, [[zero, one], [one, zero]]),
                    SmallGate::Y(q) => (q, [[zero, -i], [i, zero]]),
                    SmallGate::Z(q) => (q, [[one, zero], [zero, -one]]),
                    SmallGate::SX(q) => (
                        q,
                        [
                            [(one + i) / 2.0, (one - i) / 2.0],
                            [(one - i) / 2.0, (one + i) / 2.0],
                        ],
                    ),
                    _ => unreachable!(),
                };
                for column in 0..4 {
                    for bit in 0..2 {
                        let row = (column & !(1 << q)) | (bit << q);
                        matrix[(row, column)] = u[bit][(column >> q) & 1];
                    }
                }
            }
        }
        matrix
    }

    fn assert_dense_equal(left: &DMatrix<C>, right: &DMatrix<C>) {
        let overlap: C = left
            .iter()
            .zip(right.iter())
            .map(|(a, b)| a.conj() * b)
            .sum();
        assert!(overlap.norm() > 0.0);
        let phase = overlap / overlap.norm();
        let residual = left
            .iter()
            .zip(right.iter())
            .map(|(a, b)| (*b - phase * a).norm_sqr())
            .sum::<f64>()
            .sqrt();
        assert!(
            residual < 3e-14,
            "full matrix shared-phase residual {residual}"
        );
    }

    fn word_table(word: &[SmallGate]) -> Clifford {
        word.iter()
            .fold(Clifford::identity(), |a, gate| gate.table().compose(a))
    }

    fn parity_fixture(seed: u64, dynamic: bool) -> Circuit {
        let width = 2 + seed as usize % 5;
        let mut circuit = Circuit::new(width, width);
        circuit.register_custom_gate("opaque".into(), vec![], vec!["a".into()], vec![]);
        let mut state = seed + 1;
        for index in 0..64 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let a = ((state >> 32) as usize) % width;
            let b = (a + 1 + ((state >> 40) as usize) % (width - 1)) % width;
            let operation = match (state >> 48) % 17 {
                0 => SmallGate::Swap.operation([a, b]),
                1 => Operation::Barrier { qubits: vec![] },
                2 => Operation::Barrier { qubits: vec![a, b] },
                3 => Operation::Gate {
                    name: GateType::U,
                    qubits: vec![a],
                    params: vec![0.13, -0.27, 0.41],
                },
                4 => Operation::Gate {
                    name: GateType::Custom("opaque".into()),
                    qubits: vec![a],
                    params: vec![],
                },
                5 => SmallGate::H(0).operation([a, b]),
                6 => SmallGate::S(0).operation([a, b]),
                7 => SmallGate::Sdg(0).operation([a, b]),
                8 => SmallGate::SX(0).operation([a, b]),
                9 => SmallGate::Y(0).operation([a, b]),
                10 => Operation::Gate {
                    name: GateType::ID,
                    qubits: vec![a],
                    params: vec![],
                },
                11 => Operation::Gate {
                    name: GateType::CY,
                    qubits: vec![a, b],
                    params: vec![],
                },
                12 => Operation::Gate {
                    name: GateType::CZ,
                    qubits: vec![a, b],
                    params: vec![],
                },
                _ => SmallGate::CX(0).operation([a, b]),
            };
            circuit.operations.push(operation);
            if index % 11 == 0 {
                circuit.operations.extend([
                    SmallGate::CX(0).operation([a, b]),
                    SmallGate::S(0).operation([a, b]),
                    SmallGate::CX(0).operation([a, b]),
                ]);
            }
            if dynamic && index % 13 == 0 {
                circuit.operations.extend([
                    Operation::Measure { qubit: a, cbit: a },
                    Operation::Reset { qubit: b },
                ]);
            }
        }
        circuit
    }

    #[test]
    fn streaming_cost_matches_extraction_and_literal_depth_with_all_fences() {
        for seed in 0..96 {
            for dynamic in [false, true] {
                let circuit = parity_fixture(seed, dynamic);
                assert_eq!(cost(&circuit, false), legacy_cost(&circuit, false));
                if !dynamic {
                    assert_eq!(cost(&circuit, true), legacy_cost(&circuit, true));
                }
            }
        }
        // A conditional barrier has one literal layer under Circuit::depth,
        // unlike a top-level barrier. Preserve that existing convention.
        let mut circuit = Circuit::new(3, 1);
        circuit.operations = vec![
            SmallGate::CX(0).operation([0, 2]),
            Operation::Conditional {
                condition: crate::ir::ClassicalCondition {
                    creg: "c".into(),
                    value: 1,
                },
                op: Box::new(Operation::Barrier { qubits: vec![0, 1] }),
            },
            Operation::Barrier { qubits: vec![] },
            SmallGate::H(0).operation([1, 2]),
        ];
        assert_eq!(cost(&circuit, false), legacy_cost(&circuit, false));
    }

    #[test]
    fn streamed_edit_matches_legacy_order_and_metadata_at_every_cut() {
        for seed in 0..8 {
            let circuit = parity_fixture(seed, true);
            for start in [0, 3, 11] {
                let end = start + 8;
                for cut in start + 1..=end {
                    for replacement in [
                        vec![],
                        vec![SmallGate::H(0).operation([0, 1])],
                        vec![
                            SmallGate::Swap.operation([1, 0]),
                            SmallGate::CX(0).operation([1, 0]),
                        ],
                    ] {
                        assert_eq!(
                            folded_replacement(&circuit, start, end, cut, &replacement),
                            legacy_replacement(&circuit, start, end, cut, &replacement),
                            "seed {seed}, cut {start}/{cut}/{end}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn full_search_matches_legacy_helpers_in_outputs_and_cumulative_counters() {
        for seed in 0..96 {
            for dynamic in [false, true] {
                let circuit = parity_fixture(seed, dynamic);
                for allow_swap_replacements in [false, true] {
                    assert!(super::super::validate_circuit_structure(&circuit).is_ok());
                    let pass = CliffordSimplificationPass {
                        allow_swap_replacements,
                    };
                    let mut current = PropertySet::new();
                    let mut legacy = PropertySet::new();
                    let prefix = CliffordSimplificationStats {
                        candidates: 3,
                        rewrites: 2,
                        cx_removed: 4,
                        candidate_evaluations: 7,
                        ..CliffordSimplificationStats::default()
                    };
                    current.insert("clifford_simplification_stats", prefix.clone());
                    legacy.insert("clifford_simplification_stats", prefix);
                    let actual = pass.run(&circuit, &mut current);
                    let expected =
                        pass.run_internal(&circuit, &mut legacy, legacy_cost, legacy_replacement);
                    assert_eq!(
                        actual, expected,
                        "seed {seed}, dynamic {dynamic}, elide {allow_swap_replacements}"
                    );
                    assert_eq!(
                        current.get::<CliffordSimplificationStats>("clifford_simplification_stats"),
                        legacy.get::<CliffordSimplificationStats>("clifford_simplification_stats"),
                        "seed {seed}, dynamic {dynamic}, elide {allow_swap_replacements}"
                    );
                }
            }
        }
    }

    #[test]
    fn local_group_is_complete_and_inverse_and_synthesis_are_exact() {
        assert_eq!(frames().words.len(), 24);
        for (inverse, word) in &frames().pairs {
            let table = word_table(word);
            assert_eq!(table.compose(*inverse), Clifford::identity());
            assert_eq!(inverse.compose(table), Clifford::identity());
            assert_eq!(word_table(&local_word(table).unwrap()), table);
        }
    }

    #[test]
    fn large_rewriteable_input_has_a_deterministic_linear_work_budget() {
        let mut input = Circuit::new(2, 0);
        for _ in 0..256 {
            input.operations.push(SmallGate::CX(0).operation([0, 1]));
            input.operations.push(SmallGate::CX(0).operation([0, 1]));
        }
        let mut properties = PropertySet::new();
        let output = CliffordSimplificationPass::default().run(&input, &mut properties);
        let stats = properties
            .get::<CliffordSimplificationStats>("clifford_simplification_stats")
            .unwrap();
        assert_eq!(stats.candidate_evaluations, MAX_CANDIDATE_EVALUATIONS);
        assert!(stats.budget_exhausted);
        assert_eq!(stats.rewrites, MAX_CANDIDATE_EVALUATIONS);
        assert!(output.depth() < input.depth());
        // The residual circuit still contains an even number of identical
        // CX involutions. Check its entire operator with literal matrices.
        assert_eq!(output.operations.len() % 2, 0);
        let matrix = output
            .operations
            .iter()
            .fold(DMatrix::identity(4, 4), |m, op| {
                assert_eq!(*op, SmallGate::CX(0).operation([0, 1]));
                dense_gate(SmallGate::CX(0)) * m
            });
        assert_dense_equal(&DMatrix::identity(4, 4), &matrix);
        let repeated = CliffordSimplificationPass::default().run(&input, &mut PropertySet::new());
        assert_eq!(output.operations, repeated.operations);
    }

    #[test]
    fn every_signed_axis_endpoint_product_is_synthesized_exactly_if_reduced() {
        let axes = [(1, 0), (0, 1), (1, 1)];
        let mut recognized = 0;
        for (x0, z0) in axes {
            for (x1, z1) in axes {
                for signs in 0..4 {
                    let mut p = Pauli::hermitian(x0, z0);
                    let mut q = Pauli::hermitian(x1 << 1, z1 << 1);
                    if signs & 1 != 0 {
                        p = p.neg();
                    }
                    if signs & 2 != 0 {
                        q = q.neg();
                    }
                    for direction in 0..2 {
                        let target = SmallGate::CX(direction)
                            .table()
                            .compose(Clifford::controlled(p, q));
                        let mut found = false;
                        for swap in [false, true] {
                            let reduced = if swap {
                                target.compose(Clifford::swap())
                            } else {
                                target
                            };
                            if let Some(mut word) = synthesize(reduced) {
                                assert!(
                                    word.iter()
                                        .filter(|g| matches!(g, SmallGate::CX(_)))
                                        .count()
                                        <= 1
                                );
                                if swap {
                                    word.insert(0, SmallGate::Swap);
                                }
                                assert_eq!(word_table(&word), target);
                                let pm = dense_pauli(p);
                                let qm = dense_pauli(q);
                                let controlled = (DMatrix::identity(4, 4) + &pm + &qm - &pm * &qm)
                                    * C::new(0.5, 0.0);
                                let expected = dense_gate(SmallGate::CX(direction)) * controlled;
                                let actual = word
                                    .iter()
                                    .fold(DMatrix::identity(4, 4), |m, g| dense_gate(*g) * m);
                                assert_dense_equal(&expected, &actual);
                                found = true;
                            }
                        }
                        assert!(found);
                        recognized += 1;
                    }
                }
            }
        }
        assert_eq!(recognized, 72);
    }

    #[test]
    fn backward_local_conjugations_use_exact_adjoint_actions() {
        for q in 0..2 {
            for gate in [
                SmallGate::H(q),
                SmallGate::S(q),
                SmallGate::Sdg(q),
                SmallGate::X(q),
                SmallGate::Y(q),
                SmallGate::Z(q),
                SmallGate::SX(q),
            ] {
                let inverse = gate.inverse_table();
                assert_eq!(inverse.compose(gate.table()), Clifford::identity());
                assert_eq!(gate.table().compose(inverse), Clifford::identity());
                let u = dense_gate(gate);
                for p in Clifford::identity().0 {
                    let expected = u.adjoint() * dense_pauli(p) * &u;
                    let actual = dense_pauli(inverse.apply(p));
                    // Pauli signs are observable here: do not align a phase.
                    assert!((&expected - &actual).norm() < 1e-14);
                }
            }
        }
    }

    #[test]
    fn every_pair_of_transported_signed_interactions_has_exact_replacement() {
        let mut interactions = Vec::new();
        for (x0, z0) in [(1, 0), (0, 1), (1, 1)] {
            for (x1, z1) in [(1, 0), (0, 1), (1, 1)] {
                for signs in 0..4 {
                    let mut p = Pauli::hermitian(x0, z0);
                    let mut q = Pauli::hermitian(x1 << 1, z1 << 1);
                    if signs & 1 != 0 {
                        p = p.neg();
                    }
                    if signs & 2 != 0 {
                        q = q.neg();
                    }
                    let pm = dense_pauli(p);
                    let qm = dense_pauli(q);
                    let matrix =
                        (DMatrix::identity(4, 4) + &pm + &qm - &pm * &qm) * C::new(0.5, 0.0);
                    interactions.push((Clifford::controlled(p, q), matrix));
                }
            }
        }
        assert_eq!(interactions.len(), 36);
        let mut verified = 0;
        for (left, left_matrix) in &interactions {
            for (right, right_matrix) in &interactions {
                let product = right.compose(*left);
                let expected = right_matrix * left_matrix;
                let mut found = false;
                for swap in [false, true] {
                    let target = if swap {
                        product.compose(Clifford::swap())
                    } else {
                        product
                    };
                    if let Some(mut word) = synthesize(target) {
                        if swap {
                            word.insert(0, SmallGate::Swap);
                        }
                        assert_eq!(word_table(&word), product);
                        let actual = word.iter().fold(DMatrix::identity(4, 4), |matrix, gate| {
                            dense_gate(*gate) * matrix
                        });
                        assert_dense_equal(&expected, &actual);
                        found = true;
                    }
                }
                assert!(
                    found,
                    "a two-interaction Clifford must reduce with SWAP freedom"
                );
                verified += 1;
            }
        }
        assert_eq!(verified, 1296);
    }
}
