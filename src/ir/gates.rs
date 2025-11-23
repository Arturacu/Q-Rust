/// Quantum Gate Types
///
/// This enum represents the set of supported quantum gates.
/// It includes standard single-qubit gates (H, X, Y, Z),
/// two-qubit gates (CX), and parameterized rotation gates (RX, RY, RZ).
///
/// # Examples
///
/// ```
/// use q_rust::ir::GateType;
/// let h_gate = GateType::H;
/// let rx_gate = GateType::RX(1.57);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub enum GateType {
    /// Hadamard gate
    H,
    /// Pauli-X gate (NOT)
    X,
    /// Pauli-Y gate
    Y,
    /// Pauli-Z gate
    Z,
    /// Controlled-NOT gate
    CX,
    /// Rotation around X-axis with angle theta
    RX(f64),
    /// Rotation around Y-axis with angle theta
    RY(f64),
    /// Rotation around Z-axis with angle theta
    RZ(f64),
    /// General unitary gate U(theta, phi, lambda)
    U(f64, f64, f64),
    /// Identity gate (wait)
    ID,
    /// S gate (sqrt(Z))
    S,
    /// S-dagger gate (inverse of S)
    Sdg,
    /// T gate (sqrt(S))
    T,
    /// T-dagger gate (inverse of T)
    Tdg,
    /// Swap gate
    SWAP,
    /// Toffoli gate (CCX)
    CCX,
    // Add more as needed
    /// Custom user-defined gate
    Custom(String),
}
