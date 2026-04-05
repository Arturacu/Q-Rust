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
/// let rx_gate = GateType::RX; // Parameters are now stored in Operation::params
/// ```
use std::str::FromStr;

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
    RX,
    /// Rotation around Y-axis with angle theta
    RY,
    /// Rotation around Z-axis with angle theta
    RZ,
    /// General unitary gate U(theta, phi, lambda)
    U,
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
    /// Controlled-Z gate
    CZ,
    /// Controlled-Y gate
    CY,
    /// Controlled-RX gate
    CRX,
    /// Controlled-RY gate
    CRY,
    /// Controlled-RZ gate
    CRZ,
    /// Controlled-Hadamard gate
    CH,
    /// Controlled-√X gate
    CSX,
    /// Ising XX interaction gate
    RXX,
    /// Ising YY interaction gate
    RYY,
    /// Ising ZZ interaction gate
    RZZ,
    /// Custom user-defined gate
    Custom(String),
}

impl FromStr for GateType {
    type Err = ();

    fn from_str(name: &str) -> Result<Self, Self::Err> {
        Ok(match name {
            "h" => GateType::H,
            "x" => GateType::X,
            "y" => GateType::Y,
            "z" => GateType::Z,
            "cx" => GateType::CX,
            "rx" => GateType::RX,
            "ry" => GateType::RY,
            "rz" => GateType::RZ,
            "u1" => GateType::RZ, // u1(lambda) = RZ(lambda)
            "u2" => GateType::U,  // u2(phi, lambda) = U(pi/2, phi, lambda)
            "u3" | "U" => GateType::U,
            "id" => GateType::ID,
            "s" => GateType::S,
            "sdg" => GateType::Sdg,
            "t" => GateType::T,
            "tdg" => GateType::Tdg,
            "swap" => GateType::SWAP,
            "ccx" => GateType::CCX,
            "cz" => GateType::CZ,
            "cy" => GateType::CY,
            "crx" => GateType::CRX,
            "cry" => GateType::CRY,
            "crz" => GateType::CRZ,
            "ch" => GateType::CH,
            "csx" => GateType::CSX,
            "rxx" => GateType::RXX,
            "ryy" => GateType::RYY,
            "rzz" => GateType::RZZ,
            _ => GateType::Custom(name.to_string()),
        })
    }
}
