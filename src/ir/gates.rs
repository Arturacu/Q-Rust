use std::f64::consts::PI;
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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
            "u" | "u3" | "U" => GateType::U,
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

impl GateType {
    /// Returns the U(theta, phi, lambda) parameters for single-qubit gates, or None if not applicable.
    pub fn single_qubit_u_params(&self, params: &[f64]) -> Option<(f64, f64, f64)> {
        match self {
            GateType::ID => Some((0.0, 0.0, 0.0)),
            GateType::X => Some((PI, 0.0, PI)),
            GateType::Y => Some((PI, PI / 2.0, PI / 2.0)),
            GateType::Z => Some((0.0, 0.0, PI)),
            GateType::H => Some((PI / 2.0, 0.0, PI)),
            GateType::S => Some((0.0, 0.0, PI / 2.0)),
            GateType::Sdg => Some((0.0, 0.0, -PI / 2.0)),
            GateType::T => Some((0.0, 0.0, PI / 4.0)),
            GateType::Tdg => Some((0.0, 0.0, -PI / 4.0)),
            GateType::RX => Some((params[0], -PI / 2.0, PI / 2.0)),
            GateType::RY => Some((params[0], 0.0, 0.0)),
            GateType::RZ => Some((0.0, 0.0, params[0])),
            _ => None,
        }
    }
    /// Returns the OpenQASM 2.0 standard string representation for the gate type.
    pub fn to_qasm_name(&self) -> String {
        match self {
            GateType::H => "h".to_string(),
            GateType::X => "x".to_string(),
            GateType::Y => "y".to_string(),
            GateType::Z => "z".to_string(),
            GateType::CX => "cx".to_string(),
            GateType::RX => "rx".to_string(),
            GateType::RY => "ry".to_string(),
            GateType::RZ => "rz".to_string(),
            GateType::U => "u".to_string(),
            GateType::ID => "id".to_string(),
            GateType::S => "s".to_string(),
            GateType::Sdg => "sdg".to_string(),
            GateType::T => "t".to_string(),
            GateType::Tdg => "tdg".to_string(),
            GateType::SWAP => "swap".to_string(),
            GateType::CCX => "ccx".to_string(),
            GateType::CZ => "cz".to_string(),
            GateType::CY => "cy".to_string(),
            GateType::CRX => "crx".to_string(),
            GateType::CRY => "cry".to_string(),
            GateType::CRZ => "crz".to_string(),
            GateType::CH => "ch".to_string(),
            GateType::CSX => "csx".to_string(),
            GateType::RXX => "rxx".to_string(),
            GateType::RYY => "ryy".to_string(),
            GateType::RZZ => "rzz".to_string(),
            GateType::Custom(name) => name.clone(),
        }
    }
}
