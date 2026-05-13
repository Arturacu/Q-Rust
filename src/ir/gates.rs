//! Gate type enumeration.

use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde-ir", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum GateType {
    H,
    X,
    Y,
    Z,
    S,
    Sdg,
    T,
    Tdg,
    ID,
    RX,
    RY,
    RZ,
    U,
    CX,
    CY,
    CZ,
    CH,
    CSX,
    CRX,
    CRY,
    CRZ,
    RXX,
    RYY,
    RZZ,
    SWAP,
    CCX,
    /// IBM-Heron native two-qubit echoed cross-resonance gate.
    /// `ECR = (1/√2)(IX − XY)`. Reference: Sheldon et al. 2016, PRA 93, 060302.
    ECR,
    /// Imaginary SWAP — native on IonQ Forte and superconducting fluxonium devices.
    /// In the {00, 01, 10, 11} basis: `iSWAP = diag(1, [[0, i], [i, 0]], 1)`.
    /// Reference: Schuch & Siewert 2003, PRA 67, 032301.
    ISwap,
    /// A barrier pseudo-gate placeholder. Not used as a real gate
    /// (barriers are represented by `Operation::Barrier`), but exists so
    /// `count_ops` and similar tools can report barriers uniformly.
    Barrier,
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
            "u1" => GateType::RZ,
            "u2" | "u" | "u3" | "U" => GateType::U,
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
            "ecr" => GateType::ECR,
            "iswap" => GateType::ISwap,
            "barrier" => GateType::Barrier,
            other => GateType::Custom(other.to_string()),
        })
    }
}

impl GateType {
    pub fn to_qasm_name(&self) -> &str {
        match self {
            GateType::H => "h",
            GateType::X => "x",
            GateType::Y => "y",
            GateType::Z => "z",
            GateType::CX => "cx",
            GateType::RX => "rx",
            GateType::RY => "ry",
            GateType::RZ => "rz",
            GateType::U => "u",
            GateType::ID => "id",
            GateType::S => "s",
            GateType::Sdg => "sdg",
            GateType::T => "t",
            GateType::Tdg => "tdg",
            GateType::SWAP => "swap",
            GateType::CCX => "ccx",
            GateType::CZ => "cz",
            GateType::CY => "cy",
            GateType::CRX => "crx",
            GateType::CRY => "cry",
            GateType::CRZ => "crz",
            GateType::CH => "ch",
            GateType::CSX => "csx",
            GateType::RXX => "rxx",
            GateType::RYY => "ryy",
            GateType::RZZ => "rzz",
            GateType::ECR => "ecr",
            GateType::ISwap => "iswap",
            GateType::Barrier => "barrier",
            GateType::Custom(n) => n.as_str(),
        }
    }

    pub fn static_qasm_name(&self) -> Option<&'static str> {
        Some(match self {
            GateType::H => "h",
            GateType::X => "x",
            GateType::Y => "y",
            GateType::Z => "z",
            GateType::CX => "cx",
            GateType::RX => "rx",
            GateType::RY => "ry",
            GateType::RZ => "rz",
            GateType::U => "u",
            GateType::ID => "id",
            GateType::S => "s",
            GateType::Sdg => "sdg",
            GateType::T => "t",
            GateType::Tdg => "tdg",
            GateType::SWAP => "swap",
            GateType::CCX => "ccx",
            GateType::CZ => "cz",
            GateType::CY => "cy",
            GateType::CRX => "crx",
            GateType::CRY => "cry",
            GateType::CRZ => "crz",
            GateType::CH => "ch",
            GateType::CSX => "csx",
            GateType::RXX => "rxx",
            GateType::RYY => "ryy",
            GateType::RZZ => "rzz",
            GateType::ECR => "ecr",
            GateType::ISwap => "iswap",
            GateType::Barrier => "barrier",
            GateType::Custom(_) => return None,
        })
    }
}

impl fmt::Display for GateType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.to_qasm_name())
    }
}
