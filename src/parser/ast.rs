/// Internal AST for parsed statements
#[derive(Debug, PartialEq, Clone)]
#[allow(dead_code)]
pub enum ParsedStatement {
    QReg(String, usize),
    CReg(String, usize),
    Gate(String, Vec<(String, usize)>, Vec<f64>),
    Measure((String, usize), (String, usize)),
    Barrier(Vec<(String, usize)>),
    Reset((String, usize)),
    Ignore,
}
