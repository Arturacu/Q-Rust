/// Internal AST for parsed statements
#[derive(Debug, PartialEq, Clone)]
#[allow(dead_code)]
pub enum Expr {
    Float(f64),
    Var(String),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
}

impl Expr {
    pub fn evaluate(&self) -> Result<f64, String> {
        match self {
            Expr::Float(val) => Ok(*val),
            Expr::Var(name) => {
                if name == "pi" {
                    Ok(std::f64::consts::PI)
                } else {
                    Err(format!("Unknown variable: {}", name))
                }
            }
            Expr::Add(lhs, rhs) => Ok(lhs.evaluate()? + rhs.evaluate()?),
            Expr::Sub(lhs, rhs) => Ok(lhs.evaluate()? - rhs.evaluate()?),
            Expr::Mul(lhs, rhs) => Ok(lhs.evaluate()? * rhs.evaluate()?),
            Expr::Div(lhs, rhs) => {
                let denom = rhs.evaluate()?;
                if denom == 0.0 {
                    Err("Division by zero".to_string())
                } else {
                    Ok(lhs.evaluate()? / denom)
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ParsedStatement {
    QReg(String, usize),
    CReg(String, usize),
    Gate(String, Vec<(String, Option<usize>)>, Vec<Expr>), // Name, Qubits, Params
    Measure((String, Option<usize>), (String, Option<usize>)), // Qubit -> Cbit
    Include(String),                                       // Filename
    Barrier(Vec<(String, Option<usize>)>),                 // Qubits
    GateDef(String, Vec<String>, Vec<String>, Vec<ParsedStatement>), // Name, Params, Qubits, Body
    If(String, usize, Box<ParsedStatement>),               // CReg, Val, Op
    Ignore,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expr_float_literal() {
        let expr = Expr::Float(3.14);
        assert_eq!(expr.evaluate(), Ok(3.14));
    }

    #[test]
    fn test_expr_pi_constant() {
        let expr = Expr::Var("pi".to_string());
        assert_eq!(expr.evaluate(), Ok(std::f64::consts::PI));
    }

    #[test]
    fn test_expr_unknown_variable() {
        let expr = Expr::Var("theta".to_string());
        assert!(expr.evaluate().is_err());
        assert!(expr.evaluate().unwrap_err().contains("Unknown variable"));
    }

    #[test]
    fn test_expr_addition() {
        let expr = Expr::Add(Box::new(Expr::Float(2.0)), Box::new(Expr::Float(3.0)));
        assert_eq!(expr.evaluate(), Ok(5.0));
    }

    #[test]
    fn test_expr_subtraction() {
        let expr = Expr::Sub(Box::new(Expr::Float(5.0)), Box::new(Expr::Float(3.0)));
        assert_eq!(expr.evaluate(), Ok(2.0));
    }

    #[test]
    fn test_expr_multiplication() {
        let expr = Expr::Mul(Box::new(Expr::Float(4.0)), Box::new(Expr::Float(2.5)));
        assert_eq!(expr.evaluate(), Ok(10.0));
    }

    #[test]
    fn test_expr_division() {
        let expr = Expr::Div(Box::new(Expr::Float(10.0)), Box::new(Expr::Float(2.0)));
        assert_eq!(expr.evaluate(), Ok(5.0));
    }

    #[test]
    fn test_expr_division_by_zero() {
        let expr = Expr::Div(Box::new(Expr::Float(10.0)), Box::new(Expr::Float(0.0)));
        assert!(expr.evaluate().is_err());
        assert!(expr.evaluate().unwrap_err().contains("Division by zero"));
    }

    #[test]
    fn test_expr_complex_expression() {
        // (pi / 2) + 1
        let expr = Expr::Add(
            Box::new(Expr::Div(
                Box::new(Expr::Var("pi".to_string())),
                Box::new(Expr::Float(2.0)),
            )),
            Box::new(Expr::Float(1.0)),
        );
        let result = expr.evaluate().unwrap();
        let expected = std::f64::consts::PI / 2.0 + 1.0;
        assert!((result - expected).abs() < 1e-10);
    }
}
