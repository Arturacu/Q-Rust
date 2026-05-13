//! Parser-internal AST for OpenQASM 2.0 statements and expressions.

use crate::error::{QRustError, Result};
use std::collections::HashMap;
use std::f64::consts::PI;

#[derive(Debug, PartialEq, Clone)]
#[non_exhaustive]
pub enum Expr {
    Float(f64),
    Var(String),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
}

impl Expr {
    pub fn evaluate_with_scope(&self, scope: &HashMap<String, f64>) -> Result<f64> {
        match self {
            Expr::Float(v) => Ok(*v),
            Expr::Var(name) => {
                if name == "pi" {
                    Ok(PI)
                } else {
                    scope
                        .get(name)
                        .copied()
                        .ok_or_else(|| QRustError::Undefined(name.clone()))
                }
            }
            Expr::Add(l, r) => Ok(l.evaluate_with_scope(scope)? + r.evaluate_with_scope(scope)?),
            Expr::Sub(l, r) => Ok(l.evaluate_with_scope(scope)? - r.evaluate_with_scope(scope)?),
            Expr::Mul(l, r) => Ok(l.evaluate_with_scope(scope)? * r.evaluate_with_scope(scope)?),
            Expr::Div(l, r) => {
                let d = r.evaluate_with_scope(scope)?;
                if d == 0.0 {
                    Err(QRustError::Arithmetic("division by zero".into()))
                } else {
                    Ok(l.evaluate_with_scope(scope)? / d)
                }
            }
        }
    }

    #[inline]
    pub fn evaluate(&self) -> Result<f64> {
        self.evaluate_with_scope(&HashMap::new())
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ParsedStatement {
    QReg(String, usize),
    CReg(String, usize),
    Gate(String, Vec<(String, Option<usize>)>, Vec<Expr>),
    Measure((String, Option<usize>), (String, Option<usize>)),
    Include(String),
    Barrier(Vec<(String, Option<usize>)>),
    GateDef(String, Vec<String>, Vec<String>, Vec<ParsedStatement>),
    If(String, usize, Box<ParsedStatement>),
    Ignore,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expr_float_literal() {
        assert_eq!(Expr::Float(3.14).evaluate().unwrap(), 3.14);
    }

    #[test]
    fn test_expr_pi_constant() {
        assert_eq!(Expr::Var("pi".into()).evaluate().unwrap(), PI);
    }

    #[test]
    fn test_expr_unknown_variable() {
        assert!(matches!(
            Expr::Var("theta".into()).evaluate(),
            Err(QRustError::Undefined(_))
        ));
    }

    #[test]
    fn test_expr_addition() {
        let e = Expr::Add(Box::new(Expr::Float(2.0)), Box::new(Expr::Float(3.0)));
        assert_eq!(e.evaluate().unwrap(), 5.0);
    }

    #[test]
    fn test_expr_division_by_zero() {
        let e = Expr::Div(Box::new(Expr::Float(10.0)), Box::new(Expr::Float(0.0)));
        assert!(matches!(e.evaluate(), Err(QRustError::Arithmetic(_))));
    }

    #[test]
    fn test_expr_complex_expression() {
        let e = Expr::Add(
            Box::new(Expr::Div(
                Box::new(Expr::Var("pi".into())),
                Box::new(Expr::Float(2.0)),
            )),
            Box::new(Expr::Float(1.0)),
        );
        let got = e.evaluate().unwrap();
        assert!((got - (PI / 2.0 + 1.0)).abs() < 1e-12);
    }
}