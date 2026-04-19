//! Nom grammar rules for the OpenQASM 2.0 parser.

use crate::ir::ast::{Expr, ParsedStatement};
use nom::{
    branch::alt,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::{alpha1, alphanumeric1, char, digit1, space0, space1},
    combinator::{map, map_res, opt, recognize, value},
    multi::{many0, separated_list0},
    number::complete::recognize_float,
    sequence::{delimited, pair, tuple},
    IResult,
};

fn identifier(input: &str) -> IResult<&str, String> {
    map(
        recognize(pair(
            alt((alpha1, tag("_"))),
            many0(alt((alphanumeric1, tag("_")))),
        )),
        str::to_string,
    )(input)
}

fn usize_parser(input: &str) -> IResult<&str, usize> {
    map_res(digit1, str::parse::<usize>)(input)
}

fn parse_f64(input: &str) -> IResult<&str, f64> {
    map_res(recognize_float, str::parse::<f64>)(input)
}

pub fn comment(input: &str) -> IResult<&str, ()> {
    value((), pair(tag("//"), take_while(|c| c != '\n')))(input)
}

pub fn openqasm_version(input: &str) -> IResult<&str, String> {
    map(
        tuple((
            tag("OPENQASM"),
            space1,
            take_while1(|c: char| c != ';'),
            space0,
            tag(";"),
        )),
        |(_, _, version, _, _): (&str, &str, &str, &str, &str)| version.to_string(),
    )(input)
}

pub fn include(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("include"),
            space1::<&str, nom::error::Error<&str>>,
            delimited(char('"'), take_while1(|c| c != '"'), char('"')),
            space0::<&str, nom::error::Error<&str>>,
            tag(";"),
        )),
        |(_, _, filename, _, _)| ParsedStatement::Include(filename.to_string()),
    )(input)
}

pub fn qreg(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("qreg"),
            space1,
            identifier,
            space0,
            delimited(char('['), usize_parser, char(']')),
            space0,
            tag(";"),
        )),
        |(_, _, name, _, size, _, _)| ParsedStatement::QReg(name, size),
    )(input)
}

pub fn creg(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("creg"),
            space1,
            identifier,
            space0,
            delimited(char('['), usize_parser, char(']')),
            space0,
            tag(";"),
        )),
        |(_, _, name, _, size, _, _)| ParsedStatement::CReg(name, size),
    )(input)
}

fn qubit_ref(input: &str) -> IResult<&str, (String, Option<usize>)> {
    pair(
        identifier,
        opt(delimited(
            tuple((space0, char('['), space0)),
            usize_parser,
            tuple((space0, char(']'), space0)),
        )),
    )(input)
}

fn term(input: &str) -> IResult<&str, Expr> {
    let (input, init) = factor(input)?;
    let (input, res) = many0(pair(
        delimited(space0, alt((char('*'), char('/'))), space0),
        factor,
    ))(input)?;
    Ok((
        input,
        res.into_iter().fold(init, |acc, (op, val)| match op {
            '*' => Expr::Mul(Box::new(acc), Box::new(val)),
            '/' => Expr::Div(Box::new(acc), Box::new(val)),
            _ => unreachable!(),
        }),
    ))
}

fn factor(input: &str) -> IResult<&str, Expr> {
    alt((
        map(
            pair(delimited(space0, char('-'), space0), factor),
            |(_, f)| Expr::Sub(Box::new(Expr::Float(0.0)), Box::new(f)),
        ),
        delimited(
            tuple((space0, char('('), space0)),
            expr,
            tuple((space0, char(')'), space0)),
        ),
        map(parse_f64, Expr::Float),
        map(identifier, Expr::Var),
    ))(input)
}

pub fn expr(input: &str) -> IResult<&str, Expr> {
    let (input, init) = term(input)?;
    let (input, res) = many0(pair(
        delimited(space0, alt((char('+'), char('-'))), space0),
        term,
    ))(input)?;
    Ok((
        input,
        res.into_iter().fold(init, |acc, (op, val)| match op {
            '+' => Expr::Add(Box::new(acc), Box::new(val)),
            '-' => Expr::Sub(Box::new(acc), Box::new(val)),
            _ => unreachable!(),
        }),
    ))
}

pub fn gate_call(input: &str) -> IResult<&str, ParsedStatement> {
    let (input, name) = identifier(input)?;
    let (input, params) = opt(delimited(
        tuple((space0, char('('), space0)),
        separated_list0(tuple((space0, char(','), space0)), expr),
        tuple((space0, char(')'), space0)),
    ))(input)?;

    let input = if params.is_some() {
        let (input, _) = space0(input)?;
        input
    } else {
        let (input, _) = space1(input)?;
        input
    };

    let (input, qubits) = separated_list0(tuple((space0, char(','), space0)), qubit_ref)(input)?;
    let (input, _) = space0(input)?;
    let (input, _) = tag(";")(input)?;

    Ok((
        input,
        ParsedStatement::Gate(name, qubits, params.unwrap_or_default()),
    ))
}

pub fn measure(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("measure"),
            space1,
            qubit_ref,
            space0,
            tag("->"),
            space0,
            qubit_ref,
            space0,
            tag(";"),
        )),
        |(_, _, q, _, _, _, c, _, _)| ParsedStatement::Measure(q, c),
    )(input)
}

pub fn barrier(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("barrier"),
            space0,
            separated_list0(tuple((space0, char(','), space0)), qubit_ref),
            space0,
            tag(";"),
        )),
        |(_, _, qubits, _, _)| ParsedStatement::Barrier(qubits),
    )(input)
}

/// Parses `reset q[i];` or register-wide `reset q;`, producing a synthetic
/// gate-call statement (the upper parser will expand it into
/// `Operation::Reset` by intercepting the "reset" name).
///
/// We re-use `ParsedStatement::Gate` with the special name `"__reset__"` —
/// but cleaner is to introduce a dedicated parser that the top-level parser
/// dispatches on. Here we emit a `Gate("reset", ...)` call, which the
/// existing parser machinery turns into an unknown gate error. To avoid
/// churning the AST, we instead expose a dedicated helper.
pub fn reset(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("reset"),
            space1,
            separated_list0(tuple((space0, char(','), space0)), qubit_ref),
            space0,
            tag(";"),
        )),
        |(_, _, qubits, _, _)| {
            // Represent reset as a Gate with a sentinel name. The parser
            // `handle_statement` must be aware of this — but to avoid
            // touching the AST enum shape, we piggyback on Gate with a
            // well-known sentinel identifier.
            ParsedStatement::Gate("__reset__".to_string(), qubits, Vec::new())
        },
    )(input)
}

pub fn if_stmt(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("if"),
            space0,
            char('('),
            space0,
            identifier,
            space0,
            tag("=="),
            space0,
            usize_parser,
            space0,
            char(')'),
            space0,
            alt((measure, reset, gate_call, barrier)),
        )),
        |(_, _, _, _, creg, _, _, _, val, _, _, _, op)| {
            ParsedStatement::If(creg, val, Box::new(op))
        },
    )(input)
}

fn gate_body_stmt(input: &str) -> IResult<&str, ParsedStatement> {
    alt((barrier, gate_call))(input)
}

pub fn gate_def(input: &str) -> IResult<&str, ParsedStatement> {
    let (input, _) = tag("gate")(input)?;
    let (input, _) = space1(input)?;
    let (input, name) = identifier(input)?;
    let (input, params) = opt(delimited(
        tuple((space0, char('('), space0)),
        separated_list0(tuple((space0, char(','), space0)), identifier),
        tuple((space0, char(')'), space0)),
    ))(input)?;
    let (input, _) = space0(input)?;
    let (input, qubits) = separated_list0(tuple((space0, char(','), space0)), identifier)(input)?;
    let (input, _) = space0(input)?;

    let (input, body) = delimited(
        tuple((space0, char('{'), space0)),
        many0(delimited(space0, gate_body_stmt, space0)),
        tuple((space0, char('}'), space0)),
    )(input)?;

    Ok((
        input,
        ParsedStatement::GateDef(name, params.unwrap_or_default(), qubits, body),
    ))
}