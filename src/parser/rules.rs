use super::ast::ParsedStatement;
use nom::{
    IResult,
    branch::alt,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::{alpha1, alphanumeric1, char, digit1, space0, space1},
    combinator::{map, map_res, opt, recognize, value},
    multi::{many0, separated_list0},
    sequence::{delimited, pair, tuple},
};

// --- Helpers ---

fn identifier(input: &str) -> IResult<&str, String> {
    map(
        recognize(pair(
            alt((alpha1, tag("_"))),
            many0(alt((alphanumeric1, tag("_")))),
        )),
        |s: &str| s.to_string(),
    )(input)
}

fn usize_parser(input: &str) -> IResult<&str, usize> {
    map_res(digit1, |s: &str| s.parse::<usize>())(input)
}

fn float_parser(input: &str) -> IResult<&str, f64> {
    map_res(
        recognize(tuple((
            opt(tag("-")),
            digit1,
            opt(tuple((tag("."), digit1))),
        ))),
        |s: &str| s.parse::<f64>(),
    )(input)
}

pub fn comment(input: &str) -> IResult<&str, ()> {
    value((), pair(tag("//"), take_while(|c| c != '\n')))(input)
}

// --- QASM Parsers ---

pub fn openqasm_version(input: &str) -> IResult<&str, String> {
    map(
        tuple((
            tag("OPENQASM"),
            space1,
            take_while1(|c: char| c != ';'),
            tag(";"),
        )),
        |(_, _, version, _): (&str, &str, &str, &str)| version.to_string(),
    )(input)
}

pub fn include(input: &str) -> IResult<&str, ParsedStatement> {
    value(
        ParsedStatement::Ignore,
        tuple((
            tag("include"),
            space1,
            delimited(char('"'), take_while1(|c| c != '"'), char('"')),
            tag(";"),
        )),
    )(input)
}

pub fn qreg(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("qreg"),
            space1,
            identifier,
            delimited(char('['), usize_parser, char(']')),
            tag(";"),
        )),
        |(_, _, name, size, _)| ParsedStatement::QReg(name, size),
    )(input)
}

pub fn creg(input: &str) -> IResult<&str, ParsedStatement> {
    map(
        tuple((
            tag("creg"),
            space1,
            identifier,
            delimited(char('['), usize_parser, char(']')),
            tag(";"),
        )),
        |(_, _, name, size, _)| ParsedStatement::CReg(name, size),
    )(input)
}

fn qubit_ref(input: &str) -> IResult<&str, (String, usize)> {
    pair(identifier, delimited(char('['), usize_parser, char(']')))(input)
}

pub fn gate_call(input: &str) -> IResult<&str, ParsedStatement> {
    let (input, name) = identifier(input)?;
    let (input, params) = opt(delimited(
        char('('),
        separated_list0(tuple((space0, char(','), space0)), float_parser),
        char(')'),
    ))(input)?;
    let (input, _) = space1(input)?;
    let (input, qubits) = separated_list0(tuple((space0, char(','), space0)), qubit_ref)(input)?;
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
            tag(";"),
        )),
        |(_, _, q, _, _, _, c, _)| ParsedStatement::Measure(q, c),
    )(input)
}
