//! [E2E-NEW-FEATURE] Smoke tests for the `qrust` CLI binary.
//!
//! These tests are `#[ignore]` by default — they invoke `cargo run`
//! recursively, which is fragile in some CI environments. Run with
//! `cargo test -- --ignored` to exercise the full CLI surface.

use std::process::Command;

fn cargo_run(args: &[&str]) -> (i32, String, String) {
    let mut cmd = Command::new(env!("CARGO"));
    cmd.args(["run", "--quiet", "--bin", "qrust", "--"]);
    cmd.args(args);
    let out = cmd.output().expect("cargo run");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
    )
}

#[test]
#[ignore]
fn test_cli_help() {
    let (code, stdout, _) = cargo_run(&["--help"]);
    assert_eq!(code, 0);
    assert!(stdout.contains("Q-Rust"));
    assert!(stdout.contains("USAGE"));
}

#[test]
#[ignore]
fn test_cli_transpile_bell_state() {
    let (code, stdout, stderr) = cargo_run(&[
        "tests/fixtures/bell_state.qasm",
        "--opt",
        "1",
        "--no-decompose",
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    assert!(stdout.contains("OPENQASM 2.0;"), "stdout: {stdout}");
    assert!(stdout.contains("qreg q[2];"));
}

#[test]
#[ignore]
fn test_cli_with_report() {
    let (code, _stdout, stderr) = cargo_run(&[
        "tests/fixtures/ghz_3.qasm",
        "--opt",
        "2",
        "--backend",
        "linear-3",
        "--report",
    ]);
    assert_eq!(code, 0, "stderr: {stderr}");
    assert!(stderr.contains("Transpilation Report"));
    assert!(stderr.contains("parsed"));
}
