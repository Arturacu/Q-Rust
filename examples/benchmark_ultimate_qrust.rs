use std::process::Command;

fn run_example(name: &str) {
    println!("\n=============================================");
    println!("Running: {}", name);
    println!("=============================================");
    let status = Command::new("cargo")
        .args(["run", "--release", "--example", name])
        .status()
        .expect("Failed to execute example");

    if !status.success() {
        eprintln!("{} failed!", name);
    }
}

fn main() {
    println!("Starting Quantum Decathlon...");
    run_example("benchmark_qrust");
    run_example("benchmark_chain_qrust");
    run_example("benchmark_structure_qrust");
    run_example("benchmark_massive_qrust");
    println!("\nDecathlon Complete!");
}