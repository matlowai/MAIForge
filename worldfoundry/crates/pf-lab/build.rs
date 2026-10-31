// Records the rustc version into the binary for ExecutionReceipt.rustc
// (a real measurement, never a guess). Runs offline.
use std::process::Command;

fn main() {
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    let version = Command::new(rustc)
        .arg("--version")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "rustc (version unavailable)".into());
    println!("cargo:rustc-env=PF_LAB_RUSTC_VERSION={version}");
    println!("cargo:rerun-if-changed=build.rs");
}
