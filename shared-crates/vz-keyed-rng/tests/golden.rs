//! Golden vectors: in-process output, a separate process (run twice), and the
//! checked-in file must agree byte for byte.
use std::process::Command;

const FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/golden/keyed_rng_v1.golden.json");

fn file() -> String { std::fs::read_to_string(FILE).expect("golden file") }

#[test]
fn in_process_matches_golden_file() {
    assert_eq!(vz_keyed_rng::golden::golden_bytes(), file());
}

#[test]
fn separate_processes_are_bit_equal_to_golden_file() {
    let want = file();
    for run in 0..2 {
        let out = Command::new(env!("CARGO_BIN_EXE_vz-keyed-rng-golden")).output().expect("spawn");
        assert!(out.status.success(), "run {run}");
        assert_eq!(String::from_utf8(out.stdout).unwrap(), want, "process run {run} differs");
    }
}

#[test]
fn golden_pins_version_and_table() {
    let f = file();
    assert!(f.contains("\"keyed_rng\":\"vz.keyed-rng/v1\""));
    assert!(f.contains(&format!("\"digest\":\"{}\"", vz_keyed_rng::normal_table_digest())));
    assert_eq!(vz_keyed_rng::normal_table_lens(), (4095, 4097));
}
