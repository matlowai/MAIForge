//! Golden vectors: in-process output, a separate process (run twice), and the
//! checked-in file (which the JS twin also reproduces) agree byte for byte.
use std::process::Command;
const FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/golden/tilekey_v2.golden.json");

#[test]
fn in_process_matches_golden_file() {
    assert_eq!(vz_tilekey::golden::golden_bytes(), std::fs::read_to_string(FILE).unwrap());
}

#[test]
fn separate_processes_are_bit_equal_to_golden_file() {
    let want = std::fs::read_to_string(FILE).unwrap();
    for run in 0..2 {
        let out = Command::new(env!("CARGO_BIN_EXE_vz-tilekey-golden")).output().unwrap();
        assert!(out.status.success());
        assert_eq!(String::from_utf8(out.stdout).unwrap(), want, "run {run}");
    }
}

#[test]
fn ladder_declaration_is_canonical_and_matches_golden() {
    let s = vz_tilekey::lod_ladder_v1_canonical_json();
    assert_eq!(s, vz_tilekey::golden::to_string(&vz_tilekey::golden::ladder_j()));
    let v: serde_json::Value = serde_json::from_str(&s).unwrap();
    assert_eq!(v["named_lods"][2]["tile_side_mm"], 16_384_000);
    assert_eq!(v["named_lods"][4]["post_spacing_mm"], 2_048_000);
}
