use std::process::Command;
#[test]
fn disabled_key_never_reads_inputs() {
    for args in [
        vec![],
        vec![
            "coast=0",
            "windows=/nonexistent-pf-coast-input",
            "out=/nonexistent-pf-coast-output",
        ],
    ] {
        let o = Command::new(env!("CARGO_BIN_EXE_pf-coast"))
            .args(args)
            .output()
            .unwrap();
        assert!(o.status.success());
        assert!(String::from_utf8(o.stdout)
            .unwrap()
            .contains("no input reads or writes"));
    }
}
#[test]
fn unknown_duplicate_and_missing_tide_fail_loudly() {
    for args in [
        vec!["coast=2"],
        vec!["coast=1", "typo=1"],
        vec!["coast=0", "coast=1"],
        vec!["coast=1"],
    ] {
        let o = Command::new(env!("CARGO_BIN_EXE_pf-coast"))
            .args(args)
            .output()
            .unwrap();
        assert!(!o.status.success());
    }
}
