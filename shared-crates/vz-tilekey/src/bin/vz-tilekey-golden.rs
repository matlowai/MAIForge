//! Prints the vz.tilekey/v2 golden document (canonical JSON + newline).
use std::io::Write;
fn main() {
    std::io::stdout().write_all(vz_tilekey::golden::golden_bytes().as_bytes()).unwrap();
}
