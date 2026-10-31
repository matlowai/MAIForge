//! Prints the vz.keyed-rng/v1 golden document (canonical JSON + newline).
//! `--bench N` instead times N plain draws and N normal draws (not golden).

use std::io::Write;
use vz_keyed_rng::Scope;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 3 && args[1] == "--bench" {
        let n: u64 = args[2].parse().expect("N");
        let sc = Scope::root(1).law_set("bench").unwrap().chart(0).unwrap().tile(0, 0, 0).unwrap();
        let t = std::time::Instant::now();
        let mut acc = 0u64;
        for c in 0..n { acc ^= sc.draw("bench.u", c, 0).unwrap().u64(); }
        let du = t.elapsed();
        let t = std::time::Instant::now();
        for c in 0..n { acc ^= sc.draw("bench.n", c, 0).unwrap().normal_q32() as u64; }
        let dn = t.elapsed();
        let t = std::time::Instant::now();
        for c in 0..n { acc ^= sc.address(&format!("a/{c}")).unwrap().draw("bench.a", 0, 0).unwrap().u64(); }
        let da = t.elapsed();
        println!("draws={n} u64_ns_per_draw={:.1} normal_ns_per_draw={:.1} address_plus_draw_ns={:.1} xor={acc:016x}",
            du.as_nanos() as f64 / n as f64, dn.as_nanos() as f64 / n as f64, da.as_nanos() as f64 / n as f64);
        return;
    }
    let out = vz_keyed_rng::golden::golden_bytes();
    std::io::stdout().write_all(out.as_bytes()).unwrap();
}
