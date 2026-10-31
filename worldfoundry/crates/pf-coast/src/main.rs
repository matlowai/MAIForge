use std::{collections::BTreeMap, path::Path};
fn main() {
    if let Err(e) = run() {
        eprintln!("pf-coast: {e}");
        std::process::exit(1);
    }
}
fn run() -> Result<(), String> {
    let mut args = BTreeMap::new();
    for a in std::env::args().skip(1) {
        let (k, v) = a.split_once('=').ok_or(
            "usage: pf-coast coast=1 windows=DIR t0=DIR climate=DIR out=NEW_DIR tide_range_m=2",
        )?;
        if !["coast", "windows", "t0", "climate", "out", "tide_range_m"].contains(&k)
            || args.insert(k.to_string(), v.to_string()).is_some()
        {
            return Err(format!("unknown/duplicate key {k}"));
        }
    }
    match args.get("coast").map(String::as_str).unwrap_or("0") {
        "0" => {
            println!("coast off: no input reads or writes");
            return Ok(());
        }
        "1" => (),
        _ => return Err("coast must be 0 or 1".into()),
    }
    let arg = |k: &str| {
        args.get(k)
            .map(String::as_str)
            .ok_or(format!("missing {k}"))
    };
    let config =
        pf_coast::Config::new(arg("tide_range_m")?.parse().map_err(|_| "bad tide range")?)?;
    let result = pf_coast::publish::run(
        Path::new(arg("windows")?),
        Path::new(arg("t0")?),
        Path::new(arg("climate")?),
        Path::new(arg("out")?),
        &config,
    )?;
    println!("{}", serde_json::to_string(&result).unwrap());
    Ok(())
}
