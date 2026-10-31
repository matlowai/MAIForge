//! Native artifact adapter: observed C7 snapshots -> historical climate book.
//! Usage: cargo run -p pf-surface --example bed_climate_archive -- WORLD SNAPSHOTS OUT.json
//! No simulation advance, stock mutation, authored place or alternate climate law.
use pf_contracts::{
    artifact::{digest_canonical, verify_artifact, ArtifactOnDisk},
    terrain::bed_history::{HistoricalClimate, CLIMATE_SCHEMA, LABEL},
};
use pf_grid::GridSpec;
use pf_surface::{bed_climate::reconstruct_frame, climate::ClimateDomainV0};
use serde_json::Value;
use std::{
    fs,
    io::{BufReader, BufWriter, Write},
    path::PathBuf,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: bed_climate_archive WORLD SNAPSHOTS OUT.json".into());
    }
    let world = PathBuf::from(&args[1]);
    let status: Value =
        serde_json::from_reader(BufReader::new(fs::File::open(world.join("status.json"))?))?;
    let artifact: ArtifactOnDisk<ClimateDomainV0> = serde_json::from_reader(BufReader::new(
        fs::File::open(world.join("climate/artifact.json"))?,
    ))?;
    let verified = verify_artifact::<ClimateDomainV0>(&world.join("climate"))?;
    if status["climate_digest"] != verified {
        return Err("world/climate artifact mismatch".into());
    }
    let domain = artifact.header.domain;
    if domain.recipe.digest()? != domain.recipe_digest {
        return Err("climate recipe digest mismatch".into());
    }
    let recipe: Value =
        serde_json::from_reader(BufReader::new(fs::File::open(world.join("recipe.json"))?))?;
    if digest_canonical(&recipe)? != status["recipe_digest"] {
        return Err("world recipe/status digest mismatch".into());
    }
    // C7's archive uses the fixed output datum after the volume-primary solve.
    let profile = &recipe["inputs"]["profile"];
    let radius = profile["planet"]["radius_m"]
        .as_f64()
        .ok_or("planet radius")?;
    let expected_sea = profile["pinned_universe"]["constants"]["sea_level_datum_m"]["value"]
        .as_f64()
        .ok_or("pinned output datum")?;
    if domain.inputs.sea_level_m != expected_sea {
        return Err("snapshot datum differs from climate datum".into());
    }
    let mut snapshots = Vec::new();
    for entry in fs::read_dir(&args[2])? {
        let p = entry?.path();
        if !p.is_dir() {
            continue;
        }
        let meta: Value =
            serde_json::from_reader(BufReader::new(fs::File::open(p.join("meta.json"))?))?;
        snapshots.push((meta["t_my"].as_f64().ok_or("snapshot time")?, p, meta));
    }
    snapshots.sort_by(|a, b| a.0.total_cmp(&b.0));
    let first = &snapshots.first().ok_or("no retained snapshots")?.2;
    let nlat = first["nlat"].as_u64().ok_or("nlat")? as usize;
    let nlon = first["nlon"].as_u64().ok_or("nlon")? as usize;
    let grid = GridSpec::new(nlat, nlon, radius)?;
    let mut frames = Vec::new();
    for (time, path, meta) in snapshots {
        if meta["nlat"] != nlat
            || meta["nlon"] != nlon
            || meta["seed"] != recipe["seed_path"]["seed"]
        {
            return Err("snapshot grid/seed mismatch".into());
        }
        let bytes = fs::read(path.join("elevation_m.f32"))?;
        if bytes.len() != grid.cells() * 4 {
            return Err("snapshot elevation length".into());
        }
        let z: Vec<f64> = bytes
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes(b.try_into().unwrap()) as f64)
            .collect();
        let pin = digest_canonical(&z)?;
        frames.push(reconstruct_frame(
            &grid,
            time,
            &z,
            domain.inputs.sea_level_m,
            &domain.recipe,
            pin,
        )?);
        eprintln!(
            "reconstructed climate frame t={time} My ({} cells)",
            grid.cells()
        );
    }
    // Pin the final snapshot to the retained world elevation. Earlier driver
    // snapshots have no original signed manifest; their identity limitation
    // is preserved in assumptions instead of fabricating a receipt.
    let final_bytes = fs::read(world.join("terrain/fields/elevation_m.f32"))?;
    let final_z: Vec<f64> = final_bytes
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes(b.try_into().unwrap()) as f64)
        .collect();
    if frames.last().ok_or("no frames")?.elevation_source_digest != digest_canonical(&final_z)? {
        return Err("final snapshot does not match retained world terrain".into());
    }
    let history=HistoricalClimate {
        schema:CLIMATE_SCHEMA.into(),label:LABEL.into(),
        physical_recipe_digest:status["recipe_digest"].as_str().ok_or("recipe digest")?.into(),
        climate_recipe:domain.recipe,climate_recipe_digest:domain.recipe_digest,
        nlat,nlon,sea_level_m:domain.inputs.sea_level_m,frames,
        assumptions:vec!["Same planet climate law and fixed output sea datum applied to retained past terrain; never fed back into C7 deposition".into(),
            "Nearest time snapshot inside its recorded span only; no extrapolation before the first frame or after the last".into(),
            "Nearest grid cell at each recorded carrier position; climate topography is coarse and can differ from subcell native facies".into(),
            "Snapshots are unmanifested C7 driver outputs: seed/grid/time checked; each elevation vector is pinned here. No original per-snapshot hash receipt exists".into()],
    };
    history.validate()?;
    let mut out = BufWriter::new(
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&args[3])?,
    );
    serde_json::to_writer(&mut out, &history)?;
    out.flush()?;
    Ok(())
}
