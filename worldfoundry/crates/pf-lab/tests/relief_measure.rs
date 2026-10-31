//! SEED-RELIEF-LAW (lane RL, 2026-10-12): ignored measurement driver, a copy
//! of lane OG's `orogen_measure.rs` with the relief-law modes added (OG's
//! modes kept so both are measured by the same binary). Writes
//! raw fields and metadata for one (seed, fidelity, mode, dt) run; the lane
//! analysis scripts compute the rubric rows from them. Not part of the
//! normal suite. Env: OG_SEED, OG_FIDELITY (preview | standard), OG_MODE,
//! OG_OUT, optional OG_DT (step length, My).
//!
//! Lane RL modes: `orogen_relief` (orogen_v1 with stream-power-subgrid-
//! nonlinear-v1 in place of stream-power-subgrid-v1, PreviewOptions::
//! orogen_relief_v1), `relief_only` (bodies + the relief law, legacy
//! thickening), and LABELLED SENSITIVITY modes `orogen_relief_kref<v>` /
//! `orogen_relief_n<v>` (never a proposal or a tuning step).
//!
//! Modes: `bodies` (old laws: two-column + continental bodies),
//! `orogen` (bodies + rate-limited-belt-v1 + stream-power-subgrid-v1),
//! `orogen_belt` (bodies + belt only), `orogen_subgrid` (bodies + subgrid
//! only), `legacy` (default chain), and the DIAGNOSTIC counterfactuals
//! `orogen_noerosion` / `bodies_noerosion` (K = D = 0; never a proposal).

use pf_contracts::params::{
    FidelityPreset, NonlinearSubgridErosionBlock, OrogenyBlock, ScenarioPreset, SurfaceBlock,
};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::kinematic::CrustState;
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewOptions};
use pf_lab::preview_module::outputs_digest;
use pf_surface::water::largest_ocean;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn write_f32(path: &Path, v: &[f64]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    for x in v {
        f.write_all(&(*x as f32).to_le_bytes()).unwrap();
    }
}

fn options(mode: &str) -> PreviewOptions {
    let no_erosion = || {
        let mut s = SurfaceBlock::default();
        s.erodibility = 0.0;
        s.hillslope_diffusivity_m2_per_yr = 0.0;
        s
    };
    match mode {
        "legacy" => PreviewOptions::legacy(),
        "bodies" => PreviewOptions::crust_bodies_v0(),
        "orogen" => PreviewOptions::orogen_v1(),
        "orogen_belt" => PreviewOptions {
            orogeny: Some(OrogenyBlock {
                subgrid_erosion: None,
                ..OrogenyBlock::v1()
            }),
            ..PreviewOptions::crust_bodies_v0()
        },
        "orogen_subgrid" => PreviewOptions {
            orogeny: Some(OrogenyBlock {
                rate_limited_belt: None,
                ..OrogenyBlock::v1()
            }),
            ..PreviewOptions::crust_bodies_v0()
        },
        // SENSITIVITY ONLY (labelled; never a proposal or a tuning step):
        // the pre-registered v1 profile with a different channel viscosity.
        m if m.starts_with("orogen_mu") => {
            let mu: f64 = m
                .trim_start_matches("orogen_mu")
                .parse()
                .expect("orogen_mu<value>");
            let mut o = OrogenyBlock::v1();
            o.rate_limited_belt
                .as_mut()
                .unwrap()
                .lower_crust_channel_viscosity_pa_s = mu;
            PreviewOptions {
                orogeny: Some(o),
                ..PreviewOptions::crust_bodies_v0()
            }
        }
        "orogen_relief" => PreviewOptions::orogen_relief_v1(),
        "relief_only" => PreviewOptions {
            orogeny: Some(OrogenyBlock {
                rate_limited_belt: None,
                ..OrogenyBlock::relief_v1()
            }),
            ..PreviewOptions::crust_bodies_v0()
        },
        // SENSITIVITY ONLY (labelled; never a proposal or a tuning step).
        m if m.starts_with("orogen_relief_kref") || m.starts_with("orogen_relief_n") => {
            let mut b = NonlinearSubgridErosionBlock::v1();
            if let Some(v) = m.strip_prefix("orogen_relief_kref") {
                b.reference_steepness_m = v.parse().expect("orogen_relief_kref<value>");
            } else {
                let v = m.strip_prefix("orogen_relief_n").unwrap();
                b.slope_exponent_n = v.parse().expect("orogen_relief_n<value>");
            }
            PreviewOptions {
                orogeny: Some(OrogenyBlock {
                    nonlinear_subgrid_erosion: Some(b),
                    ..OrogenyBlock::relief_v1()
                }),
                ..PreviewOptions::crust_bodies_v0()
            }
        }
        "orogen_relief_noerosion" => PreviewOptions {
            surface_override: Some(no_erosion()),
            ..PreviewOptions::orogen_relief_v1()
        },
        "orogen_noerosion" => PreviewOptions {
            surface_override: Some(no_erosion()),
            ..PreviewOptions::orogen_v1()
        },
        "bodies_noerosion" => PreviewOptions {
            surface_override: Some(no_erosion()),
            ..PreviewOptions::crust_bodies_v0()
        },
        other => panic!("unknown mode {other}"),
    }
}

#[test]
#[ignore]
fn relief_measure() {
    let seed: u64 = std::env::var("OG_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("OG_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("OG_MODE").unwrap();
    let dt: Option<f64> = std::env::var("OG_DT").ok().map(|v| v.parse().unwrap());
    let opts = options(&mode);
    let out_dir = PathBuf::from(std::env::var("OG_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my: dt,
    };
    let t = Instant::now();
    let out = match simulate_with(&cfg, &opts) {
        Ok(o) => o,
        Err(e) => {
            std::fs::write(out_dir.join("ERROR.txt"), format!("{e}")).unwrap();
            return;
        }
    };
    let secs = t.elapsed().as_secs_f64();
    let grid = &out.grid;
    let sea = out.params.planet.sea_level_m;
    let seeding = seed_plates(grid, &out.params, &cfg.seed_path()).unwrap();
    let prior = UpliftPrior::preview_v0();
    let state0 = match &out.params.continental_bodies {
        None => CrustState::seed(grid, &seeding, &prior).unwrap(),
        Some(b) => {
            CrustState::seed_bodies(
                grid,
                &seeding,
                &prior,
                b.resolved_area_fraction(&out.params.tectonics),
                b.rifted_margin_width_m,
            )
            .unwrap()
            .0
        }
    };
    write_f32(&out_dir.join("elevation_t0_m.f32"), &out.elevation_t0_m);
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_t0_m.f32"), &state0.thickness_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    write_f32(
        &out_dir.join("eroded_total_m.f32"),
        &out.state.eroded_total_m,
    );
    write_f32(
        &out_dir.join("shortening_added_m.f32"),
        &out.state.shortening_added_m,
    );
    write_f32(&out_dir.join("arc_added_m.f32"), &out.state.arc_added_m);
    write_f32(
        &out_dir.join("uplift_rate_m_per_my.f32"),
        &out.uplift_rate_m_per_my,
    );
    if let Some(r) = out
        .orogen
        .as_ref()
        .and_then(|o| o.subgrid_relief_m.as_ref())
    {
        write_f32(&out_dir.join("subgrid_relief_m.f32"), r);
    }
    let code = |k: &[CrustKindCodeV0]| k.iter().map(|x| x.code()).collect::<Vec<u8>>();
    std::fs::write(out_dir.join("kind.u8"), code(&out.state.kind)).unwrap();
    std::fs::write(out_dir.join("kind_t0.u8"), code(&state0.kind)).unwrap();
    std::fs::write(out_dir.join("boundary_class.u8"), &out.state.boundary_class).unwrap();
    let plate: Vec<u8> = out
        .state
        .plate
        .iter()
        .flat_map(|p| p.0.to_le_bytes())
        .collect();
    std::fs::write(out_dir.join("plate_id.u32"), plate).unwrap();
    let ocean = largest_ocean(grid, &out.elevation_m, sea).unwrap();
    std::fs::write(
        out_dir.join("ocean_connected.u8"),
        ocean.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    )
    .unwrap();
    let ocean0 = largest_ocean(grid, &out.elevation_t0_m, sea).unwrap();
    std::fs::write(
        out_dir.join("ocean_connected_t0.u8"),
        ocean0.iter().map(|&b| b as u8).collect::<Vec<_>>(),
    )
    .unwrap();
    // Contributing area and flow receiver of the final routing (drainage
    // rubric rows; R5 is not this node's exit, recorded for context).
    write_f32(
        &out_dir.join("contributing_area_m2.f32"),
        &out.routing.contributing_area_m2,
    );
    let row_areas: Vec<f64> = (0..grid.nlat())
        .map(|i| grid.cell_area_m2(i).unwrap())
        .collect();
    let digest = outputs_digest(&out).unwrap();
    let meta = serde_json::json!({
        "seed": seed,
        "mode": mode,
        "dt_my": out.params.tectonics.step_my,
        "nlat": grid.nlat(),
        "nlon": grid.nlon(),
        "row_area_m2": row_areas,
        "total_area_m2": row_areas.iter().sum::<f64>() * grid.nlon() as f64,
        "sea_level_datum_m": sea,
        "ocean_fraction_param": out.params.tectonics.ocean_fraction,
        "simulate_seconds": secs,
        "outputs_digest": digest,
        "sea_level_t0": out.sea_level_t0,
        "steps": out.steps,
        "orogen_steps": out.orogen.as_ref().map(|o| &o.steps),
        "summary": out.terrain.summary,
        "crust_volume": out.terrain.crust_volume,
        "sediment_balance": out.terrain.sediment_balance,
        "approximations": out.terrain.approximations,
        "assumptions": out.terrain.assumptions,
        "uplift_prior": out.terrain.uplift_prior,
        "erosion": out.terrain.erosion,
        "params": out.params,
    });
    std::fs::write(
        out_dir.join("meta.json"),
        serde_json::to_vec_pretty(&meta).unwrap(),
    )
    .unwrap();
}
