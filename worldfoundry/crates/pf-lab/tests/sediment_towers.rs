//! SEDIMENT-TOWERS (lane ST, 2026-10-12): ignored per-step trace driver.
//! Steps one run (seed, fidelity, mode) and writes a per-step CSV of the
//! sediment ledger (advance: belt / stacked / subducted / residual; erosion:
//! eroded, land deposition, export) and of the thick-sediment population
//! (max, counts above 10 / 20 km, how many thick cells stand in a closed
//! basin on the post-step surface). Writes the final fields like
//! `carry_measure.rs`. Not part of the normal suite.
//! Env: ST_SEED, ST_FIDELITY (preview | standard), ST_MODE, ST_OUT.
//! Modes: orogen, orogen_carry, orogen_carry_noerosion (DIAGNOSTIC: K = D
//! = 0, never a proposal), bodies, legacy, and the corrected profiles
//! orogen_v2 / orogen_carry_v2 (`accreted-to-belt-crust-v1`).
//!
//! The normal test runs the corrected profiles on lane OG's mini grid
//! (36 x 72, 20 My in 2 My steps) through the genesis module (every step's
//! control-volume closures are checked there) and the batch path: module
//! replay equals batch, the crust ledger closes, and three runs agree.

use pf_contracts::params::{FidelityPreset, ScenarioPreset, SurfaceBlock};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_contracts::sim_module::{ModuleRun, RunBinding};
use pf_lab::preview::{simulate_with, PreviewConfig, PreviewEnv, PreviewOptions, PreviewRun};
use pf_lab::preview_module::{outputs_digest, preview_recipe, PreviewGenesisModule as M};
use pf_surface::routing::route;
use pf_surface::water::largest_ocean;
use std::io::Write;
use std::path::{Path, PathBuf};

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
        "orogen_carry" => PreviewOptions::orogen_carry_v1(),
        "orogen_carry_v2" => PreviewOptions::orogen_carry_v2(),
        "orogen_v2" => PreviewOptions::orogen_v2(),
        "orogen_carry_noerosion" => PreviewOptions {
            surface_override: Some(no_erosion()),
            ..PreviewOptions::orogen_carry_v1()
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
fn sediment_tower_trace() {
    let seed: u64 = std::env::var("ST_SEED").unwrap().parse().unwrap();
    let fidelity = match std::env::var("ST_FIDELITY").unwrap().as_str() {
        "preview" => FidelityPreset::Preview,
        "standard" => FidelityPreset::Standard,
        other => panic!("unknown fidelity {other}"),
    };
    let mode = std::env::var("ST_MODE").unwrap();
    let out_dir = PathBuf::from(std::env::var("ST_OUT").unwrap());
    std::fs::create_dir_all(&out_dir).unwrap();
    let cfg = PreviewConfig {
        scenario: ScenarioPreset::toy(),
        fidelity,
        grid: None,
        seed,
        duration_my: None,
        step_my: None,
    };
    let params = cfg.profile_with(&options(&mode)).unwrap();
    let env = PreviewEnv::new(params, cfg.grid, cfg.seed_path()).unwrap();
    let mut run = PreviewRun::init(env).unwrap();
    let grid = run.env().grid().clone();
    let areas = run.env().areas().to_vec();
    let n = grid.cells();
    let sea = run.env().params().planet.sea_level_m;
    let rain = vec![1.0f64; n];
    let mut csv = std::io::BufWriter::new(std::fs::File::create(out_dir.join("steps.csv")).unwrap());
    writeln!(
        csv,
        "t_my,sed_before_m3,sed_after_adv_m3,sed_final_m3,belt_sed_m3,stacked_m3,subducted_m3,adv_residual_m3,eroded_m3,sed_eroded_m3,deposited_m3,exported_m3,max_sed_m,max_sed_cont_m,n_gt10km,n_gt20km,n_gt10km_wet,max_sed_cell,max_sed_cell_elev,max_sed_cell_fill,max_sed_cell_thick,max_sed_cell_kind,max_elev_m,sed_at_max_elev_m,accreted_m3,moved_max_m,moved_max_cell,sed_at_moved_cell,elev_at_moved_cell"
    )
    .unwrap();
    while run.steps_taken() < run.steps_scheduled() {
        let r = run.step().unwrap();
        let st = run.state();
        let z = run.elevation_m();
        let ocean = largest_ocean(&grid, z, sea).unwrap();
        let routing = route(&grid, z, &rain, &ocean, sea).unwrap();
        let wet = |c: usize| !ocean[c] && routing.fill_level_m[c] > z[c] + 1.0e-9;
        let (mut mx, mut mxc, mut ic) = (0.0f64, 0.0f64, 0usize);
        let (mut g10, mut g20, mut g10w) = (0u64, 0u64, 0u64);
        for c in 0..n {
            let s = st.sediment_m[c];
            if s > mx {
                mx = s;
                ic = c;
            }
            if st.kind[c] == CrustKindCodeV0::Continental {
                mxc = mxc.max(s);
            }
            if s > 1.0e4 {
                g10 += 1;
                if wet(c) {
                    g10w += 1;
                }
            }
            if s > 2.0e4 {
                g20 += 1;
            }
        }
        let iz = (0..n).max_by(|&a, &b| z[a].total_cmp(&z[b])).unwrap();
        let a = &r.advance;
        writeln!(
            csv,
            "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.1},{:.1},{},{},{},{},{:.1},{:.1},{:.1},{},{:.1},{:.1},{:.6e},{:.1},{},{:.1},{:.1}",
            r.t_my,
            r.sediment_before_m3,
            r.sediment_after_advance_m3,
            r.sediment_final_m3,
            a.belt_sediment_m3,
            a.sediment_stacked_m3,
            a.sediment_subducted_m3,
            a.sediment_residual_m3(),
            r.eroded_m3,
            r.sediment_eroded_m3,
            r.deposited_m3,
            r.exported_m3,
            mx,
            mxc,
            g10,
            g20,
            g10w,
            ic,
            z[ic],
            routing.fill_level_m[ic],
            st.thickness_m[ic],
            st.kind[ic] as u8,
            z[iz],
            st.sediment_m[iz],
            a.sediment_accreted_m3,
            a.sediment_moved_max_m,
            a.sediment_moved_max_cell,
            st.sediment_m[a.sediment_moved_max_cell as usize],
            z[a.sediment_moved_max_cell as usize],
        )
        .unwrap();
    }
    csv.flush().unwrap();
    let out = run.finish().unwrap();
    write_f32(&out_dir.join("elevation_m.f32"), &out.elevation_m);
    write_f32(&out_dir.join("thickness_m.f32"), &out.state.thickness_m);
    write_f32(&out_dir.join("sediment_m.f32"), &out.state.sediment_m);
    write_f32(&out_dir.join("eroded_total_m.f32"), &out.state.eroded_total_m);
    write_f32(&out_dir.join("shortening_added_m.f32"), &out.state.shortening_added_m);
    let kind: Vec<u8> = out.state.kind.iter().map(|k| *k as u8).collect();
    std::fs::write(out_dir.join("kind.u8"), kind).unwrap();
    let meta = serde_json::json!({
        "seed": seed, "mode": mode, "nlat": grid.nlat(), "nlon": grid.nlon(),
        "row_area_m2": (0..grid.nlat()).map(|i| areas[i * grid.nlon()]).collect::<Vec<_>>(),
    });
    std::fs::write(out_dir.join("meta.json"), meta.to_string()).unwrap();
    let _ = &areas;
}

#[test]
fn accreted_profiles_close_replay_and_repeat() {
    for (seed, opts) in [
        (240907u64, PreviewOptions::orogen_carry_v2()),
        (3, PreviewOptions::orogen_v2()),
    ] {
        let cfg = PreviewConfig {
            scenario: ScenarioPreset::toy(),
            fidelity: FidelityPreset::Preview,
            grid: Some((36, 72)),
            seed,
            duration_my: Some(20.0),
            step_my: Some(2.0),
        };
        let a = simulate_with(&cfg, &opts).unwrap();
        let da = outputs_digest(&a).unwrap();
        let b = simulate_with(&cfg, &opts).unwrap();
        assert_eq!(outputs_digest(&b).unwrap(), da, "seed {seed}");
        let cv = &a.terrain.crust_volume;
        assert!(
            cv.closure_plug_m3.abs() <= 1e-9 * cv.initial_m3,
            "seed {seed}: closure plug {} of {}",
            cv.closure_plug_m3,
            cv.initial_m3
        );
        assert!(a
            .terrain
            .assumptions
            .iter()
            .any(|x| x.id == "A-COLLISION-SEDIMENT-ACCRETED"));
        // The law acts on this mini world (accretion happened).
        let env = PreviewEnv::new(cfg.profile_with(&opts).unwrap(), cfg.grid, cfg.seed_path()).unwrap();
        let mut stepper = PreviewRun::init(env).unwrap();
        let mut accreted = 0.0;
        while stepper.steps_taken() < stepper.steps_scheduled() {
            accreted += stepper.step().unwrap().advance.sediment_accreted_m3;
        }
        let stepped = stepper.finish().unwrap();
        assert_eq!(outputs_digest(&stepped).unwrap(), da, "seed {seed}");
        eprintln!("seed {seed}: accreted {accreted:.3e} m3");
        // Only consumed losers (front-rescue-carry-v1) accrete; under the
        // kinematic-flux supply the loser keeps its sediment (no copy).
        if seed == 240907 {
            assert!(accreted > 0.0);
        } else {
            assert_eq!(accreted, 0.0);
        }
        let recipe = preview_recipe(&cfg, &opts, false).unwrap();
        let mut run = ModuleRun::<M>::start(recipe, RunBinding::GENESIS_SIM_TIME).unwrap();
        run.run_to_end(None).unwrap();
        let out = M::finish(run.into_state()).unwrap();
        assert_eq!(outputs_digest(&out).unwrap(), da, "seed {seed}: module != batch");
    }
}
