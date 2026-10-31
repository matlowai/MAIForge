use pf_surface::river_water::{WaterCell, WaterGrid, WaterSource};
fn grid(nx: usize, ny: usize, dx: f64) -> WaterGrid {
    WaterGrid {
        source_identity: pf_contracts::sim_module::digest_bytes(
            b"numeric hydraulic control, not a world",
        ),
        nx,
        ny,
        spacing_m: dx,
        gravity_m_s2: 9.81,
        friction_cf: 0.,
        bed_m: vec![0.; nx * ny],
        water: vec![WaterCell::default(); nx * ny],
        seconds: 0.,
        steps: 0,
    }
}
fn run(s: &mut WaterGrid, until: f64) -> usize {
    let mut cut = 0;
    while s.seconds < until {
        let dt = s.timestep_limit().unwrap().min(until - s.seconds);
        let q = s.step(dt, &[]).unwrap();
        cut += q.moving_wet_faces_truncated_by_bed;
        assert!(q.residual_m3.abs() < 1e-7);
    }
    cut
}
#[test]
fn lake_at_rest_with_islands_and_dry_banks() {
    let mut s = grid(31, 19, 2.);
    for y in 0..s.ny {
        for x in 0..s.nx {
            let k = y * s.nx + x;
            s.bed_m[k] = 0.2 * (x as f64 - 15.).abs() + 0.1 * (y as f64 - 9.).abs() - 1.;
            s.water[k].depth_m = (-s.bed_m[k]).max(0.);
        }
    }
    let original = s.clone();
    run(&mut s, 10.);
    for (a, b) in s.water.iter().zip(&original.water) {
        assert!((a.depth_m - b.depth_m).abs() < 1e-12);
        assert!(a.flux_x_m2_s.abs() + a.flux_y_m2_s.abs() < 1e-12);
    }
    assert!((s.volume_m3() - original.volume_m3()).abs() < 1e-10);
}
fn dam_error(nx: usize) -> f64 {
    let dx = 160. / nx as f64;
    let mut s = grid(nx, 4, dx);
    for y in 0..s.ny {
        for x in 0..nx / 2 {
            s.water[y * nx + x].depth_m = 1.;
        }
    }
    let initial = s.volume_m3();
    run(&mut s, 5.);
    assert!((s.volume_m3() - initial).abs() < 1e-8);
    let c = 9.81f64.sqrt();
    let mut err = 0.;
    for x in 0..nx {
        let position = (x as f64 + 0.5) * dx - 80.;
        let exact = if position <= -c * 5. {
            1.
        } else if position >= 2. * c * 5. {
            0.
        } else {
            (2. * c - position / 5.).powi(2) / (9. * 9.81)
        };
        err += (s.water[nx + x].depth_m - exact).abs() * dx;
        for y in 1..s.ny {
            assert!((s.water[y * nx + x].depth_m - s.water[x].depth_m).abs() < 1e-12);
        }
    }
    err / 160.
}
#[test]
fn dry_dam_break_converges_to_ritter_rarefaction() {
    let coarse = dam_error(160);
    let fine = dam_error(320);
    eprintln!("Ritter mean absolute depth error: coarse={coarse}, fine={fine}");
    assert!(coarse < 0.025);
    assert!(fine < 0.85 * coarse);
}
#[test]
fn sources_and_refusal_are_atomic() {
    let mut s = grid(5, 5, 10.);
    let src = WaterSource {
        cell: 12,
        volume_m3_s: 4.,
        velocity_x_m_s: 1.,
        velocity_y_m_s: 0.,
    };
    let r = s.step(0.1, &[src]).unwrap();
    assert!((r.final_m3 - 0.4).abs() < 1e-12);
    assert_eq!(r.supplied_m3, 0.4);
    let old = s.clone();
    assert!(s.step(1e6, &[src]).is_err());
    assert_eq!(s, old);
    assert!(s.step(0.1, &[WaterSource { cell: 25, ..src }]).is_err());
    assert_eq!(s, old);
    assert!(s
        .step(
            0.1,
            &[WaterSource {
                volume_m3_s: -1.,
                ..src
            }]
        )
        .is_err());
    assert_eq!(s, old);
}
#[test]
fn bed_step_does_not_leak_at_diagonal_wet_corner() {
    let mut s = grid(4, 4, 1.);
    s.bed_m.fill(2.);
    s.bed_m[5] = -1.;
    s.bed_m[10] = -1.;
    s.water[5].depth_m = 1.;
    s.water[10].depth_m = 1.;
    let initial = s.water.clone();
    run(&mut s, 1.);
    assert_eq!(s.water, initial);
}
#[test]
fn drag_dissipates_interior_momentum_without_water_loss() {
    let mut s = grid(9, 9, 10.);
    s.friction_cf = 0.02;
    for c in &mut s.water {
        c.depth_m = 2.;
        c.flux_x_m2_s = 4.;
    }
    let before = s.volume_m3();
    s.step(0.01, &[]).unwrap();
    assert!((s.volume_m3() - before).abs() < 1e-10);
    assert!((s.water[40].flux_x_m2_s - 4. / (1. + 0.01 * 0.02 * 2. / 2.)).abs() < 1e-12);
}
#[test]
fn bad_state_and_timestep_are_refused() {
    let mut s = grid(2, 2, 1.);
    s.water[0].flux_x_m2_s = 1.;
    assert!(s.validate().is_err());
    s.water[0] = WaterCell::default();
    assert!(s.step(f64::NAN, &[]).is_err());
    s.bed_m[0] = f64::NAN;
    assert!(s.validate().is_err());
}

/// Actual game terrain, conditional water experiment; never a delivered world state.
#[test]
#[ignore = "requires retained native detail block and a new lane output directory"]
fn native_outlet_closed_domain_comparison() {
    use serde_json::{json, Value};
    use std::{fs, path::PathBuf, time::Instant};
    let src = PathBuf::from(std::env::var("VZ_RIVER_WATER_DETAIL").unwrap());
    let out = PathBuf::from(std::env::var("VZ_RIVER_WATER_OUT").unwrap());
    assert!(!out.exists());
    let doc: Value = serde_json::from_slice(&fs::read(src.join("block.json")).unwrap()).unwrap();
    assert_eq!(doc["block_size_posts"], json!([257, 257]));
    let bytes = fs::read(src.join("block_height_mm.i32")).unwrap();
    assert_eq!(bytes.len(), 257 * 257 * 4);
    assert_eq!(
        pf_contracts::sim_module::digest_bytes(&bytes),
        doc["block_digests"]["height"].as_str().unwrap()
    );
    let z: Vec<f64> = bytes
        .chunks_exact(4)
        .map(|v| i32::from_le_bytes(v.try_into().unwrap()) as f64 / 1000.)
        .collect();
    let mut s = grid(256, 256, 2.);
    s.source_identity = pf_contracts::sim_module::digest_bytes(&bytes);
    s.friction_cf = 0.008;
    for y in 0..256 {
        for x in 0..256 {
            let i = y * 256 + x;
            s.bed_m[i] = (z[y * 257 + x]
                + z[y * 257 + x + 1]
                + z[(y + 1) * 257 + x]
                + z[(y + 1) * 257 + x + 1])
                * 0.25;
            s.water[i].depth_m = (-s.bed_m[i]).max(0.);
        }
    }
    let initial = s.volume_m3();
    let still = s.clone();
    let t = Instant::now();
    run(&mut s, 2.);
    let rest_seconds = t.elapsed().as_secs_f64();
    let rest_error = s
        .water
        .iter()
        .zip(&still.water)
        .map(|(a, b)| (a.depth_m - b.depth_m).abs())
        .fold(0., f64::max);
    assert!(rest_error < 1e-10);
    // Deliberate documented impulse, not native river forcing or sediment displacement.
    let k = s.water.iter().position(|c| c.depth_m > 0.8).unwrap();
    s.water[k].depth_m += 0.05;
    let injected = 0.05 * 4.;
    let mut resumed = s.clone();
    let t = Instant::now();
    let cuts = run(&mut s, 12.);
    let seconds = t.elapsed().as_secs_f64();
    let final_volume = s.volume_m3();
    assert!((initial + injected - final_volume).abs() < 1e-7);
    let changed = s
        .water
        .iter()
        .zip(&still.water)
        .filter(|(a, b)| (a.depth_m - b.depth_m).abs() > 1e-8)
        .count();
    assert!(changed > 1);
    for _ in 0..50 {
        let dt = resumed.timestep_limit().unwrap().min(12. - resumed.seconds);
        resumed.step(dt, &[]).unwrap();
    }
    let checkpoint = resumed.encode_state().unwrap();
    resumed = WaterGrid::decode_state(&checkpoint, &s.source_identity).unwrap();
    run(&mut resumed, 12.);
    assert_eq!(s.encode_state().unwrap(), resumed.encode_state().unwrap());
    fs::create_dir(&out).unwrap();
    fs::write(out.join("checkpoint.rwc"), checkpoint).unwrap();
    for (name, values) in [
        (
            "depth_m.f64",
            s.water.iter().map(|c| c.depth_m).collect::<Vec<_>>(),
        ),
        (
            "flux_x_m2_s.f64",
            s.water.iter().map(|c| c.flux_x_m2_s).collect(),
        ),
        (
            "flux_y_m2_s.f64",
            s.water.iter().map(|c| c.flux_y_m2_s).collect(),
        ),
    ] {
        fs::write(
            out.join(name),
            values
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect::<Vec<_>>(),
        )
        .unwrap();
    }
    let result = json!({"source_detail":src,"source_height_digest":pf_contracts::sim_module::digest_bytes(&bytes),"source_recipe_digest":doc["recipe_digest"],"cells":65536,"spacing_m":2,"model":"first-order hydrostatic reconstruction/Rusanov, reflective walls", "water_initialization":"conditional sea-stage0 in every subzero depression; not accepted actual water state", "geometry":"2m cell bed = average of four native ground posts; no native terrain edits", "rest_physical_seconds":2,"rest_max_depth_error_m":rest_error,"rest_wall_seconds":rest_seconds,"perturbation_volume_m3":injected,"perturbation_cell":k,"perturbation_physical_seconds":10,"perturbation_wall_seconds":seconds,"total_steps":s.steps,"actual_source_checkpoint_resume_exact":true,"changed_depth_cells":changed,"moving_wet_faces_truncated_by_bed":cuts,"water_residual_m3":initial+injected-final_volume,"world_mutated":false,"boundary_accepted":false,"river_evolution_accepted":false});
    fs::write(
        out.join("RESULT.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    eprintln!("{result}");
}

#[test]
fn checkpoint_resume_is_exact_and_corruption_or_other_source_refuses() {
    let mut a = grid(19, 19, 2.);
    for c in &mut a.water {
        c.depth_m = 1.;
    }
    a.water[180].depth_m = 1.1;
    let mut b = a.clone();
    for _ in 0..12 {
        a.step(a.timestep_limit().unwrap(), &[]).unwrap();
    }
    for _ in 0..5 {
        b.step(b.timestep_limit().unwrap(), &[]).unwrap();
    }
    let bytes = b.encode_state().unwrap();
    let source = b.source_identity.clone();
    b = WaterGrid::decode_state(&bytes, &source).unwrap();
    for _ in 0..7 {
        b.step(b.timestep_limit().unwrap(), &[]).unwrap();
    }
    assert_eq!(a.encode_state().unwrap(), b.encode_state().unwrap());
    let mut bad = bytes.clone();
    let k = bad.len() - 1;
    bad[k] ^= 1;
    assert!(WaterGrid::decode_state(&bad, &source).is_err());
    assert!(WaterGrid::decode_state(
        &bytes,
        &pf_contracts::sim_module::digest_bytes(b"other seed")
    )
    .is_err());
    assert!(WaterGrid::decode_state(&bytes[..bytes.len() - 1], &source).is_err());
}
#[test]
fn two_dimensional_pulse_preserves_axis_symmetry() {
    let mut s = grid(21, 21, 1.);
    for c in &mut s.water {
        c.depth_m = 1.;
    }
    s.water[10 * 21 + 10].depth_m = 1.2;
    run(&mut s, 0.4);
    for y in 0..21 {
        for x in 0..21 {
            let a = s.water[y * 21 + x];
            let b = s.water[x * 21 + y];
            assert!((a.depth_m - b.depth_m).abs() < 1e-12);
            assert!((a.flux_x_m2_s - b.flux_y_m2_s).abs() < 1e-12);
        }
    }
}

#[test]
fn open_boundaries_account_for_uniform_through_flow_per_reservoir() {
    use pf_surface::river_water::{Edge, WaterBoundary};
    let mut s = grid(8, 4, 3.);
    for c in &mut s.water {
        c.depth_m = 2.;
        c.flux_x_m2_s = 2.;
    }
    let before = s.water.clone();
    let west = pf_contracts::sim_module::digest_bytes(b"upstream reservoir");
    let east = pf_contracts::sim_module::digest_bytes(b"downstream reservoir");
    let mut boundaries = Vec::new();
    for index in 0..4 {
        for (edge, id) in [(Edge::West, &west), (Edge::East, &east)] {
            boundaries.push(WaterBoundary {
                edge,
                index,
                reservoir_identity: id.clone(),
                exterior_bed_m: 0.,
                exterior_water: WaterCell {
                    depth_m: 2.,
                    flux_x_m2_s: 2.,
                    flux_y_m2_s: 0.,
                },
            });
        }
    }
    let report = s.step_open(0.01, &[], &boundaries).unwrap();
    assert_eq!(s.water, before);
    assert!((report.boundary_in_m3 - 0.24).abs() < 1e-12);
    assert!((report.boundary_out_m3 - 0.24).abs() < 1e-12);
    assert!(report.residual_m3.abs() < 1e-12);
    assert_eq!(report.boundary_volumes.len(), 2);
    let upstream = report
        .boundary_volumes
        .iter()
        .find(|r| r.reservoir_identity == west)
        .unwrap();
    assert!((upstream.incoming_m3 - 0.24).abs() < 1e-12);
    assert_eq!(upstream.outgoing_m3, 0.);
}
#[test]
fn exterior_wetting_and_cfl_are_conservative() {
    use pf_surface::river_water::{Edge, WaterBoundary};
    let mut s = grid(8, 4, 1.);
    let b = WaterBoundary {
        edge: Edge::West,
        index: 1,
        reservoir_identity: pf_contracts::sim_module::digest_bytes(b"wet reservoir"),
        exterior_bed_m: 0.,
        exterior_water: WaterCell {
            depth_m: 1.,
            ..WaterCell::default()
        },
    };
    let original = s.clone();
    assert!(s.step_open(1., &[], &[b.clone()]).is_err());
    assert_eq!(s, original);
    let dt = s.timestep_limit_open(&[b.clone()]).unwrap();
    let report = s.step_open(dt, &[], &[b.clone()]).unwrap();
    assert!(report.boundary_in_m3 > 0.);
    assert!((s.volume_m3() - report.boundary_in_m3).abs() < 1e-12);
    let original = s.clone();
    assert!(s.step_open(dt, &[], &[b.clone(), b.clone()]).is_err());
    assert_eq!(s, original);
    let mut bad = b;
    bad.index = 4;
    assert!(s.timestep_limit_open(&[bad]).is_err());
}
#[test]
fn open_reservoir_at_equal_stage_preserves_rest() {
    use pf_surface::river_water::{Edge, WaterBoundary};
    let mut s = grid(8, 4, 1.);
    for c in &mut s.water {
        c.depth_m = 1.;
    }
    let boundaries: Vec<_> = (0..4)
        .map(|index| WaterBoundary {
            edge: Edge::West,
            index,
            reservoir_identity: pf_contracts::sim_module::digest_bytes(b"sea"),
            exterior_bed_m: -1.,
            exterior_water: WaterCell {
                depth_m: 2.,
                ..WaterCell::default()
            },
        })
        .collect();
    let initial = s.clone();
    let q = s.step_open(0.01, &[], &boundaries).unwrap();
    assert_eq!(s.water, initial.water);
    assert_eq!(q.boundary_in_m3, 0.);
    assert_eq!(q.boundary_out_m3, 0.);
}
