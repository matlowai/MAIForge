use pf_surface::{
    river_bed::{capacity_per_width_m2_s, Gravel},
    river_bend::Hydraulics,
    river_flow::FlowLaw,
    river_reach::*,
};
fn reach(curvature: f64) -> FixedBankReach {
    let slope = 0.008 / (9.81 * 3.0);
    let h = Hydraulics {
        discharge_m3_s: 240.0,
        width_m: 80.0,
        depth_m: 3.0,
        friction_cf: 0.008,
        scour_factor: 3.0,
        gravity_m_s2: 9.81,
        water_density_kg_m3: 1000.0,
    };
    let gravel = Gravel {
        diameter_m: 0.002,
        solid_density_kg_m3: 2650.0,
        porosity: 0.35,
        grain_shear_fraction: 1.0,
    };
    FixedBankReach {
        source_identity: "engineering-reference:fixed-bank-v0".into(),
        flow: FlowLaw {
            discharge_m3_s: 240.0,
            friction_cf: 0.008,
            gravity_m_s2: 9.81,
            density_kg_m3: 1000.0,
        },
        gravel,
        downstream_stage_m: 3.0,
        upstream_solid_m3_s: 80.0 * capacity_per_width_m2_s(h, gravel, 8.0).unwrap(),
        slope_factor: 0.88,
        von_karman: 0.4,
        rows: (0..21)
            .map(|i| BedRow {
                distance_m: 50.0 * i as f64,
                strip_width_m: 5.0,
                curvature_per_m: curvature,
                basement_m: vec![slope * (1000.0 - 50.0 * i as f64) - 1.0; 16],
                mobile_m: vec![1.0; 16],
            })
            .collect(),
        sim_seconds: 0.0,
        steps: 0,
        cumulative_in_solid_m3: 0.0,
        cumulative_out_solid_m3: 0.0,
        pending_boundary_water_m3: 0.0,
    }
}
#[test]
fn coupled_steps_change_bed_then_recompute_flow_with_solid_closure() {
    let mut r = reach(0.0001);
    let before = r.hydraulics().unwrap();
    for _ in 0..4 {
        let receipt = r.step(3600.0).unwrap();
        assert!(receipt.solid_residual_m3.abs() < 1e-7);
    }
    let after = r.hydraulics().unwrap();
    assert!(r.rows[20].mobile_m[0] > 1.0);
    assert!(r.rows[20].mobile_m[15] < 1.0);
    assert!(after[20].velocities_m_s[0] < before[20].velocities_m_s[0]);
    for f in after {
        assert!((f.discharges_m3_s.iter().sum::<f64>() - 240.0).abs() < 1e-9);
    }
    assert_eq!(r.sim_seconds, 14400.0);
    assert_eq!(r.steps, 4);
}
#[test]
fn checkpoint_resume_matches_the_same_native_schedule_exactly() {
    let mut uninterrupted = reach(0.0001);
    let mut split = uninterrupted.clone();
    for _ in 0..8 {
        uninterrupted.step(1800.0).unwrap();
    }
    for _ in 0..3 {
        split.step(1800.0).unwrap();
    }
    let checkpoint = split.encode_state().unwrap();
    let mut restored = FixedBankReach::decode_state(&checkpoint, &split.source_identity).unwrap();
    assert_eq!(restored.encode_state().unwrap(), checkpoint);
    for _ in 3..8 {
        restored.step(1800.0).unwrap();
    }
    assert_eq!(
        uninterrupted.encode_state().unwrap(),
        restored.encode_state().unwrap()
    );
    assert!(FixedBankReach::decode_state(&checkpoint, "wrong-world").is_err());
    assert!(FixedBankReach::decode_state(
        &checkpoint[..checkpoint.len() - 1],
        &split.source_identity
    )
    .is_err());
}
#[test]
fn rejected_steps_are_atomic_and_zero_flow_does_not_create_sediment() {
    let mut r = reach(0.0001);
    let before = r.encode_state().unwrap();
    assert!(r.step(1e12).is_err());
    assert_eq!(r.encode_state().unwrap(), before);
    assert!(r.step(f64::NAN).is_err());
    assert_eq!(r.encode_state().unwrap(), before);
    r.flow.discharge_m3_s = 0.0;
    r.upstream_solid_m3_s = 0.0;
    let initial = r.rows.clone();
    let receipt = r.step(3600.0).unwrap();
    assert_eq!(receipt.solid_in_m3, 0.0);
    assert_eq!(receipt.solid_out_m3, 0.0);
    for (a, b) in initial.iter().zip(&r.rows) {
        assert_eq!(a.mobile_m, b.mobile_m);
    }
}

#[test]
fn coupled_timestep_refinement_reduces_bed_error() {
    let simulate = |dt: f64| {
        let mut r = reach(0.001);
        for _ in 0..(3600.0 / dt) as usize {
            r.step(dt).unwrap();
        }
        r
    };
    let a = simulate(900.0);
    let b = simulate(450.0);
    let c = simulate(225.0);
    let difference = |a: &FixedBankReach, b: &FixedBankReach| {
        a.rows
            .iter()
            .zip(&b.rows)
            .flat_map(|(a, b)| {
                a.mobile_m
                    .iter()
                    .zip(&b.mobile_m)
                    .map(|(a, b)| (a - b).powi(2))
            })
            .sum::<f64>()
            .sqrt()
    };
    let ab = difference(&a, &b);
    let bc = difference(&b, &c);
    assert!(ab > 0.0 && ab / bc > 1.7, "coupled refinement {ab} {bc}");
}

#[test]
#[ignore = "explicit actual-world fixed-bank coupled-state probe"]
fn native_reach_checkpoint_probe() {
    use serde_json::{json, Value};
    use std::{fs, path::Path, time::Instant};
    let path = std::env::var("VZ_RIVER_REACH_REQUEST").unwrap();
    let request: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let source: Value =
        serde_json::from_slice(&fs::read(request["flow_report"].as_str().unwrap()).unwrap())
            .unwrap();
    let root = Path::new(request["world_root"].as_str().unwrap());
    let status: Value =
        serde_json::from_slice(&fs::read(root.join("status.json")).unwrap()).unwrap();
    assert_eq!(source["source_state"], status["state_digest"]);
    let z = fs::read(root.join("terrain/fields/elevation_m.f32")).unwrap();
    let q = fs::read(root.join("hydrology/fields/discharge_proxy_m3_s.f32")).unwrap();
    let mut cases = Vec::new();
    for c in source["results"].as_array().unwrap() {
        if c["status"] != "flow_and_erosion_proposal_only" {
            continue;
        }
        let h: Hydraulics = serde_json::from_value(c["hydraulics"].clone()).unwrap();
        let cell = c["cell"].as_u64().unwrap() as usize;
        let z0 = f32::from_le_bytes(z[4 * cell..4 * cell + 4].try_into().unwrap()) as f64;
        assert_eq!(
            h.discharge_m3_s,
            f32::from_le_bytes(q[4 * cell..4 * cell + 4].try_into().unwrap()) as f64
        );
        let slope = c["native_coarse_slope"].as_f64().unwrap();
        let mut r = reach(0.0);
        r.source_identity = format!("{}:cell{cell}", status["state_digest"].as_str().unwrap());
        r.flow.discharge_m3_s = h.discharge_m3_s;
        r.downstream_stage_m = z0 - slope * 20000.0 + h.depth_m;
        r.upstream_solid_m3_s = h.width_m
            * capacity_per_width_m2_s(
                h,
                r.gravel,
                h.water_density_kg_m3 * h.friction_cf * h.speed_m_s().powi(2),
            )
            .unwrap();
        r.rows = (0..101)
            .map(|i| BedRow {
                distance_m: 200.0 * i as f64,
                strip_width_m: 5.0,
                curvature_per_m: 0.0,
                basement_m: vec![z0 - slope * 200.0 * i as f64 - 0.5; 16],
                mobile_m: vec![0.5; 16],
            })
            .collect();
        let mut full = r.clone();
        let start = Instant::now();
        let mut residual = 0.0_f64;
        for _ in 0..8 {
            residual = residual.max(full.step(1800.0).unwrap().solid_residual_m3.abs());
        }
        let elapsed = start.elapsed().as_secs_f64();
        for _ in 0..3 {
            r.step(1800.0).unwrap();
        }
        let checkpoint = r.encode_state().unwrap();
        let mut resumed = FixedBankReach::decode_state(&checkpoint, &r.source_identity).unwrap();
        for _ in 3..8 {
            resumed.step(1800.0).unwrap();
        }
        assert_eq!(
            full.encode_state().unwrap(),
            resumed.encode_state().unwrap()
        );
        cases.push(json!({"cell":cell,"source_identity":r.source_identity,"steps":8,"physical_seconds":full.sim_seconds,
            "sections":101,"strips":16,"full_schedule_wall_seconds":elapsed,"checkpoint_bytes":checkpoint.len(),
            "resume_byte_exact":true,"max_solid_residual_m3":residual,"pending_boundary_water_m3":full.pending_boundary_water_m3,
            "assumptions":"fixed banks, straight coarse reach,80m width,0.5m mobile uniform2mm gravel; equilibrium inlet sediment supply; not fine native alluvial calibration"}));
    }
    let output = Path::new(request["output"].as_str().unwrap());
    assert!(!output.exists());
    fs::write(output,serde_json::to_vec_pretty(&json!({"cases":cases,"fixed_bank_control_only":true,"full_meander_history_complete":false,"game_integration":false})).unwrap()).unwrap();
}

#[test]
fn coastal_refinement_replaces_local_runoff_without_duplicating_parent_flow() {
    use pf_surface::river_reach::{partition_parent_runoff, RegionalRunoffArea};
    let p = partition_parent_runoff(
        1000.,
        30.,
        0.01,
        &[
            RegionalRunoffArea {
                region_xy: [0, 0],
                river_area_m2: 100.,
                direct_marine_area_m2: 200.,
            },
            RegionalRunoffArea {
                region_xy: [-1, 0],
                river_area_m2: 150.,
                direct_marine_area_m2: 50.,
            },
        ],
    )
    .unwrap();
    assert_eq!(p.upstream_m3_s, 20.);
    assert_eq!(p.unrefined_local_m3_s, 5.);
    assert_eq!(p.refined_river_m3_s, 2.5);
    assert_eq!(p.direct_marine_m3_s, 2.5);
    assert_eq!(p.river_out_m3_s, 27.5);
    assert_eq!(p.residual_m3_s, 0.);
    let all = partition_parent_runoff(
        1000.,
        30.,
        0.01,
        &[RegionalRunoffArea {
            region_xy: [0, 0],
            river_area_m2: 0.,
            direct_marine_area_m2: 1000.,
        }],
    )
    .unwrap();
    assert_eq!(
        all.river_out_m3_s, 20.,
        "upstream river remains even when local cell drains straight to coast"
    );
    let unchanged = partition_parent_runoff(1000., 30., 0.01, &[]).unwrap();
    assert_eq!(unchanged.river_out_m3_s, 30.);
}
#[test]
fn refinement_partition_rejects_duplicate_coverage_and_unfunded_parent_flow() {
    use pf_surface::river_reach::{partition_parent_runoff, RegionalRunoffArea};
    let r = RegionalRunoffArea {
        region_xy: [0, 0],
        river_area_m2: 200.,
        direct_marine_area_m2: 200.,
    };
    assert!(partition_parent_runoff(1000., 30., 0.01, &[r.clone(), r.clone()]).is_err());
    assert!(partition_parent_runoff(300., 30., 0.01, &[r.clone()]).is_err());
    assert!(partition_parent_runoff(1000., 9., 0.01, &[r.clone()]).is_err());
    assert!(partition_parent_runoff(1000., 30., f64::NAN, &[r.clone()]).is_err());
    let mut invalid = r;
    invalid.direct_marine_area_m2 = -1.;
    assert!(partition_parent_runoff(1000., 30., 0.01, &[invalid]).is_err());
}
