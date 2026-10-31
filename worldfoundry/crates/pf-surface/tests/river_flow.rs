use pf_surface::river_flow::*;
fn law() -> FlowLaw {
    FlowLaw {
        discharge_m3_s: 240.0,
        friction_cf: 0.008,
        gravity_m_s2: 9.81,
        density_kg_m3: 1000.0,
    }
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t, "{a} != {b}");
}
fn channel(n: usize) -> Vec<Section> {
    let slope = 0.008 / (9.81 * 3.0); // Q240, W80, H3 -> U1
    (0..n)
        .map(|i| {
            let x = 1000.0 * i as f64 / (n - 1) as f64;
            Section {
                distance_m: x,
                strip_width_m: 5.0,
                bed_m: vec![slope * (1000.0 - x); 16],
            }
        })
        .collect()
}
#[test]
fn normal_flow_profile_matches_analytic_depth_and_discharge() {
    let s = channel(51);
    let p = profile(law(), &s, 3.0).unwrap();
    for (s, f) in s.iter().zip(p) {
        close(f.water_surface_m - s.bed_m[0], 3.0, 1e-10);
        close(f.discharges_m3_s.iter().sum(), 240.0, 1e-10);
        for u in f.velocities_m_s {
            close(u, 1.0, 1e-10);
        }
    }
}
#[test]
fn bed_shape_changes_flow_distribution_and_energy_closes() {
    let mut s = channel(21);
    for row in &mut s {
        for j in 0..8 {
            row.bed_m[j] += 0.1;
        }
    }
    let p = profile(law(), &s, 3.0).unwrap();
    for f in &p {
        assert!(f.velocities_m_s[0] < f.velocities_m_s[15]);
        close(f.discharges_m3_s.iter().sum(), 240.0, 1e-10);
    }
    for i in 0..p.len() - 1 {
        close(
            p[i].energy_head_m - p[i + 1].energy_head_m,
            0.5 * (p[i].friction_slope + p[i + 1].friction_slope)
                * (s[i + 1].distance_m - s[i].distance_m),
            1e-10,
        );
    }
}
#[test]
fn downstream_stage_drives_backwater_and_static_water_is_level() {
    let s = channel(21);
    let p = profile(law(), &s, 3.4).unwrap();
    assert!(p[0].water_surface_m > s[0].bed_m[0] + 3.0);
    assert!(p[0].depths_m[0] < 3.4);
    let mut dry = law();
    dry.discharge_m3_s = 0.0;
    for f in profile(dry, &s, 3.0).unwrap() {
        assert_eq!(f.water_surface_m, 3.0);
        assert!(f.discharges_m3_s.iter().all(|q| *q == 0.0));
    }
}
#[test]
fn refuses_dry_critical_or_invalid_boundaries() {
    let s = channel(21);
    assert!(profile(law(), &s, 0.0).is_err());
    assert!(profile(law(), &s, 0.1).is_err());
    assert!(profile(law(), &s, f64::NAN).is_err());
    let mut bad = s.clone();
    bad[2].distance_m = bad[1].distance_m;
    assert!(profile(law(), &bad, 3.0).is_err());
}

#[test]
fn standard_step_converges_against_independent_gvf_integration() {
    // Independent RK4 integration of dh/dx=(S0-Sf)/(1-Fr²), backwards
    // from downstream h3.4m in the rectangular constant-width case.
    let l = law();
    let s0 = 0.008 / (9.81 * 3.0);
    let derivative = |depth: f64| {
        let fr2 = l.discharge_m3_s.powi(2) / (l.gravity_m_s2 * 80.0_f64.powi(2) * depth.powi(3));
        (s0 - l.friction_cf * fr2) / (1.0 - fr2)
    };
    let mut reference = 3.4;
    let dx = -1.0;
    for _ in 0..1000 {
        let a = derivative(reference);
        let b = derivative(reference + dx * a / 2.0);
        let c = derivative(reference + dx * b / 2.0);
        let d = derivative(reference + dx * c);
        reference += dx * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for n in [11, 21, 41] {
        let s = channel(n);
        let p = profile(l, &s, 3.4).unwrap();
        errors.push((p[0].depths_m[0] - reference).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
}

#[test]
#[ignore = "explicit source-bound bank/bed to hydraulic-feedback probe"]
fn native_bed_to_flow_feedback() {
    use pf_surface::{
        river_bed::{self, BedCell, Gravel},
        river_bend::{BankMaterial, Hydraulics, Station},
    };
    use serde_json::{json, Value};
    use std::{fs, path::Path, time::Instant};
    let input = std::env::var("VZ_RIVER_FLOW_REQUEST").unwrap();
    let request: Value = serde_json::from_slice(&fs::read(input).unwrap()).unwrap();
    let flow_report: Value =
        serde_json::from_slice(&fs::read(request["flow_report"].as_str().unwrap()).unwrap())
            .unwrap();
    let root = Path::new(request["world_root"].as_str().unwrap());
    let native_z = fs::read(root.join("terrain/fields/elevation_m.f32")).unwrap();
    let status: Value =
        serde_json::from_slice(&fs::read(root.join("status.json")).unwrap()).unwrap();
    assert_eq!(status["state_digest"], flow_report["source_state"]);
    let mut cases = Vec::new();
    for c in flow_report["results"].as_array().unwrap() {
        if c["status"] != "flow_and_erosion_proposal_only" {
            continue;
        }
        let hydro: Hydraulics = serde_json::from_value(c["hydraulics"].clone()).unwrap();
        let cell = c["cell"].as_u64().unwrap() as usize;
        let z_up = f32::from_le_bytes(native_z[4 * cell..4 * cell + 4].try_into().unwrap()) as f64;
        let slope = c["native_coarse_slope"].as_f64().unwrap();
        let stations: Vec<_> = (0..1001)
            .map(|i| Station {
                distance_m: 20.0 * i as f64,
                curvature_per_m: 0.0,
            })
            .collect();
        let bed: Vec<_> = (0..1000)
            .map(|i| BedCell {
                bed_elevation_m: z_up - slope * (20.0 * i as f64 + 10.0),
                mobile_thickness_m: 0.0,
            })
            .collect();
        let banks = vec![
            BankMaterial {
                erodibility_m_pa_s: 1e-9,
                critical_shear_pa: 4.0,
                bank_height_m: hydro.depth_m + 1.0,
                porosity: 0.35
            };
            1001
        ];
        let sediment = river_bed::propose_step(
            hydro,
            &stations,
            0.0,
            &banks,
            &banks,
            Gravel {
                diameter_m: 0.029,
                solid_density_kg_m3: 2650.0,
                porosity: 0.35,
                grain_shear_fraction: 1.0,
            },
            &bed,
            0.0,
            86400.0,
        )
        .unwrap();
        let sections = |bed: &[BedCell]| {
            bed.iter()
                .enumerate()
                .map(|(i, c)| Section {
                    distance_m: 20.0 * i as f64 + 10.0,
                    strip_width_m: 5.0,
                    bed_m: vec![c.bed_elevation_m; 16],
                })
                .collect::<Vec<_>>()
        };
        let before = sections(&bed);
        let after = sections(&sediment.next_bed);
        let boundary = bed.last().unwrap().bed_elevation_m + hydro.depth_m;
        let fluid = FlowLaw {
            discharge_m3_s: hydro.discharge_m3_s,
            friction_cf: hydro.friction_cf,
            gravity_m_s2: hydro.gravity_m_s2,
            density_kg_m3: hydro.water_density_kg_m3,
        };
        let start = Instant::now();
        let old = profile(fluid, &before, boundary).unwrap();
        let new = profile(fluid, &after, boundary).unwrap();
        let wall = start.elapsed().as_secs_f64();
        let volume = |rows: &[SectionFlow]| {
            rows.windows(2)
                .map(|p| 0.5 * (p[0].area_m2 + p[1].area_m2) * 20.0)
                .sum::<f64>()
        };
        let delta_u = old
            .iter()
            .zip(&new)
            .map(|(a, b)| (a.velocities_m_s[0] - b.velocities_m_s[0]).abs())
            .fold(0.0, f64::max);
        assert!(delta_u > 0.0);
        for row in &new {
            close(row.discharges_m3_s.iter().sum(), hydro.discharge_m3_s, 1e-9);
        }
        cases.push(json!({"cell":cell,"bank_solid_m3":sediment.bank.removed_solid_m3,"bed_solid_residual_m3":sediment.residual_solid_m3,
            "profile_pair_wall_seconds":wall,"max_velocity_response_m_s":delta_u,"water_storage_change_interior_m3":volume(&new)-volume(&old),
            "hydraulic_span_m":19980.0,"sections":1000,"strips":16,"all_section_Q_checks_pass":true,
            "assumptions":"native coarse grade/Q;80mwide; declared29mm uniform gravel, bank properties and normal downstream stage; one-day candidate, fixed banks"}));
    }
    let out = Path::new(request["output"].as_str().unwrap());
    assert!(!out.exists());
    fs::write(out,serde_json::to_vec_pretty(&json!({"source_state":status["state_digest"],"cases":cases,
        "water_storage_transfer_not_committed":true,"bank_geometry_unmodified":true,"recent_history_complete":false,"game_integration":false})).unwrap()).unwrap();
}

#[test]
fn connected_dry_margins_do_not_receive_water_or_discharge() {
    let mut s = channel(21);
    for row in &mut s {
        let bank = row.bed_m[0] + 5.0;
        let mut bed = vec![bank; 4];
        bed.extend(&row.bed_m);
        bed.extend(vec![bank; 4]);
        row.bed_m = bed;
    }
    let p = profile(law(), &s, 3.0).unwrap();
    for f in p {
        assert_eq!(f.wetted_width_m, 80.0);
        close(f.area_m2, 240.0, 1e-8);
        assert!(f.discharges_m3_s[..4].iter().all(|q| *q == 0.0));
        assert!(f.discharges_m3_s[20..].iter().all(|q| *q == 0.0));
        close(f.discharges_m3_s.iter().sum(), 240.0, 1e-10);
    }
}
#[test]
fn disconnected_channels_are_refused_instead_of_bridging_a_dry_bar() {
    let mut s = channel(2);
    for row in &mut s {
        row.bed_m[8] += 5.0;
    }
    assert!(profile(law(), &s, 3.0)
        .unwrap_err()
        .contains("disconnected"));
}
