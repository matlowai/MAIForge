use pf_surface::{
    river_bed::*,
    river_bend::{BankMaterial, Hydraulics, Station},
};

fn h() -> Hydraulics {
    Hydraulics {
        discharge_m3_s: 240.0,
        width_m: 80.0,
        depth_m: 3.0,
        friction_cf: 0.008,
        scour_factor: 3.0,
        gravity_m_s2: 9.81,
        water_density_kg_m3: 1000.0,
    }
}
fn gravel(d: f64) -> Gravel {
    Gravel {
        diameter_m: d,
        solid_density_kg_m3: 2650.0,
        porosity: 0.35,
        grain_shear_fraction: 1.0,
    }
}
fn bank(k: f64) -> BankMaterial {
    BankMaterial {
        erodibility_m_pa_s: k,
        critical_shear_pa: 4.0,
        bank_height_m: 4.0,
        porosity: 0.35,
    }
}
fn stations(n: usize) -> Vec<Station> {
    (0..n)
        .map(|i| Station {
            distance_m: 1000.0 * i as f64 / (n - 1) as f64,
            curvature_per_m: 0.0,
        })
        .collect()
}
fn bed(n: usize, thickness: f64) -> Vec<BedCell> {
    vec![
        BedCell {
            bed_elevation_m: 100.0,
            mobile_thickness_m: thickness
        };
        n - 1
    ]
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{a} != {b}");
}

#[test]
fn capacity_matches_published_equation_24_and_skin_threshold() {
    let mut flow = h();
    flow.gravity_m_s2 = 10.0;
    let g = Gravel {
        diameter_m: 0.01,
        solid_density_kg_m3: 2000.0,
        porosity: 0.4,
        grain_shear_fraction: 1.0,
    };
    // Shields=0.0595, excess=0.01; volumetric capacity, not bulk or kg/s.
    close(
        capacity_per_width_m2_s(flow, g, 5.95).unwrap(),
        0.000012554242310868467,
        1e-18,
    );
    assert_eq!(capacity_per_width_m2_s(flow, g, 4.94).unwrap(), 0.0);
    assert_eq!(
        capacity_per_width_m2_s(
            flow,
            Gravel {
                grain_shear_fraction: 0.5,
                ..g
            },
            5.95
        )
        .unwrap(),
        0.0
    );
}
#[test]
fn below_grain_threshold_bank_material_deposits_with_correct_porosity() {
    let s = stations(11);
    let b = vec![bank(1e-9); 11];
    let initial = bed(11, 0.0);
    let p = propose_step(h(), &s, 0.0, &b, &b, gravel(0.029), &initial, 0.0, 86400.0).unwrap();
    assert!(p.bank.removed_solid_m3 > 0.0);
    assert_eq!(p.downstream_solid_m3, 0.0);
    close(p.final_bed_solid_m3, p.bank.removed_solid_m3, 1e-12);
    assert_eq!(p.bank.unrouted_solid_m3, 0.0);
    assert_eq!(
        p.bank.segment_solid_m3.iter().sum::<f64>(),
        p.bank.removed_solid_m3
    );
    // Two banks, rate=4e-9m/s; full bank height4m; deposit over80m bed.
    let expected_dz = 2.0 * 4e-9 * 86400.0 * 4.0 / 80.0;
    for c in &p.next_bed {
        close(c.bed_elevation_m - 100.0, expected_dz, 1e-13);
    }
}
#[test]
fn available_bed_bounds_export_and_no_negative_reservoir() {
    let s = stations(11);
    let b = vec![bank(0.0); 11];
    let initial = bed(11, 0.001);
    let p = propose_step(h(), &s, 0.0, &b, &b, gravel(0.002), &initial, 0.0, 7200.0).unwrap();
    close(p.initial_bed_solid_m3, 52.0, 1e-10);
    close(p.downstream_solid_m3, 52.0, 1e-10);
    assert!(p.next_bed.iter().all(|b| b.mobile_thickness_m == 0.0));
    close(p.residual_solid_m3, 0.0, 1e-10);
    assert_eq!(initial, bed(11, 0.001)); // candidate did not mutate input
}
#[test]
fn equilibrium_flux_preserves_bed_and_supersaturation_deposits() {
    let s = stations(11);
    let b = vec![bank(0.0); 11];
    let initial = bed(11, 0.1);
    let g = gravel(0.006);
    let cap = 80.0 * capacity_per_width_m2_s(h(), g, 8.0).unwrap();
    let p = propose_step(h(), &s, 0.0, &b, &b, g, &initial, cap, 3600.0).unwrap();
    for c in &p.next_bed {
        close(c.bed_elevation_m, 100.0, 1e-12);
    }
    let more = propose_step(h(), &s, 0.0, &b, &b, g, &initial, cap + 1e-4, 3600.0).unwrap();
    close(more.deposited_solid_m3[0], 0.36, 1e-12);
    close(more.deposited_solid_m3.iter().sum::<f64>(), 0.36, 1e-12);
}
#[test]
fn fixed_flow_split_timestep_and_spacing_preserve_uniform_deposition() {
    let mut heights = Vec::new();
    for n in [11, 21, 41] {
        let s = stations(n);
        let b = vec![bank(1e-9); n];
        let initial = bed(n, 0.0);
        let full =
            propose_step(h(), &s, 0.0, &b, &b, gravel(0.029), &initial, 0.0, 86400.0).unwrap();
        let half =
            propose_step(h(), &s, 0.0, &b, &b, gravel(0.029), &initial, 0.0, 43200.0).unwrap();
        let half2 = propose_step(
            h(),
            &s,
            0.0,
            &b,
            &b,
            gravel(0.029),
            &half.next_bed,
            0.0,
            43200.0,
        )
        .unwrap();
        close(full.final_bed_solid_m3, half2.final_bed_solid_m3, 1e-12);
        for (a, b) in full.next_bed.iter().zip(half2.next_bed) {
            close(a.bed_elevation_m, b.bed_elevation_m, 2e-14);
        }
        heights.push(full.next_bed[0].bed_elevation_m);
    }
    close(heights[0], heights[2], 1e-12);
}
#[test]
fn dry_and_invalid_candidates_never_mutate_bed() {
    let s = stations(11);
    let b = vec![bank(1e-9); 11];
    let initial = bed(11, 0.1);
    let mut dry = h();
    dry.discharge_m3_s = 0.0;
    let p = propose_step(dry, &s, 0.0, &b, &b, gravel(0.006), &initial, 0.0, 86400.0).unwrap();
    assert_eq!(p.next_bed, initial);
    assert!(propose_step(dry, &s, 0.0, &b, &b, gravel(0.006), &initial, 0.1, 86400.0).is_err());
    assert!(propose_step(h(), &s, 0.0, &b, &b, gravel(0.00001), &initial, 0.0, 1.0).is_err());
    assert!(propose_step(
        h(),
        &s,
        0.0,
        &b,
        &b,
        gravel(0.006),
        &initial,
        100.0,
        86400.0
    )
    .is_err());
    assert!(propose_step(h(), &s, 0.0, &b, &b, gravel(0.006), &initial, -1.0, 1.0).is_err());
    assert_eq!(initial, bed(11, 0.1));
}

#[test]
#[ignore = "explicit source-bound sediment component measurement"]
fn native_bank_to_bed_probe() {
    use serde_json::{json, Value};
    use std::{fs, path::Path, time::Instant};
    let request_path = std::env::var("VZ_RIVER_SEDIMENT_REQUEST").unwrap();
    let request: Value = serde_json::from_slice(&fs::read(request_path).unwrap()).unwrap();
    let flow_report: Value =
        serde_json::from_slice(&fs::read(request["flow_report"].as_str().unwrap()).unwrap())
            .unwrap();
    let status: Value = serde_json::from_slice(
        &fs::read(Path::new(request["world_root"].as_str().unwrap()).join("status.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(flow_report["source_state"], status["state_digest"]);
    let native_z = fs::read(
        Path::new(request["world_root"].as_str().unwrap()).join("terrain/fields/elevation_m.f32"),
    )
    .unwrap();
    let mut cases = Vec::new();
    for c in flow_report["results"].as_array().unwrap() {
        if c["status"] != "flow_and_erosion_proposal_only" {
            continue;
        }
        let hydro: Hydraulics = serde_json::from_value(c["hydraulics"].clone()).unwrap();
        let s: Vec<_> = (0..1001)
            .map(|i| Station {
                distance_m: 20.0 * i as f64,
                curvature_per_m: 0.0,
            })
            .collect();
        let banks = vec![
            BankMaterial {
                bank_height_m: hydro.depth_m + 1.0,
                ..bank(1e-9)
            };
            s.len()
        ];
        // Explicit grain and bed-cover scenarios; no inferred native gravel map.
        let mut outputs = Vec::new();
        for diameter in [0.006, 0.029] {
            let cell = c["cell"].as_u64().unwrap() as usize;
            let z_up =
                f32::from_le_bytes(native_z[4 * cell..4 * cell + 4].try_into().unwrap()) as f64;
            let initial: Vec<_> = (0..1000)
                .map(|i| BedCell {
                    bed_elevation_m: z_up
                        - c["native_coarse_slope"].as_f64().unwrap() * (20.0 * i as f64 + 10.0),
                    mobile_thickness_m: 0.0,
                })
                .collect();
            let start = Instant::now();
            let p = propose_step(
                hydro,
                &s,
                0.0,
                &banks,
                &banks,
                gravel(diameter),
                &initial,
                0.0,
                86400.0,
            )
            .unwrap();
            outputs.push(json!({"assumed_grain_m":diameter,"assumed_mobile_bed_m":0.0,
                "bank_solid_m3":p.bank.removed_solid_m3,"deposited_solid_m3":p.deposited_solid_m3.iter().sum::<f64>(),
                "downstream_solid_m3":p.downstream_solid_m3,"remaining_unrouted_bank_m3":p.bank.unrouted_solid_m3,"bed_water_space_displaced_m3":p.bed_displacement_bulk_m3,
                "residual_solid_m3":p.residual_solid_m3,"max_bed_change_m":p.next_bed.iter().zip(&initial).map(|(c,b)|(c.bed_elevation_m-b.bed_elevation_m).abs()).fold(0.0,f64::max),
                "wall_seconds":start.elapsed().as_secs_f64()}));
        }
        cases.push(json!({"cell":c["cell"],"scenarios":outputs}));
    }
    let output = Path::new(request["output"].as_str().unwrap());
    assert!(!output.exists());
    fs::write(output,serde_json::to_vec_pretty(&json!({"model":MODEL,"source_state":status["state_digest"],
        "cases":cases,"morphology_history_evolved":false,"terrain_published":false,
        "scope":"single-day conditional bed candidates; native coarse grade, static banks and assumed uniform grains; water displacement not committed"})).unwrap()).unwrap();
    println!("{}", output.display());
}

fn lateral_force(ratio: f64) -> LateralForcing {
    LateralForcing {
        total_shear_pa: 8.0,
        near_bed_velocity_ratio: ratio,
        slope_factor: 0.5,
    }
}
#[test]
fn lateral_drift_builds_inner_deposit_and_mirrors_without_losing_solid() {
    let initial = vec![
        BedCell {
            bed_elevation_m: 100.0,
            mobile_thickness_m: 1.0
        };
        16
    ];
    let left = propose_lateral_step(
        h(),
        gravel(0.002),
        &initial,
        100.0,
        lateral_force(-0.02),
        21600.0,
    )
    .unwrap();
    let right = propose_lateral_step(
        h(),
        gravel(0.002),
        &initial,
        100.0,
        lateral_force(0.02),
        21600.0,
    )
    .unwrap();
    assert!(left.next_bed[0].bed_elevation_m > 100.0);
    assert!(left.next_bed[15].bed_elevation_m < 100.0);
    for (a, b) in left.next_bed.iter().zip(right.next_bed.iter().rev()) {
        close(a.bed_elevation_m, b.bed_elevation_m, 1e-12);
    }
    close(left.residual_solid_m3, 0.0, 1e-8);
    assert_eq!(left.face_transfer_solid_m3[0], 0.0);
    assert_eq!(left.face_transfer_solid_m3[16], 0.0);
}
#[test]
fn lateral_equilibrium_is_force_balance_not_forced_bar_volume() {
    let g = gravel(0.002);
    let f = lateral_force(-0.02);
    let theta: f64 = 8.0 / ((2650.0 - 1000.0) * 9.81 * 0.002);
    let slope = f.near_bed_velocity_ratio * theta.sqrt() / f.slope_factor;
    let n = 16;
    let initial: Vec<_> = (0..n)
        .map(|j| BedCell {
            bed_elevation_m: 100.0 + slope * (j as f64 + 0.5) * 80.0 / n as f64,
            mobile_thickness_m: 1.0,
        })
        .collect();
    let p = propose_lateral_step(h(), g, &initial, 100.0, f, 86400.0 * 365.0).unwrap();
    for (a, b) in initial.iter().zip(p.next_bed) {
        close(a.bed_elevation_m, b.bed_elevation_m, 2e-12);
    }
    assert!(p.face_transfer_solid_m3.iter().all(|f| f.abs() < 1e-7));
}
#[test]
fn lateral_zero_drive_and_insufficient_supply_are_honest() {
    let initial = vec![
        BedCell {
            bed_elevation_m: 100.0,
            mobile_thickness_m: 0.0
        };
        16
    ];
    let p = propose_lateral_step(
        h(),
        gravel(0.002),
        &initial,
        100.0,
        lateral_force(0.0),
        86400.0,
    )
    .unwrap();
    assert_eq!(p.next_bed, initial);
    assert!(propose_lateral_step(
        h(),
        gravel(0.002),
        &initial,
        100.0,
        lateral_force(-0.02),
        86400.0
    )
    .is_err());
    let low = propose_lateral_step(
        h(),
        gravel(0.029),
        &initial,
        100.0,
        lateral_force(-0.02),
        86400.0,
    )
    .unwrap();
    assert_eq!(low.next_bed, initial);
}

fn lateral_mode_error(n: usize, steps: usize, continuous_time: bool) -> f64 {
    let flow = h();
    let g = gravel(0.002);
    let f = lateral_force(0.0);
    let total = 864000.0;
    let dt = total / steps as f64;
    let theta: f64 = 8.0 / ((2650.0 - 1000.0) * 9.81 * 0.002);
    let diffusivity = capacity_per_width_m2_s(flow, g, 8.0).unwrap() * f.slope_factor
        / theta.sqrt()
        / (1.0 - g.porosity);
    let wave = std::f64::consts::PI / flow.width_m;
    let factor = if continuous_time {
        (-diffusivity * wave * wave * total).exp()
    } else {
        (1.0 + diffusivity * wave * wave * dt).powf(-(steps as f64))
    };
    let mut state: Vec<_> = (0..n)
        .map(|j| BedCell {
            bed_elevation_m: 0.01 * (wave * (j as f64 + 0.5) * flow.width_m / n as f64).cos(),
            mobile_thickness_m: 1.0,
        })
        .collect();
    for _ in 0..steps {
        state = propose_lateral_step(flow, g, &state, 100.0, f, dt)
            .unwrap()
            .next_bed;
    }
    (state
        .iter()
        .enumerate()
        .map(|(j, c)| {
            let exact = 0.01 * factor * (wave * (j as f64 + 0.5) * flow.width_m / n as f64).cos();
            (c.bed_elevation_m - exact).powi(2)
        })
        .sum::<f64>()
        / n as f64)
        .sqrt()
}
#[test]
fn lateral_slope_relaxation_converges_to_independent_diffusion_mode() {
    let a = lateral_mode_error(16, 10, false);
    let b = lateral_mode_error(32, 10, false);
    let c = lateral_mode_error(64, 10, false);
    assert!(a / b > 3.9 && b / c > 3.9, "spatial {a} {b} {c}");
    let a = lateral_mode_error(512, 10, true);
    let b = lateral_mode_error(512, 20, true);
    let c = lateral_mode_error(512, 40, true);
    assert!(a / b > 1.9 && b / c > 1.9, "temporal {a} {b} {c}");
}

#[test]
fn native_bend_hydraulics_drive_lateral_material_transfer() {
    let mut s = stations(11);
    for point in &mut s {
        point.curvature_per_m = 0.0001;
    }
    let rows = vec![
        vec![
            BedCell {
                bed_elevation_m: 100.0,
                mobile_thickness_m: 1.0
            };
            16
        ];
        10
    ];
    let step =
        propose_lateral_reach(h(), &s, 0.0, 0.4, gravel(0.002), &rows, 0.88, 21600.0).unwrap();
    assert!(step[9].next_bed[0].bed_elevation_m > 100.0);
    assert!(step[9].next_bed[15].bed_elevation_m < 100.0);
    for row in &step {
        close(row.residual_solid_m3, 0.0, 1e-8);
    }
    let straight = propose_lateral_reach(
        h(),
        &stations(11),
        0.0,
        0.4,
        gravel(0.002),
        &rows,
        0.88,
        21600.0,
    )
    .unwrap();
    for (before, after) in rows.iter().zip(straight) {
        assert_eq!(*before, after.next_bed);
    }
}

#[test]
fn bank_failure_conserves_material_and_lowers_potential_to_analytic_slope() {
    let b = vec![
        BedCell {
            bed_elevation_m: 2.0,
            mobile_thickness_m: 2.0,
        },
        BedCell {
            bed_elevation_m: 0.0,
            mobile_thickness_m: 0.0,
        },
    ];
    let p = relax_mobile_slope(&b, 1.0, 1.0, gravel(0.002), 9.81, 1.0, 1e-10, 8).unwrap();
    close(p.next_bed[0].bed_elevation_m, 1.5, 1e-12);
    close(p.next_bed[1].bed_elevation_m, 0.5, 1e-12);
    close(p.transferred_solid_m3, 0.5 * 0.65, 1e-12);
    close(p.solid_residual_m3, 0.0, 1e-12);
    assert!(p.potential_energy_change_j < 0.0);
    assert_eq!(b[0].bed_elevation_m, 2.0);
}
#[test]
fn rock_substrate_is_not_invented_as_mobile_bank_supply() {
    let b = vec![
        BedCell {
            bed_elevation_m: 4.0,
            mobile_thickness_m: 0.1,
        },
        BedCell {
            bed_elevation_m: 0.0,
            mobile_thickness_m: 0.0,
        },
    ];
    let p = relax_mobile_slope(&b, 1.0, 1.0, gravel(0.002), 9.81, 0.7, 1e-10, 8).unwrap();
    close(p.next_bed[0].bed_elevation_m, 3.9, 1e-12);
    close(p.next_bed[1].bed_elevation_m, 0.1, 1e-12);
    assert_eq!(p.next_bed[0].mobile_thickness_m, 0.0);
    close(p.solid_residual_m3, 0.0, 1e-12);
}
#[test]
fn bank_relaxation_budget_fails_without_mutating_the_input() {
    let b = vec![
        BedCell {
            bed_elevation_m: 10.0,
            mobile_thickness_m: 10.0,
        },
        BedCell {
            bed_elevation_m: 0.0,
            mobile_thickness_m: 0.0,
        },
        BedCell {
            bed_elevation_m: 0.0,
            mobile_thickness_m: 0.0,
        },
    ];
    assert!(relax_mobile_slope(&b, 1.0, 1.0, gravel(0.002), 9.81, 0.7, 1e-12, 1).is_err());
    let p = relax_mobile_slope(&b, 1.0, 1.0, gravel(0.002), 9.81, 0.7, 1e-10, 128).unwrap();
    assert!(p.max_mobile_excess_m <= 1e-10);
    assert_eq!(b[0].bed_elevation_m, 10.0);
}
