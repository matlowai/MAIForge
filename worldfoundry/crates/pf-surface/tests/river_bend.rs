use pf_surface::river_bend::*;
use std::f64::consts::PI;

fn hydraulics() -> Hydraulics {
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
fn stations(n: usize, length: f64, f: impl Fn(f64) -> f64) -> Vec<Station> {
    (0..n)
        .map(|i| {
            let s = length * i as f64 / (n - 1) as f64;
            Station {
                distance_m: s,
                curvature_per_m: f(s),
            }
        })
        .collect()
}
fn bank() -> BankMaterial {
    BankMaterial {
        erodibility_m_pa_s: 1e-9,
        critical_shear_pa: 4.0,
        bank_height_m: 4.0,
        porosity: 0.35,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}, tol={tol}");
}

#[test]
fn constant_bend_matches_analytic_relaxation_and_correct_bank() {
    let h = hydraulics();
    let k = 1e-4;
    let s = stations(401, 4000.0, |_| k);
    let f = solve(h, &s, 0.0).unwrap();
    let steady = 0.5 * h.width_m * h.speed_m_s() * h.alpha() * k;
    for (i, p) in s.iter().enumerate() {
        close(
            f.velocity_perturbation_m_s[i],
            steady * (1.0 - (-p.distance_m / h.response_length_m()).exp()),
            1e-14,
        );
        assert!(f.right_shear_pa[i] >= f.left_shear_pa[i]);
    }
    let straight = stations(401, 4000.0, |_| 0.0);
    let g = solve(h, &straight, 0.01).unwrap();
    for (i, p) in straight.iter().enumerate() {
        close(
            g.velocity_perturbation_m_s[i],
            0.01 * (-p.distance_m / h.response_length_m()).exp(),
            1e-14,
        );
    }
}

fn harmonic_error(n: usize) -> f64 {
    let h = hydraulics();
    let k = 0.0001;
    let w = 2.0 * PI / 2000.0;
    let l = 1.0 / h.response_length_m();
    let s = stations(n, 10000.0, |s| k * (w * s).cos());
    let f = solve(h, &s, 0.0).unwrap();
    let scale = 0.5 * h.width_m * h.speed_m_s() * k;
    let re = scale * (h.alpha() * l * l - w * w) / (l * l + w * w);
    let im = -scale * (h.alpha() + 1.0) * l * w / (l * l + w * w);
    (s.iter()
        .enumerate()
        .map(|(i, p)| {
            let exact = re * (w * p.distance_m).cos()
                - im * (w * p.distance_m).sin()
                - re * (-l * p.distance_m).exp();
            (f.velocity_perturbation_m_s[i] - exact).powi(2)
        })
        .sum::<f64>()
        / n as f64)
        .sqrt()
}
#[test]
fn harmonic_reference_converges_quadratically_with_station_spacing() {
    let coarse = harmonic_error(101);
    let fine = harmonic_error(201);
    let finest = harmonic_error(401);
    assert!(
        coarse / fine > 3.8 && fine / finest > 3.8,
        "{coarse} {fine} {finest}"
    );
    assert!(finest < 1e-5);
}
#[test]
fn upstream_geometry_matters_and_mirror_reverses_response() {
    let h = hydraulics();
    let a = stations(201, 2000.0, |s| if s < 1000.0 { 1e-4 } else { 0.0 });
    let mirror: Vec<_> = a
        .iter()
        .map(|s| Station {
            distance_m: s.distance_m,
            curvature_per_m: -s.curvature_per_m,
        })
        .collect();
    let f = solve(h, &a, 0.0).unwrap();
    let m = solve(h, &mirror, 0.0).unwrap();
    assert!(f.velocity_perturbation_m_s[110].abs() > 1e-6); // local curvature is zero
    for i in 0..a.len() {
        close(
            f.velocity_perturbation_m_s[i],
            -m.velocity_perturbation_m_s[i],
            1e-14,
        );
    }
}
#[test]
fn bank_prisms_have_independent_solid_budget_and_no_fake_deposit() {
    let h = hydraulics();
    let s = stations(101, 1000.0, |_| 0.0);
    let a = vec![bank(); s.len()];
    let mut b = a.clone();
    for x in &mut b {
        x.porosity = 0.5;
        x.bank_height_m = 2.0;
    }
    let dt = 86400.0;
    let p = propose_bank_erosion(h, &s, 0.0, &a, &b, dt).unwrap();
    // Straight reach tau=rho Cf U²=8 Pa, erosion rate=4e-9 m/s.
    let retreat = 4e-9 * dt;
    close(p.removed_bulk_m3, 1000.0 * retreat * 6.0, 1e-12);
    close(
        p.removed_solid_m3,
        1000.0 * retreat * (4.0 * 0.65 + 2.0 * 0.5),
        1e-12,
    );
    assert_eq!(p.removed_solid_m3, p.unrouted_solid_m3);
    let half = propose_bank_erosion(h, &s, 0.0, &a, &b, dt / 2.0).unwrap();
    assert_eq!(p.removed_solid_m3, 2.0 * half.removed_solid_m3);
}
#[test]
fn dry_threshold_and_resistant_banks_do_not_erode() {
    let mut h = hydraulics();
    let s = stations(30, 1000.0, |_| 1e-4);
    let mut b = vec![bank(); s.len()];
    h.discharge_m3_s = 0.0;
    let p = propose_bank_erosion(h, &s, 0.0, &b, &b, 86400.0).unwrap();
    assert_eq!(p.removed_solid_m3, 0.0);
    h = hydraulics();
    for x in &mut b {
        x.critical_shear_pa = 1e5;
    }
    assert_eq!(
        propose_bank_erosion(h, &s, 0.0, &b, &b, 86400.0)
            .unwrap()
            .removed_solid_m3,
        0.0
    );
    for x in &mut b {
        x.critical_shear_pa = 0.0;
        x.erodibility_m_pa_s = 0.0;
    }
    assert_eq!(
        propose_bank_erosion(h, &s, 0.0, &b, &b, 86400.0)
            .unwrap()
            .removed_solid_m3,
        0.0
    );
}
#[test]
fn refuses_invalid_and_outside_domain_instead_of_clipping() {
    let h = hydraulics();
    let s = stations(20, 1000.0, |_| 0.0);
    let b = vec![bank(); s.len()];
    let mut bad = s.clone();
    bad[2].distance_m = bad[1].distance_m;
    assert!(solve(h, &bad, 0.0).is_err());
    bad = s.clone();
    bad[2].curvature_per_m = f64::NAN;
    assert!(solve(h, &bad, 0.0).is_err());
    assert!(solve(h, &stations(20, 1000.0, |_| 0.01), 0.0).is_err());
    assert!(solve(h, &s, h.speed_m_s()).is_err());
    let mut overflow = h;
    overflow.width_m = f64::MAX;
    assert!(solve(overflow, &s, 0.0).is_err());
    let mut fast = h;
    fast.discharge_m3_s = 1e6;
    assert!(solve(fast, &s, 0.0).is_err());
    let mut wrong = b.clone();
    wrong[0].porosity = 1.0;
    assert!(propose_bank_erosion(h, &s, 0.0, &wrong, &b, 1.0).is_err());
    assert!(propose_bank_erosion(h, &s, 0.0, &b, &b, 2000.0 * 31557600.0).is_err());
    assert!(propose_bank_erosion(h, &s, 0.0, &b, &b, f64::NAN).is_err());
}
#[test]
fn forcing_and_material_control_actual_erosion_rates() {
    let h = hydraulics();
    let s = stations(200, 10000.0, |_| 1e-4);
    let b = vec![bank(); s.len()];
    let p = propose_bank_erosion(h, &s, 0.0, &b, &b, 86400.0).unwrap();
    let mut more = h;
    more.discharge_m3_s *= 1.2;
    let q = propose_bank_erosion(more, &s, 0.0, &b, &b, 86400.0).unwrap();
    assert!(q.removed_solid_m3 > p.removed_solid_m3);
    let mut resistant = b.clone();
    for a in &mut resistant {
        a.erodibility_m_pa_s = 0.0;
    }
    let r = propose_bank_erosion(h, &s, 0.0, &resistant, &b, 86400.0).unwrap();
    assert!(r.left_retreat_m.iter().all(|x| *x == 0.0));
    assert_eq!(r.right_retreat_m, p.right_retreat_m);
}

/// Engineering capture of this library on native world fields, not coursework
/// or an accepted terrain publication. Input request is source-pinned by the
/// lane manifest; this test independently loads Q/elevation/receiver/reach.
#[test]
#[ignore = "explicit bounded actual-world component probe; VZ_RIVER_PROBE_REQUEST required"]
fn native_source_probe() {
    use serde_json::{json, Value};
    use std::{fs, path::Path, time::Instant};
    let request = std::env::var("VZ_RIVER_PROBE_REQUEST").unwrap();
    let input: Value = serde_json::from_slice(&fs::read(&request).unwrap()).unwrap();
    let root = Path::new(input["world_root"].as_str().unwrap());
    let read = |path: &str| fs::read(root.join(path)).unwrap();
    let qb = read("hydrology/fields/discharge_proxy_m3_s.f32");
    let zb = read("terrain/fields/elevation_m.f32");
    let tb = read("hydrology/fields/flow_to_element.u32");
    let rb = read("hydrology/fields/reach_id.u32");
    let wb = read("hydrology/fields/water_body_id.u32");
    let f32_at =
        |b: &[u8], i: usize| f32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap()) as f64;
    let u32_at = |b: &[u8], i: usize| u32::from_le_bytes(b[4 * i..4 * i + 4].try_into().unwrap());
    let status: Value = serde_json::from_slice(&read("status.json")).unwrap();
    assert_eq!(status["state_digest"], input["state_digest"]);
    let mut results = Vec::new();
    for c in input["cases"].as_array().unwrap() {
        let cell = c["cell"].as_u64().unwrap() as usize;
        let q = f32_at(&qb, cell);
        let z = f32_at(&zb, cell);
        let receiver = u32_at(&tb, cell) as usize;
        assert_eq!(q, c["native_q_m3_s"].as_f64().unwrap());
        assert_eq!(receiver as u64, c["receiver_cell"].as_u64().unwrap());
        let reach = u32_at(&rb, cell);
        if receiver >= qb.len() / 4
            || u32_at(&wb, cell) != u32::MAX
            || u32_at(&wb, receiver) != u32::MAX
        {
            results.push(json!({"cell":cell,"status":"ineligible", "reason":"lake/marine/terminal boundary requires water-surface boundary solve; seafloor drop is not channel slope"}));
            continue;
        }
        let slope = (z - f32_at(&zb, receiver)) / c["edge_length_m"].as_f64().unwrap();
        let width = c["assumed_width_m"].as_f64().unwrap();
        let cf = 0.008;
        let g = 9.81;
        let depth = (cf * (q / width).powi(2) / (g * slope)).cbrt();
        let h = Hydraulics {
            discharge_m3_s: q,
            width_m: width,
            depth_m: depth,
            friction_cf: cf,
            scour_factor: 3.0,
            gravity_m_s2: g,
            water_density_kg_m3: 1000.0,
        };
        // Actual coarse native edge has no resolved fine bend: zero curvature.
        // Do not present an invented sine curve as native river geometry.
        let s = stations(1001, 20000.0, |_| 0.0);
        let a = vec![
            BankMaterial {
                bank_height_m: depth + 1.0,
                ..bank()
            };
            s.len()
        ];
        let start = Instant::now();
        let outcome = propose_bank_erosion(h, &s, 0.0, &a, &a, 86400.0);
        match outcome {
            Ok(proposal) => {
                let mut checksum = 0.0;
                let rounds = 2000;
                for _ in 0..rounds {
                    checksum += std::hint::black_box(
                        propose_bank_erosion(h, &s, 0.0, &a, &a, 86400.0).unwrap(),
                    )
                    .removed_solid_m3;
                }
                results.push(json!({"cell":cell,"receiver":receiver,"reach":reach,"hydraulics":h,"native_coarse_slope":slope,
                    "status":"flow_and_erosion_proposal_only","proposal":proposal,"benchmark_wall_seconds":start.elapsed().as_secs_f64(),
                    "repeated_fixed_geometry_evaluations":rounds,"stations":s.len(),"checksum":checksum,
                    "morphological_years_evolved":0,"geometry":"straight subreach of coarse native receiver edge; fine planform unresolved"}));
            }
            Err(reason) => results
                .push(json!({"cell":cell,"hydraulics":h,"status":"ineligible","reason":reason})),
        }
    }
    let output = json!({"model":MODEL,"source_state":status["state_digest"],"results":results,
        "native_component_only":true,"terrain_integration":false,"history_years_evolved":0});
    let out = Path::new(input["output"].as_str().unwrap());
    assert!(!out.exists(), "preserve prior probe");
    fs::write(out, serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    println!("native river component report: {}", out.display());
}

#[test]
fn secondary_flow_relaxes_to_curvature_forcing_and_retains_upstream_memory() {
    let h = hydraulics();
    let s = stations(201, 1000.0, |_| 0.0001);
    let p = solve_secondary(h, &s, 0.0, 0.4).unwrap();
    for (i, point) in s.iter().enumerate() {
        close(
            p.intensity_m_s[i],
            h.depth_m
                * h.speed_m_s()
                * 0.0001
                * (1.0 - (-point.distance_m / p.adaptation_length_m).exp()),
            1e-15,
        );
    }
    assert!(p.near_bed_velocity_ratio[200] < 0.0); // inward, toward left bank
    let straight = stations(201, 1000.0, |_| 0.0);
    let tail = solve_secondary(h, &straight, p.intensity_m_s[200], 0.4).unwrap();
    assert!(tail.intensity_m_s[1] > 0.0);
    assert!(tail.intensity_m_s[200] < 1e-12);
    let mut dry = h;
    dry.discharge_m3_s = 0.0;
    assert!(solve_secondary(dry, &s, 0.0, 0.4)
        .unwrap()
        .near_bed_velocity_ratio
        .iter()
        .all(|v| *v == 0.0));
}
