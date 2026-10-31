//! SEED-UPLANDS (lane UP, 2026-10-13): synthetic belt-decay experiment for
//! `critical-slope-envelope-v1` (uplands.rs) with the subgrid stream-power
//! laws (relief.rs), no tectonics. PREREG.txt P4 in
//! `generated-test-data/`.
//!
//! Setup (declared): L1 grid (288 x 576, Earth radius); one continent
//! (lat -30..30, lon 0..90 deg) of 35 km crust at +300 m; a 4-cell-wide
//! north-south belt of crust thick enough for +6 km in its middle; ocean
//! elsewhere held at -4 km. Elevation is Airy: z = 300 + 0.16 (t - 35 km)
//! (0.16 m per m of crust, measured from extras_v2 seed 42 L1 output), so
//! every eroded metre returns 0.84 m of rock uplift. Erosion runs through
//! the shared stage (`erode_step_bounded`, toy K = 2e-6/yr, m = 0.5,
//! D = 0.1 m2/yr) with the subgrid effective areas, exactly as the preview
//! tier does. Simplification: eroded material leaves (no sediment kept).
//!
//! The quick test checks the claims on 40 My; the ignored driver writes
//! the full 400 My curves (UP_OUT=<dir>).

use pf_grid::GridSpec;
use pf_surface::erosion::{erode_step_bounded, DiffusionStencil, ErosionParams};
use pf_surface::relief::{subgrid_step_with_floor, NonlinearSubgridParams, SubgridInputs};
use pf_surface::routing::route;
use pf_surface::uplands::{adjacency, threshold_lowering_isostatic, CriticalSlopeParams};
use pf_surface::water::largest_ocean;

const NLAT: usize = 288;
const NLON: usize = 576;
const AIRY: f64 = 0.16;

struct Variant {
    name: &'static str,
    n: f64,
    envelope: bool,
}

struct Curve {
    t_my: Vec<f64>,
    belt_mean_m: Vec<f64>,
    belt_max_m: Vec<f64>,
    flank_mean_m: Vec<f64>,
    interior_mean_m: Vec<f64>,
    belt_crust_km: Vec<f64>,
}

fn run(v: &Variant, duration_my: f64, profile_at: &[f64]) -> (Curve, Vec<(f64, Vec<f64>)>) {
    let grid = GridSpec::new(NLAT, NLON, 6_371_000.0).unwrap();
    let n = grid.cells();
    let areas: Vec<f64> = (0..n)
        .map(|c| grid.cell_area_m2(c / NLON).unwrap())
        .collect();
    let cell_size: Vec<f64> = areas.iter().map(|a| a.sqrt()).collect();
    let stencil = DiffusionStencil::build(&grid).unwrap();
    let adj = adjacency(&grid).unwrap();
    let ep = ErosionParams {
        k_per_my: 2.0,
        m: 0.5,
        n: 1.0,
        diffusivity_m2_per_my: 0.1e6,
        stability_limit: 0.5,
    };
    let rp = NonlinearSubgridParams {
        k_per_my: 2.0,
        slope_exponent_n: v.n,
        reference_steepness_m: 500.0,
        channel_head_length_m: 100.0,
    };
    let cs = CriticalSlopeParams {
        critical_slope: 0.04,
        knickpoint_celerity_m_per_my: 70_000.0,
    };
    // Continent rows 96..192 (lat 30..-30), cols 0..144 (lon 0..90).
    let cont = |c: usize| {
        let (i, j) = (c / NLON, c % NLON);
        (96..192).contains(&i) && j < 144
    };
    let belt = |c: usize| {
        let (i, j) = (c / NLON, c % NLON);
        (100..188).contains(&i) && (70..74).contains(&j)
    };
    let flank = |c: usize| {
        let (i, j) = (c / NLON, c % NLON);
        (100..188).contains(&i) && ((66..70).contains(&j) || (74..78).contains(&j))
    };
    let interior = |c: usize| {
        let (i, j) = (c / NLON, c % NLON);
        (110..178).contains(&i) && ((20..40).contains(&j) || (104..124).contains(&j))
    };
    let mut crust: Vec<f64> = (0..n)
        .map(|c| {
            if belt(c) {
                35_000.0 + 5_700.0 / AIRY
            } else {
                35_000.0
            }
        })
        .collect();
    let elev = |crust: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|c| if cont(c) { 300.0 + AIRY * (crust[c] - 35_000.0) } else { -4000.0 })
            .collect()
    };
    let mut relief = vec![0.0f64; n];
    let rain = vec![1.0e6; n];
    let sed0 = vec![0.0f64; n];
    let dt = 2.0;
    let steps = (duration_my / dt).round() as usize;
    let mut cv = Curve {
        t_my: vec![],
        belt_mean_m: vec![],
        belt_max_m: vec![],
        flank_mean_m: vec![],
        interior_mean_m: vec![],
        belt_crust_km: vec![],
    };
    let mut profiles = Vec::new();
    let row = 144usize;
    let record = |cv: &mut Curve, t: f64, z: &[f64], crust: &[f64]| {
        let mean = |f: &dyn Fn(usize) -> bool, x: &[f64]| {
            let (mut s, mut w) = (0.0, 0.0);
            for c in 0..n {
                if f(c) {
                    s += x[c] * areas[c];
                    w += areas[c];
                }
            }
            s / w
        };
        cv.t_my.push(t);
        cv.belt_mean_m.push(mean(&belt, z));
        cv.belt_max_m.push((0..n).filter(|&c| belt(c)).map(|c| z[c]).fold(f64::MIN, f64::max));
        cv.flank_mean_m.push(mean(&flank, z));
        cv.interior_mean_m.push(mean(&interior, z));
        cv.belt_crust_km.push(mean(&belt, crust) / 1000.0);
    };
    let z0 = elev(&crust);
    record(&mut cv, 0.0, &z0, &crust);
    for s in 1..=steps {
        let z = elev(&crust);
        let ocean = largest_ocean(&grid, &z, 0.0).unwrap();
        let routing = route(&grid, &z, &rain, &ocean, 0.0).unwrap();
        let wet: Vec<bool> = (0..n)
            .map(|i| ocean[i] || routing.fill_level_m[i] > z[i] + 1.0e-9)
            .collect();
        let flow_to: Vec<Option<usize>> = routing.flow_to.iter().map(|t| t.map(|e| e.0 as usize)).collect();
        let mut reach = vec![1.0f64; n];
        for i in 0..n {
            if let Some(j) = flow_to[i] {
                reach[i] = grid
                    .center_distance_m(i / NLON, i % NLON, j / NLON, j % NLON)
                    .unwrap()
                    .max(1.0);
            }
        }
        let floor = if v.envelope {
            Some(
                threshold_lowering_isostatic(&adj, &z, &wet, &routing.fill_level_m, &cell_size, &cs, AIRY, dt)
                    .unwrap()
                    .lowering_m,
            )
        } else {
            None
        };
        let sg = subgrid_step_with_floor(
            &SubgridInputs {
                receiver_first: &routing.receiver_first,
                flow_to: &flow_to,
                contributing_area_m2: &routing.contributing_area_m2,
                reach_length_m: &reach,
                wet: &wet,
                fill_level_m: &routing.fill_level_m,
                cell_size_m: &cell_size,
            },
            &z,
            &relief,
            &rp,
            ep.k_per_my,
            ep.m,
            dt,
            floor.as_deref(),
        )
        .unwrap();
        let mut eff = routing.clone();
        for i in 0..n {
            if !wet[i] && flow_to[i].is_some() {
                eff.contributing_area_m2[i] = sg.effective_area_m2[i];
            }
        }
        let st = erode_step_bounded(&grid, &stencil, &eff, &z, &sed0, &ocean, 0.0, dt, &ep, Some(&crust))
            .unwrap();
        for c in 0..n {
            if cont(c) {
                crust[c] -= st.bedrock_removed_m[c];
            }
            relief[c] = if wet[c] { 0.0 } else { ((z[c] + st.dz_m[c]) - sg.trunk_m[c]).max(0.0) };
        }
        let za = elev(&crust);
        let t = s as f64 * dt;
        record(&mut cv, t, &za, &crust);
        if profile_at.iter().any(|&p| (p - t).abs() < 1e-9) {
            profiles.push((t, (0..144).map(|j| za[row * NLON + j]).collect()));
        }
    }
    (cv, profiles)
}

fn variants() -> Vec<Variant> {
    vec![
        Variant { name: "n1", n: 1.0, envelope: false },
        Variant { name: "n1_env", n: 1.0, envelope: true },
        Variant { name: "n2", n: 2.0, envelope: false },
        Variant { name: "n2_env", n: 2.0, envelope: true },
    ]
}

/// Quick check (40 My): the envelope lowers the young belt edge faster and
/// leaves a taper; n = 2 keeps the old belt higher than n = 1.
#[test]
fn envelope_tapers_and_n2_slows_decay() {
    let res: Vec<Curve> = variants().iter().map(|v| run(v, 40.0, &[]).0).collect();
    let last = |c: &Curve| *c.belt_mean_m.last().unwrap();
    let (n1, n1e, n2, n2e) = (last(&res[0]), last(&res[1]), last(&res[2]), last(&res[3]));
    eprintln!("belt mean at 40 My: n1 {n1:.0} n1_env {n1e:.0} n2 {n2:.0} n2_env {n2e:.0}");
    assert!(n2 > n1, "n = 2 should keep the belt higher than n = 1");
    assert!(n2e <= n2 + 1e-6, "the envelope never raises the belt");
    assert!(n1e <= n1 + 1e-6);
    // Rebound: belt crust removed / surface lowered ~ 1 / AIRY.
    let c = &res[3];
    let dz = c.belt_mean_m[0] - last(c);
    let dc = (c.belt_crust_km[0] - c.belt_crust_km.last().unwrap()) * 1000.0;
    assert!(dc > 0.0 && (dz / dc - AIRY).abs() < 1e-6, "dz {dz} dc {dc}");
}

#[test]
#[ignore]
fn uplands_decay_curves() {
    let out = std::path::PathBuf::from(std::env::var("UP_OUT").unwrap());
    std::fs::create_dir_all(&out).unwrap();
    let dur: f64 = std::env::var("UP_DURATION_MY").ok().map(|s| s.parse().unwrap()).unwrap_or(400.0);
    for v in variants() {
        let (c, prof) = run(&v, dur, &[0.0, 2.0, 20.0, 50.0, 100.0, 200.0, 300.0, 400.0]);
        let mut s = String::from("t_my,belt_mean_m,belt_max_m,flank_mean_m,interior_mean_m,belt_crust_km\n");
        for k in 0..c.t_my.len() {
            s += &format!(
                "{},{:.1},{:.1},{:.1},{:.1},{:.3}\n",
                c.t_my[k], c.belt_mean_m[k], c.belt_max_m[k], c.flank_mean_m[k], c.interior_mean_m[k], c.belt_crust_km[k]
            );
        }
        std::fs::write(out.join(format!("decay_{}.csv", v.name)), s).unwrap();
        let mut p = String::from("t_my,col0..col143 (row 144, lon 0..90 deg)\n");
        for (t, z) in prof {
            p += &format!("{t},{}\n", z.iter().map(|x| format!("{x:.1}")).collect::<Vec<_>>().join(","));
        }
        std::fs::write(out.join(format!("profile_{}.csv", v.name)), p).unwrap();
        eprintln!("{} done", v.name);
    }
}
