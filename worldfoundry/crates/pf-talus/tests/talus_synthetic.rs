//! Lane TALUS-1: synthetic cliff tests (determinism, mass balance, fall sorting, repose, clean faces, tile seams).

use pf_talus::model::{run, Rect, RunOut};
use pf_talus::output::{block_record, field_bytes};
use pf_talus::params::TalusParams;
use pf_talus::surface::Surface;
use pf_talus::{domain_for, full_rect};

/// 160 x 140 posts at 2 m: plateau (north) at 150 m, a 70 deg cliff down to 60 m, a 12 deg ramp, then flat ground
/// at ~20 m; basalt, regime unknown; the plateau edge wiggles so chutes and cones have something to follow.
fn synth() -> Surface {
    let (nx, ny) = (160usize, 140usize);
    let dx = 2.0;
    let mut z = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let x = i as f64 * dx;
            let y = j as f64 * dx;
            let edge = 200.0 + 8.0 * (x / 37.0).sin();
            let cliff_w = (150.0 - 60.0) / 70f64.to_radians().tan();
            let h = if y >= edge { 150.0 } else if y >= edge - cliff_w { 60.0 + (y - (edge - cliff_w)) * 70f64.to_radians().tan() } else {
                let d = (edge - cliff_w) - y;
                (60.0 - d * 12f64.to_radians().tan()).max(20.0)
            };
            z[j * nx + i] = h;
        }
    }
    let n = nx * ny;
    Surface { dx, gx0: 1000, gy0: 2000, nx, ny, z, litho: vec![3; n], facies: vec![0; n], regime: vec![255; n], sink: vec![0; n], unit: vec![3 << 8; n], digests: vec![] }
}

fn params() -> TalusParams {
    let mut p = TalusParams::default();
    p.l_max_m = 150.0;
    p.chute_len_max_m = 60.0;
    p.chute_spacing_m = 40.0;
    p.relax_levels = vec![(8, 6), (4, 8)];
    p.relax_iters = 24;
    p.threads = 4;
    p
}

fn bytes(o: &RunOut) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    (field_bytes(o).into_iter().map(|f| f.2).collect(), o.blocks.iter().map(block_record).collect())
}

#[test]
fn deterministic_and_thread_independent() {
    let s = synth();
    let mut p = params();
    let a = run(&s, &p, full_rect(&s)).unwrap();
    p.threads = 1;
    let b = run(&s, &p, full_rect(&s)).unwrap();
    assert_eq!(bytes(&a), bytes(&b));
    assert!(!a.blocks.is_empty());
}

#[test]
fn field_volume_balances_exactly() {
    let s = synth();
    let o = run(&s, &params(), full_rect(&s)).unwrap();
    for b in o.budget.values() {
        assert_eq!(b.field_bulk_u, b.field_dep_u + b.field_water_u + b.field_off_u);
        assert!(b.rock_u > 0 && b.weathered_u > 0);
    }
}

fn mean_src(o: &RunOut, c: usize) -> f64 {
    let h = &o.classes[c].h_src;
    let tot: i64 = h.iter().sum();
    h.iter().enumerate().map(|(i, &v)| (i as f64 + 0.5) * pf_talus::model::HBIN_M * v as f64).sum::<f64>() / tot.max(1) as f64
}

fn q_src(o: &RunOut, c: usize, q: f64) -> f64 {
    let h = &o.classes[c].h_src;
    let tot: i64 = h.iter().sum();
    let mut cum = 0i64;
    for (i, &v) in h.iter().enumerate() {
        cum += v;
        if cum as f64 >= q * tot as f64 { return (i as f64 + 0.5) * pf_talus::model::HBIN_M; }
    }
    f64::NAN
}

/// Fall sorting: the distal tail (p90 of runout from the source, volume weighted) grows with size from fines through
/// the mid field class to the largest block class, and the largest blocks also have the longest mean runout.
#[test]
fn fall_sorting_big_blocks_run_farthest() {
    let s = synth();
    let p = params();
    let o = run(&s, &p, full_rect(&s)).unwrap();
    let big = (0..p.n_classes()).filter(|&c| p.is_block_class(c) && o.classes[c].rest_u > 0).max().unwrap();
    let (f90, m90, b90) = (q_src(&o, 0, 0.9), q_src(&o, 4, 0.9), q_src(&o, big, 0.9));
    assert!(f90 < m90 && m90 < b90, "p90 fine {f90} mid {m90} big {b90}");
    let (fm, mm, bm) = (mean_src(&o, 0), mean_src(&o, 4), mean_src(&o, big));
    assert!(fm < mm && mm < bm, "mean fine {fm} mid {mm} big {bm}");
}

#[test]
fn debris_rests_near_repose_and_faces_stay_clean() {
    let s = synth();
    let p = params();
    let o = run(&s, &p, full_rect(&s)).unwrap();
    let (nx, ny) = (o.nx, o.ny);
    let h: Vec<f64> = (0..nx * ny).map(|k| s.z[k] - o.chute_mm[k] as f64 / 1000.0 + o.thickness_mm[k] as f64 / 1000.0).collect();
    let mut deep = 0;
    let mut steep = 0;
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let k = j * nx + i;
            if o.thickness_mm[k] < 1000 { continue; }
            deep += 1;
            let gx = (h[k + 1] - h[k - 1]) / 4.0;
            let gy = (h[k + nx] - h[k - nx]) / 4.0;
            if (gx * gx + gy * gy).sqrt() > 40f64.to_radians().tan() { steep += 1; }
        }
    }
    assert!(deep > 100);
    assert!((steep as f64) < 0.1 * deep as f64, "{steep} of {deep} debris posts steeper than 40 deg");
    // visible blocks do not sit on the cliff
    let on_face = o.blocks.iter().filter(|b| b.visible).filter(|b| {
        let i = ((b.x / 2.0).round() as i64 - s.gx0).clamp(1, nx as i64 - 2) as usize;
        let j = ((b.y / 2.0).round() as i64 - s.gy0).clamp(1, ny as i64 - 2) as usize;
        let k = j * nx + i;
        let gx = (h[k + 1] - h[k - 1]) / 4.0;
        let gy = (h[k + nx] - h[k - nx]) / 4.0;
        (gx * gx + gy * gy).sqrt() > 45f64.to_radians().tan()
    }).count();
    let vis = o.blocks.iter().filter(|b| b.visible).count();
    assert!(vis > 0 && (on_face as f64) < 0.02 * vis as f64, "{on_face} of {vis} visible blocks on > 45 deg");
}

#[test]
fn tiles_with_halo_match_the_whole_domain() {
    let s = synth();
    let p = params();
    let full = run(&s, &p, full_rect(&s)).unwrap();
    let (fb, _) = bytes(&full);
    let mid = s.gx0 + 80;
    let halo = p.halo_m(s.dx);
    for core in [Rect { gx0: s.gx0, gy0: s.gy0, gx1: mid, gy1: s.gy0 + s.ny as i64 - 1 }, Rect { gx0: mid, gy0: s.gy0, gx1: s.gx0 + s.nx as i64 - 1, gy1: s.gy0 + s.ny as i64 - 1 }] {
        let t = run(&domain_for(&s, core, halo), &p, core).unwrap();
        let (tb, tblocks) = bytes(&t);
        for (f, (a, b)) in fb.iter().zip(tb.iter()).enumerate() {
            let bpp = a.len() / (full.nx * full.ny);
            for r in 0..t.ny {
                let kf = (r * full.nx + (core.gx0 - s.gx0) as usize) * bpp;
                let kt = r * t.nx * bpp;
                assert_eq!(&a[kf..kf + t.nx * bpp], &b[kt..kt + t.nx * bpp], "field {f} row {r}");
            }
        }
        let want: Vec<Vec<u8>> = full.blocks.iter().filter(|b| b.x >= core.gx0 as f64 * 2.0 && b.x < core.gx1 as f64 * 2.0 && b.y >= core.gy0 as f64 * 2.0 && b.y < core.gy1 as f64 * 2.0).map(block_record).collect();
        assert_eq!(want, tblocks);
    }
}

#[test]
fn throughput_checks_pass_on_a_real_run() {
    let s = synth();
    let p = params();
    let o = run(&s, &p, full_rect(&s)).unwrap();
    for b in o.budget.values() {
        // rounding allowance: one unit per parcel / block rounding
        let rounding = (b.sources as i64) * (p.n_classes() as i64 * p.sub_parcels as i64 + 4);
        pf_talus::budget::check_unit(b, p.porosity, 1e-9, rounding).unwrap();
    }
    pf_talus::budget::check_relax(&o.stats.relax_throughput, 1e-9).unwrap();
    assert!(o.stats.relax_throughput.moved_m3 > 0.0);
}

/// Integrator request (2026-10-21): a 10 % transfer mismatch must fail the throughput check, also when it is small
/// next to the stock.
#[test]
fn ten_percent_transfer_mismatch_fails_the_throughput_check() {
    let s = synth();
    let p = params();
    let o = run(&s, &p, full_rect(&s)).unwrap();
    let b0 = o.budget.values().next().unwrap().clone();
    // field: 10 % of what moved vanishes between release and deposit
    let mut b = b0.clone();
    b.field_dep_u -= b.field_bulk_u / 10;
    let e = pf_talus::budget::check_unit(&b, p.porosity, 1e-3, 1000).unwrap_err();
    assert!(e.contains("gross transfer"), "{e}");
    // relaxation: residual = 10 % of the volume moved but only 0.1 % of the stock -> must fail
    let r = pf_talus::budget::relax_throughput(10_000.0, 10_000.0 - 1.0 - 9.0, 9.0, 1.0);
    assert!(r.rel_stock < 1e-3 && (r.rel_moved - 0.1).abs() < 1e-12);
    assert!(pf_talus::budget::check_relax(&r, 1e-3).is_err());
    // and a stock-only check would have passed it
    assert!(r.rel_stock <= 1e-3);
}
