//! T1b regional relief windows (ADR 043, node TER-REGION; lane TR).
//! Labelled pre-signoff; SimTime t1b-relaxation-not-history.
//!
//! Usage:
//!   pf-region-relief run <t0_run_dir> <subgrid_relief.f32> <site.json> <t1a_skeleton_dir> <out_dir> [key=value ...] <rx,ry> [<rx,ry> ...]
//!       keys: realization, gain_permille, simtime_kyr, dt_kyr, rough_ppm, rough_cap_mm, halo_mm, base_mode, inflow_mode,
//!             litho=1 (exposed lithology + threshold-slope failure, generator 0.6.0; lane T3)
//!             v4=1 (needs guides=: guide-surface ring, nonlinear creep, curved guides, zero-contour coast; generator 0.7.0; lane T3b)
//!             coast=0|1|2 with v4=1: 0 T0 cell edge, 1 zero contour only, 2 zero contour + reclaimed plain + mouth extensions (default)
//!       windows are computed in the order given, in one process
//!       guides=<guides_dir>: chart-wide tributary guides (generator 0.5.0, guide_mode 1; lane TW)
//!   pf-region-relief guides <t0_run_dir> <subgrid_relief.f32> <site.json> <t1a_skeleton_dir> <out_dir> [threshold_m2=N] [simtime_kyr=N] [gain_permille=N] [floor_mode=0|1] [confine_mode=0|1]
//!       builds the chart-wide t1b.tributary_guides layer once per chart (512 m, lane TW)
//!       preset=v3 (first argument after out_dir): guides 0.4.0, 16 km2 threshold + lithology law (lane T3); later keys override
//!   pf-region-relief net <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> [key=value ...] <tx,ty> [...]
//!       T1b 0.8.0 stage 1 (lane NO): one 64 m lattice per lod-8 tile; writes <out_dir>/t_<tx>_<ty>/net.json + fields
//!       keys: halo_mm, dt_kyr, warm, every, updates, final, gain0, a_ch, realization, fill_halo_mm, fill_steps, fill_dt_kyr
//!   pf-region-relief gains <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir>      (PREREG REVISION 1 stage 0, chart-wide per-cell gains)
//!   pf-region-relief net ... rev=1 gains=<gains_dir> <tx,ty>...                       (REVISION 1 tiles)
//!   pf-region-relief fill <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> net=<net_out_dir> <rx,ry> [...]
//!       T1b 0.8.0 stage 2: 32 m local fill of region (rx, ry) over its owning tile (params read from the tile)
//!   pf-region-relief rivernet <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> [strata=<final_dir>] [box=x0,y0,x1,y1 (m)]... [key=value ...]
//!       T1b 0.9.0 chart layer (lane RV): river network (512 m), E solve (128 m), publish <out_dir>/rivernet.json + fields
//!   pf-region-relief rnwin <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> net=<layer_dir> [strata=<final_dir>] [halo_mm=N] [margin_mm=N] [key=value ...] <rx,ry> [...]
//!       T1b 0.9.0 windows from the layer (params must match the layer's recipe)
//!       T1b 0.9.1 (lane CO) on both: coast=1 (shelf, per-cell flooding offset for F4, fans, wave cut, E taper; keys shelf_mm,
//!       taper_mm, taper_min_permille, flood_init_mm, flood_max_mm, flood_gain_mm, fan_max_mm, cliff_retreat_mm), volcano=1 with
//!       volcanoes=<volcanoes.json> (edifices by kind and record age; key clear_mm), f6_ref=2 (T0 cell mean on coastal cells)
//!       T1b 0.9.2 (lane CF) on both: cf=1 (with coast=1 volcano=1): scaled cone with foot blend, flows, gullies, foot clearance,
//!       attached fans, T1a shelf grade, damped E/F solve, cliff profile; single items off/overridden by cf_cone_hw, cf_clear_foot,
//!       cf_flows, cf_gullies, cf_fan_nocut, cf_shelf, cf_damp, cf_cliff, cf_softk_mm, cf_taper, cf_bury_inner, cf_bury_outer
//!       T1b 0.9.3 (lane TD) on both: fx=1 (with cf=1): capture routing, soft floors, hillslope-length mean over hillslope posts,
//!       contact step, soft-rock bluff, over / volcano_id publication; single items fx_lh, fx_capture, fx_capture_passes,
//!       fx_kappa_permille, fx_floor, fx_contact, fx_bluff, fx_publish. diag_dump=<dir> (rnwin) writes unpublished per-stage arrays.
//!       T1b 0.9.4 (lane TA) on both: ls=1 (with fx=1): consistent level solve (replaces TD2 capture), T1a bank ramp, Broyden
//!       E update; single items ls_solve, ls_passes, ls_stop_ppm, ls_bank, ls_e. Env VZ_SOLVE_TRACE / VZ_SOLVE_PROBE=c1,c2
//!       (rivernet) print per-evaluation solve rows and an end-state E response probe (stdout only).
//!   pf-region-relief seam <window_dir_a> <window_dir_b> <out_json> [<blend_out.i32>]
//!   pf-region-relief f12 <window_dir_a> <window_dir_b> <out_json>
//!   pf-region-relief blend <rx> <ry> <out.i32> <window_dir>...   (region core posts; raw where a neighbour is missing, reported)

use pf_region::t0::{SiteRecord, T0};
use pf_region::t1b::*;
use pf_region::t1b_checks::{self, check_window, relief_stats, seam_compare};
use pf_region::t1b_view::SkeletonView;
use pf_region::t1b_window::{self, run_window};
use std::path::PathBuf;
use std::time::Instant;

fn peak_rss_kib() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/self/status").ok()?;
    s.lines().find(|l| l.starts_with("VmHWM:"))?.split_whitespace().nth(1)?.parse().ok()
}

fn die(msg: String) -> ! { eprintln!("{msg}"); std::process::exit(1) }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(String::as_str) {
        Some("run") => cmd_run(&a[2..]),
        Some("seam") => cmd_seam(&a[2..]),
        Some("f12") => cmd_f12(&a[2..]),
        Some("blend") => cmd_blend(&a[2..]),
        Some("guides") => cmd_guides(&a[2..]),
        Some("net") => cmd_net(&a[2..]),
        Some("fill") => cmd_fill(&a[2..]),
        Some("gains") => cmd_gains(&a[2..]),
        Some("rivernet") => cmd_rivernet(&a[2..]),
        Some("rnwin") => cmd_rnwin(&a[2..]),
        Some("surface-build") if a.get(2).map(String::as_str)==Some("hash") => {
            if a.len()!=4 {die("surface-build hash <artifact>".into());}
            let bytes=std::fs::read(&a[3]).unwrap_or_else(|e|die(e.to_string()));println!("{}",blake3::hash(&bytes).to_hex());
        },
        Some("surface-build") if a.get(2).map(String::as_str)==Some("material") => cmd_rnwin_impl(&a[3..],true,false,true),
        Some("surface-build") if a.get(2).map(String::as_str)==Some("native") => cmd_rnwin_impl(&a[3..],true,false,false),
        Some("surface-build") if a.get(2).map(String::as_str)==Some("channels") => cmd_rnwin_impl(&a[3..],true,true,false),
        Some("surface-build") => {
            if a.len()!=4 {die("surface-build <content-bound-stage-request.json> <new-candidate-directory>".into());}
            pf_region::t1b_surface_domain::build_from_request(std::path::Path::new(&a[2]),std::path::Path::new(&a[3])).unwrap_or_else(|e|die(e.to_string()));
        },
        _ => { eprintln!("usage: see the file header of src/bin/pf-region-relief.rs"); std::process::exit(2) }
    }
}

fn cmd_run(a: &[String]) {
    if a.len() < 6 { die("run: <t0_run_dir> <subgrid> <site> <skeleton_dir> <out_dir> [k=v] <rx,ry>...".into()); }
    let (run, sg, site, skd, out) = (PathBuf::from(&a[0]), PathBuf::from(&a[1]), PathBuf::from(&a[2]), PathBuf::from(&a[3]), PathBuf::from(&a[4]));
    let mut p = ReliefParamsV1::default();
    let mut windows = Vec::new();
    let mut guides_dir: Option<PathBuf> = None;
    let mut litho = false;
    let mut v4 = false;
    let mut coast: Option<u32> = None;
    for s in &a[5..] {
        if let Some(d) = s.strip_prefix("guides=") { guides_dir = Some(PathBuf::from(d)); continue; }
        if let Some((k, v)) = s.split_once('=') {
            let iv: i64 = v.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            match k {
                "realization" => p.roughness_realization = iv as u32,
                "gain_permille" => p.uplift_gain_permille = iv,
                "simtime_kyr" => p.sim_time_kyr = iv,
                "dt_kyr" => p.dt_kyr = iv,
                "rough_ppm" => p.roughness_rms_fraction_ppm = iv,
                "rough_cap_mm" => p.roughness_rms_cap_mm = iv,
                "halo_mm" => p.halo_mm = iv,
                "base_mode" => p.initial_base_mode = iv as u32,
                "inflow_mode" => p.inflow_mode = iv as u32,
                "litho" => litho = iv == 1,
                "d_ppm" => p.hillslope_d_ppm_of_t0 = Some(iv),
                "v4" => v4 = iv == 1,
                "coast" => coast = Some(iv as u32),
                _ => die(format!("unknown key {k}")),
            }
        } else {
            let (x, y) = s.split_once(',').unwrap_or_else(|| die(format!("bad window {s}")));
            windows.push((x.parse::<i64>().unwrap_or_else(|_| die(format!("bad rx {s}"))), y.parse::<i64>().unwrap_or_else(|_| die(format!("bad ry {s}")))));
        }
    }
    let t_all = Instant::now();
    let site = SiteRecord::load(&site).unwrap_or_else(|e| die(format!("site: {e}")));
    let t = Instant::now();
    let t0 = T0::load(&run, &sg, &site, p.threshold_m3_s()).unwrap_or_else(|e| die(format!("t0: {e}")));
    let extras = T0Extras::load(&run, &t0).unwrap_or_else(|e| die(format!("t0 extras: {e}")));
    let skel = SkeletonView::load(&skd).unwrap_or_else(|e| die(format!("t1a: {e}")));
    if skel.climate_digest.is_some() { die("the T1a skeleton carries climate runoff; only rivernet / rnwin read it (lane CLIMATE-RUNOFF)".into()); }
    let gv = guides_dir.as_ref().map(|d| GuideView::load(d).unwrap_or_else(|e| die(format!("guides: {e}"))));
    if gv.is_some() { p = p.with_guides(); }
    if v4 { if gv.is_none() { die("v4=1 needs guides=<dir>".into()); } p = p.with_v4(); if let Some(c) = coast { p = p.with_coast(c); } }
    else if litho { p = p.with_lithology(); }
    let load_s = t.elapsed().as_secs_f64();
    let inp = ReliefInputs { t0: &t0, extras: &extras, site: &site, skel: &skel, guides: gv.as_ref(), params: p.clone() };
    let recipe = recipe_for(&inp).unwrap_or_else(|e| die(format!("recipe: {e}")));
    let rd = recipe.digest().unwrap_or_else(|e| die(e));
    std::fs::create_dir_all(&out).unwrap_or_else(|e| die(e.to_string()));
    std::fs::write(out.join("recipe.json"), serde_json::to_vec_pretty(&serde_json::json!({"recipe": recipe, "digest": rd})).unwrap()).unwrap();
    let chart = vz_tilekey::chart::Chart::new(recipe.chart.clone()).unwrap();
    let t = Instant::now();
    let zc = t1b_checks::t0_zero_contour(&t0, &chart).unwrap_or_else(|e| die(e));
    let zc_s = t.elapsed().as_secs_f64();
    let mut exec = serde_json::Map::new();
    exec.insert("load_t0_t1a_verify_seconds".into(), load_s.into());
    exec.insert("t0_zero_contour_seconds".into(), zc_s.into());
    for (rx, ry) in windows {
        let wr = run_window(&inp, &recipe, &rd, rx, ry).unwrap_or_else(|e| die(format!("window {rx},{ry}: {e}")));
        let wdir = out.join(format!("w_{rx}_{ry}"));
        t1b_window::publish(&wr, &wdir).unwrap_or_else(|e| die(e));
        let t = Instant::now();
        let ch = check_window(&wr.raw, &wr.prep, &wr.header.ledger, &t0, &skel, &p, &recipe.chart, Some(&zc)).unwrap_or_else(|e| die(e));
        let check_s = t.elapsed().as_secs_f64();
        std::fs::write(wdir.join("checks.json"), serde_json::to_vec_pretty(&ch).unwrap()).unwrap();
        println!("window {} digest {} recipe {}", wr.raw.lat.key(), wr.raw.digest, rd);
        println!("  F1 strict {} adr {} (off-channel run {} posts) (posts {}, leave {}, rise {}, wrong outlet {}, {:?})  F2 {} (max {:.2} m3/s)  F3 {}  F4 applicable {} c1 {} c2 {} ({:.1} km)  F5 pits {} (>1 m {}, max {:.2} m)  F6 {}  F7 {} ({:.2e}, vol {:.2e}, flux {:.2e})",
            ch.f1_pass, ch.f1_adr_reading_pass, ch.f1_max_offchannel_run_posts, ch.f1_core_channel_posts, ch.f1_paths_leaving_channel, ch.f1_bed_rises, ch.f1_paths_to_wrong_outlet, ch.f1_outlets,
            ch.f2_pass, ch.f2_max_free_discharge_m3_s, ch.f3_pass, ch.f4_applicable, ch.f4_clause1_pass, ch.f4_clause2_pass,
            ch.f4_max_coast_to_t0_zero_contour_km, ch.f5_core_pit_posts, ch.f5_core_pit_posts_deeper_1m, ch.f5_core_pit_max_depth_m,
            ch.f6_pass, ch.f7_pass, ch.f7_relative_residual, ch.f7_volume_relative_residual, ch.f7_diffusion_flux_relative_residual);
        for c in &ch.cells {
            println!("  cell {} {} posts {} T1b mean {:.1} bilinear {:.1} diff {:.1} tol {:.0} F6 {} | T0 mean {:.1} literal diff {:.1} | F3 {:.1} permille (strict {:.1}) | coastal {} F4 whole-cell {:.3} vs {:.3}",
                c.cell, c.water, c.core_posts, c.t1b_mean_m, c.bilinear_mean_m, c.f6_bilinear_diff_m, c.f6_tolerance_m, c.f6_pass,
                c.t0_mean_m, c.f6_literal_cell_mean_diff_m, c.f3_share_permille, c.f3_strict_share_permille, c.coastal,
                c.f4_whole_cell_t1b_land_fraction, c.f4_whole_cell_bilinear_land_fraction);
        }
        let s = &ch.stats;
        println!("  stats mean {:.1} std {:.1} p05/50/95 {:.0}/{:.0}/{:.0} HI {:.3} slope mean {:.1} deg p90 {:.1} >30 {:.3} >45 {:.3} local relief {:.0} m valley {:.4}",
            s.mean_m, s.std_m, s.p05_m, s.p50_m, s.p95_m, s.hypsometric_integral, s.slope_mean_deg, s.slope_p90_deg,
            s.slope_frac_gt_30deg, s.slope_frac_gt_45deg, s.local_relief_mean_m, s.valley_post_fraction);
        let ps = &wr.header.prep;
        if let Some(l) = &ps.lithology { println!("  lithology: basement code counts {:?} cover posts {} cover max {:.0} m Sc {:.3}-{:.3} k {:.3}-{:.3}; failure exported {:.3e} m3", l.basement_code_counts, l.cover_posts, l.cover_thickness_max_m, l.basement_sc_min, l.basement_sc_max, l.basement_k_min, l.basement_k_max, wr.header.ledger.failure_eroded_m3.unwrap_or(0.0)); }
        if let Some(g) = &ps.guides { println!("  guides: links {} posts {} conflicts {} (max {:.2} m) under channel/water {}", g.guide_links, g.guide_posts, g.guide_post_conflicts, g.guide_post_conflict_max_m, g.guide_posts_under_channel_or_water); }
        println!("  prep: channel posts {} conflicts {} soft axes {} marine {} free {} ring(free) {} inflow edges {} by side W/E/S/N km2 {:?} unplaced {:.1} | U core {:.1} m/My (T0 uplift {:.1})",
            ps.channel_posts, ps.channel_post_conflicts, ps.soft_axis_nodes_in_window, ps.marine_posts, ps.free_posts, ps.ring_posts, ps.inflow_edges,
            ps.inflow_area_by_side_km2.iter().map(|v| v.round()).collect::<Vec<_>>(), ps.inflow_area_unplaced_km2, ps.uplift_core_mean_m_per_my, ps.t0_uplift_core_mean_m_per_my);
        let last = wr.header.series.last().unwrap();
        println!("  relax {:.1} s prepare {:.1} s checks {:.1} s; last series t {:.1} My core mean {:.1} mean|rate| {:.2} m/My max|rate| {:.1}",
            wr.seconds_relax, wr.seconds_prepare, check_s, last.t_relax_my, last.core_free_mean_m, last.core_mean_abs_rate_m_per_my, last.core_max_abs_rate_m_per_my);
        exec.insert(format!("window_{rx}_{ry}"), serde_json::json!({"prepare_seconds": wr.seconds_prepare, "relax_seconds": wr.seconds_relax, "checks_seconds": check_s}));
    }
    exec.insert("total_seconds".into(), t_all.elapsed().as_secs_f64().into());
    exec.insert("peak_rss_kib".into(), peak_rss_kib().into());
    exec.insert("not_digested".into(), "execution receipt; not part of any window digest".into());
    std::fs::write(out.join("execution.json"), serde_json::to_vec_pretty(&serde_json::Value::Object(exec.clone())).unwrap()).unwrap();
    println!("{}", serde_json::Value::Object(exec));
}

fn cmd_seam(a: &[String]) {
    if a.len() < 3 { die("seam: <dir_a> <dir_b> <out_json> [blend_out.i32]".into()); }
    let (wa, ha) = t1b_window::load(&PathBuf::from(&a[0])).unwrap_or_else(|e| die(e));
    let (wb, _) = t1b_window::load(&PathBuf::from(&a[1])).unwrap_or_else(|e| die(e));
    let hw = ha.pointer("/recipe/params/seam_half_width_mm").and_then(|v| v.as_i64()).unwrap_or(512_000) / T1B_SPACING_MM;
    let (s, blended) = seam_compare(&wa, &wb, hw).unwrap_or_else(|e| die(e));
    std::fs::write(&a[2], serde_json::to_vec_pretty(&s).unwrap()).unwrap();
    if let Some(p) = a.get(3) {
        let mut b = Vec::with_capacity(blended.len() * 20);
        for (x, y, v) in blended { b.extend_from_slice(&x.to_le_bytes()); b.extend_from_slice(&y.to_le_bytes()); b.extend_from_slice(&v.to_le_bytes()); }
        std::fs::write(p, b).unwrap();
    }
    println!("{}", serde_json::to_string_pretty(&s).unwrap());
}

fn cmd_f12(a: &[String]) {
    if a.len() < 3 { die("f12: <dir_a> <dir_b> <out_json>".into()); }
    let (wa, _) = t1b_window::load(&PathBuf::from(&a[0])).unwrap_or_else(|e| die(e));
    let (wb, _) = t1b_window::load(&PathBuf::from(&a[1])).unwrap_or_else(|e| die(e));
    let p = ReliefParamsV1::default();
    let (sa, sb) = (relief_stats(&wa, &p), relief_stats(&wb, &p));
    let r = t1b_checks::f12_compare(&sa, &sb, &wa, &wb);
    let doc = serde_json::json!({"a": a[0], "b": a[1], "a_digest": wa.digest, "b_digest": wb.digest, "stats_a": sa, "stats_b": sb, "f12": r});
    std::fs::write(&a[2], serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
    for row in &r.rows { println!("{:24} a {:12.4} b {:12.4} diff {:10.4} {:18} {}", row.stat, row.a, row.b, row.diff, row.tolerance, row.pass); }
    println!("channels fixed identical {}  F12 pass {}", r.channels_fixed_identical, r.pass);
}

fn cmd_blend(a: &[String]) {
    if a.len() < 4 { die("blend: <rx> <ry> <out.i32> <window_dir>...".into()); }
    let rx: i64 = a[0].parse().unwrap();
    let ry: i64 = a[1].parse().unwrap();
    let mut ws = Vec::new();
    let mut hw = 16;
    for d in &a[3..] {
        let (w, h) = t1b_window::load(&PathBuf::from(d)).unwrap_or_else(|e| die(e));
        hw = h.pointer("/recipe/params/seam_half_width_mm").and_then(|v| v.as_i64()).unwrap_or(512_000) / T1B_SPACING_MM;
        ws.push(w);
    }
    let me = ws.iter().find(|w| w.lat.rx == rx && w.lat.ry == ry).unwrap_or_else(|| die("own window missing".into()));
    let (x0, x1, y0, y1) = me.lat.core();
    let refs: Vec<&t1b_window::RawWindow> = ws.iter().collect();
    let mut out = Vec::new();
    let mut missing = 0u64;
    for gy in y0..=y1 { for gx in x0..=x1 {
        // A post needs every window whose seam band covers it.
        let needs: Vec<(i64, i64)> = [-1i64, 0, 1].iter().flat_map(|dy| [-1i64, 0, 1].iter().map(move |dx| (rx + dx, ry + dy)))
            .filter(|&(qx, qy)| seam_weight(&WinLat::new(qx, qy, me.lat.halo_posts), gx, gy, hw) > 0.0).collect();
        let have = needs.iter().all(|k| ws.iter().any(|w| (w.lat.rx, w.lat.ry) == *k));
        let v = if have { t1b_checks::blend_post(&refs, gx, gy, hw).unwrap_or_else(|e| die(e)).unwrap() } else {
            missing += 1; me.height_mm[me.lat.idx(gx, gy).unwrap()] };
        out.extend_from_slice(&v.to_le_bytes());
    } }
    std::fs::write(&a[2], out).unwrap();
    println!("blended core {}x{} posts; {} posts kept raw (a neighbour window was not supplied)", x1 - x0 + 1, y1 - y0 + 1, missing);
}

use pf_region::t1b_checks::seam_weight;

fn cmd_guides(a: &[String]) {
    use pf_region::t1b_guides::*;
    if a.len() < 5 { die("guides: <t0_run_dir> <subgrid> <site> <skeleton_dir> <out_dir> [threshold_m2=N] [simtime_kyr=N]".into()); }
    let (run, sg, site, skd, out) = (PathBuf::from(&a[0]), PathBuf::from(&a[1]), PathBuf::from(&a[2]), PathBuf::from(&a[3]), PathBuf::from(&a[4]));
    let mut gp = GuideParamsV1::default();
    for s in &a[5..] {
        if s == "preset=v3" { gp = GuideParamsV1::v3(); continue; }
        let (k, v) = s.split_once('=').unwrap_or_else(|| die(format!("bad arg {s}")));
        let iv: i64 = v.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
        match k { "threshold_m2" => gp.guide_area_threshold_m2 = iv, "simtime_kyr" => gp.sim_time_kyr = iv, "gain_permille" => gp.uplift_gain_permille = iv,
            "floor_mode" => gp.floor_mode = iv as u32, "confine_mode" => gp.confine_mode = iv as u32, _ => die(format!("unknown key {k}")) }
    }
    let t = Instant::now();
    let site = SiteRecord::load(&site).unwrap_or_else(|e| die(format!("site: {e}")));
    let t0 = T0::load(&run, &sg, &site, ReliefParamsV1::default().threshold_m3_s()).unwrap_or_else(|e| die(format!("t0: {e}")));
    let extras = T0Extras::load(&run, &t0).unwrap_or_else(|e| die(format!("t0 extras: {e}")));
    let skel = SkeletonView::load(&skd).unwrap_or_else(|e| die(format!("t1a: {e}")));
    if skel.climate_digest.is_some() { die("the T1a skeleton carries climate runoff; only rivernet / rnwin read it (lane CLIMATE-RUNOFF)".into()); }
    let g = build_guides(&t0, &extras, &site, &skel, gp).unwrap_or_else(|e| die(e));
    publish_guides(&g, &out).unwrap_or_else(|e| die(e));
    let back = GuideView::load(&out).unwrap_or_else(|e| die(format!("re-read: {e}")));
    println!("guides digest {} recipe {} (re-read ok: {})", g.digest, g.recipe_digest, back.digest == g.digest);
    println!("{}", serde_json::to_string_pretty(&g.stats).unwrap());
    let last = g.series.last().unwrap();
    println!("ledger f7 {:.2e} vol {:.2e}; last series t {:.1} My mean {:.1} mean|rate| {:.2} m/My; seconds {:.1}; peak rss {:?} KiB",
        g.ledger.f7_relative_residual, g.ledger.volume_relative_residual, last.t_relax_my, last.core_free_mean_m, last.core_mean_abs_rate_m_per_my,
        t.elapsed().as_secs_f64(), peak_rss_kib());
}

use pf_region::t1b_guides::GuideView;

fn cpu_seconds() -> Option<f64> {
    // utime + stime of this process (clock ticks, assumed 100/s as on Linux x86-64).
    let s = std::fs::read_to_string("/proc/self/stat").ok()?;
    let rest = &s[s.rfind(')')? + 2..];
    let f: Vec<&str> = rest.split_whitespace().collect();
    Some((f.get(11)?.parse::<f64>().ok()? + f.get(12)?.parse::<f64>().ok()?) / 100.0)
}

fn net_inputs(a: &[String]) -> (T0, T0Extras, SiteRecord, SkeletonView, f64) {
    let t = Instant::now();
    let site = SiteRecord::load(&PathBuf::from(&a[2])).unwrap_or_else(|e| die(format!("site: {e}")));
    let mut t0 = T0::load(&PathBuf::from(&a[0]), &PathBuf::from(&a[1]), &site, ReliefParamsV1::default().threshold_m3_s()).unwrap_or_else(|e| die(format!("t0: {e}")));
    // Lane CLIMATE-RUNOFF: attach the climate runoff field when climate=<dir> is given (rivernet / rnwin only).
    if let Some(d) = a[5.min(a.len())..].iter().find_map(|s| s.strip_prefix("climate=")) {
        let mut cp = pf_region::climate_runoff::ClimateRunoffParamsV1::v1();
        if let Some(f) = a[5..].iter().find_map(|s| s.strip_prefix("climate_floor_permille=")) { cp.floor_permille = f.parse().unwrap_or_else(|_| die(format!("bad climate_floor_permille {f}"))); }
        t0.attach_climate(&PathBuf::from(d), cp.clone()).unwrap_or_else(|e| die(format!("climate: {e}")));
        println!("lane CLIMATE-RUNOFF: climate {} tier {} floor {} permille, r_ref {} m/yr", t0.climate.as_ref().unwrap().digest, t0.climate.as_ref().unwrap().tier, cp.floor_permille, t0.runoff_m_per_yr);
    }
    let extras = T0Extras::load(&PathBuf::from(&a[0]), &t0).unwrap_or_else(|e| die(format!("t0 extras: {e}")));
    let skel = SkeletonView::load(&PathBuf::from(&a[3])).unwrap_or_else(|e| die(format!("t1a: {e}")));
    if skel.climate_digest.is_some() && t0.climate.is_none() { die("the T1a skeleton carries climate runoff: give the same climate=<dir> (rivernet / rnwin only)".into()); }
    (t0, extras, site, skel, t.elapsed().as_secs_f64())
}

fn parse_xy(s: &str) -> (i64, i64) {
    let (x, y) = s.split_once(',').unwrap_or_else(|| die(format!("bad key {s}")));
    (x.parse().unwrap_or_else(|_| die(format!("bad {s}"))), y.parse().unwrap_or_else(|_| die(format!("bad {s}"))))
}

fn cmd_net(a: &[String]) {
    use pf_region::t1b_net::*;
    if a.len() < 6 { die("net: <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> [k=v] <tx,ty>...".into()); }
    let t_all = Instant::now();
    let mut np = NetParamsV1::default();
    let mut keys = Vec::new();
    let mut gains_dir = None;
    if a[5..].iter().any(|s| s == "rev=1") { np = NetParamsV1::rev1(); }
    for s in &a[5..] {
        if s == "rev=1" { continue; }
        if let Some(d) = s.strip_prefix("gains=") { gains_dir = Some(PathBuf::from(d)); continue; }
        if let Some((k, v)) = s.split_once('=') {
            let iv: i64 = v.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            match k {
                "halo_mm" => np.halo_mm = iv, "dt_kyr" => np.dt_kyr = iv, "warm" => np.warm_steps = iv as u32, "every" => np.gain_every_steps = iv as u32,
                "updates" => np.gain_updates_max = iv as u32, "final" => np.final_steps = iv as u32, "gain0" => np.gain0_permille = iv,
                "a_ch" => np.channel_head_area_m2 = iv, "realization" => np.roughness_realization = iv as u32, "fill_halo_mm" => np.fill_halo_mm = iv,
                "fill_steps" => np.fill_steps = iv as u32, "fill_dt_kyr" => np.fill_dt_kyr = iv, "creep_solver" => np.creep_solver = iv as u32,
                _ => die(format!("unknown key {k}")),
            }
        } else { keys.push(parse_xy(s)); }
    }
    let (t0, extras, site, skel, load_s) = net_inputs(a);
    let out = PathBuf::from(&a[4]);
    let mut exec = serde_json::Map::new();
    exec.insert("load_seconds".into(), load_s.into());
    for (tx, ty) in keys {
        let g = gains_dir.as_ref().map(|d| load_gains(d).unwrap_or_else(|e| die(format!("gains: {e}"))));
        let t = build_tile_g(&t0, &extras, &site, &skel, np.clone(), tx, ty, g.as_ref()).unwrap_or_else(|e| die(format!("tile {tx},{ty}: {e}")));
        let dir = out.join(format!("t_{tx}_{ty}"));
        publish_tile(&t, &dir).unwrap_or_else(|e| die(e));
        let s = &t.stats;
        println!("tile {} digest {} recipe {}", t.lat.key(), t.digest, t.recipe_digest);
        println!("  steps {} ({:.1} My) gain updates {} converged {} history max|diff| {:?}", s.steps, s.sim_time_my, s.gain_updates_done, s.gains_converged,
            s.gain_history_max_abs_diff_m.iter().map(|v| v.round()).collect::<Vec<_>>());
        for c in &s.cells { if c.land { println!("  cell {} solved {} mass {:.2} gain {:.3}{} z {:.0} bil {:.0} base {:.0} diff {:.0}", c.cell, c.solved, c.hat_mass_fraction, c.gain, if c.at_bound { " (BOUND)" } else { "" }, c.hat_mean_z_m, c.hat_mean_bilinear_m, c.hat_mean_base_m, c.diff_m); } }
        println!("  core mean {:.0} bilinear {:.0}; net channel posts {} mfd posts {} pits {} unconfined {}; F3@64m min {:.1} permille ({} cells, {} < 950); newton max {} unconverged {}; F7 {:.2e} vol {:.2e}",
            s.core_mean_m, s.core_bilinear_mean_m, s.net_channel_posts_core, s.mfd_posts_last_step_core, s.pits_core_final, s.unconfined_posts_final,
            s.f3_min_share_permille, s.f3_cells_checked, s.f3_cells_below_950, s.creep_newton_iters_max, s.creep_lines_not_converged,
            t.header["ledger"]["f7_relative_residual"].as_f64().unwrap_or(-1.0), t.header["ledger"]["volume_relative_residual"].as_f64().unwrap_or(-1.0));
        println!("  seconds {:?}", t.seconds);
        println!("  newton iterations total {} per post-step-sweep {:.2}", t.newton_iters_total, t.newton_iters_total as f64 / (2.0 * s.steps as f64 * s.free_posts as f64).max(1.0));
        exec.insert(format!("tile_{tx}_{ty}"), serde_json::json!({"seconds": t.seconds}));
    }
    exec.insert("total_seconds".into(), t_all.elapsed().as_secs_f64().into());
    exec.insert("cpu_seconds".into(), cpu_seconds().into());
    exec.insert("peak_rss_kib".into(), peak_rss_kib().into());
    exec.insert("not_digested".into(), "execution receipt".into());
    std::fs::create_dir_all(&out).unwrap();
    let name = format!("execution_net_{}.json", std::process::id());
    std::fs::write(out.join(name), serde_json::to_vec_pretty(&serde_json::Value::Object(exec.clone())).unwrap()).unwrap();
    println!("{}", serde_json::Value::Object(exec));
}

fn cmd_fill(a: &[String]) {
    use pf_region::t1b_net::*;
    if a.len() < 7 { die("fill: <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> net=<dir> <rx,ry>...".into()); }
    let t_all = Instant::now();
    let mut net_dir = None;
    let mut keys = Vec::new();
    for s in &a[5..] { if let Some(d) = s.strip_prefix("net=") { net_dir = Some(PathBuf::from(d)); } else { keys.push(parse_xy(s)); } }
    let net_dir = net_dir.unwrap_or_else(|| die("fill needs net=<dir>".into()));
    let (t0, extras, site, skel, load_s) = net_inputs(a);
    let out = PathBuf::from(&a[4]);
    std::fs::create_dir_all(&out).unwrap();
    let mut exec = serde_json::Map::new();
    exec.insert("load_seconds".into(), load_s.into());
    let chart_zc = std::cell::OnceCell::new();
    for (rx, ry) in keys {
        let (tx, ty) = NetLat::owner_of_region(rx, ry);
        let t = Instant::now();
        let nv = NetView::load(&net_dir.join(format!("t_{tx}_{ty}"))).unwrap_or_else(|e| die(format!("net tile {tx},{ty}: {e}")));
        let load_net_s = t.elapsed().as_secs_f64();
        let p = fill_params(&nv.params);
        let inp = ReliefInputs { t0: &t0, extras: &extras, site: &site, skel: &skel, guides: None, params: p.clone() };
        let wr = run_fill(&inp, &nv, &nv.recipe_digest, rx, ry).unwrap_or_else(|e| die(format!("fill {rx},{ry}: {e}")));
        let wdir = out.join(format!("w_{rx}_{ry}"));
        t1b_window::publish(&wr, &wdir).unwrap_or_else(|e| die(e));
        let zc = chart_zc.get_or_init(|| { let chart = vz_tilekey::chart::Chart::new(wr.header.recipe.chart.clone()).unwrap(); t1b_checks::t0_zero_contour(&t0, &chart).unwrap_or_else(|e| die(e)) });
        let ch = check_window(&wr.raw, &wr.prep, &wr.header.ledger, &t0, &skel, &p, &wr.header.recipe.chart, Some(zc)).unwrap_or_else(|e| die(e));
        std::fs::write(wdir.join("checks.json"), serde_json::to_vec_pretty(&ch).unwrap()).unwrap();
        println!("window {} digest {} recipe {} (tile {})", wr.raw.lat.key(), wr.raw.digest, wr.raw.recipe_digest, nv.lat.key());
        println!("  F1 strict {} adr {} (posts {}, leave {}, rise {}, wrong {})  F2 {} ({:.2})  F3 {}  F4 app {} c1 {} c2 {}  F5 pits {} (>1 m {}, max {:.2})  F6 {}  F7 {} ({:.2e}, vol {:.2e})",
            ch.f1_pass, ch.f1_adr_reading_pass, ch.f1_core_channel_posts, ch.f1_paths_leaving_channel, ch.f1_bed_rises, ch.f1_paths_to_wrong_outlet,
            ch.f2_pass, ch.f2_max_free_discharge_m3_s, ch.f3_pass, ch.f4_applicable, ch.f4_clause1_pass, ch.f4_clause2_pass,
            ch.f5_core_pit_posts, ch.f5_core_pit_posts_deeper_1m, ch.f5_core_pit_max_depth_m, ch.f6_pass, ch.f7_pass, ch.f7_relative_residual, ch.f7_volume_relative_residual);
        for c in &ch.cells { println!("  cell {} {} posts {} mean {:.0} bilinear {:.0} diff {:.0} tol {:.0} F6 {} F3 {:.1}", c.cell, c.water, c.core_posts, c.t1b_mean_m, c.bilinear_mean_m, c.f6_bilinear_diff_m, c.f6_tolerance_m, c.f6_pass, c.f3_share_permille); }
        let s = &ch.stats;
        println!("  stats mean {:.0} slope mean {:.1} p90 {:.1} >30 {:.3} >45 {:.3} local relief {:.0} valley {:.4}", s.mean_m, s.slope_mean_deg, s.slope_p90_deg, s.slope_frac_gt_30deg, s.slope_frac_gt_45deg, s.local_relief_mean_m, s.valley_post_fraction);
        println!("  net {}", wr.header.prep.net.as_ref().map(|v| v.to_string()).unwrap_or_default());
        println!("  seconds: load net {:.1} prepare {:.1} relax {:.1}", load_net_s, wr.seconds_prepare, wr.seconds_relax);
        exec.insert(format!("window_{rx}_{ry}"), serde_json::json!({"load_net_seconds": load_net_s, "prepare_seconds": wr.seconds_prepare, "relax_seconds": wr.seconds_relax}));
    }
    exec.insert("total_seconds".into(), t_all.elapsed().as_secs_f64().into());
    exec.insert("cpu_seconds".into(), cpu_seconds().into());
    exec.insert("peak_rss_kib".into(), peak_rss_kib().into());
    exec.insert("not_digested".into(), "execution receipt".into());
    let name = format!("execution_fill_{}.json", std::process::id());
    std::fs::write(out.join(name), serde_json::to_vec_pretty(&serde_json::Value::Object(exec.clone())).unwrap()).unwrap();
    println!("{}", serde_json::Value::Object(exec));
}

fn cmd_gains(a: &[String]) {
    use pf_region::t1b_net::*;
    if a.len() < 5 { die("gains: <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir>".into()); }
    let t_all = Instant::now();
    let (t0, extras, site, skel, load_s) = net_inputs(a);
    let g = solve_gains(&t0, &extras, &site, &skel, NetParamsV1::rev1()).unwrap_or_else(|e| die(e));
    let out = PathBuf::from(&a[4]);
    publish_gains(&g, &out).unwrap_or_else(|e| die(e));
    println!("gains digest {} recipe {} updates {} converged {} history {} mean solved {}", g.digest, g.header["recipe_digest"], g.header["updates"], g.header["converged"], g.header["history_max_abs_diff_m"], g.header["mean_solved_gain"]);
    for r in &g.rows { if r.land && r.hat_mass_fraction > 0.3 { println!("  cell {} solved {} mass {:.2} gain {:.3}{} z {:.0} bil {:.0} base {:.0} diff {:.0}", r.cell, r.solved, r.hat_mass_fraction, r.gain, if r.at_bound { " (BOUND)" } else { "" }, r.hat_mean_z_m, r.hat_mean_bilinear_m, r.hat_mean_base_m, r.diff_m); } }
    let exec = serde_json::json!({"load_seconds": load_s, "solve_seconds": g.seconds, "total_seconds": t_all.elapsed().as_secs_f64(), "cpu_seconds": cpu_seconds(), "peak_rss_kib": peak_rss_kib(), "not_digested": "execution receipt"});
    std::fs::write(out.join("execution_gains.json"), serde_json::to_vec_pretty(&exec).unwrap()).unwrap();
    println!("{exec}");
}

// ------------------------------------------------------------------ T1b 0.9.0 RIVERNET (lane RV)

struct RnArgs { p: pf_region::t1b_rivernet::RiverNetParamsV1, volcanoes: Option<PathBuf>, strata: Option<PathBuf>, net: Option<PathBuf>, boxes: Vec<[f64; 4]>, halo_mm: i64, margin_mm: Option<i64>, diag_dx_mm: Option<i64>, keys: Vec<(i64, i64)>,
    /// Lane PV: prov=1 publishes prov_id + header.provenance (rnwin only); prov_strata= a StratigraphyDomainV0 JSON when the strata dir has no env/age dumps.
    unit_contacts: bool,
    prov: bool, prov_json: Option<PathBuf>,
    /// Lane CLIMATE-RUNOFF: climate=<seed-climate-v1 artifact dir> (opt-in; T1a must be built with the same key).
    climate: Option<PathBuf>,
    /// Lane GLACIAL-T1: glacial=<t0.glacial dir> (opt-in, layer and windows; params set after T0 is loaded).
    glacial: Option<PathBuf>,
    /// Lane MS-CEIL-2: ceil_ped_ref=<key-off layer dir> (rnwin only; needs ceil and ped_base): the cone's local ground G is read
    /// from that layer's lane_ped instead of this layer's, so the ceiling keys do not move the built cone (a declared coupling).
    ped_ref: Option<PathBuf>,
    /// Lane GLACIER-1: gl_now=1 (with glacial=<present-ice-only t0.glacial>): present-day ice (NOW state, icefield sheet, window ice fields).
    gl_now: bool,
    /// Lane GLACIER-2: gl_bal=1 (with gl_now=1): the balance-flux icefield, calving and the 32 m margin (ice_window_v2).
    gl_bal: bool,
    /// Lane GLACIER-3: gl_geo=1 (with gl_bal=1): FMM plastic sheet, D-infinity routing, aligned marine claims, converged passes, ice_sub_mm.
    gl_geo: bool }

fn rn_args(a: &[String]) -> RnArgs {
    use pf_region::t1b_rivernet::RiverNetParamsV1;
    let strata = a[5..].iter().find_map(|s| s.strip_prefix("strata=")).map(PathBuf::from);
    let mut r = RnArgs { p: RiverNetParamsV1::v1(strata.is_some()), volcanoes: None, strata, net: None, boxes: Vec::new(), halo_mm: 2_048_000, margin_mm: None, diag_dx_mm: None, keys: Vec::new(), unit_contacts: false, prov: false, prov_json: None, climate: None, glacial: None, gl_now: false, gl_bal: false, ped_ref: None, gl_geo: false };
    let mut climate_floor: Option<i64> = None;
    // Lane CO: coast=1 / volcano=1 switch the 0.9.1 additions on (v1 constants); their key overrides apply afterwards.
    let mut co_over: Vec<(String, i64)> = Vec::new();
    // Lane CF: cf=1 switches the 0.9.2 additions on (v1 constants, needs coast=1 volcano=1); cf_* keys override single items.
    let mut cf_over: Vec<(String, i64)> = Vec::new();
    // Lane TD: fx=1 switches the 0.9.3 fixes on (v1 constants, needs cf=1); fx_* keys override single items.
    let mut fx_over: Vec<(String, i64)> = Vec::new();
    // Lane TA: ls=1 switches the 0.9.4 additions on (v1 constants, needs fx=1); ls_* keys override single items.
    let mut ls_over: Vec<(String, i64)> = Vec::new();
    // Lane VALLEY-FLOOR-FIX: vf_floor=1 (levee-backswamp floor, footslope on unless vf_foot=0), vf_foot=0|1, vf_fan=1.
    let mut vf_over: Vec<(String, i64)> = Vec::new();
    // Lane HS-PASS: hs=1 (hillslope transient, v0 constants) with hs_* overrides; hs_volcano=1|2 (CF young gullies off).
    let mut hs_over: Vec<(String, i64)> = Vec::new();
    for s in &a[5..] {
        if s.starts_with("strata=") { continue; }
        if let Some(d) = s.strip_prefix("net=") { r.net = Some(PathBuf::from(d)); continue; }
        if let Some(d) = s.strip_prefix("volcanoes=") { r.volcanoes = Some(PathBuf::from(d)); continue; }
        if let Some(d) = s.strip_prefix("unit_contacts=") { r.unit_contacts = match d { "0" => false, "1" => true, _ => die(format!("bad value {s}: expected 0 or 1")) }; continue; }
        if let Some(d) = s.strip_prefix("prov=") { r.prov = d == "1"; continue; }
        if let Some(d) = s.strip_prefix("prov_strata=") { r.prov_json = Some(PathBuf::from(d)); continue; }
        if let Some(d) = s.strip_prefix("climate=") { r.climate = Some(PathBuf::from(d)); continue; }
        if let Some(d) = s.strip_prefix("glacial=") { r.glacial = Some(PathBuf::from(d)); continue; }
        if let Some(d) = s.strip_prefix("gl_now=") { r.gl_now = d == "1"; continue; } // Lane GLACIER-1
        if let Some(d) = s.strip_prefix("gl_bal=") { r.gl_bal = d == "1"; continue; } // Lane GLACIER-2
        if let Some(d) = s.strip_prefix("gl_geo=") { r.gl_geo = d == "1"; continue; } // Lane GLACIER-3
        if s.starts_with("surface_owner=") { continue; }
        if s.starts_with("ve_flank2=") || s.starts_with("ve_flank2_ref=") || s.starts_with("ve_flank2_cache=") { continue; }
        if s.starts_with("ve_flank=") { continue; } // CX-CONE-STARBURST: bound with ve_erec after T0 load.
        if s.starts_with("ve_erec=") { continue; } // CX-EDIFICE-EREC: bind after T0 load in both commands.
        if s.starts_with("ring_plain=") || s.starts_with("ve_wind=") { continue; } // Lane VH-I (applied after the loop)
        // Lane T0-IMPRINT-PROD: conf=cell|basin|none|warp (+ conf_warp_mm=, conf_warp_spacing_mm=; layer AND windows) and outlet_checks=1 (checks only; ADR 045).
        if pf_region::t1b_outlet::parse_key(s, &mut r.p).unwrap_or_else(|e| die(e)) { continue; }
        // Lane CHANNEL-CENTRELINE: cl=1 (T1a-curve centreline), cl_mn=1 (declared decorative meander belt; needs cl=1),
        // cl_k_min / cl_k_max (sinuosity permille), cl_belt (belt fill permille), cl_lambda (wavelength widths permille).
        if let Some(d) = s.strip_prefix("cl=") {
            match d { "1" => { if r.p.cl.is_none() { r.p.cl = Some(pf_region::t1b_cl::ClParamsV1::v1()); } } "0" => r.p.cl = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| k.starts_with("cl_")) {
            let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            let Some(c) = r.p.cl.as_mut() else { die(format!("{k} needs cl=1 given before it")) };
            match k {
                "cl_mn" => { if v == 1 { *c = c.clone().with_meander(); } else if v != 0 { die("cl_mn is 0 or 1".into()); } }
                "cl_k_min" => c.sinuosity_min_permille = v, "cl_k_max" => c.sinuosity_max_permille = v,
                "cl_belt" => c.belt_fill_permille = v, "cl_lambda" => c.wavelength_widths_permille = v,
                k if k.starts_with("cl_mq") => pf_region::t1b_mq::parse_key(c, k, v).unwrap_or_else(|e| die(e)), // lane MEANDER-Q
                _ => die(format!("unknown key {k}")),
            }
            continue;
        }
        // Lane MS-CEIL: ceil=<c_h permille> (ceiling-min v1, MS-NEXT; layer AND windows): free non-floor land posts take min(z_ss, T0 bilinear + c_h x subgrid relief).
        if let Some(d) = s.strip_prefix("ceil=") { let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}"))); r.p.ceil = Some(pf_region::t1b_rivernet::CeilParamsV1::v1(v)); continue; }
        // Lane MS-CEIL-2: ceil_bed=1 (ceiling_bed_grade_v1) / ceil_ped=1 (ceiling_free_pedestal_v1), sub-keys of ceil (give ceil= before them; layer AND windows).
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| *k == "ceil_bed" || *k == "ceil_ped") { let v: u32 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}"))); let Some(c) = r.p.ceil.as_mut() else { die(format!("{k} needs ceil= given before it")) }; if k == "ceil_bed" { c.set_bed(v) } else { c.set_ped(v) }; continue; }
        if let Some(d) = s.strip_prefix("ceil_ped_ref=") { r.ped_ref = Some(PathBuf::from(d)); continue; }
        // Lane E-SOLVE-CONVERGE: esc=1 (needs es=1 before it; layer AND windows): convergent E-solve driver (t1b_esc.rs); esc_<item>=<int> overrides.
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| *k == "esc" || k.starts_with("esc_")) { let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}"))); let Some(c) = r.p.es.as_mut() else { die(format!("{k} needs es=1 given before it")) }; pf_region::t1b_esc::parse_key(c, k, v).unwrap_or_else(|e| die(e)); continue; }
        // Lane E-SOLVE-GENERALITY: es=1 (secant E update, chart-rim cells, bound-KKT stop; layer AND windows); es_update,
        // es_cov (coverage_min permille), es_rim (rim_mode), es_kkt, es_gmin / es_gmax / es_step (permille) override single items (need es=1 before);
        // es_emin sets the E lower bound e_min_m_per_my (m/My; a plain recipe param).
        if let Some(d) = s.strip_prefix("es=") {
            match d { "1" => { if r.p.es.is_none() { r.p.es = Some(pf_region::t1b_es::EsParamsV1::v1()); } } "0" => r.p.es = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| k.starts_with("es_")) {
            let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            if k == "es_emin" { if r.p.es.is_none() { die("es_emin needs es=1 given before it".into()); } r.p.e_min_m_per_my = v; continue; }
            let Some(c) = r.p.es.as_mut() else { die(format!("{k} needs es=1 given before it")) };
            match k {
                "es_update" => c.update = v as u32, "es_cov" => c.coverage_min_permille = v, "es_rim" => c.rim_mode = v as u32, "es_target" => c.target = v as u32, "es_kkt" => c.kkt = v as u32,
                "es_gmin" => c.g_min_permille = v, "es_gmax" => c.g_max_permille = v, "es_step" => c.step_max_permille = v,
                _ => die(format!("unknown key {k}")),
            }
            continue;
        }
        // Lane T1B-PEDESTAL: ped_base=1 (layer AND windows): composite cones stand on the layer-published local ground.
        if let Some(d) = s.strip_prefix("ped_base=") {
            let q = r.p.pedestal.take().unwrap_or_else(pf_region::t1b_ped::PedestalParamsV1::v1);
            match d { "1" => r.p.pedestal = Some(q.with_base(true)), "0" => { let q = q.with_base(false); r.p.pedestal = if q.any() { Some(q) } else { None }; } _ => die(format!("bad value {s}")) }
            continue;
        }
        // Lane VOLCANIC-SPLIT: ped_relief=1 (after ped_base=1; layer AND windows): the cone keeps the record relief over G.
        if let Some(d) = s.strip_prefix("ped_relief=") {
            let Some(q) = r.p.pedestal.take() else { die("ped_relief needs ped_base=1 given before it".into()) };
            match d { "1" => r.p.pedestal = Some(q.with_relief(true)), "0" => r.p.pedestal = Some(q.with_relief(false)), _ => die(format!("bad value {s}")) }
            continue;
        }
        // Lane T1B-PEDESTAL: ped_solve=1 (layer AND windows): the E solve counts edifice posts at their pre-edifice ground.
        if let Some(d) = s.strip_prefix("ped_solve=") {
            let q = r.p.pedestal.take().unwrap_or_else(pf_region::t1b_ped::PedestalParamsV1::v1);
            match d { "1" => r.p.pedestal = Some(q.with_solve(true)), "0" => { let q = q.with_solve(false); r.p.pedestal = if q.any() { Some(q) } else { None }; } _ => die(format!("bad value {s}")) }
            continue;
        }
        // Lane VOLCANIC-SPLIT: vs=1 (layer AND windows): lava / volcaniclastic packages inside lithology-3 strata layers
        // (strata= source only); sub-keys after vs=1: vs_scope=1|2, vs_sc= (rock-mass tangent permille), vs_face= (lava face
        // tangent), vs_k= / vs_lava_k=, vs_riser_mm=, vs_contrast=, vs_lens_mid= / vs_lens_amp= (permille), vs_lens_mm=lo,hi,
        // vs_lava_t= / vs_vc_t=median_mm,sigma_permille,min_mm,max_mm, vs_start=, vs_ct=0|1.
        if let Some(d) = s.strip_prefix("vs=") {
            match d { "1" => r.p.lithology.volcanic_split = Some(pf_region::t1b_volc::VolcanicSplitV1::v1()), "0" => r.p.lithology.volcanic_split = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| k.starts_with("vs_")) {
            let Some(v) = r.p.lithology.volcanic_split.as_mut() else { die(format!("{k} needs vs=1 given before it")) };
            let pair = |d: &str| -> [i64; 2] { let w: Vec<i64> = d.split(',').map(|x| x.parse().unwrap_or_else(|_| die(format!("bad value {s}")))).collect(); if w.len() != 2 { die(format!("bad value {s}")) } [w[0], w[1]] };
            let one = |d: &str| -> i64 { d.parse().unwrap_or_else(|_| die(format!("bad value {s}"))) };
            match k {
                "vs_scope" => v.scope = one(d) as u32, "vs_sc" => v.mass_sc_permille = one(d), "vs_face" => v.lava_face_sc_permille = one(d),
                "vs_k" => v.vc_k_permille = one(d), "vs_lava_k" => v.lava_k_permille = one(d), "vs_start" => v.lava_start_permille = one(d),
                "vs_riser_mm" => v.riser_min_mm = one(d), "vs_contrast" => v.contrast_min_permille = one(d),
                "vs_lens_mid" => v.lens_mid_permille = one(d), "vs_lens_amp" => v.lens_amp_permille = one(d), "vs_lens_mm" => v.lens_lambda_mm = pair(d),
                "vs_lava_t" | "vs_vc_t" => { let w: Vec<i64> = d.split(',').map(|x| x.parse().unwrap_or_else(|_| die(format!("bad value {s}")))).collect(); if w.len() != 4 { die(format!("bad value {s}")) } let t = [w[0], w[1], w[2], w[3]]; if k == "vs_lava_t" { v.lava_t = t } else { v.vc_t = t } }
                "vs_ct" => { v.ct = one(d) as u32; v.ct_rule = if v.ct == 1 { Some(pf_region::t1b_volc::CT_RULE.into()) } else { None }; }
                // Lane SPLIT-T1B: vs_ct_bed=1 / vs_ct_blend=1 (ct remap, windows) and vs_aa=1 (layer AND windows), vs_band_mm=, vs_aa_dx_mm=.
                "vs_ct_bed" => v.ct_bed = if one(d) == 1 { Some(pf_region::t1b_volc::CT_BED_RULE.into()) } else { None },
                "vs_ct_blend" => v.ct_blend = if one(d) == 1 { Some(pf_region::t1b_volc::CT_BLEND_RULE.into()) } else { None },
                "vs_aa" => v.aa = if one(d) == 1 { Some(pf_region::t1b_volc::AA_RULE.into()) } else { None },
                "vs_band_mm" => v.band_mm = Some(one(d)), "vs_aa_dx_mm" => v.aa_dx_mm = Some(one(d)),
                _ => die(format!("unknown key {k}")),
            }
            continue;
        }
        // E3: content-pinned source envelope; env=0 is exactly absent.
        if let Some(d) = s.strip_prefix("env=") {
            r.p.lithology.envelope = if d == "0" { None } else { Some(pf_region::t1b_litho::EnvelopeConfig::from_path(d).unwrap_or_else(|e| die(e))) };
            continue;
        }
        // Lane LITHO-CONTINUITY: lcont=1 (layer AND windows; strata= source only): planet-wide rock-body identity and continuous
        // contacts across T0 cells (t1b_lcont); lcont_min_mm= (after lcont=1) the absorbed-unit thickness.
        if let Some(d) = s.strip_prefix("lcont=") {
            match d { "1" => r.p.lithology.continuity = Some(pf_region::t1b_lcont::LithoContinuityV1::v1()), "0" => r.p.lithology.continuity = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some(d) = s.strip_prefix("lcont_min_mm=") {
            let Some(c) = r.p.lithology.continuity.as_mut() else { die("lcont_min_mm needs lcont=1 given before it".into()) };
            c.min_unit_mm = d.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            continue;
        }
        // Lane LITHO-AGE-GATE: lcont_record=<stratigraphy.json> (after lcont=1; layer AND windows): age / environment gate on the
        // body identity from the T0 record (digest in params); lcont_age_tol_kyr=, lcont_env_tol= (after lcont_record=).
        if let Some(d) = s.strip_prefix("lcont_record=") {
            let Some(c) = r.p.lithology.continuity.as_mut() else { die("lcont_record needs lcont=1 given before it".into()) };
            c.record = Some(pf_region::t1b_lcont::LcontRecordV1::from_path(d).unwrap_or_else(|e| die(e)));
            continue;
        }
        // Lane CX-LITHO-M2O: conservative ordered subdivision, requires the record gate (layer AND windows).
        if let Some(d) = s.strip_prefix("lcont_m2o=") {
            let Some(g) = r.p.lithology.continuity.as_mut().and_then(|c| c.record.as_mut()) else { die("lcont_m2o needs lcont_record= given before it".into()) };
            g.m2o_rule = match d { "1" => Some(pf_region::t1b_lcont::M2O_RULE.into()), "0" => None, _ => die(format!("bad value {s}")) };
            continue;
        }
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| *k == "lcont_age_tol_kyr" || *k == "lcont_env_tol") {
            let Some(g) = r.p.lithology.continuity.as_mut().and_then(|c| c.record.as_mut()) else { die(format!("{k} needs lcont_record= given before it")) };
            let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            if k == "lcont_age_tol_kyr" { g.age_tol_kyr = Some(v) } else { g.env_tol = u32::try_from(v).unwrap_or_else(|_| die(format!("bad value {s}"))) }
            continue;
        }
        // Lane DAMMED-VALLEY: dam_lake=1 (window-only edifice-dammed lakes), dam_min_depth_mm=<mm> (after dam_lake=1).
        if let Some(d) = s.strip_prefix("dam_lake=") {
            match d { "1" => r.p.dam_lake = Some(pf_region::t1b_dam::DamLakeParamsV1::v1()), "0" => r.p.dam_lake = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        // Lane EDIFICE-FOOT: dam_crater=1 (after dam_lake=1) crater lakes; ef_cap=1 (window-only foot drape cap),
        // ef_k_permille=, ef_fill_mm= (after ef_cap=1).
        if let Some(d) = s.strip_prefix("dam_crater=") {
            let Some(c) = r.p.dam_lake.take() else { die("dam_crater needs dam_lake=1 given before it".into()) };
            match d { "1" => r.p.dam_lake = Some(c.with_crater()), "0" => r.p.dam_lake = Some(pf_region::t1b_dam::DamLakeParamsV1 { crater: 0, crater_rule: None, ..c }), _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some(d) = s.strip_prefix("ef_cap=") {
            match d { "1" => r.p.edifice_foot = Some(pf_region::t1b_ef::EdificeFootParamsV1::v1()), "0" => r.p.edifice_foot = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some((k, d)) = s.split_once('=').filter(|(k, _)| *k == "ef_k_permille" || *k == "ef_fill_mm") {
            let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            let Some(c) = r.p.edifice_foot.as_mut() else { die(format!("{k} needs ef_cap=1 given before it")) };
            if k == "ef_k_permille" { c.k_permille = v; } else { c.fill_mm = v; }
            continue;
        }
        if let Some(d) = s.strip_prefix("dam_min_depth_mm=") {
            let v: i64 = d.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            let Some(c) = r.p.dam_lake.as_mut() else { die("dam_min_depth_mm needs dam_lake=1 given before it".into()) };
            c.min_depth_mm = v;
            continue;
        }
        // Lane SEAM-CERT: seam_cert=1 (window-only pinned-only fills, evaluation-lattice routing, certificate).
        if let Some(d) = s.strip_prefix("seam_cert=") {
            match d { "1" => r.p.seam_cert = Some(pf_region::t1b_seam::SeamCertParamsV1::v1()), "0" => r.p.seam_cert = None, _ => die(format!("bad value {s}")) }
            continue;
        }
        if let Some(d) = s.strip_prefix("climate_floor_permille=") { climate_floor = Some(d.parse().unwrap_or_else(|_| die(format!("bad value {s}")))); continue; }
        // Lane TD: diagnostic crop of per-stage arrays (not published, not digested).
        if let Some(d) = s.strip_prefix("diag_dump=") {
            pf_region::t1b_rivernet::DIAG.store(true, std::sync::atomic::Ordering::Relaxed);
            *pf_region::t1b_rivernet::DIAG_DIR.lock().unwrap() = Some(PathBuf::from(d));
            continue;
        }
        if s.starts_with("ct_escarp=") || s.starts_with("ct_contrast=") || s.starts_with("ct_contrast_n=") { continue; }
        if let Some((k, v)) = s.split_once('=') {
            if matches!(k, "hs" | "hs_n" | "hs_dt_yr" | "hs_d" | "hs_k_permille" | "hs_a1_m2" | "hs_th" | "hs_p_permille" | "hs_talus_permille" | "hs_route_every" | "hs_uplift" | "hs_runout" | "hs_volcano") {
                hs_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}")))));
                continue;
            }
            if matches!(k, "vf_floor" | "vf_foot" | "vf_fan" | "vf_window") {
                vf_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}")))));
                continue;
            }
            if matches!(k, "cf" | "cf_cone_hw" | "cf_clear_foot" | "cf_flows" | "cf_gullies" | "cf_fan_nocut" | "cf_shelf" | "cf_damp" | "cf_cliff" | "cf_softk_mm" | "cf_taper" | "cf_bury_inner" | "cf_bury_outer") {
                cf_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}")))));
                continue;
            }
            if matches!(k, "fx" | "fx_lh" | "fx_capture" | "fx_capture_passes" | "fx_kappa_permille" | "fx_floor" | "fx_contact" | "fx_bluff" | "fx_publish" | "fx_steep" | "fx_gamma_permille" | "fx_t1a_keep") {
                fx_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}")))));
                continue;
            }
            // Lane GUIDE-SPIKES: ls_undercut=1 (guide_undercut_v1; needs ls=1 ls_steep=2), applied with the ls_* overrides.
            if k == "ls_undercut" { ls_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}"))))); continue; }
            if matches!(k, "ls" | "ls_solve" | "ls_passes" | "ls_stop_ppm" | "ls_bank" | "ls_e" | "ls_hyst_mm" | "ls_t1a" | "ls_steep" | "ls_r4" | "ls_floor" | "ls_clip" | "ls_r4n" | "ls_pert_ppm" | "ls_pert_mode" | "ls_freeze" | "ls_curve" | "ls_cell" | "ls_cell_sweeps" | "ls_cell_margin_mm" | "tf" | "tf_report" | "tf_bluff" | "tf_subpost" | "tf_expo" | "tf_fetch_dirs" | "tf_fetch_mm" | "tf_cone" | "tg" | "tg_shore" | "tg_shore_cap_mm" | "tg_sill" | "tg_smooth" | "tg_fdb" | "tg_fdb_permille" | "th" | "th_expo" | "th_expo_len_mm" | "th_uncut" | "ve" | "ve_tmap_mm" | "ve_max_units" | "ve2" | "ve2_steer" | "ve2_erode" | "ve2_dome" | "ve2_taper" | "ve3" | "ve3_rough" | "ve3_pond" | "ve3_aa" | "ve4" | "ve4_lobes" | "ve4_lahar_pond" | "ve4_front" | "ve5" | "ve5_exposure" | "ve5_summit" | "ct" | "ct_step_mm" | "ct_crest" | "ct_apron" | "ct_rockmass" | "ct_apply" | "ct_window" | "ct_cap_permille" | "ct_var" | "ct_base" | "ct_local" | "ct_tread" | "ct_tread_deg" | "ct_tread_eta" | "ct_tread_full_deg" | "ct_tread_dip" | "ct_bed_median" | "ct_bed_sigma" | "ct_bed_min" | "ct_bed_max" | "ct_bed_eta_lo" | "ct_bed_eta_hi" | "ct_ends") {
                ls_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}")))));
                continue;
            }
            if matches!(k, "coast" | "volcano" | "shelf_mm" | "taper_mm" | "taper_min_permille" | "flood_init_mm" | "flood_max_mm" | "flood_gain_mm" | "clear_mm" | "fan_max_mm" | "cliff_retreat_mm") {
                co_over.push((k.to_string(), v.parse().unwrap_or_else(|_| die(format!("bad value {s}")))));
                continue;
            }
        }
        if let Some(b) = s.strip_prefix("box=") {
            let v: Vec<f64> = b.split(',').map(|x| x.parse().unwrap_or_else(|_| die(format!("bad box {s}")))).collect();
            if v.len() != 4 { die(format!("bad box {s}")); }
            r.boxes.push([v[0], v[1], v[2], v[3]]);
            continue;
        }
        if let Some((k, v)) = s.split_once('=') {
            let iv: i64 = v.parse().unwrap_or_else(|_| die(format!("bad value {s}")));
            let p = &mut r.p;
            match k {
                "realization" => p.realization = iv as u32, "net_mode" => p.network_mode = iv as u32, "f6_ref" => p.f6_reference = iv as u32, "rockmass_sigma_permille" => p.rockmass_sigma_permille = iv, "jitter_permille" => p.jitter_permille = iv, "meander_permille" => p.meander_permille = iv,
                "chaikin" => p.chaikin_passes = iv as u32, "a_ch" => p.channel_head_area_m2 = iv, "bed_cap_permille" => p.bed_slope_cap_permille_of_sc = iv,
                "floor_coeff_ppm" => p.floor_width_coeff_ppm = iv, "floor_max_mm" => p.floor_width_max_mm = iv, "floor_slope_ppm" => p.floor_slope_ppm = iv,
                "lref_mm" => p.hill_ref_length_mm = iv, "smin_ppm" => p.hill_slope_min_ppm = iv, "crest_max_mm" => p.crest_radius_max_mm = iv,
                "e0_gain_permille" => p.e0_gain_permille = iv, "solve_spacing_mm" => p.solve_spacing_mm = iv, "solve_tol_mm" => p.solve_tolerance_mm = iv,
                "solve_updates" => p.solve_updates_max = iv as u32, "margin_param_mm" => p.margin_mm = iv,
                "halo_mm" => r.halo_mm = iv, "margin_mm" => r.margin_mm = Some(iv), "diag_dx_mm" => r.diag_dx_mm = Some(iv),
                _ => die(format!("unknown key {k}")),
            }
        } else { r.keys.push(parse_xy(s)); }
    }
    // Lane CLIMATE-RUNOFF: the climate key sets params.climate_runoff (the recipe digest changes only when on).
    if r.climate.is_some() {
        let mut cp = pf_region::climate_runoff::ClimateRunoffParamsV1::v1();
        if let Some(f) = climate_floor { cp.floor_permille = f; }
        r.p.climate_runoff = Some(cp);
    } else if climate_floor.is_some() { die("climate_floor_permille needs climate=<dir>".into()); }
    for (k, v) in &co_over {
        match k.as_str() {
            "coast" => r.p.coast = if *v == 1 { Some(pf_region::t1b_coast::CoastParamsV1::v1()) } else { None },
            "volcano" => r.p.volcano = if *v == 1 { Some(pf_region::t1b_coast::VolcanoParamsV1::v1()) } else { None },
            _ => {}
        }
    }
    for (k, v) in &co_over {
        let c = r.p.coast.as_mut();
        match (k.as_str(), c) {
            ("coast" | "volcano", _) => {}
            ("shelf_mm", Some(c)) => c.shelf_depth_mm = *v, ("taper_mm", Some(c)) => c.taper_height_mm = *v, ("taper_min_permille", Some(c)) => c.taper_min_permille = *v,
            ("flood_init_mm", Some(c)) => c.flood_init_mm = *v, ("flood_max_mm", Some(c)) => c.flood_max_mm = *v, ("flood_gain_mm", Some(c)) => c.flood_gain_mm = *v,
            ("fan_max_mm", Some(c)) => c.fan_radius_max_mm = *v, ("cliff_retreat_mm", Some(c)) => c.cliff_retreat_mm = *v,
            ("clear_mm", _) => match r.p.volcano.as_mut() { Some(vp) => vp.clear_t1a_mm = *v, None => die("clear_mm needs volcano=1".into()) },
            _ => die(format!("{k} needs coast=1")),
        }
    }
    for (k, v) in &cf_over {
        if k == "cf" { r.p.cf = if *v == 1 { Some(pf_region::t1b_coast::CfParamsV1::v1()) } else { None }; }
    }
    for (k, v) in &cf_over {
        if k == "cf" { continue; }
        let Some(c) = r.p.cf.as_mut() else { die(format!("{k} needs cf=1")) };
        match k.as_str() {
            "cf_cone_hw" => c.cone_hw_permille = *v, "cf_clear_foot" => c.clear_foot_permille = *v, "cf_flows" => c.flow_count = *v as u32,
            "cf_gullies" => c.gully_count = *v as u32, "cf_fan_nocut" => c.fan_nocut = *v as u32, "cf_shelf" => c.t1a_shelf_grade = *v as u32,
            "cf_damp" => c.solve_damp = *v as u32, "cf_cliff" => c.cliff_profile = *v as u32, "cf_softk_mm" => c.softmax_k_mm = *v,
            "cf_taper" => c.foot_taper_permille = *v, "cf_bury_inner" => c.bury_inner_permille = *v, "cf_bury_outer" => c.bury_outer_permille = *v,
            _ => die(format!("unknown key {k}")),
        }
    }
    for (k, v) in &fx_over {
        if k == "fx" { r.p.fx = if *v == 1 { Some(pf_region::t1b_fix::FxParamsV1::v1()) } else { None }; }
    }
    for (k, v) in &fx_over {
        if k == "fx" { continue; }
        let Some(f) = r.p.fx.as_mut() else { die(format!("{k} needs fx=1")) };
        match k.as_str() {
            "fx_lh" => f.lh_free_mean = *v as u32, "fx_capture" => f.capture = *v as u32, "fx_capture_passes" => f.capture_passes_max = *v as u32,
            "fx_kappa_permille" => f.capture_kappa_permille = *v, "fx_floor" => f.soft_floor = *v as u32, "fx_contact" => f.contact_step = *v as u32,
            "fx_bluff" => f.bluff = *v as u32, "fx_publish" => f.publish = *v as u32, "fx_steep" => f.steep_cap = *v as u32, "fx_gamma_permille" => f.steep_gamma_permille = *v, "fx_t1a_keep" => f.t1a_keep = *v as u32,
            _ => die(format!("unknown key {k}")),
        }
    }
    for (k, v) in &ls_over {
        if k == "ls" { r.p.ls = if *v == 1 { Some(pf_region::t1b_ls::LsParamsV1::v1()) } else { None }; }
    }
    for (k, v) in &ls_over {
        if k == "ls" { continue; }
        let Some(l) = r.p.ls.as_mut() else { die(format!("{k} needs ls=1")) };
        match k.as_str() {
            "ls_solve" => l.solve = *v as u32, "ls_passes" => l.passes_max = *v as u32, "ls_stop_ppm" => l.stop_ppm = *v,
            "ls_bank" => l.bank_ramp = *v as u32, "ls_e" => l.e_update = *v as u32, "ls_hyst_mm" => l.hysteresis_mm = *v,
            // Lane TB options (the sub-struct appears in the recipe only when one of them is given).
            "ls_t1a" => l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).t1a_levels = *v as u32,
            "ls_steep" => l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).steep_eik = *v as u32,
            "ls_undercut" => l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).set_undercut(*v as u32),
            "ls_r4" => l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).r4_plain_permille = *v,
            // Lane TC options (sub-struct params.ls.tb.tc appears only when one of them is given).
            "ls_floor" | "ls_clip" | "ls_r4n" | "ls_pert_ppm" | "ls_pert_mode" | "ls_freeze" | "ls_curve" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off);
                match k.as_str() { "ls_floor" => c.floor = *v as u32, "ls_clip" => c.clip = *v as u32, "ls_r4n" => c.r4n_permille = *v, "ls_pert_ppm" => c.pert_ppm = *v,
                    "ls_pert_mode" => c.pert_mode = *v as u32, "ls_freeze" => c.freeze_from = *v as u32, _ => c.curve = *v as u32 }
            }
            // Lane TE options (sub-struct params.ls.tb.tc.te appears only when one of them is given).
            "ls_cell" | "ls_cell_sweeps" | "ls_cell_margin_mm" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).te.get_or_insert_with(pf_region::t1b_ls::TeParamsV1::v1);
                match k.as_str() { "ls_cell" => c.cell_solve = *v as u32, "ls_cell_sweeps" => c.sweeps_max = *v as u32, _ => c.sub_margin_mm = *v }
            }
            // Lane TF options (sub-struct params.ls.tb.tc.tf appears only when one of them is given; tf=1 all on, tf=0 all off,
            // then single switches).
            "tf" | "tf_report" | "tf_bluff" | "tf_subpost" | "tf_expo" | "tf_fetch_dirs" | "tf_fetch_mm" | "tf_cone" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).tf.get_or_insert_with(pf_region::t1b_tf::TfParamsV1::off);
                match k.as_str() {
                    "tf" => { let f = if *v == 1 { pf_region::t1b_tf::TfParamsV1::v1() } else { pf_region::t1b_tf::TfParamsV1::off() }; c.report = f.report; c.bluff_v2 = f.bluff_v2; c.subpost = f.subpost; c.exposure = f.exposure; c.cone_c1 = f.cone_c1; }
                    "tf_report" => c.report = *v as u32, "tf_bluff" => c.bluff_v2 = *v as u32, "tf_subpost" => c.subpost = *v as u32, "tf_expo" => c.exposure = *v as u32,
                    "tf_fetch_dirs" => c.fetch_dirs = *v as u32, "tf_fetch_mm" => c.fetch_cap_mm = *v, _ => c.cone_c1 = *v as u32,
                }
            }
            // Lane TG options (sub-struct params.ls.tb.tc.tg appears only when one of them is given; tg=1 = shore + sill +
            // smooth, tg=0 all off, then single switches; tg_fdb (F deadband) is never part of tg=1).
            "tg" | "tg_shore" | "tg_shore_cap_mm" | "tg_sill" | "tg_smooth" | "tg_fdb" | "tg_fdb_permille" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).tg.get_or_insert_with(pf_region::t1b_tg::TgParamsV1::off);
                match k.as_str() {
                    "tg" => { let f = if *v == 1 { pf_region::t1b_tg::TgParamsV1::v1() } else { pf_region::t1b_tg::TgParamsV1::off() }; c.shore = f.shore; c.sill_grade = f.sill_grade; c.smooth_bil = f.smooth_bil; }
                    "tg_shore" => c.shore = *v as u32, "tg_shore_cap_mm" => c.shore_cap_mm = *v, "tg_sill" => c.sill_grade = *v as u32, "tg_smooth" => c.smooth_bil = *v as u32,
                    "tg_fdb" => c.f_deadband = *v as u32, _ => c.f_deadband_permille = *v,
                }
            }
            // Lane TH options (sub-struct params.ls.tb.tc.th appears only when one of them is given; th=1 = expo smoothing +
            // uncut_mm, th=0 all off, then single switches).
            "th" | "th_expo" | "th_expo_len_mm" | "th_uncut" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).th.get_or_insert_with(pf_region::t1b_th::ThParamsV1::off);
                match k.as_str() {
                    "th" => { let f = if *v == 1 { pf_region::t1b_th::ThParamsV1::v1() } else { pf_region::t1b_th::ThParamsV1::off() }; c.expo_smooth = f.expo_smooth; c.uncut = f.uncut; }
                    "th_expo" => c.expo_smooth = *v as u32, "th_expo_len_mm" => c.expo_len_mm = *v, _ => c.uncut = *v as u32,
                }
            }
            // Lane VE options (sub-struct params.ls.tb.tc.ve appears only when one of them is given; ve=1 on, ve=0 off).
            "ve" | "ve_tmap_mm" | "ve_max_units" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).ve.get_or_insert_with(pf_region::t1b_ve::VeParamsV1::v1);
                match k.as_str() { "ve" => c.on = *v as u32, "ve_tmap_mm" => c.t_map_mm = *v, _ => c.max_units = *v as u32 }
            }
            // Lane VE2 options (params.ls.tb.tc.ve.ve2 appears only when one of them is given; ve2=1 all on, ve2=0 all off).
            "ve2" | "ve2_steer" | "ve2_erode" | "ve2_dome" | "ve2_taper" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).ve.get_or_insert_with(pf_region::t1b_ve::VeParamsV1::v1).ve2.get_or_insert_with(pf_region::t1b_ve::Ve2ParamsV1::off);
                match k.as_str() {
                    "ve2" => { let f = if *v == 1 { pf_region::t1b_ve::Ve2ParamsV1::v1() } else { pf_region::t1b_ve::Ve2ParamsV1::off() }; c.steer = f.steer; c.erode = f.erode; c.dome = f.dome; c.taper = f.taper; }
                    "ve2_steer" => c.steer = *v as u32, "ve2_erode" => c.erode = *v as u32, "ve2_dome" => c.dome = *v as u32, _ => c.taper = *v as u32,
                }
            }
            // Lane VE3 options (params.ls.tb.tc.ve.ve3; ve3=1 all on, ve3=0 all off; needs ve2 steer).
            "ve3" | "ve3_rough" | "ve3_pond" | "ve3_aa" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).ve.get_or_insert_with(pf_region::t1b_ve::VeParamsV1::v1).ve3.get_or_insert_with(pf_region::t1b_ve::Ve3ParamsV1::off);
                match k.as_str() {
                    "ve3" => { let f = if *v == 1 { pf_region::t1b_ve::Ve3ParamsV1::v1() } else { pf_region::t1b_ve::Ve3ParamsV1::off() }; c.rough = f.rough; c.pond = f.pond; c.aa = f.aa; }
                    "ve3_rough" => c.rough = *v as u32, "ve3_pond" => c.pond = *v as u32, _ => c.aa = *v as u32,
                }
            }
            // Lane VE4 options (params.ls.tb.tc.ve.ve4; ve4=1 all on, ve4=0 all off; needs ve3 rough).
            "ve4" | "ve4_lobes" | "ve4_lahar_pond" | "ve4_front" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).ve.get_or_insert_with(pf_region::t1b_ve::VeParamsV1::v1).ve4.get_or_insert_with(pf_region::t1b_ve::Ve4ParamsV1::off);
                match k.as_str() {
                    "ve4" => { let f = if *v == 1 { pf_region::t1b_ve::Ve4ParamsV1::v1() } else { pf_region::t1b_ve::Ve4ParamsV1::off() }; c.lobes = f.lobes; c.lahar_pond = f.lahar_pond; c.front = f.front; }
                    "ve4_lobes" => c.lobes = *v as u32, "ve4_lahar_pond" => c.lahar_pond = *v as u32, _ => c.front = *v as u32,
                }
            }
            // Lane VE5 options (params.ls.tb.tc.ve.ve5; ve5=1 all on, ve5=0 all off; needs ve4 lobes).
            "ve5" | "ve5_exposure" | "ve5_summit" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).ve.get_or_insert_with(pf_region::t1b_ve::VeParamsV1::v1).ve5.get_or_insert_with(pf_region::t1b_ve::Ve5ParamsV1::off);
                match k.as_str() {
                    "ve5" => { let f = if *v == 1 { pf_region::t1b_ve::Ve5ParamsV1::v1() } else { pf_region::t1b_ve::Ve5ParamsV1::off() }; c.exposure = f.exposure; c.summit = f.summit; }
                    "ve5_exposure" => c.exposure = *v as u32, _ => c.summit = *v as u32,
                }
            }
            // Lane STEP options (params.ls.tb.tc.ct; ct=1 the cliff-and-talus rule with its v1 table; absent = byte-identical).
            "ct" | "ct_step_mm" | "ct_crest" | "ct_apron" | "ct_rockmass" | "ct_apply" | "ct_window" | "ct_cap_permille" | "ct_var" | "ct_base" | "ct_local" | "ct_tread" | "ct_tread_deg" | "ct_tread_eta" | "ct_tread_full_deg" | "ct_tread_dip" | "ct_bed_median" | "ct_bed_sigma" | "ct_bed_min" | "ct_bed_max" | "ct_bed_eta_lo" | "ct_bed_eta_hi" | "ct_ends" => {
                let c = l.tb.get_or_insert_with(pf_region::t1b_ls::TbParamsV1::off).tc.get_or_insert_with(pf_region::t1b_ls::TcParamsV1::off).ct.get_or_insert_with(pf_region::t1b_ct::CtParamsV1::v1);
                // Lane STEP2 sub-keys (params.ls.tb.tc.ct.v2; absent = lane STEP byte-identical): ct_apply=2, ct_var, ct_base.
                if (k == "ct_apply" && *v == 2) || k == "ct_var" || k == "ct_base" || k == "ct_local" { c.v2.get_or_insert_with(pf_region::t1b_ct::CtV2ParamsV1::v1); }
                if k == "ct_var" { let w = c.v2.as_mut().unwrap(); w.var = *v as u32; if *v == 3 { w.bed.get_or_insert_with(pf_region::t1b_ct::CtBedParamsV1::v1); } continue; }
                // Lane BENCH-VAR sub-keys (params.ls.tb.tc.ct.v2.bed; need ct_var=3, any order; absent = lane TREAD2 byte-identical).
                if k.starts_with("ct_bed_") {
                    let b = c.v2.get_or_insert_with(pf_region::t1b_ct::CtV2ParamsV1::v1).bed.get_or_insert_with(pf_region::t1b_ct::CtBedParamsV1::v1);
                    match k.as_str() { "ct_bed_median" => b.median_permille = *v, "ct_bed_sigma" => b.sigma_permille = *v, "ct_bed_min" => b.min_permille = *v, "ct_bed_max" => b.max_permille = *v, "ct_bed_eta_lo" => b.eta_lo_permille = *v, _ => b.eta_hi_permille = *v }
                    continue;
                }
                if k == "ct_base" { c.v2.as_mut().unwrap().base = *v as u32; continue; }
                // Lane PERIOD-ENDS (params.ls.tb.tc.ct.v2.ends; needs ct_apply=2; absent = byte-identical): publish period ends.
                if k == "ct_ends" {
                    if *v != 1 { die("ct_ends: only ct_ends=1 (omit the key for off)".into()); }
                    let w = c.v2.get_or_insert_with(pf_region::t1b_ct::CtV2ParamsV1::v1);
                    w.ends = 1; w.ends_rule = pf_region::t1b_ct::ENDS_RULE.into();
                    continue;
                }
                // Lane TREAD2 sub-keys (params.ls.tb.tc.ct.v2.tread; absent = lane STEP2 byte-identical).
                if k.starts_with("ct_tread") {
                    let w = c.v2.get_or_insert_with(pf_region::t1b_ct::CtV2ParamsV1::v1);
                    if k == "ct_tread" { if *v != 1 { die("ct_tread: only ct_tread=1 (omit the key for off)".into()); } w.tread.get_or_insert_with(pf_region::t1b_ct::CtTreadParamsV1::v1); continue; }
                    let Some(t) = w.tread.as_mut() else { die(format!("{k} needs ct_tread=1")) };
                    let tanp = |d: i64| ((d as f64).to_radians().tan() * 1000.0).round() as i64;
                    match k.as_str() { "ct_tread_deg" => t.tread_tan_permille = tanp(*v), "ct_tread_eta" => t.eta_permille = *v, "ct_tread_full_deg" => t.full_tan_permille = tanp(*v), _ => t.dip = *v as u32 }
                    continue;
                }
                if k == "ct_local" { let w = c.v2.as_mut().unwrap(); w.local = *v as u32; if w.local == 1 { w.local_rule = pf_region::t1b_ct::LOCAL_RULE.into(); } continue; }
                match k.as_str() {
                    "ct" => { if *v != 1 { die("ct: only ct=1 (omit the key for off)".into()); } }
                    "ct_step_mm" => c.step_mm = *v, "ct_crest" => c.crest_exempt = *v as u32, "ct_apron" => c.apron = *v as u32, "ct_apply" => c.apply = *v as u32, "ct_window" => c.window_only = *v as u32, "ct_cap_permille" => c.shift_cap_permille = *v, _ => c.rockmass_cliff = *v as u32,
                }
            }
            _ => die(format!("unknown key {k}")),
        }
    }
    if !vf_over.is_empty() {
        let mut vp = pf_region::t1b_vf::VfParamsV1::off();
        let mut foot_given = false;
        for (k, v) in &vf_over {
            match k.as_str() { "vf_floor" => vp.floor = *v as u32, "vf_foot" => { vp.foot = *v as u32; foot_given = true; } "vf_window" => vp.window_only = *v as u32, _ => vp.fan = *v as u32 }
        }
        if !foot_given && vp.floor == 1 { vp.foot = 1; }
        if vp.floor == 1 || vp.fan == 1 { r.p.vf = Some(vp); }
    }
    if !hs_over.is_empty() {
        let mut hp = pf_region::t1b_hs::HsParamsV1::v0();
        hp.on = 0;
        for (k, v) in &hs_over {
            match k.as_str() {
                "hs" => hp.on = *v as u32, "hs_n" => hp.steps = *v as u32, "hs_dt_yr" => hp.dt_yr = *v, "hs_d" => hp.d_m2_per_my = *v,
                "hs_k_permille" => hp.k_permille = *v, "hs_a1_m2" => hp.a1_m2 = *v, "hs_th" => hp.threshold_m_per_my = *v, "hs_p_permille" => hp.mfd_p_permille = *v,
                "hs_talus_permille" => hp.talus_tan_permille = *v, "hs_route_every" => hp.route_every = *v as u32, "hs_uplift" => hp.uplift = *v as u32,
                "hs_runout" => hp.runout_max_posts = *v as u32, _ => hp.volcano = *v as u32,
            }
        }
        if hp.on == 0 && hs_over.iter().any(|(k, _)| k != "hs" && k != "hs_volcano") { die("hs_* rate keys need hs=1".into()); }
        if hp.on == 1 || hp.volcano >= 1 { r.p.hs = Some(hp); }
    }
    // E5 uses the existing window-only hillslope stage, after all construction.
    let mut contrast_on = false;
    let mut escarp_mode = 0;
    let mut contrast_n = None;
    for s in &a[5..] {
        if let Some(v) = s.strip_prefix("ct_escarp=") { escarp_mode=match v {"0"=>0,"1"=>1,"2"=>2,"3"=>3,"4"=>4,"6"=>6,_=>die("ct_escarp expects 0|1|2|3|4|6".into())}; }
        if let Some(v) = s.strip_prefix("ct_contrast=") { contrast_on = match v { "0" => false, "1" => true, _ => die("ct_contrast expects 0|1".into()) }; }
        if let Some(v) = s.strip_prefix("ct_contrast_n=") { contrast_n = Some(v.parse::<u32>().unwrap_or_else(|_| die("bad ct_contrast_n".into()))); }
    }
    if escarp_mode != 0 && !contrast_on {die("ct_escarp requires ct_contrast=1".into());}
    if escarp_mode >= 2 {
        let ct = r.p.ls.as_mut().and_then(|l| l.tb.as_mut()).and_then(|t| t.tc.as_mut()).and_then(|t| t.ct.as_mut()).unwrap_or_else(|| die("ct_escarp=2|3|4|6 needs ct=1".into()));
        ct.escarp_reconstruct = Some(escarp_mode);
    }
    if contrast_n.is_some() && !contrast_on { die("ct_contrast_n needs ct_contrast=1".into()); }
    if contrast_on {
        if r.p.lithology.envelope.is_none() { die("ct_contrast=1 needs env=<source envelope>".into()); }
        if r.p.hs.as_ref().is_some_and(|h| h.runs()) { die("ct_contrast is a separate transient; do not combine hs=1".into()); }
        let mut c = pf_region::t1b_ct::ContrastParams::v1();
        if let Some(n) = contrast_n { c.steps = n; }
        if escarp_mode == 1 {c.subgrid_spacing_mm=Some(2000); c.rule=format!("ESCARP-1 carrier v1: sparse2m geometry and1000yr fine steps for the same interval; coarse removal restriction; requires the escarp-v1 reader. Inherited conditional law: {}",c.rule); }
        c.validate().unwrap_or_else(|e| die(e));
        let h = r.p.hs.get_or_insert_with(pf_region::t1b_hs::HsParamsV1::v0);
        h.on = 1; h.steps = c.steps; h.dt_yr = c.dt_yr;
        h.label = pf_region::t1b_ct::CONTRAST_LABEL.into();
        h.creep_rule = "No soil creep on weathering-limited bare caprock in the E5 transient".into();
        h.slide_rule = c.rule.clone(); h.incision_rule = c.rule.clone();
        h.fill_rule = "No fill; every lowering retains a descending neighbour, Jacobi updates".into();
        h.contrast = Some(c);
    }
    // Lane VH-I (window-only): ring_plain=1 (needs cf=1) ring-plain tail; ve_wind=1 (needs ve=1 and climate=<dir>) tephra on the T0 wind.
    for s in &a[5..] {
        if let Some(d) = s.strip_prefix("ring_plain=") { match (d, r.p.cf.as_mut()) { ("1", Some(c)) => c.ring_plain = Some(pf_region::t1b_coast::RingPlainV1::v1()), ("0", _) => {}, _ => die(format!("{s}: needs cf=1 (values 0|1)")) } }
        if let Some(d) = s.strip_prefix("ve_wind=") { if d == "1" { let w = pf_region::t1b_ve::VeWindParamsV1::load(r.climate.as_ref().unwrap_or_else(|| die("ve_wind=1 needs climate=<dir>".into()))).unwrap_or_else(|e| die(format!("ve_wind: {e}"))); match r.p.ls.as_mut().and_then(|l| l.tb.as_mut()).and_then(|t| t.tc.as_mut()).and_then(|t| t.ve.as_mut()) { Some(v) => v.wind = Some(w), None => die("ve_wind=1 needs ve=1".into()) } } else if d != "0" { die(format!("bad value {s}")) } }
    }
    if let Err(e) = r.p.validate() { die(format!("params: {e}")); }
    r
}

/// CX-EDIFICE-EREC: ve_erec=0|1 (layer AND windows), bound once before recipe creation.
fn ve_erec_bind(a: &[String], ra: &mut RnArgs, t0: &pf_region::t0::T0) {
    // Keep zero absent from the serialized recipe, so explicit zero has the old bytes.
    let mut flank = false;
    for s in a { if let Some(v) = s.strip_prefix("ve_flank=") { flank = match v { "0" => false, "1" => true, _ => die(format!("bad value {s}: expected 0|1")) }; } }
    if flank {
        if ra.p.volcano.is_none() || ra.p.cf.is_none() || pf_region::t1b_ve::ve(&ra.p).is_none() { die("ve_flank=1 needs volcano=1 cf=1 ve=1".into()); }
        ra.p.ls.as_mut().unwrap().tb.as_mut().unwrap().tc.as_mut().unwrap().ve.as_mut().unwrap().flank = Some(pf_region::t1b_ve::FlankParamsV1::v1());
        ra.p.validate().unwrap_or_else(|e| die(format!("ve_flank: {e}")));
    }
    let mut mode = "0";
    let mut reference = None;
    let mut cache = None;
    for s in a {
        if let Some(v) = s.strip_prefix("ve_flank2=") { mode=v; }
        if let Some(v) = s.strip_prefix("ve_flank2_ref=") { reference=Some(v); }
        if let Some(v) = s.strip_prefix("ve_flank2_cache=") { cache=Some(v); }
    }
    let mut on = false;
    for s in &a[5..] {
        if let Some(v) = s.strip_prefix("ve_erec=") { match v { "0" => on = false, "1" => on = true, _ => die(format!("bad value {s}: expected 0|1")) } }
    }
    if on {
        if ra.p.volcano.is_none() || pf_region::t1b_ve::ve(&ra.p).is_none() { die("ve_erec=1 needs volcano=1 and ve=1".into()); }
        let record = pf_region::t0_erate::ErecRecordV1::load(&PathBuf::from(&a[0]), t0, 200, 5).unwrap_or_else(|e| die(format!("ve_erec: {e}")));
        let v = ra.p.ls.as_mut().and_then(|l| l.tb.as_mut()).and_then(|t| t.tc.as_mut()).and_then(|t| t.ve.as_mut()).unwrap();
        v.erec = Some(pf_region::t1b_ve::VeErecParamsV1 { rule: pf_region::t1b_ve::EREC_AGE_RULE.into(), record });
        ra.p.validate().unwrap_or_else(|e| die(format!("ve_erec: {e}")));
    }
    if mode != "0" {
        if !flank { die("ve_flank2 needs ve_flank=1".into()); }
        let fp=pf_region::t1b_ve::Flank2Params::bind(mode,reference,t0,std::path::Path::new(&a[2]),ra.p.realization,&serde_json::to_value(&ra.p).unwrap(),cache).unwrap_or_else(|e|die(format!("ve_flank2: {e}")));
        ra.p.ls.as_mut().unwrap().tb.as_mut().unwrap().tc.as_mut().unwrap().ve.as_mut().unwrap().flank2=Some(fp);
    }
}

/// Lane GLACIAL-T1: bind the t0.glacial layer to T0 and set params.glacial (its digest is part of the recipe).
fn gl_bind(ra: &mut RnArgs, t0: &pf_region::t0::T0) -> Option<pf_region::t1b_glacial::T0Glacial> {
    if ra.gl_now && ra.glacial.is_none() { die("gl_now=1 needs glacial=<t0.glacial dir>".into()); } // Lane GLACIER-1
    if ra.gl_bal && !ra.gl_now { die("gl_bal=1 needs gl_now=1".into()); } // Lane GLACIER-2
    if ra.gl_geo && !ra.gl_bal { die("gl_geo=1 needs gl_bal=1".into()); } // Lane GLACIER-3
    let d = ra.glacial.as_ref()?;
    let g = pf_region::t1b_glacial::T0Glacial::load(d, t0).unwrap_or_else(|e| die(format!("glacial: {e}")));
    ra.p.glacial = Some(pf_region::t1b_glacial::GlacialParamsV1::for_layer(&g, ra.gl_now).and_then(|p| p.with_bal(ra.gl_bal)).and_then(|p| p.with_geo(ra.gl_geo)).unwrap_or_else(|e| die(format!("glacial: {e}")))); // Lane GLACIER-1: v1 + scenario + gl_now; Lane GLACIER-2: + gl_bal; Lane GLACIER-3: + gl_geo
    if let Err(e) = ra.p.validate() { die(format!("params: {e}")); }
    println!("lane GLACIAL-T1: t0.glacial {} (K_g/K_f {}, Hf - He {} m, phase f_P {:?}, char f_P {:.3}, max f_P {})", g.digest, g.ratio, g.full_m, g.phase_fp, g.char_fp, g.max_fp);
    Some(g)
}

/// Lane CO: the volcano records (when given) and the recipe upstream entry for them.
/// Explicit canonical native owner, shared by rnwin and the resident server.
/// Default/none retains the established product and cache identity.
fn surface_owner_mode(a:&[String])->u8 {
    let mut mode=0;
    for s in a {if let Some(v)=s.strip_prefix("surface_owner="){mode=match v {
        "none"=>0,"native"=>1,"channels"=>2,_=>die(format!("bad {s}: expected none|native|channels")),
    };}}
    mode
}

fn co_records(ra: &RnArgs, recipe: &mut pf_region::t1b_rivernet::RiverNetRecipeV1) -> Option<pf_region::t1b_coast::VolcanoRecords> {
    let v = ra.volcanoes.as_ref().map(|p| pf_region::t1b_coast::VolcanoRecords::load(p).unwrap_or_else(|e| die(format!("volcanoes: {e}"))));
    if ra.p.volcano.is_some() && v.is_none() { die("volcano=1 needs volcanoes=<volcanoes.json>".into()); }
    if let (Some(v), Some(_)) = (&v, &ra.p.volcano) { recipe.upstream.push(("t0.volcano_features".to_string(), v.digest.clone())); recipe.upstream.sort(); }
    v
}

fn cmd_rivernet(a: &[String]) {
    use pf_region::t1b_rivernet::*;
    if a.len() < 5 { die("rivernet: <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> [strata=] [box=] [k=v]".into()); }
    let t_all = Instant::now();
    let mut ra = rn_args(a);
    if pf_region::t1b_ct::ct(&ra.p).map(|c| c.window_only == 1).unwrap_or(false) { die("ct_window=1 is a window-only measurement switch: build the layer without ct".into()); }
    if ra.p.vf.as_ref().map(|v| v.window_only == 1).unwrap_or(false) { die("vf_window=1 is a window-only switch: build the layer without vf keys".into()); }
    if ra.p.hs.is_some() { die("hs / hs_* / hs_volcano are window-only keys (rnwin): build the layer without them".into()); }
    if ra.p.seam_cert.is_some() { die("seam_cert is a window-only key (rnwin): build the layer without it".into()); }
    if ra.p.dam_lake.is_some() { die("dam_lake / dam_* are window-only keys (rnwin): build the layer without them".into()); }
    if ra.p.edifice_foot.is_some() { die("ef_cap / ef_* are window-only keys (rnwin): build the layer without them".into()); }
    if ra.p.cl.is_some() { die("cl / cl_* are window-only keys (rnwin): build the layer without them".into()); }
    if ra.p.cf.as_ref().map(|c| c.ring_plain.is_some()).unwrap_or(false) || pf_region::t1b_ve::ve(&ra.p).map(|v| v.wind.is_some()).unwrap_or(false) { die("ring_plain / ve_wind are window-only keys (rnwin): build the layer without them".into()); }
    let (t0, extras, site, skel, load_s) = net_inputs(a);
    ve_erec_bind(a, &mut ra, &t0);
    let gl0 = gl_bind(&mut ra, &t0);
    let strata = ra.strata.as_ref().map(|d| StrataCols::load(d).unwrap_or_else(|e| die(format!("strata: {e}"))));
    let mut recipe = rivernet_recipe(&t0, &extras, &site, &skel, strata.as_ref(), &ra.p).unwrap_or_else(|e| die(e));
    let vrec = co_records(&ra, &mut recipe);
    let rd = recipe.digest().unwrap_or_else(|e| die(e));
    if pf_region::t1b_esc_local::reuse_unchanged(&ra.p, &rd, std::path::Path::new(&a[4])).unwrap_or_else(|e| die(e)) { println!("esc_ref: unchanged inputs, reference bytes reused"); return; }
    let drd = recipe.draw_digest().unwrap_or_else(|e| die(e));
    let chart = vz_tilekey::chart::Chart::new(recipe.chart.clone()).unwrap();
    let anchor = (site.site_in_chart_mm[0] as f64 / 1000.0, site.site_in_chart_mm[1] as f64 / 1000.0);
    let t = Instant::now();
    let net = build_network(&t0, &skel, &chart, site.seed, &drd, &ra.p).unwrap_or_else(|e| die(e));
    let net_s = t.elapsed().as_secs_f64();
    let co = CoInputs { volcanoes: if ra.p.volcano.is_some() { vrec.as_ref() } else { None }, seed: site.seed, recipe_digest: drd.clone(), realization: ra.p.realization, glacial: gl0.as_ref(), ped_ground: None };
    if drd != rd { println!("lane TB paired keying: draws under {drd} (recipe without params.ls.tb); recipe {rd}"); }
    if ra.p.coast.is_some() { println!("coast: conduits {} (nodes {})", net.conduit_paths.len(), net.conduit.iter().filter(|&&c| c).count()); }
    let nctx = node_ctx(&net, &chart, &t0, &extras, &ra.p, strata.as_ref(), anchor).unwrap_or_else(|e| die(e));
    println!("network: land nodes {} unreached by the confined pass {} ({:.1} s)", net.nodes_land, net.unreached_confined, net_s);
    let sol = solve_e(&t0, &extras, &skel, &chart, &net, &nctx, &ra.p, strata.as_ref(), anchor, &co, &mut |s| println!("  {s}")).unwrap_or_else(|e| die(e));
    let mut beds = node_beds_warm_f(&net, sol.net.as_ref(), &nctx, &t0, &extras.laws, &sol.table, &sol.flood, &ra.p, sol.frozen);
    // Lane GLACIAL-T1: the published glacier table is built from the published (fluvial) beds; compared with the last solve evaluation's.
    if ra.p.glacial.is_some() {
        let tb = sol.cx.gl.clone().unwrap_or_else(|| die("glacial: the solve returned no glacier table".into()));
        let bn = beds.net(&net);
        let land_rcv_diff = (0..bn.kind.len()).filter(|&n| bn.kind[n] == pf_region::t1b_rivernet::K_LAND && bn.recv[n] != tb.rcv[n] && !(bn.recv[n] == u32::MAX - 1 && tb.rcv[n] == u32::MAX)).count();
        if land_rcv_diff > 0 { die(format!("glacial: published network differs from the glacier table's in {land_rcv_diff} land receivers")); }
        let kind = bn.kind.clone();
        tb.apply_beds(&kind, &mut beds.bed);
        println!("lane GLACIAL-T1: published glacier table (land receivers equal to the published network): {}", tb.stats);
    }
    // Lane GLACIER-1 (gl_now=1): the icefield sheet on the published layer surface (windows' volcano sites).
    let sol = { let mut sol = sol; if ra.p.glacial.as_ref().map(|g| g.now.is_some()).unwrap_or(false) { let m = pf_region::t1b_rivernet::gl_now_sheet(&t0, &extras, &skel, &chart, &net, &nctx, &ra.p, strata.as_ref(), anchor, &co, &mut sol, gl0.as_ref().unwrap()).unwrap_or_else(|e| die(format!("glacial now: {e}"))); println!("lane GLACIER-1: icefield sheet: {m}"); } sol };
    if let Some(c) = &beds.capture { println!("capture (TD2, published network): {}", serde_json::to_string(c).unwrap()); }
    if let Some(c) = &beds.ls { println!("level solve (LS1, published network): {}", serde_json::to_string(c).unwrap()); }
    let mut boxes = ra.boxes.clone();
    boxes.insert(0, [-131072.0, -131072.0, 131072.0, 131072.0]);
    let stats: Vec<NetStats> = boxes.iter().map(|b| net_stats(&net, &beds, &t0, ReliefParamsV1::default().threshold_m3_s(), *b)).collect();
    let out = PathBuf::from(&a[4]);
    // Lane MS-CEIL (ceil only): the lane_ceil header record from one evaluation of the published layer.
    let ceil_hdr: Vec<(&str, serde_json::Value)> = if ra.p.ceil.is_some() { let v = ceil_record(&t0, &extras, &skel, &chart, &net, &nctx, &ra.p, strata.as_ref(), anchor, &co, &sol).unwrap_or_else(|e| die(format!("ceil: {e}"))); println!("lane MS-CEIL: ceiling posts {} (solve lattice), counted held at e_max {} at e_min {}, counted F6-literal failures {}", v["ceil_posts"], v["counted_held_at_e_max"], v["counted_held_at_e_min"], v["counted_f6_literal_failures"]); vec![("lane_ceil", v)] } else { Vec::new() };
    let layer = publish_layer_ext(&recipe, &rd, &net, &beds, &sol, &stats, &out, &ceil_hdr).unwrap_or_else(|e| die(e));
    println!("rivernet layer digest {} recipe {} solve converged {} evals {} history {:?} ({:.1} s)", layer.digest, rd, sol.converged, sol.evals,
        sol.history_max_abs_diff_m.iter().map(|v| v.round()).collect::<Vec<_>>(), sol.seconds);
    for r in &sol.rows { if r.solved && r.land_hat_fraction > 0.3 && r.hat_mean_bilinear_m.abs() > 0.0 {
        let near = boxes.len() > 1;
        let _ = near;
        println!("  cell {} land {:.2} E0 {:.0} E {:.0}{} z {:.0} bil {:.0} trunk {:.0} diff {:.0}", r.cell, r.land_hat_fraction, r.e0_m_per_my, r.e_m_per_my, if r.at_bound { " (BOUND)" } else { "" }, r.hat_mean_z_m, r.hat_mean_bilinear_m, r.hat_mean_trunk_m, r.diff_m);
    } }
    for s in &stats { println!("  stats {}", serde_json::to_string(s).unwrap()); }
    // Lane T0-IMPRINT-PROD outlet_checks=1: ADR 045 per-T0-outlet F3 of the published network (own file; rivernet.json unchanged).
    if ra.p.outlet_checks == 1 {
        let f3 = pf_region::t1b_outlet::layer_f3_outlet(beds.net(&net), &skel, &t0, &chart).unwrap_or_else(|e| die(e));
        std::fs::write(out.join("checks_outlet.json"), serde_json::to_vec_pretty(&serde_json::json!({"layer_digest": layer.digest, "f3_outlet": f3})).unwrap()).unwrap();
        println!("  lane T0-IMPRINT-PROD F3 per T0 outlet pass {} | cells {} | basins {}", f3.pass,
            f3.cells.iter().map(|c| format!("{}:{:.0}/{:.0}{}", c.cell, c.share_permille, c.strict_share_permille, if c.gated { "" } else { "(rim)" })).collect::<Vec<_>>().join(" "),
            f3.basins.iter().map(|b| format!("{}:{:+.0}", b.basin, b.rel_err_permille)).collect::<Vec<_>>().join(" "));
    }
    for r in &sol.rows { if r.coastal == Some(true) && r.solved {
        println!("  coastal cell {} F {:.1} m land {:.3} bilinear land {:.3} | plain mean z {:.0} ref {:.0} bil {:.0} posts {} | hat diff {:.0} E {:.0}", r.cell, r.flood_offset_m.unwrap_or(0.0), r.land_fraction.unwrap_or(0.0), r.land_fraction_bilinear.unwrap_or(0.0),
            r.cell_mean_z_m.unwrap_or(0.0), r.cell_mean_ref_m.unwrap_or(0.0), r.cell_mean_bilinear_m.unwrap_or(0.0), r.cell_land_posts.unwrap_or(0), r.diff_m, r.e_m_per_my);
    } }
    for v in &sol.cx.volcanoes { if v.vent_xy_m[0].abs() < 60_000.0 && v.vent_xy_m[1].abs() < 60_000.0 { println!("  volcano {}", serde_json::to_string(v).unwrap()); } }
    if ra.p.coast.is_some() { println!("  fans {} (R >= 1 km: {}); F4 history {:?}", sol.cx.fans.len(), sol.cx.fans.iter().filter(|f| f.radius_m >= 1000.0).count(), sol.history_max_abs_f4.iter().map(|v| (v * 1000.0).round() / 1000.0).collect::<Vec<_>>()); }
    let exec = serde_json::json!({"load_seconds": load_s, "network_seconds": net_s, "solve_seconds": sol.seconds, "total_seconds": t_all.elapsed().as_secs_f64(), "cpu_seconds": cpu_seconds(), "peak_rss_kib": peak_rss_kib(), "not_digested": "execution receipt"});
    std::fs::write(out.join("execution_rivernet.json"), serde_json::to_vec_pretty(&exec).unwrap()).unwrap();
    println!("{exec}");
}

fn cmd_rnwin(a: &[String]) { let mode=surface_owner_mode(a); cmd_rnwin_impl(a,mode>0,mode==2,false); }

fn cmd_rnwin_impl(a: &[String], native_owner:bool, native_channels:bool, material_only:bool) {
    use pf_region::t1b_rivernet::*;
    if a.len() < 7 { die("rnwin: <t0_run_dir> <subgrid> <site> <t1a_dir> <out_dir> net=<dir> [strata=] [halo_mm=] [margin_mm=] <rx,ry>...".into()); }
    let t_all = Instant::now();
    let mut ra = rn_args(a);
    let net_dir = ra.net.clone().unwrap_or_else(|| die("rnwin needs net=<dir>".into()));
    let (t0, extras, site, skel, load_s) = net_inputs(a);
    ve_erec_bind(a, &mut ra, &t0);
    let _gl0 = gl_bind(&mut ra, &t0);
    let strata = ra.strata.as_ref().map(|d| StrataCols::load(d).unwrap_or_else(|e| die(format!("strata: {e}"))));
    let mut recipe = rivernet_recipe(&t0, &extras, &site, &skel, strata.as_ref(), &ra.p).unwrap_or_else(|e| die(e));
    let vrec = co_records(&ra, &mut recipe);
    let rd = recipe.digest().unwrap_or_else(|e| die(e));
    let (ld, lrd, table, flood) = load_layer_co(&net_dir).unwrap_or_else(|e| die(e));
    // Lane T1B-PEDESTAL ped_base=1: windows read the layer's published local ground (never recompute it).
    let ped_ground = if pf_region::t1b_ped::local_base(&ra.p) {
        // Lane MS-CEIL-2 ceil_ped_ref=<dir>: G from the named (key-off) layer instead (needs ceil).
        if ra.ped_ref.is_some() && ra.p.ceil.is_none() { die("ceil_ped_ref needs ceil=".into()); }
        let ped_dir = ra.ped_ref.clone().unwrap_or_else(|| net_dir.clone());
        if ra.ped_ref.is_some() { println!("lane MS-CEIL-2: local ground read from the reference layer {}", ped_dir.display()); }
        let m = load_layer_ped(&ped_dir).unwrap_or_else(|e| die(e)).unwrap_or_else(|| die("ped_base=1: the layer has no lane_ped.ground (build the layer with ped_base=1)".into()));
        println!("lane T1B-PEDESTAL: local ground read from the layer: {:?}", m);
        Some(m)
    } else { None };
    // Lane STEP ct_window=1 (measurement): the layer was solved without ct; check it against these params minus ct.
    let ct_win = pf_region::t1b_ct::ct(&ra.p).map(|c| c.window_only == 1).unwrap_or(false);
    let vf_win = ra.p.vf.as_ref().map(|v| v.window_only == 1).unwrap_or(false);
    let hs_win = ra.p.hs.is_some();
    let seam_win = ra.p.seam_cert.is_some();
    let cl_win = ra.p.cl.is_some();
    let dam_win = ra.p.dam_lake.is_some();
    let ef_win = ra.p.edifice_foot.is_some();
    let vh_win = ra.p.cf.as_ref().map(|c| c.ring_plain.is_some()).unwrap_or(false) || pf_region::t1b_ve::ve(&ra.p).map(|v| v.wind.is_some()).unwrap_or(false);
    let rd_layer = if ct_win || vf_win || hs_win || seam_win || cl_win || dam_win || ef_win || vh_win {
        let mut p0 = ra.p.clone();
        // Lane VH-I: ring_plain and ve_wind are always window-only.
        if vh_win { if let Some(c) = p0.cf.as_mut() { c.ring_plain = None; } if let Some(v) = p0.ls.as_mut().and_then(|l| l.tb.as_mut()).and_then(|t| t.tc.as_mut()).and_then(|t| t.ve.as_mut()) { v.wind = None; } println!("lane VH-I: layer checked against these params without ring_plain / ve_wind"); }
        // Lane EDIFICE-FOOT: edifice_foot is always window-only.
        if ef_win { p0.edifice_foot = None; println!("lane EDIFICE-FOOT: layer checked against these params without edifice_foot"); }
        // Lane DAMMED-VALLEY: dam_lake is always window-only.
        if dam_win { p0.dam_lake = None; println!("lane DAMMED-VALLEY: layer checked against these params without dam_lake"); }
        // Lane CHANNEL-CENTRELINE: cl is always window-only.
        if cl_win { p0.cl = None; println!("lane CHANNEL-CENTRELINE: layer checked against these params without cl"); }
        // Lane SEAM-CERT: seam_cert is always window-only.
        if seam_win { p0.seam_cert = None; println!("lane SEAM-CERT: layer checked against these params without seam_cert"); }
        // Lane HS-PASS: hs is always window-only.
        if hs_win { p0.hs = None; println!("lane HS-PASS: layer checked against these params without hs"); }
        // Lane VALLEY-FLOOR-FIX vf_window=1: the layer was solved without vf.
        if vf_win { p0.vf = None; println!("lane VALLEY-FLOOR-FIX vf_window=1: layer checked against these params without vf"); }
        if ct_win { if let Some(tc) = p0.ls.as_mut().and_then(|l| l.tb.as_mut()).and_then(|t| t.tc.as_mut()) { tc.ct = None; } }
        // Lane SPLIT-T1B: vs_ct_bed / vs_ct_blend act in the ct remap only, so a ct_window=1 layer is checked without them.
        if ct_win { if let Some(v) = p0.lithology.volcanic_split.as_mut() { *v = v.without_ct_keys(); } }
        let mut r0 = rivernet_recipe(&t0, &extras, &site, &skel, strata.as_ref(), &p0).unwrap_or_else(|e| die(e));
        let _ = co_records(&ra, &mut r0);
        let d0 = r0.digest().unwrap_or_else(|e| die(e));
        println!("lane STEP ct_window=1: window recipe {rd} (with ct), layer checked against {d0} (without ct)");
        d0
    } else { rd.clone() };
    if lrd != rd_layer { die(format!("layer recipe {lrd} differs from these params' recipe {rd}; refusing")); }
    let chart = vz_tilekey::chart::Chart::new(recipe.chart.clone()).unwrap();
    let anchor = (site.site_in_chart_mm[0] as f64 / 1000.0, site.site_in_chart_mm[1] as f64 / 1000.0);
    let t = Instant::now();
    let drd = recipe.draw_digest().unwrap_or_else(|e| die(e));
    // Lane VH-I: ring_plain / ve_wind change the window edifice surface and the apron axes, not draws (key-on shares the key-off draws).
    let drd = if vh_win { let mut p1 = ra.p.clone(); if let Some(c) = p1.cf.as_mut() { c.ring_plain = None; } if let Some(v) = p1.ls.as_mut().and_then(|l| l.tb.as_mut()).and_then(|t| t.tc.as_mut()).and_then(|t| t.ve.as_mut()) { v.wind = None; } let mut r1 = rivernet_recipe(&t0, &extras, &site, &skel, strata.as_ref(), &p1).unwrap_or_else(|e| die(e)); let _ = co_records(&ra, &mut r1); r1.draw_digest().unwrap_or_else(|e| die(e)) } else { drd };
    let net = build_network(&t0, &skel, &chart, site.seed, &drd, &ra.p).unwrap_or_else(|e| die(e));
    let nctx = node_ctx(&net, &chart, &t0, &extras, &ra.p, strata.as_ref(), anchor).unwrap_or_else(|e| die(e));
    // Lane TD2: windows take the layer's published (captured) receivers; no capture in the window.
    let net = if ra.p.ls.as_ref().map(|l| l.solve == 1).unwrap_or(false) { pf_region::t1b_ls::with_published_network(&net, &net_dir, &ra.p).unwrap_or_else(|e| die(e)) }
        else if ra.p.rerouted() { pf_region::t1b_fix::with_published_receivers(&net, &net_dir, &ra.p).unwrap_or_else(|e| die(e)) } else { net };
    let mut beds = if ra.p.rerouted() { node_beds_core(&net, &nctx, &t0, &extras.laws, &table, &flood, &ra.p) } else { node_beds(&net, &nctx, &t0, &extras.laws, &table, &flood, &ra.p) };
    // Lane GLACIAL-T1: windows read the layer's glacier table (never recompute it) and apply its filled glacial beds.
    // MEASUREMENT ONLY (env VZ_GL_DIAG_NOCARVE=1, stdout-announced): the key-on layer's E table and network without the
    // glacial beds and carve, to isolate the landform change from the E re-solve. Such windows are not publishable.
    let gl_nocarve = std::env::var_os("VZ_GL_DIAG_NOCARVE").is_some();
    if gl_nocarve { println!("lane GLACIAL-T1 MEASUREMENT: VZ_GL_DIAG_NOCARVE set: glacial beds and carve suppressed (not a publishable window)"); }
    let gl_table = if ra.p.glacial.is_some() && !gl_nocarve {
        let tb = pf_region::t1b_glacial::GlacialTable::load(&net_dir).unwrap_or_else(|e| die(format!("glacial table: {e}")));
        tb.apply_beds(&net.kind, &mut beds.bed);
        println!("lane GLACIAL-T1: glacier table read from the layer: {} cirques", tb.cirques.len());
        Some(std::sync::Arc::new(tb))
    } else { None };
    pf_region::t1b_mq::attach_bank_rock(&mut beds, &net, &nctx, &skel, &ra.p); // lane MEANDER-Q (no-op without cl_mq=1)
    let (mut cv, mut cl_stats) = curves_with_stats(&net, &beds, &skel, &t0, &ra.p, None);
    // Lane CHANNEL-CENTRELINE cl_mn=1: the belt is confined by volcano edifices (sited on the unmeandered `site` curves, so
    // the vents do not depend on the bends): second pass with each vent's influence disc (+ the 600 m max floor half-width).
    if pf_region::t1b_cl::cl(&ra.p).map(|c| c.meander_on()).unwrap_or(false) && ra.p.volcano.is_some() {
        let co0 = CoInputs { volcanoes: vrec.as_ref(), seed: site.seed, recipe_digest: drd.clone(), realization: ra.p.realization, glacial: None, ped_ground: ped_ground.clone() };
        let c0 = co0.ctx(&ra.p, &t0, &chart, &cv, &table, &flood).unwrap_or_else(|e| die(e));
        let mi = ra.p.volcano.as_ref().unwrap().mask_inner_mm as f64 / 1000.0;
        let pad = ra.p.floor_width_max_mm as f64 / 2000.0;
        let discs: Vec<([f64; 2], f64)> = c0.volcanoes.iter().map(|v| (v.vent_xy_m, v.r_influence_m + pad)).collect();
        let conf = |q: [f64; 2]| -> f64 { if discs.iter().any(|(c, r)| (q[0] - c[0]).hypot(q[1] - c[1]) < *r) { mi } else { f64::INFINITY } };
        let (cv2, st2) = curves_with_stats(&net, &beds, &skel, &t0, &ra.p, Some(&conf));
        println!("lane CHANNEL-CENTRELINE: meander belt confined near {} volcano edifices (mask_inner {:.0} m)", discs.len(), mi);
        cv = cv2; cl_stats = st2;
    }
    if let Some(c) = &cl_stats { println!("lane CHANNEL-CENTRELINE: T1a links {} on curve {} junction {} non-monotone {} no segment {}; node shift mean {:.0} max {:.0} m; connector max {:.0} m; meander segments {}",
        c.t1a_links, c.links_on_curve, c.junction_links, c.nonmonotone_links_straight, c.links_without_segment, c.node_shift_mean_m, c.node_shift_max_m, c.connector_max_m, c.meander.len()); }
    // Lane CHANNEL-CENTRELINE: optional dump of the T1a link polylines for measurement (env var, never digested, no bytes).
    if let Ok(path) = std::env::var("VZ_CL_DUMP") {
        let v: Vec<serde_json::Value> = cv.iter().filter(|c| c.link != u32::MAX && (c.code == pf_region::t1b::CHANNEL) && c.label == c.link && skel.kind.get(c.link as usize).map(|&k| k == 1 || k == 3).unwrap_or(false))
            .map(|c| serde_json::json!({"link": c.link, "seg": skel.segment[c.link as usize], "pts": c.pts.iter().map(|q| [q[0], q[1], q[2]]).collect::<Vec<_>>(),
                "carrier": c.carrier.as_ref().map(|(p, _)| p.iter().map(|q| [q[0], q[1]]).collect::<Vec<_>>())})).collect();
        std::fs::write(&path, serde_json::to_vec(&v).unwrap()).unwrap_or_else(|e| die(format!("VZ_CL_DUMP: {e}")));
        println!("lane CHANNEL-CENTRELINE: dumped {} T1a channel links to {path} (measurement only)", v.len());
    }
    let co = CoInputs { volcanoes: if ra.p.volcano.is_some() { vrec.as_ref() } else { None }, seed: site.seed, recipe_digest: drd.clone(), realization: ra.p.realization, glacial: None, ped_ground: ped_ground.clone() };
    let mut coast = co.ctx(&ra.p, &t0, &chart, &cv, &table, &flood).unwrap_or_else(|e| die(e));
    coast.gl = gl_table;
    // Lane VE: units are built on the published window surface only (the layer's E / F solve sees the body, not the units).
    coast.ve_units = coast.ve.is_some();
    let net_s = t.elapsed().as_secs_f64();
    let margin = ra.margin_mm.unwrap_or(ra.p.margin_mm);
    // Lane TE (stdout only, no bytes): the TD5 bluff reads the sea distance up to 8 km; a strip post sees halo + margin - seam
    // half width of its own lattice on the far side, and the fills need room beyond that (lane TE measured seams).
    if ra.p.fx.as_ref().map(|f| f.bluff == 1).unwrap_or(false) && ra.halo_mm + margin - 512_000 < 9_000_000 {
        println!("lane TE warning: halo {} + margin {} - seam half width 512000 mm < 9000000 mm (8 km bluff EDT + 1 km): seams are not guaranteed exact (use margin_mm=8192000)", ra.halo_mm, margin);
    }
    let vdig = if ra.p.volcano.is_some() { vrec.as_ref().map(|v| v.digest.clone()) } else { None };
    let prov = if ra.prov {
        let st = strata.as_ref().unwrap_or_else(|| die("prov=1 needs strata=".into()));
        let sp = pf_region::strata_prov::StrataProv::load(ra.strata.as_ref().unwrap(), ra.prov_json.as_deref(), &st.off, &st.lith, &st.th).unwrap_or_else(|e| die(format!("prov: {e}")));
        let volc = ra.volcanoes.as_ref().map(|p| pf_region::strata_prov::volcano_records(p).unwrap_or_else(|e| die(format!("prov: {e}"))));
        println!("lane PV prov: source {} present t {} My", sp.source, sp.present_t_my);
        Some(RnProv { strata: sp, volc })
    } else { None };
    let ctx = RnContext { t0: &t0, extras: &extras, site: &site, skel: &skel, strata: strata.as_ref(), params: &ra.p, chart: chart.clone(), net, table, layer_digest: ld.clone(), rn_recipe_digest: rd.clone(), curves: cv, anchor, coast, volcano_digest: vdig, prov, cl_stats };
    pf_region::t1b_esc_deps::bind_layer(&net_dir, &ld, &ra.p).unwrap_or_else(|e| die(e));
    let out = PathBuf::from(&a[4]);
    std::fs::create_dir_all(&out).unwrap();
    let mut exec = serde_json::Map::new();
    exec.insert("load_seconds".into(), load_s.into());
    exec.insert("network_seconds".into(), net_s.into());
    let zc = std::cell::OnceCell::new();
    let owner = if native_owner { Some(pf_region::t1b_surface_stages::NativeOwner::new(&ctx,&std::env::var_os("VZ_SURFACE_CACHE").map(PathBuf::from).unwrap_or_else(||out.join("surface-cache"))).unwrap_or_else(|e|die(e))) } else {None};
    if material_only {
        let native=owner.as_ref().unwrap_or_else(||die("material source requires canonical owner".into()));
        for (tx,ty) in &ra.keys {
            let tile=pf_region::t1b_ve::flank3_material_tile(&ctx,native,*tx,*ty).unwrap_or_else(|e|die(format!("material {tx},{ty}: {e}")));
            let path=out.join(format!("t_{tx}_{ty}"));
            pf_region::t1b_ve::flank3_publish_material_tile(&tile,&path).unwrap_or_else(|e|die(e));
            println!("material {tx},{ty} {}",tile.header["digest"]);
        }
        return;
    }
    let area_owner=owner.as_ref().map(|o|pf_region::t1b_surface_area::NativeAreaOwner::new(o,&out.join("surface-area"),native_channels)).transpose().unwrap_or_else(|e|die(e));
    // Lane T0-IMPRINT-PROD outlet_checks=1: the layer network's per-T0-outlet F3, once (ADR 045).
    let outlet_f3 = (ra.p.outlet_checks == 1).then(|| pf_region::t1b_outlet::layer_f3_outlet(&ctx.net, &skel, &t0, &chart).unwrap_or_else(|e| die(e)));
    for (rx, ry) in ra.keys.clone() {
        let wr = pf_region::t1b_rivernet::run_window_with_contacts(&ctx, ra.halo_mm, margin, rx, ry, owner.as_ref().map(|o|o as &dyn pf_region::t1b_surface_stages::Conditioner), area_owner.as_ref(), ra.unit_contacts).unwrap_or_else(|e| die(format!("window {rx},{ry}: {e}")));
        let wdir = out.join(format!("w_{rx}_{ry}"));
        if !pf_region::t1b_ve::publish_flank2_surface(&ra.p, &wr, &wdir).unwrap_or_else(|e| die(e)) && !pf_region::t1b_esc_deps::publish_evidence(&wr, &wdir).unwrap_or_else(|e| die(e)) { t1b_window::publish(&wr, &wdir).unwrap_or_else(|e| die(e)); }
        if ra.halo_mm == 0 {
            println!("window {} digest {} (halo 0: published only; the falsifier checks need a halo of at least one post)", wr.raw.lat.key(), wr.raw.digest);
            if let Some(sc) = wr.prep.stats.net.as_ref().and_then(|v| v.get("seam_cert")) { println!("  seam cert pass {} fill-tainted published {} route-bad published free {} lake extent {} donor extent {} ring {}", sc["pass"], sc["fill_tainted_published_posts"], sc["route_bad_published_free_posts"], sc["fill_lake_extent_published_posts"], sc["route_up_extent_published_posts"], sc["ring_distance_posts"]); }
            continue;
        }
        let t = Instant::now();
        let zc = zc.get_or_init(|| t1b_checks::t0_zero_contour(&t0, &chart).unwrap_or_else(|e| die(e)));
        let p = window_params(&ra.p, ra.halo_mm);
        let ch = check_window(&wr.raw, &wr.prep, &wr.header.ledger, &t0, &skel, &p, &wr.header.recipe.chart, Some(zc)).unwrap_or_else(|e| die(e));
        let check_s = t.elapsed().as_secs_f64();
        // Lane SEAM-CERT: the window's certificate is appended to checks.json (key on only).
        let seam = wr.prep.stats.net.as_ref().and_then(|v| v.get("seam_cert")).cloned();
        // Lane T0-IMPRINT-PROD: the ADR 045 gates block (outlet_checks=1 only).
        let ob = outlet_f3.as_ref().map(|lf| { let wf = t1b_checks::window_f3_outlet(&wr.raw, &wr.prep, &t0, &skel).unwrap_or_else(|e| die(e)); pf_region::t1b_outlet::window_block(lf, &net_dir, &wr.raw, &ch, &wf, &t0).unwrap_or_else(|e| die(e)) });
        if let Some(b) = &ob { println!("  lane T0-IMPRINT-PROD gates {}", b["gates"]); }
        // Lane DRAIN-FIX WALL-CHECK (wall_checks=1 only): the one-post wall census of the window core (reported, not gated).
        let wall = (ra.p.wall_checks == 1).then(|| t1b_checks::wall_census(&wr.raw));
        if let Some(w) = &wall { println!("  lane DRAIN-FIX walls (X {:?} m): all {:?} channel {:?} (guide {:?} t1a {:?}) pits {:?}; max wall {:.0} m, channel {:.0} m", w.thresholds_m, w.walls_all, w.walls_channel, w.walls_guide, w.walls_t1a, w.pits_all, w.max_wall_all_m, w.max_wall_channel_m); }
        if seam.is_some() || ob.is_some() || wall.is_some() {
            let mut v = serde_json::to_value(&ch).unwrap();
            if let serde_json::Value::Object(m) = &mut v { if let Some(sc) = &seam { m.insert("seam_cert".into(), sc.clone()); } if let Some(b) = &ob { m.insert("outlet_checks".into(), b.clone()); } if let Some(w) = &wall { m.insert("wall_check".into(), serde_json::to_value(w).unwrap()); } }
            std::fs::write(wdir.join("checks.json"), serde_json::to_vec_pretty(&v).unwrap()).unwrap();
        } else {
            std::fs::write(wdir.join("checks.json"), serde_json::to_vec_pretty(&ch).unwrap()).unwrap();
        }
        println!("window {} digest {} recipe {} (layer {})", wr.raw.lat.key(), wr.raw.digest, wr.raw.recipe_digest, ld);
        if let Some(sc) = &seam {
            println!("  seam cert pass {} (fill {} route {} core {}): fill-tainted published {} core {} (taint min dist out {}), lake extent {} posts, route-bad published free {} core {}, donor extent pub {} core {} posts, ring {} posts",
                sc["pass"], sc["fill_pass"], sc["route_pass"], sc["core_pass"], sc["fill_tainted_published_posts"], sc["fill_tainted_core_posts"], sc["fill_taint_min_dist_out_posts"], sc["fill_lake_extent_published_posts"],
                sc["route_bad_published_free_posts"], sc["route_bad_core_free_posts"], sc["route_up_extent_published_posts"], sc["route_up_extent_core_posts"], sc["ring_distance_posts"]);
        }
        println!("  F1 strict {} adr {} (core channel posts {}, T1a nodes not fixed {}, leave {}, rise {}, wrong {})  F2 {} ({:.2})  F3 {}  F4 app {} c1 {} c2 {}  F5 pits {} (>1 m {}, max {:.2})  F6 {}  F7 n/a",
            ch.f1_pass, ch.f1_adr_reading_pass, ch.f1_core_channel_posts, ch.f1_t1a_nodes_not_fixed, ch.f1_paths_leaving_channel, ch.f1_bed_rises, ch.f1_paths_to_wrong_outlet,
            ch.f2_pass, ch.f2_max_free_discharge_m3_s, ch.f3_pass, ch.f4_applicable, ch.f4_clause1_pass, ch.f4_clause2_pass,
            ch.f5_core_pit_posts, ch.f5_core_pit_posts_deeper_1m, ch.f5_core_pit_max_depth_m, ch.f6_pass);
        for c in &ch.cells { println!("  cell {} {} posts {} mean {:.0} bilinear {:.0} diff {:.0} tol {:.0} F6 {} F3 {:.1}", c.cell, c.water, c.core_posts, c.t1b_mean_m, c.bilinear_mean_m, c.f6_bilinear_diff_m, c.f6_tolerance_m, c.f6_pass, c.f3_share_permille); }
        let s = &ch.stats;
        println!("  stats mean {:.0} max {:.0} slope mean {:.1} p50 {:.1} p90 {:.1} >30 {:.3} >45 {:.3} local relief {:.0} valley {:.4}", s.mean_m, s.max_m, s.slope_mean_deg, s.slope_p50_deg, s.slope_p90_deg, s.slope_frac_gt_30deg, s.slope_frac_gt_45deg, s.local_relief_mean_m, s.valley_post_fraction);
        println!("  prep {}", wr.header.prep.net.as_ref().map(|v| v.to_string()).unwrap_or_default());
        println!("  seconds: prepare {:.1} evaluate {:.1} checks {:.1}", wr.seconds_prepare, wr.seconds_relax, check_s);
        exec.insert(format!("window_{rx}_{ry}"), serde_json::json!({"prepare_seconds": wr.seconds_prepare, "evaluate_seconds": wr.seconds_relax, "checks_seconds": check_s}));
        if let Some(dd) = ra.diag_dx_mm {
            // Resolution diagnostic: the same core evaluated at another spacing (not published).
            let mp = margin / dd;
            let (c0x, c0y) = (rx * 16_384_000 / dd, ry * 16_384_000 / dd);
            let el = EvalLat { gx0: c0x - mp, gy0: c0y - mp, n: (16_384_000 / dd + 2 * mp + 1) as usize, dx_mm: dd };
            let lp = lat_prep(el, &ctx.chart, &t0, &extras, &ra.p, strata.as_ref(), anchor).unwrap_or_else(|e| die(e));
            let sd = evaluate(&lp, &ctx.net, &ctx.curves, &t0, &extras.laws, &ctx.table, &ra.p, &ctx.coast);
            let mut acc: std::collections::BTreeMap<u32, (f64, f64, u64)> = std::collections::BTreeMap::new();
            for i in 0..el.nodes() {
                let (gx, gy) = el.g(i);
                if gx < c0x || gy < c0y || gx > c0x + 16_384_000 / dd || gy > c0y + 16_384_000 / dd { continue; }
                let e = acc.entry(lp.cell[i]).or_default();
                e.0 += sd.z[i]; e.1 += lp.bil[i] as f64; e.2 += 1;
            }
            for (c, (sz, sb, k)) in acc { println!("  diag dx {} mm: cell {c} mean {:.0} bilinear {:.0} posts {k}", dd, sz / k as f64, sb / k as f64); }
        }
    }
    exec.insert("total_seconds".into(), t_all.elapsed().as_secs_f64().into());
    exec.insert("cpu_seconds".into(), cpu_seconds().into());
    exec.insert("peak_rss_kib".into(), peak_rss_kib().into());
    exec.insert("not_digested".into(), "execution receipt".into());
    std::fs::write(out.join(format!("execution_rnwin_{}.json", std::process::id())), serde_json::to_vec_pretty(&serde_json::Value::Object(exec.clone())).unwrap()).unwrap();
    println!("{}", serde_json::Value::Object(exec));
}
