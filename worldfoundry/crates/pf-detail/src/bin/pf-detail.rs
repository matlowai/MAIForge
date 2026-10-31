//! pf-detail: T2 walk tiles (lane T2, TER-DETAIL). Labelled pre-signoff, Preview, derived, never canon.
//!
//! Usage:
//!   pf-detail block <config.json> <out_dir> <ix0> <iy0> <nx> <ny> [order=rev] [key=value ...]
//!     computes lod-0 tiles ix0..ix0+nx-1, iy0..iy0+ny-1 (TileKeyV2 chart 0), writes per-tile fields and checks,
//!     plus a block mosaic (height_mm.i32, litho.u8, base_mm.i32) and block.json.
//!     keys: remap, relax, channels, mean, fill (0/1), closure_rms_cap_m, steps, k_sp, d_creep, remap_window_m, all_off=1
//!   pf-detail compare <dir_a> <dir_b>      byte-compare every tile field of two block outputs

use pf_detail::checks::check_tile;
use pf_detail::input::*;
use pf_detail::tile::{assemble, Ctx, TileOut};
use pf_detail::*;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn die(m: String) -> ! { eprintln!("{m}"); std::process::exit(1) }

fn bytes_i32(v: &[i32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn bytes_f32(v: &[f32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn bytes_u16(v: &[u16]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }
fn bytes_u32(v: &[u32]) -> Vec<u8> { v.iter().flat_map(|x| x.to_le_bytes()).collect() }

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(String::as_str) {
        Some("block") => block(&a[2..]),
        Some("compare") => compare(&a[2..]),
        Some("aux") => aux(&a[2..]),
        Some("probe") => probe(&a[2..]),
        Some("column") => column(&a[2..]),
        _ => die("usage: see src/bin/pf-detail.rs".into()),
    }
}

fn write_tile(dir: &Path, t: &TileOut, header: &serde_json::Value) -> serde_json::Value {
    std::fs::create_dir_all(dir).unwrap();
    let mut files: Vec<(&str, &str, Vec<u8>)> = vec![
        ("height_mm", "i32 mm, 129x129 row-major, row 0 south", bytes_i32(&t.height_mm)),
        ("base_mm", "i32 mm, B1 bicubic T1b base (the F6/F8 reference)", bytes_i32(&t.base_mm)),
        ("litho", "u8 LithologyCodeV1 exposed at the final surface (255 none)", t.litho.clone()),
        ("layer", "u8 index of the exposed layer in the T0 column (top = 0)", t.layer.clone()),
        ("fixed", "u8 0 free, 1 channel interior (carved), 2 channel bank/corridor, 3 fixed sea/lake/low", t.fixed.clone()),
        ("closure_mm", "i32 mm, B4 keyed closure applied before relaxation (0 on fixed posts)", bytes_i32(&t.closure_mm)),
        ("cut_minus_base_m", "f32 m, B3 remap displacement after the channel taper (diagnostic)", bytes_f32(&t.cut_minus_base_m)),
    ];
    if !t.bed.is_empty() {
        files.push(("bed", "u16 v2: index of the exposed bed (hard/soft pairs; strata column or edifice sequence)", bytes_u16(&t.bed)));
        files.push(("over", "u8 v2: 0 strata, 1 edifice (T1b over-code), 2 deposit (fan/beach)", t.over.clone()));
        files.push(("remap_mm", "i32 mm v2: bounded stratigraphic remap displacement applied after the relaxation", bytes_i32(&t.remap_mm)));
    }
    if !t.prov_id.is_empty() {
        files.push(("prov_id", pf_detail::strata_prov::ID_NOTE, bytes_u32(&t.prov_id)));
    }
    if !t.depth_in_unit_m.is_empty() {
        files.push(("depth_in_unit_m", "f32 m lane BS unit_depth=1: depth of the surface below the TOP contact of the exposed unit (strata layer k >= 1: depth - top of layer k, vertical, unwarped); NaN where the unit has no recorded top contact (layer 0: its top is the present land surface; volcanic edifice; deposit) or water", bytes_f32(&t.depth_in_unit_m)));
        files.push(("height_above_unit_base_m", "f32 m lane BS unit_depth=1: height of the surface above the BASE contact of the exposed unit (strata: base of layer k; edifice: above the T1b pre-construct surface base_mm, bicubic); NaN for deposits, water, or an edifice without base_mm", bytes_f32(&t.height_above_unit_base_m)));
    }
    if !t.bed_style_id.is_empty() {
        files.push(("bed_style_id", pf_detail::bedstyle::FACIES_NOTE, t.bed_style_id.clone()));
        files.push(("flow_unit", pf_detail::bedstyle::UNIT_NOTE, bytes_u16(&t.flow_unit)));
    }
    if !t.pyro_sub.is_empty() {
        files.push(("pyro_sub", pf_detail::pyro::PYRO_FIELD_NOTE, t.pyro_sub.clone()));
    }
    if !t.de_mm.is_empty() {
        files.push(("de_talus", "u8 lane DE de=1: 0 none, 1 talus mantle (slope below a bare hard face carrying its debris; talus transport), 2 talus deposit >= 0.25 m from the threshold-failure step (seam-blended weights)", t.de_talus.clone()));
        files.push(("de_mm", "i32 mm lane DE de=1: the differential-erosion change (layered hillslope equilibrium minus the relaxed surface) after the seam blend (diagnostic; already in height_mm)", bytes_i32(&t.de_mm)));
    }
    if !t.sharp_mm.is_empty() {
        files.push(("sharp_mm", pf_detail::sharp::SHARP_MM_NOTE, bytes_i32(&t.sharp_mm)));
        files.push(("ct_regime", if header.get("ct_sharp_tread_note").is_some() { pf_detail::sharp::SHARP_REGIME_TREAD_NOTE } else { pf_detail::sharp::SHARP_REGIME_NOTE }, t.ct_regime.clone()));
    }
    let mut fl = Vec::new();
    for (name, note, b) in files {
        let ext = if name == "litho" || name == "layer" || name == "fixed" || name == "over" || name == "bed_style_id" || name == "de_talus" || name == "pyro_sub" || name == "ct_regime" { "u8" } else if name == "cut_minus_base_m" || name == "depth_in_unit_m" || name == "height_above_unit_base_m" { "f32" } else if name == "bed" || name == "flow_unit" { "u16" } else if name == "prov_id" { "u32" } else { "i32" };
        let f = format!("{name}.{ext}");
        std::fs::write(dir.join(&f), &b).unwrap();
        fl.push(serde_json::json!({"name": name, "file": f, "note": note, "digest": blake(&b)}));
    }
    let doc = serde_json::json!({"header": header, "tile": format!("0.{}.{}", t.ix, t.iy), "fields": fl,
        "base_mean_m": t.base_mean_m, "mean_shift_c_m": t.mean_shift_c_m, "fill_raised_posts": t.fill_raised_posts, "fill_max_raise_mm": t.fill_max_raise_mm});
    let mut doc = doc;
    if header.get("bedding_note").is_some() { doc["mean_deficit_m"] = t.mean_deficit_m.into(); doc["shift_max_m"] = t.shift_max_m.into(); }
    if t.lip_posts > 0 { doc["lip_posts"] = t.lip_posts.into(); doc["lip_raise_mm"] = t.lip_raise_mm.into(); }
    if let Some(v) = t.levee_mean_m { doc["levee_mean_m"] = v.into(); }
    if let Some(v) = t.carve_mean_m { doc["carve_mean_m"] = v.into(); }
    if let Some(v) = t.b8_cap_infeasible { doc["b8_cap_infeasible"] = v.into(); }
    if let Some(d) = &t.de_stats { doc["de"] = d.clone(); }
    if let Some(q) = &t.pyro_stats { doc["pyro_sub_counts"] = q.clone(); }
    if let Some(q) = &t.sharp_stats { doc["ct_sharp"] = q.clone(); }
    if !t.prov_id.is_empty() { doc["provenance"] = serde_json::json!({"field": "prov_id", "rows": t.prov_rows}); }
    std::fs::write(dir.join("tile.json"), serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
    serde_json::json!({"relax_seconds_own": t.relax_seconds_own, "assemble_seconds": t.assemble_seconds})
}

fn block(a: &[String]) {
    if a.len() < 6 { die("block <config.json> <out_dir> <ix0> <iy0> <nx> <ny> [k=v]".into()); }
    let t_all = Instant::now();
    let cfg: RunConfig = serde_json::from_slice(&std::fs::read(&a[0]).unwrap_or_else(|e| die(e.to_string()))).unwrap_or_else(|e| die(e.to_string()));
    let out = PathBuf::from(&a[1]);
    let p_i = |s: &str| s.parse::<i64>().unwrap_or_else(|e| die(format!("{s}: {e}")));
    let (ix0, iy0, nx, ny) = (p_i(&a[2]), p_i(&a[3]), p_i(&a[4]), p_i(&a[5]));
    let mut p = DetailParamsV1::v1();
    let mut rev = false;
    // lane PV: prov=1 (on), prov_strata=<StratigraphyDomainV0 JSON> (when the T0 dir has no env/age dumps), prov_volcanoes=<T0 volcanoes.json>
    let (mut prov_on, mut prov_json, mut prov_volc) = (false, None::<PathBuf>, None::<PathBuf>);
    for kv in &a[6..] {
        let (k, v) = kv.split_once('=').unwrap_or_else(|| die(format!("bad arg {kv}")));
        if k == "order" { rev = v == "rev"; continue; }
        if k == "prov" { prov_on = v == "1"; continue; }
        if k == "prov_strata" { prov_json = Some(PathBuf::from(v)); continue; }
        if k == "prov_volcanoes" { prov_volc = Some(PathBuf::from(v)); continue; }
        if k == "version" { if v == "2" { p = DetailParamsV1::v2(); } else if v == "3" { p = DetailParamsV1::v3(); } else if v != "1" { die(format!("version {v}")); } continue; }
        p.set(k, v).unwrap_or_else(|e| die(e));
    }
    let t_load = Instant::now();
    let w = T1bWindow::load(Path::new(&cfg.t1b_window_dir)).unwrap_or_else(|e| die(e));
    let t1a = T1a::load(Path::new(&cfg.t1a_dir)).unwrap_or_else(|e| die(e));
    let t0 = T0Strata::load(Path::new(&cfg.t0_run_root)).unwrap_or_else(|e| die(e));
    let (anchor, site_digest) = site_anchor_m(Path::new(&cfg.site_record)).unwrap_or_else(|e| die(e));
    if prov_on { p.prov = Some(ProvParamsV1::v1()); }
    // lane PV records; lane BS (beds=2) reads the same records (loaded once)
    let beds_on = p.beds.is_some();
    let rec = if prov_on || beds_on {
        let sp = pf_detail::strata_prov::StrataProv::load(&Path::new(&cfg.t0_run_root).join("final"), prov_json.as_deref(), &t0.off, &t0.lith, &t0.th).unwrap_or_else(|e| die(format!("prov: {e}")));
        let volc = prov_volc.as_ref().map(|q| pf_detail::strata_prov::volcano_records(q).unwrap_or_else(|e| die(format!("prov: {e}"))));
        println!("lane PV/BS records: source {} present t {} My; T1b over {} volcano_id {} volcanoes {}", sp.source, sp.present_t_my, w.over.is_some(), w.volcano_id.is_some(), w.volcanoes.len());
        Some((std::sync::Arc::new(sp), volc))
    } else { None };
    let prov = if prov_on { rec.as_ref().map(|(sp, volc)| pf_detail::tile::DetailProv { strata: sp.clone(), volc: volc.clone() }) } else { None };
    let bs_rec = if beds_on { rec.as_ref().map(|(sp, volc)| pf_detail::tile::BsRecords { strata: sp.clone(), volc: volc.as_ref().map(|v| v.0.clone()) }) } else { None };
    let load_s = t_load.elapsed().as_secs_f64();
    let params_json = serde_json::to_value(&p).unwrap();
    let params_digest = pf_contracts::artifact::digest_canonical(&params_json).unwrap();
    let mut upstream: Vec<(String, String)> = vec![
        ("t1b.relief".into(), w.digest.clone()), ("t1b.recipe".into(), w.recipe_digest.clone()),
        ("t1a.skeleton".into(), t1a.digest.clone()), ("site".into(), site_digest.clone()),
    ];
    upstream.extend(t0.digests.iter().cloned());
    if let (Some((sp, volc)), true) = (&rec, prov_on || beds_on) {
        let mut extra = sp.digests.clone();
        if let Some((_, d)) = volc { extra.push(("t0.volcano_features".into(), d.clone())); }
        for e in extra { if !upstream.contains(&e) { upstream.push(e); } }
    }
    upstream.sort();
    let recipe = serde_json::json!({
        "schema": "vz.t2-ground-recipe/v1 (lane T2 proposal; not a contract)",
        "layer": LAYER, "generator": {"name": GENERATOR_NAME, "version": p.generator_version(), "tier": "Preview", "reproducibility": "exact"},
        "params": params_json, "params_digest": params_digest, "upstream": upstream,
        "seed_scope": {"keyed_rng": "vz.keyed-rng/v1", "world_seed": cfg.world_seed, "layer": format!("{GENERATOR_NAME}/{}", p.generator_version()), "chart": cfg.chart,
            "tile": "lod 0 (ix, iy) of the post's owning tile", "address": "p/<px>/<py> local post", "purpose": PURPOSE_ROUGHNESS},
        "anchor_chart_m": [anchor.0, anchor.1],
        "label": LABEL, "simtime": SIMTIME_LABEL, "approximations": [CLOSURE_LABEL, "strata-profile-remap-v1 (slope proportional to rock critical slope; declared parametrization of differential erosion, PREREG B3)",
            "column anchored at the site point (70 km column has no sub-cell position)", "detachment-limited relaxation: eroded material leaves the tile (no fans or talus deposition)",
            "lithology table LithologyV1 standin_v1 copied from pf-region t1b_litho.rs"],
    });
    let mut recipe = recipe;
    if p.v2.is_some() {
        recipe["approximations"].as_array_mut().unwrap().extend([
            serde_json::json!("v2 channels_v2: curves traced on the T1b channel posts (codes 1, 5) and smoothed (binomial + Chaikin); valley reprojected onto the curves; banks graded at the alluvium / bed critical slope (PREREG T2B V1-V3)"),
            serde_json::json!("v2 over-codes: edifice and deposit posts recovered from T1b exposed_litho against our own cut (T1b does not publish its over field); one edifice bed sequence per chart (no volcano id in the window)"),
            serde_json::json!("v2 bedding: hard/soft bed pairs conserving the class rock-mass Sc (sub-resolution-closure-v1, amplitude bound min(0.04 R, 8 m), R = T1b local relief)"),
        ]);
    }
    if p.v3p().is_some() {
        recipe["approximations"].as_array_mut().unwrap().extend([
            serde_json::json!("v3 bedding packages and lateral pinch-and-swell: a STATISTICAL sub-resolution presentation of the modeled T0 stratigraphy (keyed draws; class thickness-weighted mean of 1/Sc conserved per layer; class contacts, dip and column order are the modeled data and are not moved). Not simulated deposition: these bands are not newly simulated depositional events (design rule 2026-10-15) (sub-resolution-closure-v1, PREREG T2C W1-W2)"),
            serde_json::json!("v3 soil mantle: the stratigraphic remap fades out where the base slope is well below the window critical slope (W3)"),
            serde_json::json!("v3 relaxation reads remap-window Sc and K; channel heads follow the area-slope law A_c = A_1 / S^2 (Montgomery and Dietrich 1988; A_1 = 300 m2) (W4, W8)"),
            serde_json::json!("v3 local mass balance: the relaxation's and the remap's change are re-centred by a 20 m Gaussian low-pass over free posts (eroded volume put back within about 40 m; not a sediment-transport law); the residual tile-mean shift is tapered off at channel corridors and weighted toward a footslope apron near channels (W6, W7)"),
        ]);
    }
    if beds_on {
        recipe["approximations"].as_array_mut().unwrap().extend([
            serde_json::json!(format!("lane BS (opt-in beds=2): bed styles from the record replace the v2 hard/soft pairs and v3 packages. {}", pf_detail::bedstyle::BED_STYLE_NOTE)),
            serde_json::json!("lane BS: beds thinner than rmin_m are merged into the more similar neighbour (thickness-weighted mean of 1/m; the 2 m grid shows bedsets, not thinner beds); style numbers are project assumptions from recalled textbook ranges (lane PV rock_look_map.json sources), unverified"),
            serde_json::json!("lane BS: edifice bedding per volcano (T1b volcano_id, or the nearest listed vent within r_influence_m when the window has none); interflow weathered-horizon probability from cone height / (last eruption - onset) of the T0 volcano record (unknown span: a declared default)"),
        ]);
    }
    if p.beds.as_ref().and_then(|b| b.csets.as_ref()).is_some() {
        recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!("lane BS2 (opt-in bs_csets=1): shelf limestone cycles grouped in keyed 10-30 m sets; a set's marliness moves the marly-base fraction and the clean-bed and marl strength together (statistical sub-resolution closure; amplitude by recorded env, contrast by recorded burial; layer mean 1/Sc conserved)"));
    }
    if let Some(d) = p.de.as_ref() {
        recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!(format!("{} Guard bound |eta| <= B tanh(eta/B), B = min({} R, {} m) (ADR 043 condition 7 style, declared). Talus {}; stream power by bed {}; DE in W6 {}.",
            pf_detail::de::DE_NOTE, d.amp_frac, d.amp_cap_m, if d.talus.is_some() { "on (reduced-order: talus properties below bare hard faces, no explicit debris volume ledger)" } else { "off" },
            if d.sp_area_m2.is_some() { "on" } else { "off" }, if d.mb == Some(true) { "included" } else { "excluded" })));
        if d.lobes.is_some() { recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!(pf_detail::lobes::LOBE_NOTE)); }
    }
    if let Some(py) = p.pyro.as_ref() {
        recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!(pf_detail::pyro::PYRO_SUB_NOTE));
        if py.hard.is_some() { recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!(pf_detail::pyro::PYRO_HARD_NOTE)); }
    }
    if p.unit_depth == Some(true) {
        recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!("lane BS (opt-in unit_depth=1): depth_in_unit_m / height_above_unit_base_m copy the cut's own layer geometry (no new process, no draws)"));
    }
    if prov.is_some() {
        recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!("lane PV (opt-in prov=1): prov_id names the exposed unit per post and tile.json provenance copies its upstream record (T0 strata layer, volcano record, or T1b present-day deposit); no new process, no draws; unrecorded facts are null"));
    }
    if p.ct_sharp.is_some() {
        recipe["approximations"].as_array_mut().unwrap().push(serde_json::json!(format!("{} T1b's per-post cliff tangent (rock-mass multiplier) and period boundaries are recovered from the published fields, not re-run; a construction, no rockfall process", pf_detail::sharp::SHARP_NOTE)));
    }
    let recipe_digest = pf_contracts::artifact::digest_canonical(&recipe).unwrap();
    let layer_digest = format!("{GENERATOR_NAME}/{}", p.generator_version());
    let mut tiles: Vec<(i64, i64)> = Vec::new();
    for iy in iy0..iy0 + ny { for ix in ix0..ix0 + nx { tiles.push((ix, iy)); } }
    if rev { tiles.reverse(); }
    let t_ctx = Instant::now();
    let mut ctx = Ctx::new_with(&w, &t1a, &t0, p.clone(), anchor, &layer_digest, cfg.world_seed, cfg.chart, &tiles, bs_rec.as_ref()).unwrap_or_else(|e| die(e));
    ctx.prov = prov;
    let ctx = ctx;
    let ctx_s = t_ctx.elapsed().as_secs_f64();
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(out.join("recipe.json"), serde_json::to_vec_pretty(&serde_json::json!({"recipe": recipe, "digest": recipe_digest})).unwrap()).unwrap();
    let mut header = serde_json::json!({"recipe_digest": recipe_digest, "label": LABEL});
    if p.v3p().is_some() {
        header["bedding_note"] = serde_json::json!(pf_detail::v3::BEDDING_NOTE);
    }
    if beds_on { header["bed_style_note"] = serde_json::json!(pf_detail::bedstyle::BED_STYLE_NOTE); }
    if p.ct_sharp.is_some() { header["ct_sharp_note"] = serde_json::json!(pf_detail::sharp::SHARP_NOTE); }
    // lane SHARP2: the 3-limb (tread) re-profile is on when the window carries pf-region ct_tread's fields
    if ctx.sharp.as_ref().map(|s| s.tread).unwrap_or(false) { header["ct_sharp_tread_note"] = serde_json::json!(pf_detail::sharp::SHARP_TREAD_NOTE); }
    if let Some(d) = p.de.as_ref() {
        header["de_note"] = serde_json::json!(pf_detail::de::DE_NOTE);
        if d.lobes.is_some() { header["de_lobes_note"] = serde_json::json!(pf_detail::lobes::LOBE_NOTE); }
    }
    if let Some(bp) = p.beds.as_ref() {
        if bp.bed_meta == Some(true) { header["bed_meta"] = serde_json::json!({"bed_style_id": pf_detail::bedstyle::FACIES_NOTE, "flow_unit": pf_detail::bedstyle::UNIT_NOTE}); }
        if let Some(py) = p.pyro.as_ref() {
            header["pyro_sub_note"] = serde_json::json!(pf_detail::pyro::PYRO_SUB_NOTE);
            header["pyro_sub"] = serde_json::json!({"file": "pyro_sub.u8", "note": pf_detail::pyro::PYRO_FIELD_NOTE, "names": pf_detail::pyro::SUB_NAMES});
            if py.hard.is_some() { header["pyro_sub_hard_note"] = serde_json::json!(pf_detail::pyro::PYRO_HARD_NOTE); }
        }
        if bp.ed_contacts == Some(true) { header["ed_contacts_note"] = serde_json::json!("lane BS2 ed_contacts=1: on volcanic edifice ground depth_in_unit_m / height_above_unit_base_m are measured to the top / base of the GENERATED flow unit or pyroclastic package containing the point (bed coordinates of the per-volcano beds=2 sequence; keyed closure, not recorded eruptions). Strata posts keep the recorded layer contacts."); }
        if bp.csets.is_some() { header["csets_note"] = serde_json::json!("lane BS2 bs_csets=1: shelf limestone cycles grouped in keyed sets (10-30 m) whose marliness moves marl fraction and clean-bed / marl strength together (statistical closure; set amplitude by recorded env, contrast by recorded burial)"); }
    }
    if let Some(pv) = &ctx.prov {
        let mut h = pv.strata.header();
        h["field"] = serde_json::json!("prov_id");
        h["volcano_records"] = match &pv.volc { Some((_, d)) => serde_json::json!(d), None => serde_json::Value::Null };
        h["t1b_window"] = serde_json::json!(w.digest);
        header["provenance"] = h;
    }
    let mut timing = serde_json::Map::new();
    let mut checks = Vec::new();
    let mut outs: Vec<TileOut> = Vec::new();
    for &(ix, iy) in &tiles {
        let tt = Instant::now();
        let t = assemble(&ctx, ix, iy);
        let tm = write_tile(&out.join(format!("tiles/t_{ix}_{iy}")), &t, &header);
        let c = check_tile(&t, 60.0);
        std::fs::write(out.join(format!("tiles/t_{ix}_{iy}/checks.json")), serde_json::to_vec_pretty(&c).unwrap()).unwrap();
        println!("tile 0.{ix}.{iy}: F6 {:+.4} m {} | slope mean {:.1} (base {:.1}) max bin {:.3} at {} deg (base {:.3}) peaks {:?} bimodal {} | F10 {:+.1} pts | pits {} fill {} (max {} mm) | F1 {} ({} posts, leave {}, pit {}) | F8c p99 {:.1} max {:.1} m | litho {:?} layers {} | {:.2} s",
            c.f6_diff_m, c.f6_pass, c.slope_mean_deg, c.base_slope_mean_deg, c.slope_max_bin_share, c.slope_max_bin_deg, c.base_slope_max_bin_share, c.slope_peaks_deg, c.slope_bimodal,
            c.f10_points_over_base, c.f5_interior_pits, t.fill_raised_posts, t.fill_max_raise_mm, c.f1_pass, c.f1_t1a_interior_posts, c.f1_paths_leaving_corridor, c.f1_paths_ending_in_pit,
            c.f8c_abs_diff_p99_m, c.f8c_abs_diff_max_m, c.litho_share.iter().map(|(k, v)| (*k, (v * 100.0).round() / 100.0)).collect::<Vec<_>>(), c.layers_exposed, tt.elapsed().as_secs_f64());
        let mut tmo = tm.as_object().unwrap().clone();
        tmo.insert("tile_total_seconds".into(), tt.elapsed().as_secs_f64().into());
        timing.insert(format!("0.{ix}.{iy}"), serde_json::Value::Object(tmo));
        checks.push(c);
        outs.push(t);
    }
    // block mosaic (edge-inclusive): (nx*128+1)^2
    let (mw, mh) = ((nx * 128 + 1) as usize, (ny * 128 + 1) as usize);
    let mut mz = vec![0i32; mw * mh];
    let mut mb = vec![0i32; mw * mh];
    let mut ml = vec![0u8; mw * mh];
    let mut mf = vec![0u8; mw * mh];
    let mut mp: Vec<u32> = if ctx.prov.is_some() { vec![u32::MAX; mw * mh] } else { Vec::new() };
    // lane PYRO2: block mosaic of pyro_sub (additive file, pyro_sub=1 only)
    let mut mps: Vec<u8> = if outs.iter().any(|t| !t.pyro_sub.is_empty()) { vec![pf_detail::pyro::SUB_NONE; mw * mh] } else { Vec::new() };
    // lane SHARP: block mosaic of the 2 m ct_regime (additive file, ct_sharp=1 only)
    let mut mrg: Vec<u8> = if outs.iter().any(|t| !t.ct_regime.is_empty()) { vec![0u8; mw * mh] } else { Vec::new() };
    let mut prov_rows: std::collections::BTreeMap<u64, serde_json::Value> = std::collections::BTreeMap::new();
    let mut edge_mismatch = 0u64;
    let mut written = vec![false; mw * mh];
    for t in &outs {
        let (bx, by) = (((t.ix - ix0) * 128) as usize, ((t.iy - iy0) * 128) as usize);
        for j in 0..TILE_POSTS { for i in 0..TILE_POSTS {
            let k = (by + j) * mw + bx + i;
            let s = j * TILE_POSTS + i;
            if written[k] && mz[k] != t.height_mm[s] { edge_mismatch += 1; }
            written[k] = true;
            mz[k] = t.height_mm[s]; mb[k] = t.base_mm[s]; ml[k] = t.litho[s]; mf[k] = t.fixed[s];
            if !mp.is_empty() { mp[k] = t.prov_id[s]; }
            if !mps.is_empty() { mps[k] = t.pyro_sub[s]; }
            if !mrg.is_empty() { mrg[k] = t.ct_regime[s]; }
        } }
        for r in &t.prov_rows { if let Some(id) = r.get("id").and_then(serde_json::Value::as_u64) { prov_rows.entry(id).or_insert_with(|| r.clone()); } }
    }
    if !mp.is_empty() { std::fs::write(out.join("block_prov_id.u32"), bytes_u32(&mp)).unwrap(); }
    if !mps.is_empty() { std::fs::write(out.join("block_pyro_sub.u8"), &mps).unwrap(); }
    if !mrg.is_empty() { std::fs::write(out.join("block_ct_regime.u8"), &mrg).unwrap(); }
    std::fs::write(out.join("block_height_mm.i32"), bytes_i32(&mz)).unwrap();
    std::fs::write(out.join("block_base_mm.i32"), bytes_i32(&mb)).unwrap();
    std::fs::write(out.join("block_litho.u8"), &ml).unwrap();
    std::fs::write(out.join("block_fixed.u8"), &mf).unwrap();
    let strata: Vec<serde_json::Value> = ctx.strata.cells.values().map(|c| serde_json::json!({"cell": c.cell, "top_m": c.top_m, "grad_raw": [c.g_raw.0, c.g_raw.1], "grad_clamped": [c.g.0, c.g.1],
        "layers": c.col.lith.iter().zip(c.col.base_depth.iter()).take(40).map(|(l, d)| (*l, *d)).collect::<Vec<_>>(), "remap_max_drift_m": c.max_drift_m})).collect();
    let segs = ctx.ch.segs.iter().filter(|s| s.code == 1).count();
    let v2info = ctx.v2.as_ref().map(|v2| serde_json::json!({"curves": v2.curves.stats, "aux": v2.aux.counts,
        "edifice_beds_top": v2.ed.col.bed_base.iter().zip(v2.ed.col.bed_sc.iter()).take(12).map(|(b, s)| (*b, *s)).collect::<Vec<_>>()}));
    let doc = serde_json::json!({
        "recipe_digest": recipe_digest, "tiles": tiles, "order": if rev { "reverse" } else { "forward" }, "block_size_posts": [mw, mh],
        "shared_edge_mismatch_posts": edge_mismatch, "strata_cells": strata, "v2": v2info, "channel_segments": ctx.ch.segs.len(), "t1a_segments": segs,
        "timing": {"load_seconds": load_s, "context_seconds": ctx_s, "total_seconds": t_all.elapsed().as_secs_f64(), "per_tile": timing},
        "checks": checks,
        "block_digests": {"height": blake(&bytes_i32(&mz)), "litho": blake(&ml)},
    });
    let mut doc = doc;
    if beds_on {
        let st: Vec<serde_json::Value> = ctx.strata.bs_styles.iter().map(|((c, k), s)| serde_json::json!({"cell": c, "layer": k, "style": s})).collect();
        let (ed, src) = match ctx.v2.as_ref().and_then(|v| v.ed_bs.as_ref()) { Some(e) => (serde_json::json!(e.info), serde_json::json!(e.source)), None => (serde_json::Value::Null, serde_json::Value::Null) };
        doc["bed_styles"] = serde_json::json!({"note": pf_detail::bedstyle::BED_STYLE_NOTE, "params": p.beds, "strata_layers": st, "edifices": ed, "edifice_volcano_source": src});
        if let Some(lw) = ctx.v2.as_ref().and_then(|v| v.lobes.as_ref()) { doc["de_lobes"] = serde_json::json!({"note": pf_detail::lobes::LOBE_NOTE, "volcanoes": lw.info}); }
        if let Some(pw) = ctx.v2.as_ref().and_then(|v| v.pyro.as_ref()) {
            // block code counts by source (summed tile counts; tile-edge posts are counted in each tile that holds them)
            let mut tot: std::collections::BTreeMap<String, std::collections::BTreeMap<String, u64>> = Default::default();
            for t in &outs { if let Some(q) = t.pyro_stats.as_ref().and_then(|v| v.as_object()) { for (src, cs) in q { for (c, n) in cs.as_object().unwrap() { *tot.entry(src.clone()).or_default().entry(c.clone()).or_default() += n.as_u64().unwrap_or(0); } } } }
            doc["pyro_sub"] = serde_json::json!({"note": pf_detail::pyro::PYRO_SUB_NOTE, "params": p.pyro, "records": pw.recs, "names": pf_detail::pyro::SUB_NAMES, "counts_by_source": tot,
                "block_field": {"file": "block_pyro_sub.u8", "dtype": "u8", "size": [mw, mh], "layout": "row-major, row 0 south, edge-inclusive (nx*128+1) x (ny*128+1)", "digest": blake(&mps), "note": pf_detail::pyro::PYRO_FIELD_NOTE}});
        }
    }
    if let Some(sf) = &ctx.sharp {
        let mut st = std::collections::BTreeMap::<&str, f64>::new();
        for t in &outs { if let Some(q) = t.sharp_stats.as_ref() { for k in ["touched_posts", "lowered_posts", "raised_posts", "refill_raised_posts"] { *st.entry(k).or_default() += q[k].as_f64().unwrap_or(0.0); } } }
        doc["ct_sharp"] = serde_json::json!({"note": pf_detail::sharp::SHARP_NOTE, "params": p.ct_sharp, "window_build": sf.stats, "tile_totals": st,
            "block_field": {"file": "block_ct_regime.u8", "dtype": "u8", "size": [mw, mh], "digest": blake(&mrg), "note": if sf.tread { pf_detail::sharp::SHARP_REGIME_TREAD_NOTE } else { pf_detail::sharp::SHARP_REGIME_NOTE }}});
    }
    if !mp.is_empty() {
        let mut h = header["provenance"].clone();
        h["rows"] = serde_json::Value::Array(prov_rows.into_values().collect());
        h["block_field"] = serde_json::json!({"file": "block_prov_id.u32", "dtype": "u32", "digest": blake(&bytes_u32(&mp)), "note": pf_detail::strata_prov::ID_NOTE});
        doc["provenance"] = h;
    }
    std::fs::write(out.join("block.json"), serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
    println!("block {}x{} tiles; shared-edge mismatches {}; height digest {}; total {:.1} s (load {:.1}, ctx {:.1})", nx, ny, edge_mismatch, blake(&bytes_i32(&mz)), t_all.elapsed().as_secs_f64(), load_s, ctx_s);
}

fn compare(a: &[String]) {
    let (da, db) = (PathBuf::from(&a[0]), PathBuf::from(&a[1]));
    let mut n = 0; let mut diff = 0;
    let ta = da.join("tiles");
    let mut names: Vec<_> = std::fs::read_dir(&ta).unwrap().map(|e| e.unwrap().file_name()).collect();
    names.sort();
    for t in names {
        let tdir = ta.join(&t);
        if !db.join("tiles").join(&t).exists() { continue; }
        let mut files: Vec<_> = std::fs::read_dir(&tdir).unwrap().map(|e| e.unwrap().file_name()).collect();
        files.sort();
        for f in files {
            let fs = f.to_string_lossy();
            if fs.ends_with(".json") { continue; }
            n += 1;
            let x = std::fs::read(tdir.join(&f)).unwrap();
            let y = std::fs::read(db.join("tiles").join(&t).join(&f)).unwrap_or_default();
            if x != y { diff += 1; println!("DIFF {}/{}", t.to_string_lossy(), fs); }
        }
    }
    println!("compared {n} tile fields: {diff} differ");
}

/// aux <config.json> <out_dir> <ix0> <iy0> <nx> <ny>: v2 curve and over-code diagnostics over a tile box (no tiles computed):
/// aux.json (curve stats, over-code counts vs the strata cut) and curves.json (smoothed polylines, chart m, with bed).
fn aux(a: &[String]) {
    let cfg: RunConfig = serde_json::from_slice(&std::fs::read(&a[0]).unwrap_or_else(|e| die(e.to_string()))).unwrap_or_else(|e| die(e.to_string()));
    let out = PathBuf::from(&a[1]);
    let p_i = |s: &str| s.parse::<i64>().unwrap_or_else(|e| die(format!("{s}: {e}")));
    let (ix0, iy0, nx, ny) = (p_i(&a[2]), p_i(&a[3]), p_i(&a[4]), p_i(&a[5]));
    let p = DetailParamsV1::v2();
    let w = T1bWindow::load(Path::new(&cfg.t1b_window_dir)).unwrap_or_else(|e| die(e));
    let t1a = T1a::load(Path::new(&cfg.t1a_dir)).unwrap_or_else(|e| die(e));
    let t0 = T0Strata::load(Path::new(&cfg.t0_run_root)).unwrap_or_else(|e| die(e));
    let (anchor, _) = site_anchor_m(Path::new(&cfg.site_record)).unwrap_or_else(|e| die(e));
    let mut tiles = Vec::new();
    for iy in iy0..iy0 + ny { for ix in ix0..ix0 + nx { tiles.push((ix, iy)); } }
    let t = Instant::now();
    let ctx = Ctx::new(&w, &t1a, &t0, p.clone(), anchor, &format!("{GENERATOR_NAME}/{}", p.generator_version()), cfg.world_seed, cfg.chart, &tiles).unwrap_or_else(|e| die(e));
    let v2 = ctx.v2.as_ref().unwrap();
    std::fs::create_dir_all(&out).unwrap();
    let pl: Vec<&Vec<[f64; 3]>> = v2.curves.polylines.iter().collect();
    std::fs::write(out.join("curves.json"), serde_json::to_vec(&pl).unwrap()).unwrap();
    std::fs::write(out.join("raw_lines.json"), serde_json::to_vec(&v2.curves.raw_lines).unwrap()).unwrap();
    let ax = &v2.aux;
    std::fs::write(out.join("aux_over.u8"), &ax.over).unwrap();
    std::fs::write(out.join("aux_zmod_mm.i32"), bytes_i32(&ax.zmod.iter().map(|&v| to_mm(v)).collect::<Vec<_>>())).unwrap();
    let doc = serde_json::json!({"curves": v2.curves.stats, "aux": ax.counts, "aux_lattice": {"gx0": ax.gx0, "gy0": ax.gy0, "nx": ax.nx, "ny": ax.ny},
        "seconds": t.elapsed().as_secs_f64()});
    std::fs::write(out.join("aux.json"), serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
    println!("{}", serde_json::to_string(&doc).unwrap());
}

/// probe <config.json> <x_m> <y_m> [k=v]: v2 internals at one chart position (diagnostic).
fn probe(a: &[String]) {
    let cfg: RunConfig = serde_json::from_slice(&std::fs::read(&a[0]).unwrap()).unwrap();
    let (x, y): (f64, f64) = (a[1].parse().unwrap(), a[2].parse().unwrap());
    let mut p = DetailParamsV1::v2();
    for kv in &a[3..] { let (k, v) = kv.split_once('=').unwrap(); p.set(k, v).unwrap_or_else(|e| die(e)); }
    let w = T1bWindow::load(Path::new(&cfg.t1b_window_dir)).unwrap_or_else(|e| die(e));
    let t1a = T1a::load(Path::new(&cfg.t1a_dir)).unwrap_or_else(|e| die(e));
    let t0 = T0Strata::load(Path::new(&cfg.t0_run_root)).unwrap_or_else(|e| die(e));
    let (anchor, _) = site_anchor_m(Path::new(&cfg.site_record)).unwrap_or_else(|e| die(e));
    let (ix, iy) = ((x / 256.0).floor() as i64, (y / 256.0).floor() as i64);
    let ctx = Ctx::new(&w, &t1a, &t0, p.clone(), anchor, &format!("{GENERATOR_NAME}/{}", p.generator_version()), cfg.world_seed, cfg.chart, &[(ix, iy)]).unwrap_or_else(|e| die(e));
    let v2 = ctx.v2.as_ref().unwrap();
    let ax = &v2.aux;
    let pv = ctx.post((x / 2.0).round() as i64, (y / 2.0).round() as i64);
    println!("base {:.2} cut(valley) {:.2} zmod-bicubic {:.2} G {:?} fixed_sea {} corridor {} interior {:?}", pv.base, pv.cut, ax.bicubic(&ax.zmod, x, y), ax.bilinear(&ax.g, x, y), pv.fixed_sea, pv.hit.corridor, pv.hit.interior);
    let mut n = 0;
    v2.curves.near(x, y, 96.0, &mut |q| { if n < 12 { println!("  seg r {:.1} bed {:.2} hw {:.2} depth {:.2}", q.r, q.bed, q.hw, q.depth); } n += 1; });
    for dy in -2i64..=2 { let gy = (y / 32.0).round() as i64 + dy; let mut line = String::new();
        for dx in -2i64..=2 { let gx = (x / 32.0).round() as i64 + dx; let i = ax.idx(gx, gy);
            line += &format!(" [{:7.1} {} g{:5.2} zm{:7.1}]", w.h(gx, gy), w.idx(gx, gy).map(|k| w.fixed[k]).unwrap_or(9), i.map(|k| ax.g[k]).unwrap_or(f64::NAN), i.map(|k| ax.zmod[k]).unwrap_or(f64::NAN)); }
        println!("{line}"); }
}

/// column <config.json> <cell> <out.json> [version=2|3] [beds=2] [prov_strata=..] [k=v ...]: lane BS diagnostic: the placed column
/// of one T0 cell exactly as a block run builds it (kept layers, bed base depths, Sc, K; up-beds), as JSON.
fn column(a: &[String]) {
    let cfg: RunConfig = serde_json::from_slice(&std::fs::read(&a[0]).unwrap_or_else(|e| die(e.to_string()))).unwrap_or_else(|e| die(e.to_string()));
    let cell: u32 = a[1].parse().unwrap_or_else(|e| die(format!("cell: {e}")));
    let mut p = DetailParamsV1::v1();
    let mut prov_json = None::<PathBuf>;
    for kv in &a[3..] {
        let (k, v) = kv.split_once('=').unwrap_or_else(|| die(format!("bad arg {kv}")));
        if k == "prov_strata" { prov_json = Some(PathBuf::from(v)); continue; }
        if k == "version" { if v == "2" { p = DetailParamsV1::v2(); } else if v == "3" { p = DetailParamsV1::v3(); } continue; }
        p.set(k, v).unwrap_or_else(|e| die(e));
    }
    let t0 = T0Strata::load(Path::new(&cfg.t0_run_root)).unwrap_or_else(|e| die(e));
    let (anchor, _) = site_anchor_m(Path::new(&cfg.site_record)).unwrap_or_else(|e| die(e));
    let scope = vz_keyed_rng::Scope::root(cfg.world_seed).layer(&format!("{GENERATOR_NAME}/{}", p.generator_version())).unwrap().chart(cfg.chart).unwrap();
    let mut m = pf_detail::strata::StrataModel::new(anchor);
    if p.beds.is_some() {
        m.bs = Some(std::sync::Arc::new(pf_detail::strata_prov::StrataProv::load(&Path::new(&cfg.t0_run_root).join("final"), prov_json.as_deref(), &t0.off, &t0.lith, &t0.th).unwrap_or_else(|e| die(e))));
    }
    m.ensure(&t0, cell, &p, &scope).unwrap_or_else(|e| die(e));
    let c = &m.cells[&cell].col;
    let doc = serde_json::json!({"cell": cell, "lith": c.lith, "layer_base_m": c.base_depth, "bed_base_m": c.bed_base, "bed_sc": c.bed_sc, "bed_k": c.bed_k,
        "n_down": c.n_down, "up_top_m": c.up_top, "bed_fac": c.bed_fac, "bed_unit": c.bed_unit, "unit_top": c.unit_top, "styles": m.bs_styles.iter().map(|((_, k), s)| (k.to_string(), *s)).collect::<std::collections::BTreeMap<_, _>>()});
    std::fs::write(&a[2], serde_json::to_vec(&doc).unwrap()).unwrap();
    println!("cell {cell}: {} layers, {} down beds, {} up beds", c.lith.len(), c.n_down, c.up_top.len());
}
