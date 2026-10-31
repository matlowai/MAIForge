//! Exposed lithology for T1b 0.6.0 / guides 0.4.0 (lane T3, TER-REGION v3,
//! 2026-10-13). LABEL: pre-signoff, fidelity tier Preview, derived, never canon.
//!
//! design 2026-10-13: "feels like the steep slope thing should be a function of
//! rock type ... some geology lends itself to cliffs and gorges some lends
//! itsself to slopes" and "the world needs stratified layers of rock type 70km
//! cells isnt going to work for that at all is it". A real per-parcel layer
//! stack (SEED-STRATIGRAPHY) is queued; until it exists this module is an
//! honest STAND-IN behind a small interface:
//!
//! - [`ExposedLithology`]: what rock is exposed at post i when the surface is
//!   at z, as a mixture of exposed-lithology CODES (a 32 m post can sit between
//!   T0 cells, so the basement is a land-weighted mixture, never a noise draw).
//! - [`LithologyV1`]: the declared table code -> (critical slope, K multiplier).
//! - [`CoverOverBasementV1`]: the stand-in. Sedimentary cover occupies the band
//!   between the T0 surface (land-weighted bilinear) and that surface minus the
//!   T0 sediment thickness (land-weighted bilinear). A post whose current
//!   surface is at or above the base of cover exposes cover; below it, basement
//!   (class from T0 rock_family; for a SedimentaryCover cell, from crust kind).
//!   Deep valleys expose basement, interfluves keep cover. No draws.
//!
//! A real stack plugs in by implementing `ExposedLithology` (surface cut through
//! the deformed stack) and reusing the same table.

use crate::t0::{CellWater, T0};
use crate::t1b::bilinear;
use crate::*;
use serde::Serialize;

/// Exposed-lithology codes (published as `exposed_litho` u8). The 12 classes of the design review
/// stratigraphy design (`unpublished-design-record` section 1, LithologyV1),
/// so a real stack (SEED-STRATIGRAPHY stage B) plugs into the same table.
pub const LITHO_CRYSTALLINE: u8 = 0;
pub const LITHO_GRANITE: u8 = 1;
pub const LITHO_METAMORPHIC: u8 = 2;
pub const LITHO_BASALT: u8 = 3;
pub const LITHO_CARBONATE: u8 = 4;
pub const LITHO_SHALE: u8 = 5;
pub const LITHO_SANDSTONE: u8 = 6;
pub const LITHO_CONGLOMERATE: u8 = 7;
pub const LITHO_EVAPORITE: u8 = 8;
pub const LITHO_COAL: u8 = 9;
pub const LITHO_UNCONSOLIDATED: u8 = 10;
pub const LITHO_OPHIOLITE: u8 = 11;
pub const LITHO_CODES: usize = 12;
/// The class the stand-in uses for undifferentiated sedimentary cover.
pub const LITHO_COVER: u8 = LITHO_SHALE;

#[derive(Clone, Debug, Serialize, serde::Deserialize, PartialEq, Eq)]
pub struct LithoClassV1 {
    pub code: u8,
    pub name: String,
    /// Critical (threshold) slope as a D8 pair gradient, permille.
    pub critical_slope_permille: i64,
    /// Stream-power erodibility multiplier on the T0 header K, permille.
    pub k_mult_permille: i64,
    /// E5 conditional bare-rock lowering (mm/kyr). Absent means unknown, not zero.
    /// Set only by ct_contrast; the scenario evidence is serialized in hs.contrast.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bare_weathering_mm_per_kyr: Option<i64>,
    pub source: String,
}

/// The declared code table plus the rules (integers and strings only: exact canonical JSON).
#[derive(Clone, Debug, Serialize, serde::Deserialize, PartialEq, Eq)]
pub struct LithologyV1 {
    pub law: String,
    pub classes: Vec<LithoClassV1>,
    pub basement_rule: String,
    pub cover_rule: String,
    pub cap_rule: String,
    pub erodibility_rule: String,
    /// Cover thinner than this is treated as absent [default 1 m].
    pub min_cover_thickness_mm: i64,
    /// Lane VOLCANIC-SPLIT (opt-in `vs=1`): lava vs volcaniclastic facies inside lithology-3 strata layers (t1b_volc);
    /// absent = unchanged (not serialized).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub volcanic_split: Option<crate::t1b_volc::VolcanicSplitV1>,
    /// Lane LITHO-CONTINUITY (opt-in `lcont=1`): planet-wide rock-body identity and continuous contacts across T0 cells
    /// (t1b_lcont); absent = unchanged (not serialized).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub continuity: Option<crate::t1b_lcont::LithoContinuityV1>,
    /// E3: source-record geometry supersedes private continuity/split geometry.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub envelope: Option<EnvelopeConfig>,
}

const SRC_RECALLED: &str = "project assumption: angle and K multiplier from the design review stratigraphy design table (2026-10-13, section 4, 'project assumptions until cited'); Sc = tan(angle) as a D8 pair gradient (lane T3). Background recalled, unverified: threshold hillslopes 30-40 deg modal in rapidly eroding ranges (Burbank et al. 1996 Nature 379:505; Montgomery and Brandon 2002 EPSL 201:481); rock-mass strength sets the limiting slope (Selby 1980 Z. Geomorph. 24:31; Schmidt and Montgomery 1995 Science 270:617); erodibility falls with rock tensile strength (Sklar and Dietrich 2001 Geology 29:1087) and spans orders of magnitude across lithologies (Stock and Montgomery 1999 JGR 104:4983); multipliers kept within 0.5-2 so T0's calibration still holds";

impl LithologyV1 {
    pub fn standin_v1() -> LithologyV1 {
        let c = |code: u8, name: &str, sc: i64, k: i64| LithoClassV1 { code, name: name.into(), critical_slope_permille: sc, k_mult_permille: k, bare_weathering_mm_per_kyr: None, source: SRC_RECALLED.into() };
        LithologyV1 {
            law: "threshold-slope-lithology-v0: stream power n = 1 with K_post = K_T0 x k_mult(exposed lithology); linear diffusion D_T0 below the critical slope; after diffusion each step, every free post is lowered so that z_i <= z_j + Sc_i x d_ij for all 8 neighbours j (threshold failure; exact min-plus closure by a Dijkstra pass from the lower neighbours of violating posts, key (z, index)); failed volume is exported along the flow path (detachment-limited, like stream power)".into(),
            classes: vec![
                c(LITHO_CRYSTALLINE, "crystalline basement (gneiss, undifferentiated); 38 deg", 781, 500),
                c(LITHO_GRANITE, "granite pluton; 45 deg", 1000, 400),
                c(LITHO_METAMORPHIC, "metamorphic (slate, schist, marble, quartzite); 35 deg (isotropic here: no fabric yet)", 700, 800),
                c(LITHO_BASALT, "basalt / volcanic (lava, arc volcanics, LIP), caprock; 60 deg", 1732, 600),
                c(LITHO_CARBONATE, "marine carbonate; 70 deg", 2747, 500),
                c(LITHO_SHALE, "marine shale / mudstone; 22 deg", 404, 2000),
                c(LITHO_SANDSTONE, "sandstone; 55 deg", 1428, 800),
                c(LITHO_CONGLOMERATE, "conglomerate; 45 deg", 1000, 900),
                c(LITHO_EVAPORITE, "evaporite; 25 deg", 466, 3000),
                c(LITHO_COAL, "coal-bearing shale/sandstone; 25 deg", 466, 1500),
                c(LITHO_UNCONSOLIDATED, "unconsolidated (alluvium, regolith); 32 deg", 625, 3000),
                c(LITHO_OPHIOLITE, "oceanic igneous / ophiolite; 40 deg", 839, 700),
            ],
            basement_rule: "cover_over_basement_v1 basement: per T0 cell, code from rock_family (0 oceanic_igneous -> 11 ophiolite, 1 continental_undifferentiated -> 0 crystalline basement, 2 arc_magmatic -> 1 granite pluton (the intrusive 0.7 of an arc; no stack says where the volcanic 0.3 lies), 3 collisional_thickened -> 2 metamorphic, 4 sedimentary_cover -> by crust kind: continental -> 0 crystalline, oceanic -> 11 ophiolite); at a post, critical slope and k_mult are the land-weighted bilinear mixture of the 4 surrounding T0 cells' class values (the post's own cell's class if none is land); the dominant code is the post's own T0 cell's".into(),
            cover_rule: "cover_over_basement_v1 cover: thickness = land-weighted bilinear T0 sediment_thickness_m; base of cover = land-weighted bilinear T0 elevation_m minus thickness; a post exposes cover (class 5 shale/mudstone: undifferentiated cover read as the commonest sedimentary rock, recalled) iff thickness >= min_cover_thickness and its current surface z >= base of cover, else basement (no draws; a real stratigraphy stack replaces this through the ExposedLithology interface)".into(),
            cap_rule: "the critical slope of the post being lowered (the upper post of the pair), as a D8 pair gradient (axis distance 32 m, diagonal 32 sqrt2 m); exposure evaluated on the surface after the source term, once per step".into(),
            erodibility_rule: "the stream-power K of post i is K_T0 x k_mult(exposure of i); the relaxation source term U keeps T0's uniform K (T0's closure rate), so hard rock stands higher and soft rock lower (differential erosion)".into(),
            min_cover_thickness_mm: 1_000,
            volcanic_split: None,
            continuity: None, envelope: None,
        }
    }
    pub fn class(&self, code: u8) -> &LithoClassV1 { &self.classes[code as usize] }
    pub fn sc(&self, code: u8) -> f64 { self.class(code).critical_slope_permille as f64 / 1000.0 }
    pub fn k(&self, code: u8) -> f64 { self.class(code).k_mult_permille as f64 / 1000.0 }
    pub fn validate(&self) -> R<()> {
        if self.classes.len() != LITHO_CODES { return Err("lithology table must list all 12 codes".into()); }
        for (i, c) in self.classes.iter().enumerate() {
            if c.bare_weathering_mm_per_kyr.is_some_and(|v| v < 0) { return Err("negative bare weathering rate".into()); }
            if c.code as usize != i { return Err("lithology table out of code order".into()); }
            if c.critical_slope_permille <= 0 || c.k_mult_permille <= 0 { return Err("lithology values must be positive".into()); }
        }
        if let Some(v) = &self.volcanic_split { v.validate()?; }
        if let Some(c) = &self.continuity { c.validate()?; }
        Ok(())
    }
}

/// Basement code of a T0 cell.
pub fn basement_code(rock_family: u8, crust_kind: u8) -> u8 {
    match rock_family {
        0 => LITHO_OPHIOLITE,
        1 => LITHO_CRYSTALLINE,
        2 => LITHO_GRANITE,
        3 => LITHO_METAMORPHIC,
        4 => if crust_kind == 0 { LITHO_OPHIOLITE } else { LITHO_CRYSTALLINE },
        _ => LITHO_CRYSTALLINE,
    }
}

/// What rock a post exposes at surface height z.
pub trait ExposedLithology {
    /// (critical slope gradient, K multiplier, dominant exposed code) of post i at surface z.
    fn at(&self, i: usize, z: f64) -> (f64, f64, u8);
}

/// The stand-in: sedimentary cover over a basement mixture (see module docs).
pub struct CoverOverBasementV1 {
    pub basement_sc: Vec<f64>,
    pub basement_k: Vec<f64>,
    pub basement_code: Vec<u8>,
    pub cover_thickness_m: Vec<f64>,
    pub cover_base_m: Vec<f64>,
    pub cover_sc: f64,
    pub cover_k: f64,
    pub min_cover_m: f64,
}

impl ExposedLithology for CoverOverBasementV1 {
    #[inline]
    fn at(&self, i: usize, z: f64) -> (f64, f64, u8) {
        if self.cover_thickness_m[i] >= self.min_cover_m && z >= self.cover_base_m[i] {
            (self.cover_sc, self.cover_k, LITHO_COVER)
        } else {
            (self.basement_sc[i], self.basement_k[i], self.basement_code[i])
        }
    }
}

/// T0 fields the stand-in reads (from the T0 terrain artifact already named upstream).
pub struct T0Rock<'a> { pub rock_family: &'a [u8], pub crust_kind: &'a [u8], pub sediment_m: &'a [f32] }

/// Build the stand-in for posts at (lat, lon) degrees; `pos(i)` gives a post's position.
pub fn build_standin(table: &LithologyV1, t0: &T0, rock: &T0Rock, n_posts: usize, pos: &dyn Fn(usize) -> R<(f64, f64)>) -> R<CoverOverBasementV1> {
    table.validate()?;
    let code_of = |c: u32| basement_code(rock.rock_family[c as usize], rock.crust_kind[c as usize]);
    let sc_of = |c: u32| table.sc(code_of(c));
    let k_of = |c: u32| table.k(code_of(c));
    let sed = |c: u32| (rock.sediment_m[c as usize] as f64).max(0.0);
    let elev = |c: u32| t0.elevation_m[c as usize] as f64;
    let mut s = CoverOverBasementV1 {
        basement_sc: vec![0.0; n_posts], basement_k: vec![0.0; n_posts], basement_code: vec![0; n_posts],
        cover_thickness_m: vec![0.0; n_posts], cover_base_m: vec![0.0; n_posts],
        cover_sc: table.sc(LITHO_COVER), cover_k: table.k(LITHO_COVER),
        min_cover_m: table.min_cover_thickness_mm as f64 / 1000.0,
    };
    let any_land = |la: f64, lo: f64| bilinear(t0, &|_c| 1.0, la, lo, true) > 0.5;
    for i in 0..n_posts {
        let (la, lo) = pos(i)?;
        let own = t0.cell_of_deg(la, lo);
        s.basement_code[i] = code_of(own);
        if any_land(la, lo) {
            s.basement_sc[i] = bilinear(t0, &sc_of, la, lo, true);
            s.basement_k[i] = bilinear(t0, &k_of, la, lo, true);
            let th = bilinear(t0, &sed, la, lo, true);
            s.cover_thickness_m[i] = th;
            s.cover_base_m[i] = bilinear(t0, &elev, la, lo, true) - th;
        } else {
            s.basement_sc[i] = sc_of(own);
            s.basement_k[i] = k_of(own);
            s.cover_thickness_m[i] = if t0.water(own) == CellWater::Land { sed(own) } else { 0.0 };
            s.cover_base_m[i] = elev(own) - s.cover_thickness_m[i];
        }
    }
    Ok(s)
}

/// Summary of a stand-in over a set of posts (for the window header).
#[derive(Clone, Debug, Default, Serialize)]
pub struct LithoPrepStats {
    pub basement_code_counts: [u64; LITHO_CODES],
    pub cover_posts: u64,
    pub cover_thickness_max_m: f64,
    /// Cover thickness percentiles over the counted posts (design review risk: compare with T1b local relief).
    pub cover_thickness_p50_m: f64,
    pub cover_thickness_p90_m: f64,
    pub basement_sc_min: f64,
    pub basement_sc_max: f64,
    pub basement_k_min: f64,
    pub basement_k_max: f64,
}

pub fn prep_stats(s: &CoverOverBasementV1, free: &dyn Fn(usize) -> bool) -> LithoPrepStats {
    let mut th = Vec::new();
    let mut st = LithoPrepStats { basement_sc_min: f64::MAX, basement_sc_max: f64::MIN, basement_k_min: f64::MAX, basement_k_max: f64::MIN, ..Default::default() };
    for i in 0..s.basement_sc.len() {
        if !free(i) { continue; }
        st.basement_code_counts[s.basement_code[i] as usize] += 1;
        if s.cover_thickness_m[i] >= s.min_cover_m { st.cover_posts += 1; }
        st.cover_thickness_max_m = st.cover_thickness_max_m.max(s.cover_thickness_m[i]);
        th.push(s.cover_thickness_m[i]);
        st.basement_sc_min = st.basement_sc_min.min(s.basement_sc[i]);
        st.basement_sc_max = st.basement_sc_max.max(s.basement_sc[i]);
        st.basement_k_min = st.basement_k_min.min(s.basement_k[i]);
        st.basement_k_max = st.basement_k_max.max(s.basement_k[i]);
    }
    th.sort_by(|a, b| a.total_cmp(b));
    if !th.is_empty() {
        st.cover_thickness_p50_m = th[(th.len() - 1) / 2];
        st.cover_thickness_p90_m = th[((th.len() - 1) as f64 * 0.9).round() as usize];
    }
    st
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cover_is_exposed_above_its_base_and_basement_below() {
        let t = LithologyV1::standin_v1();
        let s = CoverOverBasementV1 { basement_sc: vec![0.75], basement_k: vec![0.7], basement_code: vec![LITHO_METAMORPHIC],
            cover_thickness_m: vec![300.0], cover_base_m: vec![1000.0], cover_sc: t.sc(LITHO_COVER), cover_k: t.k(LITHO_COVER), min_cover_m: 1.0 };
        assert_eq!(s.at(0, 1200.0).2, LITHO_COVER);
        assert_eq!(s.at(0, 999.0), (0.75, 0.7, LITHO_METAMORPHIC));
        let none = CoverOverBasementV1 { cover_thickness_m: vec![0.0], ..s };
        assert_eq!(none.at(0, 5000.0).2, LITHO_METAMORPHIC);
        assert!(t.validate().is_ok());
    }
}

// E3 source-envelope consumer. This is a derived geometry, not stock admission.
pub const ENVELOPE_RULE: &str = "source_envelope_cell_v1: native spherical cell owns its published ordered records, including every positive thin layer; no identity correlation, base interpolation, smoothing, dip clamp or lcont/vs subunit construction. Bases use the source mean anchor and raw uplift gradient at the cell centre. Top record extends up to constructed surface only as named zero-stock construction excess above its retained top. Below the last source base is unsupported and fails accounting. Nonland cells are outside E1 support and retain their actual full strata column, without E1 stock claims. Half-open contacts expose the lower positive unit at its top. Materials remain unresolved ADR047. Native source budget remains parent; geometric occupancy never admits stock.";

#[derive(Clone, Debug, Serialize, serde::Deserialize, PartialEq, Eq)]
pub struct EnvelopeConfig {
    pub digest: String,
    pub rule: String,
    #[serde(skip)]
    pub path: Option<String>,
}
impl EnvelopeConfig {
    pub fn from_path(path: &str) -> R<Self> {
        let mut f = std::fs::File::open(path).map_err(|e| format!("env {path}: {e}"))?;
        let mut h = blake3::Hasher::new();
        std::io::copy(&mut f, &mut h).map_err(|e| e.to_string())?;
        Ok(Self { digest: format!("blake3:{}", h.finalize().to_hex()), rule: ENVELOPE_RULE.into(), path: Some(path.into()) })
    }
    pub fn load(&self) -> R<std::sync::Arc<EnvelopeSet>> {
        use std::io::BufRead;
        use std::sync::{Arc, Mutex, OnceLock};
        static CACHE: OnceLock<Mutex<std::collections::BTreeMap<String, Arc<EnvelopeSet>>>> = OnceLock::new();
        if self.rule != ENVELOPE_RULE { return Err("unknown envelope rule".into()); }
        let path = self.path.as_ref().ok_or("env path missing; supply env=<envelopes.jsonl>")?;
        let cache = CACHE.get_or_init(Default::default);
        let mut cache = cache.lock().map_err(|e| e.to_string())?;
        // A context holds immutable digest-pinned input after its first load.
        let key = format!("{}:{path}", self.digest);
        if let Some(v) = cache.get(&key) { return Ok(v.clone()); }
        let mut set = EnvelopeSet { cells: Default::default() };
        let f = std::fs::File::open(path).map_err(|e| format!("env {path}: {e}"))?;
        let mut reader = std::io::BufReader::new(f);
        let mut line = Vec::new();
        let mut hash = blake3::Hasher::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line).map_err(|e| e.to_string())? == 0 { break; }
            hash.update(&line);
            let e: pf_contracts::column_envelope::SourceColumnEnvelope = serde_json::from_slice(&line).map_err(|e| format!("env record: {e}"))?;
            e.validate()?;
            let cell = match e.planet_key { pf_contracts::planet::PlanetKey::Cell { cell } => cell, _ => return Err("env noncell".into()) };
            if e.layers.is_empty() || e.layers.len() > 256 { return Err(format!("env cell {cell}: unsupported layer count for existing 8-bit published provenance index")); }
            let c = EnvelopeCell { shape: e.support.shape, radius_m: e.support.radius_m, sea_m: e.support.sea_level_m,
                mean_m: e.moments.mean_z_m, relief_m: e.moments.relief_r_m, grad: e.moments.uplift_gradient_east_north,
                area_m2: e.support.whole_cell_area_m2, receiver: e.hydraulics.published_receiver.map(|p| match p { pf_contracts::planet::PlanetKey::Cell { cell } => cell, _ => u32::MAX }),
                layers: e.layers.into_iter().map(|l| EnvelopeLayer { id: l.record_id, alias: l.native_unit_key, h: l.h_retained_m as f64, base: l.base_at_cell_centre_m,
                    depth: l.depth_base_exact.metres(), code: l.lithology_code, index: l.column_index as u16, csr: l.source_csr_index }).collect() };
            if !c.layers.iter().any(|l| l.h > 0.) { return Err(format!("env cell {cell}: no positive unit")); }
            if set.cells.insert(cell, c).is_some() { return Err("duplicate envelope cell".into()); }
        }
        if format!("blake3:{}", hash.finalize().to_hex()) != self.digest { return Err("env digest mismatch".into()); }
        let v = Arc::new(set); cache.insert(key, v.clone()); Ok(v)
    }
}
#[derive(Debug)]
pub struct EnvelopeLayer { pub id: String, pub alias: Option<String>, pub h: f64, pub base: f64, pub depth: f64, pub code: u8, pub index: u16, pub csr: u32 }
#[derive(Debug)]
pub struct EnvelopeCell {
    pub shape: [u32; 2], pub radius_m: f64, pub sea_m: f64, pub mean_m: f64, pub relief_m: f64, pub grad: [f64; 2], pub area_m2: f64,
    pub receiver: Option<u32>, pub layers: Vec<EnvelopeLayer>,
}
#[derive(Debug)]
pub struct EnvelopeSet { pub cells: std::collections::BTreeMap<u32, EnvelopeCell> }
impl EnvelopeCell {
    pub fn tilt(&self, cell: u32, lat: f64, lon: f64) -> f64 {
        let la = -90. + (cell / self.shape[1]) as f64 * 180. / self.shape[0] as f64 + 90. / self.shape[0] as f64;
        let lo = -180. + (cell % self.shape[1]) as f64 * 360. / self.shape[1] as f64 + 180. / self.shape[1] as f64;
        let east = self.radius_m * la.to_radians().cos() * ((lon-lo+540.).rem_euclid(360.)-180.).to_radians();
        let north = self.radius_m * (lat-la).to_radians();
        self.grad[0]*east + self.grad[1]*north
    }
    /// Positive layers only; thin records stay. Top extension is construction excess, never source stock.
    pub fn exposed(&self, tilt: f64, z: f64) -> Option<&EnvelopeLayer> {
        let depth = self.mean_m + tilt - z;
        self.layers.iter().find(|l| l.h > 0. && depth < l.depth)
    }
    /// Finite retained interval clipped above the declared sea datum. No top extrapolation in source volume.
    pub fn occupancy(&self, l: &EnvelopeLayer, tilt: f64, surface: f64) -> f64 {
        (surface.min(l.base + tilt + l.h) - (l.base + tilt).max(self.sea_m)).max(0.)
    }
}
impl EnvelopeSet {
    pub fn validate_inputs(&self, t0: &T0, cols: &crate::t1b_rivernet::StrataCols) -> R<()> {
        for c in 0..t0.elevation_m.len() as u32 {
            if t0.water(c) != CellWater::Land { continue; }
            let e = self.cells.get(&c).ok_or_else(|| format!("env missing land cell {c}"))?;
            let i = c as usize;
            if e.shape != [t0.nlat as u32,t0.nlon as u32] || e.radius_m != t0.radius_m || e.sea_m != t0.sea_level_m
                || e.mean_m != t0.elevation_m[i] as f64 || e.relief_m != t0.subgrid_relief_m[i] as f64
                || e.grad != [cols.grad_e[i] as f64, cols.grad_n[i] as f64] || e.receiver != t0.flow(c)
                || e.layers.len() != (cols.off[i+1]-cols.off[i]) as usize {
                return Err(format!("env differs from native T0/strata cell {c}"));
            }
            for l in &e.layers {
                let k = cols.off[i] as usize + l.index as usize;
                if l.csr as usize != k || l.code != cols.lith[k] || l.h != cols.th[k] as f64 { return Err(format!("env layer mismatch {c}:{}", l.index)); }
            }
        }
        Ok(())
    }
}
