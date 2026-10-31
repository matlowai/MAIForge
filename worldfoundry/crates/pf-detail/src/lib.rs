//! T2 walkable ground (ADR 043 terrain tier T2, node TER-DETAIL; lane T2, 2026-10-13).
//!
//! LABEL: pre-signoff; fidelity tier Preview; derived, never canon. SimTime of the
//! relaxation is `t2-relaxation-not-history` (a numerical device, not world time).
//!
//! One walk tile = TileKeyV2 lod 0 (256 m, 129 x 129 edge-inclusive posts at 2 m,
//! ADR 044 section 3). Pipeline (PREREG part B in `unpublished-design-record`):
//! B1 Catmull-Rom bicubic of the T1b relief; B2 stratigraphy cut through the T0
//! column (ADR 043 2026-10-13 amendment: T1b/T2 may read the stack); B3
//! differential-resistance remap (slope proportional to rock strength); B4 keyed
//! roughness closure (`sub-resolution-closure-v1`, ADR 044 keyed draws); B5 tile-local
//! gully relaxation (MFD area, implicit D8 stream power, lithology nonlinear creep);
//! B6 channel cross-sections from declared hydraulic geometry; B7 seam blend of
//! neighbouring tiles' relaxations; B8 tile-mean preservation; B9 no new sinks;
//! B10 integer-millimetre quantization.
//!
//! Every stage except B5 is pointwise in (chart position, upstream data), so a post
//! shared by two tiles gets the same bytes from either side (F9); B5 is made
//! consistent by B7 (each tile's relaxation is a pure function of its key).

pub mod base_interp;
pub mod bedstyle;
pub mod chain;
pub mod channels;
pub mod chn;
pub mod checks;
pub mod de;
pub mod input;
pub mod lobes;
pub mod pyro;
pub mod relax;
pub mod relief;
pub mod sharp;
pub mod strata;
pub mod tile;
pub mod v2;
pub mod v3;
/// Lane PV (ROCK-LOOK-PROVENANCE, 2026-10-15): the exposed-layer provenance core, the same module as pf-region's so both
/// tiers encode the same id for the same layer. Lane BS (2026-10-15): a crate dependency instead of the `#[path]` share.
pub use pf_region::strata_prov;

use serde::Serialize;

pub type R<T> = Result<T, String>;

pub const GENERATOR_NAME: &str = "t2-ground-detail";
pub const GENERATOR_VERSION: &str = "0.1.0";
pub const LAYER: &str = "t2.ground";
pub const SIMTIME_LABEL: &str = "t2-relaxation-not-history";
pub const CLOSURE_LABEL: &str = "sub-resolution-closure-v1";
pub const LABEL: &str = "pre-signoff; T0 substrate pre-signoff (world_strata_v2 opt-in profile); fidelity tier Preview; derived, never canon; SimTime t2-relaxation-not-history";
pub const PURPOSE_ROUGHNESS: &str = "t2.closure.roughness";
pub const PURPOSE_BEDDING: &str = "t2.closure.bedding";
/// Bed median thickness (m) per LithologyCodeV1 class; 0 = unbedded (Part B revision 1 R1; project assumptions).
pub const BED_MEDIAN_M: [f64; 12] = [0.0, 0.0, 0.0, 6.0, 3.0, 1.5, 2.0, 3.0, 2.0, 1.5, 0.0, 0.0];

/// Posts per tile side (edge-inclusive) and cells per tile side.
pub const TILE_POSTS: usize = 129;
pub const TILE_CELLS: i64 = 128;
pub const POST_M: f64 = 2.0;
pub const TILE_M: f64 = 256.0;
/// T1b post spacing (lod 4).
pub const T1B_M: f64 = 32.0;
/// T1a node spacing (lod 8) in T1b posts.
pub const L8_TO_L4: i64 = 16;
pub const SECONDS_PER_YEAR: f64 = 31_557_600.0;

/// LithologyV1 `standin_v1` (copied verbatim from pf-region `t1b_litho.rs`, lane T3, 2026-10-13;
/// project assumptions until cited; see that file for the sources). Code order = pf-contracts
/// `LithologyCodeV1`.
pub const LITHO_NAMES: [&str; 12] = [
    "crystalline basement", "granite pluton", "metamorphic", "basalt / volcanic", "marine carbonate", "shale / mudstone",
    "sandstone", "conglomerate", "evaporite", "coal-bearing", "unconsolidated", "oceanic igneous / ophiolite",
];
pub const LITHO_SC_PERMILLE: [i64; 12] = [781, 1000, 700, 1732, 2747, 404, 1428, 1000, 466, 466, 625, 839];
pub const LITHO_K_PERMILLE: [i64; 12] = [500, 400, 800, 600, 500, 2000, 800, 900, 3000, 1500, 3000, 700];
pub fn sc_of(code: u8) -> f64 { LITHO_SC_PERMILLE[code as usize] as f64 / 1000.0 }
pub fn k_of(code: u8) -> f64 { LITHO_K_PERMILLE[code as usize] as f64 / 1000.0 }

/// Declared parameters (PREREG part B). Stage switches exist for F8 (zero invention) tests.
#[derive(Clone, Debug, Serialize)]
pub struct DetailParamsV1 {
    /// Opt-in bounded C1 reconstruction of process-evolved T1b parents.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub base_interp: Option<base_interp::Params>,
    pub halo_posts: usize,
    pub seam_half_width_m: f64,
    pub mean_taper_m: f64,
    // B2
    pub min_layer_m: f64,
    pub dip_clamp_tan: f64,
    // B3
    pub remap_on: bool,
    /// Part B revision 1 R1: keyed bedding closure.
    pub bedding_on: bool,
    pub bedding_sigma_ln: f64,
    pub bedding_mult_sigma_ln: f64,
    pub bedding_max_depth_m: f64,
    pub remap_window_m: f64,
    pub remap_drift_window_m: f64,
    pub remap_table_step_m: f64,
    pub remap_channel_taper_m: f64,
    // B4
    pub closure_rms_cap_m: f64,
    pub closure_rms_fraction_of_relief: f64,
    // B5
    pub relax_on: bool,
    pub steps: usize,
    pub mfd_exponent: f64,
    pub channel_head_area_m2: f64,
    pub k_sp_per_step: f64,
    pub d_creep_m2_per_step: f64,
    pub creep_ratio_cap: f64,
    // B6
    pub channels_on: bool,
    pub hg_width_coeff: f64,
    pub hg_width_exp: f64,
    pub hg_depth_coeff: f64,
    pub hg_depth_exp: f64,
    pub min_width_m: f64,
    pub min_depth_m: f64,
    pub bank_width_m: f64,
    pub bank_slope: f64,
    pub runoff_m_per_yr: f64,
    // B8, B9
    pub mean_on: bool,
    pub fill_on: bool,
    pub fixed_above_sea_m: f64,
    /// v2 (lane T2B, generator 0.2.0): None = v1 (its params JSON, and so its recipe digest, are unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub v2: Option<v2::V2Params>,
    /// Lane PV (opt-in `prov=1`): publish prov_id + the provenance table. None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub prov: Option<ProvParamsV1>,
    /// Lane BS (opt-in `beds=2`, needs version >= 2): record-driven bed styles replace the v2 pairs / v3 packages. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub beds: Option<bedstyle::BedStyleParams>,
    /// Lane BS (opt-in `unit_depth=1`): publish depth_in_unit_m / height_above_unit_base_m per post (output only). None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit_depth: Option<bool>,
    /// Lane DE (opt-in `de=1`, needs version=3 and beds=2): differential erosion by hardness as a process. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub de: Option<de::DeParams>,
    /// Lane PYRO (opt-in `pyro_sub=1`, needs beds=2 and bed_meta=1): publish the pyroclastic subtype of facies 5 / 6 beds;
    /// `pyro_sub_hard=1` also sets the lobe columns' pyroclastic strengths from it. None = off (params JSON unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pyro: Option<pyro::PyroParams>,
    /// Lane SHARP (opt-in `ct_sharp=1`, needs a T1b window with the ct fields and ct_base_mm): 2 m re-profile of T1b's
    /// cliff-and-talus staircase (sharp.rs). None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ct_sharp: Option<sharp::SharpParams>,
    /// Lane CLIMATE-RUNOFF (opt-in `climate_runoff=1`, needs a T1b window built with pf-region climate=<dir>): every B5
    /// post's own drainage area is multiplied by the window's runoff_weight R(x) / r_ref at the post (bilinear between T1b
    /// posts), so the 2 m rill / gully areas are effective areas Q / r_ref like the T1b areas they join. None = off (params
    /// JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub climate_runoff: Option<ClimateDetailParamsV1>,
    /// Lane VALLEY-FLOOR-FIX (opt-in `b8_capfix=1`, `b8_levee=1`): B8 apron-cap fallback and carved-volume levees.
    /// None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vff: Option<VffDetailParamsV1>,
    /// Lane INTS2n4 (opt-in `arm_cl=1`, `ringflow_ch=1`): marine arm on a smooth centreline; ring inflow without channel area
    /// (chn.rs). None = off (params JSON and every byte unchanged).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chn: Option<chn::ChnParams>,
}

/// Lane VALLEY-FLOOR-FIX (2026-10-20, EXPLORATORY; constants are project assumptions awaiting review).
#[derive(Clone, Debug, Serialize)]
pub struct VffDetailParamsV1 {
    /// b8_capfix_v1: see `B8_CAPFIX_RULE`.
    pub capfix: bool,
    /// b8_levee_v1: see `B8_LEVEE_RULE`.
    pub levee: bool,
    /// Levee width E as a multiple of the bankfull width, and its cap (m).
    pub levee_width_widths: f64,
    pub levee_width_max_m: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub capfix_rule: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub levee_rule: Option<String>,
}

pub const B8_CAPFIX_RULE: &str = "b8_capfix_v1 (lane VALLEY-FLOOR-FIX): the T2D X3 apron cap (raise <= apron_cap x |deficit|, weight ceiling by bisection) is infeasible when even a uniform raise over every free post with weight > 0 exceeds the cap (free posts < np / apron_cap: tiles mostly carved channel, measured 57.9 % free in tiles 36_26, 37_27, 44_27, 45_26 of period-ends2 r_pe); the bisection then drove the ceiling to ~2^-60 of the largest weight, so min(phi, ceiling) was the same tiny number on every free post: a uniform 1.15-1.3 m pad with vertical edges and mean_shift_c_m ~1.2e18. With the fix an infeasible cap switches that tile to the most spread-out raise that still fades at the tile edges: weights = the B8 edge taper alone (no channel taper or apron peak, no ceiling), also in the post-fill re-correction; the tile mean is exact as before (same c formula). Uncapped apron weights were tried first and peaked at 3.9-4.5 m beside the channel in those tiles";
pub const B8_LEVEE_RULE: &str = "b8_levee_v1 (lane VALLEY-FLOOR-FIX): after the B6 carve and before the B8 mean shift, the carved cross-section of each channel line is returned beside it as natural levees: per line (nearest point, curves near_lines), A_c = (4/3) hw D + 2 hw h_b (the B6 parabolic section of half-width hw and depth D below the bed, plus the T1b surface standing h_b above the bed at the banks, h_b = mean of the T1b bicubic at the two bank points - bed, >= 0), levee width E = min(c_E 2 hw, E_max), thickness t(e) = H (1 - e / E)^2 for 0 <= e = r - hw <= E with H = 3 A_c / (2 E) (each side holds A_c / 2), C1 to zero at E; a post takes the largest t over the lines near it (no double count at junctions) and only free or bank posts (fixed 0 / 2) are raised. Pure function of position, so seamless across tiles. Mean preservation moves from the 256 m tile to the channel segment: carve and levee balance per segment, so B8 (same weights, same c formula) restores each tile's mean for every OTHER stage (target = B1 base mean - (carved tile mean - levee tile mean)) and leaves the carve-plus-levee imbalance of a tile, which is volume carried across the tile edge along the same segment, in place instead of turning it into a per-tile pad; the tile stats record carve_mean_m and levee_mean_m so the balance can be summed over any set of tiles (a mismatch of the analytic A_c against the real carve shows there)";

impl VffDetailParamsV1 {
    pub fn off() -> Self { VffDetailParamsV1 { capfix: false, levee: false, levee_width_widths: 3.0, levee_width_max_m: 300.0, capfix_rule: None, levee_rule: None } }
    /// The levee width and crest thickness of a line of half-width hw and depth D.
    pub fn levee_shape(&self, hw: f64, depth: f64, bank_m: f64) -> (f64, f64) {
        let e = (self.levee_width_widths * 2.0 * hw).min(self.levee_width_max_m).max(1.0);
        let ac = 4.0 / 3.0 * hw * depth + 2.0 * hw * bank_m;
        (e, 1.5 * ac / e)
    }
}

/// Lane CLIMATE-RUNOFF: the declared T2 climate runoff rule (no tunable constants; the runoff field is T1b's).
#[derive(Clone, Debug, Serialize)]
pub struct ClimateDetailParamsV1 { pub rule: String }

pub const CLIMATE_DETAIL_RULE: &str = "climate_runoff_t2_v1 (lane CLIMATE-RUNOFF): B5 own area of each 2 m post = cell area x runoff factor (de_runoff / de_talus_runoff facies factor, 1 without them) x w, w = the T1b window's runoff_weight (R / r_ref, pf_region climate_runoff.rs) bilinear at the post; transmission losses unchanged; channel widths (B6) read T1b / T1a discharge areas, which are already effective areas Q / r_ref, with runoff_m_per_yr = r_ref";

/// Lane PV: the declared provenance rule (no tunable constants; it only copies upstream records).
#[derive(Clone, Debug, Serialize)]
pub struct ProvParamsV1 {
    pub rule: String,
    /// Search radius (T1b posts) for the T1b volcano_id / over code of an edifice / deposit post whose nearest T1b post has none.
    pub t1b_search_posts: i64,
}

impl ProvParamsV1 {
    pub fn v1() -> Self {
        ProvParamsV1 { rule: format!("{} T2: the exposed layer is the one the B2/v2 cut reports (kept-layer index mapped back to the full-column index by skipping layers thinner than min_layer_m, as the cut does); edifice posts take the volcano index and deposit posts the fan/beach kind of the nearest T1b post that has one within t1b_search_posts", strata_prov::PROV_RULE), t1b_search_posts: 2 }
    }
}

impl DetailParamsV1 {
    pub fn v1() -> Self {
        DetailParamsV1 {
            base_interp: None,
            halo_posts: 32,
            seam_half_width_m: 24.0,
            mean_taper_m: 64.0,
            min_layer_m: 1.0,
            dip_clamp_tan: 0.087_488_663_525_924, // tan 5 deg
            remap_on: true,
            bedding_on: true,
            bedding_sigma_ln: 0.6,
            bedding_mult_sigma_ln: 0.35,
            bedding_max_depth_m: 7000.0,
            remap_window_m: 30.0,
            remap_drift_window_m: 120.0,
            remap_table_step_m: 0.05,
            remap_channel_taper_m: 48.0,
            closure_rms_cap_m: 0.30,
            closure_rms_fraction_of_relief: 0.002,
            relax_on: true,
            steps: 300,
            mfd_exponent: 1.1,
            channel_head_area_m2: 300.0,
            k_sp_per_step: 0.0007,
            d_creep_m2_per_step: 0.04,
            creep_ratio_cap: 0.95,
            channels_on: true,
            hg_width_coeff: 3.5,
            hg_width_exp: 0.5,
            hg_depth_coeff: 0.27,
            hg_depth_exp: 0.4,
            min_width_m: 2.0,
            min_depth_m: 0.2,
            bank_width_m: 6.0,
            bank_slope: 0.5,
            runoff_m_per_yr: 1.0,
            mean_on: true,
            fill_on: true,
            fixed_above_sea_m: 1.0,
            v2: None,
            prov: None,
            beds: None,
            unit_depth: None,
            de: None,
            pyro: None,
            ct_sharp: None,
            climate_runoff: None,
            vff: None,
            chn: None,
        }
    }
    /// v2 (lane T2B PREREG): v1 plus the V1-V7 stages; remap window 60 m of section, drift low-pass 240 m.
    pub fn v2() -> Self {
        let mut p = Self::v1();
        p.remap_window_m = 60.0;
        p.remap_drift_window_m = 240.0;
        p.v2 = Some(v2::V2Params::default());
        p
    }
    /// v3 (lane T2C PREREG): v2 plus the W1-W7 stages.
    pub fn v3() -> Self {
        let mut p = Self::v2();
        p.v2.as_mut().unwrap().v3 = Some(v3::V3Params::default());
        p
    }
    pub fn generator_version(&self) -> &'static str {
        match &self.v2 { Some(v) if v.v3.is_some() => v3::GENERATOR_VERSION_V3, Some(_) => v2::GENERATOR_VERSION_V2, None => GENERATOR_VERSION }
    }
    pub fn v3p(&self) -> Option<&v3::V3Params> { self.v2.as_ref().and_then(|v| v.v3.as_ref()) }
    /// Every inventive stage off: the output must be the quantized bicubic base (F8a).
    pub fn all_off(mut self) -> Self {
        self.remap_on = false;
        self.relax_on = false;
        self.channels_on = false;
        self.mean_on = false;
        self.fill_on = false;
        self.closure_rms_cap_m = 0.0;
        self.bedding_on = false;
        if let Some(v) = self.v2.as_mut() { v.valley_on = false; }
        self
    }
    pub fn set(&mut self, k: &str, v: &str) -> R<()> {
        let f = || v.parse::<f64>().map_err(|e| format!("{k}={v}: {e}"));
        let b = || Ok::<bool, String>(v == "1" || v == "true");
        match k {
            "base_interp" => self.base_interp = match v {
                "c1" => Some(base_interp::Params::default()),
                "off" => None,
                _ => return Err(format!("base_interp={v}: expected c1 or off")),
            },
            "remap" => self.remap_on = b()?,
            "climate_runoff" => self.climate_runoff = if b()? { Some(ClimateDetailParamsV1 { rule: CLIMATE_DETAIL_RULE.into() }) } else { None },
            "b8_capfix" | "b8_levee" => {
                let on = b()?;
                let mut q = self.vff.take().unwrap_or_else(VffDetailParamsV1::off);
                if k == "b8_capfix" { q.capfix = on; q.capfix_rule = if on { Some(B8_CAPFIX_RULE.into()) } else { None }; }
                else { if on && self.v2.is_none() { return Err("b8_levee needs version=2 or 3 (curves)".into()); } q.levee = on; q.levee_rule = if on { Some(B8_LEVEE_RULE.into()) } else { None }; }
                self.vff = if q.capfix || q.levee { Some(q) } else { None };
            }
            "relax" => self.relax_on = b()?,
            "channels" => self.channels_on = b()?,
            "mean" => self.mean_on = b()?,
            "fill" => self.fill_on = b()?,
            "closure_rms_cap_m" => self.closure_rms_cap_m = f()?,
            "bedding" => self.bedding_on = b()?,
            "beds" => {
                match v {
                    "2" => { if self.v2.is_none() { return Err("beds=2 needs version=2 or 3 (give version first)".into()); } self.beds = Some(bedstyle::BedStyleParams::default()); }
                    "0" | "1" => self.beds = None,
                    _ => return Err(format!("beds={v}: 0/1 (off) or 2")),
                }
            }
            "bs_csets" => { let bp = self.beds.as_mut().ok_or_else(|| "bs_csets needs beds=2 (give it first)".to_string())?; bp.csets = if b()? { Some(bedstyle::CsetParams::default()) } else { None }; }
            "bed_meta" => { let bp = self.beds.as_mut().ok_or_else(|| "bed_meta needs beds=2 (give it first)".to_string())?; bp.bed_meta = if b()? { Some(true) } else { None }; }
            "ed_datum" => { let bp = self.beds.as_mut().ok_or_else(|| "ed_datum needs beds=2 (give it first)".to_string())?; bp.ed_datum = if b()? { Some(crate::bedstyle::ED_DATUM_NOTE.to_string()) } else { None }; }
            // lane ARCH RIDGE: the datum from the interfluve envelope (closing along the regional strike over +-L T1b posts); 0 = off
            "ed_env" => { let bp = self.beds.as_mut().ok_or_else(|| "ed_env needs beds=2 (give it first)".to_string())?;
                if bp.ed_datum.is_none() { return Err("ed_env needs ed_datum=1 (give it first)".into()); }
                let h: i64 = v.parse().map_err(|_| format!("ed_env={v}: posts"))?;
                if !(0..=64).contains(&h) { return Err(format!("ed_env={v}: 0..64")); }
                bp.ed_env = if h > 0 { Some(crate::bedstyle::EdEnvParams::new(h)) } else { None }; }
            "ed_env_strike" => { let e = self.beds.as_mut().and_then(|b| b.ed_env.as_mut()).ok_or_else(|| "ed_env_strike needs ed_env (give it first)".to_string())?;
                let t: i64 = v.parse().map_err(|_| format!("ed_env_strike={v}: posts"))?;
                if !(1..=128).contains(&t) { return Err(format!("ed_env_strike={v}: 1..128")); }
                e.strike_half_posts = t; }
            "ed_env_route" => { let e = self.beds.as_mut().and_then(|b| b.ed_env.as_mut()).ok_or_else(|| "ed_env_route needs ed_env (give it first)".to_string())?;
                e.route = if b()? { Some(true) } else { None }; }
            "ed_contacts" => { let bp = self.beds.as_mut().ok_or_else(|| "ed_contacts needs beds=2 (give it first)".to_string())?; bp.ed_contacts = if b()? { Some(true) } else { None }; }
            // lane INTS2n4 (chn.rs)
            "arm_cl" | "ringflow_ch" => {
                let on = b()?;
                let mut c = self.chn.take().unwrap_or(chn::ChnParams { arm: None, ringflow_ch: None });
                if k == "arm_cl" {
                    if on && self.v2.is_none() { self.chn = Some(c); return Err("arm_cl needs version=2 or 3 (curves; give version first)".into()); }
                    c.arm = if on { Some(c.arm.take().unwrap_or_else(chn::ArmParams::new)) } else { None };
                } else {
                    if on && self.de.as_ref().and_then(|d| d.ringflow) != Some(1) { self.chn = Some(c); return Err("ringflow_ch needs de_ringflow=1 (give it first)".into()); }
                    c.ringflow_ch = if on { Some(chn::RINGFLOW_CH_RULE.into()) } else { None };
                }
                self.chn = if c.arm.is_some() || c.ringflow_ch.is_some() { Some(c) } else { None };
            }
            "arm_cl_passes" | "arm_cl_b_m" | "arm_cl_hw_max_m" => {
                let a = self.chn.as_mut().and_then(|c| c.arm.as_mut()).ok_or_else(|| format!("{k} needs arm_cl=1 (give it first)"))?;
                let x = f()?;
                match k {
                    "arm_cl_passes" => { if !(0.0..=64.0).contains(&x) || x.fract() != 0.0 { return Err(format!("{k}={v}: 0..64")); } a.smooth_passes = x as u32; }
                    "arm_cl_b_m" => { if !(x > 0.0) { return Err(format!("{k}={v}: > 0")); } a.b_m = x; }
                    _ => { if !(x > 0.0) { return Err(format!("{k}={v}: > 0")); } a.hw_max_m = x; }
                }
            }
            "bs_rmin_m" => { self.beds.as_mut().ok_or_else(|| "bs_rmin_m needs beds=2".to_string())?.rmin_m = f()?; }
            "de" => {
                if b()? {
                    if self.v3p().is_none() || self.beds.is_none() { return Err("de=1 needs version=3 and beds=2 (give them first)".into()); }
                    self.de = Some(de::DeParams::default());
                } else { self.de = None; }
            }
            "de_talus" | "de_sp" | "de_mb" | "de_lobes" | "de_amp_frac" | "de_amp_cap_m" | "de_outer" | "de_outer2" | "de_damp" | "de_screen_m" | "de_cg_tol" | "de_window_m" | "de_fs" | "de_repose" | "de_ve" | "de_fs_screen_m" | "de_fs_ring" | "de_fs_cap" | "de_fs_lh" | "de_runoff" | "de_ro_top" | "de_ro_core" | "de_ro_base" | "de_ro_pyro" | "de_ro_hard" | "de_ro_ltop" | "de_ro_lcore" | "de_ro_lbase" | "de_ro_lhard" | "de_lobate" | "de_lobate_amp" | "de_ringflow" | "de_stack" | "de_stack_blur" | "de_route" | "de_lens" | "de_seq" | "de_seq_sigma" | "de_seq_pond" | "de_seq_fill" | "de_tephra_rec" | "de_seq_level" | "de_seq_stretch" | "de_lobe_relief" | "de_lobe_relief_kappa" | "de_b5_chain" | "de_b5_chain_cell" | "de_b5_chain_dt" | "de_b5_chain_loss" | "de_b5_chain_snap" | "de_b5_chain_ring" | "de_b5_chain_b8" | "de_b5_chain_interp" | "de_b5_chain_head" | "de_b5_chain_dump" | "de_talus_mantle" | "de_talus_runoff" | "de_talus_f" | "de_talus_l_m" | "de_fail_settle" | "de_fail_settle_cap" => {
                let d = self.de.as_mut().ok_or_else(|| format!("{k} needs de=1 (give it first)"))?;
                match k {
                    "de_talus" => d.talus = if b()? { Some(de::TalusParams::default()) } else { None },
                    "de_sp" => d.sp_area_m2 = if b()? { Some([2000.0, 20000.0]) } else { None },
                    "de_mb" => d.mb = if b()? { Some(true) } else { None },
                    "de_lobes" => d.lobes = if b()? { Some(lobes::LobeParams::default()) } else { None },
                    "de_amp_frac" => d.amp_frac = f()?,
                    "de_amp_cap_m" => d.amp_cap_m = f()?,
                    "de_outer2" => d.outer[1] = f()? as usize,
                    "de_damp" => d.damp = f()?,
                    "de_screen_m" => d.screen_m = f()?,
                    "de_cg_tol" => d.cg_tol = f()?,
                    "de_window_m" => d.window_m = Some(f()?),
                    "de_fs" => d.fs = if b()? { Some(de::FsParams::default()) } else { None },
                    "de_fs_screen_m" | "de_fs_ring" | "de_fs_cap" | "de_fs_lh" => {
                        let fs = d.fs.as_mut().ok_or_else(|| format!("{k} needs de_fs=1 (give it first)"))?;
                        match k { "de_fs_screen_m" => fs.screen_m = f()?, "de_fs_ring" => fs.natural_ring = b()?, "de_fs_cap" => fs.fail_cap = f()? as usize, _ => fs.l_h_m = f()? }
                    }
                    // lane DECAP: the failure step settles past its operation budget (a compute budget, not physics)
                    "de_fail_settle" => { let t = d.talus.as_mut().ok_or_else(|| "de_fail_settle needs the talus step (de_talus=1)".to_string())?; t.settle = if b()? { Some(de::SettleParams::new()) } else { None }; }
                    "de_fail_settle_cap" => { let t = d.talus.as_mut().filter(|t| t.settle.is_some()).ok_or_else(|| "de_fail_settle_cap needs de_fail_settle=1 (give it first)".to_string())?;
                        let x = f()?; if !(x >= t.sweeps as f64 && x <= 1.0e6) { return Err(format!("de_fail_settle_cap={v}: {}..1e6 ops per free post", t.sweeps)); }
                        t.settle.as_mut().unwrap().cap_posts = x as usize; }
                    "de_repose" => { let t = d.talus.as_mut().ok_or_else(|| "de_repose needs the talus step (de_talus=1)".to_string())?; t.repose = if b()? { Some(true) } else { None }; }
                    "de_ve" => { let l = d.lobes.as_mut().ok_or_else(|| "de_ve needs de_lobes=1 (give it first)".to_string())?; l.ve = if b()? { Some(true) } else { None }; }
                    // lane RIB: lobate lens margins (and smoothed VE paths)
                    "de_lobate" => { let l = d.lobes.as_mut().ok_or_else(|| "de_lobate needs de_lobes=1 (give it first)".to_string())?; l.lobate = if b()? { Some(lobes::LobateParams::default()) } else { None }; }
                    // lane DATUM: level stacking of the lobe column; paths on the box-mean (constructional proxy) surface
                    "de_stack" => { let l = d.lobes.as_mut().ok_or_else(|| "de_stack needs de_lobes=1 (give it first)".to_string())?; l.stack = if b()? { Some(lobes::StackParams::default()) } else { None }; }
                    "de_stack_blur" => { d.lobes.as_mut().and_then(|l| l.stack.as_mut()).ok_or_else(|| "de_stack_blur needs de_stack=1 (give it first)".to_string())?.blur_hw_mult = f()?; }
                    "de_route" => { let l = d.lobes.as_mut().ok_or_else(|| "de_route needs de_lobes=1 (give it first)".to_string())?; l.route = if b()? { Some(lobes::RouteParams::default()) } else { None }; }
                    // lane DATUM2: sequential (compensational) emplacement (1 routing + bottom-up column, 2 routing only)
                    "de_seq" => { let l = d.lobes.as_mut().ok_or_else(|| "de_seq needs de_lobes=1 (give it first)".to_string())?;
                        l.seq = match v { "0" => None, "1" => Some(lobes::SeqParams::new(1)), "2" => Some(lobes::SeqParams::new(2)), _ => return Err(format!("de_seq={v}: 0, 1 or 2")) }; }
                    // lane PYRO2: fill the gap above the keyed pile with the youngest thick keyed unit (mode 1 column)
                    "de_seq_fill" => { let q = d.lobes.as_mut().and_then(|l| l.seq.as_mut()).ok_or_else(|| "de_seq_fill needs de_seq (give it first)".to_string())?;
                        if q.mode != 1 { return Err("de_seq_fill needs de_seq=1".into()); }
                        q.fill = if b()? { Some(lobes::SeqFillParams::default()) } else { None }; }
                    // lane PYRO2: tephra mantles from the recorded explosive eruptions
                    "de_tephra_rec" => { if self.pyro.is_none() { return Err("de_tephra_rec needs pyro_sub=1 (give it first)".into()); }
                        let l = d.lobes.as_mut().ok_or_else(|| "de_tephra_rec needs de_lobes=1 (give it first)".to_string())?;
                        l.tephra_rec = if b()? { Some(pyro::TEPHRA_REC_NOTE.to_string()) } else { None }; }
                    // lane ARCH: levelled keyed units in the de_seq=1 column (1 elliptic lava lens excess, 2 flat-topped lobe)
                    "de_seq_level" => { let q = d.lobes.as_mut().and_then(|l| l.seq.as_mut()).ok_or_else(|| "de_seq_level needs de_seq (give it first)".to_string())?;
                        if q.mode != 1 { return Err("de_seq_level needs de_seq=1".into()); }
                        q.level = match v { "0" => None, "1" => Some(lobes::SeqLevelParams::new(1)), "2" => Some(lobes::SeqLevelParams::new(2)), _ => return Err(format!("de_seq_level={v}: 0, 1 or 2")) }; }
                    // lane ARCH: close the levelled pile's gap by stretching the whole keyed pile
                    "de_seq_stretch" => { let q = d.lobes.as_mut().and_then(|l| l.seq.as_mut()).and_then(|q| q.level.as_mut()).ok_or_else(|| "de_seq_stretch needs de_seq_level (give it first)".to_string())?;
                        q.stretch = if b()? { Some(lobes::SEQ_STRETCH_NOTE.to_string()) } else { None }; }
                    // lane LOBE: constructional relief of lane VE's lava flows (relief.rs RELIEF_NOTE)
                    "de_lobe_relief" => { let l = d.lobes.as_mut().ok_or_else(|| "de_lobe_relief needs de_lobes=1 (give it first)".to_string())?;
                        if l.ve != Some(true) { return Err("de_lobe_relief needs de_ve=1 (give it first)".into()); }
                        l.relief = match v { "0" => None, "1" | "2" => { let mut r = relief::ReliefParams::default(); r.mode = if v == "2" { 2 } else { 1 }; Some(r) } _ => return Err(format!("de_lobe_relief={v}: 0, 1 or 2")) }; }
                    "de_lobe_relief_kappa" => { d.lobes.as_mut().and_then(|l| l.relief.as_mut()).ok_or_else(|| "de_lobe_relief_kappa needs de_lobe_relief=1 (give it first)".to_string())?.kappa_m2_yr = f()?; }
                    "de_seq_pond" => { d.lobes.as_mut().and_then(|l| l.seq.as_mut()).ok_or_else(|| "de_seq_pond needs de_seq (give it first)".to_string())?.pond_stop = if b()? { Some(true) } else { None }; }
                    "de_seq_sigma" => { d.lobes.as_mut().and_then(|l| l.seq.as_mut()).ok_or_else(|| "de_seq_sigma needs de_seq (give it first)".to_string())?.base_sigma_m = f()?; }
                    // lane DATUM2: flat-based lava lenses in the level stack (1 elliptic dome, 2 flat-topped lobe profile)
                    "de_lens" => { let st = d.lobes.as_mut().and_then(|l| l.stack.as_mut()).ok_or_else(|| "de_lens needs de_stack=1 (give it first)".to_string())?;
                        st.lens = match v { "0" => None, "1" => Some(lobes::LensParams::new(1)), "2" => Some(lobes::LensParams::new(2)), _ => return Err(format!("de_lens={v}: 0, 1 or 2")) }; }
                    "de_lobate_amp" => { d.lobes.as_mut().and_then(|l| l.lobate.as_mut()).ok_or_else(|| "de_lobate_amp needs de_lobate=1 (give it first)".to_string())?.amp = f()?; }
                    // lane RIB2: domain-ring inflow along T1b flow (1 inward flux, 2 T1b D8 flow-line entry points)
                    "de_ringflow" => d.ringflow = match v { "0" => None, "1" => Some(1), "2" => Some(2), _ => return Err(format!("de_ringflow={v}: 0, 1 or 2")) },
                    // lane RILL: coarse window-wide B5 pass conditions the tile domains (1 ring inflow, 2 + initial surface)
                    "de_b5_chain" => d.chain = match v { "0" => None, "1" => Some(chain::ChainParams::new(1)), "2" => Some(chain::ChainParams::new(2)), "3" => Some(chain::ChainParams::new(3)), _ => return Err(format!("de_b5_chain={v}: 0, 1, 2 or 3")) },
                    "de_b5_chain_interp" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain=2|3 (give it first)"))?; if c.mode < 2 { return Err("de_b5_chain_interp needs de_b5_chain=2|3".into()); }
                        c.interp = match v { "1" => None, "3" => Some(3), "4" => Some(4), _ => return Err(format!("{k}={v}: 1 (bilinear), 3 (Catmull-Rom) or 4 (cubic B-spline)")) }; }
                    "de_b5_chain_head" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain (give it first)"))?; c.head_scale = if b()? { Some(true) } else { None }; }
                    "de_b5_chain_dump" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain (give it first)"))?; c.dump = Some(v.to_string()); }
                    "de_b5_chain_b8" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain (give it first)"))?; c.b8_taper_after_cap = if b()? { Some(true) } else { None }; }
                    "de_b5_chain_ring" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain (give it first)"))?; c.natural_ring = if b()? { Some(true) } else { None }; }
                    "de_b5_chain_snap" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain=2|3 (give it first)"))?; if c.mode < 2 { return Err("de_b5_chain_snap needs de_b5_chain=2|3".into()); } c.snap = if b()? { Some(true) } else { None }; }
                    "de_b5_chain_loss" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain=1|2 (give it first)"))?; c.no_loss = if b()? { None } else { Some(true) }; }
                    "de_b5_chain_cell" | "de_b5_chain_dt" => { let c = d.chain.as_mut().ok_or_else(|| format!("{k} needs de_b5_chain=1|2 (give it first)"))?;
                        let u = f()?; if u < 1.0 || u.fract() != 0.0 { return Err(format!("{k}={v}: a positive integer")); }
                        if k == "de_b5_chain_cell" { c.cell_posts = u as usize } else { c.dt = u as usize } }
                    // lane RIB: runoff-weighted rill drainage area on lobe columns
                    // lane DE-TALUS (need ct_sharp=1 ct_sharp_stage=pre; checked when the context is built)
                    "de_talus_mantle" => { let m = match v { "0" => 0, "1" => 1, "2" => 2, _ => return Err(format!("{k}={v}: 0, 1 or 2")) };
                        if m == 0 { if let Some(r) = d.regime.as_mut() { r.mantle = 0; if r.runoff.is_none() { d.regime = None; } } } else { d.regime.get_or_insert_with(de::RegimeParams::new).mantle = m; } }
                    "de_talus_runoff" => { if b()? { d.regime.get_or_insert_with(de::RegimeParams::new).runoff = Some(de::RegimeParams::RUNOFF_TALUS); }
                        else if let Some(r) = d.regime.as_mut() { r.runoff = None; if r.mantle == 0 { d.regime = None; } } }
                    "de_talus_f" | "de_talus_l_m" => { let r = d.regime.as_mut().and_then(|r| r.runoff.as_mut()).ok_or_else(|| format!("{k} needs de_talus_runoff=1 (give it first)"))?;
                        let x = f()?; if !(x > 0.0) { return Err(format!("{k}={v}: > 0")); } if k == "de_talus_f" { r[0] = x; } else { r[1] = x; } }
                    "de_runoff" => { if d.lobes.is_none() { return Err("de_runoff needs de_lobes=1 (give it first)".into()); } d.runoff = if b()? { Some(de::RunoffParams::default()) } else { None }; }
                    "de_ro_top" | "de_ro_core" | "de_ro_base" | "de_ro_pyro" | "de_ro_hard" | "de_ro_ltop" | "de_ro_lcore" | "de_ro_lbase" | "de_ro_lhard" => {
                        let r = d.runoff.as_mut().ok_or_else(|| format!("{k} needs de_runoff=1 (give it first)"))?;
                        match k { "de_ro_top" => r.f_lava_top = f()?, "de_ro_core" => r.f_lava_core = f()?, "de_ro_base" => r.f_lava_base = f()?, "de_ro_pyro" => r.f_pyro = f()?,
                            "de_ro_ltop" => r.l_lava_top = f()?, "de_ro_lcore" => r.l_lava_core = f()?, "de_ro_lbase" => r.l_lava_base = f()?, "de_ro_lhard" => r.l_pyro_hard = f()?, _ => r.f_pyro_hard = f()? }
                    }
                    _ => d.outer[0] = f()? as usize,
                }
            }
            // lane PYRO: pyroclastic subtype (look) and, opt-in, its strength (terrain)
            "pyro_sub" => {
                if b()? {
                    if self.beds.as_ref().and_then(|b| b.bed_meta) != Some(true) { return Err("pyro_sub=1 needs beds=2 and bed_meta=1 (give them first)".into()); }
                    if self.pyro.is_none() { self.pyro = Some(pyro::PyroParams::default()); }
                } else { self.pyro = None; }
            }
            "pyro_downwind" => {
                let py = self.pyro.as_mut().ok_or_else(|| "pyro_downwind needs pyro_sub=1 (give it first)".to_string())?;
                py.downwind = if b()? { Some(pyro::PYRO_DOWNWIND_NOTE.to_string()) } else { None };
            }
            // lane VH-I: fall dispersal axes on the T0 wind published by the T1b window (ve_wind)
            "pyro_wind" => {
                let py = self.pyro.as_mut().ok_or_else(|| "pyro_wind needs pyro_sub=1 (give it first)".to_string())?;
                py.wind = if b()? { Some(pyro::PyroWind::default()) } else { None };
            }
            "pyro_sub_hard" => {
                let py = self.pyro.as_mut().ok_or_else(|| "pyro_sub_hard needs pyro_sub=1 (give it first)".to_string())?;
                py.hard = if b()? { Some(pyro::PyroHard::default()) } else { None };
            }
            // lane SHARP: 2 m re-profile of the T1b cliff-and-talus staircase
            "ct_sharp" => self.ct_sharp = if b()? { Some(sharp::SharpParams::default()) } else { None },
            "ct_sharp_tol_m" | "ct_sharp_rad" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if k == "ct_sharp_tol_m" { sp.merge_tol_m = f()?; } else { let r = f()?; if !(1.0..=16.0).contains(&r) || r.fract() != 0.0 { return Err(format!("{k}={v}: 1..16")); } sp.radius_posts = r as i64; }
            }
            // lane TREAD: where the 2 m staircase is applied (post: after quantization, the default; pre: before B5 / DE with the
            // riser posts protected)
            "ct_sharp_consensus" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.t3_consensus = b()?;
            }
            // lane PHASE: how the cell's four T1b corners combine at a 2 m post (f: blend the F values, the default; phase: one F
            // from the blended period knots; nearest: the dominant corner only, diagnostic)
            "ct_sharp_blend" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.blend = match v { "f" => None, "phase" => Some(sharp::SHARP_BLEND_PHASE_RULE.to_string()), "nearest" => Some(sharp::SHARP_BLEND_NEAREST_RULE.to_string()), _ => return Err(format!("{k}={v}: f | phase | nearest")) };
            }
            // lane SPLIT-2M: in the phase blend, knots are lent only between corners on the same period (sharp::SHARP_BLEND_PERIOD_RULE)
            "ct_sharp_blend_period" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.blend.as_deref() != Some(sharp::SHARP_BLEND_PHASE_RULE) { return Err(format!("{k} needs ct_sharp_blend=phase (give it first)")); }
                sp.blend_period = if b()? { Some(sharp::SHARP_BLEND_PERIOD_RULE.to_string()) } else { None };
            }
            // lane TREADSOIL: benches next to a soil-mantled T1b corner stay whole up to the mantle edge (ws = 0.5)
            "ct_sharp_soil_edge" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.soil_edge = if b()? { Some(sharp::SHARP_SOIL_EDGE_RULE.to_string()) } else { None };
            }
            "ct_sharp_soil_edge_tw" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.soil_edge.is_none() { return Err(format!("{k} needs ct_sharp_soil_edge=1 (give it first)")); }
                let t = f()?; if !(t >= 0.0 && t <= 10.0) { return Err(format!("{k}={v}: 0..10 (0 = no limit)")); }
                sp.soil_edge_tw = Some(t);
            }
            // lane STEEP-TALUS: talus posts steeper than talus can stand at 8 m are rock (pre stage only; checked at build)
            "ct_sharp_talus_rock" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.talus_rock = if b()? { Some(sharp::SHARP_TALUS_ROCK_RULE.to_string()) } else { None };
                if sp.talus_rock.is_none() { sp.talus_rock_deg = None; }
            }
            "ct_sharp_talus_rock_deg" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.talus_rock.is_none() { return Err(format!("{k} needs ct_sharp_talus_rock=1 (give it first)")); }
                let t = f()?; if !(t > 30.0 && t < 80.0) { return Err(format!("{k}={v}: 30 < deg < 80")); }
                sp.talus_rock_deg = Some(t);
            }
            // lane PERIOD-ENDS-2: exact per-post F from pf-region ct_ends (window fields checked at build)
            "ct_sharp_ends" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.ends = if b()? { Some(sharp::SHARP_ENDS_RULE.to_string()) } else { None };
                if sp.ends.is_none() { sp.ends_label = None; sp.cell_riser = None; sp.cell_riser_m = None; sp.cell_riser_dh = None; sp.cell_riser_donor = None; sp.cell_riser_cap = None; }
            }
            // lane CELL-RISER: same-period posts whose cliff foot fans out re-profile on a neighbour's knots with a pseudo-base shift
            "ct_sharp_cell_riser" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.ends.is_none() || sp.blend.as_deref() != Some(sharp::SHARP_BLEND_PHASE_RULE) { return Err(format!("{k} needs ct_sharp_blend=phase and ct_sharp_ends=1 (give them first)")); }
                sp.cell_riser = if b()? { Some(sharp::SHARP_CELL_RISER_RULE.to_string()) } else { None };
                if sp.cell_riser.is_none() { sp.cell_riser_m = None; sp.cell_riser_dh = None; sp.cell_riser_donor = None; sp.cell_riser_cap = None; }
            }
            "ct_sharp_cell_riser_m" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.cell_riser.is_none() { return Err(format!("{k} needs ct_sharp_cell_riser=1 (give it first)")); }
                let t = f()?; if !(t > 0.0 && t <= 200.0) { return Err(format!("{k}={v}: 0 < m <= 200")); }
                sp.cell_riser_m = Some(t);
            }
            "ct_sharp_cell_riser_dh" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.cell_riser.is_none() { return Err(format!("{k} needs ct_sharp_cell_riser=1 (give it first)")); }
                if !matches!(v, "bilinear" | "bicubic" | "smooth" | "monotone") { return Err(format!("{k}={v}: bilinear | bicubic | smooth | monotone")); }
                sp.cell_riser_dh = if v == sharp::CELL_RISER_DH_DEFAULT { None } else { Some(v.to_string()) };
            }
            // lane CELL-RISER-2: robust donor and shift cap (sharp::CELL_RISER2_RULE)
            "ct_sharp_cell_riser_donor" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.cell_riser.is_none() { return Err(format!("{k} needs ct_sharp_cell_riser=1 (give it first)")); }
                if !matches!(v, "min" | "median" | "pair") { return Err(format!("{k}={v}: min | median | pair")); }
                sp.cell_riser_donor = if v == "min" { None } else { Some(v.to_string()) };
            }
            "ct_sharp_cell_riser_cap" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.cell_riser.is_none() { return Err(format!("{k} needs ct_sharp_cell_riser=1 (give it first)")); }
                let t = f()?; if !(t > 0.0 && t <= 8.0) { return Err(format!("{k}={v}: 0 < cells <= 8")); }
                sp.cell_riser_cap = Some(t);
            }
            // lane CARRY: held riser posts carried by a change low-passed along the riser contour (pre stage only; checked at build)
            "ct_sharp_carry_strike" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.carry_strike = if b()? { Some(sharp::SHARP_CARRY_STRIKE_RULE.to_string()) } else { None };
                if sp.carry_strike.is_none() { sp.carry_strike_m = None; }
            }
            "ct_sharp_carry_strike_m" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.carry_strike.is_none() { return Err(format!("{k} needs ct_sharp_carry_strike=1 (give it first)")); }
                let t = f()?; if !(t >= 2.0 && t <= 32.0) { return Err(format!("{k}={v}: 2 <= m <= 32")); }
                sp.carry_strike_m = Some(t);
            }
            // lane CARRY: cliff-limb posts read the B4 closure low-passed along the T1b contour (pre stage only; checked at build)
            "ct_sharp_limb_lp" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.limb_lp = if b()? { Some(sharp::SHARP_LIMB_LP_RULE.to_string()) } else { None };
                if sp.limb_lp.is_none() { sp.limb_lp_m = None; }
            }
            "ct_sharp_limb_lp_m" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.limb_lp.is_none() { return Err(format!("{k} needs ct_sharp_limb_lp=1 (give it first)")); }
                let t = f()?; if !(t >= 2.0 && t <= 32.0) { return Err(format!("{k}={v}: 2 <= m <= 32")); }
                sp.limb_lp_m = Some(t);
            }
            // lane CARRY2: each staircase post is the mean of the staircase over its 2 m footprint (pre stage only; checked at build)
            "ct_sharp_aa" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.aa = if b()? { Some(sharp::SHARP_AA_RULE.to_string()) } else { None };
                if sp.aa.is_none() { sp.aa_w = None; }
            }
            "ct_sharp_aa_w" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.aa.is_none() { return Err(format!("{k} needs ct_sharp_aa=1 (give it first)")); }
                let t = f()?; if !(t >= 0.5 && t <= 2.0) { return Err(format!("{k}={v}: 0.5 <= posts <= 2")); }
                sp.aa_w = Some(t);
            }
            // lane CARRY2: B5 creep edges at held riser posts carry no flux (pre stage only; checked at build)
            "ct_sharp_riser_noflux" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.riser_noflux = if b()? { Some(sharp::SHARP_RISER_NOFLUX_RULE.to_string()) } else { None };
            }
            // lane CARRY2: the riser fringe's B5 change low-passed along the contour (needs ct_sharp_carry_strike)
            "ct_sharp_carry_fringe" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.carry_strike.is_none() { return Err(format!("{k} needs ct_sharp_carry_strike=1 (give it first)")); }
                sp.carry_fringe = if b()? { Some(sharp::SHARP_CARRY_FRINGE_RULE.to_string()) } else { None };
            }
            "ct_sharp_ends_label" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                if sp.ends.is_none() { return Err(format!("{k} needs ct_sharp_ends=1 (give it first)")); }
                sp.ends_label = if b()? { Some(sharp::SHARP_ENDS_LABEL_RULE.to_string()) } else { None };
            }
            "ct_sharp_stage" => {
                let sp = self.ct_sharp.as_mut().ok_or_else(|| format!("{k} needs ct_sharp=1 (give it first)"))?;
                sp.stage = match v { "post" => None, "pre" => Some(sharp::SHARP_PRE_RULE.to_string()), _ => return Err(format!("{k}={v}: pre | post")) };
            }
            "unit_depth" => { if self.v2.is_none() { return Err("unit_depth needs version=2 or 3".into()); } self.unit_depth = if b()? { Some(true) } else { None }; }
            "zero_closure" => { if b()? { self.closure_rms_cap_m = 0.0; self.bedding_on = false; } }
            "steps" => self.steps = f()? as usize,
            "k_sp" => self.k_sp_per_step = f()?,
            "d_creep" => self.d_creep_m2_per_step = f()?,
            "remap_window_m" => self.remap_window_m = f()?,
            "all_off" => { if b()? { *self = self.clone().all_off(); } }
            "valley" => { self.v2.as_mut().ok_or_else(|| "valley needs version=2".to_string())?.valley_on = b()?; }
            "ed_support" => { self.v2.as_mut().ok_or_else(|| "ed_support needs version=2 or version=3 (give it first)".to_string())?.edifice_support = if b()? { Some(true) } else { None }; }
            "overcodes" => { self.v2.as_mut().ok_or_else(|| "overcodes needs version=2".to_string())?.overcodes_on = b()?; }
            "slope_clamp" => { self.v2.as_mut().ok_or_else(|| "slope_clamp needs version=2".to_string())?.slope_clamp_on = b()?; }
            "amp_cap_m" => { self.v2.as_mut().ok_or_else(|| "amp_cap_m needs version=2".to_string())?.amp_cap_m = f()?; }
            "deperch" => { self.v2.as_mut().ok_or_else(|| "deperch needs version=2".to_string())?.deperch_on = b()?; }
            "amp_frac" => { self.v2.as_mut().ok_or_else(|| "amp_frac needs version=2".to_string())?.amp_frac_of_relief = f()?; }
            "pkg" | "pinch" | "expose" | "relax_win" | "rough_corr" | "massbal" | "massbal_remap" | "taper_c" | "massbal_sigma_m" | "expose_lo" | "pinch_amp" | "apron" | "head_law" | "clamp_h_m" | "bank_grade" | "apron_cap" | "bed_sig_v2" | "thalweg_hw_m" | "bank_fin_m" | "b8_taper_after" => {
                let w = self.v2.as_mut().and_then(|v| v.v3.as_mut()).ok_or_else(|| format!("{k} needs version=3"))?;
                match k {
                    "pkg" => w.pkg_on = b()?,
                    "pinch" => w.pinch_on = b()?,
                    "expose" => w.expose_on = b()?,
                    "relax_win" => w.relax_win_on = b()?,
                    "rough_corr" => w.rough_corr_on = b()?,
                    "massbal" => w.massbal_on = b()?,
                    "massbal_remap" => w.massbal_remap_on = b()?,
                    "taper_c" => w.taper_c_on = b()?,
                    "massbal_sigma_m" => w.massbal_sigma_m = f()?,
                    "expose_lo" => w.expose_lo = f()?,
                    "apron" => w.apron_on = b()?,
                    "head_law" => w.head_law_on = b()?,
                    "clamp_h_m" => w.clamp_h_m = Some(f()?),
                    "bank_grade" => w.bank_grade_max_m = if b()? { Some(0.5) } else { None },
                    "apron_cap" => w.apron_cap = Some(f()?),
                    "thalweg_hw_m" => w.thalweg_hw_m = Some(f()?),
                    "bank_fin_m" => w.bank_fin_m = Some(f()?),
                    "bed_sig_v2" => w.bed_sig_v2 = if b()? { Some(true) } else { None },
                    "b8_taper_after" => w.b8_taper_after = if b()? { Some(true) } else { None },
                    _ => w.pinch_amp = f()?,
                }
            }
            _ => return Err(format!("unknown key {k}")),
        }
        Ok(())
    }
}

#[inline]
pub fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Round half away from zero to integer millimetres.
#[inline]
pub fn to_mm(z_m: f64) -> i32 {
    let v = z_m * 1000.0;
    (if v >= 0.0 { (v + 0.5).floor() } else { (v - 0.5).ceil() }) as i32
}

pub mod contrast;

pub mod escarpment;
