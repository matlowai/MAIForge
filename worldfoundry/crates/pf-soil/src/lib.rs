//! `pf-soil` (lane SOIL-VEG-1, 2026-10-22): a vegetation zonation field built only from fields that already exist.
//! LABEL: pre-signoff, fidelity tier Preview, derived, never canon; an ecology / presentation proxy (the kernel never
//! reads it for gameplay). Opt-in: nothing else in the chain reads it until a consumer is switched over.
//!
//! PHYSICS (per post; constants and their source classes in `params.rs`):
//! 1. Temperature. The canonical T1b climate (seed-climate-v1) reduced to sea level per T0 cell (pf-flight's
//!    climate_look_proxy_v1 convention), bilinear between cell centres, lapsed to the post height at the artifact's own
//!    lapse rate, plus an aspect term k_sun (R_sun - 1). The seasonal cycle is the artifact's one harmonic
//!    T(t) = T_mean + half_range cos(2 pi t).
//! 2. Tree line (Paulsen and Koerner 2014): trees where the mean temperature of the growing season (days above 0.9 C)
//!    is >= 6.4 C and the season lasts >= 94 days; `treeline_margin_m` is the post's height below its own tree line.
//! 3. Above the tree line: an ecotone band, then the alpine / tundra cover from the summer warmth index (sum of the
//!    positive monthly means) on the CAVM bioclimate-subzone cover curve (Walker et al. 2005).
//! 4. Snow line (Ohmura et al. 1992): the glacier equilibrium line where P = 645 + 296 T_JJA + 9 T_JJA^2;
//!    above it the post is nival (perennial snow / ice or frost barren), with a soft +-50 m ramp.
//! 5. Moisture: the UNEP aridity index P / PET, raised on converging ground by the topographic wetness index
//!    ln(a / tan b) (Beven and Kirkby 1979; a from T1b's climate-key effective drainage area) and lowered on sunny
//!    slopes (PET grows with insolation); mires on high-TWI flats and lake fringes, riparian strips on channels.
//! 6. Aspect: warm-half-year clear-sky direct-beam ratio slope / flat at the post's latitude (numerical day and
//!    declination quadrature, no horizon shading).
//! 7. Substrate and disturbance: the cliff / talus / tread regime (T1b ct_regime at 32 m, pf-detail's 2 m regime),
//!    a slope taper toward the soil-mantle threshold gradient, pf-talus debris cover (when present), the age of the
//!    exposed volcanic unit (lane VE records: lava, cinder, tephra, lahar, dome) through a Weibull colonisation law
//!    whose time scale grows in the cold (Q10) and the dry, a vent ring on active volcanoes and tephra burial
//!    zones around explosive events (recovering with age).
//! 8. Coast: intertidal foreshore and salt marsh bands by height above sea next to the T1b shoreline (no tides are
//!    modelled upstream: tide_m is declared), beaches from the T1b over-code.
//!
//! Determinism: no keyed draws at all; every post is a pure function of world-keyed inputs (chart position, published
//! heights, codes and records); f64 in a fixed order; outputs quantized. Seams: see `tier32` and `tier2`.
//! Not available (documented, not modelled): soil depth beyond an age / slope proxy (no regolith model), fire,
//! grazing, wind exposure, horizon shading, per-eruption tephra isopachs, glacier ice extent (no glacier model in T1b).

pub mod climate;
pub mod output;
pub mod params;
pub mod snowpack;
pub mod tier2;
pub mod tier32;
pub mod window;
pub mod zonation;

pub const GENERATOR: &str = "pf-soil";
pub const GENERATOR_VERSION: &str = "0.1.0";
pub const SCHEMA: &str = "vz.soil-veg/v1";
pub const LABEL: &str = "pre-signoff; fidelity tier Preview; derived, never canon; vegetation zonation proxy from existing fields (lane SOIL-VEG-1)";

pub const RULES: [(&str, &str); 10] = [
    ("temperature", "t_sl_lapse_v1: T_sl(c) = temperature_mean_c(c) + lapse max(max(elevation, fill_level)(c) - sea, 0) per T0 cell (pf-flight climate_look_proxy_v1); at a post T_mean = bilinear(T_sl) - lapse max(h - sea, 0) + clamp(k_sun (R_sun - 1), +-dt_sun_max); half range bilinear; climate tier seed-climate-v1 only"),
    ("treeline", "paulsen_koerner_v1: on T(t) = T_mean + A cos(2 pi t) the season is the part above t_season_c; trees iff season >= lgs_min_days and its mean >= t_gs_min_c; treeline_margin = (T_mean - T_mean*(A)) / lapse with T_mean* the bisected threshold"),
    ("alpine", "cavm_swi_v1: above the tree line f_T = CAVM cover curve of SWI = sum_m max(0, T_mean + A cos(2 pi m / 12)); ecotone class within ecotone_m above the tree line; alpine meadow class where SWI >= swi_meadow, fell-field below"),
    ("snowline", "ohmura_ela_v1: T_JJA = T_mean + jja_factor A; ELA where T_JJA = T*(P) solving P = a + b T + c T^2 with P = bilinear annual precipitation; nival_margin = (T_JJA - T*) / lapse; f_snow = smoothstep(-ramp, +ramp, nival_margin); nival class where nival_margin < 0"),
    ("moisture", "aridity_twi_v1: AI_loc = AI (1 + twi_wet_gain max(TWI - twi_wet_ref, 0)) / (1 + pet_sun_gain (R_sun - 1)); f_W = clamp((AI_loc - ai_lo) / (ai_hi - ai_lo), 0, 1); TWI = mean over the 3 x 3 posts of ln(max(A_eff, 1)) - ln(32 m max(|grad z|, twi_tan_min)) at 32 m (geometric-mean area: a multiple-flow-direction stand-in; central differences over +-1 post)"),
    ("insolation", "direct_beam_ratio_v1: sum over insol_decl_n declinations of the warm half year (+-tilt cos u, u uniform over (-pi/2, pi/2), hemisphere of the post) and insol_hours_n hour angles of max(0, s . n) with the sun up, divided by the same sum of max(0, s . up); gradient: 32 m central differences over +-1 post, 2 m over +-insol_baseline_2m_posts"),
    ("substrate", "regime_slope_v1: f_S = min(regime factor (soil 1, talus f_talus, cliff f_cliff, tread f_tread), 1 - (1 - slope_floor) smoothstep(slope_soil_lo, slope_soil_hi, slope)) x (1 - debris_gain debris_cover); 2 m tier: the regime factor and the debris cover are box means over +-substrate_smooth_2m_posts / +-debris_smooth_2m_posts (classes keep the post's own regime)"),
    ("volcanic", "unit_age_weibull_v1: on edifice posts (over 1) the exposed lane VE unit (edifice_unit; unit 0 = composite body aged at the record's history window T_h) gives age t and style; f_age = 1 - exp(-(t / tau')^k), tau' = tau(style) q10^((t_ref - max(T_gs, 0)) / 10) max(1, ai_hi / AI); kill: active volcano (youngest unit < active_age_yr) f *= smoothstep(r_crater, r_crater + vent_ring, r); each unit with VEI >= 3: f *= 1 - (1 - f_age(t, tau_tephra)) (1 - smoothstep(r_k / 2, r_k, r)), r_k by VEI"),
    ("wet_and_coast", "wet_coast_v1: intertidal = h - sea <= tide_m and shore <= intertidal_shore_m (cover 0); salt marsh = tide_m < h - sea <= tide_m + saltmarsh_band_m, slope <= saltmarsh_slope_deg, shore <= saltmarsh_shore_m; marsh = lake fringe (0 < h - lake stage <= lake_marsh_m; stage = highest lake post within lake_marsh_posts at 32 m, max of the 4 bilinear corners at 2 m) or TWI >= marsh_twi on slope <= marsh_slope_deg; riparian = channel post (32 m), bank / within riparian_2m_m of a carved channel (2 m) or TWI >= riparian_twi; wet cover = max(base, wet_cover f_T f_snow f_age f_kill)"),
    ("class_priority", "class_priority_v1: water > beach > intertidal > nival > kill zone (f_kill < 0.5) > fresh volcanic (f_age < 0.5) > bare rock (cliff or slope factor < 0.3, not wet) > scree (talus with cover < 0.5 or debris >= 0.5, not wet) > salt marsh > marsh > riparian > below tree line (closed f_W >= 0.8, woodland >= 0.4, grassland >= 0.2, dry sparse) > ecotone > alpine meadow > fell-field"),
];

pub fn rules_json() -> serde_json::Value {
    serde_json::Value::Object(RULES.iter().map(|(k, v)| (k.to_string(), serde_json::Value::String(v.to_string()))).collect())
}

pub fn classes_json() -> serde_json::Value {
    serde_json::Value::Array(zonation::CLASS_NAMES.iter().map(|(c, n, d)| serde_json::json!({"code": c, "name": n, "rule": d})).collect())
}

pub fn recipe_digest(v: &serde_json::Value) -> String {
    format!("blake3:{}", blake3::hash(&serde_json::to_vec(v).expect("json")).to_hex())
}
