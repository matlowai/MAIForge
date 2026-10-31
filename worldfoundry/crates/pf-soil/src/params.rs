//! Declared constants of the zonation rules, each with its unit and source class:
//! [V] verified here (read from a record or exact maths), [R] recalled literature (not re-checked in this lane),
//! [A] project assumption awaiting review (a value picked inside a recalled range, or a presentation choice).
//! Every constant is serialized with its class and basis into the output JSON and the recipe digest.

use serde_json::{json, Value};

macro_rules! params {
    ($( $name:ident : $val:expr, $unit:expr, $class:expr, $basis:expr; )*) => {
        /// Zonation constants (all f64; integers are stored as exact small floats).
        #[derive(Clone, Debug, PartialEq)]
        pub struct SoilParams { $(pub $name: f64,)* }
        impl SoilParams {
            pub fn v1() -> Self { SoilParams { $($name: $val,)* } }
            /// Override one constant by name (CLI `k=v`); unknown names are an error.
            pub fn set(&mut self, k: &str, v: f64) -> Result<(), String> {
                if !v.is_finite() { return Err(format!("param {k}: non-finite value")); }
                match k { $(stringify!($name) => { self.$name = v; Ok(()) })* _ => Err(format!("unknown param {k}")) }
            }
            /// The constants with unit, class and basis (published in the output JSON).
            pub fn table(&self) -> Value {
                let mut v: Vec<Value> = Vec::new();
                $( v.push(json!({"name": stringify!($name), "value": self.$name, "unit": $unit, "class": $class, "basis": $basis})); )*
                Value::Array(v)
            }
        }
    };
}

params! {
    // ---- temperature limits ----
    t_gs_min_c: 6.4, "degC", "R", "Paulsen and Koerner 2014 (Alpine Botany 124, global treeline model): trees need a mean growing-season air temperature >= 6.4 C (recalled, not re-checked)";
    lgs_min_days: 94.0, "d", "R", "Paulsen and Koerner 2014: and a growing season of >= 94 days (recalled, not re-checked)";
    t_season_c: 0.9, "degC", "R", "Paulsen and Koerner 2014: a growing-season day has a daily mean above 0.9 C (recalled); applied here to the climate's one-harmonic seasonal cycle, with no daily or synoptic variance [A]";
    ecotone_m: 150.0, "m", "A", "treeline ecotone (krummholz, tree islands) band above the tree line; recalled order 100-200 m (Koerner 2012, Alpine Treelines), value assumed";
    swi_meadow: 12.0, "degC month", "R", "CAVM bioclimate subzone C/D boundary (summer warmth index 12; Walker et al. 2005, recalled): above it closed (>= 50 %) dwarf-shrub / graminoid tundra, the alpine meadow class here";
    cav_swi_a: 6.0, "degC month", "R", "CAVM subzone A (polar desert, < 5 % vascular cover) below SWI 6 (Walker et al. 2005; CAVM Team 2003; recalled)";
    cav_cov_a: 0.05, "1", "R", "cover at the subzone A/B boundary (upper end of the < 5 % range)";
    cav_swi_b: 9.0, "degC month", "R", "CAVM subzone B (5-25 % cover) to SWI 9 (recalled)";
    cav_cov_b: 0.25, "1", "R", "cover at the B/C boundary (upper end of 5-25 %)";
    cav_swi_c: 12.0, "degC month", "R", "CAVM subzone C (5-50 %) to SWI 12 (recalled)";
    cav_cov_c: 0.5, "1", "R", "cover at the C/D boundary (upper end of 5-50 %)";
    cav_swi_d: 20.0, "degC month", "R", "CAVM subzone D (50-80 %) to SWI 20 (recalled)";
    cav_cov_d: 0.8, "1", "R", "cover at the D/E boundary (upper end of 50-80 %)";
    cav_swi_e: 35.0, "degC month", "R", "CAVM subzone E (80-100 %) to SWI 35, the forest-tundra boundary (recalled); piecewise-linear between the points [A]";
    // ---- nival limit (perennial snow line) ----
    ohm_a: 645.0, "mm/yr", "R", "Ohmura, Kasser and Funk 1992 (J. Glaciol. 38): at the glacier equilibrium-line altitude P = a + b T + c T^2, P annual precipitation (mm), T mean JJA air temperature (C) (recalled, not re-checked)";
    ohm_b: 296.0, "mm/yr/K", "R", "Ohmura et al. 1992, linear coefficient (recalled)";
    ohm_c: 9.0, "mm/yr/K2", "R", "Ohmura et al. 1992, quadratic coefficient (recalled)";
    jja_factor: 0.900316316157106, "1", "V", "mean of cos over the warmest quarter (+-45 deg) of the climate's sinusoidal cycle: (2 / pi) x 2 sin(pi / 4); T_JJA = T_mean + jja_factor x half range (exact for the artifact's declared cycle)";
    nival_ramp_m: 50.0, "m", "A", "half width of the soft transition across the snow line (patchy late-lying snow)";
    // ---- moisture ----
    ai_lo: 0.05, "1", "R", "UNEP aridity index P / PET: hyper-arid below 0.05 (Middleton and Thomas 1997, recalled); cover ramp start";
    ai_hi: 0.65, "1", "R", "UNEP humid above 0.65 (recalled); the same 0.65 as pf-flight climate_look_proxy_v1 [V project consistency]; linear ramp between [A]";
    twi_wet_ref: 8.0, "ln(m)", "A", "topographic wetness above which the post gains local water (converging, gentle ground)";
    twi_wet_gain: 0.25, "1/ln(m)", "A", "local aridity index gain per TWI unit above twi_wet_ref";
    pet_sun_gain: 0.5, "1", "A", "PET scales as 1 + gain (R_sun - 1): sun-facing slopes dry faster";
    marsh_twi: 13.0, "ln(m)", "A", "mire / marsh where TWI >= this on gentle ground (32 m TWI; measured distributions in the lane report)";
    marsh_slope_deg: 3.0, "deg", "A", "mire / marsh only on ground gentler than this";
    riparian_twi: 11.0, "ln(m)", "A", "riparian (valley-floor wet meadow / gallery) where TWI >= this";
    wet_cover: 0.9, "1", "A", "cover of wet classes (marsh, riparian, salt marsh) before the temperature, volcanic and kill factors";
    twi_tan_min: 0.001, "1", "A", "floor of tan(beta) in TWI (flats)";
    // ---- aspect / insolation ----
    k_sun_k: 1.0, "K", "A", "warming per unit of (R_sun - 1); treeline elevations differ by order 100 m between sunny and shaded slopes (recalled, Koerner 2012), value assumed";
    dt_sun_max_k: 2.0, "K", "A", "clamp of the aspect temperature term";
    insol_decl_n: 12.0, "1", "V", "numerical quadrature: declination samples over the warm half year (exact quadrature count, not physics)";
    insol_hours_n: 48.0, "1", "V", "numerical quadrature: hour-angle samples per day";
    tilt_deg: 23.44, "deg", "A", "axial tilt used for the warm-half-year declinations; Earth value; the climate artifact's own axial_tilt_deg is recorded alongside when present";
    // ---- substrate ----
    slope_soil_lo_deg: 38.0, "deg", "A", "soil / vegetation mantle starts to thin above this slope; basis: soil-mantled hillslopes approach a threshold gradient S_c ~ 1.0-1.25 (45-51 deg; Roering et al. 1999, recalled)";
    slope_soil_hi_deg: 55.0, "deg", "A", "mantle gone above this slope (bare rock except crevice plants)";
    slope_floor: 0.05, "1", "A", "cover floor on the steepest ground (crevice plants, lichen not counted)";
    f_cliff: 0.05, "1", "A", "cover factor on cliff / riser regime posts (ct_regime 2)";
    f_talus: 0.7, "1", "A", "cover factor on talus / bench-at-repose posts (ct_regime 1): T1b's regime 1 marks repose-angle hillslopes, not active debris; stable old screes carry sward or forest in humid climates, active ones are bare (recalled qualitative); active rockfall debris enters through pf-talus debris_cover at 2 m. 0.4 gave 31 % scree class on C7 site 1 (first run, reported)";
    f_tread: 0.8, "1", "A", "cover factor on 2 m structural treads (ct_regime 3)";
    debris_gain: 0.9, "1", "A", "cover loss per unit of pf-talus debris_cover (fresh rockfall debris)";
    // ---- volcanic surfaces (unit age from the window's VE records) ----
    tau_lava_yr: 1000.0, "yr", "A", "colonisation time scale of lava / composite-body surfaces at t_ref; recalled chronosequences: closed vascular cover on lava in centuries in warm wet climates (Hawaii) to millennia in cool or dry ones (Etna, Iceland, Mt Fuji); value assumed";
    tau_cinder_yr: 600.0, "yr", "A", "cinder / scoria and tuff surfaces (porous, hold water; recalled faster than blocky lava)";
    tau_tephra_yr: 150.0, "yr", "A", "fine tephra / pyroclastic apron / crater floor surfaces: revegetate in decades (recalled: Mount St Helens 1980 pumice plains, Krakatau)";
    tau_lahar_yr: 100.0, "yr", "A", "lahar deposits (mixed fines, water): decades (recalled)";
    tau_dome_yr: 2000.0, "yr", "A", "lava dome, crater wall and neck surfaces (massive, steep)";
    age_weibull_k: 2.0, "1", "A", "f_age = 1 - exp(-(t / tau)^k): k = 2 gives the lag of primary succession (lichens and fines before vascular cover)";
    q10: 2.0, "1", "R", "temperature sensitivity of biological rates, Q10 ~ 2 (recalled general value)";
    t_ref_c: 15.0, "degC", "A", "growing-season temperature at which the tau values apply";
    active_age_yr: 11700.0, "yr", "R", "a volcano with a unit younger than the Holocene (11.7 ky) counts as active (Smithsonian GVP convention, recalled)";
    vent_ring_m: 1000.0, "m", "A", "vent kill zone of an active volcano: cover rises from 0 at the crater rim to full over this distance (persistent degassing and small ash falls; recalled qualitative: Arenal, Villarrica, Sakurajima summits)";
    kill_r_vei3_m: 2000.0, "m", "A", "radius inside which a VEI 3 fall buries ground vegetation (fall thickness >~ 0.3 m; recalled exponential isopach thinning, Pyle 1989; values assumed)";
    kill_r_vei4_m: 5000.0, "m", "A", "same, VEI 4";
    kill_r_vei5_m: 12000.0, "m", "A", "same, VEI 5";
    kill_r_vei6_m: 30000.0, "m", "A", "same, VEI >= 6";
    // ---- coast and lakes ----
    tide_m: 1.0, "m", "A", "half tidal range (T1b models no tides); land posts up to this height above sea next to the shore are intertidal foreshore";
    saltmarsh_band_m: 1.5, "m", "A", "salt marsh from tide_m up to tide_m + this (between mean high water and storm / spring high water)";
    saltmarsh_slope_deg: 2.0, "deg", "A", "salt marsh only on flats";
    saltmarsh_shore_m: 500.0, "m", "A", "salt marsh only this close to the sea shoreline (T1b shore_mm)";
    intertidal_shore_m: 200.0, "m", "A", "intertidal foreshore only this close to the shoreline";
    lake_marsh_m: 1.0, "m", "A", "fresh marsh on lake fringes up to this height above the lake stage";
    lake_marsh_posts: 2.0, "post", "A", "lake fringe search radius at 32 m (Chebyshev posts)";
    riparian_2m_m: 10.0, "m", "A", "2 m tier: riparian band beside carved channel posts";
    // ---- soil depth proxy (not a regolith model) ----
    soil_dmax_m: 1.0, "m", "A", "soil depth proxy ceiling on stable gentle ground";
    tau_soil_yr: 30000.0, "yr", "V", "project constant: pf-detail bed_styles params tau_soil_ky = 30 (itself a project assumption there); soil on a volcanic unit grows as 1 - exp(-t / tau)";
    soil_slope_lo_deg: 25.0, "deg", "A", "soil proxy thins above this slope";
    soil_slope_hi_deg: 50.0, "deg", "A", "soil proxy zero above this slope";
    // ---- 2 m tier ----
    ring_2m_posts: 8.0, "post", "A", "2 m tier: stencil radius (16 m insolation baseline, riparian search); posts within this many posts of the block edge are published as no-data";
    substrate_smooth_2m_posts: 4.0, "post", "A", "2 m tier: the regime cover factor is the box mean over +-this many posts (8 m): vegetation does not follow pf-detail's 32 m-grid-aligned regime edges post by post (first 2 m run showed straight cover steps there); classes keep the post's own regime";
    debris_smooth_2m_posts: 2.0, "post", "A", "2 m tier: pf-talus debris cover box-averaged over +-this many posts (4 m) before the cover factor and the scree test (first run: salt-and-pepper scree speckle at the debris toe)";
    insol_baseline_2m_posts: 8.0, "post", "A", "2 m tier: half baseline of the insolation / climate slope (central difference over +-8 posts = 32 m, the 32 m tier's scale is +-1 post = 64 m)";
}

/// Climate-zone class thresholds below the tree line (on f_W): [A] presentation classes.
pub const CLOSED_FW: f64 = 0.8;
pub const WOODLAND_FW: f64 = 0.4;
pub const OPEN_FW: f64 = 0.2;
