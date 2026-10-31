//! WEATHER-V0: the hourly weather state at a planet point.
//!
//! state = f(T0 climate fields, world seed, world hour, planet point, elevation)
//! and nothing else: no stored state, no chart, no wall clock. See the crate
//! docs for the model and `constants.rs` for every number.

use serde::Serialize;

use crate::climate::{ClimateGrid, ClimateSample};
use crate::constants::WeatherConstantsV0;
use crate::disturb::{DisturbanceField, Hit, ProcessSpec};
use crate::fronts::{FrontField, FrontsConstantsV1};
use crate::sun::{solar_clock, SolarClock};

pub const LAW_ID: &str = "pf-weather-v0";
pub const SCHEMA: &str = "pf.weather-state/v0";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum PrecipType {
    None,
    Drizzle,
    Rain,
    Mixed,
    Snow,
}

impl PrecipType {
    pub fn code(self) -> u8 {
        match self {
            PrecipType::None => 0,
            PrecipType::Drizzle => 1,
            PrecipType::Rain => 2,
            PrecipType::Mixed => 3,
            PrecipType::Snow => 4,
        }
    }
}

/// Extinction budget (1/km); V = koschmieder_ln / total.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Extinction {
    pub rayleigh: f64,
    pub haze: f64,
    pub precip: f64,
    pub dust: f64,
    pub fog: f64,
}

impl Extinction {
    pub fn total(&self) -> f64 { self.rayleigh + self.haze + self.precip + self.dust + self.fog }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct WeatherState {
    pub t_h: f64,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub elevation_m: f64,
    pub t0_cell: u32,
    pub solar_hour: f64,
    pub sun_elevation_deg: f64,
    pub cloud_cover: f64,
    pub precip_type: PrecipType,
    pub precip_rate_mm_h: f64,
    pub convective: bool,
    pub wind_u_m_s: f64,
    pub wind_v_m_s: f64,
    pub wind_speed_m_s: f64,
    pub air_temperature_c: f64,
    pub dewpoint_c: f64,
    pub relative_humidity: f64,
    /// Meteorological visibility of the air mass (km), excluding radiation
    /// fog: what the streamer and the far-view wash use above any pool.
    pub visibility_km: f64,
    /// Visibility on open flat ground including the open-ground radiation-fog
    /// term (km). On relief, use COLDPOOL-1 per post instead.
    pub visibility_open_ground_km: f64,
    pub extinction_per_km: Extinction,
    /// 0..1: how favourable this hour is for radiation cooling (clear, calm, dark).
    pub radiative_cooling: f64,
    /// 0..1: open-ground radiation fog factor included in the fog extinction.
    pub flat_fog: f64,
    /// 0..1: synoptic activity (any system's halo, rain or not).
    pub synoptic_activity: f64,
    /// Expected (climatological) values at this point and hour, for checks.
    pub wet_hour_fraction: f64,
    pub conditional_rate_mm_h: f64,
}

/// Magnus saturation vapour pressure over water (hPa), Alduchov and
/// Eskridge (1996) constants 6.1094, 17.625, 243.04 (recalled).
pub fn es_hpa(t_c: f64) -> f64 { 6.1094 * (17.625 * t_c / (t_c + 243.04)).exp() }
pub fn dewpoint_c(e_hpa: f64) -> f64 {
    let g = (e_hpa / 6.1094).ln();
    243.04 * g / (17.625 - g)
}

pub(crate) fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

pub struct WeatherModel<'a> {
    pub climate: &'a ClimateGrid,
    pub k: WeatherConstantsV0,
    pub world_seed: u64,
    pub field: DisturbanceField,
    pub strat: ProcessSpec,
    pub conv: ProcessSpec,
    /// WEATHER-V1 (lane WEATHER-V1, opt-in): frontal bands and convective
    /// clusters replace the V0 discs when set; None is WEATHER-V0 unchanged.
    pub fronts: Option<FrontField>,
    /// WEATHER-V1 orographic modulation of frontal rain (only with fronts).
    pub oro: bool,
}

impl<'a> WeatherModel<'a> {
    /// The world's own constants (axial tilt, lapse) override the defaults.
    pub fn new(climate: &'a ClimateGrid, world_seed: u64, planet_r_km: Option<f64>) -> Result<WeatherModel<'a>, String> {
        let mut k = WeatherConstantsV0::default();
        if let Some(t) = climate.axial_tilt_deg { k.axial_tilt_deg = t; }
        if let Some(l) = climate.lapse_rate_k_per_m { k.lapse_rate_k_per_m = l; }
        if let Some(r) = planet_r_km { k.planet_radius_km = r; }
        Self::with_constants(climate, world_seed, k)
    }

    pub fn with_constants(climate: &'a ClimateGrid, world_seed: u64, k: WeatherConstantsV0) -> Result<WeatherModel<'a>, String> {
        let strat = ProcessSpec {
            tag: "strat", id: 1, lattice_km: k.strat_lattice_km, slot_h: k.strat_slot_h,
            r_min_km: k.strat_r_min_km, r_max_km: k.strat_r_max_km, d_min_h: k.strat_d_min_h, d_max_h: k.strat_d_max_h,
            mean_cover: k.strat_mean_cover, halo: k.strat_halo, intensity_sigma: k.strat_intensity_sigma,
        };
        let conv = ProcessSpec {
            tag: "conv", id: 2, lattice_km: k.conv_lattice_km, slot_h: k.conv_slot_h,
            r_min_km: k.conv_r_min_km, r_max_km: k.conv_r_max_km, d_min_h: k.conv_d_min_h, d_max_h: k.conv_d_max_h,
            mean_cover: k.conv_mean_cover, halo: k.conv_halo, intensity_sigma: k.conv_intensity_sigma,
        };
        strat.check()?;
        conv.check()?;
        let lmax = -(1.0 - k.wet_hour_fraction_max).ln();
        if k.strat_mean_cover < lmax {
            return Err(format!("strat_mean_cover {} < {lmax}: wet fraction cap unreachable", k.strat_mean_cover));
        }
        if k.conv_mean_cover < k.conv_share_max * (1.0 + k.conv_diurnal_amp) * lmax {
            return Err("conv_mean_cover too small for the convective share".into());
        }
        let field = DisturbanceField::new(world_seed, LAW_ID, k.planet_radius_km, k.steer_tropical_m_s, k.steer_midlat_m_s);
        Ok(WeatherModel { climate, k, world_seed, field, strat, conv, fronts: None, oro: false })
    }

    /// Switch this model to WEATHER-V1 fronts (and optionally the orographic
    /// term). Keyed draws of the new processes use new purposes
    /// (`weather.front.*`, `weather.convcl.*`) under the same law-set scope, so
    /// the V0 hourly jitter draws are shared.
    pub fn enable_fronts(&mut self, kf: FrontsConstantsV1, oro: bool) -> Result<(), String> {
        let lmax = -(1.0 - self.k.wet_hour_fraction_max).ln();
        if kf.front_mean_cover < lmax {
            return Err(format!("front_mean_cover {} < {lmax}", kf.front_mean_cover));
        }
        if kf.cl_mean_cover < self.k.conv_share_max * (1.0 + self.k.conv_diurnal_amp) * lmax {
            return Err("cl_mean_cover too small for the convective share".into());
        }
        let law = vz_keyed_rng::Scope::root(self.world_seed).law_set(LAW_ID).map_err(|e| format!("{e:?}"))?;
        self.fronts = Some(FrontField::new(law, self.k.planet_radius_km, self.k.steer_tropical_m_s, self.k.steer_midlat_m_s, kf)?);
        self.oro = oro;
        Ok(())
    }

    /// Clear the memo of keyed disturbances (no effect on any output).
    pub fn clear_memo(&self) {
        self.field.clear_memo();
        if let Some(f) = &self.fronts { f.clear_memo(); }
    }

    /// Seasonal-mean temperature (no diurnal, no synoptic) at the T0 elevation.
    pub(crate) fn seasonal_t(&self, c: &ClimateSample, sc: &SolarClock, lat_deg: f64) -> f64 {
        let peak = if lat_deg >= 0.0 { self.k.season_peak_day_north } else { self.k.season_peak_day_north - 182.5 };
        let ph = 2.0 * std::f64::consts::PI * (sc.doy - peak) / 365.0;
        c.temperature_mean_c + c.temperature_half_range_k * ph.cos()
    }

    /// Local occurrence: (lambda = -ln(1 - p_wet), convective share) from the
    /// non-orographic part of the T0 precipitation and the conditional-rate law.
    pub(crate) fn occurrence_local(&self, c: &ClimateSample, sc: &SolarClock, lat_deg: f64) -> (f64, f64) {
        let k = &self.k;
        let r_cc = (k.cond_rate_ref_mm_h * (1.0 + k.cond_rate_cc_per_k).powf(c.temperature_mean_c - k.cond_rate_ref_t_c))
            .clamp(k.cond_rate_min_mm_h, k.cond_rate_max_mm_h);
        let p_oro = c.precipitation_orographic_mm_per_yr.min(c.precipitation_mm_per_yr);
        let occ_rate = (c.precipitation_mm_per_yr - p_oro + k.orographic_frequency_share * p_oro).max(0.0) / k.hours_per_year;
        let p = (occ_rate / r_cc).min(k.wet_hour_fraction_max);
        let t_season = self.seasonal_t(c, sc, lat_deg);
        let share = k.conv_share_max * smoothstep(k.conv_share_t0_c, k.conv_share_t1_c, t_season);
        (-(1.0 - p).ln(), share)
    }

    /// Stratiform occurrence rate (1 - share) * lambda, averaged over the point
    /// and 8 points on a ring of `occurrence_smoothing_km` (weights 0.2 and 0.1).
    pub fn lam_strat_smooth(&self, lat_deg: f64, lon_deg: f64, t_h: f64) -> f64 {
        let k = &self.k;
        let d = k.occurrence_smoothing_km / k.planet_radius_km;
        let (p1, l1) = (lat_deg * crate::sun::DEG, lon_deg * crate::sun::DEG);
        let mut acc = 0.0;
        for i in 0..9 {
            let (la, lo, w) = if i == 0 {
                (lat_deg, lon_deg, 0.2)
            } else {
                let b = (i - 1) as f64 * std::f64::consts::PI / 4.0;
                let p2 = (p1.sin() * d.cos() + p1.cos() * d.sin() * b.cos()).clamp(-1.0, 1.0).asin();
                let l2 = l1 + (b.sin() * d.sin() * p1.cos()).atan2(d.cos() - p1.sin() * p2.sin());
                (p2 / crate::sun::DEG, l2 / crate::sun::DEG, 0.1)
            };
            let c = self.climate.sample(la, lo);
            let sc = solar_clock(t_h, la, lo, k.axial_tilt_deg);
            let (lam, share) = self.occurrence_local(&c, &sc, la);
            acc += w * (1.0 - share) * lam;
        }
        acc
    }

    /// Hour jitter: keyed standard normal per (T0 cell, integer hour, purpose),
    /// blended bilinearly over the 4 nearest cell centres (smooth in space).
    pub(crate) fn jitter(&self, lat_deg: f64, lon_deg: f64, hour: u64, what: &str) -> f64 {
        let w = self.climate.bilinear(lat_deg, lon_deg);
        let mut s = 0.0;
        let mut w2 = 0.0;
        for (cell, wt) in w {
            if wt == 0.0 { continue; }
            let scope = self.field.law_scope().address(&format!("cell/{cell}")).expect("address");
            s += wt * scope.draw(&format!("weather.hour.{what}"), hour, 0).expect("purpose").normal_f64();
            w2 += wt * wt;
        }
        // Rescale to unit variance (the blend of independent normals has
        // variance sum w^2), so lognormal factors keep mean 1 everywhere.
        if w2 > 0.0 { s / w2.sqrt() } else { 0.0 }
    }

    /// The weather state at (lat, lon) in degrees, at world hour `t_h`
    /// (f64 hours; integer = top of the hour), at `elevation_m` (None = the T0
    /// surface elevation there).
    pub fn state_at(&self, lat_deg: f64, lon_deg: f64, t_h: f64, elevation_m: Option<f64>) -> WeatherState {
        if let Some(ff) = &self.fronts {
            return self.state_at_v1(ff, lat_deg, lon_deg, t_h, elevation_m);
        }
        let k = &self.k;
        let c = self.climate.sample(lat_deg, lon_deg);
        let z = elevation_m.unwrap_or(c.t0_elevation_m);
        let sc = solar_clock(t_h, lat_deg, lon_deg, k.axial_tilt_deg);
        let hour = t_h.max(0.0).floor() as u64;
        let t_season = self.seasonal_t(&c, &sc, lat_deg);

        // --- precipitation climatology -> occurrence and rates
        // Synoptic occurrence uses a ring-smoothed climate (300 km), so the keep
        // threshold varies slowly and a system's rain area is not cut at a T0
        // coastline (keeping the test at the point keeps the void-probability
        // calibration exact); the local T0 precipitation, orographic part
        // included, sets the rate (seeder-feeder). Convective occurrence is local.
        let mean_rate = c.precipitation_mm_per_yr / k.hours_per_year; // mm/h, long-run mean
        let (lam_loc, share) = self.occurrence_local(&c, &sc, lat_deg);
        let lam_s = self.lam_strat_smooth(lat_deg, lon_deg, t_h);
        let theta_s = lam_s / self.strat.mean_cover;
        let conv_diurnal = 1.0 + k.conv_diurnal_amp * (2.0 * std::f64::consts::PI * (sc.solar_hour - k.conv_peak_solar_hour) / 24.0).cos();
        let lam_c_mean = share * lam_loc;
        let lam_c = lam_c_mean * conv_diurnal;
        let theta_c = lam_c / self.conv.mean_cover;
        // E[sum of kept marks] = lambda (Campbell; marks have mean 1), so with
        // r_c = ratio * r_s the expected rate is r_s (lam_s + ratio lam_c) = mean_rate.
        let denom = lam_s + k.conv_rate_ratio * lam_c_mean;
        let r_s = if denom > 0.0 { mean_rate / denom } else { 0.0 };
        let r_c = k.conv_rate_ratio * r_s;
        let p_wet = 1.0 - (-(lam_s + lam_c_mean)).exp();
        let r_cond = if p_wet > 0.0 { mean_rate / p_wet } else { 0.0 };

        // --- synoptic discs
        let hs = self.field.hits(&self.strat, lat_deg, lon_deg, t_h);
        let mut sum_s = 0.0;
        let mut n_s = 0u32;
        let mut cloud_s: f64 = 0.0;
        let mut activity: f64 = 0.0;
        let mut storm_en = [0.0f64, 0.0];
        let mut t_anom = 0.0;
        for h in &hs {
            let env = (std::f64::consts::PI * h.tau).sin();
            let halo_prof = 1.0 - smoothstep(1.0, k.strat_halo, h.d_frac);
            activity = activity.max(halo_prof * env.max(0.25));
            let wind_prof = (std::f64::consts::PI * h.d_frac / k.strat_halo).sin().max(0.0) * env;
            let v = k.strat_storm_wind_m_s * h.intensity.sqrt() * wind_prof;
            storm_en[0] += v * h.cyclonic_en[0];
            storm_en[1] += v * h.cyclonic_en[1];
            // Tropical air masses barely differ: the swing ramps in from 10 to 40 deg.
            let lat_w = 0.3 + 0.7 * smoothstep(10.0, 40.0, lat_deg.abs());
            t_anom += k.strat_temp_anomaly_k * lat_w * halo_prof * (std::f64::consts::PI * h.tau).cos();
            let kept = h.strength < theta_s;
            // Cloud lingers: ramp over the first and last 15 % of the life.
            let life = smoothstep(0.0, 0.15, h.tau) * (1.0 - smoothstep(0.85, 1.0, h.tau));
            let cl = if kept { halo_prof.max(0.0) } else { k.dry_front_cloud * halo_prof };
            cloud_s = cloud_s.max(cl * life.max(0.3));
            if kept && h.d_frac < 1.0 {
                let g = 2.0 * (1.0 - h.d_frac * h.d_frac);
                let e = 0.5 * std::f64::consts::PI * env;
                sum_s += h.intensity * g * e;
                n_s += 1;
            }
        }
        let t_anom = t_anom.clamp(-2.0 * k.strat_temp_anomaly_k, 2.0 * k.strat_temp_anomaly_k);

        // --- convective discs
        let hc = if lam_c > 0.0 { self.field.hits(&self.conv, lat_deg, lon_deg, t_h) } else { Vec::new() };
        let mut sum_c = 0.0;
        let mut n_c = 0u32;
        let mut cloud_c: f64 = 0.0;
        let mut gust = 0.0;
        for h in &hc {
            if h.strength >= theta_c { continue; }
            let prof = 1.0 - smoothstep(1.0, k.conv_halo, h.d_frac);
            cloud_c = cloud_c.max(0.75 * prof);
            if h.d_frac < 1.0 {
                let g = 2.0 * (1.0 - h.d_frac * h.d_frac);
                let e = 0.5 * std::f64::consts::PI * (std::f64::consts::PI * h.tau).sin();
                sum_c += h.intensity * g * e;
                n_c += 1;
                gust = k.conv_gust_m_s * h.intensity.sqrt();
            }
        }
        let jit = (k.hourly_jitter_sigma * self.jitter(lat_deg, lon_deg, hour, "rate") - 0.5 * k.hourly_jitter_sigma * k.hourly_jitter_sigma).exp();
        let rate_s = if n_s > 0 { r_s * sum_s } else { 0.0 };
        let rate_c = if n_c > 0 { r_c * sum_c } else { 0.0 };
        let rate = (rate_s + rate_c) * jit;
        let convective = rate_c > rate_s;

        // --- fair-weather cloud (daytime cumulus in humid air), keyed jitter
        let rh_mean = k.rh_mean_floor + k.rh_mean_span * c.aridity_index / (c.aridity_index + k.rh_mean_half_ai);
        let oro_share = if c.precipitation_mm_per_yr > 0.0 { (c.precipitation_orographic_mm_per_yr / c.precipitation_mm_per_yr).clamp(0.0, 1.0) } else { 0.0 };
        let day = smoothstep(-2.0, 10.0, sc.elevation_deg);
        let cj = self.jitter(lat_deg, lon_deg, hour, "cloud");
        let fair = ((rh_mean - 0.35) * 0.8 * day + 0.25 * oro_share + 0.08 * cj).clamp(0.0, 0.6);
        let cloud = 1.0 - (1.0 - cloud_s) * (1.0 - cloud_c) * (1.0 - fair);
        let cloud = if rate > 0.0 { cloud.max(0.9) } else { cloud }.clamp(0.0, 1.0);

        // --- temperature
        let t_cell_day = t_season + t_anom - k.lapse_rate_k_per_m * (z - c.t0_elevation_m);
        let amp = (k.diurnal_half_amp_base_k + k.diurnal_half_amp_arid_k * 0.5 / (c.aridity_index + 0.5)) * (1.0 - k.diurnal_cloud_damping * cloud);
        let (sr, hm) = (sc.sunrise_h.clamp(0.5, 11.0), k.diurnal_max_solar_hour);
        let hsol = sc.solar_hour;
        let shape = if hsol >= sr && hsol <= hm {
            -((std::f64::consts::PI * (hsol - sr) / (hm - sr)).cos())
        } else {
            let hh = if hsol > hm { hsol } else { hsol + 24.0 };
            (std::f64::consts::PI * (hh - hm) / (sr + 24.0 - hm)).cos()
        };
        let t_air = t_cell_day + amp * shape;

        // --- humidity: vapour pressure fixed through the day at the daily-mean RH
        let e = rh_mean * es_hpa(t_cell_day);
        let mut rh = (e / es_hpa(t_air)).min(1.0);
        if rate > 0.0 { rh = rh.max(k.rh_in_precip); }
        let td = dewpoint_c(rh * es_hpa(t_air));

        // --- wind
        let u0 = c.wind_u_m_s * (k.wind_quiet_factor + (1.0 - k.wind_quiet_factor) * activity);
        let v0 = c.wind_v_m_s * (k.wind_quiet_factor + (1.0 - k.wind_quiet_factor) * activity);
        let sun_up = sc.elevation_deg > 0.0;
        let clear = 1.0 - cloud;
        let mix = if sun_up {
            1.0 + k.wind_day_mixing * clear * (sc.elevation_deg * crate::sun::DEG).sin()
        } else {
            1.0 - k.wind_night_decoupling * clear
        };
        let ju = k.wind_jitter_sigma_m_s * self.jitter(lat_deg, lon_deg, hour, "wind-u");
        let jv = k.wind_jitter_sigma_m_s * self.jitter(lat_deg, lon_deg, hour, "wind-v");
        let mut wu = (u0 + storm_en[0]) * mix + ju * mix;
        let mut wv = (v0 + storm_en[1]) * mix + jv * mix;
        let mut speed = (wu * wu + wv * wv).sqrt();
        if gust > 0.0 && speed > 0.0 {
            let f = (speed + gust) / speed;
            wu *= f;
            wv *= f;
            speed += gust;
        }

        // --- precipitation type
        let precip_type = if rate <= 0.0 {
            PrecipType::None
        } else if t_air <= k.snow_max_t_c {
            PrecipType::Snow
        } else if t_air < k.rain_min_t_c {
            PrecipType::Mixed
        } else if rate < k.drizzle_max_mm_h && !convective {
            PrecipType::Drizzle
        } else {
            PrecipType::Rain
        };

        // --- radiative cooling and open-ground fog
        let night = 1.0 - smoothstep(-6.0, 3.0, sc.elevation_deg);
        let calm = 1.0 - smoothstep(k.calm_full_m_s, k.calm_none_m_s, speed);
        let cooling = night * clear.powf(k.cooling_cloud_exponent) * calm * if rate > 0.0 { 0.0 } else { 1.0 };
        // Open-ground fog may persist into the first hours after sunrise.
        let fog_time = 1.0 - smoothstep(3.0, 12.0, sc.elevation_deg);
        let flat_fog = fog_time * clear.powf(k.cooling_cloud_exponent) * calm * smoothstep(k.flat_fog_rh_onset, 1.0, rh)
            * if rate > 0.0 { 0.0 } else { 1.0 };

        // --- extinction
        let zk = z.max(0.0) / 1000.0;
        let ln = k.koschmieder_ln;
        let rayleigh = k.rayleigh_beta_per_km * (-zk / k.rayleigh_scale_height_km).exp();
        let humid = c.aridity_index / (c.aridity_index + 0.5);
        let v_dry = k.dry_visibility_arid_km + (k.dry_visibility_humid_km - k.dry_visibility_arid_km) * humid;
        let rh_g = rh.min(k.haze_rh_cap);
        let growth = if rh_g > k.haze_rh_ref { ((1.0 - rh_g) / (1.0 - k.haze_rh_ref)).powf(-k.haze_growth_gamma) } else { 1.0 };
        let haze = ln / v_dry * growth * (-zk / k.aerosol_scale_height_km).exp().max(k.aerosol_floor);
        let conc = |r: f64, v: f64| r / (3.6 * v);
        let beta_rain = |r: f64, v: f64| k.rain_beta_coef * conc(r, v).powf(k.rain_beta_exp);
        let beta_snow = |r: f64| k.snow_beta_coef * conc(r, k.snow_fall_speed_m_s).powf(k.snow_beta_exp);
        let precip = match precip_type {
            PrecipType::None => 0.0,
            PrecipType::Drizzle => beta_rain(rate, k.drizzle_fall_speed_m_s),
            PrecipType::Rain => beta_rain(rate, k.rain_fall_speed_m_s),
            PrecipType::Mixed => 0.5 * (beta_rain(rate, k.rain_fall_speed_m_s) + beta_snow(rate)),
            PrecipType::Snow => beta_snow(rate),
        };
        let dryness = (1.0 - c.aridity_index / k.dust_aridity_max).max(0.0);
        let dust = if rate > 0.0 { 0.0 } else { k.dust_coef * dryness * (speed - k.dust_threshold_m_s).max(0.0).powi(3) };
        let fog = flat_fog * ln / k.flat_fog_visibility_km;
        let ext = Extinction { rayleigh, haze, precip, dust, fog };
        let visibility_km = ln / (ext.total() - ext.fog);
        let visibility_open_ground_km = ln / ext.total();

        WeatherState {
            t_h,
            lat_deg,
            lon_deg,
            elevation_m: z,
            t0_cell: self.climate.cell_of(lat_deg, lon_deg),
            solar_hour: sc.solar_hour,
            sun_elevation_deg: sc.elevation_deg,
            cloud_cover: cloud,
            precip_type,
            precip_rate_mm_h: rate,
            convective: rate > 0.0 && convective,
            wind_u_m_s: wu,
            wind_v_m_s: wv,
            wind_speed_m_s: speed,
            air_temperature_c: t_air,
            dewpoint_c: td,
            relative_humidity: rh,
            visibility_km,
            visibility_open_ground_km,
            extinction_per_km: ext,
            radiative_cooling: cooling,
            flat_fog,
            synoptic_activity: activity,
            wet_hour_fraction: p_wet,
            conditional_rate_mm_h: r_cond,
        }
    }

    /// Strat hits at a point (exposed for diagnostics and tests).
    pub fn strat_hits(&self, lat_deg: f64, lon_deg: f64, t_h: f64) -> Vec<Hit> {
        self.field.hits(&self.strat, lat_deg, lon_deg, t_h)
    }
}
