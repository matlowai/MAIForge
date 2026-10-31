//! WEATHER-V1 hourly state (lane WEATHER-V1, 2026-10-21; opt-in `fronts=1`).
//!
//! Same law as `state::state_at` (WEATHER-V0) for climatology, occurrence and
//! rate calibration, temperature, humidity, wind and extinction; only the
//! precipitation geometry changes: frontal band families and convective
//! clusters (`fronts`) instead of discs, plus the optional orographic
//! modulation of frontal rain (`oro`). The V0 function is untouched; this is a
//! separate code path selected by `WeatherModel::fronts`.

use crate::fronts::FrontField;
use crate::state::{dewpoint_c, es_hpa, smoothstep, Extinction, PrecipType, WeatherModel, WeatherState};
use crate::sun::{solar_clock, DEG};

impl<'a> WeatherModel<'a> {
    /// T0 elevation gradient (dz/dx east, dz/dy north; m/m) at a point,
    /// central differences over +-`step_cells` T0 cells of the bilinear surface.
    pub fn t0_gradient(&self, lat_deg: f64, lon_deg: f64, step_cells: f64) -> [f64; 2] {
        let dla = step_cells * self.climate.dlat_deg();
        let dlo = step_cells * self.climate.dlon_deg();
        let z = |la: f64, lo: f64| self.climate.sample(la, lo).t0_elevation_m;
        let r_m = self.k.planet_radius_km * 1000.0;
        let dx = 2.0 * dlo * DEG * r_m * (lat_deg * DEG).cos().max(0.05);
        let dy = 2.0 * dla * DEG * r_m;
        [(z(lat_deg, lon_deg + dlo) - z(lat_deg, lon_deg - dlo)) / dx, (z(lat_deg + dla, lon_deg) - z(lat_deg - dla, lon_deg)) / dy]
    }

    pub(crate) fn state_at_v1(&self, ff: &FrontField, lat_deg: f64, lon_deg: f64, t_h: f64, elevation_m: Option<f64>) -> WeatherState {
        let k = &self.k;
        let kf = &ff.k;
        let c = self.climate.sample(lat_deg, lon_deg);
        let z = elevation_m.unwrap_or(c.t0_elevation_m);
        let sc = solar_clock(t_h, lat_deg, lon_deg, k.axial_tilt_deg);
        let hour = t_h.max(0.0).floor() as u64;
        let t_season = self.seasonal_t(&c, &sc, lat_deg);

        // --- precipitation climatology (identical to V0)
        let mean_rate = c.precipitation_mm_per_yr / k.hours_per_year;
        let (lam_loc, share) = self.occurrence_local(&c, &sc, lat_deg);
        let lam_s = self.lam_strat_smooth(lat_deg, lon_deg, t_h);
        let theta_s = lam_s / kf.front_mean_cover;
        let conv_diurnal = 1.0 + k.conv_diurnal_amp * (2.0 * std::f64::consts::PI * (sc.solar_hour - k.conv_peak_solar_hour) / 24.0).cos();
        let lam_c_mean = share * lam_loc;
        let lam_c = lam_c_mean * conv_diurnal;
        let theta_c = lam_c / kf.cl_mean_cover;
        let denom = lam_s + k.conv_rate_ratio * lam_c_mean;
        let r_s = if denom > 0.0 { mean_rate / denom } else { 0.0 };
        let r_c = k.conv_rate_ratio * r_s;
        let p_wet = 1.0 - (-(lam_s + lam_c_mean)).exp();
        let r_cond = if p_wet > 0.0 { mean_rate / p_wet } else { 0.0 };

        // --- frontal band families and convective clusters
        let fa = ff.eval_fronts(lat_deg, lon_deg, t_h, theta_s);
        let ca = if lam_c > 0.0 { ff.eval_clusters(lat_deg, lon_deg, t_h, theta_c) } else { Default::default() };
        let cloud_s = fa.cloud;
        let cloud_c = ca.cloud;
        // Wind, synoptic activity and the warm/cold swing stay those of the V0
        // synoptic systems (the same keyed discs, used here for their wind
        // field only, not for rain or cloud): this keeps the V0 calm-night and
        // wind climate that COLDPOOL depends on. Measured alternative (lane
        // notes): band-tied activity and wind cut calm clear nights at site 1
        // from 442 to 34-287 h/yr.
        let hs = self.field.hits(&self.strat, lat_deg, lon_deg, t_h);
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
            let lat_w = 0.3 + 0.7 * smoothstep(10.0, 40.0, lat_deg.abs());
            t_anom += k.strat_temp_anomaly_k * lat_w * halo_prof * (std::f64::consts::PI * h.tau).cos();
        }
        let t_anom = t_anom.clamp(-2.0 * k.strat_temp_anomaly_k, 2.0 * k.strat_temp_anomaly_k);
        // Band cloud thins under high pressure (low V0 activity), so clear and
        // calm stay correlated as in V0 (rain still forces cloud >= 0.9 below).
        let cloud_s = cloud_s * (kf.cloud_activity_floor + (1.0 - kf.cloud_activity_floor) * activity);
        let gust = if ca.n_rain > 0 { k.conv_gust_m_s * ca.gust_intensity.sqrt() } else { 0.0 };
        let jit = (k.hourly_jitter_sigma * self.jitter(lat_deg, lon_deg, hour, "rate") - 0.5 * k.hourly_jitter_sigma * k.hourly_jitter_sigma).exp();
        let rate_s0 = if fa.n_rain > 0 { r_s * fa.sum_rain } else { 0.0 };
        let rate_c = if ca.n_rain > 0 { r_c * ca.sum_rain } else { 0.0 };
        let raining = rate_s0 + rate_c > 0.0;
        // Soft edges: V0 forces cloud >= 0.9 and RH >= 0.95 at any rate > 0,
        // which draws a hard visibility edge at every rain footprint; V1 ramps
        // both in with the (pre-orographic) rate.
        let wet_soft = smoothstep(0.0, kf.rain_soft_mm_h, (rate_s0 + rate_c) * jit);

        // --- fair-weather cloud (as V0)
        let rh_mean = k.rh_mean_floor + k.rh_mean_span * c.aridity_index / (c.aridity_index + k.rh_mean_half_ai);
        let oro_share = if c.precipitation_mm_per_yr > 0.0 { (c.precipitation_orographic_mm_per_yr / c.precipitation_mm_per_yr).clamp(0.0, 1.0) } else { 0.0 };
        let day = smoothstep(-2.0, 10.0, sc.elevation_deg);
        let cj = self.jitter(lat_deg, lon_deg, hour, "cloud");
        let fair = ((rh_mean - 0.35) * 0.8 * day + 0.25 * oro_share + 0.08 * cj).clamp(0.0, 0.6);
        let cloud = 1.0 - (1.0 - cloud_s) * (1.0 - cloud_c) * (1.0 - fair);
        let cloud = cloud.max(0.9 * wet_soft).clamp(0.0, 1.0);

        // --- temperature (as V0)
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

        // --- humidity (as V0)
        let e = rh_mean * es_hpa(t_cell_day);
        let mut rh = (e / es_hpa(t_air)).min(1.0);
        if raining { rh += (k.rh_in_precip - rh).max(0.0) * wet_soft; }
        let td = dewpoint_c(rh * es_hpa(t_air));

        // --- wind (as V0, frontal wind from the most active system)
        let u0 = c.wind_u_m_s * (k.wind_quiet_factor + (1.0 - k.wind_quiet_factor) * activity);
        let v0 = c.wind_v_m_s * (k.wind_quiet_factor + (1.0 - k.wind_quiet_factor) * activity);
        let sun_up = sc.elevation_deg > 0.0;
        let clear = 1.0 - cloud;
        let mix = if sun_up {
            1.0 + k.wind_day_mixing * clear * (sc.elevation_deg * DEG).sin()
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

        // --- orographic modulation of frontal rain (oro=1)
        let oro_f = if self.oro && rate_s0 > 0.0 {
            let g = self.t0_gradient(lat_deg, lon_deg, kf.oro_grad_step_cells);
            let w = ((storm_en[0] + ju) * g[0] + (storm_en[1] + jv) * g[1]) * mix;
            1.0 + kf.oro_amp * (w / kf.oro_w_scale_m_s).tanh()
        } else {
            1.0
        };
        let rate_s = rate_s0 * oro_f;
        let rate = (rate_s + rate_c) * jit;
        let convective = rate_c > rate_s;

        // --- precipitation type (as V0)
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

        // --- radiative cooling and open-ground fog (as V0)
        let night = 1.0 - smoothstep(-6.0, 3.0, sc.elevation_deg);
        let calm = 1.0 - smoothstep(k.calm_full_m_s, k.calm_none_m_s, speed);
        let cooling = night * clear.powf(k.cooling_cloud_exponent) * calm * if rate > 0.0 { 0.0 } else { 1.0 };
        let fog_time = 1.0 - smoothstep(3.0, 12.0, sc.elevation_deg);
        let flat_fog = fog_time * clear.powf(k.cooling_cloud_exponent) * calm * smoothstep(k.flat_fog_rh_onset, 1.0, rh)
            * if rate > 0.0 { 0.0 } else { 1.0 };

        // --- extinction (as V0)
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
}
