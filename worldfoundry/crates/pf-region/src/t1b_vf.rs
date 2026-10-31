//! Lane VALLEY-FLOOR-FIX (2026-10-20): opt-in valley-floor profile and tributary fans for T1b RIVERNET.
//! LABEL: EXPLORATORY, every constant a project assumption awaiting review. Absent (`params.vf` None) = byte-identical.
//!
//! Why (lane VALLEY-FLOOR diagnosis, `unpublished-design-record`): nearest_channel_floor_v1 gives every
//! floor post the bed of its nearest channel post + max(0.5 %, bed gradient) x distance, a straight V falling TOWARD the
//! channel from both walls, with a hard kink at the wall foot that pf-detail's bicubic (B1) undershoots into a 1.5-3 m moat.
//! Real floodplains do the opposite near the river (a natural levee falling AWAY from the channel to a backswamp) and meet
//! the valley side in a concave colluvial footslope; tributaries leave their valleys onto the trunk floor in fans.
//!
//! Three switches (keys of pf-region-relief: `vf_floor=1`, `vf_foot=0|1`, `vf_fan=1`):
//! * floor (levee_backswamp_floor_v1): see `FLOOR_RULE`.
//! * foot (colluvial_footslope_v2, needs floor): see `FOOT_RULE`.
//! * fan (tributary_fan_v1): see `FAN_RULE`.
//!
//! The keyed draws never see these params (`RiverNetRecipeV1::draw_digest` strips them), so key on and key off share the
//! network, meanders, rock-mass and volcano draws; the E table is re-solved by the layer with the key (F6 sees the new floor).

use serde::{Deserialize, Serialize};

pub const FLOOR_RULE: &str = "levee_backswamp_floor_v1 (replaces nearest_channel_floor_v1's value, same W/2 band and nearest-post choice): a floor post at distance d from its nearest channel post (bed z_c, drainage area A, bed gradient S_b) takes z_c + max(s_x d + L(d), S_b min(d, g dx)); L(d) = H exp(-max(d - w/2, 0) / lambda) is the natural levee: Q = A r_ref / yr (A is the effective area Q / r_ref under the climate key), bankfull width w = 3.5 Q^0.5 and depth D = 0.27 Q^0.4 (the pf-detail B6 hydraulic geometry, DetailParamsV1 v1), crest H = min(f_l D, H_max) above the backswamp, decay length lambda = c_lambda w (overbank deposition thinning away from the bank, recalled form of Pizzuto 1987 / Adams et al. 2004; coefficients are project assumptions); s_x is a small backswamp drainage slope toward the river (backswamps drain to the river through tie channels and tributaries); the S_b min(d, g dx) guard keeps the posts within g posts of a channel above the channel's next post downstream (F1, as REVISION 1 of nearest_channel_floor_v1, now only near the channel instead of across the whole floor). Network channel posts inside a T1a floor (lane TC C1, ls_floor) take min(own bed, z_c + max(s_x d, S_b min(d, g dx))): no levee, so a tributary crossing the floor never climbs the levee";
pub const FOOT_RULE: &str = "colluvial_footslope_v2 (v1, a pre-eikonal parabola claimed as soft floor, was dropped: the claimed sources shortened every hillslope length L_H and moved hillslopes 1.5 km away by up to 90 m): on the finished hillslope surface (after crest rounding and the remap, before the fans, the bank rule and the fill), the concave foot of every wide floor's wall is filled by the depositional part of linear hillslope diffusion: z <- z + w max(0, G_sigma * z - z), G_sigma = a 1-D Gaussian of sigma = L_f / 4 along the across-corner direction (the gradient of the signed distance to the nearest edge; a tributary valley crossing the zone is smoothed along its axis, not filled across), bilinear samples of the free posts (channel and water posts carry no weight) of the surface limited to the rock's limiting slope from the edges (z_lim = min(z, z_edge + min(Sc, Sc_ref) r), Sc_ref = sc_ref 0.7, about the repose slope of colluvium; a cliff feeds the fill only as a wall at Sc would, so a cliff foot gets the apron of a critical-slope wall and the cliff face above it is untouched), so the gradient rises from the floor slope to the wall's own gradient over about L_f (erf profile, 2 % to 98 % of the change over 4.1 sigma), C1 into both, and planar or convex ground (the wall face, cliffs, crests) is untouched; edges = floor posts of a floor of half-width W/2 >= L_f (narrower floors would be filled wall to wall) with a free non-floor 8-neighbour, standing at least 1 m above sea level (a drowned wall foot belongs to the coast stage); zone = free posts within L_f of an edge, weight w = clamp(2 (1 - r / L_f), 0, 1) (r = distance to the nearest edge); mass is added (colluvium shed from the wall above; the wall itself is not lowered)";
pub const FAN_RULE: &str = "tributary_fan_v1: for every channel curve, walking downstream, the first curve point whose post lies in the floor band of a channel post with drainage area >= r_a x the point's own area is the fan apex (one per curve); H = the curve's own bed there - the larger channel's floor value there (no fan if H < H_min); fan slope s_f = clamp(c_f (A / 1 km2)^e_f, s_min, s_max) (recalled form of the Bull 1964 fan-gradient / drainage-area relation, coefficients project assumptions); supply V = E(apex) x A x T (lane CO bayhead_fan_v1 supply law, T = its fan_time 7 ky); R = min(R_mass, R_geom, R_trim): R_mass from lane CO's half-disc law R = sqrt(2 V / (pi h)) with h = the mean thickness s_f R / 3 of a cone, i.e. R = (6 V / (pi s_f))^(1/3); R_geom = H / s_f (the fan fills the slope break and no more); R_trim = the apex distance to the larger channel (the trunk river trims the fan toe); apex height z_a = floor + s_f R; after crest rounding and the remap, before the bank rule and the fill, every post within R except T1a channel and water posts is raised to max(z, z_a - s_f r) (deposit only raises; tributary network channel posts on the fan ride its surface, which falls away from the apex); raised posts by >= 0.1 m expose the fan deposit (over class OVER_FAN, lithology 10)";
pub const DELTA_RULE: &str = "delta_surface_v1 (part of vf_floor, coast on only): after lane CO's coast stage, land posts its bayhead fan / fan delta covers (over OVER_FAN, not tributary fans) take the levee L(d) of their nearest channel-curve post (rasterized channel or estuary arm, its drainage area; within 4 decay lengths + w/2, at most 600 m) on top of the fan surface, and (vf_foot) the colluvial_footslope_v2 fill at the fan's walls (edges = fan posts with a free non-fan neighbour); then the priority-flood fill again with the fixed posts as outlets";
pub const LABEL: &str = "EXPLORATORY lane VALLEY-FLOOR-FIX 2026-10-20; every constant a project assumption awaiting review";

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VfParamsV1 {
    /// levee_backswamp_floor_v1 on (1) / off (0: nearest_channel_floor_v1 values).
    pub floor: u32,
    /// colluvial_footslope_v2 on (needs floor = 1).
    pub foot: u32,
    /// tributary_fan_v1 on.
    pub fan: u32,
    /// 1 = window-only (key vf_window=1): the window applies vf on the E table and network of a layer solved WITHOUT vf (the
    /// layer's level solve reroutes and re-solves E from the surface, so a vf layer changes hillslopes far from any floor);
    /// the layer refuses it. 0 = the layer is solved with vf too.
    #[serde(default, skip_serializing_if = "crate::t1b::is_zero_u32_pub")]
    pub window_only: u32,
    pub runoff_ref_mm_per_yr: i64,
    pub hg_width_coeff_permille: i64,
    pub hg_width_exp_permille: i64,
    pub hg_depth_coeff_permille: i64,
    pub hg_depth_exp_permille: i64,
    pub levee_depth_permille: i64,
    pub levee_decay_widths_permille: i64,
    pub levee_max_mm: i64,
    pub backswamp_slope_ppm: i64,
    pub guard_posts: u32,
    pub foot_length_mm: i64,
    pub fan_area_ratio_permille: i64,
    pub fan_slope_coeff_ppm: i64,
    pub fan_slope_exp_permille: i64,
    pub fan_slope_min_ppm: i64,
    pub fan_slope_max_ppm: i64,
    pub fan_time_yr: i64,
    pub fan_height_min_mm: i64,
    pub floor_rule: String,
    pub foot_rule: String,
    pub fan_rule: String,
    #[serde(default)]
    pub delta_rule: String,
    pub label: String,
}

impl VfParamsV1 {
    /// v1 constants with every switch off (keys switch them on).
    pub fn off() -> Self {
        VfParamsV1 {
            floor: 0, foot: 0, fan: 0, window_only: 0,
            runoff_ref_mm_per_yr: 1000,
            hg_width_coeff_permille: 3500, hg_width_exp_permille: 500, hg_depth_coeff_permille: 270, hg_depth_exp_permille: 400,
            levee_depth_permille: 400, levee_decay_widths_permille: 1500, levee_max_mm: 2000,
            backswamp_slope_ppm: 500, guard_posts: 2,
            foot_length_mm: 150_000,
            fan_area_ratio_permille: 3000, fan_slope_coeff_ppm: 50_000, fan_slope_exp_permille: -150, fan_slope_min_ppm: 10_000, fan_slope_max_ppm: 120_000,
            fan_time_yr: 7_000, fan_height_min_mm: 500,
            floor_rule: FLOOR_RULE.into(), foot_rule: FOOT_RULE.into(), fan_rule: FAN_RULE.into(), delta_rule: DELTA_RULE.into(), label: LABEL.into(),
        }
    }
    pub fn validate(&self) -> Result<(), String> {
        if self.floor > 1 || self.foot > 1 || self.fan > 1 || self.window_only > 1 { return Err("vf switches are 0 or 1".into()); }
        if self.foot == 1 && self.floor != 1 { return Err("vf_foot needs vf_floor=1".into()); }
        if self.floor == 0 && self.fan == 0 { return Err("vf: no switch on (vf_floor=1 and/or vf_fan=1)".into()); }
        if !(50_000..=200_000).contains(&self.foot_length_mm) { return Err("vf foot length must be 50-200 m".into()); }
        if self.levee_max_mm < 0 || self.levee_decay_widths_permille <= 0 || self.fan_slope_min_ppm <= 0 || self.fan_slope_max_ppm < self.fan_slope_min_ppm { return Err("vf constants".into()); }
        Ok(())
    }
    pub fn floor_on(&self) -> bool { self.floor == 1 }
    pub fn foot_on(&self) -> bool { self.foot == 1 }
    pub fn fan_on(&self) -> bool { self.fan == 1 }
    fn q(&self, area_m2: f64) -> f64 { area_m2.max(0.0) * self.runoff_ref_mm_per_yr as f64 / 1000.0 / crate::SECONDS_PER_YEAR }
    /// Bankfull width and depth (m) of a channel of drainage area A (pf-detail B6 hydraulic geometry).
    pub fn hydraulic(&self, area_m2: f64) -> (f64, f64) {
        let q = self.q(area_m2);
        (self.hg_width_coeff_permille as f64 / 1000.0 * q.powf(self.hg_width_exp_permille as f64 / 1000.0),
         self.hg_depth_coeff_permille as f64 / 1000.0 * q.powf(self.hg_depth_exp_permille as f64 / 1000.0))
    }
    /// Natural levee height (m) at distance d (m) from the channel centreline.
    pub fn levee(&self, area_m2: f64, d: f64) -> f64 {
        let (w, dep) = self.hydraulic(area_m2);
        let h = (self.levee_depth_permille as f64 / 1000.0 * dep).min(self.levee_max_mm as f64 / 1000.0);
        let lam = self.levee_decay_widths_permille as f64 / 1000.0 * w;
        if h <= 0.0 || lam <= 1e-6 { return 0.0; }
        h * (-(d - 0.5 * w).max(0.0) / lam).exp()
    }
    fn guard(&self, s_b: f64, d: f64, dx: f64) -> f64 { s_b.max(0.0) * d.min(self.guard_posts as f64 * dx) }
    /// Floor height above the nearest channel post's bed at distance d (levee + backswamp + F1 guard).
    pub fn floor_rel(&self, area_m2: f64, s_b: f64, d: f64, dx: f64) -> f64 {
        (self.backswamp_slope_ppm as f64 / 1e6 * d + self.levee(area_m2, d)).max(self.guard(s_b, d, dx))
    }
    /// Lane TC C1 value with the floor rule on: no levee (a tributary does not climb the levee).
    pub fn tc_rel(&self, s_b: f64, d: f64, dx: f64) -> f64 { (self.backswamp_slope_ppm as f64 / 1e6 * d).max(self.guard(s_b, d, dx)) }
    pub fn fan_slope(&self, area_m2: f64) -> f64 {
        let a = (area_m2 / 1e6).max(1e-6);
        (self.fan_slope_coeff_ppm as f64 / 1e6 * a.powf(self.fan_slope_exp_permille as f64 / 1000.0))
            .clamp(self.fan_slope_min_ppm as f64 / 1e6, self.fan_slope_max_ppm as f64 / 1e6)
    }
}

/// Stats (published in the window header's surface stats only when vf is on).
#[derive(Clone, Debug, Default, Serialize)]
pub struct VfStats {
    pub floor_posts_levee_gt_0p1m: u32,
    pub floor_mean_change_m: f64,
    pub foot_edge_posts: u32,
    pub foot_raised_posts: u32,
    pub foot_raise_max_m: f64,
    pub foot_volume_m3: f64,
    pub fans: Vec<FanRec>,
    pub fan_raised_posts: u32,
    pub fan_volume_m3: f64,
    /// delta_surface_v1 (coast on): levee and footslope on lane CO's bayhead fan surface, and the refill.
    pub delta_levee_posts: u32,
    pub delta_levee_volume_m3: f64,
    pub delta_foot_raised_posts: u32,
    pub delta_foot_volume_m3: f64,
    pub delta_fill_raised_posts: u32,
    pub delta_fill_raised_max_m: f64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct FanRec { pub x_m: f64, pub y_m: f64, pub area_m2: f64, pub step_m: f64, pub slope: f64, pub supply_m3: f64, pub r_mass_m: f64, pub r_geom_m: f64, pub r_trim_m: f64, pub radius_m: f64, pub apex_m: f64 }

/// colluvial_footslope_v2 on a square lattice of n x n posts, spacing dx (see `FOOT_RULE`). `free[i]` = fixed FREE; `floor`
/// the floor posts; `wide[i]` = floor post of a floor with half-width >= L_f; z is modified on free posts only.
pub fn foot_fill(n: usize, dx: f64, free: &[bool], floor: &[bool], wide: &[bool], z: &mut [f64], sc: &dyn Fn(usize, f64) -> f64, sc_cap: f64, z_min: f64, lf: f64, st: &mut VfStats) {
    let nn = n * n;
    const D8: [(i64, i64); 8] = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)];
    let mut edges: Vec<usize> = Vec::new();
    for i in 0..nn {
        if !floor[i] || !wide[i] || !free[i] || !z[i].is_finite() || z[i] < z_min { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        if D8.iter().any(|&(ox, oy)| {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { return false; }
            let j = yy as usize * n + xx as usize;
            free[j] && !floor[j]
        }) { edges.push(i); }
    }
    st.foot_edge_posts = edges.len() as u32;
    if edges.is_empty() { return; }
    // distance to the nearest edge (brute force within L_f)
    let r = (lf / dx).ceil() as i64;
    let mut dist = vec![f64::INFINITY; nn];
    let scv: Vec<f64> = (0..nn).map(|i| if z[i].is_finite() { sc(i, z[i]).min(sc_cap) } else { sc_cap }).collect();
    // the wall limited to the rock's limiting slope from the edges (cliffs do not feed the colluvial fill)
    let mut lim = vec![f64::INFINITY; nn];
    for &e in &edges {
        let (x, y) = ((e % n) as i64, (e / n) as i64);
        for oy in -r..=r { for ox in -r..=r {
            let (xx, yy) = (x + ox, y + oy);
            if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let d = ((ox * ox + oy * oy) as f64).sqrt() * dx;
            let k = yy as usize * n + xx as usize;
            if d < dist[k] { dist[k] = d; }
            let c = z[e] + scv[k] * d;
            if c < lim[k] { lim[k] = c; }
        } }
    }
    // signed distance (+ hillslope side, - floor side) and its gradient = the across-corner direction
    let sd: Vec<f64> = (0..nn).map(|i| if floor[i] { -dist[i] } else { dist[i] }).collect();
    let zl: Vec<f64> = (0..nn).map(|i| if free[i] && z[i].is_finite() { z[i].min(lim[i]) } else { f64::NAN }).collect();
    // bilinear sample of zl at lattice coordinates (fx, fy); NaN posts carry no weight
    let samp = |fx: f64, fy: f64| -> Option<f64> {
        let (x0, y0) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = (fx - x0 as f64, fy - y0 as f64);
        let (mut a, mut b) = (0.0, 0.0);
        for (ox, oy, w) in [(0i64, 0i64, (1.0 - tx) * (1.0 - ty)), (1, 0, tx * (1.0 - ty)), (0, 1, (1.0 - tx) * ty), (1, 1, tx * ty)] {
            let (xx, yy) = (x0 + ox, y0 + oy);
            if w <= 0.0 || xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 { continue; }
            let v = zl[yy as usize * n + xx as usize];
            if v.is_finite() { a += w * v; b += w; }
        }
        if b > 1e-9 { Some(a / b) } else { None }
    };
    // 1-D Gaussian (sigma = L_f / 4, radius 3 sigma, steps of dx) along the across-corner direction
    let sig = lf / 4.0;
    let kr = (3.0 * sig / dx).ceil() as i64;
    let z0: Vec<f64> = z.to_vec();
    for i in 0..nn {
        if dist[i] > lf || !free[i] || !z0[i].is_finite() { continue; }
        let (x, y) = ((i % n) as i64, (i / n) as i64);
        if x == 0 || y == 0 || x == n as i64 - 1 || y == n as i64 - 1 { continue; }
        let gx = sd[i + 1] - sd[i - 1];
        let gy = sd[i + n] - sd[i - n];
        let gl = gx.hypot(gy);
        if !(gl > 1e-9) || !gl.is_finite() { continue; }
        let (ux, uy) = (gx / gl, gy / gl);
        let (mut a, mut b) = (0.0, 0.0);
        for o in -kr..=kr {
            let t = o as f64;
            let w = (-0.5 * (t * dx / sig).powi(2)).exp();
            if let Some(v) = samp(x as f64 + t * ux, y as f64 + t * uy) { a += w * v; b += w; }
        }
        if b <= 1e-12 { continue; }
        let g = a / b;
        let w = (2.0 * (1.0 - dist[i] / lf)).clamp(0.0, 1.0);
        let dz = w * (g - z0[i]).max(0.0);
        if dz > 0.0 { z[i] = z0[i] + dz; st.foot_raised_posts += 1; st.foot_raise_max_m = st.foot_raise_max_m.max(dz); st.foot_volume_m3 += dz * dx * dx; }
    }
}

/// Fan apexes from the channel curves: (post index, curve point) of the first point (downstream) inside a larger channel's
/// floor band. `tr_area[i]` = the largest drainage area of a channel post whose floor band covers post i.
pub fn fan_apexes(curves: &[crate::t1b_rivernet::Curve], idx: &dyn Fn(f64, f64) -> Option<usize>, tr_area: &[f64], ratio: f64) -> Vec<(usize, [f64; 4])> {
    let mut out = Vec::new();
    for c in curves {
        if c.pts.len() < 2 { continue; }
        // walk downstream (drainage area grows downstream)
        let fwd = c.pts[c.pts.len() - 1][3] >= c.pts[0][3];
        let n = c.pts.len();
        for k in 0..n {
            let q = if fwd { c.pts[k] } else { c.pts[n - 1 - k] };
            let Some(i) = idx(q[0], q[1]) else { continue };
            if q[3] > 0.0 && tr_area[i] >= ratio * q[3] { out.push((i, q)); break; }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn levee_falls_away_from_the_channel_and_is_bounded() {
        let p = VfParamsV1 { floor: 1, ..VfParamsV1::off() };
        let a = 7.385e9; // the S2j trunk (lane VALLEY-FLOOR): Q ~ 234 m3/s
        let (w, d) = p.hydraulic(a);
        assert!((w - 53.5).abs() < 1.0 && (d - 2.39).abs() < 0.05, "w {w} d {d}");
        let l32 = p.levee(a, 32.0);
        let l600 = p.levee(a, 600.0);
        assert!(l32 > 0.5 && l32 < 1.0 && l600 < 0.01, "{l32} {l600}");
        // floor near the channel above the F1 guard and falling away from it
        let s_b = 0.001;
        assert!(p.floor_rel(a, s_b, 32.0, 32.0) > p.floor_rel(a, s_b, 160.0, 32.0));
        assert!(p.floor_rel(a, s_b, 32.0, 32.0) >= s_b * 32.0);
        assert!(p.levee(1e14, 0.0) <= 2.0 + 1e-12);
    }
    #[test]
    fn fan_slope_bounds() {
        let p = VfParamsV1 { fan: 1, ..VfParamsV1::off() };
        assert!((p.fan_slope(1e6) - 0.05).abs() < 1e-12);
        assert!(p.fan_slope(1e3) <= 0.12 && p.fan_slope(1e13) >= 0.01);
    }
    #[test]
    fn foot_fill_fills_the_corner_and_leaves_planes() {
        // 1-D strip: floor x < 20 (z = 0), wall z = 0.3 (x - 19) dx beyond
        let (n, dx, lf) = (40usize, 32.0, 150.0);
        let nn = n * n;
        let free = vec![true; nn];
        let floor: Vec<bool> = (0..nn).map(|i| i % n < 20).collect();
        let wide = floor.clone();
        let mut z: Vec<f64> = (0..nn).map(|i| { let x = (i % n) as f64; if x < 20.0 { 0.0 } else { 0.3 * (x - 19.0) * dx } }).collect();
        let z0 = z.clone();
        let mut st = VfStats::default();
        foot_fill(n, dx, &free, &floor, &wide, &mut z, &|_, _| 0.7, 0.7, f64::NEG_INFINITY, lf, &mut st);
        let row = 20 * n;
        assert!(z[row + 19] > z0[row + 19] + 1.0, "corner filled: {}", z[row + 19]);
        assert!((z[row + 5] - z0[row + 5]).abs() < 1e-9 && (z[row + 35] - z0[row + 35]).abs() < 1e-6, "far floor and wall face untouched");
        for x in 1..n { assert!(z[row + x] >= z[row + x - 1] - 1e-9); }
        // slope rises monotonically across the footslope
        let g: Vec<f64> = (1..n).map(|x| z[row + x] - z[row + x - 1]).collect();
        for k in 12..24 { assert!(g[k] >= g[k - 1] - 1e-6, "k {k}: {:?}", &g[10..26]); }
        assert!(st.foot_raised_posts > 0 && st.foot_volume_m3 > 0.0);
    }
}
