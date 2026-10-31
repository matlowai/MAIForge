//! Lane PYRO (2026-10-16, opt-in `pyro_sub=1`, needs beds=2 and bed_meta=1): pyroclastic subtype of the exposed bed where the
//! BS / DE closures record only facies 5 (weak pyroclastic) or 6 (indurated pyroclastic).
//!
//! design 2026-10-15/16: rock appearance attached to provenance; "geology-derived generation". Real stratovolcano and cinder-cone
//! sections show dark scoria bands, pale pumice / ash, massive ignimbrite sheets and lahar diamicts in proportions set by the
//! eruption style. The subtype is chosen from the VOLCANO RECORD (lane VE's per-volcano event counts in the T1b window header:
//! cinder cones, explosive eruptions by VEI, domes, lahars, maars; else a declared default by volcano kind), the magma family
//! (volcano kind), the distance of the post from the vent, and the unit's own geometry (local deposited thickness, position
//! inside the unit). Keyed draws (purpose `t2.closure.pyro_sub.v1`, a separate stream: no existing draw changes) pick which
//! recorded event class an unrecorded (keyed) unit belongs to, with the record's own frequencies, and the sub-resolution
//! sequence inside a unit (thickness of a co-ignimbrite / ash-cloud ash top, of a basal Plinian fall layer, of the non-welded
//! top and base zones of a welded ignimbrite). Opt-in `pyro_sub_hard=1` (needs pyro_sub=1) also sets the strength of the
//! pyroclastic parts of the lobe columns from the subtype (de_lobes); without it the subtype is presentation only.

use serde::Serialize;
use vz_keyed_rng::Scope;
use crate::R;

pub const PURPOSE_PYRO: &str = "t2.closure.pyro_sub.v1";

pub const SUB_SCORIA: u8 = 1;
pub const SUB_ASH: u8 = 2;
pub const SUB_PUMICE: u8 = 3;
pub const SUB_IGN_NW: u8 = 4;
pub const SUB_IGN_W: u8 = 5;
pub const SUB_LAHAR: u8 = 6;
pub const SUB_BAF: u8 = 7;
pub const SUB_UNRESOLVED: u8 = 254;
pub const SUB_NONE: u8 = 255;
/// Names by code (index 1-7); the first six are lane CLIFF2 facies_contrast_v1.json `subtypes_not_recorded` keys.
pub const SUB_NAMES: [&str; 8] = ["", "scoria_fall", "ash_fall", "pumice_fall", "ignimbrite_nonwelded", "ignimbrite_welded", "lahar", "block_and_ash"];

pub const PYRO_FIELD_NOTE: &str = "u8 lane PYRO pyro_sub=1: pyroclastic subtype of the exposed bed where bed_style_id is 5 or 6 (a declared closure from the volcano record, see pyro_sub_note): 1 scoria_fall (scoria / spatter lapilli: Strombolian to Vulcanian fall, cinder cones, spatter ramparts), 2 ash_fall (fine tephra beyond the lapilli limit, co-ignimbrite and ash-cloud ash tops, phreatomagmatic maar ash), 3 pumice_fall (Plinian pumice lapilli of evolved arc magma, VEI >= 4), 4 ignimbrite_nonwelded, 5 ignimbrite_welded (welded core of a thick ignimbrite of a VEI >= 5 eruption), 6 lahar (volcanic debris-flow / hyperconcentrated-flow deposit), 7 block_and_ash (small-volume pyroclastic density current: dome-collapse block-and-ash or VEI 3 scoria flow; no lane CLIFF2 entry), 254 pyroclastic but no volcano record to choose from (strata lava-stack package), 255 not a pyroclastic bed (other facies, deposit, water)";

pub const PYRO_SUB_NOTE: &str = "lane PYRO pyro_sub=1: the pyroclastic subtype of every facies 5 / 6 bed is a declared closure from the volcano record, not a simulated eruption history. Record: lane VE's event counts per volcano (T1b window header prep.net.volcanoes[k].ve.counts: cinder_cone = a Strombolian fall event of VEI 2; explosive_veiN = a fall event and a pyroclastic-density-current event of VEI N; dome = a dome-collapse block-and-ash event; lahar = a lahar; maar = a phreatomagmatic ash event), else a declared default by volcano kind (kind_default). Magma family by kind: continental / island arc = evolved (andesite-dacite; Plinian lapilli are pumice), hotspot / rift / LIP = basaltic (all lapilli are scoria). Rules: recorded lane VE units keep their kind (lahar paths = lahar; cinder cones and fissure ramparts = scoria_fall); a keyed pyroclastic-density-current / lahar lobe is a lahar or a PDC event with the record's frequencies (keyed draw); a PDC of VEI >= ign_vei_min is an ignimbrite (Plinian column collapse), smaller ones block_and_ash; a keyed tephra-fall mantle is one fall event with the record's VEI frequencies; fall grain size by distance r from the vent: lapilli (scoria, or pumice for evolved magma at VEI >= 4) within r_lap = r_lap_bb_mult x bb(VEI), ash beyond, with bb the thickness e-folding length of lane VE's apron fall rule (bb 2 / 4 / 8 km at VEI 4 / 5 / 6, declared 0.5 / 1 km at VEI 2 / 3); a pyroclastic package of the per-volcano 1D sequence (BS beds) is one event drawn from all event classes. Where the column rule extends a unit above its deposited top (DATUM2 de_seq: the youngest present keyed unit filled up to the recorded units / D), the extension takes the unit's body code. Sub-resolution sequence inside a unit (keyed): block_and_ash carries an ash-cloud ash top (lognormal median 0.3 m, at most 20 % of the local thickness); an ignimbrite carries a co-ignimbrite ash top (median 0.4 m, at most 15 %) and a basal Plinian fall layer (median 0.8 m, at most 30 %), and where the eruption is VEI >= weld_vei_min and the local deposited thickness is at least weld_min_m its body has a welded core between non-welded top (15-35 % of the body) and base (10-25 %) zones (cooling-unit zonation). Thickness, extent, paths, ages and hardness of the units are unchanged (presentation only) unless pyro_sub_hard=1. Numbers are project assumptions from recalled volcanology (walker1973, sparks1976, smith1960, pyle1989, vessell_davies1981, smith1991), unverified.";

pub const PYRO_DOWNWIND_NOTE: &str = "lane PYRO2 pyro_downwind=1: the lapilli limit of a fall is elongated downwind as lane VE's fall isopachs (x2.5 downwind, x0.6 upwind of the vent): along the direction at angle a from the fall's dispersal axis it is r_lap (1.55 + 0.95 cos a). The axis is the tephra unit's own keyed wind direction (lobe tephra mantles) or a keyed direction per package (1D sequence; separate draw at the end of the fall branch, so no other draw moves). Ignimbrite basal falls keep the circular limit.";

pub const TEPHRA_REC_NOTE: &str = "lane PYRO2 de_tephra_rec=1 (needs pyro_sub=1 and de_lobes=1): the volcano-wide tephra-fall mantles of the lobe column are the RECORDED explosive eruptions (lane VE counts explosive_veiN, N >= 3; cinder cones and maars are local and stay as VE units) instead of a count from the style share: one mantle per recorded eruption, of that eruption's VEI, with lane VE's fall rule for its size (bulk volume 10^(VEI + 4.5) m3, 20 % as fall, thickness e-folding length bb(VEI) = 1 / 2 / 4 / 8 km at VEI 3 / 4 / 5 / 6 (VE's apron rule; VEI 3 declared), vent thickness t0 = 0.2 V / (1.55 2 pi bb^2)); age and dispersal direction stay keyed (the same draws as before for the first mantles). The fall's subtype uses the recorded VEI (no draw). Without lane VE counts the kind default event list is used.";

pub const PYRO_HARD_NOTE: &str = "lane PYRO pyro_sub_hard=1 (needs pyro_sub=1; acts on the de_lobes columns): the strength multiplier of each pyroclastic part of a lobe column is set by its subtype (ln m mean per code, recalled order of strengths: welded ignimbrite > lahar / block-and-ash > non-welded ignimbrite / scoria > ash > pumice) plus a keyed per-unit spread, and the parts follow the unit's sub-resolution sequence (a welded core between non-welded zones erodes as a ledge over a slope). Per volcano and material (pdc / lahar lobes, tephra falls) the multipliers are rescaled so the expected mean 1/Sc over the record's event frequencies equals the BS / DE draws' expectation (thin ash tops and fall layers ignored; the lapilli share of falls taken as the area share of the r_lap disc within the cone radius). Facies 5 / 6 of those parts follow the new strength (ln m > 0 = indurated). The per-volcano 1D sequence below the lobe stack keeps its BS strengths.";

#[derive(Clone, Debug, Serialize)]
pub struct PyroHard {
    pub version: u32,
    /// ln m mean by subtype code (index 1-7; 0 unused), relative to the class Sc as the DE lobe parts.
    pub lnm: [f64; 8],
    /// Keyed per-unit spread (sigma of ln m).
    pub spread: f64,
    pub note: String,
}
impl Default for PyroHard {
    fn default() -> Self {
        PyroHard { version: 1, lnm: [0.0, -0.6, -0.8, -0.9, -0.7, 0.2, -0.3, -0.4], spread: 0.25, note: PYRO_HARD_NOTE.into() }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct PyroParams {
    pub version: u32,
    /// Thickness e-folding length of fall by VEI (index 0-6, m) and the lapilli limit as a multiple of it.
    pub bb_m: [f64; 7],
    pub r_lap_bb_mult: f64,
    /// PDC events of VEI >= this are ignimbrites; ignimbrites of VEI >= weld_vei_min weld where locally thicker than weld_min_m.
    pub ign_vei_min: u8,
    pub weld_vei_min: u8,
    pub weld_min_m: f64,
    pub note: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hard: Option<PyroHard>,
    /// Lane PYRO2 (opt-in `pyro_downwind=1`): the lapilli limit of a fall follows lane VE's elongated isopach shape. None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub downwind: Option<String>,
    /// Lane VH-I (opt-in `pyro_wind=1`): fall dispersal axes on the T0 wind published by T1b (lane VE ve_wind). None = off.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind: Option<PyroWind>,
}

pub const PYRO_WIND_NOTE: &str = "lane VH-I pyro_wind=1 (interim of the eruption-history design I2): a fall's dispersal axis is the volcano's T0 downwind azimuth as T1b publishes it (prep.net.volcanoes[k].ve.wind.chart_azimuth_rad, lane VE ve_wind=1: the climate= artifact's surface wind at the record's T0 cell, in the chart frame) plus a wrapped-Cauchy deviation 2 atan(((1 - rho) / (1 + rho)) tan(pi (u - 1/2))) of the SAME keyed uniform u that drew the uniform axis (no draw count changes); rho = 0.698 (von Mises kappa 2 mean resultant length) [H]. Volcanoes whose window lists no wind (ve_wind off, calm) keep the uniform axis.";

/// Lane VH-I pyro_wind parameters.
#[derive(Clone, Debug, Serialize)]
pub struct PyroWind { pub version: u32, pub rho_permille: i64, pub note: String }
impl Default for PyroWind { fn default() -> Self { PyroWind { version: 1, rho_permille: 698, note: PYRO_WIND_NOTE.into() } } }

/// Lane VH-I: the dispersal axis of a fall from the keyed uniform u: the T0 downwind azimuth plus a wrapped-Cauchy deviation
/// (pyro_wind on and the record carries a wind), else u x 2 pi (the axis before; same draw either way).
pub fn wind_axis(rec: &PyroRec, p: &PyroParams, u: f64) -> f64 {
    match (p.wind.as_ref(), rec.wind_rad) {
        (Some(w), Some(az)) => { let rho = w.rho_permille as f64 / 1000.0; az + 2.0 * (((1.0 - rho) / (1.0 + rho)) * (std::f64::consts::PI * (u - 0.5)).tan()).atan() }
        _ => u * std::f64::consts::TAU,
    }
}
impl Default for PyroParams {
    fn default() -> Self {
        PyroParams { version: 1, bb_m: [250.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0], r_lap_bb_mult: 3.0, ign_vei_min: 4, weld_vei_min: 5, weld_min_m: 8.0,
                     note: PYRO_SUB_NOTE.into(), hard: None, downwind: None, wind: None }
    }
}

/// One recorded (or default) event class: (VEI, weight, lapilli code for falls).
#[derive(Clone, Debug, Serialize)]
pub struct Ev { pub vei: u8, pub w: f64, pub lap: u8 }

/// The pyroclastic record of one volcano.
#[derive(Clone, Debug, Serialize)]
pub struct PyroRec {
    pub volcano: String,
    pub kind: String,
    pub evolved: bool,
    pub source: &'static str,
    pub fall: Vec<Ev>,
    pub pdc: Vec<Ev>,
    pub lahar: f64,
    pub vent: [f64; 2],
    pub r_cone_m: f64,
    /// pyro_sub_hard rescale factors (pdc / lahar lobes, tephra falls); 1 when off.
    pub k_pdc: f64,
    pub k_fall: f64,
    /// Lane VH-I (pyro_wind=1 only): the T0 downwind azimuth (chart frame, rad) published by T1b; None = uniform axis.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_rad: Option<f64>,
}

impl PyroRec {
    /// From the T1b volcano kind, lane VE's counts (header ve.counts) when present, the vent and the cone radius.
    pub fn new(id: &str, kind: &str, counts: Option<&serde_json::Value>, vent: [f64; 2], r_cone_m: f64, p: &PyroParams) -> PyroRec {
        let evolved = matches!(kind, "continental_arc" | "island_arc");
        let lap_of = |vei: u8| if vei >= 4 && evolved { SUB_PUMICE } else { SUB_SCORIA };
        let (mut fall, mut pdc) = (Vec::new(), Vec::new());
        let lahar;
        let add = |v: &mut Vec<Ev>, vei: u8, w: f64, lap: u8| { if w <= 0.0 { return; } match v.iter_mut().find(|e| e.vei == vei && e.lap == lap) { Some(e) => e.w += w, None => v.push(Ev { vei, w, lap }) } };
        let source;
        match counts.and_then(|c| c.as_object()) {
            Some(c) => {
                source = "ve_counts";
                let n = |k: &str| c.get(k).and_then(|v| v.as_f64()).unwrap_or(0.0);
                add(&mut fall, 2, n("cinder_cone"), SUB_SCORIA);
                add(&mut fall, 2, n("maar"), SUB_ASH);
                for vei in 3u8..=6 { let k = n(&format!("explosive_vei{vei}")); add(&mut fall, vei, k, lap_of(vei)); add(&mut pdc, vei, k, 0); }
                add(&mut pdc, 3, n("dome"), 0);
                lahar = n("lahar");
            }
            None => {
                source = "kind_default";
                match kind {
                    "continental_arc" | "island_arc" => {
                        add(&mut fall, 2, 1.0, SUB_SCORIA); add(&mut fall, 3, 4.0, SUB_SCORIA); add(&mut fall, 4, 2.0, lap_of(4)); add(&mut fall, 5, 1.0, lap_of(5));
                        add(&mut pdc, 3, 6.0, 0); add(&mut pdc, 4, 2.0, 0); add(&mut pdc, 5, 1.0, 0); lahar = 8.0;
                    }
                    "hotspot" | "rift" => { add(&mut fall, 2, 4.0, SUB_SCORIA); add(&mut fall, 3, 1.0, SUB_SCORIA); add(&mut pdc, 3, 0.5, 0); lahar = 1.0; }
                    _ => { add(&mut fall, 2, 1.0, SUB_SCORIA); add(&mut pdc, 3, 0.1, 0); lahar = 0.2; }
                }
            }
        }
        if fall.is_empty() { add(&mut fall, 2, 1.0, if evolved { SUB_ASH } else { SUB_SCORIA }); }
        let mut r = PyroRec { volcano: id.into(), kind: kind.into(), evolved, source, fall, pdc, lahar, vent, r_cone_m, k_pdc: 1.0, k_fall: 1.0, wind_rad: None };
        if let Some(h) = p.hard.as_ref() { let (a, b) = r.hard_factors(p, h); r.k_pdc = a; r.k_fall = b; }
        r
    }

    pub fn r_lap(&self, p: &PyroParams, vei: u8) -> f64 { p.r_lap_bb_mult * p.bb_m[(vei as usize).min(6)] }
    /// Lapilli code of a fall of this VEI (pumice for evolved magma at VEI >= 4, else scoria).
    pub fn lap_of(&self, vei: u8) -> u8 { if vei >= 4 && self.evolved { SUB_PUMICE } else { SUB_SCORIA } }
    /// Lane PYRO2 de_tephra_rec: the recorded explosive eruptions (VEI >= 3) as a list of VEI, largest first (counts rounded).
    pub fn explosive_events(&self) -> Vec<u8> {
        let mut v: Vec<u8> = Vec::new();
        let mut ev: Vec<&Ev> = self.fall.iter().filter(|e| e.vei >= 3).collect();
        ev.sort_by(|a, b| b.vei.cmp(&a.vei));
        for e in ev { for _ in 0..(e.w.round().max(0.0) as usize) { v.push(e.vei); } }
        v
    }

    /// pyro_sub_hard: rescale factors so E[1/m] over the record equals the BS / DE draws' expectation (PYRO_HARD_NOTE).
    fn hard_factors(&self, p: &PyroParams, h: &PyroHard) -> (f64, f64) {
        let s2 = 0.5 * h.spread * h.spread;
        let inv = |c: u8| (-h.lnm[c as usize] + s2).exp();
        // DE unit_parts: pdc weak (85 %) ln m -0.55 + 0.4 z, indurated (15 %) 0.1 + 0.2 z; tephra -0.7 + 0.3 z
        let old_pdc = 0.85 * (0.55f64 + 0.08).exp() + 0.15 * (-0.1f64 + 0.02).exp();
        let old_fall = (0.7f64 + 0.045).exp();
        let mut w = self.lahar; let mut e = self.lahar * inv(SUB_LAHAR);
        for ev in &self.pdc {
            w += ev.w;
            e += ev.w * if ev.vei >= p.ign_vei_min {
                if ev.vei >= p.weld_vei_min { 0.575 * inv(SUB_IGN_W) + 0.425 * inv(SUB_IGN_NW) } else { inv(SUB_IGN_NW) }
            } else { inv(SUB_BAF) };
        }
        let k_pdc = if w > 0.0 { (e / w) / old_pdc } else { 1.0 };
        let (mut wf, mut ef) = (0.0, 0.0);
        for ev in &self.fall {
            let a = (self.r_lap(p, ev.vei) / self.r_cone_m.max(1.0)).powi(2).min(1.0);
            wf += ev.w; ef += ev.w * (a * inv(ev.lap) + (1.0 - a) * inv(SUB_ASH));
        }
        let k_fall = if wf > 0.0 { (ef / wf) / old_fall } else { 1.0 };
        (k_pdc, k_fall)
    }
}

/// Lane PYRO2 de_tephra_rec: lane VE's fall size for an eruption of this VEI: (vent thickness t0 m, e-folding length bb m).
pub fn ve_fall(vei: u8, p: &PyroParams) -> (f64, f64) {
    let bb = p.bb_m[(vei as usize).min(6)];
    let v = 10f64.powf(vei as f64 + 4.5);
    (0.2 * v / (1.55 * std::f64::consts::TAU * bb * bb), bb)
}

/// The subtype assignment of one unit (lobe unit or 1D-sequence package).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum UnitSub {
    /// Fall of VEI `vei`: `lap` (scoria / pumice / ash) within r_lap of the vent, ash beyond.
    Fall { lap: u8, r_lap: f64, wind: f64 },
    /// Recorded cinder cone or spatter rampart.
    Scoria,
    Lahar,
    /// Block-and-ash / scoria flow with an ash-cloud ash top (m).
    Baf { ash: f64 },
    /// Ignimbrite: co-ignimbrite ash top (m), basal Plinian fall (m, its lapilli code), welded-capable, non-welded top / base
    /// fractions of the body.
    Ign { ash: f64, pf: f64, lap: u8, weld: bool, nwt: f64, nwb: f64 },
}

/// Which pool of recorded events an unrecorded unit draws from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Pool { PdcLahar, Fall, All }

struct Dr { sk: Scope, n: u64 }
impl Dr {
    fn u(&mut self) -> R<f64> { self.n += 1; Ok(self.sk.draw(PURPOSE_PYRO, self.n, 0).map_err(|e| e.0)?.uniform_f64()) }
    fn z(&mut self) -> R<f64> { self.n += 1; Ok(self.sk.draw(PURPOSE_PYRO, self.n, 0).map_err(|e| e.0)?.normal_f64()) }
    fn ln(&mut self, med: f64, sig: f64, lo: f64, hi: f64) -> R<f64> { Ok((med * (sig * self.z()?).exp()).clamp(lo, hi)) }
}

fn pick(evs: &[Ev], u: f64) -> Option<&Ev> {
    let tot: f64 = evs.iter().map(|e| e.w).sum();
    if tot <= 0.0 { return None; }
    let mut acc = 0.0;
    for e in evs { acc += e.w / tot; if u < acc { return Some(e); } }
    evs.last()
}

/// Draw the subtype of an unrecorded unit from the record (keyed by `key`). Also returns the keyed per-unit strength z.
pub fn draw_sub(rec: &PyroRec, p: &PyroParams, scope: &Scope, key: &str, pool: Pool) -> R<(UnitSub, f64)> {
    let mut dr = Dr { sk: scope.address(&format!("{key}/pyro")).map_err(|e| e.0)?, n: 0 };
    let (u0, u1, zs) = (dr.u()?, dr.u()?, dr.z()?);
    let wf: f64 = rec.fall.iter().map(|e| e.w).sum();
    let wp: f64 = rec.pdc.iter().map(|e| e.w).sum();
    // class: fall, pdc or lahar with the record's frequencies (within the pool)
    let (cf, cp, cl) = match pool { Pool::Fall => (wf.max(1e-12), 0.0, 0.0), Pool::PdcLahar => (0.0, wp, rec.lahar), Pool::All => (wf, wp, rec.lahar) };
    let tot = cf + cp + cl;
    let s = if tot <= 0.0 { UnitSub::Lahar } else {
        let x = u0 * tot;
        if x < cf {
            let e = pick(&rec.fall, u1).cloned().unwrap_or(Ev { vei: 2, w: 1.0, lap: SUB_SCORIA });
            // the dispersal direction is the last draw of this branch (pyro_downwind; no other draw moves)
            UnitSub::Fall { lap: e.lap, r_lap: rec.r_lap(p, e.vei), wind: wind_axis(rec, p, dr.u()?) }
        } else if x < cf + cp {
            let e = pick(&rec.pdc, u1).cloned().unwrap_or(Ev { vei: 3, w: 1.0, lap: 0 });
            if e.vei >= p.ign_vei_min {
                let lap = if rec.evolved { SUB_PUMICE } else { SUB_SCORIA };
                let (ash, pf, nwt, nwb) = (dr.ln(0.4, 0.5, 0.1, 1.5)?, dr.ln(0.8, 0.6, 0.1, 3.0)?, 0.15 + 0.2 * dr.u()?, 0.10 + 0.15 * dr.u()?);
                UnitSub::Ign { ash, pf, lap, weld: e.vei >= p.weld_vei_min, nwt, nwb }
            } else { UnitSub::Baf { ash: dr.ln(0.3, 0.5, 0.05, 1.0)? } }
        } else { UnitSub::Lahar }
    };
    Ok((s, zs))
}

impl UnitSub {
    /// Zones top-down as (top fraction, base fraction, code) of the local deposited thickness `tu` at distance `r` from the vent.
    pub fn zones(&self, tu: f64, r: f64, p: &PyroParams) -> Vec<(f64, f64, u8)> {
        let tu = tu.max(1e-6);
        match *self {
            UnitSub::Fall { lap, r_lap, .. } => vec![(0.0, 1.0, if r <= r_lap { lap } else { SUB_ASH })],
            UnitSub::Scoria => vec![(0.0, 1.0, SUB_SCORIA)],
            UnitSub::Lahar => vec![(0.0, 1.0, SUB_LAHAR)],
            UnitSub::Baf { ash } => { let a = ash.min(0.2 * tu) / tu; vec![(0.0, a, SUB_ASH), (a, 1.0, SUB_BAF)] }
            UnitSub::Ign { ash, pf, lap, weld, nwt, nwb } => {
                let a = ash.min(0.15 * tu) / tu;
                let f = pf.min(0.3 * tu) / tu;
                let fall = if r <= p.r_lap_bb_mult * p.bb_m[4] { lap } else { SUB_ASH };
                let body = 1.0 - a - f;
                let mut z = vec![(0.0, a, SUB_ASH)];
                if weld && tu >= p.weld_min_m {
                    let (t1, t2) = (a + nwt * body, 1.0 - f - nwb * body);
                    z.push((a, t1, SUB_IGN_NW)); z.push((t1, t2, SUB_IGN_W)); z.push((t2, 1.0 - f, SUB_IGN_NW));
                } else { z.push((a, 1.0 - f, SUB_IGN_NW)); }
                z.push((1.0 - f, 1.0, fall));
                z
            }
        }
    }
    /// Lane PYRO2: the distance a code lookup uses for a post at (dx, dy) from the vent: the true distance, or with
    /// pyro_downwind=1 for a fall the distance scaled by lane VE's isopach elongation along the fall's axis (PYRO_DOWNWIND_NOTE).
    pub fn r_eff(&self, dx: f64, dy: f64, p: &PyroParams) -> f64 {
        let r = dx.hypot(dy);
        match (*self, p.downwind.is_some()) {
            (UnitSub::Fall { wind, .. }, true) if r > 0.0 => r / (1.55 + 0.95 * (dy.atan2(dx) - wind).cos()),
            _ => r,
        }
    }
    /// The body code: the thickest zone (the unit's main material).
    pub fn body(&self, tu: f64, r: f64, p: &PyroParams) -> u8 {
        self.zones(tu, r, p).into_iter().fold((-1.0f64, SUB_UNRESOLVED), |a, q| if q.1 - q.0 > a.0 { (q.1 - q.0, q.2) } else { a }).1
    }
    /// Code at fraction `frac` (0 = the deposited top, 1 = the base) of the local deposited thickness. Above the deposited top
    /// (frac < 0: the DATUM2 column rule extends the youngest present keyed unit up to the recorded units / D) the body code:
    /// the extension is the unit's material, not its thin top layer.
    pub fn code_at(&self, frac: f64, tu: f64, r: f64, p: &PyroParams) -> u8 {
        if frac < -1e-9 { return self.body(tu, r, p); }
        let f = frac.clamp(0.0, 1.0);
        let z = self.zones(tu, r, p);
        for q in &z { if f < q.1 { return q.2; } }
        z.last().map(|q| q.2).unwrap_or(SUB_UNRESOLVED)
    }
}

/// The records of every listed volcano of a T1b window, and a cache of the 1D-sequence package subtypes.
pub struct PyroWorld {
    pub p: PyroParams,
    pub recs: Vec<PyroRec>,
    cache: std::sync::Mutex<std::collections::HashMap<(u8, u16), UnitSub>>,
}

impl PyroWorld {
    pub fn build(p: &PyroParams, w: &crate::input::T1bWindow) -> PyroWorld {
        let recs = w.volcanoes.iter().enumerate().map(|(i, v)| {
            let g = w.volcano_geo.get(i).copied().unwrap_or([f64::NAN; 4]);
            let r_cone = w.volcano_r_body.get(i).copied().flatten().unwrap_or(g[2]);
            let mut r = PyroRec::new(&v.id, &v.kind, w.volcano_ve_counts.get(i).and_then(|c| c.as_ref()), [g[0], g[1]], r_cone, p);
            if p.wind.is_some() { r.wind_rad = w.volcano_ve_wind.get(i).copied().flatten(); }
            r
        }).collect();
        PyroWorld { p: p.clone(), recs, cache: Default::default() }
    }
    /// Code of a facies 5 / 6 bed of volcano v's 1D sequence (BS beds): package `unit` (one event drawn from all classes,
    /// keyed by the volcano id and the unit), fraction `frac` of its thickness `tu`, post at (x, y).
    pub fn ed_code(&self, scope: &Scope, v: u8, unit: u16, frac: f64, tu: f64, x: f64, y: f64) -> u8 {
        let Some(rec) = self.recs.get(v as usize) else { return SUB_UNRESOLVED };
        if unit == u16::MAX { return SUB_UNRESOLVED; }
        let s = { self.cache.lock().unwrap().get(&(v, unit)).copied() };
        let s = match s {
            Some(s) => s,
            None => match draw_sub(rec, &self.p, scope, &format!("edifice/{}/pkg/{unit}", rec.volcano), Pool::All) {
                Ok((s, _)) => { self.cache.lock().unwrap().insert((v, unit), s); s }
                Err(_) => return SUB_UNRESOLVED,
            },
        };
        let r = s.r_eff(x - rec.vent[0], y - rec.vent[1], &self.p);
        s.code_at(frac, tu, if r.is_finite() { r } else { f64::INFINITY }, &self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn scope() -> Scope { Scope::root(42).layer("t2-ground-detail/test").unwrap().chart(0).unwrap() }
    fn counts() -> serde_json::Value { serde_json::json!({"cinder_cone": 3, "dome": 3, "explosive_vei3": 8, "explosive_vei5": 1, "lahar": 11, "lava_flow": 18}) }
    #[test]
    fn record_frequencies_and_distance_set_the_subtypes() {
        let p = PyroParams::default();
        let rec = PyroRec::new("v", "continental_arc", Some(&counts()), [0.0, 0.0], 10_000.0, &p);
        assert_eq!(rec.source, "ve_counts");
        assert!((rec.lahar - 11.0).abs() < 1e-12);
        let s = scope();
        // keyed pdc / lahar lobes: the lahar share follows the record (11 / 23) within sampling error; ignimbrites only from VEI 5
        let (mut nl, mut ni, mut nb, n) = (0, 0, 0, 3000);
        for i in 0..n {
            match draw_sub(&rec, &p, &s, &format!("u/{i}"), Pool::PdcLahar).unwrap().0 { UnitSub::Lahar => nl += 1, UnitSub::Ign { weld, .. } => { assert!(weld); ni += 1 } UnitSub::Baf { .. } => nb += 1, o => panic!("{o:?}") }
        }
        let fl = nl as f64 / n as f64;
        assert!((fl - 11.0 / 23.0).abs() < 0.04, "lahar share {fl}");
        assert!(ni > 50 && ni < 220 && nb > 1000, "ign {ni} baf {nb}");
        // the same key gives the same subtype (keyed), and the draw does not depend on anything else
        assert_eq!(draw_sub(&rec, &p, &s, "u/7", Pool::PdcLahar).unwrap(), draw_sub(&rec, &p, &s, "u/7", Pool::PdcLahar).unwrap());
        // falls: lapilli near the vent, ash beyond the lapilli limit; Plinian (VEI 5) lapilli of an arc are pumice
        let f5 = UnitSub::Fall { lap: SUB_PUMICE, r_lap: rec.r_lap(&p, 5), wind: 0.0 };
        assert_eq!(f5.code_at(0.5, 2.0, 3000.0, &p), SUB_PUMICE);
        assert_eq!(f5.code_at(0.5, 2.0, 30_000.0, &p), SUB_ASH);
        let f2 = UnitSub::Fall { lap: SUB_SCORIA, r_lap: rec.r_lap(&p, 2), wind: 0.0 };
        assert_eq!(f2.code_at(0.5, 2.0, 1000.0, &p), SUB_SCORIA);
        assert_eq!(f2.code_at(0.5, 2.0, 3000.0, &p), SUB_ASH);
        // a basaltic shield has no pumice
        let sh = PyroRec::new("h", "hotspot", None, [0.0, 0.0], 10_000.0, &p);
        assert!(sh.fall.iter().all(|e| e.lap == SUB_SCORIA), "{:?}", sh.fall);
    }
    #[test]
    fn downwind_limit_and_recorded_falls() {
        let mut p = PyroParams::default();
        let f = UnitSub::Fall { lap: SUB_PUMICE, r_lap: 12_000.0, wind: 0.0 };
        // off: circular
        assert_eq!(f.r_eff(-15_000.0, 0.0, &p), 15_000.0);
        p.downwind = Some(PYRO_DOWNWIND_NOTE.into());
        // downwind (x 2.5) lapilli reach 30 km, upwind (x 0.6) 7.2 km
        assert_eq!(f.code_at(0.5, 1.0, f.r_eff(25_000.0, 0.0, &p), &p), SUB_PUMICE);
        assert_eq!(f.code_at(0.5, 1.0, f.r_eff(-8_000.0, 0.0, &p), &p), SUB_ASH);
        assert_eq!(f.code_at(0.5, 1.0, f.r_eff(-7_000.0, 0.0, &p), &p), SUB_PUMICE);
        // recorded explosive eruptions, largest first; VE's fall size grows with VEI
        let rec = PyroRec::new("v", "continental_arc", Some(&counts()), [0.0, 0.0], 10_000.0, &p);
        assert_eq!(rec.explosive_events(), vec![5, 3, 3, 3, 3, 3, 3, 3, 3]);
        let (t3, b3) = ve_fall(3, &p); let (t5, b5) = ve_fall(5, &p);
        assert!(t5 > 5.0 * t3 && b5 == 4000.0 && b3 == 1000.0 && (t5 - 4.06).abs() < 0.05, "{t3} {t5}");
    }
    #[test]
    fn wind_axis_uses_the_same_draw() {
        // Lane VH-I: off (or no published wind) the axis is u x 2 pi as before; on, the T0 azimuth plus a deviation of u.
        let mut p = PyroParams::default();
        let mut rec = PyroRec::new("v", "continental_arc", Some(&counts()), [0.0, 0.0], 10_000.0, &p);
        assert!(!serde_json::to_string(&p).unwrap().contains("wind\"") && rec.wind_rad.is_none());
        assert_eq!(wind_axis(&rec, &p, 0.3), 0.3 * std::f64::consts::TAU);
        p.wind = Some(PyroWind::default());
        assert_eq!(wind_axis(&rec, &p, 0.3), 0.3 * std::f64::consts::TAU, "no published wind: uniform axis");
        rec.wind_rad = Some(1.0);
        assert!((wind_axis(&rec, &p, 0.5) - 1.0).abs() < 1e-12 && wind_axis(&rec, &p, 0.6) > 1.0 && wind_axis(&rec, &p, 0.4) < 1.0);
    }
    #[test]
    fn ignimbrite_zonation_needs_thickness() {
        let p = PyroParams::default();
        let ig = UnitSub::Ign { ash: 0.4, pf: 0.8, lap: SUB_PUMICE, weld: true, nwt: 0.25, nwb: 0.15 };
        // thick: ash top, non-welded top, welded core, non-welded base, basal pumice fall; contiguous
        let z = ig.zones(20.0, 2000.0, &p);
        let codes: Vec<u8> = z.iter().map(|q| q.2).collect();
        assert_eq!(codes, vec![SUB_ASH, SUB_IGN_NW, SUB_IGN_W, SUB_IGN_NW, SUB_PUMICE]);
        assert!(z[0].0 == 0.0 && (z.last().unwrap().1 - 1.0).abs() < 1e-12);
        for w in z.windows(2) { assert!((w[0].1 - w[1].0).abs() < 1e-12 && w[1].1 >= w[1].0); }
        assert_eq!(ig.code_at(0.5, 20.0, 2000.0, &p), SUB_IGN_W);
        // thin: no welded core
        assert!(!ig.zones(5.0, 2000.0, &p).iter().any(|q| q.2 == SUB_IGN_W));
        // the basal fall is ash beyond the VEI 4 lapilli limit
        assert_eq!(ig.code_at(0.999, 20.0, 50_000.0, &p), SUB_ASH);
        // above the deposited top (column extension): the body, not the ash top
        assert_eq!(ig.code_at(-0.5, 5.0, 2000.0, &p), SUB_IGN_NW);
        assert_eq!(UnitSub::Baf { ash: 0.3 }.code_at(-0.2, 4.0, 2000.0, &p), SUB_BAF);
        assert_eq!(UnitSub::Baf { ash: 0.3 }.code_at(0.0, 4.0, 2000.0, &p), SUB_ASH);
    }
    #[test]
    fn keys_are_opt_in() {
        let mut p = crate::DetailParamsV1::v3();
        let j0 = serde_json::to_string(&p).unwrap();
        assert!(!j0.contains("pyro"));
        assert!(p.set("pyro_sub", "1").is_err(), "needs beds=2 and bed_meta=1");
        p.set("beds", "2").unwrap();
        assert!(p.set("pyro_sub", "1").is_err(), "needs bed_meta=1");
        p.set("bed_meta", "1").unwrap();
        assert!(p.set("pyro_sub_hard", "1").is_err(), "needs pyro_sub=1");
        let j1 = serde_json::to_string(&p).unwrap();
        assert!(!j1.contains("pyro"));
        p.set("pyro_sub", "1").unwrap();
        assert!(p.pyro.as_ref().unwrap().hard.is_none());
        assert!(!serde_json::to_string(&p).unwrap().contains("\"hard\""));
        p.set("pyro_sub_hard", "1").unwrap();
        assert!(p.pyro.as_ref().unwrap().hard.is_some());
        p.set("pyro_sub", "0").unwrap();
        assert!(p.pyro.is_none() && serde_json::to_string(&p).unwrap() == j1);
    }
    #[test]
    fn hard_factors_keep_the_expected_strength() {
        let mut p = PyroParams::default();
        p.hard = Some(PyroHard::default());
        let rec = PyroRec::new("v", "continental_arc", Some(&counts()), [0.0, 0.0], 10_000.0, &p);
        assert!(rec.k_pdc > 0.2 && rec.k_pdc < 5.0 && rec.k_fall > 0.2 && rec.k_fall < 5.0, "{} {}", rec.k_pdc, rec.k_fall);
        let off = PyroRec::new("v", "continental_arc", Some(&counts()), [0.0, 0.0], 10_000.0, &PyroParams::default());
        assert_eq!((off.k_pdc, off.k_fall), (1.0, 1.0));
    }
}
