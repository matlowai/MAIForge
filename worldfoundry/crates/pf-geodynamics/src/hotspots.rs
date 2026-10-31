//! Fixed mantle hotspots and plate-crossing detection (brief K6;
//! acceptance: seeded fixed points, crossing events with age-progressive
//! chains; NO parcel transport).
//!
//! Hotspots are fixed in the mantle frame (the datum's inertial frame):
//! seeded once from the SeedPath's independent "hotspots" stream, they
//! never move. A crossing happens when the PLATE under a hotspot changes
//! between two consecutive plate assignments (the plate's motion carried
//! new crust over the plume). Each crossing mints a HotspotCrossing event
//! into the caller's event ledger — this module detects and proposes; the
//! ledger owns the records.

use crate::events::EventLedger;
use crate::plates::{seed_from_path_tagged, SplitMix64};
use pf_contracts::domain::{EventKindV0, PlateId};
use pf_contracts::params::SeedPath;
use pf_contracts::spatial::ElementId;

#[derive(Debug, Clone, PartialEq)]
pub enum HotspotsError {
    NonFinite { what: &'static str, value: f64 },
    CountMismatch { expected: usize, assignment: usize },
    EmptyAssignment,
    NoElementForHotspot { hotspot: usize },
}

impl std::fmt::Display for HotspotsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HotspotsError::NonFinite { what, value } => write!(f, "nonfinite {what}: {value:?}"),
            HotspotsError::CountMismatch {
                expected,
                assignment,
            } => {
                write!(
                    f,
                    "plate assignment covers {assignment} of {expected} elements"
                )
            }
            HotspotsError::EmptyAssignment => write!(f, "empty plate assignment"),
            HotspotsError::NoElementForHotspot { hotspot } => {
                write!(
                    f,
                    "no containing element found for hotspot {hotspot} (grid topology broken)"
                )
            }
        }
    }
}

impl std::error::Error for HotspotsError {}

/// One fixed mantle hotspot (unit vector in the datum frame + its id).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hotspot {
    pub id: usize,
    pub position: [f64; 3],
}

/// Seed `count` hotspots uniformly on the sphere from the SeedPath's
/// independent hotspot stream (draw order fixed: position per hotspot in
/// id order).
pub fn seed_hotspots(count: u32, seed_path: &SeedPath) -> Vec<Hotspot> {
    let mut rng = SplitMix64::new(seed_from_path_tagged(seed_path, "hotspots"));
    (0..count as usize)
        .map(|id| {
            let z = 2.0 * rng.next_f64() - 1.0;
            let phi = 2.0 * std::f64::consts::PI * rng.next_f64();
            let r = (1.0 - z * z).max(0.0).sqrt();
            Hotspot {
                id,
                position: [r * phi.cos(), r * phi.sin(), z],
            }
        })
        .collect()
}

/// The element whose centroid contains a point: nearest centroid on the
/// grid (declared rule; ties break to the lowest element id — same
/// convention as plate seeding).
pub fn containing_element(
    point: [f64; 3],
    centroids: &[[f64; 3]],
) -> Result<ElementId, HotspotsError> {
    if centroids.is_empty() {
        return Err(HotspotsError::EmptyAssignment);
    }
    let mut best_d = f64::INFINITY;
    let mut best_k = 0usize;
    for (k, c) in centroids.iter().enumerate() {
        let d =
            crate::rotation::great_circle_rad(point, *c).map_err(|_| HotspotsError::NonFinite {
                what: "centroid distance",
                value: 0.0,
            })?;
        if d < best_d {
            best_d = d;
            best_k = k;
        }
    }
    Ok(ElementId(best_k as u64))
}

/// State of the plates moving under the hotspots between assignments.
#[derive(Debug, Clone, PartialEq)]
pub struct HotspotTracker {
    hotspots: Vec<Hotspot>,
    /// Current plate under each hotspot (index-aligned with hotspots).
    current_plate: Vec<Option<PlateId>>,
}

impl HotspotTracker {
    pub fn new(hotspots: Vec<Hotspot>) -> Self {
        let n = hotspots.len();
        Self {
            hotspots,
            current_plate: vec![None; n],
        }
    }

    pub fn hotspots(&self) -> &[Hotspot] {
        &self.hotspots
    }

    /// Feed one plate assignment (per element, same indexing as the grid
    /// flat order) plus the element centroids. On the FIRST feed the state
    /// initializes; afterwards, every plate change under a hotspot mints a
    /// HotspotCrossing event (K6) into `ledger` and returns the minted
    /// event ids.
    pub fn observe(
        &mut self,
        t_my: f64,
        assignment: &[PlateId],
        centroids: &[[f64; 3]],
        ledger: &mut EventLedger,
    ) -> Result<Vec<pf_contracts::domain::EventId>, HotspotsError> {
        if centroids.len() != assignment.len() {
            return Err(HotspotsError::CountMismatch {
                expected: assignment.len(),
                assignment: centroids.len(),
            });
        }
        let mut minted = Vec::new();
        for (k, hs) in self.hotspots.iter().enumerate() {
            let element = containing_element(hs.position, centroids)?;
            let plate = assignment[element.0 as usize];
            match self.current_plate[k] {
                None => self.current_plate[k] = Some(plate),
                Some(prev) if prev != plate => {
                    let id = ledger
                        .push(
                            t_my,
                            EventKindV0::HotspotCrossing,
                            vec![plate, prev],
                            vec![element],
                            vec![],
                            format!(
                                "hotspot {} crossed by plate {:?} -> {:?}",
                                hs.id, prev, plate
                            ),
                        )
                        .map_err(|_| HotspotsError::NonFinite {
                            what: "event push",
                            value: t_my,
                        })?;
                    minted.push(id);
                    self.current_plate[k] = Some(plate);
                }
                Some(_) => {}
            }
        }
        Ok(minted)
    }
}

/// Age-progressive chain helper (test/inspection): crossings on ONE plate,
/// ordered by event time. A plate sweeping over a fixed hotspot yields
/// strictly increasing times (K6's island-chain signature).
pub fn chain_on_plate(
    events: &[pf_contracts::domain::TectonicEventV0],
    plate: PlateId,
) -> Vec<f64> {
    let mut times: Vec<f64> = events
        .iter()
        .filter(|e| e.kind == EventKindV0::HotspotCrossing && e.plates.contains(&plate))
        .map(|e| e.t_my)
        .collect();
    times.sort_by(|a, b| {
        a.partial_cmp(b)
            .expect("event times are finite (ledger invariant)")
    });
    times
}

// ---------------------------------------------------------------------------
// SEED-EXTRAS (lane EX, 2026-10-12): `fixed-plume-chains-v1`.
//
// The K6 machinery above detects crossings and adds no crust. The plume law
// below adds magmatic crust: each plume is fixed in the mantle (datum) frame
// (positions from the same "hotspots" stream as `seed_hotspots`), switches on
// at a seeded onset time with a HEAD (a large igneous province spread over a
// wide disc), and from then on a TAIL adds a steady flux along the
// age-progressive track: the crust that sat over the plume at time t inside
// the step has been carried by its plate to rotate(p, pole, t_end - t) at the
// step end, so successive steps leave a chain whose age increases away from
// the plume (Hawaii-Emperor style) on oceanic plates and a volcanic track or
// plateau on continental ones.
// ---------------------------------------------------------------------------

/// Parameters of `fixed-plume-chains-v1` (pf-contracts laws.rs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlumeParams {
    pub head_volume_m3: f64,
    pub head_radius_m: f64,
    pub tail_flux_m3_per_my: f64,
    pub edifice_radius_m: f64,
}

/// One fixed mantle plume (datum frame) with its seeded onset time and the
/// azimuth of the cut used if its head splits a plate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plume {
    pub id: usize,
    pub position: [f64; 3],
    pub onset_my: f64,
    pub split_azimuth_rad: f64,
}

/// Seed `count` plumes: positions are `seed_hotspots(count, ..)` (the K6
/// stream); onsets uniform in [-lead_my, duration_my) and cut azimuths
/// uniform in [0, 2 pi) come from an independent "plume-onset" stream (draw
/// order: onset then azimuth per plume, in id order).
pub fn seed_plumes(count: u32, lead_my: f64, duration_my: f64, seed_path: &SeedPath) -> Vec<Plume> {
    let mut rng = SplitMix64::new(seed_from_path_tagged(seed_path, "plume-onset"));
    seed_hotspots(count, seed_path)
        .into_iter()
        .map(|h| {
            let onset_my = -lead_my + (duration_my + lead_my) * rng.next_f64();
            let split_azimuth_rad = 2.0 * std::f64::consts::PI * rng.next_f64();
            Plume {
                id: h.id,
                position: h.position,
                onset_my,
                split_azimuth_rad,
            }
        })
        .collect()
}

/// What one plume pass added (m^3 and per-cell increments). `cells` is in
/// ascending cell order: (cell, added thickness m, latest eruption time My).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PlumeReport {
    pub head_m3: f64,
    pub tail_m3: f64,
    pub heads: u64,
    pub tail_samples: u64,
    pub cells: Vec<(u32, f64, f64)>,
}

/// Add the plume heads whose onset falls in (t_old, t_new] and every active
/// plume's tail over the step to `thickness_m` (the post-carry crust of the
/// step end; `plate` is its plate field). Deterministic: plumes in id order,
/// samples in time order, disc cells ascending.
#[allow(clippy::too_many_arguments)]
pub fn plume_pass(
    grid: &pf_grid::GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    plumes: &[Plume],
    params: &PlumeParams,
    t_old: f64,
    t_new: f64,
    plate: &[PlateId],
    poles: &std::collections::BTreeMap<PlateId, pf_contracts::domain::EulerPoleV0>,
    thickness_m: &mut [f64],
) -> Result<PlumeReport, String> {
    use crate::kinematic::nearest_cell;
    use crate::orogen::disc_weights;
    let n = grid.cells();
    if cents.len() != n || areas.len() != n || plate.len() != n || thickness_m.len() != n {
        return Err("plume_pass: field sizes differ from the grid".into());
    }
    let radius = grid.radius_m();
    let dx = radius * std::f64::consts::PI / grid.nlat() as f64;
    let mut rep = PlumeReport::default();
    let mut added: std::collections::BTreeMap<usize, (f64, f64)> = std::collections::BTreeMap::new();
    let mut put = |weights: &[(usize, f64)], v: f64, t: f64, th: &mut [f64]| -> f64 {
        let norm: f64 = weights.iter().map(|(c, w)| w * areas[*c]).sum();
        if !(norm > 0.0) || v == 0.0 {
            return 0.0;
        }
        let mut placed = 0.0;
        for &(c, w) in weights {
            let dh = v * w / norm;
            th[c] += dh;
            placed += dh * areas[c];
            let e = added.entry(c).or_insert((0.0, t));
            e.0 += dh;
            if t > e.1 {
                e.1 = t;
            }
        }
        placed
    };
    for pl in plumes {
        if !(pl.onset_my < t_new) {
            continue;
        }
        let t_start = pl.onset_my.max(t_old);
        let active = t_new - t_start;
        if !(active > 0.0) {
            continue;
        }
        let c0 = nearest_cell(grid, pl.position).map_err(|e| e.to_string())?;
        if pl.onset_my > t_old {
            // Head: once, over a wide disc, any plate.
            let w = disc_weights(grid, cents, c0, params.head_radius_m, &|_| true)?;
            rep.head_m3 += put(&w, params.head_volume_m3, pl.onset_my, thickness_m);
            rep.heads += 1;
            if crate::volcanoes::supply_capture::enabled() && params.head_volume_m3 > 0.0 {
                crate::volcanoes::supply_capture::record(vec![crate::volcanoes::supply_capture::Sample::placed(
                    "plume_head", format!("plume/{}/head", pl.id), c0 as u32,
                    t_old, t_new, params.head_volume_m3, &w, areas,
                )?])?;
            }
        }
        // Tail along the step's age-progressive track on the overlying plate.
        let pid = plate[c0];
        let pole = poles
            .get(&pid)
            .ok_or_else(|| format!("plume {}: no pole for plate {}", pl.id, pid.0))?;
        let w_vec = crate::boundaries::angular_velocity(pole).map_err(|e| e.to_string())?;
        let p = pl.position;
        let cr = [
            w_vec[1] * p[2] - w_vec[2] * p[1],
            w_vec[2] * p[0] - w_vec[0] * p[2],
            w_vec[0] * p[1] - w_vec[1] * p[0],
        ];
        let speed = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt() * radius;
        let track_m = speed * active;
        let m = ((track_m / (0.5 * dx)).ceil() as u64).clamp(1, 100_000);
        let v = params.tail_flux_m3_per_my * active / m as f64;
        for j in 0..m {
            let t_j = t_start + (j as f64 + 0.5) / m as f64 * active;
            let x = crate::rotation::rotate_by_euler(p, pole, t_new - t_j).map_err(|e| e.to_string())?;
            let c = nearest_cell(grid, x).map_err(|e| e.to_string())?;
            let w = disc_weights(grid, cents, c, params.edifice_radius_m, &|k| plate[k] == pid)?;
            rep.tail_m3 += put(&w, v, t_j, thickness_m);
            rep.tail_samples += 1;
            if crate::volcanoes::supply_capture::enabled() && v > 0.0 {
                crate::volcanoes::supply_capture::record(vec![crate::volcanoes::supply_capture::Sample::placed(
                    "plume_tail", format!("plume/{}/tail/{j}", pl.id), c as u32,
                    t_old, t_new, v, &w, areas,
                )?])?;
            }
        }
    }
    rep.cells = added
        .into_iter()
        .map(|(c, (dh, t))| (c as u32, dh, t))
        .collect();
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::events::DEFAULT_MAX_EVENTS;

    fn centroids_ring(n: usize) -> Vec<[f64; 3]> {
        // n points around the equator, plus poles for realism.
        let mut c = Vec::new();
        for k in 0..n / 2 {
            let lon = (k as f64) * 360.0 / (n as f64 / 2.0) - 180.0;
            c.push(crate::rotation::lat_lon_to_unit(0.0, lon).unwrap());
        }
        c.push([0.0, 0.0, 1.0]);
        c.push([0.0, 0.0, -1.0]);
        c
    }

    #[test]
    fn hotspots_seed_deterministically_and_independently_of_plates() {
        let sp = SeedPath {
            seed: 5,
            derivation: "d".into(),
            salt: "s".into(),
        };
        let a = seed_hotspots(4, &sp);
        let b = seed_hotspots(4, &sp);
        assert_eq!(a, b);
        assert_eq!(a.len(), 4);
        for h in &a {
            let n = (h.position[0] * h.position[0]
                + h.position[1] * h.position[1]
                + h.position[2] * h.position[2])
                .sqrt();
            assert!((n - 1.0).abs() < 1.0e-12);
        }
        // Different tag -> independent stream (with overwhelming
        // probability different positions for count > 1).
        let plates_stream = crate::plates::seed_from_path(&sp);
        let hotspots_stream = seed_from_path_tagged(&sp, "hotspots");
        assert_ne!(plates_stream, hotspots_stream);
    }

    #[test]
    fn plate_change_under_hotspot_mints_crossing_events() {
        let centroids = centroids_ring(12);
        let n = centroids.len();
        let hotspots = vec![Hotspot {
            id: 0,
            position: centroids[0],
        }];
        let mut tracker = HotspotTracker::new(hotspots);
        let mut ledger = EventLedger::new(DEFAULT_MAX_EVENTS);

        // First observe initializes without events.
        let assignment = vec![PlateId(1); n];
        assert!(tracker
            .observe(1.0, &assignment, &centroids, &mut ledger)
            .unwrap()
            .is_empty());
        assert_eq!(ledger.len(), 0);

        // Same plate again: still silent.
        assert!(tracker
            .observe(2.0, &assignment, &centroids, &mut ledger)
            .unwrap()
            .is_empty());

        // Plate under the hotspot changes: one event, correct plates.
        let mut moved = vec![PlateId(1); n];
        moved[0] = PlateId(2);
        let minted = tracker
            .observe(3.0, &moved, &centroids, &mut ledger)
            .unwrap();
        assert_eq!(minted.len(), 1);
        let ev = &ledger.events()[0];
        assert_eq!(ev.kind, EventKindV0::HotspotCrossing);
        assert!(ev.plates.contains(&PlateId(1)) && ev.plates.contains(&PlateId(2)));
        assert_eq!(ev.footprint_elements, vec![ElementId(0)]);

        // Plate 2's involvement: entered AND left -> two crossings, both
        // Plate 2's involvement in THIS test: exactly the one 1->2
        // crossing (it never leaves here).
        let chain = chain_on_plate(ledger.events(), PlateId(2));
        assert_eq!(chain, vec![3.0]);
    }

    #[test]
    fn sweep_produces_age_ordered_chain() {
        // A hotspot sits on element 0; three successive plates pass over it.
        let centroids = centroids_ring(8);
        let n = centroids.len();
        let mut tracker = HotspotTracker::new(vec![Hotspot {
            id: 0,
            position: centroids[0],
        }]);
        let mut ledger = EventLedger::new(DEFAULT_MAX_EVENTS);
        let mut t = 0.0;
        for plate in [PlateId(1), PlateId(2), PlateId(3), PlateId(4)] {
            let mut assignment = vec![PlateId(9); n];
            assignment[0] = plate;
            tracker
                .observe(t, &assignment, &centroids, &mut ledger)
                .unwrap();
            t += 5.0;
        }
        // Three crossings total (1->2, 2->3, 3->4), strictly increasing
        // times: the age-progressive signature. Plate 2 is involved in the
        // first two (entered, then left).
        let all: Vec<f64> = ledger.events().iter().map(|e| e.t_my).collect();
        assert_eq!(all, vec![5.0, 10.0, 15.0]);
        assert!(all.windows(2).all(|w| w[0] < w[1]), "age-progressive");
        assert_eq!(chain_on_plate(ledger.events(), PlateId(2)), vec![5.0, 10.0]);
    }

    #[test]
    fn mismatched_assignment_is_loud() {
        let centroids = centroids_ring(4);
        let mut tracker = HotspotTracker::new(vec![Hotspot {
            id: 0,
            position: centroids[0],
        }]);
        let mut ledger = EventLedger::new(DEFAULT_MAX_EVENTS);
        let short = vec![PlateId(1); centroids.len() - 1];
        let err = tracker
            .observe(1.0, &short, &centroids, &mut ledger)
            .unwrap_err();
        assert!(matches!(err, HotspotsError::CountMismatch { .. }), "{err}");
    }

    #[test]
    fn containing_element_ties_break_to_lowest_id() {
        // Two identical centroids: element 0 wins by declaration.
        let c = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert_eq!(
            containing_element([1.0, 0.0, 0.0], &c).unwrap(),
            ElementId(0)
        );
        // Nearest wins in general.
        let c = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert_eq!(
            containing_element([0.0, 0.9, 0.1], &c).unwrap(),
            ElementId(1)
        );
    }
}

// ---------------------------------------------------------------------------
// SEED-MAGMATISM (lane MG, 2026-10-13): plume sources of
// `slab-arc-backarc-llsvp-v1` (`unpublished-design-record`).
//
// Plumes rise preferentially from the margins of the two antipodal deep-mantle
// provinces (LLSVPs, the "plume generation zones"), and drift slowly in the
// mantle frame so chains can bend. The plume list itself (ids, onsets, split
// azimuths) is unchanged; only positions are re-drawn and then drifted.
// ---------------------------------------------------------------------------

/// Plume-source constants of `slab-arc-backarc-llsvp-v1` (laws.rs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlumeSourceParams {
    pub edge_deg: f64,
    pub edge_sigma_deg: f64,
    pub antipode_jitter_deg: f64,
    pub background_fraction: f64,
    pub drift_rate_m_per_my: f64,
}

/// One plume's seeded drift: from `origin` at its onset along the great
/// circle in the unit tangent `dir` at `rate_rad_per_my`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlumeDrift {
    pub id: usize,
    pub origin: [f64; 3],
    pub dir: [f64; 3],
    pub rate_rad_per_my: f64,
    pub onset_my: f64,
    /// Which LLSVP edge the plume was drawn at (None = background).
    pub llsvp: Option<usize>,
}

fn unit_cross(a: [f64; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let n = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
    if n > 1.0e-9 {
        Some([c[0] / n, c[1] / n, c[2] / n])
    } else {
        None
    }
}

/// A tangent basis (e1, e2) at unit vector `c` (e1 east-like, e2 = c x e1).
fn tangent_basis(c: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let e1 = unit_cross([0.0, 0.0, 1.0], c).unwrap_or([0.0, 1.0, 0.0]);
    let e2 = unit_cross(c, e1).unwrap_or([1.0, 0.0, 0.0]);
    (e1, e2)
}

/// The point at angular distance `theta` (rad) and azimuth `az` (rad) from
/// unit vector `c`.
pub fn offset_point(c: [f64; 3], theta: f64, az: f64) -> [f64; 3] {
    let (e1, e2) = tangent_basis(c);
    let (s, co) = theta.sin_cos();
    let (sa, ca) = az.sin_cos();
    let t = [ca * e1[0] + sa * e2[0], ca * e1[1] + sa * e2[1], ca * e1[2] + sa * e2[2]];
    let p = [co * c[0] + s * t[0], co * c[1] + s * t[1], co * c[2] + s * t[2]];
    let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    [p[0] / n, p[1] / n, p[2] / n]
}

fn uniform_unit(rng: &mut SplitMix64) -> [f64; 3] {
    let z = 2.0 * rng.next_f64() - 1.0;
    let phi = 2.0 * std::f64::consts::PI * rng.next_f64();
    let r = (1.0 - z * z).max(0.0).sqrt();
    [r * phi.cos(), r * phi.sin(), z]
}

/// The two LLSVP centres: one uniform, the other its antipode moved by
/// U(0, jitter) in a uniform direction (stream "magmatism/llsvp", four
/// draws after the first centre's two).
pub fn llsvp_centres(p: &PlumeSourceParams, seed_path: &SeedPath) -> [[f64; 3]; 2] {
    let mut rng = SplitMix64::new(seed_from_path_tagged(seed_path, "magmatism/llsvp"));
    let a = uniform_unit(&mut rng);
    let jitter = p.antipode_jitter_deg.to_radians() * rng.next_f64();
    let az = 2.0 * std::f64::consts::PI * rng.next_f64();
    let b = offset_point([-a[0], -a[1], -a[2]], jitter, az);
    [a, b]
}

/// Re-draw every plume's position at an LLSVP edge (or, with probability
/// `background_fraction`, uniformly) and give it a drift. Plumes are taken
/// in slice order (ids ascending in every caller); each consumes exactly
/// nine draws of the stream "magmatism/plumes" in id order, so a plume's
/// position does not depend on how many plumes the run has.
pub fn llsvp_edge_plumes(
    plumes: &mut [Plume],
    p: &PlumeSourceParams,
    seed_path: &SeedPath,
    radius_m: f64,
) -> Result<Vec<PlumeDrift>, String> {
    if !(p.edge_deg > 0.0
        && p.edge_deg < 180.0
        && p.edge_sigma_deg >= 0.0
        && p.antipode_jitter_deg >= 0.0
        && (0.0..=1.0).contains(&p.background_fraction)
        && p.drift_rate_m_per_my >= 0.0
        && radius_m > 0.0)
    {
        return Err(format!("plume source parameters out of range: {p:?}"));
    }
    for w in plumes.windows(2) {
        if w[1].id <= w[0].id {
            return Err("llsvp_edge_plumes: plumes must be in ascending id order".into());
        }
    }
    let centres = llsvp_centres(p, seed_path);
    let mut rng = SplitMix64::new(seed_from_path_tagged(seed_path, "magmatism/plumes"));
    let mut next_id = 0usize;
    let mut out = Vec::with_capacity(plumes.len());
    for pl in plumes.iter_mut() {
        // Skip the draws of ids absent from the slice (keeps prefix stability).
        let mut d = [0.0f64; 9];
        while next_id <= pl.id {
            for x in d.iter_mut() {
                *x = rng.next_f64();
            }
            next_id += 1;
        }
        let tau = 2.0 * std::f64::consts::PI;
        let (pos, which) = if d[0] < p.background_fraction {
            let z = 2.0 * d[1] - 1.0;
            let phi = tau * d[2];
            let r = (1.0 - z * z).max(0.0).sqrt();
            ([r * phi.cos(), r * phi.sin(), z], None)
        } else {
            let k = if d[3] < 0.5 { 0 } else { 1 };
            // Box-Muller from (d[4], d[5]).
            let g = (-2.0 * (1.0 - d[4]).ln()).sqrt() * (tau * d[5]).cos();
            let theta = (p.edge_deg + p.edge_sigma_deg * g).clamp(0.0, 180.0).to_radians();
            (offset_point(centres[k], theta, tau * d[6]), Some(k))
        };
        pl.position = pos;
        let (e1, e2) = tangent_basis(pos);
        let (sa, ca) = (tau * d[7]).sin_cos();
        let dir = [ca * e1[0] + sa * e2[0], ca * e1[1] + sa * e2[1], ca * e1[2] + sa * e2[2]];
        out.push(PlumeDrift {
            id: pl.id,
            origin: pos,
            dir,
            rate_rad_per_my: p.drift_rate_m_per_my * (0.5 + d[8]) / radius_m,
            onset_my: pl.onset_my,
            llsvp: which,
        });
    }
    Ok(out)
}

/// A drifting plume's position at `t_my` (its origin before the onset).
pub fn drifted_position(d: &PlumeDrift, t_my: f64) -> [f64; 3] {
    let theta = d.rate_rad_per_my * (t_my - d.onset_my).max(0.0);
    let (s, c) = theta.sin_cos();
    let p = [
        c * d.origin[0] + s * d.dir[0],
        c * d.origin[1] + s * d.dir[1],
        c * d.origin[2] + s * d.dir[2],
    ];
    let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    [p[0] / n, p[1] / n, p[2] / n]
}

/// Move every plume with a drift record (matched by id) to its position at
/// `t_my`.
pub fn drift_plumes(plumes: &mut [Plume], drifts: &[PlumeDrift], t_my: f64) {
    for pl in plumes.iter_mut() {
        if let Ok(k) = drifts.binary_search_by_key(&pl.id, |d| d.id) {
            pl.position = drifted_position(&drifts[k], t_my);
        }
    }
}

#[cfg(test)]
mod llsvp_tests {
    use super::*;

    fn params() -> PlumeSourceParams {
        PlumeSourceParams {
            edge_deg: 35.0,
            edge_sigma_deg: 5.0,
            antipode_jitter_deg: 15.0,
            background_fraction: 0.15,
            drift_rate_m_per_my: 10_000.0,
        }
    }

    fn plumes(n: usize) -> Vec<Plume> {
        (0..n)
            .map(|id| Plume {
                id,
                position: [1.0, 0.0, 0.0],
                onset_my: id as f64 * 10.0,
                split_azimuth_rad: 0.0,
            })
            .collect()
    }

    #[test]
    fn llsvp_plumes_cluster_at_edges_and_are_prefix_stable() {
        let sp = SeedPath { seed: 42, derivation: "d".into(), salt: "s".into() };
        let p = params();
        let mut a = plumes(200);
        let da = llsvp_edge_plumes(&mut a, &p, &sp, 6.371e6).unwrap();
        let mut b = plumes(50);
        let db = llsvp_edge_plumes(&mut b, &p, &sp, 6.371e6).unwrap();
        assert_eq!(&a[..50], &b[..]);
        assert_eq!(&da[..50], &db[..]);
        let c = llsvp_centres(&p, &sp);
        let mut near = 0;
        for pl in &a {
            let d = c
                .iter()
                .map(|x| crate::rotation::great_circle_rad(*x, pl.position).unwrap().to_degrees())
                .map(|g| (g - 35.0).abs())
                .fold(f64::INFINITY, f64::min);
            if d <= 15.0 {
                near += 1;
            }
        }
        assert!(near as f64 >= 0.8 * a.len() as f64, "{near} of {}", a.len());
        // Antipodes within the jitter.
        let anti = crate::rotation::great_circle_rad(c[0], [-c[1][0], -c[1][1], -c[1][2]]).unwrap();
        assert!(anti.to_degrees() <= 15.0 + 1e-9);
    }

    #[test]
    fn drift_moves_at_the_declared_rate_after_onset() {
        let d = PlumeDrift {
            id: 0,
            origin: [1.0, 0.0, 0.0],
            dir: [0.0, 1.0, 0.0],
            rate_rad_per_my: 1.0e-3,
            onset_my: 10.0,
            llsvp: Some(0),
        };
        assert_eq!(drifted_position(&d, 5.0), [1.0, 0.0, 0.0]);
        let q = drifted_position(&d, 20.0);
        let g = crate::rotation::great_circle_rad([1.0, 0.0, 0.0], q).unwrap();
        assert!((g - 0.01).abs() < 1e-12);
    }
}
