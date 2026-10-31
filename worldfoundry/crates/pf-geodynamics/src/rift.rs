//! SEED-EXTRAS (lane EX, 2026-10-12): continental rifting and plume-head
//! plate splits.
//!
//! `stretched-rift-v1` (`rift_pass`): in the Eulerian carry a gap between
//! two diverging plates is a RIDGE BIRTH, legacy-filled with new oceanic
//! crust at once, so two continents drifting apart break up instantly with
//! a step-sharp margin. Here a birth whose nearest older (non-birth) cells
//! of two diverging plates are both continental is a RIFT birth instead:
//! the continental crust of those two plates within the rift radius
//! stretches into the gap. Exactly the volume taken from the flanks fills
//! the birth (no creation), the flanks thin most near the axis (triangular
//! weights), and the birth is as thick as its thinned flanks (closed-form
//! uniform-stretch share; no flank is stretched below the breakup
//! thickness). Thinner crust floats lower, so the rift is a
//! valley that deepens into a narrow sea as it thins. A rift birth that
//! would be thinner than the breakup thickness is an ordinary oceanic birth:
//! seafloor spreading starts between thinned conjugate margins.
//!
//! `plume-head-continental-split-v1` (`split_plate`): with fixed Euler poles
//! no new divergent boundary ever crosses a continent, so continents could
//! never split. A plume head arriving under continental crust splits the
//! plate under it along a great circle through the plume; the new plate's
//! angular velocity is the parent's plus a relative rotation about the axis
//! on the cut 90 degrees from the plume, so the halves open at the plume
//! (scissor-like, faster near the plume, closing beyond 90 degrees).
//!
//! Determinism: ascending cell order, ties to the lowest plate / cell id.

use crate::boundaries::{angular_velocity, relative_velocity_at};
use crate::rotation::{great_circle_rad, unit_to_lat_lon};
use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_grid::GridSpec;
use std::collections::BTreeMap;

/// Parameters of `stretched-rift-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RiftParams {
    pub half_width_m: f64,
    pub breakup_thickness_m: f64,
}

/// What one rift pass did. `births` are cells turned into rift crust (in
/// ascending order); `thinning` the per-cell thickness taken from flanks.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RiftReport {
    pub births: Vec<u32>,
    /// Rift candidates (continental on two diverging plates) left oceanic
    /// because they would be thinner than the breakup thickness.
    pub breakup_births: u64,
    /// Candidates with continental crust on two plates that converge (not a
    /// rift; carry churn): left oceanic as before.
    pub convergent_skipped: u64,
    pub moved_m3: f64,
    pub thinning: Vec<(u32, f64)>,
    /// SEED-RIFT-SHOULDERS (lane RS): per birth (parallel to `births`), the
    /// weighted mean pre-stretch thickness of its flanks (m), so the birth's
    /// stretching factor is this over its thickness. Diagnostic for every
    /// caller; read only by `shoulder-uplift-thermal-v1`.
    pub birth_pre_thickness_m: Vec<f64>,
    /// SEED-RIFT-SHOULDERS (lane RS): per birth (parallel to `births`), its
    /// flank cells and their weights (triangular weight x area), so a
    /// caller holding per-cell stretching factors can recover the birth's
    /// unstretched thickness.
    pub birth_flanks: Vec<Vec<(u32, f64)>>,
}

/// Cells within great-circle distance `radius_m` of cell `centre`
/// (breadth-first over 8-neighbours) with their distances (m), the centre
/// first at 0, then ascending cell order.
pub(crate) fn disc_with_distance(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    centre: usize,
    radius_m: f64,
) -> Result<Vec<(usize, f64)>, String> {
    let r = grid.radius_m();
    let mut out = vec![(centre, 0.0)];
    let mut seen = std::collections::BTreeSet::new();
    seen.insert(centre);
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(centre);
    while let Some(c) = queue.pop_front() {
        // SEED-PERF (lane SP): neighbour ids only (the distances that
        // `GridSpec::neighbors` computed were unused), same set and order.
        grid.ij(c as u64).map_err(|e| e.to_string())?;
        let mut next = [0usize; 8];
        let nn = crate::kinematic::neighbor_ids(grid.nlat(), grid.nlon(), c, &mut next);
        let next = &mut next[..nn];
        next.sort_unstable();
        for &nb in next.iter() {
            if !seen.insert(nb) {
                continue;
            }
            let d = great_circle_rad(cents[centre], cents[nb]).map_err(|e| e.to_string())? * r;
            if d <= radius_m {
                queue.push_back(nb);
                out.push((nb, d));
            }
        }
    }
    out[1..].sort_by_key(|x| x.0);
    Ok(out)
}

/// The rift search / stretching radius on this grid: the declared half
/// width, but never less than 1.5 latitude spacings (so the flanking ring
/// is always reached).
pub fn rift_radius_m(grid: &GridSpec, params: &RiftParams) -> f64 {
    let dx = grid.radius_m() * std::f64::consts::PI / grid.nlat() as f64;
    params.half_width_m.max(1.5 * dx)
}

/// Apply `stretched-rift-v1` to the births of one step. `birth[c]` marks the
/// step's ridge births (legacy-filled oceanic cells); `kind`, `thickness_m`
/// and `plate` are the post-carry fields, modified in place. Returns the
/// report; the caller re-books the converted births out of the ridge term.
#[allow(clippy::too_many_arguments)]
pub fn rift_pass(
    grid: &GridSpec,
    cents: &[[f64; 3]],
    areas: &[f64],
    params: &RiftParams,
    birth: &[bool],
    plate: &[PlateId],
    poles: &BTreeMap<PlateId, EulerPoleV0>,
    kind: &mut [CrustKindCodeV0],
    thickness_m: &mut [f64],
) -> Result<RiftReport, String> {
    let n = grid.cells();
    if [cents.len(), areas.len(), birth.len(), plate.len(), kind.len(), thickness_m.len()]
        .iter()
        .any(|&l| l != n)
    {
        return Err("rift_pass: field sizes differ from the grid".into());
    }
    let radius = rift_radius_m(grid, params);
    let mut rep = RiftReport::default();
    let mut thinned: BTreeMap<usize, f64> = BTreeMap::new();
    // SEED-PERF (lane SP, 2026-10-13): the births' discs depend only on
    // the grid, so they are computed up front in parallel and consumed in
    // ascending birth order (bit for bit; the first error where the serial
    // loop met it).
    let births: Vec<usize> = (0..n).filter(|&b| birth[b]).collect();
    let mut discs = crate::kinematic::par::map_heavy(births.len(), |k| {
        disc_with_distance(grid, cents, births[k], radius)
    })
    .into_iter();
    for &b in &births {
        let disc = discs.next().expect("one disc per birth")?;
        // Nearest non-birth cell of every plate in reach.
        let mut nearest: BTreeMap<PlateId, (f64, usize)> = BTreeMap::new();
        for &(c, d) in &disc {
            if birth[c] {
                continue;
            }
            let e = nearest.entry(plate[c]).or_insert((d, c));
            if d < e.0 || (d == e.0 && c < e.1) {
                *e = (d, c);
            }
        }
        let mut conj: Vec<(f64, PlateId, usize)> = nearest
            .iter()
            .filter(|(_, (_, c))| kind[*c] == CrustKindCodeV0::Continental)
            .map(|(p, (d, c))| (*d, *p, *c))
            .collect();
        if conj.len() < 2 {
            continue;
        }
        conj.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let (_, pa, ca) = conj[0];
        let (_, pb, cb) = conj[1];
        let (Some(wa), Some(wb)) = (poles.get(&pa), poles.get(&pb)) else {
            continue;
        };
        // Velocity of plate a relative to plate b at the birth, against the
        // direction from b's flank to a's flank: positive = opening.
        let v = relative_velocity_at(cents[b], wa, wb).map_err(|e| e.to_string())?;
        let d = [
            cents[ca][0] - cents[cb][0],
            cents[ca][1] - cents[cb][1],
            cents[ca][2] - cents[cb][2],
        ];
        if !(v[0] * d[0] + v[1] * d[1] + v[2] * d[2] > 0.0) {
            rep.convergent_skipped += 1;
            continue;
        }
        // Flanks: continental non-birth cells of the two plates in reach.
        let flanks: Vec<(usize, f64)> = disc
            .iter()
            .filter(|(c, _)| {
                !birth[*c]
                    && kind[*c] == CrustKindCodeV0::Continental
                    && (plate[*c] == pa || plate[*c] == pb)
            })
            .map(|&(c, dist)| (c, 1.0 - dist / radius))
            .filter(|(_, w)| *w > 0.0)
            .collect();
        let (mut s1, mut s2, mut wsum) = (0.0, 0.0, 0.0);
        for &(c, w) in &flanks {
            s1 += w * thickness_m[c] * areas[c];
            s2 += w * w * thickness_m[c] * areas[c];
            wsum += w * areas[c];
        }
        if !(s1 > 0.0) {
            continue;
        }
        let ab = areas[b];
        // Share alpha so the birth ends as thick as the weighted mean of its
        // thinned flanks: alpha s1 / A_b = (s1 - alpha s2) / W. No flank is
        // stretched below the breakup thickness (crust that thin has broken
        // up; the stretching moves on), so a flank gives at most its excess.
        let alpha = s1 * ab / (s1 * wsum + s2 * ab);
        let draws: Vec<(usize, f64)> = flanks
            .iter()
            .map(|&(c, w)| {
                let excess = (thickness_m[c] - params.breakup_thickness_m).max(0.0);
                (c, (alpha * w * thickness_m[c]).min(excess))
            })
            .collect();
        let moved: f64 = draws.iter().map(|&(c, dh)| dh * areas[c]).sum();
        if moved / ab < params.breakup_thickness_m {
            rep.breakup_births += 1;
            continue;
        }
        for &(c, dh) in &draws {
            thickness_m[c] -= dh;
            if dh > 0.0 {
                *thinned.entry(c).or_insert(0.0) += dh;
            }
        }
        kind[b] = CrustKindCodeV0::Continental;
        thickness_m[b] = moved / ab;
        rep.moved_m3 += moved;
        rep.births.push(b as u32);
        rep.birth_pre_thickness_m.push(s1 / wsum);
        rep.birth_flanks
            .push(flanks.iter().map(|&(c, w)| (c as u32, w * areas[c])).collect());
    }
    rep.thinning = thinned.into_iter().map(|(c, dh)| (c as u32, dh)).collect();
    Ok(rep)
}

/// One plume-head plate split (`plume-head-continental-split-v1`).
#[derive(Debug, Clone, PartialEq)]
pub struct PlateSplit {
    pub t_my: f64,
    pub plume: u64,
    pub parent: PlateId,
    pub child: PlateId,
    /// Unit normal of the cut plane (the child is on its positive side).
    pub cut_normal: [f64; 3],
    pub cells_moved: u64,
    pub child_pole: EulerPoleV0,
}

/// Unit normal of the great-circle cut through `p` with azimuth `az` (from
/// local north toward east).
pub fn cut_normal(p: [f64; 3], az: f64) -> [f64; 3] {
    let z = if p[2].abs() < 0.999_999 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let unit = |v: [f64; 3]| {
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / l, v[1] / l, v[2] / l]
    };
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let east = unit(cross(z, p));
    let north = cross(p, east);
    let (s, c) = az.sin_cos();
    let t = [
        c * north[0] + s * east[0],
        c * north[1] + s * east[1],
        c * north[2] + s * east[2],
    ];
    unit(cross(p, t))
}

/// Split plate `parent` along the great circle through `p` with cut
/// azimuth `az`: its cells (by material position) on the positive side of
/// the cut normal join a new plate (id = max id + 1, lineage `parent`), whose
/// angular velocity is the parent's plus `opening_rate_deg_per_my` about
/// `p x n` (so at `p` the child moves away from the parent along `n`).
/// Returns `None` (nothing changed) when one side would be empty.
#[allow(clippy::too_many_arguments)]
pub fn split_plate(
    plates: &mut Vec<PlateRecordV0>,
    plate_field: &mut [PlateId],
    material_pos: &[[f64; 3]],
    parent: PlateId,
    p: [f64; 3],
    az: f64,
    opening_rate_deg_per_my: f64,
    t_my: f64,
    plume: u64,
) -> Result<Option<PlateSplit>, String> {
    let n = cut_normal(p, az);
    let side: Vec<usize> = (0..plate_field.len())
        .filter(|&c| plate_field[c] == parent)
        .filter(|&c| {
            let x = material_pos[c];
            x[0] * n[0] + x[1] * n[1] + x[2] * n[2] > 0.0
        })
        .collect();
    let total = plate_field.iter().filter(|&&q| q == parent).count();
    if side.is_empty() || side.len() == total {
        return Ok(None);
    }
    let rec = plates
        .iter()
        .find(|r| r.plate_id == parent)
        .ok_or_else(|| format!("split: unknown plate {}", parent.0))?;
    let pole = rec
        .euler_pole
        .ok_or_else(|| format!("split: plate {} has no pole", parent.0))?;
    let w = angular_velocity(&pole).map_err(|e| e.to_string())?;
    let u = [
        p[1] * n[2] - p[2] * n[1],
        p[2] * n[0] - p[0] * n[2],
        p[0] * n[1] - p[1] * n[0],
    ];
    let om = opening_rate_deg_per_my.to_radians();
    let wc = [w[0] + om * u[0], w[1] + om * u[1], w[2] + om * u[2]];
    let mag = (wc[0] * wc[0] + wc[1] * wc[1] + wc[2] * wc[2]).sqrt();
    let child_pole = if mag > 1.0e-15 {
        let (lat, lon) =
            unit_to_lat_lon([wc[0] / mag, wc[1] / mag, wc[2] / mag]).map_err(|e| e.to_string())?;
        EulerPoleV0 {
            lat_deg: lat,
            lon_deg: lon,
            rate_deg_per_my: mag.to_degrees(),
        }
    } else {
        EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: 0.0,
        }
    };
    let child = PlateId(plates.iter().map(|r| r.plate_id.0).max().unwrap_or(0) + 1);
    for &c in &side {
        plate_field[c] = child;
    }
    plates.push(PlateRecordV0 {
        plate_id: child,
        parent_plate: Some(parent),
        birth_t_my: t_my,
        euler_pole: Some(child_pole),
        retirement: None,
    });
    Ok(Some(PlateSplit {
        t_my,
        plume,
        parent,
        child,
        cut_normal: n,
        cells_moved: side.len() as u64,
        child_pole,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundaries::relative_velocity_at;
    use crate::rotation::lat_lon_to_unit;

    #[test]
    fn cut_normal_is_perpendicular_to_the_plume() {
        for (lat, lon, az) in [(0.0, 0.0, 0.3), (45.0, 120.0, 2.0), (89.9999999, 0.0, 1.0)] {
            let p = lat_lon_to_unit(lat, lon).unwrap();
            let n = cut_normal(p, az);
            let dot = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(dot.abs() < 1e-12 && (len - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn split_opens_at_the_plume() {
        // Ten cells along the equator, all plate 0; plume at lon 0.
        let pos: Vec<[f64; 3]> = (0..10)
            .map(|k| lat_lon_to_unit(0.0, -45.0 + 10.0 * k as f64).unwrap())
            .collect();
        let mut field = vec![PlateId(0); 10];
        let mut plates = vec![PlateRecordV0 {
            plate_id: PlateId(0),
            parent_plate: None,
            birth_t_my: 0.0,
            euler_pole: Some(EulerPoleV0 {
                lat_deg: 10.0,
                lon_deg: 30.0,
                rate_deg_per_my: 0.7,
            }),
            retirement: None,
        }];
        let p = lat_lon_to_unit(0.0, 0.0).unwrap();
        // Azimuth 0: the cut runs north-south, normal points east or west.
        let s = split_plate(&mut plates, &mut field, &pos, PlateId(0), p, 0.0, 0.25, 5.0, 3)
            .unwrap()
            .unwrap();
        assert_eq!(s.child, PlateId(1));
        assert_eq!(plates.len(), 2);
        assert_eq!(plates[1].parent_plate, Some(PlateId(0)));
        assert_eq!(s.cells_moved, field.iter().filter(|&&q| q == PlateId(1)).count() as u64);
        assert!(s.cells_moved > 0 && s.cells_moved < 10);
        // Child moves away from the parent along n at p, at 0.25 deg/My.
        let v = relative_velocity_at(p, &s.child_pole, &plates[0].euler_pole.unwrap()).unwrap();
        let vn = v[0] * s.cut_normal[0] + v[1] * s.cut_normal[1] + v[2] * s.cut_normal[2];
        assert!((vn - 0.25f64.to_radians()).abs() < 1e-12, "{vn}");
        // Every moved cell is on the positive side of the cut.
        for (c, q) in field.iter().enumerate() {
            let x = pos[c];
            let side = x[0] * s.cut_normal[0] + x[1] * s.cut_normal[1] + x[2] * s.cut_normal[2];
            assert_eq!(*q == PlateId(1), side > 0.0);
        }
    }
}
