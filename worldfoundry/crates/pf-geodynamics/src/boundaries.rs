//! Boundary classification with stable IDs (brief K3/K5; acceptance:
//! relative-velocity classification; stable boundary IDs; no parcel
//! transport).
//!
//! Motion inputs are PLATE-LEVEL Euler poles (from plates.rs) plus the two
//! plates' reference points (their Voronoi seeds). At a candidate boundary
//! point r, the relative angular velocity produces a tangential relative
//! velocity v = (w_a - w_b) x r; the component along the local a->b
//! direction n̂ (tangent at r pointing from plate a's side toward plate b's)
//! is the NORMAL rate: positive = convergent (b approaching a's side
//! closing), negative = divergent. Classification thresholds are declared
//! constants in rad/My scaled by radius — physical rates, not grid units.
//!
//! Convergent boundaries split by crust typing (K5): any oceanic crust in
//! the pair -> Subduction; continent-continent -> Collision. The caller
//! supplies the typing (from the parcel ledger's plate flags — this module
//! never touches parcels).

use pf_contracts::domain::{BoundaryKindV0, BoundaryRecordV0, EulerPoleV0, PlateId};
use std::collections::BTreeMap;

/// Classification threshold on the normal rate (rad/My at unit radius):
/// |normal| <= this is transform-adjacent, not convergent/divergent.
/// Declared constant; tests pin it.
pub const NORMAL_RATE_THRESHOLD_RAD_PER_MY: f64 = 1.0e-4;

#[derive(Debug, Clone, PartialEq)]
pub enum BoundariesError {
    NonFinite {
        what: &'static str,
        value: f64,
    },
    UnknownPlate {
        plate: PlateId,
    },
    BadPair {
        a: PlateId,
        b: PlateId,
        why: &'static str,
    },
    EmptyFootprint {
        boundary: u64,
    },
    NonMonotonicTime {
        last: f64,
        got: f64,
    },
}

impl std::fmt::Display for BoundariesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BoundariesError::NonFinite { what, value } => write!(f, "nonfinite {what}: {value:?}"),
            BoundariesError::UnknownPlate { plate } => write!(f, "unknown plate {plate:?}"),
            BoundariesError::BadPair { a, b, why } => write!(f, "plate pair ({a:?}, {b:?}): {why}"),
            BoundariesError::EmptyFootprint { boundary } => {
                write!(f, "boundary {boundary} with empty footprint")
            }
            BoundariesError::NonMonotonicTime { last, got } => {
                write!(
                    f,
                    "boundary update at {got:?} My not after previous {last:?} My"
                )
            }
        }
    }
}

impl std::error::Error for BoundariesError {}

/// Angular velocity vector of an Euler pole: axis * rate in rad/My.
pub fn angular_velocity(pole: &EulerPoleV0) -> Result<[f64; 3], BoundariesError> {
    pole.validate().map_err(|_| BoundariesError::NonFinite {
        what: "euler pole",
        value: pole.lat_deg,
    })?;
    let axis = crate::rotation::euler_pole_axis(pole).map_err(|_| BoundariesError::NonFinite {
        what: "euler pole axis",
        value: pole.rate_deg_per_my,
    })?;
    let rate_rad = pole.rate_deg_per_my * std::f64::consts::PI / 180.0;
    Ok([axis[0] * rate_rad, axis[1] * rate_rad, axis[2] * rate_rad])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Tangential relative velocity at unit point `r` between two plates'
/// Euler rotations: (w_a - w_b) x r. Units: rad/My (angular displacement
/// of material at r per My, tangential).
pub fn relative_velocity_at(
    r: [f64; 3],
    a: &EulerPoleV0,
    b: &EulerPoleV0,
) -> Result<[f64; 3], BoundariesError> {
    let wa = angular_velocity(a)?;
    let wb = angular_velocity(b)?;
    let dw = [wa[0] - wb[0], wa[1] - wb[1], wa[2] - wb[2]];
    Ok(cross(dw, r))
}

/// Local a->b tangent direction at boundary point r: the projection of
/// (ref_b - ref_a) onto the tangent plane at r, normalized. Fails loudly
/// when the projection is degenerate (r on the ref-antipode great circle).
pub fn boundary_normal(
    r: [f64; 3],
    ref_a: [f64; 3],
    ref_b: [f64; 3],
) -> Result<[f64; 3], BoundariesError> {
    let d = [
        ref_b[0] - ref_a[0],
        ref_b[1] - ref_a[1],
        ref_b[2] - ref_a[2],
    ];
    let dr = dot(d, r);
    let t = [d[0] - dr * r[0], d[1] - dr * r[1], d[2] - dr * r[2]];
    let n = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
    if !n.is_finite() || n < 1.0e-9 {
        return Err(BoundariesError::NonFinite {
            what: "boundary normal norm",
            value: n,
        });
    }
    Ok([t[0] / n, t[1] / n, t[2] / n])
}

/// Signed normal rate at r: positive = convergent (closing), negative =
/// divergent (opening). Transform shear is the tangential remainder.
pub fn normal_rate_at(
    r: [f64; 3],
    ref_a: [f64; 3],
    ref_b: [f64; 3],
    a: &EulerPoleV0,
    b: &EulerPoleV0,
) -> Result<f64, BoundariesError> {
    let n = boundary_normal(r, ref_a, ref_b)?;
    let v = relative_velocity_at(r, a, b)?;
    Ok(dot(v, n))
}

/// K3/K5 classification from the normal rate and crust typing.
/// Divergent -> Ridge; Convergent + any oceanic crust -> Subduction;
/// Convergent + continent-continent -> Collision; |normal| below the
/// threshold with motion -> Transform.
pub fn classify(
    normal_rate_rad_per_my: f64,
    both_continental: bool,
) -> Result<BoundaryKindV0, BoundariesError> {
    if !normal_rate_rad_per_my.is_finite() {
        return Err(BoundariesError::NonFinite {
            what: "normal rate",
            value: normal_rate_rad_per_my,
        });
    }
    let thr = NORMAL_RATE_THRESHOLD_RAD_PER_MY;
    if normal_rate_rad_per_my < -thr {
        Ok(BoundaryKindV0::Ridge)
    } else if normal_rate_rad_per_my > thr {
        if both_continental {
            Ok(BoundaryKindV0::Collision)
        } else {
            Ok(BoundaryKindV0::Subduction)
        }
    } else {
        Ok(BoundaryKindV0::Transform)
    }
}

/// SEED-TRANSFORMS (lane TF, 2026-10-13): `smoothed-normal-transform-v1`.
///
/// Boundary classification v2. The boundary normal at a boundary cell is
/// the first moment of the two plates' indicators over a Gaussian window
/// (sigma `normal_smoothing_m`, cut at 3 sigma, cell-area weighted):
/// n ~ sum_k w_k (chi_b(k) - chi_a(k)) (r_k - (r_k . r_c) r_c), so it
/// follows the boundary trace instead of the grid staircase (lane VL
/// defect 1: the staircase direction flips the closing/opening sign on
/// 16-23% of moving boundary length). The relative velocity v = (w_a -
/// w_b) x r splits into a normal part v_n (positive = closing) and a
/// tangential part v_t; a boundary is a TRANSFORM when its obliquity
/// atan(|v_n| / |v_t|) is below `transform_obliquity_deg` (or |v_n| is
/// below the legacy threshold), else convergent / divergent by the sign of
/// v_n. In the carry (`kinematic.rs`), gaps on transform or convergent
/// boundaries are filled from the same side's carried material (no ridge
/// birth), and oceanic overlaps on transform or divergent boundaries
/// retire the loser without subduction (no arc sample).
pub const LAW_TRANSFORM_V1: &str = "smoothed-normal-transform-v1";

/// Window reach of the smoothed normal, in units of sigma (fixed by the
/// law version).
pub const TRANSFORM_WINDOW_SIGMAS: f64 = 3.0;
/// Rings of cells searched around a gap or overlap cell for the nearest
/// step-start boundary class (fixed by the law version).
pub const TRANSFORM_CLASS_SEARCH_CELLS: i64 = 2;
/// Gap-fill layers: cells filled in one layer can feed the next (fixed by
/// the law version).
pub const TRANSFORM_FILL_LAYERS: usize = 2;

/// PLATE-POLES (lane PLATE-POLES, 2026-10-21): `rigid-group-contacts-v1`,
/// opt-in on top of `smoothed-normal-transform-v1` and the suture law.
///
/// Why: plates welded into one suture group (`contact-length-pole-share-v1`
/// sutures, breakup sides, fragment pieces) carry bit-identical Euler poles,
/// so they are one rigid body; their mutual contacts have zero relative
/// velocity everywhere. The legacy classifier still labels those contacts
/// TRANSFORM (|v_n| below the threshold), so a fossil suture inside a
/// continent is published as an active strike-slip boundary (C6, seed 42,
/// t = 400 My: 1407 of 13732 boundary cells, all transform; lane
/// PLATE-POLES). A plate boundary needs relative motion; a contact inside
/// one rigid body is plate interior (fossil sutures keep their record in
/// the shortening and stratigraphy fields, not in `boundary_class`).
///
/// Two parts, two keys: `rigid-group-interior-v1` is part 1 alone (labels
/// and boundary records only; motion and material unchanged);
/// `rigid-group-contacts-v1` is parts 1 and 2.
/// 1. Classification (`classify_cells_v2` with `rigid_group_interior`):
///    neighbouring plates whose angular velocity is bit-identical to the
///    cell's own plate are not boundary partners; a cell whose only other
///    plates are such is Interior, and no boundary record is minted for the
///    pair. The carry's step-start remap class is NOT changed: a gap or
///    overlap beside a welded contact keeps the same-body remap rule (fill
///    from the same body, no ridge birth, no subduction), which is the
///    correct rule inside a rigid body.
/// 2. Release (`reorg::reorg_pass_release`): at each scheduled
///    reorganization event, after the fragment pass and before the
///    re-solve, member plates of a suture group that touch no other member
///    of the group's main body (largest-area 8-connected set of members)
///    leave the group as their own group(s). A weld needs contact; a piece
///    detached from its group (e.g. a `fragment-plate` piece, which the
///    legacy pass keeps in its parent's group) otherwise co-moves with it
///    at a distance forever. The re-solve in the same event then gives the
///    released group its own pole from the same seeded slab-torque draw as
///    every group (keyed SplitMix64 stream, ADR 044 determinism); no new
///    constant, no invented perturbation. Released pieces share no boundary
///    with their former group, so there is no new boundary whose opening
///    rate a perturbation would have to set.
pub const LAW_RIGID_GROUP_CONTACTS_V1: &str = "rigid-group-contacts-v1";
/// Part 1 alone (see `LAW_RIGID_GROUP_CONTACTS_V1`).
pub const LAW_RIGID_GROUP_INTERIOR_V1: &str = "rigid-group-interior-v1";

/// Constants of `smoothed-normal-transform-v1`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransformParams {
    /// Obliquity (degrees between the relative velocity and the boundary
    /// trace) below which a moving boundary is a transform.
    pub transform_obliquity_deg: f64,
    /// Gaussian sigma (m) of the boundary-normal smoothing window.
    pub normal_smoothing_m: f64,
    /// PLATE-POLES `rigid-group-interior-v1` (part 1 of
    /// `rigid-group-contacts-v1`): contacts between plates
    /// with bit-identical angular velocity are plate interior (false =
    /// legacy, bit for bit).
    pub rigid_group_interior: bool,
}

impl TransformParams {
    /// Pre-registered values (`unpublished-design-record`).
    pub fn v1() -> Self {
        Self {
            transform_obliquity_deg: 30.0,
            normal_smoothing_m: 200_000.0,
            rigid_group_interior: false,
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(self.transform_obliquity_deg.is_finite()
            && self.transform_obliquity_deg > 0.0
            && self.transform_obliquity_deg < 90.0)
        {
            return Err(format!(
                "{LAW_TRANSFORM_V1}: transform_obliquity_deg must be in (0, 90), got {}",
                self.transform_obliquity_deg
            ));
        }
        if !(self.normal_smoothing_m.is_finite() && self.normal_smoothing_m > 0.0) {
            return Err(format!(
                "{LAW_TRANSFORM_V1}: normal_smoothing_m must be finite and positive, got {}",
                self.normal_smoothing_m
            ));
        }
        Ok(())
    }
}

/// v2 classification from the normal and tangential relative rates
/// (rad/My at unit radius): transform when |v_n| is at or below the legacy
/// threshold or the obliquity atan(|v_n| / |v_t|) is below the declared
/// angle; else ridge / subduction / collision by the sign of v_n.
pub fn classify_oblique(
    normal_rate_rad_per_my: f64,
    tangential_rate_rad_per_my: f64,
    both_continental: bool,
    obliquity_deg: f64,
) -> Result<BoundaryKindV0, BoundariesError> {
    if !normal_rate_rad_per_my.is_finite() || !tangential_rate_rad_per_my.is_finite() {
        return Err(BoundariesError::NonFinite {
            what: "normal or tangential rate",
            value: normal_rate_rad_per_my,
        });
    }
    let obl = normal_rate_rad_per_my
        .abs()
        .atan2(tangential_rate_rad_per_my.abs())
        .to_degrees();
    if normal_rate_rad_per_my.abs() <= NORMAL_RATE_THRESHOLD_RAD_PER_MY || obl < obliquity_deg {
        return Ok(BoundaryKindV0::Transform);
    }
    classify(normal_rate_rad_per_my, both_continental)
}

/// Stable boundary table across steps: (ordered plate pair, kind) keeps
/// its ID; a kind change retires the old record and mints a new ID
/// (K7 lineage-friendly). Tracks the last update time.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryTable {
    next_id: u64,
    last_t_my: f64,
    /// Active boundaries keyed by ordered plate pair.
    active: BTreeMap<(PlateId, PlateId), BoundaryRecordV0>,
    /// Retired (kind-changed) records, kept for lineage.
    retired: Vec<BoundaryRecordV0>,
}

pub const BOUNDARY_ID_SEED: u64 = 1;

impl BoundaryTable {
    pub fn new(t_start_my: f64) -> Self {
        Self {
            next_id: BOUNDARY_ID_SEED,
            last_t_my: t_start_my,
            active: BTreeMap::new(),
            retired: Vec::new(),
        }
    }

    /// Record one boundary for this step. Same pair + same kind keeps the
    /// existing ID and refreshes activity/footprint; a kind change retires
    /// the old record and mints a new one. Footprint must be nonempty.
    /// Time is consumed only on success.
    pub fn upsert(
        &mut self,
        t_my: f64,
        a: PlateId,
        b: PlateId,
        kind: BoundaryKindV0,
        footprint_element_count: u64,
    ) -> Result<u64, BoundariesError> {
        let id = self.upsert_inner(t_my, a, b, kind, footprint_element_count)?;
        self.last_t_my = t_my;
        Ok(id)
    }

    fn upsert_inner(
        &mut self,
        t_my: f64,
        a: PlateId,
        b: PlateId,
        kind: BoundaryKindV0,
        footprint_element_count: u64,
    ) -> Result<u64, BoundariesError> {
        if !t_my.is_finite() || t_my < 0.0 {
            return Err(BoundariesError::NonFinite {
                what: "t_my",
                value: t_my,
            });
        }
        if t_my < self.last_t_my {
            return Err(BoundariesError::NonMonotonicTime {
                last: self.last_t_my,
                got: t_my,
            });
        }
        if footprint_element_count == 0 {
            return Err(BoundariesError::EmptyFootprint {
                boundary: self.next_id,
            });
        }
        let key = if a <= b { (a, b) } else { (b, a) };
        match self.active.get(&key) {
            Some(existing) if existing.kind == kind => {
                let id = existing.boundary_id;
                let rec = self.active.get_mut(&key).expect("checked");
                rec.last_active_t_my = t_my;
                rec.footprint_element_count = footprint_element_count;
                Ok(id)
            }
            Some(existing) => {
                let mut old = existing.clone();
                old.last_active_t_my = t_my;
                self.retired.push(old);
                let id = self.next_id;
                self.next_id += 1;
                self.active.insert(
                    key,
                    BoundaryRecordV0 {
                        boundary_id: id,
                        kind,
                        plate_a: key.0,
                        plate_b: key.1,
                        birth_t_my: t_my,
                        last_active_t_my: t_my,
                        footprint_element_count,
                    },
                );
                Ok(id)
            }
            None => {
                let id = self.next_id;
                self.next_id += 1;
                self.active.insert(
                    key,
                    BoundaryRecordV0 {
                        boundary_id: id,
                        kind,
                        plate_a: key.0,
                        plate_b: key.1,
                        birth_t_my: t_my,
                        last_active_t_my: t_my,
                        footprint_element_count,
                    },
                );
                Ok(id)
            }
        }
    }

    pub fn active(&self) -> impl Iterator<Item = &BoundaryRecordV0> {
        self.active.values()
    }

    pub fn retired(&self) -> &[BoundaryRecordV0] {
        &self.retired
    }

    pub fn get(&self, a: PlateId, b: PlateId) -> Option<&BoundaryRecordV0> {
        let key = if a <= b { (a, b) } else { (b, a) };
        self.active.get(&key)
    }
}

/// Convenience: normal rate between two plates at the midpoint direction
/// of their reference points (coarse B1 sampling point).
pub fn pair_normal_rate(
    ref_a: [f64; 3],
    ref_b: [f64; 3],
    a: &EulerPoleV0,
    b: &EulerPoleV0,
) -> Result<f64, BoundariesError> {
    let mid = [
        ref_a[0] + ref_b[0],
        ref_a[1] + ref_b[1],
        ref_a[2] + ref_b[2],
    ];
    let norm = (dot(mid, mid)).sqrt();
    if !norm.is_finite() || norm < 1.0e-12 {
        return Err(BoundariesError::BadPair {
            a: PlateId(0),
            b: PlateId(1),
            why: "antipodal references",
        });
    }
    let r = [mid[0] / norm, mid[1] / norm, mid[2] / norm];
    normal_rate_at(r, ref_a, ref_b, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rotation::lat_lon_to_unit;

    fn pole(lat: f64, lon: f64, rate: f64) -> EulerPoleV0 {
        EulerPoleV0 {
            lat_deg: lat,
            lon_deg: lon,
            rate_deg_per_my: rate,
        }
    }

    #[test]
    fn divergent_convergent_transform_classification() {
        // Clean geometry: plate a south of the boundary, plate b north;
        // boundary point on the equator at lon 0; a->b normal points +z.
        let ref_a = lat_lon_to_unit(-10.0, 0.0).unwrap();
        let ref_b = lat_lon_to_unit(10.0, 0.0).unwrap();
        let r = lat_lon_to_unit(0.0, 0.0).unwrap();
        let n = boundary_normal(r, ref_a, ref_b).unwrap();
        assert!(n[2] > 0.99, "normal must point north (+z), got {n:?}");

        // Poles along +/-y give vertical (normal) relative motion at r:
        // w = +y*rate rotates r=(1,0,0) toward -z; opposite signs diverge
        // or converge cleanly across the E-W boundary.
        let rate_rad = 0.01; // rad/My, well above threshold
        let p_y_pos = EulerPoleV0 {
            lat_deg: 0.0,
            lon_deg: 90.0,
            rate_deg_per_my: rate_rad * 180.0 / std::f64::consts::PI,
        };
        let p_y_neg = EulerPoleV0 {
            lat_deg: 0.0,
            lon_deg: 90.0,
            rate_deg_per_my: -rate_rad * 180.0 / std::f64::consts::PI,
        };

        let divergent = normal_rate_at(r, ref_a, ref_b, &p_y_pos, &p_y_neg).unwrap();
        assert!(
            divergent < -NORMAL_RATE_THRESHOLD_RAD_PER_MY,
            "expected opening, got {divergent}"
        );
        assert_eq!(classify(divergent, false).unwrap(), BoundaryKindV0::Ridge);

        let convergent = normal_rate_at(r, ref_a, ref_b, &p_y_neg, &p_y_pos).unwrap();
        assert!(
            convergent > NORMAL_RATE_THRESHOLD_RAD_PER_MY,
            "expected closing, got {convergent}"
        );
        assert_eq!(
            classify(convergent, false).unwrap(),
            BoundaryKindV0::Subduction
        );
        assert_eq!(
            classify(convergent, true).unwrap(),
            BoundaryKindV0::Collision
        );

        // Identical rotations: no relative motion -> transform band.
        let still = normal_rate_at(r, ref_a, ref_b, &p_y_pos, &p_y_pos).unwrap();
        assert!(still.abs() <= NORMAL_RATE_THRESHOLD_RAD_PER_MY);
        assert_eq!(classify(still, false).unwrap(), BoundaryKindV0::Transform);

        // The classic shear case: plates north/south of an E-W boundary
        // (normal +z) with opposite z-poles move boundary-PARALLEL (+/-y):
        // pure shear, normal rate ~0, transform.
        let ra = lat_lon_to_unit(-10.0, 0.0).unwrap();
        let rb = lat_lon_to_unit(10.0, 0.0).unwrap();
        let z_pos = EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: rate_rad * 180.0 / std::f64::consts::PI,
        };
        let z_neg = EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: -rate_rad * 180.0 / std::f64::consts::PI,
        };
        let shear_rate = normal_rate_at(r, ra, rb, &z_pos, &z_neg).unwrap();
        assert!(
            shear_rate.abs() <= NORMAL_RATE_THRESHOLD_RAD_PER_MY,
            "pure shear must classify transform, got {shear_rate}"
        );

        // Near-threshold classification is exact at the constant.
        assert_eq!(
            classify(NORMAL_RATE_THRESHOLD_RAD_PER_MY, false).unwrap(),
            BoundaryKindV0::Transform
        );
        assert_eq!(
            classify(NORMAL_RATE_THRESHOLD_RAD_PER_MY * 1.0001, false).unwrap(),
            BoundaryKindV0::Subduction
        );
        assert_eq!(
            classify(-NORMAL_RATE_THRESHOLD_RAD_PER_MY * 1.0001, false).unwrap(),
            BoundaryKindV0::Ridge
        );
        assert!(classify(f64::NAN, false).is_err());
    }

    #[test]
    fn oblique_classification_has_a_transform_class() {
        let obl = TransformParams::v1().transform_obliquity_deg;
        // Pure shear and slightly oblique motion: transform.
        assert_eq!(classify_oblique(0.0, 1.0e-2, false, obl).unwrap(), BoundaryKindV0::Transform);
        let v = 1.0e-2;
        let at = |deg: f64| (v * deg.to_radians().sin(), v * deg.to_radians().cos());
        let (n, t) = at(29.0);
        assert_eq!(classify_oblique(n, t, false, obl).unwrap(), BoundaryKindV0::Transform);
        assert_eq!(classify_oblique(-n, t, true, obl).unwrap(), BoundaryKindV0::Transform);
        // Beyond the obliquity: convergent / divergent by sign, typed by crust.
        let (n, t) = at(31.0);
        assert_eq!(classify_oblique(n, t, false, obl).unwrap(), BoundaryKindV0::Subduction);
        assert_eq!(classify_oblique(n, t, true, obl).unwrap(), BoundaryKindV0::Collision);
        assert_eq!(classify_oblique(-n, t, false, obl).unwrap(), BoundaryKindV0::Ridge);
        // Head-on motion is never a transform; no motion is (legacy band).
        assert_eq!(classify_oblique(v, 0.0, false, obl).unwrap(), BoundaryKindV0::Subduction);
        assert_eq!(classify_oblique(0.0, 0.0, false, obl).unwrap(), BoundaryKindV0::Transform);
        assert!(classify_oblique(f64::NAN, 0.0, false, obl).is_err());
        assert!(TransformParams::v1().validate().is_ok());
        assert!(TransformParams { transform_obliquity_deg: 90.0, ..TransformParams::v1() }.validate().is_err());
        assert!(TransformParams { normal_smoothing_m: 0.0, ..TransformParams::v1() }.validate().is_err());
    }

    #[test]
    fn boundary_ids_are_stable_and_kind_changes_mint_new_ids() {
        let mut t = BoundaryTable::new(0.0);
        let (a, b) = (PlateId(1), PlateId(2));
        let id1 = t.upsert(1.0, a, b, BoundaryKindV0::Ridge, 10).unwrap();
        let id2 = t.upsert(2.0, a, b, BoundaryKindV0::Ridge, 12).unwrap();
        assert_eq!(id1, id2, "same pair+kind keeps its id");
        assert_eq!(t.get(a, b).unwrap().footprint_element_count, 12);
        // Order-insensitive key.
        let id3 = t.upsert(3.0, b, a, BoundaryKindV0::Ridge, 9).unwrap();
        assert_eq!(id3, id1);
        // Kind change: new id, old retired with lineage kept.
        let id4 = t.upsert(4.0, a, b, BoundaryKindV0::Subduction, 8).unwrap();
        assert_ne!(id4, id1);
        assert_eq!(t.retired().len(), 1);
        assert_eq!(t.retired()[0].boundary_id, id1);
        assert_eq!(t.retired()[0].kind, BoundaryKindV0::Ridge);
        assert_eq!(t.get(a, b).unwrap().kind, BoundaryKindV0::Subduction);
        // Loud failures: empty footprint, non-monotonic time.
        assert!(t.upsert(5.0, a, b, BoundaryKindV0::Ridge, 0).is_err());
        assert!(t.upsert(3.5, a, b, BoundaryKindV0::Ridge, 5).is_err());
    }

    #[test]
    fn angular_velocity_matches_pole_geometry() {
        let w = angular_velocity(&pole(90.0, 0.0, 180.0 / std::f64::consts::PI)).unwrap();
        // rate = 1 rad/My about +z.
        assert!(
            (w[0] - 0.0).abs() < 1.0e-12 && (w[2] - 1.0).abs() < 1.0e-12,
            "{w:?}"
        );
        // Rigid sanity: a point's speed under one plate equals |w x r|.
        let r = lat_lon_to_unit(0.0, 90.0).unwrap();
        let v = relative_velocity_at(
            r,
            &pole(90.0, 0.0, 180.0 / std::f64::consts::PI),
            &EulerPoleV0 {
                lat_deg: 0.0,
                lon_deg: 0.0,
                rate_deg_per_my: 0.0,
            },
        )
        .unwrap();
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!(
            (speed - 1.0).abs() < 1.0e-12,
            "|w x r| must equal |w| for unit r, got {speed}"
        );
    }
}
