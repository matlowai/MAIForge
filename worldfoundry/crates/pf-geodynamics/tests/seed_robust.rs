//! SEED-ROBUST (lane R5, 2026-10-12): new plate contacts that form during
//! a step. Before this node the kernel refused any subduction overlap
//! whose pair had no step-start interface within the swept locality bound
//! ("no input interface ... event-time geometry required" or "pair
//! interface ... is nonlocal"). Every refusal measured on the toy preset
//! was two plates converging onto a cell across a third plate's sliver or
//! a detached fragment, so the pair met during the step. The kernel now
//! samples that new contact across the material that met (winner source to
//! loser source) and keeps a typed refusal for sources that are not local.
//!
//! Expected values here are built by hand from the rigid-rotation geometry,
//! never read back from the kernel's own bookkeeping.

use pf_contracts::domain::{EulerPoleV0, PlateId, PlateRecordV0};
use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath};
use pf_contracts::terrain::{BoundaryClassCodeV0, CrustKindCodeV0};
use pf_geodynamics::kinematic::{advance, ArcContactGeometry, CrustState};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::rotation::lat_lon_to_unit;
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;

const R: f64 = 6.371e6;

fn areas(grid: &GridSpec) -> Vec<f64> {
    let mut a = Vec::new();
    for i in 0..grid.nlat() {
        let x = grid.cell_area_m2(i).unwrap();
        a.extend(std::iter::repeat(x).take(grid.nlon()));
    }
    a
}

fn plate(id: u32, rate: f64) -> PlateRecordV0 {
    PlateRecordV0 {
        plate_id: PlateId(id),
        parent_plate: None,
        birth_t_my: 0.0,
        euler_pole: Some(EulerPoleV0 {
            lat_deg: 90.0,
            lon_deg: 0.0,
            rate_deg_per_my: rate,
        }),
        retirement: None,
    }
}

/// Three oceanic plates in longitude strips: plate 0 west of -5 deg, a
/// one-column sliver plate 1 on [-5, 0), plate 2 east of 0. Plates 0 and 2
/// spin about the geographic pole toward each other at `r` deg/My each; the
/// sliver is fixed. Plates 0 and 2 also touch across the antimeridian, so
/// the pair HAS a step-start interface, about 180 deg from the sliver.
fn sliver_world(r: f64) -> (GridSpec, Vec<PlateRecordV0>, CrustState) {
    let grid = GridSpec::new(36, 72, R).unwrap();
    let n = grid.cells();
    let mut st = CrustState {
        t_my: 0.0,
        plate: vec![PlateId(0); n],
        kind: vec![CrustKindCodeV0::Oceanic; n],
        birth_t_my: vec![0.0; n],
        thickness_m: vec![7_000.0; n],
        sediment_m: vec![0.0; n],
        eroded_total_m: vec![0.0; n],
        arc_added_m: vec![0.0; n],
        shortening_added_m: vec![0.0; n],
        boundary_class: vec![BoundaryClassCodeV0::Interior.code(); n],
        normal_rate_m_per_my: vec![0.0; n],
        material_pos: vec![[0.0; 3]; n],
    };
    for i in 0..grid.nlat() {
        let lat = grid.center_lat_deg(i).unwrap();
        for j in 0..grid.nlon() {
            let lon = grid.center_lon_deg(j).unwrap();
            let c = i * grid.nlon() + j;
            st.material_pos[c] = lat_lon_to_unit(lat, lon).unwrap();
            // Declared ages: plate 0 youngest, so it overrides both others.
            let (p, birth) = if lon < -5.0 {
                (0, -10.0)
            } else if lon < 0.0 {
                (1, -100.0)
            } else {
                (2, -50.0)
            };
            st.plate[c] = PlateId(p);
            st.birth_t_my[c] = birth;
        }
    }
    (grid, vec![plate(0, r), plate(1, 0.0), plate(2, -r)], st)
}

#[test]
fn sliver_closure_forms_a_new_contact_sampled_across_the_met_material() {
    let r = 3.0; // deg/My: 0.6 of a 5 deg column per 1 My step
    let dt = 1.0;
    let (grid, plates, st) = sliver_world(r);
    let prior = UpliftPrior::preview_v0();
    let a = areas(&grid);
    // Before SEED-ROBUST this step refused: the only 0/2 interface is at
    // the antimeridian, about pi rad from the sliver.
    let (out, ledger, _) = advance(&grid, &st, &plates, &prior, dt).unwrap();
    assert_eq!(ledger.arc_normals.len(), ledger.arc_contact_geometry.len());
    let nlon = grid.nlon();
    let sliver_col = 35;
    let mut bridged_rows = 0;
    let mut expected_arc_m3 = 0.0;
    let mut expected_subducted_m3 = 0.0;
    for (s, g) in ledger.arc_normals.iter().zip(&ledger.arc_contact_geometry) {
        let c = s.overlap_element.0 as usize;
        expected_arc_m3 += prior.arc_thickening_m(s.normal_rate_m_per_my.max(0.0), dt) * a[c];
        expected_subducted_m3 += 7_000.0 * a[c];
        let (i, j) = (c / nlon, c % nlon);
        let pair = (s.winner.0, s.loser.0);
        if *g == ArcContactGeometry::NewContactBridge {
            // Only the 0/2 pair, only in the sliver column.
            assert_eq!(pair, (0, 2), "unexpected bridged pair at cell {c}");
            assert_eq!(j, sliver_col, "bridge outside the sliver at cell {c}");
        }
        let lat = grid.center_lat_deg(i).unwrap();
        if pair == (0, 2) && lat.abs() < 60.0 {
            assert_eq!(*g, ArcContactGeometry::NewContactBridge, "cell {c}");
            // Sources: plate 0 one column west, plate 2 one column east.
            let rec = [
                s.input_interface[0].0 as usize,
                s.input_interface[1].0 as usize,
            ];
            assert_eq!(rec, [c - 1, c + 1], "cell {c}");
            // Hand rate: the normal across the sources is due east at their
            // midpoint, whose latitude is atan(tan(lat) / cos(5 deg)); the
            // relative angular velocity is 2r about the pole, so the closing
            // speed is 2r cos(lat_mid) R.
            let lat_mid = (lat.to_radians().tan() / 5f64.to_radians().cos()).atan();
            let want = 2.0 * r.to_radians() * lat_mid.cos() * R;
            assert!(
                (s.normal_rate_m_per_my - want).abs() < 1e-9 * want,
                "cell {c}: {} vs {want}",
                s.normal_rate_m_per_my
            );
            bridged_rows += 1;
        }
        if pair == (0, 1) && lat.abs() < 60.0 {
            // The sliver still touches plate 0: the D5 step-start rule.
            assert_eq!(*g, ArcContactGeometry::StepStartInterface, "cell {c}");
        }
    }
    // Every row with |lat| < 60 deg (24 of 36) closed the sliver this step.
    assert_eq!(bridged_rows, 24);
    // Bookkeeping: arc volume and subducted volume are exactly the sum
    // over the recorded samples (one loser column retired per sample).
    assert!((ledger.arc_created_m3 - expected_arc_m3).abs() <= 1e-9 * expected_arc_m3);
    assert!((ledger.subducted_m3 - expected_subducted_m3).abs() <= 1e-9 * expected_subducted_m3);
    // The winner holds the sliver column where it closed.
    for i in 0..grid.nlat() {
        if grid.center_lat_deg(i).unwrap().abs() < 60.0 {
            assert_eq!(out.plate[i * nlon + sliver_col], PlateId(0));
        }
    }
}

#[test]
fn sliver_closure_is_deterministic() {
    let (grid, plates, st) = sliver_world(3.0);
    let prior = UpliftPrior::preview_v0();
    let x = advance(&grid, &st, &plates, &prior, 1.0).unwrap();
    let y = advance(&grid, &st, &plates, &prior, 1.0).unwrap();
    assert_eq!(x, y);
}

/// The first measured refusal: toy preset, seed 3, L0 (144 x 288), legacy
/// crust. The tectonic-only replay reaches the same state as `pf-lab
/// simulate` (erosion changes thickness and sediment only; claims depend on
/// plate, kind, age and material position). Pre-fix: refused at t = 16 My,
/// overlap 11334, plates 3 / 6, "no input interface". Plate 1's tongue,
/// one cell wide, separated them; both reached its tip in one step.
#[test]
fn seed3_l0_tongue_closure_completes_and_closes_its_ledgers() {
    let p = ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    let grid = GridSpec::level(p.fidelity.grid_level, p.planet.radius_m).unwrap();
    let sp = SeedPath {
        seed: 3,
        derivation: "kinematic-preview".into(),
        salt: "toy".into(),
    };
    let seeding = seed_plates(&grid, &p, &sp).unwrap();
    let prior = UpliftPrior::preview_v0();
    let a = areas(&grid);
    let run = || {
        let mut st = CrustState::seed(&grid, &seeding, &prior).unwrap();
        let mut ledgers = Vec::new();
        for _ in 0..9 {
            let (next, ledger, _) =
                advance(&grid, &st, &seeding.plates, &prior, p.tectonics.step_my).unwrap();
            ledgers.push(ledger);
            st = next;
        }
        (st, ledgers)
    };
    let (st, ledgers) = run();
    assert_eq!(st.t_my, 18.0);
    let last = ledgers.last().unwrap();
    let k = last
        .arc_normals
        .iter()
        .position(|s| s.overlap_element.0 == 11334)
        .expect("the refusing overlap is sampled");
    let s = &last.arc_normals[k];
    assert_eq!((s.winner.0, s.loser.0), (3, 6));
    assert_eq!(
        last.arc_contact_geometry[k],
        ArcContactGeometry::NewContactBridge
    );
    // Recorded sources: plate 3's from 11335, plate 6's from 11332.
    assert_eq!(
        [s.input_interface[0].0, s.input_interface[1].0],
        [11335, 11332]
    );
    assert!(s.normal_rate_m_per_my.is_finite() && s.normal_rate_m_per_my > 0.0);
    let bridges: usize = ledgers
        .iter()
        .map(|l| {
            l.arc_contact_geometry
                .iter()
                .filter(|g| **g == ArcContactGeometry::NewContactBridge)
                .count()
        })
        .sum();
    assert_eq!(bridges, 1, "only the tongue closure is a new contact");
    for l in &ledgers {
        assert_eq!(l.arc_normals.len(), l.arc_contact_geometry.len());
        let arc: f64 = l
            .arc_normals
            .iter()
            .map(|s| {
                prior.arc_thickening_m(s.normal_rate_m_per_my.max(0.0), 2.0)
                    * a[s.overlap_element.0 as usize]
            })
            .sum();
        assert!((l.arc_created_m3 - arc).abs() <= 1e-9 * arc.max(1.0));
        // Wedge volume closes (advance refuses otherwise; asserted here too)
        // and no sediment exists in a tectonic-only replay.
        assert!(
            (l.wedge_in_m3 - l.wedge_stacked_m3 - l.wedge_delaminated_m3).abs()
                <= 1e-12 * l.wedge_in_m3.max(1.0)
        );
        assert_eq!(l.sediment_residual_m3(), 0.0);
    }
    let (st2, ledgers2) = run();
    assert_eq!(st, st2);
    assert_eq!(ledgers, ledgers2);
}
