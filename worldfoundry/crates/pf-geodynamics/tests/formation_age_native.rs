//! Small native solver fixture, not a world run or historical C7 recapture.
//! Seeded plates supply every location/front. Broad back-arc bands are a
//! declared low-resolution fixture law, never hand-placed geographic features.
use pf_contracts::params::{FidelityPreset, MantleBlock, ScenarioPreset, SeedPath};
use pf_contracts::terrain::CrustKindCodeV0;
use pf_geodynamics::arcs::{find_fronts, ArcParams};
use pf_geodynamics::backarc::{
    backarc_pass_episode, backarc_pass_episode_with_formation_age_v1, BackarcParams,
};
use pf_geodynamics::boundaries::TransformParams;
use pf_geodynamics::carry::CarryLaw;
use pf_geodynamics::conserve::{ParcelLayer, ParcelParams};
use pf_geodynamics::formation_age::{BirthOrigin, EventKind, Formation, FormationAgeRecorderV1};
use pf_geodynamics::kinematic::{
    advance_with, advance_with_formation_age_v1, advance_with_parcels_transform, CrustState,
    ExtrasLaws, RidgeAccretion, NO_SOURCE,
};
use pf_geodynamics::orogen::{BeltParams, CollisionSedimentLaw};
use pf_geodynamics::plates::seed_plates;
use pf_geodynamics::uplift::UpliftPrior;
use pf_grid::GridSpec;
use serde_json::json;

fn seed() -> (
    GridSpec,
    Vec<pf_contracts::domain::PlateRecordV0>,
    CrustState,
    UpliftPrior,
) {
    let grid = GridSpec::new(18, 36, 6.371e6).unwrap();
    let mut profile =
        ScenarioPreset::toy().profile(FidelityPreset::Preview, MantleBlock::default());
    profile.tectonics.plate_count = 5;
    profile.tectonics.continental_seeds = 2;
    let path = SeedPath {
        seed: 240907,
        derivation: "kinematic-test".into(),
        salt: "k".into(),
    };
    let seeding = seed_plates(&grid, &profile, &path).unwrap();
    let prior = UpliftPrior::preview_v0();
    let state = CrustState::seed(&grid, &seeding, &prior).unwrap();
    (grid, seeding.plates, state, prior)
}
fn belt() -> BeltParams {
    BeltParams {
        shear_zone_dip_deg: 30.0,
        channel_viscosity_pa_s: 1e19,
        channel_reference_thickness_m: 35_000.0,
        delamination_thickness_m: 70_000.0,
        arc_productivity_m: 500.0,
        arc_belt_width_m: 100_000.0,
        gravity_m_s2: 9.81,
        crust_density_kg_m3: 2800.0,
        mantle_density_kg_m3: 3300.0,
    }
}
fn arc() -> ArcParams {
    ArcParams {
        slab_depth_m: 100_000.0,
        dip_young_deg: 30.0,
        dip_old_deg: 60.0,
        dip_age_my: 100.0,
        productivity_m: 500.0,
        rate_floor_m_per_my: 0.0,
        rate_cap_m_per_my: 200_000.0,
        half_width_m: 100_000.0,
        continental_thickness_m: 30_000.0,
        thermal_uplift_m: 1000.0,
        thermal_sigma_m: 100_000.0,
        thermal_decay_my: 10.0,
        crust_cap_m: 70_000.0,
    }
}

#[test]
fn native_fixture_births_projection_cancelled_gaps_and_thermal_resets() {
    pf_geodynamics::kinematic::par::set_threads(2);
    let (grid, plates, mut state, prior) = seed();
    let n = grid.cells();
    let cents = pf_geodynamics::flood_basalt::cell_centroids(&grid).unwrap();
    let areas: Vec<_> = (0..n)
        .map(|c| grid.cell_area_m2(c / grid.nlon()).unwrap())
        .collect();
    let mut layer = ParcelLayer::seed(&grid, &state).unwrap();
    let mut recorder = FormationAgeRecorderV1::inherited(&grid, &state, 100_000).unwrap();
    assert!(recorder.columns().iter().all(
        |c| c.formation == Formation::UnknownInherited || c.formation == Formation::NotOceanic
    ));
    let belt = belt();
    let pp = ParcelParams {
        overlap_cell_fraction: 0.5,
        fragment_min_area_fraction: 0.001,
    };
    let tp = TransformParams::v1();
    let extras = ExtrasLaws::default();
    // Global fixture parameters. One grid spacing resolves the broad band;
    // fronts, target plates and thermal birth values are derived by native code.
    let ba = BackarcParams {
        slab_age_min_my: 40.0,
        advance_max_m_per_my: 0.0,
        rate_m_per_my: 200_000.0,
        axis_offset_m: 500_000.0,
        half_width_m: grid.radius_m() * std::f64::consts::PI / grid.nlat() as f64,
        min_crust_m: 7000.0,
    };
    let mut steps = Vec::new();
    let mut gaps = 0;
    let mut projections = 0;
    let mut cancelled = 0;
    let mut resets = 0;
    let mut known_resets = 0;
    let mut unknown_resets = 0;
    let mut gap_formation_resets = 0;
    let mut carried = 0;
    let mut key_off_reset_checked = false;
    for step in 0..30 {
        let old_columns = recorder.columns().to_vec();
        let old_events = recorder.events().len();
        let mut off_layer = layer.clone();
        let (mut next, ledger, thick, source) = advance_with_formation_age_v1(
            &grid,
            &state,
            &plates,
            &prior,
            2.0,
            Some(&belt),
            CarryLaw::ContinentalParcelsV1,
            CollisionSedimentLaw::default(),
            Some(&extras),
            Some((&mut layer, pp)),
            RidgeAccretion::SymmetricV1,
            None,
            Some(&tp),
            Some(&mut recorder),
        )
        .unwrap();
        // One key-off comparison against the original production wrapper.
        if step == 4 {
            let off = advance_with_parcels_transform(
                &grid,
                &state,
                &plates,
                &prior,
                2.0,
                Some(&belt),
                CollisionSedimentLaw::default(),
                &extras,
                &mut off_layer,
                pp,
                RidgeAccretion::SymmetricV1,
                None,
                Some(&tp),
            )
            .unwrap();
            assert_eq!(
                (&next, &ledger, &thick, &source),
                (&off.0, &off.1, &off.2, &off.3)
            );
            assert_eq!(layer, off_layer);
        }
        let receipt = recorder.last_step().clone();
        let rep = ledger.parcels.as_ref().unwrap();
        let filled =
            ledger.transform.filled_transform_cells + ledger.transform.filled_convergent_cells;
        assert_eq!(
            receipt.cancelled_gap_cells as u64,
            filled + rep.overridden_births
        );
        assert_eq!(receipt.projection_birth_cells as u64, rep.vacated_births);
        assert_eq!(
            (receipt.accepted_gap_cells + receipt.projection_birth_cells) as u64,
            ledger.ridge_cells
        );
        gaps += receipt.accepted_gap_cells;
        projections += receipt.projection_birth_cells;
        cancelled += receipt.cancelled_gap_cells;
        for c in 0..n {
            if source[c] != NO_SOURCE
                && next.kind[c] == CrustKindCodeV0::Oceanic
                && state.kind[source[c] as usize] == CrustKindCodeV0::Oceanic
            {
                assert_eq!(recorder.columns()[c], old_columns[source[c] as usize]);
                carried += usize::from(matches!(
                    recorder.columns()[c].formation,
                    Formation::Recorded(_)
                ));
            }
        }
        let mut volume = 0.0;
        for e in &recorder.events()[old_events..] {
            assert_eq!(e.interval.start_my, state.t_my);
            assert_eq!(e.interval.end_my, next.t_my);
            assert_eq!(source[e.cell as usize], NO_SOURCE);
            assert_eq!(next.kind[e.cell as usize], CrustKindCodeV0::Oceanic);
            if let EventKind::Birth { created_m3, .. } = e.kind {
                assert!(created_m3.is_finite() && created_m3 > 0.0);
                assert!(
                    (created_m3 - areas[e.cell as usize] * prior.record.oceanic_crust_thickness_m)
                        .abs()
                        <= 1e-4 * created_m3
                );
                volume += created_m3;
            } else {
                panic!("carry emitted a thermal reset")
            }
        }
        assert!((volume - receipt.created_m3).abs() <= 1e-4 * volume.max(1.0));
        // The native continental parcel accounting remains supported separately.
        assert!(
            (layer.totals().1 + rep.consumed_m3 - rep.volume_before_m3).abs()
                <= 1e-4 * rep.volume_before_m3.max(1.0)
        );
        let fronts = find_fronts(&grid, &cents, &areas, &next, &arc()).unwrap();
        let before_thermal = recorder.columns().to_vec();
        let event_start = recorder.events().len();
        let mut off_state = next.clone();
        let volume_before = next.crust_volume_m3(&areas);
        let t_new = next.t_my;
        let mut dh = vec![0.0; n];
        let mut markers = vec![-1.0; n];
        let report = backarc_pass_episode_with_formation_age_v1(
            &grid,
            &cents,
            &areas,
            &mut next,
            &plates,
            &fronts,
            &ba,
            arc().half_width_m,
            arc().crust_cap_m,
            t_new,
            2.0,
            &mut dh,
            &mut markers,
            None,
            Some(&mut recorder),
        )
        .unwrap();
        if !key_off_reset_checked && report.renewed_cells > 0 {
            let mut off_dh = vec![0.0; n];
            let mut off_markers = vec![-1.0; n];
            let off_report = backarc_pass_episode(
                &grid,
                &cents,
                &areas,
                &mut off_state,
                &plates,
                &fronts,
                &ba,
                arc().half_width_m,
                arc().crust_cap_m,
                t_new,
                2.0,
                &mut off_dh,
                &mut off_markers,
                None,
            )
            .unwrap();
            assert_eq!(
                (&next, &report, &dh, &markers),
                (&off_state, &off_report, &off_dh, &off_markers)
            );
            key_off_reset_checked = true;
        }
        assert!((next.crust_volume_m3(&areas) - volume_before).abs() <= 1e-4 * volume_before);
        assert_eq!(
            recorder.events().len() - event_start,
            report.renewed_cells as usize
        );
        for e in &recorder.events()[event_start..] {
            assert_eq!(e.interval.start_my, t_new - 2.0);
            assert_eq!(e.interval.end_my, t_new);
            let c = e.cell as usize;
            assert_eq!(recorder.columns()[c].formation, before_thermal[c].formation);
            match e.kind {
                EventKind::ThermalReset {
                    formation,
                    previous_thermal_birth_my,
                    effective_thermal_birth_my,
                    affected_m3,
                } => {
                    assert!(effective_thermal_birth_my > previous_thermal_birth_my);
                    assert!(affected_m3.is_finite() && affected_m3 > 0.0);
                    assert!(t_new - effective_thermal_birth_my >= 0.5 - 1e-9);
                    known_resets += usize::from(matches!(formation, Formation::Recorded(_)));
                    unknown_resets += usize::from(formation == Formation::UnknownInherited);
                    if let Formation::Recorded(id) = formation {
                        gap_formation_resets += usize::from(matches!(
                            recorder.events()[id as usize - 1].kind,
                            EventKind::Birth {
                                origin: BirthOrigin::KinematicGap,
                                ..
                            }
                        ));
                    }
                }
                _ => panic!("thermal reset minted material"),
            }
        }
        resets += report.renewed_cells as usize;
        steps.push(json!({"end_my":t_new, "gap_births":receipt.accepted_gap_cells,
            "projection_births":receipt.projection_birth_cells, "cancelled_gaps":receipt.cancelled_gap_cells,
            "thermal_resets":report.renewed_cells, "birth_volume_m3":receipt.created_m3,
            "native_ridge_source_m3":ledger.ridge_created_m3, "projection_copied_cells":rep.vacated_copied,
            "native_remap_residual_m3":ledger.residual_m3(),
            "native_volume_before_m3":ledger.volume_before_m3}));
        state = next;
    }
    println!("gaps={gaps} projection_births={projections} cancelled={cancelled} resets={resets} known_resets={known_resets} unknown_resets={unknown_resets} carried_known={carried}");
    assert!(gaps > 0 && projections > 0 && cancelled > 0 && resets > 0 && carried > 0);
    assert!(
        gap_formation_resets > 0 && known_resets > 0 && unknown_resets > 0 && key_off_reset_checked
    );
    assert!(recorder.events().iter().any(|e| matches!(
        e.kind,
        EventKind::Birth {
            origin: BirthOrigin::ProjectionInitialization,
            ..
        }
    )));
    if let Ok(dir) = std::env::var("VZ_FORMATION_AGE_EVIDENCE") {
        std::fs::create_dir_all(&dir).unwrap();
        let path = std::path::Path::new(&dir);
        std::fs::write(
            path.join("native_events.json"),
            serde_json::to_vec_pretty(&recorder.to_json()).unwrap(),
        )
        .unwrap();
        std::fs::write(
            path.join("fixture_summary.json"),
            serde_json::to_vec_pretty(&json!({
            "fixture":"seeded native 18x36, seed 240907, 30 steps x 2 My; not C7/C8",
            "gaps":gaps, "projection_births":projections, "cancelled_gaps":cancelled,
            "thermal_resets":resets, "known_formation_resets":known_resets,
            "unknown_formation_resets":unknown_resets, "gap_formation_resets":gap_formation_resets, "carried_known_columns":carried,
            "steps":steps}))
            .unwrap(),
        )
        .unwrap();
    }
}

#[test]
fn disabled_and_full_recorder_keep_legacy_state_and_fail_before_mutation() {
    let (grid, plates, state, prior) = seed();
    let old = advance_with(&grid, &state, &plates, &prior, 2.0, None).unwrap();
    let off = advance_with_formation_age_v1(
        &grid,
        &state,
        &plates,
        &prior,
        2.0,
        None,
        CarryLaw::default(),
        CollisionSedimentLaw::default(),
        None,
        None,
        RidgeAccretion::PreviousOccupantV0,
        None,
        None,
        None,
    )
    .unwrap();
    assert_eq!(old, off);
    let mut recorder = FormationAgeRecorderV1::inherited(&grid, &state, 1).unwrap();
    let unchanged = recorder.clone();
    let err = advance_with_formation_age_v1(
        &grid,
        &state,
        &plates,
        &prior,
        2.0,
        None,
        CarryLaw::default(),
        CollisionSedimentLaw::default(),
        None,
        None,
        RidgeAccretion::PreviousOccupantV0,
        None,
        None,
        Some(&mut recorder),
    )
    .unwrap_err();
    assert!(err.to_string().contains("event limit"));
    assert_eq!(recorder, unchanged);
}
