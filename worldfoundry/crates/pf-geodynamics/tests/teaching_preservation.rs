//! Authored preservation counterexamples through the existing public native API.
//! Fixture parameters follow the existing native unit-test constants; no production edits.
use pf_contracts::terrain::{STRATA_FLAG_MERGED, STRATA_FLAG_TRUNCATED_TOP};
use pf_geodynamics::stratigraphy::{apply_parcel, ParcelAmounts, StrataLayer,
    StrataParams, StrataStack, StrataStepLog, StrataV2Params,
    ENV_PLUTON, ENV_VOLCANIC, LITH_BASEMENT, LITH_BASALT, LITH_GRANITE};

fn params(capacity: usize) -> StrataParams {
        StrataParams {
            max_layers: capacity,
            merge_gap_my: 2.0,
            unconformity_gap_my: 5.0,
            min_layer_m: 10.0,
            trace_m: 0.01,
            arc_extrusive_fraction: 0.3,
            plume_extrusive_fraction: 0.2,
            shortening_pluton_fraction: 0.1,
            pluton_roof_m: 4_000.0,
            metamorphic_depth_m: 10_000.0,
            high_grade_depth_m: 20_000.0,
            shelf_depth_m: 200.0,
            carbonate_max_paleolat_deg: 30.0,
            carbonate_max_supply_m_per_my: 50.0,
            arid_lat_min_deg: 15.0,
            arid_lat_max_deg: 35.0,
            humid_lat_max_deg: 15.0,
            humid_midlat_min_deg: 40.0,
            humid_midlat_max_deg: 65.0,
            coal_max_elevation_m: 200.0,
            conglomerate_relief_m: 1000.0,
            river_mouth_area_m2: 1.0e11,
            publish_depth_m: 3000.0,
            basement_min_m: 5000.0,
            v2: Some(StrataV2Params {
                shelf_width_m: 80000.0, slope_foot_depth_m: 4000.0,
                shelf_majority_fraction: 0.5, platform_cover_max_m: 2000.0,
                platform_full_elev_m: 0.0, platform_zero_elev_m: 600.0,
                platform_sandstone_fraction: 0.3, platform_carbonate_fraction: 0.3,
                platform_age_base_my: -600.0, platform_age_top_my: -50.0,
            }),
        }
    }

fn history(capacity: usize, hiatus: bool) -> (StrataStack, StrataStepLog) {
    let p = params(capacity);
    p.validate().expect("valid fixture parameters");
    let mut stack = StrataStack::basement_only(35000.0, LITH_BASEMENT, -1000.0);
    // Explicit inherited cover, not a simulated prehistory. Older pairs cross gaps
    // and have higher merge cost than the subsequently deposited shale pair.
    for i in 0..14 {
        let granite = i % 2 == 0;
        let age = -300.0 + 12.0 * i as f64;
        stack.cover.push(StrataLayer {
            thickness_m: 100.0, age_base_my: age, age_top_my: age + 1.0,
            max_burial_m: 0.0, lith: if granite { LITH_GRANITE } else { LITH_BASALT },
            env: if granite { ENV_PLUTON } else { ENV_VOLCANIC }, grade: 0, flags: 0,
        });
    }
    let initial = stack.total_m();
    let mut log = StrataStepLog::default();
    let intervals = if hiatus { [(0.0, 2.0), (8.0, 10.0)] }
                    else { [(0.0, 5.0), (5.0, 10.0)] };
    for (index, (start, end)) in intervals.into_iter().enumerate() {
        // Distinct shelf/deep environments prevent ordinary same-environment merging.
        let amounts = ParcelAmounts { deposit: 100.0, class: 7 + index as u8, ..Default::default() };
        let total_after = stack.total_m() + 100.0;
        apply_parcel(&mut stack, &amounts, total_after, 1.0, &p, start, end, &mut log);
    }
    assert!((stack.total_m() - initial - 200.0).abs() < 1e-9);
    assert!((stack.total_m() - initial - log.ops_m3()).abs() < 1e-9);
    (stack, log)
}

#[test]
fn capacity_merge_is_not_an_invertible_deposition_record() {
    let (gap, gl) = history(16, true);
    let (continuous, cl) = history(16, false);
    assert_eq!(gap, continuous);
    assert_eq!(gap.cover.len(), 15);
    assert_eq!((gl.merges_cap, cl.merges_cap), (1, 1));
    assert_eq!((gl.merges_unconformity, cl.merges_unconformity), (1, 0));
    let layer = gap.cover.last().unwrap();
    assert_eq!((layer.age_base_my, layer.age_top_my, layer.thickness_m), (0.0, 10.0, 200.0));
    assert_ne!(layer.flags & STRATA_FLAG_MERGED, 0);
    println!(r#"PRES0_RESULT {{"case":"cap16_history_alias","columns_equal":true,"cover_layers":{},"gap_merge_count":{},"continuous_gap_merge_count":{},"retained_age_envelope_my":[{},{}],"retained_thickness_m":{}}}"#,
        gap.cover.len(), gl.merges_unconformity, cl.merges_unconformity,
        layer.age_base_my, layer.age_top_my, layer.thickness_m);
}

#[test]
fn larger_capacity_retains_this_fixture_distinction() {
    let (gap, gl) = history(64, true);
    let (continuous, cl) = history(64, false);
    assert_ne!(gap, continuous);
    assert_eq!(gap.cover.len(), 16);
    assert_eq!((gl.merges_cap, cl.merges_cap), (0, 0));
    let a = gap.cover[15].age_base_my - gap.cover[14].age_top_my;
    let b = continuous.cover[15].age_base_my - continuous.cover[14].age_top_my;
    assert_eq!((a, b), (6.0, 0.0));
    println!(r#"PRES0_RESULT {{"case":"cap64_retains_fixture","columns_equal":false,"cover_layers":{},"retained_gap_my":{},"continuous_gap_my":{}}}"#,
        gap.cover.len(), a, b);
}

#[test]
fn partial_erosion_marks_truncation_without_resolving_intralayer_age() {
    let (mut stack, _) = history(64, true);
    let previous = *stack.cover.last().unwrap();
    let before = stack.total_m();
    let mut log = StrataStepLog::default();
    apply_parcel(&mut stack, &ParcelAmounts { erosion: 50.0, ..Default::default() },
        before - 50.0, 1.0, &params(64), 10.0, 11.0, &mut log);
    let top = stack.cover.last().unwrap();
    assert_eq!(top.thickness_m, 50.0);
    assert_eq!((top.age_base_my, top.age_top_my), (previous.age_base_my, previous.age_top_my));
    assert_ne!(top.flags & STRATA_FLAG_TRUNCATED_TOP, 0);
    assert!((stack.total_m() - before + 50.0).abs() < 1e-9);
    assert!((stack.total_m() - before - log.ops_m3()).abs() < 1e-9);
    println!(r#"PRES0_RESULT {{"case":"partial_erosion_age_support","remaining_thickness_m":{},"truncated_flag":true,"unchanged_original_age_support_my":[{},{}],"removed_thickness_m":50}}"#,
        top.thickness_m, top.age_base_my, top.age_top_my);
}
