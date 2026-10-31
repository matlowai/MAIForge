//! WS-CONST (lane WC, 2026-10-12): legacy recipe preservation.
//!
//! GOLDEN DATA captured with the PRE-WS-CONST source (commit 6b8132d,
//! worldfoundry tree identical to 654ab03) by the ignored driver
//! pf-lab tests/ws_const_capture.rs `ws_const_capture_profiles`
//! (`generated-test-data/profiles.json`):
//! every scenario preset x fidelity as the contract-level profile and as
//! the pf-lab resolved profile under each physics option (legacy,
//! two-column, bodies). Each entry is the exact JSON the params digest was
//! taken over, and that digest.
//!
//! Asserted for every entry with the WS-CONST source:
//! - the old JSON loads (no pin, no new keys required);
//! - it re-serializes to the identical bytes, so `params_digest` is unchanged;
//! - `universe()` resolves and validates, and the resolution is the law
//!   profile the blocks imply (legacy / two-column / bodies).

use pf_contracts::artifact::digest_params;
use pf_contracts::laws::{CrustSeedingLawV1, IsostasyLawV1, SeaLevelLawV1};
use pf_contracts::params::PlanetParameterProfile;

/// (name, json, params digest) captured before WS-CONST.
const LEGACY_PROFILES: &[(&str, &str, &str)] = &[
    (
        "contract/toy/Preview",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:9c112dfe000114c9ae36d76e4de8ce7c8695bc174061d0a378337a0fdd86b543",
    ),
    (
        "preview/toy/Preview/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:5ce5cfc8eaa0ba74d785746ca4ca78a98bc66ea144fba44cc0e80a64e5ee4767",
    ),
    (
        "preview/toy/Preview/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:fd03b851fad5ae9bb688971da03c60e3b8c30969615d950c0db77c55d29ad36c",
    ),
    (
        "preview/toy/Preview/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:0e7be09b34217dbd341ead0c0f7f72476ae3785f73fb6b418731d61e82a2c83e",
    ),
    (
        "contract/toy/Standard",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:a7ffd7b4c7563f83d8dc5060bd40e5642de5fee0fe084584ad8f9b086b54234e",
    ),
    (
        "preview/toy/Standard/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:11da594f25fce5ea34c89d64765e28d3bb1b0f30bb3d770b3901973150d6afac",
    ),
    (
        "preview/toy/Standard/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:83ea533d993f928010c8ea1a417582de1bade0b588f9769f69f420fd3ddc394e",
    ),
    (
        "preview/toy/Standard/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:11a2798dc432a7e06d726c557b946495eeafe67d3cde246895deb20e4c1872c9",
    ),
    (
        "contract/toy/High",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:68c585b69bd74e06051fdc5249828a0ed0906874930f7f76f762bfbcac191565",
    ),
    (
        "preview/toy/High/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:46099b90fc8b446b8496c500bd13205c83582816fcfc87f027ca927a571e34be",
    ),
    (
        "preview/toy/High/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:4ad906a12993130a8780947751bc8b86b92a98dccca1e7bffa648cd72dffe9cc",
    ),
    (
        "preview/toy/High/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":7,"rate_scale_deg_per_my":1.0,"ocean_fraction":0.62,"continental_seeds":4,"hotspot_count":2,"duration_my":100.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:4bcd2c1367ae57676d6ae58490b77662eacb520b787e0d15cd065c3e6f1417c1",
    ),
    (
        "contract/benchmark-a/Preview",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:c1bea0b182448a800aca71d6087855c221a45a4f8fc1e4c23c8e39ba13b54358",
    ),
    (
        "preview/benchmark-a/Preview/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:794db88bae6caf719dc7854536b1b5ad61d135f14cc7b3679691639f1fa56b7c",
    ),
    (
        "preview/benchmark-a/Preview/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:adfccbd7bdb7416b68ce63f5ee083a00e5bdd1b1525739e7fa5227faf9de2a38",
    ),
    (
        "preview/benchmark-a/Preview/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:f5ff44dbb1204584f157cc777a2ba8dd00f9f9c1723b555a5a76ed7cebe836e0",
    ),
    (
        "contract/benchmark-a/Standard",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:b0ef73945ca08b912a43e008be49149ca2d1c32f4e00c508ac8401266050766c",
    ),
    (
        "preview/benchmark-a/Standard/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:ea5cfe89682ccb00cf18a102f62cf8d2429e65772998c1d9cc513095ef08aec8",
    ),
    (
        "preview/benchmark-a/Standard/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:ef8ac3dd1cb4b8c3a34efeed7dcbbed7dfd7d3b5e241fe167aeedef6451698f2",
    ),
    (
        "preview/benchmark-a/Standard/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:86acd12cd24d67c13d79407379de442941e236082c75c533eba18f2329ce60eb",
    ),
    (
        "contract/benchmark-a/High",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:c145b549082b514af82739d08da4c9745a6248dede32a336d7bef4d7f3da9167",
    ),
    (
        "preview/benchmark-a/High/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:609396df2129b4d647e0370eeea2dcfc8471ef11609c6edfb69c3f5c65fce766",
    ),
    (
        "preview/benchmark-a/High/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:fd5c05e63a873726f43ded665e5cd4d7c30c32df448f66051f322f5a90c39e0e",
    ),
    (
        "preview/benchmark-a/High/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":6,"rate_scale_deg_per_my":4.0,"ocean_fraction":0.35,"continental_seeds":5,"hotspot_count":1,"duration_my":300.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:a7fc6d4e6cdb7a8ebababe3c17f731adfdb1fad2df8a70a272db5e4f25a66ea2",
    ),
    (
        "contract/benchmark-b/Preview",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:647f02b289c609bdd4dcf2ab415e7813de4892ed43ced10a225fba2fe057dcf7",
    ),
    (
        "preview/benchmark-b/Preview/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:8b7dabeaaa1fb7bd5d35ece5729600bc7a03f1f76f8c7960166ed0075c530811",
    ),
    (
        "preview/benchmark-b/Preview/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:c0a61d9aee2fff7b77ad7ef76758065c210c367c256843092017dd00f8393149",
    ),
    (
        "preview/benchmark-b/Preview/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:654573ebb57ac3d4f946aaf73dceb8cbb87fad186a23c226bc72b861647ab249",
    ),
    (
        "contract/benchmark-b/Standard",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:a5d0d91837c697d76b54ffd1edfee07a14e4ac7de964d3a6bfd5d3aace360b36",
    ),
    (
        "preview/benchmark-b/Standard/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:8cd02022138ddf90cda36fbcd900eafc7e1516a51b432dd1f761f8903082891f",
    ),
    (
        "preview/benchmark-b/Standard/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:31764fe009afd0feb41ba3fb61edce30d839ec4731dfd71ca40a7f809d67eea2",
    ),
    (
        "preview/benchmark-b/Standard/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:3694aac96cd9ab69b0fe207da2d6abc92a347212cecfffbe09a55aa8a6627767",
    ),
    (
        "contract/benchmark-b/High",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:b090bdab6e51c9f8ba01bf7ef876514a1ac7b8e24ddb6ed809ddeab429c50dd7",
    ),
    (
        "preview/benchmark-b/High/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:31ac8ced0f53be65e68438c1ff355f773d724332e63d8d2e921289018b54250d",
    ),
    (
        "preview/benchmark-b/High/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:036ddc4f25466a3f280ed787458ef0b0c5bb757e94a521847d2bd216598eddd6",
    ),
    (
        "preview/benchmark-b/High/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":8,"rate_scale_deg_per_my":3.0,"ocean_fraction":0.8,"continental_seeds":2,"hotspot_count":3,"duration_my":200.0,"step_my":2.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:1b8a03c3449d8b685c9047bcad7013d0ef5867d3220c324b9d95c0ad23b08f3a",
    ),
    (
        "contract/benchmark-c/Preview",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:5b97be3a5a9bdd6f8e18d96f0c07893261f2f925232079fb1b059c2816494315",
    ),
    (
        "preview/benchmark-c/Preview/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:fc4507700d69955a480a27eb3a248a707c88bb0b0654c9628332c06e30ed3402",
    ),
    (
        "preview/benchmark-c/Preview/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:959f6eb27247b980fd08335468316534e4d4e7ec789399eeae5ab2156360fd40",
    ),
    (
        "preview/benchmark-c/Preview/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":0,"radial_layers":4,"substeps":1},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:65196f8e5e676c17c09150226ec105dc62121d4f2674e7f7eef7aa7c2a4a94cf",
    ),
    (
        "contract/benchmark-c/Standard",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:fdf76921c4e7268d5d0f191ab30f113c5b32d2219db8b6c202092aa05fba584d",
    ),
    (
        "preview/benchmark-c/Standard/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:b74d8037fc2f4519a507fcdec7ee7d3b8fef9faecaa2d11269f78253faa50d84",
    ),
    (
        "preview/benchmark-c/Standard/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:33f1c8d3a9ee41a292175375e5e1f7a6be53330ad92ef78bcb4f4a80e4340758",
    ),
    (
        "preview/benchmark-c/Standard/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":1,"radial_layers":8,"substeps":2},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:c3b5c16d02d6f2db47cb129a81566990d3494910d41460492d734e84bebfb049",
    ),
    (
        "contract/benchmark-c/High",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Flexure","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:a015c3d004e11986d0ab61ab90083908f273144342f099d219821a12737a5cd1",
    ),
    (
        "preview/benchmark-c/High/legacy",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0}}"#,
        "blake3:0086d554d2aaafb086b053506476a60d3cf6eb633613c2aad0eb43ba50aab25b",
    ),
    (
        "preview/benchmark-c/High/two_column",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0}}"#,
        "blake3:1417297810f7db933ceb6c4263412452c85e84dc3c09cc41f7532f3fcc8a4f9b",
    ),
    (
        "preview/benchmark-c/High/bodies",
        r#"{"planet":{"radius_m":6371000.0,"gravity_m_s2":9.81,"crust_density_kg_m3":2800.0,"mantle_density_kg_m3":3300.0,"sea_level_m":0.0},"tectonics":{"plate_count":5,"rate_scale_deg_per_my":0.5,"ocean_fraction":0.55,"continental_seeds":4,"hotspot_count":1,"duration_my":900.0,"step_my":5.0},"mantle":{"thermal_model":"Plate","compensation":"Airy","effective_elastic_thickness_km":40.0,"plume_heat_k_per_step":12.0},"surface":{"rainfall_mode":"Uniform","uniform_runoff_m_per_yr":1.0,"erodibility":2e-6},"fidelity":{"grid_level":2,"radial_layers":16,"substeps":4},"world_character":{"waterfall_wonder":0.0},"two_column_isostasy":{"oceanic_crust_density_kg_m3":2900.0,"seawater_density_kg_m3":1030.0,"compensation_depth_m":100000.0,"ocean_water_volume_m3":1.335e+18,"compensation_depth_margin_m":10000.0},"continental_bodies":{"rifted_margin_width_m":200000.0}}"#,
        "blake3:56c36efeda1d86e7a5d9fad4f96cc1f6a8e3a6190919ad517e82498a6bcd0c09",
    ),
];

#[test]
fn ws_const_old_profiles_load_reserialize_identically_and_keep_their_digest() {
    assert_eq!(LEGACY_PROFILES.len(), 48);
    for (name, json, digest) in LEGACY_PROFILES {
        let p: PlanetParameterProfile =
            serde_json::from_str(json).unwrap_or_else(|e| panic!("{name}: old JSON refused: {e}"));
        assert!(p.pinned_universe.is_none(), "{name}");
        let again = serde_json::to_string(&p).unwrap();
        assert_eq!(&again, json, "{name}: re-serialization moved");
        assert_eq!(&digest_params(&p).unwrap(), digest, "{name}: params digest moved");
        p.validate().unwrap_or_else(|e| panic!("{name}: {e}"));
    }
}

#[test]
fn ws_const_old_profiles_resolve_to_the_law_profile_their_blocks_imply() {
    let mut seen = [0usize; 3];
    for (name, json, _) in LEGACY_PROFILES {
        let p: PlanetParameterProfile = serde_json::from_str(json).unwrap();
        let u = p.universe().unwrap_or_else(|e| panic!("{name}: {e}"));
        let k = &u.constants;
        // The planet block values are read through unchanged (bitwise).
        assert_eq!(k.planet_radius_m.value.to_bits(), p.planet.radius_m.to_bits(), "{name}");
        assert_eq!(k.surface_gravity_m_s2.value.to_bits(), p.planet.gravity_m_s2.to_bits(), "{name}");
        assert_eq!(k.mantle_density_kg_m3.value.to_bits(), p.planet.mantle_density_kg_m3.to_bits(), "{name}");
        assert_eq!(k.continental_crust_density_kg_m3.value.to_bits(), p.planet.crust_density_kg_m3.to_bits(), "{name}");
        assert_eq!(k.sea_level_datum_m.value.to_bits(), p.planet.sea_level_m.to_bits(), "{name}");
        assert_eq!(k.stream_power_erodibility_per_yr.value.to_bits(), p.surface.erodibility.to_bits(), "{name}");
        if name.ends_with("/two_column") {
            seen[1] += 1;
            assert!(matches!(u.law_set.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }), "{name}");
            assert_eq!(u.law_set.sea_level, SeaLevelLawV1::VolumePrimaryGlobalFillV0, "{name}");
            assert_eq!(u.law_set.crust_seeding, CrustSeedingLawV1::WholePlateContinentsV0, "{name}");
            let b = p.two_column_isostasy.as_ref().unwrap();
            assert_eq!(k.ocean_water_volume_m3.as_ref().unwrap().value.to_bits(), b.ocean_water_volume_m3.to_bits());
        } else if name.ends_with("/bodies") {
            seen[2] += 1;
            assert!(matches!(u.law_set.isostasy, IsostasyLawV1::TwoColumnCommonDepthV0 { .. }), "{name}");
            match &u.law_set.crust_seeding {
                CrustSeedingLawV1::GeodesicDiscBodiesV0 { continental_area_fraction, .. } => {
                    let b = p.continental_bodies.as_ref().unwrap();
                    assert_eq!(continental_area_fraction.value.to_bits(), b.resolved_area_fraction(&p.tectonics).to_bits());
                }
                other => panic!("{name}: {other:?}"),
            }
        } else {
            seen[0] += 1;
            assert!(u.law_set.is_legacy(), "{name}: {:?}", u.law_set);
            assert!(k.ocean_water_volume_m3.is_none() && k.oceanic_crust_density_kg_m3.is_none(), "{name}");
        }
    }
    // 12 contract-level + 12 preview legacy, 12 two-column, 12 bodies.
    assert_eq!(seen, [24, 12, 12]);
}
