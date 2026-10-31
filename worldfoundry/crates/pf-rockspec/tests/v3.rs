use pf_rockspec::{v3::*, ConstituentFraction, MaterialRef, Property, SurfaceState};
fn c() -> Catalogue {
    Catalogue::embedded().unwrap()
}
fn context(setting: &str) -> FormationContext {
    FormationContext {
        setting: setting.into(),
        protolith: None,
        metamorphic_grade: None,
        burial_depth_m: None,
        provenance_ref: "declared-test-formation".into(),
    }
}
fn admission(key: &str, id: &str) -> BodyAdmission {
    BodyAdmission {
        material_ref: MaterialRef {
            spec_id: key.into(),
            spec_version: "pf-rockspec/3".into(),
            realization_id: id.into(),
            body_id: "test-body".into(),
        },
        kind: BodyKind::Host,
        event_id: "fixture-admission".into(),
        snapshot_id: "fixture-snapshot".into(),
        geometry_ref: "fixture-geometry".into(),
        stock_account_ref: "fixture-stock".into(),
        material_frame_ref: "fixture-frame".into(),
        host_body_id: None,
    }
}
fn body(c: &Catalogue, id: &str) -> ResolvedBody {
    let r = c.realization(id).unwrap();
    let a = admission(&r.lithofacies, id);
    c.resolve(a.material_ref.clone(), context(&r.settings[0]), a)
        .unwrap()
}
fn part(id: &str, x: f64) -> ConstituentFraction {
    ConstituentFraction {
        mineral: id.into(),
        mass_fraction: x,
    }
}
fn overlay() -> Overlay {
    let mut a = admission("granite", "granite.reference");
    a.kind = BodyKind::Vein;
    a.material_ref.body_id = "test-vein".into();
    a.host_body_id = Some("test-body".into());
    Overlay {
        id: "overlay1".into(),
        target_body_id: "test-body".into(),
        source_body: a,
        event_id: "fixture-admission".into(),
        stock_transfer_ref: "fixture-transfer".into(),
        removed: vec![part("quartz", 0.02)],
        added: vec![part("native_gold", 0.02)],
        overlay: OverlayKind::Mineralization {
            texture: MineralizationTexture::VeinFill,
            ore_mass_fractions: vec![part("native_gold", 0.02)],
        },
    }
}
fn y(b: &ResolvedBody) -> f64 {
    let s = b
        .surface_sample(
            &SurfaceState::fresh("face", "cut", [0.; 3], [0., 0., 1.]),
            2.,
        )
        .unwrap();
    s.base_reflectance_linear
        .into_iter()
        .zip([0.2126, 0.7152, 0.0722])
        .map(|(x, w)| x * w)
        .sum()
}
#[test]
fn catalogue_covers_facies_minerals_and_realizations() {
    let c = c();
    assert_eq!(c.core.rocks.len(), 49);
    assert_eq!(c.core.minerals.len(), 56);
    assert_eq!(c.realizations.len(), 67);
    for r in &c.realizations {
        let b = body(&c, &r.id);
        assert!(
            (b.bulk()
                .state
                .mineral_mode
                .iter()
                .map(|x| x.mass_fraction)
                .sum::<f64>()
                - 1.)
                .abs()
                < 1e-9
        );
        b.validate().unwrap();
        assert!(y(&b).is_finite());
    }
}
#[test]
fn unknown_incompatible_or_incomplete_body_is_rejected() {
    let c = c();
    let a = admission("granite", "granite.i_type_white");
    assert!(c
        .resolve(
            a.material_ref.clone(),
            context("continental_rift"),
            a.clone()
        )
        .is_err());
    let mut bad = a.clone();
    bad.geometry_ref.clear();
    assert!(c
        .resolve(bad.material_ref.clone(), context("continental_arc"), bad)
        .is_err());
    let mut bad = a;
    bad.material_ref.realization_id = "random-fragment".into();
    assert!(c
        .resolve(bad.material_ref.clone(), context("continental_arc"), bad)
        .is_err());
    assert!(c
        .core
        .realize_prior(
            "granite",
            admission("granite", "granite.reference").material_ref
        )
        .is_err());
}
#[test]
fn persist_and_fragment_reads_preserve_realization_and_ore() {
    let c = c();
    let mut b = body(&c, "granite.i_type_white");
    b.apply_overlay(overlay()).unwrap();
    let text = serde_json::to_string(&b).unwrap();
    let loaded = ResolvedBody::from_json(&text).unwrap();
    assert_eq!(loaded.reference(), b.reference());
    assert_eq!(loaded.bulk().state, b.bulk().state);
    for id in ["fragment-a", "fragment-b"] {
        let out = loaded
            .surface_sample(
                &SurfaceState::fresh(id, "new-cut", [0.; 3], [0., 0., 1.]),
                2.,
            )
            .unwrap();
        assert_eq!(&out.material_ref, b.reference());
    }
    assert_eq!(text, serde_json::to_string(&b).unwrap());
    let mut broken: serde_json::Value = serde_json::from_str(&text).unwrap();
    broken["bulk"]["state"]["mineral_mode"][0]["mass_fraction"] = serde_json::json!(0.99);
    assert!(ResolvedBody::from_json(&broken.to_string()).is_err());
}
#[test]
fn modal_and_optical_realizations_are_distinct() {
    let c = c();
    let pink = body(&c, "granite.i_type_pink");
    let white = body(&c, "granite.i_type_white");
    assert!(y(&white) > y(&pink));
    assert_ne!(
        body(&c, "granite.s_type_white").bulk().state.mineral_mode,
        white.bulk().state.mineral_mode
    );
    assert_ne!(
        body(&c, "andesite.arc_hornblende")
            .bulk()
            .state
            .mineral_mode,
        body(&c, "andesite.two_pyroxene").bulk().state.mineral_mode
    );
    assert!(y(&body(&c, "chert.white")) > y(&body(&c, "chert.flint")) + 0.1);
    assert!(y(&body(&c, "ignimbrite.unwelded")) > y(&body(&c, "ignimbrite.welded")) + 0.1);
    assert!(y(&body(&c, "shale.grey")) > y(&body(&c, "shale.organic_rich")));
    assert!(
        (y(&body(&c, "basalt_tholeiitic.lip")) - y(&body(&c, "basalt_alkali.island_cap"))).abs()
            > 0.005
    );
}
#[test]
fn penetrative_fabric_is_not_an_open_joint() {
    let c = c();
    for r in &c.realizations {
        if matches!(r.spec.fabric.as_str(), "cleaved" | "foliated") {
            assert!(r.fabric.penetrative);
            assert_eq!(r.fabric.aperture_m, 0.);
        }
        assert!(r
            .spec
            .joint_sets
            .iter()
            .all(|j| j.relation != "fabric_parallel"));
    }
    let mud = c.realization("mudstone.reference").unwrap();
    assert_eq!(mud.fabric.kind, FabricKind::Bedding);
    assert!(!mud.fabric.penetrative);
    assert!(mud
        .spec
        .joint_sets
        .iter()
        .any(|j| j.relation == "bedding_discontinuity"));
    let mut broken = c;
    broken.realizations[0].fabric.aperture_m = 0.0003;
    assert!(broken.validate().is_err());
}
#[test]
fn joints_are_family_specific_not_fabric_string_specific() {
    let c = c();
    let spacing = |k: &str| {
        c.realization(&format!("{k}.reference"))
            .unwrap()
            .spec
            .joint_sets[0]
            .spacing_m
            .value
    };
    assert_ne!(spacing("granite"), spacing("chert"));
    assert_ne!(spacing("chert"), spacing("rock_salt"));
    assert_ne!(spacing("rock_salt"), spacing("marble"));
}
#[test]
fn family_mechanics_fix_outliers_and_declare_creep() {
    let c = c();
    let val = |key: &str, n: &str| {
        let r = c.realization(&format!("{key}.reference")).unwrap();
        match &r.spec.mechanical[n] {
            Property::KnownPrior { estimate } => estimate.value,
            _ => panic!(),
        }
    };
    assert_eq!(val("rock_salt", "young_modulus_pa"), 25e9);
    assert_eq!(val("gypsum_rock", "young_modulus_pa"), 30e9);
    assert_eq!(val("marble", "young_modulus_pa"), 58.5e9);
    assert_eq!(val("obsidian", "young_modulus_pa"), 69.6e9);
    assert_eq!(val("obsidian", "fracture_toughness_pa_sqrt_m"), 0.8e6);
    assert!((val("shale", "friction_angle_rad").to_degrees() - 21.).abs() < 1e-9);
    for k in ["rock_salt", "gypsum_rock", "shale"] {
        assert!(
            c.realization(&format!("{k}.reference"))
                .unwrap()
                .mechanics
                .creep
        );
    }
}
#[test]
fn all_inventory_ids_map_without_silent_group_resolution() {
    let c = c();
    let ids = [
        "silicic_glass",
        "quartz",
        "k_feldspar",
        "plagioclase",
        "biotite",
        "clinopyroxene",
        "olivine",
        "amphibole",
        "calcite",
        "dolomite",
        "clay",
        "muscovite",
        "halite",
        "gypsum",
    ];
    assert_eq!(c.inventory_mapping.len(), ids.len());
    for id in ids {
        let m = c
            .inventory_mapping
            .iter()
            .find(|m| m.inventory_id == id)
            .unwrap();
        let mut t = InventoryVolumeTerm {
            inventory_id: id.into(),
            solid_volume_fraction: 1.,
            source_density_kg_m3: 2700.,
            species_mass_split: None,
        };
        assert_eq!(
            c.inventory_volume_to_mass(&[t.clone()]).is_err(),
            m.requires_species_resolution
        );
        t.species_mass_split = Some(vec![part(&m.rockspec_ids[0], 1.)]);
        assert!(c.inventory_volume_to_mass(&[t]).is_ok());
    }
}
#[test]
fn volume_to_mass_preserves_source_density_and_closure() {
    let c = c();
    let ts = vec![
        InventoryVolumeTerm {
            inventory_id: "quartz".into(),
            solid_volume_fraction: 0.5,
            source_density_kg_m3: 2650.,
            species_mass_split: None,
        },
        InventoryVolumeTerm {
            inventory_id: "clay".into(),
            solid_volume_fraction: 0.5,
            source_density_kg_m3: 2600.,
            species_mass_split: Some(vec![part("illite", 0.7), part("smectite", 0.3)]),
        },
    ];
    let r = c.inventory_volume_to_mass(&ts).unwrap();
    assert_eq!(r.source_grain_density_kg_m3, 2625.);
    assert!(
        (r.mass_mode
            .iter()
            .find(|m| m.mineral == "quartz")
            .unwrap()
            .mass_fraction
            - 2650. / 5250.)
            .abs()
            < 1e-12
    );
    let mut bad = ts;
    bad[0].solid_volume_fraction = 0.6;
    assert!(c.inventory_volume_to_mass(&bad).is_err());
    bad[0].solid_volume_fraction = 0.5;
    bad[1].species_mass_split = Some(vec![part("native_gold", 1.)]);
    assert!(c.inventory_volume_to_mass(&bad).is_err());
}
#[test]
fn all_requested_mining_minerals_and_body_kinds_serialize() {
    let c = c();
    for id in [
        "galena",
        "sphalerite",
        "cassiterite",
        "chromite",
        "native_gold",
        "uraninite",
        "fluorite",
        "barite",
        "graphite",
        "gibbsite",
        "boehmite",
        "diaspore",
        "bornite",
        "chalcocite",
        "native_silver",
        "native_copper",
        "molybdenite",
        "siderite",
        "malachite",
        "azurite",
        "scheelite",
        "beryl",
        "spodumene",
        "pyrite",
        "magnetite",
    ] {
        assert!(c.core.mineral(id).is_ok());
        assert!(c.mineral_metadata.contains_key(id));
    }
    for kind in [
        BodyKind::Vein,
        BodyKind::Dyke,
        BodyKind::Sill,
        BodyKind::Fault,
        BodyKind::ShearZone,
        BodyKind::Breccia,
        BodyKind::ContactAureole,
        BodyKind::Placer,
        BodyKind::Nodule,
    ] {
        let s = serde_json::to_string(&kind).unwrap();
        assert_eq!(kind, serde_json::from_str(&s).unwrap());
    }
}
#[test]
fn alteration_types_replace_mass_in_the_same_bulk_state() {
    let c = c();
    for alteration in [
        AlterationKind::Potassic,
        AlterationKind::Phyllic,
        AlterationKind::Argillic,
        AlterationKind::Propylitic,
        AlterationKind::Weathering,
    ] {
        let (mineral, host) = match alteration {
            AlterationKind::Potassic => ("k_feldspar", "plagioclase"),
            AlterationKind::Phyllic => ("muscovite", "k_feldspar"),
            AlterationKind::Argillic => ("kaolinite", "k_feldspar"),
            AlterationKind::Propylitic => ("chlorite", "biotite"),
            _ => ("goethite", "biotite"),
        };
        let mut b = body(&c, "granite.reference");
        let mut o = overlay();
        o.removed = vec![part(host, 0.02)];
        o.added = vec![part(mineral, 0.02)];
        o.overlay = OverlayKind::Alteration {
            alteration,
            intensity: 0.02,
        };
        b.apply_overlay(o).unwrap();
        b.validate().unwrap();
        assert_eq!(b.bulk().overlays.len(), 1);
    }
}
#[test]
fn ore_overlay_is_conservative_and_requires_provenance() {
    let c = c();
    let mut b = body(&c, "granite.reference");
    let before = y(&b);
    b.apply_overlay(overlay()).unwrap();
    assert_ne!(y(&b), before);
    b.validate().unwrap();
    assert!(b.apply_overlay(overlay()).is_err());
    for i in 0..6 {
        let mut b = body(&c, "granite.reference");
        let pristine = serde_json::to_string(&b).unwrap();
        let mut o = overlay();
        match i {
            0 => o.stock_transfer_ref.clear(),
            1 => o.source_body.geometry_ref.clear(),
            2 => o.added[0].mass_fraction = 0.03,
            3 => o.removed[0] = part("native_gold", 0.02),
            4 => o.source_body.host_body_id = None,
            _ => o.event_id = "unadmitted-event".into(),
        };
        assert!(b.apply_overlay(o).is_err());
        assert_eq!(pristine, serde_json::to_string(&b).unwrap());
    }
}
#[test]
fn sandstone_coatings_and_cement_are_subsets_not_biotite_hacks() {
    let c = c();
    for key in ["sandstone_quartz.reference", "arkose.reference"] {
        let r = c.realization(key).unwrap();
        assert!(!r.spec.mode.iter().any(|m| m.mineral == "biotite"));
        let s = r.parameters.sedimentary.as_ref().unwrap();
        assert_eq!(s.grain_coatings[0].mineral, "goethite");
    }
}
#[test]
fn f8_optics_priors_and_plume_split_are_explicit() {
    let c = c();
    let m = c.core.mineral("magnetite").unwrap();
    let f = pf_rockspec::fresnel_f0(
        m.refractive_index.value,
        m.extinction_k.as_ref().unwrap()[0].value,
        1.,
    );
    assert!((0.20..0.23).contains(&f));
    assert!(c.core.mineral("plagioclase").unwrap().diffuse_reflectance[0].value > 0.7);
    assert!(c.core.mineral("olivine").unwrap().diffuse_reflectance[1].value > 0.3);
    assert_eq!(c.plume_priors.len(), 3);
    assert!(
        c.plume_priors[0].weights["basalt_tholeiitic"].value
            > c.plume_priors[2].weights["basalt_tholeiitic"].value
    );
}
#[test]
fn populated_overlay_assemblages_are_closed_and_require_admission() {
    let c = c();
    assert_eq!(c.alteration_assemblages.len(), 5);
    assert_eq!(c.mineralization_assemblages.len(), 4);
    for a in &c.alteration_assemblages {
        let mut b = body(&c, "granite.reference");
        let mut o = overlay();
        o.added = a
            .product_mode
            .iter()
            .map(|m| part(&m.mineral, m.mass_fraction.value * 0.02))
            .collect();
        o.overlay = OverlayKind::Alteration {
            alteration: a.kind.clone(),
            intensity: 0.02,
        };
        b.apply_overlay(o).unwrap();
        b.validate().unwrap();
    }
    let mut b = body(&c, "granite.reference");
    b.apply_overlay(overlay()).unwrap();
    let mut duplicate = overlay();
    duplicate.id = "different-label-same-transfer".into();
    assert!(b.apply_overlay(duplicate).is_err());
}

#[test]
fn invalid_modes_parameters_and_mechanics_are_rejected() {
    let original = c();
    let mut bad = original.clone();
    bad.realizations[0].spec.mode[0].mass_fraction.value += 0.1;
    assert!(bad.validate().is_err());
    let mut bad = original.clone();
    let phase = bad.realizations[0].spec.mode[0].clone();
    bad.realizations[0].spec.mode.push(phase);
    assert!(bad.validate().is_err());
    let mut bad = original.clone();
    bad.realizations[0].parameters.plagioclase_an_fraction = Some(1.2);
    assert!(bad.validate().is_err());
    let mut bad = original;
    if let Property::KnownPrior { estimate } = &mut bad.realizations[0]
        .spec
        .mechanical
        .get_mut("young_modulus_pa")
        .unwrap()
    {
        estimate.unit = "MPa".into();
    }
    assert!(bad.validate().is_err());
}
