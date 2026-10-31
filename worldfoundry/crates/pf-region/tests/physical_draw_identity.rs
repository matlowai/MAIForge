//! Actual pinned-world regression. No synthetic terrain or default fixture path.
use pf_region::{skeleton::{self, Inputs}, t0::{DrawRealizationV1, SiteRecord, T0}, SkeletonParamsV1};
use serde_json::{json, Value};
use std::path::Path;

#[test]
fn physical_keys_and_verified_legacy_realization() {
    let Ok(config) = std::env::var("VZ_PHYSICAL_PAIR_JSON") else {
        eprintln!("SKIP actual-world key regression: set VZ_PHYSICAL_PAIR_JSON"); return;
    };
    let cfg: Value = serde_json::from_slice(&std::fs::read(config).unwrap()).unwrap();
    let out = Path::new(cfg["output"].as_str().unwrap());
    std::fs::create_dir_all(out).unwrap();
    let mut evidence = Vec::new();
    for pair in cfg["pairs"].as_array().unwrap() {
        let path = |key: &str| Path::new(pair[key].as_str().unwrap());
        let mut old = SiteRecord::load(path("off_site")).unwrap();
        old.physical_keying = true;
        let mut new = SiteRecord::load(path("on_site")).unwrap();
        new.physical_keying = true;
        let mut t0 = T0::load(path("off_run"), path("off_subgrid"), &old, 10.0).unwrap();
        let t1 = T0::load(path("on_run"), path("on_subgrid"), &new, 10.0).unwrap();
        let physical = skeleton::physical_upstream(&t0, &old).unwrap();
        assert_eq!(physical, skeleton::physical_upstream(&t1, &new).unwrap());
        assert_eq!(t0.subgrid_digest, t1.subgrid_digest);
        assert_ne!(old.bytes_digest, new.bytes_digest);
        assert_ne!(t0.terrain_digest, t1.terrain_digest);
        let recipe = |t, site| skeleton::recipe_for(&Inputs { t0: t, site, params: SkeletonParamsV1::default() }).unwrap();
        let rd = skeleton::draw_digest(&recipe(&t0, &old)).unwrap();
        assert_eq!(rd, skeleton::draw_digest(&recipe(&t1, &new)).unwrap());
        let mut metadata = new.clone();
        metadata.label.push_str(" metadata-only change");
        metadata.bytes_digest = format!("blake3:{}", "a".repeat(64));
        assert_eq!(rd, skeleton::draw_digest(&recipe(&t1, &metadata)).unwrap());
        let mut legacy = recipe(&t0, &old);
        legacy.upstream.retain(|(role, _)| !role.starts_with("physical."));
        let legacy_draw = skeleton::draw_digest(&legacy).unwrap();
        let mut draw = vec![("physical.site".into(), old.bytes_digest.clone()),
            ("physical.t0.terrain".into(), t0.terrain_digest.clone()),
            ("physical.t0.hydrology".into(), t0.hydrology_digest.clone())];
        draw.sort();
        let binding = DrawRealizationV1 { schema: "vz.physical-draw-realization/v1".into(), physical_upstream: physical.clone(), draw_upstream: draw };
        let mut bound = new.clone(); bound.draw_realization = Some(binding.clone());
        assert_eq!(legacy_draw, skeleton::draw_digest(&recipe(&t1, &bound)).unwrap());
        let mut changed_law = recipe(&t1, &bound);
        changed_law.params.planform_amplitude_permille += 1;
        changed_law.params_digest = pf_contracts::artifact::digest_canonical(&changed_law.params).unwrap();
        assert_ne!(legacy_draw, skeleton::draw_digest(&changed_law).unwrap());
        // A physical elevation change changes the default key and refuses a
        // retained realization. Neither behavior relies on a metadata digest.
        let c = old.t0_cell as usize;
        let original_elevation = t0.elevation_m[c];
        t0.elevation_m[c] += 1.0;
        assert_ne!(physical, skeleton::physical_upstream(&t0, &old).unwrap());
        assert!(skeleton::physical_upstream(&t0, &bound).unwrap_err().contains("physical inputs differ"));
        t0.elevation_m[c] = original_elevation;
        // Independent hydrology and site-location falsifiers.
        t0.discharge_m3_s[c] += 1.0;
        assert!(skeleton::physical_upstream(&t0, &bound).is_err());
        let mut moved = bound.clone(); moved.chart_origin_lat_e7 += 1;
        assert!(skeleton::physical_upstream(&t1, &moved).is_err());
        let mut wrong = bound.clone(); wrong.draw_realization.as_mut().unwrap().draw_upstream.pop();
        assert!(skeleton::physical_upstream(&t1, &wrong).is_err());
        let id = pair["id"].as_str().unwrap();
        for side in ["off", "on"] {
            let mut record: Value = serde_json::from_slice(&std::fs::read(path(&format!("{side}_site"))).unwrap()).unwrap();
            if pair["retain_legacy_realization"].as_bool().unwrap_or(false) { record["draw_realization"] = serde_json::to_value(&binding).unwrap(); }
            else { record["physical_keying"] = json!(true); }
            let dir = out.join(format!("{id}-{side}")); std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("site.json"), serde_json::to_vec_pretty(&record).unwrap()).unwrap();
        }
        evidence.push(json!({"id":id,"physical_upstream":physical,"default_physical_draw":rd,
            "retained_legacy_draw":legacy_draw,"binding":binding,"metadata_invariant":true,
            "elevation_discharge_location_falsifiers_refused":true,"malformed_binding_refused":true}));
    }
    std::fs::write(out.join("physical-key-proof.json"), serde_json::to_vec_pretty(&evidence).unwrap()).unwrap();
}

#[test]
fn actual_rivernet_law_provenance_is_not_a_random_seed() {
    use pf_region::{t1b::T0Extras, t1b_view::SkeletonView, t1b_rivernet::{rivernet_recipe, RiverNetParamsV1}};
    let Ok(config) = std::env::var("VZ_RIVERNET_KEY_PAIR_JSON") else {
        eprintln!("SKIP actual rivernet regression: set VZ_RIVERNET_KEY_PAIR_JSON"); return;
    };
    let cfg: Value = serde_json::from_slice(&std::fs::read(config).unwrap()).unwrap();
    let mut draw = Vec::new(); let mut full = Vec::new();
    for side in cfg["sides"].as_array().unwrap() {
        let path = |key: &str| Path::new(side[key].as_str().unwrap());
        let mut site = SiteRecord::load(path("site")).unwrap();
        site.physical_keying = true;
        let mut t0 = T0::load(path("run"), path("subgrid"), &site, 10.0).unwrap();
        let cp = pf_region::climate_runoff::ClimateRunoffParamsV1::v1();
        t0.attach_climate(path("climate"), cp.clone()).unwrap();
        let extras = T0Extras::load(path("run"), &t0).unwrap();
        let skel = SkeletonView::load(path("t1a")).unwrap();
        let mut p = RiverNetParamsV1::v1(false); p.climate_runoff = Some(cp);
        let mut r = rivernet_recipe(&t0, &extras, &site, &skel, None, &p).unwrap();
        let d = r.draw_digest().unwrap(); draw.push(d.clone()); full.push(r.digest().unwrap());
        r.laws.source.push_str(" metadata-only description correction");
        assert_eq!(d, r.draw_digest().unwrap(), "source prose must not re-key draws");
        r.laws.stream_power_k_per_my *= 2.0;
        assert_ne!(d, r.draw_digest().unwrap(), "physical law keying must remain sensitive");
    }
    assert_eq!(draw.len(), 2); assert_eq!(draw[0], draw[1]); assert_ne!(full[0], full[1]);
    std::fs::write(cfg["output"].as_str().unwrap(), serde_json::to_vec_pretty(&json!({"draw_digests":draw,"full_source_recipe_digests":full,"physical_law_falsifier":true})).unwrap()).unwrap();
}

#[test]
fn default_skeletons_match_retained_headers_curves_and_every_field() {
    let Ok(config) = std::env::var("VZ_LEGACY_KEY_FIXTURES_JSON") else {
        eprintln!("SKIP retained default-byte regression: set VZ_LEGACY_KEY_FIXTURES_JSON"); return;
    };
    let cfg: Value = serde_json::from_slice(&std::fs::read(config).unwrap()).unwrap();
    let mut proof = Vec::new();
    for f in cfg["fixtures"].as_array().unwrap() {
        let path = |key: &str| Path::new(f[key].as_str().unwrap());
        let site = SiteRecord::load(path("site")).unwrap();
        assert!(!site.physical_keying_enabled());
        let mut t0 = T0::load(path("run"), path("subgrid"), &site, 10.0).unwrap();
        let mut params = SkeletonParamsV1::default();
        if let Some(dir) = f["climate"].as_str() {
            let cp = pf_region::climate_runoff::ClimateRunoffParamsV1::v1();
            t0.attach_climate(Path::new(dir), cp.clone()).unwrap();
            params.runoff_rule = format!("climate_runoff_v1: R(x) x node area (params.climate_runoff; upstream t0.climate); r_ref = A-RAIN-UNITS {} m/yr", t0.runoff_m_per_yr);
            params.climate_runoff = Some(cp);
        }
        let sk = skeleton::build_with_twin(&Inputs { t0: &t0, site: &site, params: params.clone() }).unwrap();
        assert!(!sk.recipe.upstream.iter().any(|(role, _)| role.starts_with("physical.")));
        assert!(sk.physical_climate_twin_digest.is_none());
        if params.climate_runoff.is_none() {
            let mut reordered = sk.recipe.clone(); reordered.upstream.reverse();
            assert_eq!(skeleton::draw_digest(&reordered).unwrap(), reordered.digest().unwrap(), "legacy default must not normalize upstream ordering");
        }
        let (header, blobs) = pf_region::publish::header(&sk).unwrap();
        let digest = pf_contracts::artifact::digest_canonical(&header).unwrap();
        let published = pf_region::publish::Published { header, digest: digest.clone() };
        assert_eq!(serde_json::to_vec_pretty(&published).unwrap(), std::fs::read(path("reference").join("skeleton.json")).unwrap(), "complete published skeleton bytes: {}", f["name"]);
        for (name, bytes) in &blobs {
            assert_eq!(*bytes, std::fs::read(path("reference").join(name)).unwrap(), "field {name}: {}", f["name"]);
        }
        let checks = pf_region::checks::run(&sk, &t0).unwrap();
        if let Some(n) = f["expected_tributaries"].as_u64() {
            assert_eq!(checks.f2_strict_hillslope_nodes_ge_threshold as u64, n);
        }
        // The default regional constructor must also retain the legacy upstream
        // and climate twin, even though physical identities are available.
        let extras = pf_region::t1b::T0Extras::load(path("run"), &t0).unwrap();
        let view = pf_region::t1b_view::SkeletonView::load(path("reference")).unwrap();
        assert!(!view.physical_keying);
        let mut rp = pf_region::t1b_rivernet::RiverNetParamsV1::v1(false);
        rp.climate_runoff = params.climate_runoff;
        let recipe = pf_region::t1b_rivernet::rivernet_recipe(&t0, &extras, &site, &view, None, &rp).unwrap();
        assert!(!recipe.upstream.iter().any(|(role, _)| role.starts_with("physical.")));
        assert_eq!(recipe.t1a_draw_twin, view.climate_twin_digest);
        let mut requested = site.clone(); requested.physical_keying = true;
        assert!(pf_region::t1b_rivernet::rivernet_recipe(&t0, &extras, &requested, &view, None, &rp).unwrap_err().contains("keying modes differ"));
        proof.push(json!({"name":f["name"],"published_skeleton_digest":digest,"published_header_and_curves_bytes_equal":true,"field_count":blobs.len(),"all_field_bytes_equal":true,"tributaries":checks.f2_strict_hillslope_nodes_ge_threshold,"legacy_rivernet_draw_digest":recipe.draw_digest().unwrap(),"mixed_mode_refused":true}));
    }
    std::fs::write(cfg["output"].as_str().unwrap(), serde_json::to_vec_pretty(&proof).unwrap()).unwrap();
}
