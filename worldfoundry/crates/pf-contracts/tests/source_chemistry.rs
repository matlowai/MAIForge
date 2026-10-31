use pf_contracts::source_chemistry::*;
use pf_contracts::volcanic_supply::MassUg;
use std::collections::BTreeMap;

fn table() -> ChemistryTable {
    ChemistryTable::from_json(include_str!("../../pf-rockspec/data/ore_chemistry_v1.json")).unwrap()
}
fn mass(n: i128) -> MassUg {
    MassUg::new(n).unwrap()
}
fn request(species: &str, n: i128) -> PhaseChemistryRequest {
    PhaseChemistryRequest {
        species: species.into(),
        phase_mass: mass(n),
        realization_id: "test:declared-major-only/v1".into(),
        evidence_ref: "test:arithmetic-fixture-not-source-admission".into(),
        split: None,
    }
}
fn declared(species: &str, n: i128) -> PhaseChemistryRequest {
    let mut r = request(species, n);
    let fractions = match species {
        "clinopyroxene" => vec![("diopside", 8), ("hedenbergite", 2)],
        "olivine" => vec![("forsterite", 8), ("fayalite", 2)],
        "plagioclase" => vec![("anorthite", 6), ("albite", 4)],
        "clay" => vec![("kaolinite", 10)],
        "k_feldspar" => vec![("orthoclase", 10)],
        "quartz" => return r,
        _ => panic!("unsupported qualification species"),
    };
    r.split = Some(EndmemberSplit { basis: "mass_fraction".into(), label: "conditional ORE-1 demonstration only; not a measured B1/S2 assay; MASS not molar fractions".into(), evidence_ref: "ORE-1:explicit-declaration-no-stock-admission".into(), denominator: 10, endmembers: fractions.into_iter().map(|(e, numerator)| EndmemberFraction { endmember: e.into(), numerator }).collect() });
    r
}
#[test]
fn catalogue_formulas_close_at_zero_one_and_planet_scale() {
    let t = table();
    assert_eq!(t.endmembers.len(), 18);
    for e in &t.endmembers {
        for n in [
            0,
            1,
            2,
            17,
            23_039_943_219,
            210_210_497_978_861_646_921_009,
            i128::MAX,
        ] {
            let r = t.decompose(&request(&e.id, n)).unwrap();
            r.validate(&t).unwrap();
            r.composition.validate_mass(mass(n)).unwrap();
            let bytes = serde_json::to_vec(&r).unwrap();
            let decoded: StoichiometryReceipt = serde_json::from_slice(&bytes).unwrap();
            assert_eq!(r, decoded);
        }
    }
}
#[test]
fn groups_need_labelled_split_and_refuse_wrong_or_incomplete_members() {
    let t = table();
    for g in &t.groups {
        assert!(t
            .decompose(&request(&g.id, 1000))
            .unwrap_err()
            .contains("split"));
    }
    for species in [
        "clinopyroxene",
        "olivine",
        "plagioclase",
        "clay",
        "k_feldspar",
    ] {
        let r = declared(species, 17);
        t.decompose(&r).unwrap().validate(&t).unwrap();
        let mut bad = r.clone();
        bad.split.as_mut().unwrap().label.clear();
        assert!(t.decompose(&bad).is_err());
        let mut bad = r.clone();
        bad.split.as_mut().unwrap().basis = "mole_fraction".into();
        assert!(t.decompose(&bad).is_err());
        let mut bad = r.clone();
        bad.split.as_mut().unwrap().denominator += 1;
        assert!(t.decompose(&bad).is_err());
        let mut bad = r.clone();
        bad.split.as_mut().unwrap().endmembers[0].endmember = "chalcopyrite".into();
        assert!(t.decompose(&bad).is_err());
        let mut bad = r.clone();
        let e = bad.split.as_ref().unwrap().endmembers[0].clone();
        bad.split.as_mut().unwrap().endmembers.push(e);
        assert!(t.decompose(&bad).is_err());
    }
    assert!(t.decompose(&request("not_a_mineral", 123)).is_err());
}
#[test]
fn independent_quartz_formula_and_hydrous_sulfide_falsifiers() {
    let t = table();
    // One formula mass in micro-u gives integral elemental ratios in ug.
    let quartz = t.decompose(&request("quartz", 60_083_000)).unwrap();
    let values: BTreeMap<_, _> = quartz
        .composition
        .elements
        .iter()
        .map(|e| (e.element, e.mass.micrograms()))
        .collect();
    assert_eq!(
        values,
        BTreeMap::from([(Element::Si, 28_085_000), (Element::O, 31_998_000)])
    );
    for (name, required) in [
        ("gibbsite", Element::H),
        ("kaolinite", Element::H),
        ("goethite", Element::H),
        ("chalcopyrite", Element::S),
        ("pyrite", Element::S),
    ] {
        let r = t.decompose(&request(name, 1_000_000_000)).unwrap();
        assert!(r
            .composition
            .elements
            .iter()
            .any(|e| e.element == required && e.mass.micrograms() > 0));
        let mut bad = r.clone();
        bad.composition.elements.retain(|e| e.element != required);
        assert!(bad.validate(&t).is_err());
    }
    let organic = t
        .decompose(&request("organic_matter_declared_chons_v1", 100))
        .unwrap();
    assert_eq!(
        organic
            .composition
            .elements
            .iter()
            .find(|e| e.element == Element::C)
            .unwrap()
            .mass,
        mass(80)
    );
    assert_eq!(organic.composition.elements.len(), 5);
}
#[test]
fn receipts_detect_tampering_including_balanced_relabelling() {
    let t = table();
    let r = t.decompose(&declared("plagioclase", 123_456_789)).unwrap();
    let mut bad = r.clone();
    bad.components[0].allocation.terms[0].rounded_up ^= true;
    assert!(bad.validate(&t).is_err());
    let mut bad = r.clone();
    bad.components[0].allocation.terms[0].remainder_numerator += 1;
    assert!(bad.validate(&t).is_err());
    let mut bad = r.clone();
    bad.composition.elements[0].element = Element::Cu;
    assert!(bad.validate(&t).is_err());
    let mut changed = t.clone();
    changed.atomic_weights[0]
        .evidence
        .reference
        .push_str(" changed");
    assert!(r.validate(&changed).is_err());
    let mut bad = r.composition.clone();
    bad.elements.push(bad.elements[0].clone());
    assert!(bad.validate().is_err());
    let mut bad = r.composition.clone();
    bad.phase_mass = mass(1);
    assert!(bad.validate().is_err());
}
#[test]
fn allocator_matches_independent_small_integer_division_and_extreme_bounds() {
    for n in 0..1000i128 {
        let weights = BTreeMap::from([("a", 7u64), ("b", 11), ("c", 19)]);
        let r = allocate(mass(n), weights.clone()).unwrap();
        assert_eq!(
            r.terms.iter().map(|x| x.amount.micrograms()).sum::<i128>(),
            n
        );
        for term in &r.terms {
            let product = n * i128::from(weights[term.key]);
            assert_eq!(term.floor.micrograms(), product / 37);
            assert_eq!(term.remainder_numerator as i128, product % 37);
            assert!((term.amount.micrograms() * 37 - product).abs() < 37);
        }
    }
    let r = allocate(
        mass(i128::MAX),
        BTreeMap::from([("a", u64::MAX - 1), ("b", 1)]),
    )
    .unwrap();
    assert_eq!(
        r.terms[0].amount.checked_add(r.terms[1].amount).unwrap(),
        mass(i128::MAX)
    );
    assert!(allocate(mass(1), BTreeMap::from([("a", u64::MAX), ("b", 1)])).is_err());
    assert!(allocate(mass(1), BTreeMap::from([("a", 0)])).is_err());
    assert!(MassUg::new(-1).is_err());
    assert!(mass(i128::MAX).checked_add(mass(1)).is_err());
    let r = allocate(mass(1), BTreeMap::from([("z", 1), ("a", 1)])).unwrap();
    assert_eq!(r.terms[0].key, "a");
    assert!(r.terms[0].rounded_up);
}
#[test]
fn malformed_tables_and_oxide_basis_refuse() {
    let mut t = table();
    t.endmembers.push(t.endmembers[0].clone());
    assert!(t.validate().is_err());
    let mut t = table();
    t.atomic_weights[0].micro_u += 100;
    assert!(t.validate().is_err());
    let mut t = table();
    t.groups[0].allowed_endmembers = vec!["unknown".into()];
    assert!(t.validate().is_err());
    let oxide = OxideMass {
        formula: "Al2O3".into(),
        mass: mass(1),
        basis: OxideBasis::AnalyticalEquivalent,
        evidence_ref: "test:report-only-no-stock".into(),
    };
    oxide.validate().unwrap();
    let mut bad = oxide;
    bad.formula = "gibbsite".into();
    assert!(bad.validate().is_err());
}

fn small_json(path: &std::path::Path) -> serde_json::Value {
    assert!(std::fs::metadata(path).unwrap().len() < 200_000_000);
    serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap()
}
/// Real E4 and MINE-0 artifacts only. All conditional chemistry declarations are
/// retained in the emitted receipts; no fixture replaces a missing artifact.
#[test]
fn actual_e4_six_minerals_and_mine0_read_only_demo() {
    let Some(root) = std::env::var_os("VZ_ORE_E4_PACKETS") else {
        eprintln!("SKIP real E4/MINE-0 qualification: set VZ_ORE_E4_PACKETS, VZ_ORE_MINE_EVIDENCE, VZ_ORE_OUTPUT");
        return;
    };
    let mine_path = std::env::var_os("VZ_ORE_MINE_EVIDENCE").expect("mine evidence path required");
    let out = std::env::var_os("VZ_ORE_OUTPUT").expect("output directory required");
    let out = std::path::PathBuf::from(out);
    std::fs::create_dir_all(&out).unwrap();
    let t = table();
    let mut phase_receipts = Vec::new();
    let mut checked = 0;
    for site in ["B1", "S2"] {
        let path = std::path::Path::new(&root).join(format!("{site}.terrain-packet.json"));
        let packet = small_json(&path);
        let stock = &packet["stock"];
        assert_eq!(stock["initial"].as_object().unwrap().len(), 3);
        for (species, initial) in stock["initial"].as_object().unwrap() {
            let initial: MassUg = serde_json::from_value(initial.clone()).unwrap();
            let reserve: MassUg =
                serde_json::from_value(stock["unplaced_native_reserve"][species].clone()).unwrap();
            let mut allocated = MassUg::ZERO;
            for leaf in stock["leaves"].as_array().unwrap() {
                let n: MassUg = serde_json::from_value(leaf["stock"][species].clone()).unwrap();
                allocated = allocated.checked_add(n).unwrap();
            }
            assert_eq!(allocated.checked_add(reserve).unwrap(), initial);
            for (scope, n) in [
                ("initial", initial),
                ("reserve", reserve),
                ("allocated", allocated),
            ] {
                let mut req = declared(species, n.micrograms());
                req.realization_id = "ore-1:conditional-qualification-only/v1".into();
                req.evidence_ref = format!(
                    "{}#stock/{scope}/{species}; source_snapshot={}; NOT admitted chemistry",
                    path.display(),
                    stock["source_snapshot"]
                );
                let r = t.decompose(&req).unwrap();
                r.validate(&t).unwrap();
                r.composition.validate_mass(n).unwrap();
                phase_receipts.push(serde_json::json!({"site":site,"scope":scope,"receipt":r}));
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 6);
    let mine = small_json(std::path::Path::new(&mine_path));
    let removed = mine["after"]["edits"][0]["removed"].as_object().unwrap();
    assert_eq!(removed.len(), 3);
    let mut total = MassUg::ZERO;
    let mut elements = BTreeMap::<Element, MassUg>::new();
    let mut mine_receipts = Vec::new();
    for (species, n) in removed {
        let n: i128 = n.as_str().unwrap().parse().unwrap();
        total = total.checked_add(mass(n)).unwrap();
        let mut req = declared(species, n);
        req.realization_id = "ore-1:conditional-qualification-only/v1".into();
        req.evidence_ref = format!(
            "{}#after/edits/0/removed/{species}; NOT admitted chemistry",
            std::path::Path::new(&mine_path).display()
        );
        let r = t.decompose(&req).unwrap();
        r.validate(&t).unwrap();
        for e in &r.composition.elements {
            let entry = elements.entry(e.element).or_insert(MassUg::ZERO);
            *entry = entry.checked_add(e.mass).unwrap();
        }
        mine_receipts.push(r);
    }
    assert_eq!(total.micrograms(), 23_039_943_219);
    assert_eq!(
        elements
            .values()
            .try_fold(MassUg::ZERO, |s, m| s.checked_add(*m))
            .unwrap(),
        total
    );
    assert!(!elements.contains_key(&Element::Cu));
    assert_eq!(elements.len(), 7);
    let output = serde_json::json!({"qualification":"read_only_conditional_major_element_decomposition_not_assay_or_admission","e4_mineral_accounts":checked,"e4_phase_receipts":phase_receipts,"mine_evidence":std::path::Path::new(&mine_path),"mine_mass":total,"mine_elements":elements,"mine_receipts":mine_receipts,"rounding":"each allocation term error <1 ug; nested total per element error <2*N ug for N endmembers; no independent parent/child elemental-transfer conservation claim"});
    let path = out.join("chemistry-qualification.json");
    assert!(!path.exists(), "fresh output path required");
    std::fs::write(&path, serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    eprintln!("6/6 E4 minerals (18 initial/reserve/allocated receipts); 3/3 dig minerals; exact {} ug; {}",total.micrograms(),path.display());
}
