//! Manufactured contract truth only: exact small integer stock makes failures
//! attributable to accounting/identity, not a claim about any C7 volcano.
use pf_contracts::{emplacement::*, planet::*, volcanic_supply::*};
use serde_json::json;
fn known<T>(value: T) -> Knowledge<T> {
    Knowledge::Known { value }
}
fn unknown<T>() -> Knowledge<T> {
    Knowledge::Unknown {
        reason: "not measured in manufactured contract truth".into(),
    }
}
fn mass(v: i128) -> MassUg {
    MassUg::new(v).unwrap()
}
fn time(a: f64, b: f64) -> TimeInterval {
    TimeInterval {
        start_my: a,
        end_my: b,
    }
}
fn world() -> WorldScope {
    WorldScope {
        seed: 17,
        generation_pin: format!("blake3:{}", "1".repeat(64)),
    }
}
fn lattice() -> PlanetLattice {
    PlanetLattice {
        kind: PLANET_KIND.into(),
        rows: 180,
        cols: 360,
        subdivisions: 1,
    }
}
fn support() -> NativeSupport {
    NativeSupport {
        radius_m: 6_371_000.,
        denominator: 4,
        cells: vec![
            SupportCell {
                cell: PlanetKey::Cell { cell: 100 },
                whole_cell_area_m2: 100.,
                represented_area_m2: 50.,
                weight: 1,
            },
            SupportCell {
                cell: PlanetKey::Cell { cell: 101 },
                whole_cell_area_m2: 200.,
                represented_area_m2: 100.,
                weight: 3,
            },
        ],
        selection_rule: "full manufactured source support; weights not area ratios".into(),
        evidence_ref: "manufactured:test-support".into(),
    }
}
fn density() -> Knowledge<DensityBasis> {
    known(DensityBasis {
        id: "manufactured:legacy-2500".into(),
        kind: DensityKind::LegacyUniformDry,
        density_kg_m3: EvidenceValue {
            value: 2500.,
            unit: "kg/m3".into(),
            interval: EvidenceInterval::Bounds {
                lower: 2500.,
                upper: 2500.,
                interpretation: "exact manufactured test truth".into(),
            },
            conditions: "dry solid fixture, not a calibrated volcanic material".into(),
            reference: "manufactured:vs1-fixture".into(),
            evidence_kind: "manufactured_test_truth".into(),
        },
        intrinsic_pore_fraction: known(0.),
        interparticle_pore_fraction: known(0.),
    })
}
fn partition(n: i128) -> ProductionPartition {
    ProductionPartition {
        gross_potential: known(mass(n + 20)),
        admitted: known(mass(n)),
        capped_or_rejected: known(mass(10)),
        unextracted: known(mass(10)),
    }
}
fn fixture() -> SupplyBook {
    let identity = SourceIdentity {
        world: world(),
        anchor: PlanetKey::Body { body: 41 },
        lineage: "manufactured:arc-front-7/slab-9".into(),
        interval: time(10., 12.),
    };
    let source = MeltSupply {
        id: identity.id().unwrap(),
        identity,
        model_version: "manufactured:source/v1".into(),
        origin: SupplyOrigin::BoundaryInflux {
            law_version: "manufactured:boundary-law/v1".into(),
            evidence_ref: "manufactured:captured-admission".into(),
        },
        support: support(),
        dry_rock: partition(1000),
        volatile: partition(100),
        density_basis: density(),
        material: unknown(),
        provenance_refs: vec!["manufactured:source-event".into()],
        quantity_evidence_refs: vec!["manufactured:exact-ug-no-conversion".into()],
        constants: vec![],
        composition: None,
    };
    let account = |id: &str, component, phase, initial, closing| StockAccount {
        id: id.into(),
        source_id: source.id.clone(),
        component,
        phase,
        origin_episode_id: unknown(),
        material: unknown(),
        density_basis: density(),
        support: support(),
        initial: known(mass(initial)),
        closing: known(mass(closing)),
    };
    let mut accounts = vec![
        account(
            "reservoir.dry",
            Component::DryRock,
            StockPhase::Melt,
            1000,
            650,
        ),
        account("lava", Component::DryRock, StockPhase::SolidRock, 0, 200),
        account(
            "intrusion",
            Component::DryRock,
            StockPhase::SolidRock,
            0,
            50,
        ),
        account("ash", Component::DryRock, StockPhase::AirborneSolid, 0, 60),
        account("ocean", Component::DryRock, StockPhase::SolidRock, 0, 40),
        account(
            "reservoir.volatile",
            Component::Volatile,
            StockPhase::DissolvedVolatile,
            100,
            90,
        ),
        account(
            "atmosphere",
            Component::Volatile,
            StockPhase::ExsolvedVolatile,
            0,
            10,
        ),
    ];
    let mut episode = EruptionEpisode {
        id: String::new(),
        world: world(),
        anchor: PlanetKey::Body { body: 41 },
        reservoir_id: "reservoir".into(),
        vent_id: "vent/7".into(),
        vent: GeoPoint::new(0, 180_000_000_000).unwrap(),
        parent_event_ref: "manufactured:coarse-event".into(),
        interval: time(10.5, 11.),
        timing: TimingStatus::Conditional {
            observed_interval: time(10., 12.),
            realization_id: "manufactured:timing/v1".into(),
            assumption_ref: "manufactured:conditional-coarse-time".into(),
        },
        allocations: vec![],
        emplacement_law_version: "manufactured:emplacement/v1".into(),
        placement_support: support(),
    };
    for (component, branch, n) in [
        (Component::DryRock, EruptionBranch::Lava, 200),
        (Component::DryRock, EruptionBranch::Intrusion, 50),
        (Component::DryRock, EruptionBranch::AirborneAsh, 100),
        (Component::Volatile, EruptionBranch::VolatileRelease, 10),
    ] {
        episode.allocations.push(EpisodeAllocation {
            source_id: source.id.clone(),
            component,
            branch,
            amount: mass(n),
        });
    }
    episode.id = episode.expected_id().unwrap();
    for a in &mut accounts {
        if !a.id.starts_with("reservoir.") {
            a.origin_episode_id = known(episode.id.clone());
        }
    }
    let edge = |id: &str, from: &str, to: &str, component, amount, revision, purpose| {
        EmplacementTransfer {
            id: id.into(),
            world: world(),
            source_id: source.id.clone(),
            from_account: from.into(),
            to_account: to.into(),
            component,
            amount: mass(amount),
            interval: time(10.5, 11.),
            purpose,
            density_basis: density(),
            recipient_support: support(),
            input_revision: revision,
            committed_revision: revision + 1,
            evidence_refs: vec!["manufactured:edge".into()],
        }
    };
    let erupted = |branch| TransferPurpose::Eruption {
        episode_id: episode.id.clone(),
        branch,
    };
    let transfers = vec![
        edge(
            "t1",
            "reservoir.dry",
            "lava",
            Component::DryRock,
            200,
            0,
            erupted(EruptionBranch::Lava),
        ),
        edge(
            "t2",
            "reservoir.dry",
            "intrusion",
            Component::DryRock,
            50,
            1,
            erupted(EruptionBranch::Intrusion),
        ),
        edge(
            "t3",
            "reservoir.dry",
            "ash",
            Component::DryRock,
            100,
            2,
            erupted(EruptionBranch::AirborneAsh),
        ),
        edge(
            "t4",
            "ash",
            "ocean",
            Component::DryRock,
            40,
            3,
            TransferPurpose::Redistribution {
                evidence_ref: "manufactured:ash-settling".into(),
            },
        ),
        edge(
            "t5",
            "reservoir.volatile",
            "atmosphere",
            Component::Volatile,
            10,
            4,
            erupted(EruptionBranch::VolatileRelease),
        ),
    ];
    SupplyBook {
        schema: pf_contracts::volcanic_supply::SCHEMA.into(),
        world: world(),
        lattice: lattice(),
        initial_revision: 0,
        committed_revision: 5,
        sources: vec![source],
        reservoirs: vec![MagmaReservoir {
            id: "reservoir".into(),
            account_ids: vec!["reservoir.dry".into(), "reservoir.volatile".into()],
            temperature_k: unknown(),
            volatile_model: unknown(),
            composition: None,
        }],
        accounts,
        episodes: vec![episode],
        transfers,
    }
}
fn refuses(book: &SupplyBook, reason: &str) {
    let err = book.validate().unwrap_err();
    assert!(err.contains(reason), "expected {reason:?}, got {err:?}");
}

#[test]
fn conservative_multi_branch_receipt_and_json_round_trip() {
    let b = fixture();
    b.validate().unwrap();
    let encoded = serde_json::to_string(&b).unwrap();
    let decoded: SupplyBook = serde_json::from_str(&encoded).unwrap();
    assert_eq!(decoded, b);
    decoded.validate().unwrap();
    assert_eq!(
        b.accounts
            .iter()
            .filter(|a| a.component == Component::DryRock)
            .map(|a| a.closing.require().unwrap().micrograms())
            .sum::<i128>(),
        1000
    );
    assert_eq!(
        b.accounts
            .iter()
            .filter(|a| a.component == Component::Volatile)
            .map(|a| a.closing.require().unwrap().micrograms())
            .sum::<i128>(),
        100
    );
}
#[test]
fn exact_micrograms_are_strings_with_unit_and_version() {
    let m = mass(i128::MAX);
    let v = serde_json::to_value(m).unwrap();
    assert_eq!(
        v,
        json!({"unit":"ug","version":1,"decimal":i128::MAX.to_string()})
    );
    assert_eq!(serde_json::from_value::<MassUg>(v).unwrap(), m);
    assert!(m.checked_add(mass(1)).is_err());
    assert!(mass(0).checked_sub(mass(1)).is_err());
    for value in [
        "-1",
        "+1",
        "01",
        "1.0",
        "170141183460469231731687303715884105728",
    ] {
        assert!(
            serde_json::from_value::<MassUg>(json!({"unit":"ug","version":1,"decimal":value}))
                .is_err()
        );
    }
    for (version, unit, value) in [
        (1, "kg", json!("1")),
        (2, "ug", json!("1")),
        (1, "ug", json!(1)),
    ] {
        assert!(serde_json::from_value::<MassUg>(
            json!({"unit":unit,"version":version,"decimal":value})
        )
        .is_err());
    }
}
#[test]
fn negative_stock_refused_at_wire_and_constructor() {
    assert!(MassUg::new(-1).is_err());
    let mut v = serde_json::to_value(fixture()).unwrap();
    v["accounts"][0]["initial"]["value"]["decimal"] = json!("-1");
    assert!(serde_json::from_value::<SupplyBook>(v).is_err());
}
#[test]
fn duplicate_transfer_is_not_idempotently_spent_again() {
    let mut b = fixture();
    b.transfers[1].id = b.transfers[0].id.clone();
    refuses(&b, "duplicate transfer");
}
#[test]
fn duplicate_or_renamed_source_cannot_mint_stock() {
    let mut b = fixture();
    b.sources.push(b.sources[0].clone());
    refuses(&b, "duplicate source");
    b.sources[1].id = "renamed".into();
    refuses(&b, "source identity mismatch");
}
#[test]
fn unknown_fields_are_explicit_required_and_never_default() {
    let mut v = serde_json::to_value(fixture()).unwrap();
    v["sources"][0].as_object_mut().unwrap().remove("material");
    assert!(serde_json::from_value::<SupplyBook>(v).is_err());
    let mut v = serde_json::to_value(fixture()).unwrap();
    v["sources"][0]["material"] = serde_json::Value::Null;
    assert!(serde_json::from_value::<SupplyBook>(v).is_err());
    let mut v = serde_json::to_value(fixture()).unwrap();
    v["transfers"][0]["stock_from_mesh"] = json!(42);
    assert!(serde_json::from_value::<SupplyBook>(v).is_err());
    let mut b = fixture();
    b.sources[0].material = Knowledge::Unknown {
        reason: String::new(),
    };
    refuses(&b, "empty identity");
}
#[test]
fn wrong_world_pin_refused_at_source_episode_and_transfer() {
    for which in 0..3 {
        let mut b = fixture();
        let scope = match which {
            0 => &mut b.sources[0].identity.world,
            1 => &mut b.episodes[0].world,
            _ => &mut b.transfers[0].world,
        };
        scope.seed += 1;
        refuses(&b, "wrong world pin");
    }
    let mut b = fixture();
    b.transfers[0].world.generation_pin = format!("blake3:{}", "2".repeat(64));
    refuses(&b, "wrong world pin");
}
#[test]
fn density_change_requires_reconciliation_even_with_same_numeric_density() {
    let mut b = fixture();
    if let Knowledge::Known { value } = &mut b.accounts[1].density_basis {
        value.kind = DensityKind::DepositBulkDry;
    }
    refuses(&b, "incompatible density basis");
    let mut b = fixture();
    b.transfers[0].density_basis = unknown();
    refuses(&b, "incompatible density basis");
}
#[test]
fn zero_nan_negative_density_and_porosity_refused() {
    for n in [0., -1., f64::NAN, f64::INFINITY] {
        let mut d = density();
        if let Knowledge::Known { value } = &mut d {
            value.density_kg_m3.value = n;
        }
        assert!(d.validate_with(DensityBasis::validate).is_err());
    }
    let mut d = density();
    if let Knowledge::Known { value } = &mut d {
        value.interparticle_pore_fraction = known(1.);
    }
    assert!(d.validate_with(DensityBasis::validate).is_err());
}
#[test]
fn source_ids_are_chart_independent_and_world_lineage_time_bound() {
    let b = fixture();
    let s = &b.sources[0].identity;
    let a = Atlas {
        kind: ATLAS_KIND.into(),
        radius_mm: 6_371_000_000,
    };
    let q = GeoPoint::new(0, 0).unwrap();
    let owner = a.owner(q).unwrap();
    assert_ne!(
        a.project(owner, q).unwrap().chart,
        a.project(owner + 1, q).unwrap().chart
    );
    assert_eq!(s.id().unwrap(), s.clone().id().unwrap());
    for which in 0..4 {
        let mut other = s.clone();
        match which {
            0 => other.world.seed += 1,
            1 => other.lineage.push('x'),
            2 => other.interval.end_my += 1.,
            _ => other.anchor = PlanetKey::Body { body: 42 },
        };
        assert_ne!(s.id().unwrap(), other.id().unwrap());
    }
    let mut zero = s.clone();
    zero.interval.start_my = 0.;
    let mut neg = zero.clone();
    neg.interval.start_my = -0.;
    assert_eq!(zero.id().unwrap(), neg.id().unwrap());
    let encoded = serde_json::to_string(s).unwrap();
    assert!(!encoded.contains("chart"));
}
#[test]
fn incomplete_duplicate_or_noncell_support_refused() {
    let mut b = fixture();
    b.sources[0].support.denominator += 1;
    refuses(&b, "incomplete support");
    let mut b = fixture();
    b.sources[0].support.cells[1].cell = b.sources[0].support.cells[0].cell.clone();
    refuses(&b, "duplicate support");
    let mut b = fixture();
    b.sources[0].support.cells[0].cell = PlanetKey::Body { body: 1 };
    refuses(&b, "support cell");
}
#[test]
fn potential_and_rejected_production_cannot_fund_stock() {
    let mut b = fixture();
    b.accounts[0].initial = known(mass(1020));
    refuses(&b, "opening stock");
    let mut b = fixture();
    b.sources[0].dry_rock.gross_potential = known(mass(1000));
    refuses(&b, "production partition");
}
#[test]
fn unknown_supply_is_preserved_but_cannot_be_spent() {
    let mut b = fixture();
    b.sources[0].dry_rock = ProductionPartition {
        gross_potential: unknown(),
        admitted: unknown(),
        capped_or_rejected: unknown(),
        unextracted: unknown(),
    };
    b.accounts[0].initial = unknown();
    b.accounts[0].closing = unknown();
    b.episodes.clear();
    b.transfers.clear();
    b.committed_revision = 0;
    for a in &mut b.accounts {
        a.closing = a.initial.clone();
        a.origin_episode_id = unknown();
    }
    b.validate().unwrap();
    b.episodes = fixture().episodes;
    refuses(&b, "unknown quantity");
}
#[test]
fn all_source_accounts_and_reservoir_ownership_are_unique() {
    let mut b = fixture();
    b.accounts.push(b.accounts[0].clone());
    refuses(&b, "duplicate account");
    let mut b = fixture();
    let mut r = b.reservoirs[0].clone();
    r.id = "other".into();
    b.reservoirs.push(r);
    refuses(&b, "two reservoirs");
}
#[test]
fn retained_rock_is_not_fresh_eruptible_magma() {
    let mut b = fixture();
    b.sources[0].origin = SupplyOrigin::RetainedImport {
        body_account: "old-rock".into(),
        debit_receipt_ref: "captured".into(),
    };
    refuses(&b, "retained import");
    let mut b = fixture();
    b.accounts[0].phase = StockPhase::SolidRock;
    refuses(&b, "not eruptible");
}
#[test]
fn oversubscribed_episode_reservations_refuse() {
    let mut b = fixture();
    b.episodes[0].allocations[0].amount = mass(1001);
    refuses(&b, "allocations exceed");
}
#[test]
fn transfer_over_allocation_refuses() {
    let mut b = fixture();
    b.transfers[0].amount = mass(201);
    refuses(&b, "allocation exceeded");
}
#[test]
fn negative_balance_from_secondary_transfer_refuses_without_mutation() {
    let mut b = fixture();
    b.transfers[3].amount = mass(101);
    let before = b.clone();
    refuses(&b, "negative stock");
    assert_eq!(b, before);
}
#[test]
fn dry_and_volatile_branches_cannot_substitute_for_each_other() {
    let mut b = fixture();
    b.transfers[0].component = Component::Volatile;
    refuses(&b, "source/component");
    let mut b = fixture();
    b.episodes[0].allocations[0].component = Component::Volatile;
    refuses(&b, "episode allocation");
    let mut b = fixture();
    b.accounts[0].phase = StockPhase::DissolvedVolatile;
    refuses(&b, "phase mismatch");
}
#[test]
fn no_second_eruption_disguised_as_redistribution() {
    let mut b = fixture();
    b.transfers[0].purpose = TransferPurpose::Redistribution {
        evidence_ref: "lie".into(),
    };
    refuses(&b, "redistribution changes eruption provenance");
}
#[test]
fn event_time_bounds_and_nonfinite_times_refuse() {
    let mut b = fixture();
    b.transfers[0].interval.end_my = 11.5;
    refuses(&b, "debit/recipient/time");
    let mut b = fixture();
    b.episodes[0].timing = TimingStatus::Conditional {
        observed_interval: time(11., 12.),
        realization_id: "r".into(),
        assumption_ref: "a".into(),
    };
    refuses(&b, "outside observed");
    assert!(time(f64::NAN, 12.).validate().is_err());
    assert!(time(12., 10.).validate().is_err());
}
#[test]
fn wrong_source_and_unqualified_recipient_refuse() {
    let mut b = fixture();
    b.transfers[0].source_id = "wrong".into();
    refuses(&b, "source/component");
    let mut b = fixture();
    b.transfers[0].to_account = "missing".into();
    refuses(&b, "unknown credit");
    let mut b = fixture();
    b.transfers[0].recipient_support.evidence_ref = "different".into();
    refuses(&b, "support does not match");
}
#[test]
fn stale_replayed_revision_and_wrong_closing_stock_refuse() {
    let mut b = fixture();
    b.transfers[1].input_revision = 0;
    refuses(&b, "nonsequential");
    let mut b = fixture();
    b.accounts[1].closing = known(mass(201));
    refuses(&b, "closing account");
    let mut b = fixture();
    b.committed_revision += 1;
    refuses(&b, "closing revision");
}
#[test]
fn literature_record_requires_units_conditions_reference_and_honest_interval() {
    let mut c = EvidenceValue {
        value: 70.,
        unit: "km3".into(),
        interval: EvidenceInterval::NotReported {
            reason: "Table 2 gives no uncertainty interval".into(),
        },
        conditions: "Mayon edifice above interpreted basal surface, not erupted DRE".into(),
        reference: "https://doi.org/10.1016/j.jvolgeores.2021.107251 Table 2".into(),
        evidence_kind: "published_morphometric_estimate; record only, not adopted".into(),
    };
    c.validate().unwrap();
    c.reference.clear();
    assert!(c.validate().is_err());
    c.reference = "ref".into();
    c.interval = EvidenceInterval::Bounds {
        lower: 71.,
        upper: 80.,
        interpretation: "wrong".into(),
    };
    assert!(c.validate().is_err());
}
fn deposit(b: &SupplyBook) -> DepositBody {
    DepositBody {
        id: "lava-body".into(),
        world: world(),
        account_ids: vec!["lava".into()],
        episode_ids: vec![b.episodes[0].id.clone()],
        transfer_ids: vec!["t1".into()],
        material: unknown(),
        formation_interval: time(10., 12.),
        occupied_support: support(),
        density_basis: density(),
        geometry_ref: unknown(),
        contact_refs: unknown(),
        cooling_history_ref: unknown(),
        exposure_history_ref: unknown(),
    }
}
#[test]
fn deposit_view_inherits_provenance_without_minting_stock() {
    let b = fixture();
    let d = deposit(&b);
    d.validate(&b).unwrap();
    let copy = d.clone();
    copy.validate(&b).unwrap();
    assert_eq!(b.accounts[1].closing, known(mass(200)));
    let mut bad = d.clone();
    bad.transfer_ids = vec!["t2".into()];
    assert!(bad.validate(&b).is_err());
    let mut bad = d;
    bad.account_ids = vec!["reservoir.dry".into()];
    assert!(bad.validate(&b).is_err());
}
#[test]
fn existing_body_input_wire_contract_is_unchanged() {
    let old = BodyInput {
        schema: pf_contracts::emplacement::SCHEMA.into(),
        manifest_sha256: "a".repeat(64),
    };
    old.validate().unwrap();
    assert_eq!(
        serde_json::to_value(old).unwrap(),
        json!({"schema":"vz.native-body-source/v1","manifest_sha256":"a".repeat(64)})
    );
}

#[test]
fn transported_lots_cannot_erase_eruption_identity() {
    let mut b = fixture();
    b.accounts[4].origin_episode_id = unknown();
    refuses(&b, "changes eruption provenance");
    let mut b = fixture();
    b.accounts[1].origin_episode_id = unknown();
    refuses(&b, "unknown quantity");
    let mut b = fixture();
    b.accounts[1].initial = known(mass(1));
    refuses(&b, "must originate");
}
#[test]
fn eruption_placement_cannot_silently_truncate_or_move_support() {
    let mut b = fixture();
    b.episodes[0].placement_support.cells[1].cell = PlanetKey::Cell { cell: 102 };
    refuses(&b, "outside declared placement support");
    let mut b = fixture();
    b.accounts[1].support.radius_m *= 2.;
    refuses(&b, "radius mismatch");
}
#[test]
fn reservoir_crystallization_changes_phase_without_changing_dry_mass() {
    let mut b = fixture();
    let mut crystals = b.accounts[0].clone();
    crystals.id = "crystals".into();
    crystals.phase = StockPhase::Crystal;
    crystals.initial = known(mass(0));
    crystals.closing = known(mass(25));
    b.accounts[0].closing = known(mass(625));
    b.accounts.push(crystals);
    b.reservoirs[0].account_ids.push("crystals".into());
    let mut t = b.transfers[0].clone();
    t.id = "crystallize".into();
    t.to_account = "crystals".into();
    t.amount = mass(25);
    t.input_revision = 5;
    t.committed_revision = 6;
    t.purpose = TransferPurpose::ReservoirProcess {
        law_version: "manufactured:crystallization".into(),
    };
    b.transfers.push(t);
    b.committed_revision = 6;
    b.validate().unwrap();
}
#[test]
fn secondary_transport_cannot_happen_before_its_source_exists() {
    let mut b = fixture();
    b.transfers[3].interval = time(9., 9.5);
    refuses(&b, "precedes source");
}

#[test]
fn captured_later_eruption_can_spend_stored_magma_but_not_before_production() {
    let mut b = fixture();
    b.episodes[0].interval = time(13., 14.);
    b.episodes[0].timing = TimingStatus::Captured {
        event_ref: "manufactured:later-observed-eruption".into(),
    };
    b.episodes[0].id = b.episodes[0].expected_id().unwrap();
    let id = b.episodes[0].id.clone();
    for a in &mut b.accounts {
        if matches!(a.origin_episode_id, Knowledge::Known { .. }) {
            a.origin_episode_id = known(id.clone());
        }
    }
    for t in &mut b.transfers {
        t.interval = time(13., 14.);
        if let TransferPurpose::Eruption { episode_id, .. } = &mut t.purpose {
            *episode_id = id.clone();
        }
    }
    b.validate().unwrap();
    b.episodes[0].timing = TimingStatus::Conditional {
        observed_interval: time(13., 14.),
        realization_id: "invented-date".into(),
        assumption_ref: "wrong-coarse-interval".into(),
    };
    refuses(&b, "timing evidence outside");
    b.episodes[0].timing = TimingStatus::Captured {
        event_ref: "manufactured:impossible-early-event".into(),
    };
    b.episodes[0].interval = time(9., 10.);
    b.episodes[0].id = b.episodes[0].expected_id().unwrap();
    refuses(&b, "precedes source interval");
}
#[test]
fn degassing_can_release_volatiles_without_erupting_dry_rock() {
    let mut b = fixture();
    b.episodes[0]
        .allocations
        .retain(|a| a.component == Component::Volatile);
    b.transfers.retain(|t| t.component == Component::Volatile);
    b.transfers[0].input_revision = 0;
    b.transfers[0].committed_revision = 1;
    b.committed_revision = 1;
    for a in &mut b.accounts {
        if a.component == Component::DryRock {
            a.closing = a.initial.clone();
        }
    }
    b.validate().unwrap();
    assert_eq!(b.accounts[0].closing, known(mass(1000)));
}
#[test]
fn one_reservoir_and_episode_keep_two_source_budgets_distinct() {
    let mut b = fixture();
    let mut second = b.sources[0].clone();
    second.identity.lineage = "manufactured:plume-4".into();
    second.id = second.identity.id().unwrap();
    let source_id = second.id.clone();
    b.sources.push(second);
    let original_accounts = b.accounts.clone();
    for mut a in original_accounts {
        a.id.push_str(".second");
        a.source_id = source_id.clone();
        if a.id.starts_with("reservoir.") {
            b.reservoirs[0].account_ids.push(a.id.clone());
        }
        b.accounts.push(a);
    }
    let mut extra_allocations = b.episodes[0].allocations.clone();
    for a in &mut extra_allocations {
        a.source_id = source_id.clone();
    }
    b.episodes[0].allocations.extend(extra_allocations);
    let mut extra_transfers = b.transfers.clone();
    for t in &mut extra_transfers {
        t.id.push_str(".second");
        t.source_id = source_id.clone();
        t.from_account.push_str(".second");
        t.to_account.push_str(".second");
        t.input_revision += 5;
        t.committed_revision += 5;
    }
    b.transfers.extend(extra_transfers);
    b.committed_revision = 10;
    b.validate().unwrap();
    let mut wrong_sphere = b.clone();
    wrong_sphere.sources[1].support.radius_m *= 2.;
    refuses(&wrong_sphere, "source/planet radius mismatch");
    // Crossing source accounts must fail even if both sources have ample stock.
    b.transfers[0].to_account.push_str(".second");
    refuses(&b, "source/component");
}
#[test]
fn exact_transfer_can_preserve_unknown_density_without_a_volume_claim() {
    let mut b = fixture();
    b.sources[0].density_basis = unknown();
    for a in &mut b.accounts {
        a.density_basis = unknown();
    }
    for t in &mut b.transfers {
        t.density_basis = unknown();
    }
    b.validate().unwrap();
    assert_eq!(b.accounts[1].closing, known(mass(200)));
}

#[test]
fn optional_chemistry_snapshot_is_legacy_wire_neutral_and_mass_bound() {
    use pf_contracts::source_chemistry::{Element, ElementComposition, ElementMass};
    let mut b = fixture();
    let legacy = serde_json::to_string(&b).unwrap();
    assert!(!legacy.contains("composition"));
    let roundtrip: SupplyBook = serde_json::from_str(&legacy).unwrap();
    assert_eq!(serde_json::to_string(&roundtrip).unwrap(), legacy);
    let composition = |n| ElementComposition {
        schema: pf_contracts::source_chemistry::SCHEMA.into(),
        realization_id: "manufactured:snapshot-only".into(),
        evidence_ref: "manufactured:bulk-closure-falsifier-not-geology".into(),
        phase_mass: mass(n),
        elements: vec![ElementMass {
            element: Element::O,
            mass: mass(n),
        }],
    };
    b.sources[0].composition = Some(composition(1100));
    b.reservoirs[0].composition = Some(composition(650 + 90));
    b.validate().unwrap();
    b.sources[0].composition = Some(composition(1099));
    refuses(&b, "composition snapshot mass mismatch");
    b.sources[0].composition = Some(composition(1100));
    b.reservoirs[0].composition = Some(composition(739));
    refuses(&b, "composition snapshot mass mismatch");
    b.reservoirs[0].composition = None;
    b.sources[0].dry_rock.admitted = unknown();
    refuses(&b, "unknown quantity cannot fund transfer");
}
